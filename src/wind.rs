use sounding_base::{Profile, Sounding};

use crate::error::*;
use crate::layers::{self, Layer};
use crate::levels::height_level;
use itertools::Itertools;
use metfor::spd_dir_to_uv;

/// Calculate the mean wind in a layer.
///
/// This is not the pressure weighted mean.
///
/// Returns `(direction_deg, speed_kts)`
pub fn mean_wind(layer: &Layer, snd: &Sounding) -> Result<(f64, f64)> {
    use std::f64;

    let height = snd.get_profile(Profile::GeopotentialHeight);
    let wind_dir = snd.get_profile(Profile::WindDirection);
    let wind_spd = snd.get_profile(Profile::WindSpeed);

    let max_hgt = if layer.top.height.is_some() {
        layer.top.height.unpack()
    } else {
        return Err(AnalysisError::MissingProfile);
    };

    let min_hgt = if layer.bottom.height.is_some() {
        layer.bottom.height.unpack()
    } else {
        return Err(AnalysisError::MissingProfile);
    };

    let (_, _, _, mut iu, mut iv, dz) = izip!(height, wind_dir, wind_spd)
        .filter_map(|(hgt, dir, spd)| {
            if hgt.is_some() && dir.is_some() && spd.is_some() {
                Some((hgt.unpack(), dir.unpack(), spd.unpack()))
            } else {
                None
            }
        })
        .skip_while(|&(hgt, _, _)| hgt < min_hgt)
        .take_while(|&(hgt, _, _)| hgt <= max_hgt)
        .map(|(hgt, dir, spd)| {
            let (u, v) = metfor::spd_dir_to_uv(dir, spd);
            (hgt, u, v)
        })
        .fold((f64::MAX, 0.0, 0.0, 0.0, 0.0, 0.0), |acc, (hgt, u, v)| {
            let (old_hgt, old_u, old_v, mut iu, mut iv, mut acc_dz) = acc;

            let dz = hgt - old_hgt;

            if dz < 0.0 {
                return (hgt, u, v, iu, iv, acc_dz);
            }
            iu += (u + old_u) * dz;
            iv += (v + old_v) * dz;
            acc_dz += dz;

            (hgt, u, v, iu, iv, acc_dz)
        });

    if dz == 0.0 {
        return Err(AnalysisError::NotEnoughData);
    }

    iu /= 2.0 * dz;
    iv /= 2.0 * dz;

    Ok(metfor::uv_to_spd_dir(iu, iv))
}

/// Storm relative helicity.
#[doc(hidden)]
pub fn sr_helicity(layer: &Layer, storm_motion_uv_ms: (f64, f64), snd: &Sounding) -> Result<f64> {
    let height = snd.get_profile(Profile::GeopotentialHeight);
    let spd = snd.get_profile(Profile::WindSpeed);
    let dir = snd.get_profile(Profile::WindDirection);

    let bottom = layer.bottom.height.ok_or(AnalysisError::MissingValue)?;
    let top = layer.top.height.ok_or(AnalysisError::MissingValue)?;

    let vals: Vec<(f64, f64, f64)> = izip!(height, spd, dir)
        .filter_map(|(h, s, d)| {
            if let (Some(h), Some(s), Some(d)) = (h.into_option(), s.into_option(), d.into_option())
            {
                Some((h, s, d))
            } else {
                None
            }
        })
        .skip_while(|(h, _, _)| *h < bottom)
        .take_while(|(h, _, _)| *h < top)
        .map(|(h, s, d)| {
            let (u, v) = metfor::spd_dir_to_uv(d, s);
            (h, u - storm_motion_uv_ms.0, v - storm_motion_uv_ms.1)
        })
        .collect();

    if vals.len() < 3 {
        return Err(AnalysisError::NotEnoughData);
    }

    let dz0 = vals[1].0 - vals[0].0;
    let du0 = (vals[1].1 - vals[0].1) / dz0;
    let dv0 = (vals[1].2 - vals[0].2) / dz0;

    let last_idx = vals.len() - 1;
    let dzf = vals[last_idx].0 - vals[last_idx - 1].0;
    let duf = (vals[last_idx].1 - vals[last_idx - 1].1) / dzf;
    let dvf = (vals[last_idx].2 - vals[last_idx - 1].2) / dzf;
    let uf = vals[last_idx].1;
    let vf = vals[last_idx].2;

    let mut derivatives: Vec<(f64, f64, f64, f64, f64)> = Vec::with_capacity(vals.len());
    derivatives.push((vals[0].1, vals[0].2, du0, dv0, dz0));

    for ((h0, u0, v0), (h1, u1, v1), (h2, u2, v2)) in vals.into_iter().tuple_windows::<(_, _, _)>()
    {
        let dz = h1 - h0;
        let du = (u2 - u0) / (h2 - h0);
        let dv = (v2 - v0) / (h2 - h0);

        derivatives.push((u1, v1, du, dv, dz));
    }
    derivatives.push((uf, vf, duf, dvf, dzf));

    unimplemented!()
}

/// Calculate the super cell storm motion using the "id" method.
///
/// Returns the storm motions in m/s of the right and left mover cells.
/// ((rm_u_ms, rm_v_ms), (lm_u_ms, lm_v_ms))
pub fn bunkers_storm_motion(snd: &Sounding) -> Result<((f64, f64), (f64, f64))> {
    let layer = &layers::layer_agl(snd, 6000.0)?;
    let mean_wind = mean_wind(layer, snd)?;
    let (mean_u, mean_v) = spd_dir_to_uv(mean_wind.0, mean_wind.1);

    let (shear_u, shear_v) = bulk_shear_half_km(layer, snd)?;
    const D: f64 = 7.5; // m/s

    let scale = D / shear_u.hypot(shear_v);
    let (delta_u, delta_v) = (shear_v * scale, -shear_u * scale);

    Ok((
        (mean_u + delta_u, mean_v + delta_v),
        (mean_u - delta_u, mean_v - delta_v),
    ))
}

/// Calculate the bulk shear of a layer using winds averaged over the bottom and top half km.
///
/// When using the id method for storm motion vectors, the bulk shear was calculated with top and
/// bottom wind vectors that were averaged over top/bottom half km of the layer.
///
/// Returns `(shear_u_ms, shear_v_ms)`
pub(crate) fn bulk_shear_half_km(layer: &Layer, snd: &Sounding) -> Result<(f64, f64)> {
    let bottom = layer
        .bottom
        .height
        .into_option()
        .ok_or(AnalysisError::MissingValue)?;
    let top = layer
        .top
        .height
        .into_option()
        .ok_or(AnalysisError::MissingValue)?;

    // abort of not at least 250 meters of non-overlapping area.
    if top - bottom < 750.0 {
        return Err(AnalysisError::NotEnoughData);
    }

    let top_bottom_layer = height_level(bottom + 500.0, snd)?;
    let bottom_layer = &Layer {
        top: top_bottom_layer,
        bottom: layer.bottom,
    };
    let (bottom_mean_dir, bottom_mean_spd) = mean_wind(bottom_layer, snd)?;
    let (bottom_u, bottom_v) = spd_dir_to_uv(bottom_mean_dir, bottom_mean_spd);

    let bottom_top_layer = height_level(top - 500.0, snd)?;
    let top_layer = &Layer {
        top: layer.top,
        bottom: bottom_top_layer,
    };
    let (top_mean_dir, top_mean_spd) = mean_wind(top_layer, snd)?;
    let (top_u, top_v) = spd_dir_to_uv(top_mean_dir, top_mean_spd);

    Ok((top_u - bottom_u, top_v - bottom_v))
}
