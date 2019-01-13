use metfor::{Meters, MetersPSec, Quantity, WindUV};
use sounding_base::Sounding;

use crate::error::*;
use crate::layers::{self, Layer};
use crate::levels::height_level;
use itertools::Itertools;

/// Calculate the mean wind in a layer.
///
/// This is not the pressure weighted mean.
///
pub fn mean_wind(layer: &Layer, snd: &Sounding) -> Result<WindUV<MetersPSec>> {
    use std::f64;

    let height = snd.height_profile();
    let wind = snd.wind_profile();

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

    let (_, _, _, mut iu, mut iv, dz) = izip!(height, wind)
        .filter_map(|(hgt, wind)| {
            if hgt.is_some() && wind.is_some() {
                Some((hgt.unpack(), wind.unpack()))
            } else {
                None
            }
        })
        .skip_while(|&(hgt, _)| hgt < min_hgt)
        .take_while(|&(hgt, _)| hgt <= max_hgt)
        .map(|(hgt, wind)| {
            let WindUV { u, v } = WindUV::<MetersPSec>::from(wind);
            (hgt, u, v)
        })
        .fold(
            (
                Meters(f64::MAX),
                MetersPSec(0.0),
                MetersPSec(0.0),
                MetersPSec(0.0),
                MetersPSec(0.0),
                Meters(0.0),
            ),
            |acc, (hgt, u, v)| {
                let (old_hgt, old_u, old_v, mut iu, mut iv, mut acc_dz) = acc;

                let dz = hgt - old_hgt;

                if dz < Meters(0.0) {
                    return (hgt, u, v, iu, iv, acc_dz);
                }
                iu += (u + old_u) * dz.unpack();
                iv += (v + old_v) * dz.unpack();
                acc_dz += dz;

                (hgt, u, v, iu, iv, acc_dz)
            },
        );

    if dz.unpack().abs() < std::f64::EPSILON {
        return Err(AnalysisError::NotEnoughData);
    }

    iu /= 2.0 * dz.unpack();
    iv /= 2.0 * dz.unpack();

    Ok(WindUV { u: iu, v: iv })
}

/// Storm relative helicity.
#[doc(hidden)]
pub fn sr_helicity(
    layer: &Layer,
    storm_motion_uv_ms: (MetersPSec, MetersPSec),
    snd: &Sounding,
) -> Result<f64> {
    let height = snd.height_profile();
    let wind = snd.wind_profile();

    let bottom = layer.bottom.height.ok_or(AnalysisError::MissingValue)?;
    let top = layer.top.height.ok_or(AnalysisError::MissingValue)?;

    let vals: Vec<(f64, f64, f64)> = izip!(height, wind)
        .filter_map(|(h, w)| {
            if let (Some(h), Some(w)) = (h.into_option(), w.into_option()) {
                Some((h, w))
            } else {
                None
            }
        })
        .skip_while(|(h, _)| *h < bottom)
        .take_while(|(h, _)| *h < top)
        .map(|(h, w)| {
            let WindUV { u, v } = WindUV::<MetersPSec>::from(w);
            (
                h.unpack(),
                (u - storm_motion_uv_ms.0).unpack(),
                (v - storm_motion_uv_ms.1).unpack(),
            )
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
/// (right mover, left mover)
pub fn bunkers_storm_motion(snd: &Sounding) -> Result<(WindUV<MetersPSec>, WindUV<MetersPSec>)> {
    let layer = &layers::layer_agl(snd, Meters(6000.0))?;
    let WindUV {
        u: mean_u,
        v: mean_v,
    } = mean_wind(layer, snd)?;

    let WindUV {
        u: shear_u,
        v: shear_v,
    } = bulk_shear_half_km(layer, snd)?;
    const D: f64 = 7.5; // m/s

    let scale = D / shear_u.unpack().hypot(shear_v.unpack());
    let (delta_u, delta_v) = (shear_v * scale, -shear_u * scale);

    Ok((
        WindUV {
            u: mean_u + delta_u,
            v: mean_v + delta_v,
        },
        WindUV {
            u: mean_u - delta_u,
            v: mean_v - delta_v,
        },
    ))
}

/// Calculate the bulk shear of a layer using winds averaged over the bottom and top half km.
///
/// When using the id method for storm motion vectors, the bulk shear was calculated with top and
/// bottom wind vectors that were averaged over top/bottom half km of the layer.
///
/// Returns `(shear_u_ms, shear_v_ms)`
pub(crate) fn bulk_shear_half_km(layer: &Layer, snd: &Sounding) -> Result<WindUV<MetersPSec>> {
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
    if top - bottom < Meters(750.0) {
        return Err(AnalysisError::NotEnoughData);
    }

    let top_bottom_layer = height_level(bottom + Meters(500.0), snd)?;
    let bottom_layer = &Layer {
        top: top_bottom_layer,
        bottom: layer.bottom,
    };

    let WindUV {
        u: bottom_u,
        v: bottom_v,
    } = mean_wind(bottom_layer, snd)?;

    let bottom_top_layer = height_level(top - Meters(500.0), snd)?;
    let top_layer = &Layer {
        top: layer.top,
        bottom: bottom_top_layer,
    };

    let WindUV { u: top_u, v: top_v } = mean_wind(top_layer, snd)?;

    Ok(WindUV {
        u: top_u - bottom_u,
        v: top_v - bottom_v,
    })
}
