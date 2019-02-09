use metfor::{IntHelicityM2pS2, Knots, Meters, MetersPSec, Quantity, WindSpdDir, WindUV};
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

    let max_hgt = layer.top.height.ok_or(AnalysisError::MissingProfile)?;
    let min_hgt = layer.bottom.height.ok_or(AnalysisError::MissingProfile)?;

    let (old_hgt, old_u, old_v, mut iu, mut iv, dz) = izip!(height, wind)
        // Filter out missing values
        .filter_map(|(hgt, wind)| hgt.into_option().and_then(|h| wind.map(|w| (h,w))))
        // Skip values below the layer
        .skip_while(|&(hgt, _)| hgt < min_hgt)
        // Only take values below the top of the layer
        .take_while(|&(hgt, _)| hgt <= max_hgt)
        // Get the wind u-v components in m/s
        .map(|(hgt, wind)| {
            let WindUV { u, v } = WindUV::<MetersPSec>::from(wind);
            (hgt, u, v)
        })
        // Integration with the trapezoid rule, to find the mean value 
        .fold(
            (
                Meters(f64::MAX), // height from previous level, MAX is a sentinel
                MetersPSec(0.0),  // previous step u
                MetersPSec(0.0),  // previous step v 
                MetersPSec(0.0),  // integrated u component so far
                MetersPSec(0.0),  // integrated v component so far
                Meters(0.0),      // the total distance integrated so far
            ),
            |acc, (hgt, u, v)| {
                let (old_hgt, old_u, old_v, mut iu, mut iv, mut acc_dz) = acc;

                let dz = hgt - old_hgt;

                if dz < Meters(0.0) {  // Less than zero on the first step, initialization
                    return (hgt, u, v, iu, iv, acc_dz);
                }
                iu += (u + old_u) * dz.unpack();
                iv += (v + old_v) * dz.unpack();
                acc_dz += dz;

                (hgt, u, v, iu, iv, acc_dz)
            },
        );

    if old_hgt == Meters(f64::MAX) {
        // nothing was done, no points in layer
        return Err(AnalysisError::NotEnoughData);
    } else if dz.unpack().abs() < std::f64::EPSILON {
        // only one point, use that value. This seems imposible, but it is possible to make a layer
        // where the top and bottom are equal.
        iu = old_u;
        iv = old_v;
    } else {
        // we integrated, so divide by height and constant of 2 for trapezoid rule
        iu /= 2.0 * dz.unpack();
        iv /= 2.0 * dz.unpack();
    }

    Ok(WindUV { u: iu, v: iv })
}

/// Storm relative helicity.
pub fn sr_helicity<W>(
    layer: &Layer,
    storm_motion_uv_ms: W,
    snd: &Sounding,
) -> Result<IntHelicityM2pS2>
where
    WindUV<MetersPSec>: From<W>,
{
    let height = snd.height_profile();
    let wind = snd.wind_profile();
    let storm_motion_uv_ms = WindUV::<MetersPSec>::from(storm_motion_uv_ms);

    let bottom = layer.bottom.height.ok_or(AnalysisError::MissingValue)?;
    let top = layer.top.height.ok_or(AnalysisError::MissingValue)?;

    izip!(height, wind)
        // Filter out levels with missing values
        .filter_map(|(h, w)| {
            if let (Some(h), Some(w)) = (h.into_option(), w.into_option()) {
                Some((h, w))
            } else {
                None
            }
        })
        // Convert the wind and unpack it, subtract the storm motion.
        .map(|(h, w)| {
            let WindUV { u, v }: WindUV<MetersPSec> = From::<WindSpdDir<Knots>>::from(w);
            (h, (u - storm_motion_uv_ms.u), (v - storm_motion_uv_ms.v))
        })
        // Make windows so I can see three at a time
        .tuple_windows::<(_, _, _)>()
        // Skip levels where the middle of the window is below the bottom of the layer.
        .skip_while(|(_, (h, _, _), _)| *h < bottom)
        // Take until the middle of the window is at the top of the layer.
        .take_while(|(_, (h, _, _), _)| *h <= top)
        // Add in the derivative information
        .map(|((h0, u0, v0), (h1, u1, v1), (h2, u2, v2))| {
            let dz = Meters::from(h2 - h0).unpack();
            let du = MetersPSec::from(u2 - u0).unpack() / dz;
            let dv = MetersPSec::from(v2 - v0).unpack() / dz;
            (h1, u1, v1, du, dv)
        })
        // Make a window so I can see two at a time for integrating with the trapezoid rule
        .tuple_windows::<(_, _)>()
        // Integrate with the trapezoidal rule
        .fold(Err(AnalysisError::NotEnoughData), |acc, (lvl0, lvl1)| {
            let mut integrated_helicity: f64 = acc.unwrap_or(0.0);
            let (z0, u0, v0, du0, dv0) = lvl0;
            let (z1, u1, v1, du1, dv1) = lvl1;

            let h0 = u0 * dv0 - v0 * du0;
            let h1 = u1 * dv1 - v1 * du1;
            let h = MetersPSec::from(h0 + h1).unpack() * Meters::from(z1 - z0).unpack();
            integrated_helicity += h;

            Ok(integrated_helicity)
        })
        // Multiply by constant factors, -1 for the helicity formula, 2 for trapezoid rule,
        // and wrap in integrated helicity type
        .map(|integrated_helicity| IntHelicityM2pS2(-integrated_helicity / 2.0))
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

    // abort if not at least 250 meters of non-overlapping area.
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
