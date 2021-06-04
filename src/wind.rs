use crate::{
    error::{AnalysisError, Result},
    layers::{self, Layer},
    levels::height_level,
    sounding::Sounding,
};
use itertools::{izip, Itertools};
use metfor::{IntHelicityM2pS2, Knots, Meters, MetersPSec, Quantity, WindSpdDir, WindUV};
use std::iter::once;

/// Calculate the mean wind in a layer.
///
/// This is NOT the pressure weighted mean.
///
pub fn mean_wind(layer: &Layer, snd: &Sounding) -> Result<WindUV<MetersPSec>> {
    let height = snd.height_profile();
    let wind = snd.wind_profile();

    let max_hgt = layer.top.height.ok_or(AnalysisError::MissingValue)?;
    let min_hgt = layer.bottom.height.ok_or(AnalysisError::MissingValue)?;

    let bottom_wind = layer.bottom.wind;
    let top_wind = layer.top.wind;

    let intermediate_layers = izip!(height, wind)
        .filter_map(|(hgt, wind)| hgt.into_option().map(|h| (h, wind)))
        // Skip values below the layer
        .skip_while(|&(hgt, _)| hgt < min_hgt)
        // Only take values below the top of the layer
        .take_while(|&(hgt, _)| hgt < max_hgt);

    let (mut iu, mut iv, dz) =
        // Start at the bottom of the layer
        once((min_hgt, &bottom_wind))
        // Add in any intermediate layers
        .chain(intermediate_layers)
        // Finish with the top layer
        .chain(once((max_hgt, &top_wind)))
        // Filter out missing values
        .filter_map(|(hgt, wind)| wind.map(|w| (hgt, w)))
        // Get the wind u-v components in m/s
        .map(|(hgt, wind)| {
            let WindUV { u, v } = WindUV::<MetersPSec>::from(wind);
            (hgt, u, v)
        })
        // Make windows to see two points at a time for trapezoid rule integration
        .tuple_windows::<(_, _)>()
        // Integration with the trapezoid rule to find the mean value
        .fold(
            (
                MetersPSec(0.0), // integrated u component so far
                MetersPSec(0.0), // integrated v component so far
                Meters(0.0),     // the total distance integrated so far
            ),
            |acc, ((h0, u0, v0), (h1, u1, v1))| {
                let (mut iu, mut iv, mut acc_dz) = acc;

                let dz = h1 - h0;

                iu += (u0 + u1) * dz.unpack();
                iv += (v0 + v1) * dz.unpack();
                acc_dz += dz;

                (iu, iv, acc_dz)
            },
        );

    if dz == Meters(0.0) {
        // nothing was done, 1 or zero points in the layer
        return Err(AnalysisError::NotEnoughData);
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
        .filter(|(h, w)| h.is_some() && w.is_some())
        // Unwrap from the `Optioned` type
        .map(|(h, w)| (h.unpack(), w.unpack()))
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
            let dz = (h2 - h0).unpack();
            let du = (u2 - u0).unpack() / dz;
            let dv = (v2 - v0).unpack() / dz;
            (h1, u1, v1, du, dv)
        })
        // Calculate the helicity (not, this is not the integrated helicity yet)
        .map(|(z, u, v, du, dv)| (z, u * dv - v * du))
        // Make a window so I can see two at a time for integrating with the trapezoid rule
        .tuple_windows::<(_, _)>()
        // Integrate with the trapezoidal rule
        .fold(Err(AnalysisError::NotEnoughData), |acc, (lvl0, lvl1)| {
            let mut integrated_helicity: f64 = acc.unwrap_or(0.0);
            let (z0, h0) = lvl0;
            let (z1, h1) = lvl1;

            let h = (h0 + h1).unpack() * (z1 - z0).unpack();
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
