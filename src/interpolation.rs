use metfor::{HectoPascal, Knots, Quantity, WindSpdDir, WindUV};
use optional::{none, Optioned};
use sounding_base::{DataRow, Sounding};
use std::ops::Sub;

use crate::error::AnalysisError::*;
use crate::error::*;

/// Interpolate values from the vertical sounding using pressure as the primary coordinate.
///
/// Returns a `DataRow` struct with interpolated values.
pub fn linear_interpolate_sounding(snd: &Sounding, tgt_p: HectoPascal) -> Result<DataRow> {
    let pressure: &[Optioned<HectoPascal>] = snd.pressure_profile();
    let temperature = snd.temperature_profile();
    let wet_bulb = snd.wet_bulb_profile();
    let dew_point = snd.dew_point_profile();
    let theta_e = snd.theta_e_profile();
    let wind = snd.wind_profile();
    let omega = snd.pvv_profile();
    let height = snd.height_profile();
    let cloud_fraction = snd.cloud_fraction_profile();

    let mut result = DataRow::default();
    result.pressure = Optioned::from(tgt_p);

    let mut below_idx: usize = 0;
    let mut above_idx: usize = 0;
    let mut found_bottom: bool = false;
    for (i, p) in pressure.iter().enumerate() {
        if let Some(p) = p.into_option() {
            if p > tgt_p {
                below_idx = i;
                found_bottom = true;
            } else if p < tgt_p && found_bottom {
                above_idx = i;
                break;
            } else if (p - tgt_p).unpack().abs() <= ::std::f64::EPSILON {
                return snd.data_row(i).ok_or(AnalysisError::InvalidInput);
            } else {
                break; // leave above_idx = 0 to signal error
            }
        }
    }

    if above_idx != 0 {
        let p_below = pressure[below_idx].unwrap();
        let p_above = pressure[above_idx].unwrap();
        let run = p_above - p_below;
        let dp = tgt_p - p_below;

        result.temperature = eval_linear_interp(below_idx, above_idx, run, dp, temperature);
        result.wet_bulb = eval_linear_interp(below_idx, above_idx, run, dp, wet_bulb);
        result.dew_point = eval_linear_interp(below_idx, above_idx, run, dp, dew_point);
        result.theta_e = eval_linear_interp(below_idx, above_idx, run, dp, theta_e);

        // Special interpolation for vectors
        if wind.len() > above_idx {
            if let (Some(w_below), Some(w_above)) =
                (wind[below_idx].into_option(), wind[above_idx].into_option())
            {
                let WindUV::<Knots> {
                    u: x_below,
                    v: y_below,
                } = WindUV::from(w_below);
                let WindUV::<Knots> {
                    u: x_above,
                    v: y_above,
                } = WindUV::from(w_above);
                let dp = dp.unpack();
                let run = run.unpack();

                let rise_x = x_above - x_below;
                let rise_y = y_above - y_below;

                let x = x_below + rise_x * (dp / run);
                let y = y_below + rise_y * (dp / run);

                let interped_wind = WindSpdDir::from(WindUV { u: x, v: y });

                result.wind = interped_wind.into();
            }
        }

        result.pvv = eval_linear_interp(below_idx, above_idx, run, dp, omega);
        result.height = eval_linear_interp(below_idx, above_idx, run, dp, height);
        result.cloud_fraction = eval_linear_interp(below_idx, above_idx, run, dp, cloud_fraction);
        Ok(result)
    } else {
        // Target pressure was above or below actual pressures in the sounding.
        Err(InvalidInput)
    }
}

/// Interpolate values given two parallel vectors of data and a target value.
// FIXME: Currently assume xs are sorted in descending order, change to just assuming monotonic
#[inline]
pub fn linear_interpolate<X, Y>(xs: &[Optioned<X>], ys: &[Optioned<Y>], target_x: X) -> Optioned<Y>
where
    X: Quantity + optional::Noned + PartialOrd + Sub<X>,
    <X as Sub<X>>::Output: Quantity + optional::Noned,
    Y: Quantity + optional::Noned,
{
    debug_assert_eq!(xs.len(), ys.len());

    let mut below_idx: usize = 0;
    let mut above_idx: usize = 0;
    let mut found_bottom: bool = false;
    for (i, x) in xs.iter().enumerate() {
        if x.is_some() {
            let x = x.unpack();
            if x > target_x {
                below_idx = i;
                found_bottom = true;
            } else if x < target_x && found_bottom {
                above_idx = i;
                break;
            } else if (x - target_x).unpack().abs() <= ::std::f64::EPSILON {
                return ys[i];
            } else {
                break; // leave above_idx = 0 to signal error
            }
        }
    }

    if above_idx != 0 {
        let x_below = xs[below_idx].unwrap();
        let x_above = xs[above_idx].unwrap();
        let run = x_above - x_below;
        let dx = target_x - x_below;

        eval_linear_interp(below_idx, above_idx, run, dx, ys)
    } else {
        none()
    }
}

#[inline]
fn eval_linear_interp<QX, QY>(
    blw_idx: usize,
    abv_idx: usize,
    run: QX,
    dp: QX,
    array: &[Optioned<QY>],
) -> Optioned<QY>
where
    QX: Quantity + optional::Noned,
    QY: Quantity + optional::Noned,
{
    if array.len() > abv_idx {
        if array[blw_idx].is_some() && array[abv_idx].is_some() {
            let (val_below, val_above) = (
                array[blw_idx].unpack().unpack(),
                array[abv_idx].unpack().unpack(),
            );
            let rise: f64 = (val_above - val_below).unpack();
            let run: f64 = run.unpack();
            let dp: f64 = dp.unpack();
            Optioned::from(QY::pack(val_below + dp * rise / run))
        } else {
            Optioned::default()
        }
    } else {
        Optioned::default()
    }
}

#[inline]
pub(crate) fn linear_interp<X, Y>(x_val: X, x1: X, x2: X, y1: Y, y2: Y) -> Y
where
    X: Sub<X> + Copy,
    <X as Sub<X>>::Output: Quantity,
    Y: Quantity + Sub<Y>,
    <Y as Sub<Y>>::Output: Quantity,
{
    let run = (x2 - x1).unpack();
    let rise = (y2 - y1).unpack();
    let dx = (x_val - x1).unpack();

    Y::pack(y1.unpack() + dx * (rise / run))
}
