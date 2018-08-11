use optional::{none, Optioned};

use sounding_base::{DataRow, Sounding};

use error::AnalysisError::*;
use error::*;

// FIXME: Use interpolation error.

/// Interpolate values from the vertical sounding using pressure as the primary coordinate.
///
/// Returns a `DataRow` struct with interpolated values.
pub fn linear_interpolate_sounding(snd: &Sounding, target_p: f64) -> Result<DataRow> {
    use sounding_base::Profile::*;

    let pressure = snd.get_profile(Pressure);
    let temperature = snd.get_profile(Temperature);
    let wet_bulb = snd.get_profile(WetBulb);
    let dew_point = snd.get_profile(DewPoint);
    let theta_e = snd.get_profile(ThetaE);
    let direction = snd.get_profile(WindDirection);
    let speed = snd.get_profile(WindSpeed);
    let omega = snd.get_profile(PressureVerticalVelocity);
    let height = snd.get_profile(GeopotentialHeight);
    let cloud_fraction = snd.get_profile(CloudFraction);

    let mut result = DataRow::default();
    result.pressure = Optioned::from(target_p);

    let mut below_idx: usize = 0;
    let mut above_idx: usize = 0;
    let mut found_bottom: bool = false;
    for (i, p) in pressure.iter().enumerate() {
        if p.is_some() {
            let p = p.unpack();
            if p > target_p {
                below_idx = i;
                found_bottom = true;
            } else if p < target_p && found_bottom {
                above_idx = i;
                break;
            } else if (p - target_p).abs() <= ::std::f64::EPSILON {
                return snd.get_data_row(i).ok_or(AnalysisError::InvalidInput);
            } else {
                break; // leave above_idx = 0 to signal error
            }
        }
    }

    if above_idx != 0 {
        let p_below = pressure[below_idx].unwrap();
        let p_above = pressure[above_idx].unwrap();
        let run = p_above - p_below;
        let dp = target_p - p_below;

        result.temperature = eval_linear_interp(below_idx, above_idx, run, dp, temperature);
        result.wet_bulb = eval_linear_interp(below_idx, above_idx, run, dp, wet_bulb);
        result.dew_point = eval_linear_interp(below_idx, above_idx, run, dp, dew_point);
        result.theta_e = eval_linear_interp(below_idx, above_idx, run, dp, theta_e);

        // Special interpolation for vectors
        if direction.len() > above_idx && speed.len() > above_idx {
            if direction[below_idx].is_some()
                && speed[below_idx].is_some()
                && direction[above_idx].is_some()
                && speed[above_idx].is_some()
            {
                let (dir_below, dir_above) =
                    (direction[below_idx].unpack(), direction[above_idx].unpack());
                let (spd_below, spd_above) = (speed[below_idx].unpack(), speed[above_idx].unpack());
                let x_below = dir_below.to_radians().sin() * spd_below;
                let x_above = dir_above.to_radians().sin() * spd_above;
                let y_below = dir_below.to_radians().cos() * spd_below;
                let y_above = dir_above.to_radians().cos() * spd_above;

                let rise_x = x_above - x_below;
                let rise_y = y_above - y_below;

                let x = x_below + dp * rise_x / run;
                let y = y_below + dp * rise_y / run;

                let mut dir = x.atan2(y).to_degrees();

                while dir < 0.0 {
                    dir += 360.0;
                }
                while dir > 360.0 {
                    dir -= 360.0;
                }

                let spd = x.hypot(y);

                result.direction = dir.into();
                result.speed = spd.into();
            }
        }

        result.omega = eval_linear_interp(below_idx, above_idx, run, dp, omega);
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
pub fn linear_interpolate(
    xs: &[Optioned<f64>],
    ys: &[Optioned<f64>],
    target_x: f64,
) -> Optioned<f64> {
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
            } else if (x - target_x).abs() <= ::std::f64::EPSILON {
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

fn eval_linear_interp(
    blw_idx: usize,
    abv_idx: usize,
    run: f64,
    dp: f64,
    array: &[Optioned<f64>],
) -> Optioned<f64> {
    if array.len() > abv_idx {
        if array[blw_idx].is_some() && array[abv_idx].is_some() {
            let (val_below, val_above) = (array[blw_idx].unpack(), array[abv_idx].unpack());
            let rise = val_above - val_below;
            Optioned::from(val_below + dp * rise / run)
        } else {
            Optioned::default()
        }
    } else {
        Optioned::default()
    }
}

pub(crate) fn linear_interp(x_val: f64, x1: f64, x2: f64, y1: f64, y2: f64) -> f64 {
    let run = x2 - x1;
    let rise = y2 - y1;
    let dx = x_val - x1;

    y1 + rise / run * dx
}
