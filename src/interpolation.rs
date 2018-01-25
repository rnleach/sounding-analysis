use sounding_base::{DataRow, Sounding};

/// Interpolate values from the vertical sounding using pressure as the primary coordinate.
///
/// Returns a `DataRow` struct with interpolated values.
pub fn linear_interpolate(snd: &Sounding, target_p: f64) -> DataRow {
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
    result.pressure = Option::from(target_p);

    let mut below_idx: usize = 0;
    let mut above_idx: usize = 0;
    for (i, p) in pressure.iter().enumerate() {
        if let Some(p) = *p {
            if p > target_p {
                below_idx = i;
            }
            if p < target_p {
                above_idx = i;
                break;
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

        // Special interpolation for anlges
        if direction.len() > above_idx {
            if let (Some(dir_below), Some(dir_above)) = (direction[below_idx], direction[above_idx])
            {
                let x_below = dir_below.to_radians().sin();
                let x_above = dir_above.to_radians().sin();
                let y_below = dir_below.to_radians().cos();
                let y_above = dir_above.to_radians().cos();

                let rise_x = x_above - x_below;
                let rise_y = y_above - y_below;

                let x_dir = x_below + dp * rise_x / run;
                let y_dir = y_below + dp * rise_y / run;

                let mut dir = x_dir.atan2(y_dir).to_degrees();

                while dir < 0.0 {
                    dir += 360.0;
                }
                while dir > 360.0 {
                    dir -= 360.0;
                }

                result.direction = dir.into();
            }
        }

        result.speed = eval_linear_interp(below_idx, above_idx, run, dp, speed);
        result.omega = eval_linear_interp(below_idx, above_idx, run, dp, omega);
        result.height = eval_linear_interp(below_idx, above_idx, run, dp, height);
        result.cloud_fraction = eval_linear_interp(below_idx, above_idx, run, dp, cloud_fraction);
    }

    result
}

fn eval_linear_interp(
    blw_idx: usize,
    abv_idx: usize,
    run: f64,
    dp: f64,
    array: &[Option<f64>],
) -> Option<f64> {
    if array.len() > abv_idx {
        if let (Some(val_below), Some(val_above)) = (array[blw_idx], array[abv_idx]) {
            let rise = val_above - val_below;
            Option::from(val_below + dp * rise / run)
        } else {
            Option::default()
        }
    } else {
        Option::default()
    }
}

pub fn linear_interp(x_val: f64, x1: f64, x2: f64, y1: f64, y2: f64) -> f64 {
    let run = x2 - x1;
    let rise = y2 - y1;
    let dx = x_val - x1;

    y1 + rise / run * dx
}
