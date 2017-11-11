#![warn(missing_docs)]
//! Functions for analyzing soundings from the
//! [sounding-base](https://github.com/rnleach/sounding-base.git) crate.

extern crate sounding_base;
use sounding_base::{Sounding, DataRow, OptionVal};

/// Interpolate values from the vertical sounding using pressure as the primary coordinate.
///
/// Returns a `DataRow` struct with interpolated values.
#[cfg_attr(feature = "cargo-clippy", allow(cyclomatic_complexity))]
pub fn linear_interpolate(snd: &Sounding, target_p: f64) -> DataRow {

    macro_rules! linear_interp {
        ($res:ident, $blw_idx:ident, $abv_idx:ident,  $run:ident, $dp:ident, $array:ident) => {
            if $array.len() > $abv_idx {
                if let (Some(val_below), Some(val_above)) =
                    (
                        $array[$blw_idx].as_option(),
                        $array[$abv_idx].as_option(),
                    )
                {
                    let rise = val_above - val_below;
                    $res.$array = (val_below + $dp * rise/$run).into();
                }
            }
        };
    }

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
    result.pressure = OptionVal::from(target_p);

    let mut below_idx: usize = 0;
    let mut above_idx: usize = 0;
    for (i, p) in pressure.iter().enumerate() {
        if let Some(p) = p.as_option() {
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

        linear_interp!(result, below_idx, above_idx, run, dp, temperature);
        linear_interp!(result, below_idx, above_idx, run, dp, wet_bulb);
        linear_interp!(result, below_idx, above_idx, run, dp, dew_point);
        linear_interp!(result, below_idx, above_idx, run, dp, theta_e);

        // Special interpolation for anlges
        if direction.len() > above_idx {
            if let (Some(dir_below), Some(dir_above)) =
                (
                    direction[below_idx].as_option(),
                    direction[above_idx].as_option(),
                )
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

        linear_interp!(result, below_idx, above_idx, run, dp, speed);
        linear_interp!(result, below_idx, above_idx, run, dp, omega);
        linear_interp!(result, below_idx, above_idx, run, dp, height);
        linear_interp!(result, below_idx, above_idx, run, dp, cloud_fraction);
    }

    result
}
