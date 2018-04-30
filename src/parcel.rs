//! Functions for doing parcel analysis on a sounding, specifically related to convection.
//!

use metfor;
use sounding_base::Sounding;

use error::*;
use profile::{equivalent_potential_temperature};

/// Variables defining a parcel as used in parcel analysis.
#[derive(Debug, Clone, Copy)]
pub struct Parcel {
    /// Temperature in C
    pub temperature: f64,
    /// Pressure in hPa
    pub pressure: f64,
    /// Dew point in C
    pub dew_point: f64,
}

impl Parcel {
    /// Get the potential temperatures of the parcel
    pub fn theta(&self) -> Result<f64> {
        metfor::theta_kelvin(self.pressure, self.temperature)
            .map_err(|_| AnalysisError::InvalidInput)
    }

    /// Get the equivalent potential temperature of the parcel
    pub fn theta_e(&self) -> Result<f64> {
        metfor::theta_e_kelvin(self.temperature, self.dew_point, self.pressure)
            .map_err(|_| AnalysisError::InvalidInput)
    }

    /// Get the mixing ratio of the parcel.AnalysisError
    pub fn mixing_ratio(&self) -> Result<f64> {
        metfor::mixing_ratio(self.dew_point, self.pressure).map_err(|_| AnalysisError::InvalidInput)
    }
}

/// Create a mixed layer parcel.
///
/// The values in this parcel are the simple mean of the lowest 100 hPa of the sounding.
pub fn mixed_layer_parcel(snd: &Sounding) -> Result<Parcel> {
    use sounding_base::Profile::*;

    let press = snd.get_profile(Pressure);
    let t = snd.get_profile(Temperature);
    let dp = snd.get_profile(DewPoint);

    if press.is_empty() || t.is_empty() || dp.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    let bottom_p = press
        .iter()
        .filter_map(|p| *p)
        .nth(0)
        .ok_or(AnalysisError::NoDataProfile)?;

    let (sum_p, sum_t, sum_dp, count) = izip!(press, t, dp)
        .filter_map(|(p, t, dp)| p.and_then(|p| t.and_then(|t| dp.and_then(|dp| Some((p, t, dp))))))
        .take_while(|&(p, _, _)| p >= bottom_p - 100.0)
        .fold((0.0f64, 0.0f64, 0.0f64, 0.0f64), |acc, (p, t, dp)| {
            let (sum_p, sum_t, sum_dp, count) = acc;
            (sum_p + p, sum_t + t, sum_dp + dp, count + 1.0)
        });

    if count == 0.0 {
        return Err(AnalysisError::NotEnoughData);
    }

    Ok(Parcel {
        temperature: (sum_t / count),
        pressure: (sum_p / count),
        dew_point: (sum_dp / count),
    })
}

/// Get a surface parcel.
pub fn surface_parcel(snd: &Sounding) -> Result<Parcel> {
    use sounding_base::Surface::*;

    snd.get_surface_value(StationPressure)
        .and_then(|pressure| {
            snd.get_surface_value(Temperature).and_then(|temperature| {
                snd.get_surface_value(DewPoint).and_then(|dew_point| {
                    Some(Parcel {
                        temperature,
                        pressure,
                        dew_point,
                    })
                })
            })
        })
        .ok_or(AnalysisError::MissingValue)
}

/// Get the most unstable parcel.
///
/// This is defined as the parcel in the lowest 300 hPa of the sounding with the highest equivalent
/// potential temperature.
pub fn most_unstable_parcel(snd: &Sounding) -> Result<Parcel> {
    use sounding_base::Profile::*;

    let temp_vec: Vec<Option<f64>>;

    let theta_e = if !snd.get_profile(ThetaE).is_empty() {
        snd.get_profile(ThetaE)
    } else {
        temp_vec = equivalent_potential_temperature(snd);
        &temp_vec
    };
    let press = snd.get_profile(Pressure);

    let (idx, _) = izip!(0.., press, theta_e)
        .take_while(|&(_, press_opt, _)| {
            press_opt
                .and_then(|p| if p < 300.0 { None } else { Some(()) })
                .is_some()
        })
        .filter_map(|(idx, _, theta_e_opt)| theta_e_opt.and_then(|theta_e| Some((idx, theta_e))))
        .fold(
            (::std::usize::MAX, ::std::f64::MIN),
            |(max_idx, max_val), (i, theta_e)| {
                if theta_e > max_val {
                    (i, theta_e)
                } else {
                    (max_idx, max_val)
                }
            },
        );

    let row = snd.get_data_row(idx).ok_or(AnalysisError::NoDataProfile)?;

    row.pressure
        .and_then(|pressure| {
            row.temperature.and_then(|temperature| {
                row.dew_point.and_then(|dew_point| {
                    Some(Parcel {
                        temperature,
                        pressure,
                        dew_point,
                    })
                })
            })
        })
        .ok_or(AnalysisError::MissingValue)
}

// TODO: cape, cin, el, lfc, lcl, dcape, ncape, hail zone cape
