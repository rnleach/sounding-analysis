//! Functions for doing parcel analysis on a sounding, specifically related to convection.
//!
use optional::Optioned;

use metfor;
use sounding_base::{DataRow, Profile::*, Sounding};

use error::*;
use interpolation::{linear_interp, linear_interpolate_sounding};
use profile::equivalent_potential_temperature;

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
        Ok(metfor::theta_kelvin(self.pressure, self.temperature)?)
    }

    /// Get the equivalent potential temperature of the parcel
    pub fn theta_e(&self) -> Result<f64> {
        Ok(metfor::theta_e_kelvin(
            self.temperature,
            self.dew_point,
            self.pressure,
        )?)
    }

    /// Get the mixing ratio of the parcel.AnalysisError
    pub fn mixing_ratio(&self) -> Result<f64> {
        Ok(metfor::mixing_ratio(self.dew_point, self.pressure)?)
    }

    /// Get the virtual temperature of the parcel
    pub fn virtual_temperature_c(&self) -> Result<f64> {
        Ok(metfor::virtual_temperature_c(
            self.temperature,
            self.dew_point,
            self.pressure,
        )?)
    }

    /// Try to convert a `DataRow` to a `Parcel`.
    pub fn from_datarow(mut dr: DataRow) -> Option<Self> {
        let temperature = dr.temperature.take()?;
        let pressure = dr.pressure.take()?;
        let dew_point = dr.dew_point.take()?;

        Some(Parcel {
            temperature,
            pressure,
            dew_point,
        })
    }
}

/// Create a mixed layer parcel.
///
/// Calculated by averaging potential temperature and mixing ratio in the lowest 100 hPa of the
/// sounding and calculating a temperature and dew point at the lowest pressure level using those
/// averaged values.
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
        // Remove None values
        .filter_map(|&p| p.clone().take())
        // Get the first of the non-None values
        .nth(0)
        // Check to make sure we got one, return error if not.
        .ok_or(AnalysisError::NoDataProfile)?;

    let (sum_p, sum_t, sum_mw) = izip!(press, t, dp)
        // remove levels with missing values
        .filter_map(|(p, t, dp)| {
            if p.is_some() && t.is_some() && dp.is_some() {
                Some((p.unpack(), t.unpack(), dp.unpack()))
            } else {
                None
            }
        })
        // only go up 100 hPa above the lowest level
        .take_while(|&(p, _, _)| p >= bottom_p - 100.0)
        // convert to theta and mw
        .filter_map(|(p,t,dp)|{
            metfor::theta_kelvin(p, t).ok().and_then(|theta|{
                metfor::mixing_ratio(dp, p).ok().and_then(|mw| Some((p, theta, mw)))
            })
        })
        // calculate the sums and count needed for the average
        .fold((0.0f64, 0.0f64, 0.0f64), |acc, (p, theta, mw)| {
            let (sum_p, sum_t, sum_mw) = acc;
            (sum_p + p, sum_t + theta * p, sum_mw + mw * p)
        });

    if sum_p == 0.0 {
        return Err(AnalysisError::NotEnoughData);
    }

    // convert back to temperature and dew point at the lowest pressure level.
    let pressure = bottom_p;
    let temperature = metfor::temperature_c_from_theta(sum_t / sum_p, bottom_p)?;
    let dew_point = metfor::dew_point_from_p_and_mw(bottom_p, sum_mw / sum_p)
        .map_err(|_| AnalysisError::InvalidInput)?;

    Ok(Parcel {
        temperature,
        pressure,
        dew_point,
    })
}

/// Get a surface parcel.
pub fn surface_parcel(snd: &Sounding) -> Result<Parcel> {
    use sounding_base::Surface::*;

    let pressure = snd
        .get_surface_value(StationPressure)
        .ok_or(AnalysisError::MissingValue)?;
    let temperature = snd
        .get_surface_value(Temperature)
        .ok_or(AnalysisError::MissingValue)?;
    let dew_point = snd
        .get_surface_value(DewPoint)
        .ok_or(AnalysisError::MissingValue)?;

    Ok(Parcel {
        temperature,
        pressure,
        dew_point,
    })
}

/// Get the parcel at a specific pressure.
pub fn pressure_parcel(snd: &Sounding, pressure: f64) -> Result<Parcel> {
    let row = linear_interpolate_sounding(snd, pressure)?;

    Parcel::from_datarow(row).ok_or(AnalysisError::MissingValue)
}

/// Get the most unstable parcel.
///
/// This is defined as the parcel in the lowest 300 hPa of the sounding with the highest equivalent
/// potential temperature.
pub fn most_unstable_parcel(snd: &Sounding) -> Result<Parcel> {
    use sounding_base::Profile::*;

    let temp_vec: Vec<Optioned<f64>>;

    let theta_e = if !snd.get_profile(ThetaE).is_empty() {
        snd.get_profile(ThetaE)
    } else {
        temp_vec = equivalent_potential_temperature(snd);
        &temp_vec
    };
    let press = snd.get_profile(Pressure);

    let top_pressure = press
        .iter()
        .filter(|p| p.is_some())
        .nth(0)
        .ok_or(AnalysisError::NotEnoughData)?
        .unpack()
        - 300.0;

    let (idx, _) = izip!(0.., press, theta_e)
        .take_while(|&(_, press_opt, _)| {
            if press_opt.is_some() {
                if press_opt.unpack() >= top_pressure {
                    true
                } else {
                    false // only stop if we are sure we are above 300.0 hPa
                }
            } else {
                true // Don't stop just because one is missing or the none value
            }
        }).filter_map(|(idx, _, theta_e_opt)| {
            if theta_e_opt.is_some() {
                Some((idx, theta_e_opt.unpack()))
            } else {
                None
            }
        }).fold(
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

    Parcel::from_datarow(row).ok_or(AnalysisError::MissingValue)
}

/// Get the convective parcel - this is the surface parcel that will rise without any lifting.
///
/// This is based off the mixing ratio of the mixed layer parcel.
pub fn convective_parcel(snd: &Sounding) -> Result<Parcel> {
    let Parcel {
        dew_point: dp,
        pressure: initial_p,
        ..
    } = mixed_layer_parcel(snd)?;

    let tgt_mw = metfor::mixing_ratio(dp, initial_p)?;

    let pressure = snd.get_profile(Pressure);
    let temperature = snd.get_profile(Temperature);

    let (tgt_p, tgt_t) = izip!(pressure, temperature)
        .filter_map(|(p, t)| {
            let p = p.into_option()?;
            let t = t.into_option()?;
            let mw = metfor::mixing_ratio(t, p).ok()?;

            Some((p, t, mw))
        }).skip_while(|(_, _, mw)| *mw <= tgt_mw)
        .scan((0.0, 0.0, 0.0), |(old_p, old_t, old_mw), (p, t, mw)| {
            let result = if mw < tgt_mw && *old_mw >= tgt_mw {
                // found the crossing
                let tgt_p = linear_interp(tgt_mw, *old_mw, mw, *old_p, p);
                let tgt_t = linear_interp(tgt_mw, *old_mw, mw, *old_t, t);

                Some(Some((tgt_p, tgt_t)))
            } else {
                Some(None)
            };

            // save these as the new ones
            *old_p = p;
            *old_t = t;
            *old_mw = mw;

            result
        }).filter_map(|opt_opt| opt_opt)
        .nth(0)
        .ok_or(AnalysisError::NotEnoughData)?;

    let tgt_theta = metfor::theta_kelvin(tgt_p, tgt_t)?;
    let tgt_t = metfor::temperature_c_from_theta(tgt_theta, initial_p)?;

    Ok(Parcel {
        pressure: initial_p,
        temperature: tgt_t,
        dew_point: dp,
    })
}
