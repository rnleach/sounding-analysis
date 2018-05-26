//! Functions for doing parcel analysis on a sounding, specifically related to convection.
//!

use smallvec::SmallVec;

use metfor;
use sounding_base::{Profile::*, Sounding};

use error::*;
use layers::{pressure_layer, Layers};
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
/// Calculated by averaging potential temperature and mixing ratio in the lowest  100 hPa of the
/// sounding and calculating a temperature and dew point at the lowest pressure level  using those
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
        .filter_map(|p| *p)
        // Get the first of the non-None values
        .nth(0)
        // Check to make sure we got one, return error if not.
        .ok_or(AnalysisError::NoDataProfile)?;

    let (sum_t, sum_mw, count) = izip!(press, t, dp)
        // remove levels with missing values
        .filter_map(|(p, t, dp)| p.and_then(|p| t.and_then(|t| dp.and_then(|dp| Some((p, t, dp))))))
        // only go up 100 hPa above the lowest level
        .take_while(|&(p, _, _)| p >= bottom_p - 100.0)
        // convert to theta and mw, drop the pressure after this level
        .filter_map(|(p,t,dp)|{
            metfor::theta_kelvin(p, t).ok().and_then(|theta|{
                metfor::mixing_ratio(dp, p).ok().and_then(|mw| Some((theta, mw)))
            })
        })
        // calculate the sums and count needed for the average
        .fold((0.0f64, 0.0f64, 0.0f64), |acc, (theta, mw)| {
            let (sum_t, sum_mw, count) = acc;
            (sum_t + theta, sum_mw + mw, count + 1.0)
        });

    if count == 0.0 {
        return Err(AnalysisError::NotEnoughData);
    }

    // convert back to temperature and dew point at the lowest pressure level.
    let pressure = bottom_p;
    let temperature = metfor::temperature_c_from_theta(sum_t / count, bottom_p)
        .map_err(|_| AnalysisError::InvalidInput)?;
    let dew_point = metfor::dew_point_from_p_and_mw(bottom_p, sum_mw / count)
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

/// Get the parcel at a specific pressure.
pub fn pressure_parcel(snd: &Sounding, pressure: f64) -> Result<Parcel> {
    let row = ::interpolation::linear_interpolate_sounding(snd, pressure)?;

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

/// Get the most unstable parcel.
///
/// This is defined as the parcel in the lowest 300 hPa of the sounding with the highest equivalent
/// potential temperature.
// TODO: Change to highest cape? Much more expensive to compute cape for every parcel.
//       I could handle this by looking at the highest theta_e parcel in every 100 hPa layer and
//       just calculate cape for those. That is not guaranteed to converge on the absolute most
//       unstable parcel. But in some shallow convection situations, the highest theta_e parcel may
//       me WAY above the layer with convection. This default case will capture the most unstable
//       parcel in a deep, moist convection scenario.
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

/// Hold profiles for a parcel and it's environment.
#[derive(Debug, Clone)]
pub struct ParcelProfile {
    /// Pressure profile
    pub pressure: Vec<f64>,
    /// Height profile
    pub height: Vec<f64>,
    /// Parcel temperature profile
    pub parcel_t: Vec<f64>,
    /// Environment temperature profile
    pub environment_t: Vec<f64>,
    /// The original parcel
    pub parcel: Parcel,
}

impl ParcelProfile {}

/// Lift a parcel
pub fn lift_parcel(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let (lcl_pressure, lcl_temperature) = metfor::pressure_and_temperature_at_lcl(
        parcel.temperature,
        parcel.dew_point,
        parcel.pressure,
    ).map_err(|_| AnalysisError::InvalidInput)?;

    let lcl_temperature =
        metfor::kelvin_to_celsius(lcl_temperature).map_err(|_| AnalysisError::InvalidInput)?;

    let theta = parcel.theta()?;
    let theta_e = parcel.theta_e()?;

    let parcel_start_data = ::interpolation::linear_interpolate_sounding(snd, parcel.pressure)?;

    let lcl_env = ::interpolation::linear_interpolate_sounding(snd, lcl_pressure)?;
    let lcl_height = lcl_env.height.ok_or(AnalysisError::InvalidInput)?;
    let lcl_env_temperature = lcl_env.temperature.ok_or(AnalysisError::InvalidInput)?;

    let snd_pressure = snd.get_profile(Pressure);
    let hgt = snd.get_profile(GeopotentialHeight);
    let env_t = snd.get_profile(Temperature);

    let mut pressure = Vec::with_capacity(snd_pressure.len() + 5);
    let mut height = Vec::with_capacity(snd_pressure.len() + 5);
    let mut parcel_t = Vec::with_capacity(snd_pressure.len() + 5);
    let mut environment_t = Vec::with_capacity(snd_pressure.len() + 5);

    // Nested scope to limit closure borrows
    {
        // Helper function to calculate parcel temperature
        let calc_parcel_t = |tgt_pres| {
            if tgt_pres > lcl_pressure {
                // Dry adiabatic lifting
                metfor::temperature_c_from_theta(theta, tgt_pres)
                    .map_err(|_| AnalysisError::InvalidInput)
            } else {
                // Moist adiabatic lifting
                metfor::temperature_c_from_theta_e_saturated_and_pressure(tgt_pres, theta_e)
                    .map_err(|_| AnalysisError::InvalidInput)
            }
        };

        // Helper function to add row to parcel profile
        let mut add_row = |pp, hh, pcl_tt, env_tt| {
            pressure.push(pp);
            height.push(hh);
            parcel_t.push(pcl_tt);
            environment_t.push(env_tt);
        };

        // Start by adding the parcel level
        let mut p0 = parcel.pressure;
        let h0 = parcel_start_data.height.ok_or(AnalysisError::InvalidInput)?;
        let mut pcl_t0 = parcel.temperature;
        let mut env_t0 = parcel_start_data
            .temperature
            .ok_or(AnalysisError::InvalidInput)?;

        add_row(p0, h0, pcl_t0, env_t0);

        for (p, h, env_t) in izip!(snd_pressure, hgt, env_t) {
            // Unpack options
            let (p, h, env_t) = if let (&Some(p), &Some(h), &Some(env_t)) = (p, h, env_t) {
                (p, h, env_t)
            } else {
                continue;
            };

            // Check if this level has already been added or if we are below the parcel.
            if p0 <= p {
                continue;
            }

            // Check to see if we are passing the lcl
            if p0 > lcl_pressure && p < lcl_pressure {
                add_row(
                    lcl_pressure,
                    lcl_height,
                    lcl_temperature,
                    lcl_env_temperature,
                );
            }

            // Calculate the new parcel temperature
            let pcl_t = if let Ok(new_t) = calc_parcel_t(p) {
                new_t
            } else {
                break;
            };

            // Check to see if the parcel and environment soundings have crossed
            if (pcl_t0 < env_t0 && pcl_t > env_t) || (pcl_t0 > env_t0 && pcl_t < env_t) {
                let tgt_pres =
                    ::interpolation::linear_interp(0.0, pcl_t - env_t, pcl_t0 - env_t0, p, p0);
                let tgt_row = ::interpolation::linear_interpolate_sounding(snd, tgt_pres)?;
                let h2 = tgt_row.height.ok_or(AnalysisError::InvalidInput)?;
                let env_t2 = tgt_row.temperature.ok_or(AnalysisError::InvalidInput)?;
                add_row(tgt_pres, h2, env_t2, env_t2);
            }

            // Add the new values to the array
            add_row(p, h, pcl_t, env_t);

            // Remember them for the next iteration
            p0 = p;
            pcl_t0 = pcl_t;
            env_t0 = env_t;
        }
    }

    Ok(ParcelProfile {
        pressure,
        height,
        parcel_t,
        environment_t,
        parcel,
    })
}

/// Get the lfcs and el levels for a parcel.
pub fn cape_layers(parcel_profile: &ParcelProfile, snd: &Sounding) -> Layers {
    let mut to_ret = SmallVec::new();

    let lcl_pressure = if let Ok(pres) = metfor::pressure_hpa_at_lcl(
        parcel_profile.parcel.temperature,
        parcel_profile.parcel.dew_point,
        parcel_profile.parcel.pressure,
    ) {
        pres
    } else {
        return to_ret;
    };

    let l0 = izip!(
        &parcel_profile.pressure,
        &parcel_profile.parcel_t,
        &parcel_profile.environment_t
    );
    let l1 = izip!(&parcel_profile.parcel_t, &parcel_profile.environment_t).skip(1);

    let mut bottom = ::std::f64::MIN;
    let mut top = ::std::f64::MAX;

    if !(parcel_profile.parcel_t.is_empty() || parcel_profile.environment_t.is_empty()
        || parcel_profile.pressure.is_empty())
        && (parcel_profile.parcel_t[0] >= parcel_profile.environment_t[0])
    {
        bottom = parcel_profile.pressure[0];
    }

    izip!(l0, l1)
        .skip_while(|&(l0, _)| *l0.0 > lcl_pressure)
        .for_each(|(l0, l1)| {
            let (&p0, &pt0, &et0) = l0;
            let (&pt1, &et1) = l1;

            if pt0 <= et0 && pt1 > et1 {
                bottom = p0;
            }

            if pt0 >= et0 && pt1 < et1 {
                top = p0;
            }

            if bottom > top {
                pressure_layer(snd, bottom, top).ok().and_then(|lyr| {
                    to_ret.push(lyr);
                    Some(())
                });
                bottom = ::std::f64::MIN;
                top = ::std::f64::MAX;
            }
        });

    to_ret
}

/// Get the pressure layers with CIN in this profile.
pub fn cin_layers(parcel_profile: &ParcelProfile, snd: &Sounding) -> Layers {
    let mut to_ret = SmallVec::new();

    let l0 = izip!(
        &parcel_profile.pressure,
        &parcel_profile.parcel_t,
        &parcel_profile.environment_t
    );
    let l1 = izip!(&parcel_profile.parcel_t, &parcel_profile.environment_t).skip(1);

    let mut top = ::std::f64::MAX;
    let mut bottom = {
        let (p, pt, et) = (
            parcel_profile.pressure[0],
            parcel_profile.parcel_t[0],
            parcel_profile.environment_t[0],
        );
        if pt < et {
            p
        } else {
            ::std::f64::MIN
        }
    };

    izip!(l0, l1).for_each(|(l0, l1)| {
        let (&p0, &pt0, &et0) = l0;
        let (&pt1, &et1) = l1;

        if pt0 <= et0 && pt1 > et1 {
            top = p0;
        }

        if pt0 >= et0 && pt1 < et1 {
            bottom = p0;
        }

        if bottom > top {
            pressure_layer(snd, bottom, top).ok().and_then(|lyr| {
                to_ret.push(lyr);
                Some(())
            });
            bottom = ::std::f64::MIN;
            top = ::std::f64::MAX;
        }
    });

    to_ret
}

// TODO: descend parcel dry adiabatically
// TODO: descend parcel moist adiabatically

// TODO: cape, cin, el, lfc, lcl, dcape, ncape, hail zone cape
