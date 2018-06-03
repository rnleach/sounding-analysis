//! Functions for doing parcel analysis on a sounding, specifically related to convection.
//!
use optional::Optioned;
use smallvec::SmallVec;

use metfor;
use sounding_base::{DataRow, Profile::*, Sounding};

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
        .filter_map(|&p| p.clone().take())
        // Get the first of the non-None values
        .nth(0)
        // Check to make sure we got one, return error if not.
        .ok_or(AnalysisError::NoDataProfile)?;

    let (sum_t, sum_mw, count) = izip!(press, t, dp)
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

    let pressure = snd.get_surface_value(StationPressure)
        .ok_or(AnalysisError::MissingValue)?;
    let temperature = snd.get_surface_value(Temperature)
        .ok_or(AnalysisError::MissingValue)?;
    let dew_point = snd.get_surface_value(DewPoint)
        .ok_or(AnalysisError::MissingValue)?;

    Ok(Parcel {
        temperature,
        pressure,
        dew_point,
    })
}

/// Get the parcel at a specific pressure.
pub fn pressure_parcel(snd: &Sounding, pressure: f64) -> Result<Parcel> {
    let row = ::interpolation::linear_interpolate_sounding(snd, pressure)?;

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

    let (idx, _) = izip!(0.., press, theta_e)
        .take_while(|&(_, press_opt, _)| {
            if press_opt.is_some() {
                if press_opt.unpack() >= 300.0 {
                    true
                } else {
                    false // only stop if we are sure we are above 300.0 hPa
                }
            } else {
                true // Don't stop just because one is missing or the none value
            }
        })
        .filter_map(|(idx, _, theta_e_opt)| {
            if theta_e_opt.is_some() {
                Some((idx, theta_e_opt.unpack()))
            } else {
                None
            }
        })
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

    Parcel::from_datarow(row).ok_or(AnalysisError::MissingValue)
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
            let (p, h, env_t) = if p.is_some() && h.is_some() && env_t.is_some() {
                (p.unpack(), h.unpack(), env_t.unpack())
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

/// Descend a parcel dry adiabatically.
pub fn descend_dry_adiabatically(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let theta = parcel.theta()?;
    descend_parcel(parcel, snd, theta, ::metfor::temperature_c_from_theta)
}

/// Descend a parcel dry adiabatically
pub fn descend_moist_adiabatically(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let theta = parcel.theta_e()?;

    let theta_func = |theta_e, press| {
        ::metfor::temperature_c_from_theta_e_saturated_and_pressure(press, theta_e)
    };

    descend_parcel(parcel, snd, theta, theta_func)
}

fn descend_parcel<F>(
    parcel: Parcel,
    snd: &Sounding,
    theta: f64,
    theta_func: F,
) -> Result<ParcelProfile>
where
    F: Fn(f64, f64) -> ::std::result::Result<f64, ::metfor::MetForErr>,
{
    let mut pressure = Vec::new();
    let mut height = Vec::new();
    let mut parcel_t = Vec::new();
    let mut environment_t = Vec::new();

    // Actually start at the bottom and work up.
    let press = snd.get_profile(Pressure);
    let env_t = snd.get_profile(Temperature);
    let hght = snd.get_profile(GeopotentialHeight);

    // Nested scope to limit borrows
    {
        // Helper function to add row to parcel profile
        let mut add_row = |pp, hh, pcl_tt, env_tt| {
            pressure.push(pp);
            height.push(hh);
            parcel_t.push(pcl_tt);
            environment_t.push(env_tt);
        };

        izip!(press, hght, env_t)
            .take_while(|(p_opt, _, _)| {
                if p_opt.is_some() {
                    p_opt.unpack() >= parcel.pressure
                } else {
                    true
                }
            })
            .filter_map(|(p_opt, h_opt, e_t_opt)| {
                if p_opt.is_some() && h_opt.is_some() && e_t_opt.is_some() {
                    Some((p_opt.unpack(), h_opt.unpack(), e_t_opt.unpack()))
                } else {
                    None
                }
            })
            .for_each(|(p, h, et)| {
                theta_func(theta, p).ok().and_then(|pt| {
                    add_row(p, h, pt, et);
                    Some(())
                });
            });

        // Add the parcel layer also
        let parcel_level = ::interpolation::linear_interpolate_sounding(snd, parcel.pressure)?;
        let parcel_height = parcel_level.height.ok_or(AnalysisError::MissingValue)?;
        add_row(
            parcel.pressure,
            parcel_height,
            parcel.temperature,
            parcel.temperature,
        );
    }

    Ok(ParcelProfile {
        pressure,
        height,
        parcel_t,
        environment_t,
        parcel,
    })
}

// TODO: cape, cin, el, lfc, lcl, dcape, ncape, hail zone cape
