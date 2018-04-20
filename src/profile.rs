//! Create profiles.
//!
//! There are two kinds of profiles:
//!   - Those created from a sounding, the output will be at the same levels as the sounding and
//!     these are suitable to be set as a profile in the sounding. For example, calculating a
//!     wet bulb or relative humidity profile from a sounding with temperature and dew point. If one
//!     of the profiles required for the analysis in the sounding is missing, the result cannot be
//!     calculated and an empty vector is returned.
//!   - Those created for parcel analysis on a sounding. These require an initial parcel and a
//!     sounding. They are useful doing parcel analysis for variables such CAPE.
//!

use std::iter::once;

use metfor;
use sounding_base::Sounding;
use sounding_base::Profile::*;

use error::*;
use parcel::Parcel;

/// Given a sounding, calculate a profile of wet bulb temperature.
pub fn wet_bulb(snd: &Sounding) -> Vec<Option<f64>> {
    let p_profile = snd.get_profile(Pressure);
    let t_profile = snd.get_profile(Temperature);
    let dp_profile = snd.get_profile(DewPoint);

    if p_profile.len().min(t_profile.len()).min(dp_profile.len()) == 0 {
        return vec![];
    }

    izip!(p_profile, t_profile, dp_profile)
        .map(|(p_opt, t_opt, dp_opt)| {
            p_opt.and_then(|p| {
                t_opt.and_then(|t| {
                    dp_opt.and_then(|dp| {
                        metfor::wet_bulb_c(t, dp, p)
                        // Ignore errors, if not possible to calculate just use missing value.
                        .ok()
                    })
                })
            })
        })
        .collect()
}

/// Given a sounding, calculate a profile of relative humidity.
pub fn relative_humidity(snd: &Sounding) -> Vec<Option<f64>> {
    let t_profile = snd.get_profile(Temperature);
    let dp_profile = snd.get_profile(DewPoint);

    if t_profile.len().min(dp_profile.len()) == 0 {
        return vec![];
    }

    izip!(t_profile, dp_profile)
        .map(|(t_opt, dp_opt)| {
            t_opt.and_then(|t| {
                dp_opt.and_then(|dp| {
                    metfor::rh(t,dp)
                    // Ignore errors, if not possible to calculate just use missing value.
                    .ok()
                })
            })
        })
        .collect()
}

/// Given a sounding, calculate a profile of the potential temperature.
pub fn potential_temperature(snd: &Sounding) -> Vec<Option<f64>> {
    let p_profile = snd.get_profile(Pressure);
    let t_profile = snd.get_profile(Temperature);

    if p_profile.len().min(t_profile.len()) == 0 {
        return vec![];
    }

    izip!(p_profile, p_profile)
        .map(|(p_opt, t_opt)| {
            p_opt.and_then(|p| {
                t_opt.and_then(|t| {
                    metfor::theta_kelvin(p, t)
                    // Ignore errors, if not possible to calculate just use missing value.
                    .ok()
                })
            })
        })
        .collect()
}

/// Given a sounding, calculate a profile of the equivalent potential temperature.
pub fn equivalent_potential_temperature(snd: &Sounding) -> Vec<Option<f64>> {
    let p_profile = snd.get_profile(Pressure);
    let t_profile = snd.get_profile(Temperature);
    let dp_profile = snd.get_profile(DewPoint);

    if p_profile.len().min(t_profile.len()).min(dp_profile.len()) == 0 {
        return vec![];
    }

    izip!(p_profile, t_profile, dp_profile)
        .map(|(p_opt, t_opt, dp_opt)| {
            p_opt.and_then(|p| {
                t_opt.and_then(|t| {
                    dp_opt.and_then(|dp| {
                        metfor::theta_e_kelvin(t, dp, p)
                        // Ignore errors, if not possible to calculate just use missing value.
                        .ok()
                    })
                })
            })
        })
        .collect()
}

/// Get a profile of the lapse rate between layers in &deg;C / km.
pub fn temperature_lapse_rate(snd: &Sounding) -> Vec<Option<f64>> {
    let t_profile = snd.get_profile(Temperature).iter().cloned();
    lapse_rate(snd, t_profile)
}

/// Get the lapse rate of equivalent potential temperature in &deg;K / km.
pub fn theta_e_lapse_rate(snd: &Sounding) -> Vec<Option<f64>> {
    let theta_e = snd.get_profile(ThetaE).iter().cloned();
    lapse_rate(snd, theta_e)
}

fn lapse_rate<I: Iterator<Item = Option<f64>>>(snd: &Sounding, v_profile: I) -> Vec<Option<f64>> {
    let z_profile = snd.get_profile(GeopotentialHeight);

    izip!(z_profile, v_profile)
        .scan((None, None), |prev_pair, pair| {
            let &mut (ref mut prev_z, ref mut prev_v) = prev_pair;
            let (&z, v) = pair;

            let lapse_rate = if let (Some(ref prev_z), Some(ref prev_v), Some(ref z), Some(ref v)) =
                (*prev_z, *prev_v, z, v)
            {
                Some((v - prev_v) / (z - prev_z) * 1000.0)
            } else {
                None
            };

            *prev_z = z;
            *prev_v = v;

            Some(lapse_rate)
        })
        .collect()
}

/// Get the hydrolapse in (kg/kg)/km
pub fn hydrolapse(snd: &Sounding) -> Vec<Option<f64>> {
    let z_profile = snd.get_profile(GeopotentialHeight);
    let dp_profile = snd.get_profile(DewPoint);
    let p_profile = snd.get_profile(Pressure);

    izip!(p_profile, z_profile, dp_profile)
        .scan(
            (None, None),
            |prev_pair: &mut (Option<f64>, Option<f64>), triple| {
                let &mut (ref mut prev_z, ref mut prev_mw) = prev_pair;
                let (&p, &z, &dp) = triple;

                let mw = if let (Some(p), Some(dp)) = (p, dp) {
                    ::metfor::mixing_ratio(dp, p).ok()
                } else {
                    None
                };

                let mw_lapse_rate = if let (Some(p_z), Some(p_mw), Some(z), Some(mw)) =
                    (*prev_z, *prev_mw, z, mw)
                {
                    Some((mw - p_mw) / (z - p_z) * 1000.0)
                } else {
                    None
                };

                *prev_z = z;
                *prev_mw = mw;

                Some(mw_lapse_rate)
            },
        )
        .collect()
}

#[cfg(test)]
mod test_sounding_profiles {
    use super::*;

    fn make_test_sounding() -> Sounding {
        Sounding::new()
            .set_profile(Temperature, vec![Some(9.8), Some(0.0), Some(-5.0)])
            .set_profile(
                GeopotentialHeight,
                vec![Some(1000.0), Some(2000.0), Some(3000.0)],
            )
    }

    #[test]
    fn test_temperature_lapse_rate() {
        let snd = make_test_sounding();

        let lapse_rate = temperature_lapse_rate(&snd);
        println!("{:#?}", lapse_rate);
        assert!(lapse_rate.contains(&Some(-9.8)));
        assert!(lapse_rate.contains(&Some(-5.0)));
    }
}

/// Hold profiles for a parcel and it's environment.
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

/// Lift a parcel, make sure it has the same data points as the reference sounding.
pub fn lift_parcel(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let (lcl_pressure, lcl_temperature) = metfor::pressure_and_temperature_at_lcl(
        parcel.temperature,
        parcel.dew_point,
        parcel.pressure,
    ).map_err(|_| AnalysisError::InvalidInput)?;

    let theta = parcel.theta()?;
    let theta_e = parcel.theta_e()?;

    let lcl_env = ::interpolation::linear_interpolate(snd, lcl_pressure)?;
    let lcl_height = lcl_env.height.ok_or(AnalysisError::InvalidInput)?;
    let lcl_env_temperature = lcl_env.temperature.ok_or(AnalysisError::InvalidInput)?;

    let pressure = snd.get_profile(Pressure);
    let hgt = snd.get_profile(GeopotentialHeight);
    let env_t = snd.get_profile(Temperature);

    let dry_adiabatic = izip!(pressure, hgt, env_t)
        // Skip levels below the parcel starting point
        .skip_while(|row| {
            if let Some(p) = *row.0 {
                p > parcel.pressure
            } else {
                true
            }
        })
        // Take until we reach the lcl and become saturated
        .take_while(|row|{
            if let Some(p) = *row.0 {
                p > lcl_pressure
            } else {
                true
            }
        })
        // Remove layers with missing data and unwrap options
        .filter_map(|(p_opt, h_opt, t_opt)|{
            p_opt.and_then(|p| h_opt.and_then(|h| t_opt.and_then(|t| Some((p,h,t)))))
        })
        // Calculate the parcel temperature and add it in, if it is OK
        .filter_map(|(p,h,env_t)|{
            metfor::temperature_c_from_theta(theta, p)
                .ok()
                .and_then(|parcel_t| Some((p,h,parcel_t, env_t)))
        });

    let lcl = once((
        lcl_pressure,
        lcl_height,
        lcl_temperature,
        lcl_env_temperature,
    ));

    let moist_adiabatic = izip!(pressure, hgt, env_t)
        // Skip levels below the parcel starting point and below the lcl
        .skip_while(|row| {
            if let Some(p) = *row.0 {
                p > lcl_pressure || p > parcel.pressure
            } else {
                true
            }
        })
        // Remove layers with missing data and unwrap options
        .filter_map(|(p_opt, h_opt, t_opt)|{
            p_opt.and_then(|p| h_opt.and_then(|h| t_opt.and_then(|t| Some((p,h,t)))))
        })
        // Calculate the parcel temperature and add it in, if it is OK
        .filter_map(|(p,h,env_t)|{
            metfor::temperature_c_from_theta_e_saturated_and_pressure(p, theta_e)
                .ok()
                .and_then(|parcel_t| Some((p,h,parcel_t, env_t)))
        });

    let full_profile = dry_adiabatic.chain(lcl).chain(moist_adiabatic);

    let mut pressure = Vec::with_capacity(pressure.len() + 1);
    let mut height = Vec::with_capacity(pressure.len() + 1);
    let mut parcel_t = Vec::with_capacity(pressure.len() + 1);
    let mut environment_t = Vec::with_capacity(pressure.len() + 1);

    full_profile.for_each(|(p, h, p_t, e_t)| {
        pressure.push(p);
        height.push(h);
        parcel_t.push(p_t);
        environment_t.push(e_t);
    });

    Ok(ParcelProfile {
        pressure,
        height,
        parcel_t,
        environment_t,
        parcel,
    })
}

// TODO: descend parcel dry adiabatically
// TODO: descend parcel moist adiabatically
