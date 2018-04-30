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

use smallvec::SmallVec;

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

impl ParcelProfile {
    /// Get the lfcs and el levels for this parcel.
    pub fn cape_layers(&self) -> SmallVec<[(f64, f64); ::VEC_SIZE]> {
        let mut to_ret = SmallVec::new();

        let lcl_pressure = if let Ok(pres) = metfor::pressure_hpa_at_lcl(
            self.parcel.temperature,
            self.parcel.dew_point,
            self.parcel.pressure,
        ) {
            pres
        } else {
            return to_ret;
        };

        let l0 = izip!(&self.pressure, &self.parcel_t, &self.environment_t);
        let l1 = izip!(&self.parcel_t, &self.environment_t).skip(1);

        let mut bottom = ::std::f64::MIN;
        let mut top = ::std::f64::MAX;

        if !(self.parcel_t.is_empty() || self.environment_t.is_empty() || self.pressure.is_empty()){
            if self.parcel_t[0] >= self.environment_t[0] {
                bottom = self.pressure[0];
            }
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
                    to_ret.push((bottom, top));
                    bottom = ::std::f64::MIN;
                    top = ::std::f64::MAX;
                }
            });

        to_ret
    }

    /// Get the pressure layers with CIN in this profile.
    pub fn cin_layers(&self) -> SmallVec<[(f64, f64); ::VEC_SIZE]> {
        let mut to_ret = SmallVec::new();

        let l0 = izip!(&self.pressure, &self.parcel_t, &self.environment_t);
        let l1 = izip!(&self.parcel_t, &self.environment_t).skip(1);

        let mut top = ::std::f64::MAX;
        let mut bottom = {
            let (p, pt, et) = (self.pressure[0], self.parcel_t[0], self.environment_t[0]);
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
                to_ret.push((bottom, top));
                bottom = ::std::f64::MIN;
                top = ::std::f64::MAX;
            }
        });

        to_ret
    }
}

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

    let parcel_start_data = ::interpolation::linear_interpolate(snd, parcel.pressure)?;

    let lcl_env = ::interpolation::linear_interpolate(snd, lcl_pressure)?;
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
                let tgt_row = ::interpolation::linear_interpolate(snd, tgt_pres)?;
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

// TODO: descend parcel dry adiabatically
// TODO: descend parcel moist adiabatically
