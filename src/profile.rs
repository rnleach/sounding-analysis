//! Create profiles.
//!
//! The output will be at the same levels as the sounding and these are suitable to be set as a
//! profile in the sounding. For example, calculating a wet bulb or relative humidity profile from a
//! sounding with temperature and dew point. If one of the profiles required for the analysis in the
//! sounding is missing, the result cannot be calculated and an empty vector is returned.
//!

use metfor;
use sounding_base::Profile::{DewPoint, GeopotentialHeight, Pressure, Temperature, ThetaE};
use sounding_base::Sounding;
use sounding_base::Surface;

use parcel::mixed_layer_parcel;

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

/// Get a profile of the average lapse rate from the surface to *, or the level on the y axis.
pub fn sfc_to_level_temperature_lapse_rate(snd: &Sounding) -> Vec<Option<f64>> {
    let z_profile = snd.get_profile(GeopotentialHeight);
    let t_profile = snd.get_profile(Temperature);

    let (t_sfc, z_sfc) = if let (Some(t_sfc), Some(z_sfc)) = (
        snd.get_surface_value(Surface::Temperature),
        snd.get_station_info().elevation(),
    ) {
        (t_sfc, z_sfc)
    } else {
        return vec![];
    };

    izip!(z_profile, t_profile)
        .map(|pair| {
            if let (Some(z), Some(t)) = pair {
                if (*z - z_sfc).abs() < ::std::f64::EPSILON {
                    None
                } else {
                    Some((*t - t_sfc) / (*z - z_sfc) * 1000.0)
                }
            } else {
                None
            }
        })
        .collect()
}

/// Get a profile of the average lapse rate from the surface to *, or the level on the y axis.
pub fn ml_to_level_temperature_lapse_rate(snd: &Sounding) -> Vec<Option<f64>> {
    let z_profile = snd.get_profile(GeopotentialHeight);
    let t_profile = snd.get_profile(Temperature);

    let (t_sfc, z_sfc) = match mixed_layer_parcel(snd) {
        Ok(parcel) => {
            if let Some(z_sfc) = snd.get_station_info().elevation() {
                (parcel.temperature, z_sfc)
            } else {
                return vec![];
            }
        }
        Err(_) => return vec![],
    };

    izip!(z_profile, t_profile)
        .map(|pair| {
            if let (Some(z), Some(t)) = pair {
                if (*z - z_sfc).abs() < ::std::f64::EPSILON {
                    None
                } else {
                    Some((*t - t_sfc) / (*z - z_sfc) * 1000.0)
                }
            } else {
                None
            }
        })
        .collect()
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
