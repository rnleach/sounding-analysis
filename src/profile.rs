//! Create profiles.
//!
//! The output will be at the same levels as the sounding and these are suitable to be set as a
//! profile in the sounding. For example, calculating a wet bulb or relative humidity profile from a
//! sounding with temperature and dew point. If one of the profiles required for the analysis in the
//! sounding is missing, the result cannot be calculated and an empty vector is returned.
//!

use metfor;
use optional::{none, some, Optioned};
use sounding_base::Profile::{DewPoint, GeopotentialHeight, Pressure, Temperature, ThetaE};
use sounding_base::Sounding;
use sounding_base::Surface;

/// Given a sounding, calculate a profile of wet bulb temperature.
pub fn wet_bulb(snd: &Sounding) -> Vec<Optioned<f64>> {
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
                            .into()
                    })
                })
            })
        })
        .collect()
}

/// Given a sounding, calculate a profile of relative humidity.
pub fn relative_humidity(snd: &Sounding) -> Vec<Optioned<f64>> {
    let t_profile = snd.get_profile(Temperature);
    let dp_profile = snd.get_profile(DewPoint);

    if t_profile.len().min(dp_profile.len()) == 0 {
        return vec![];
    }

    izip!(t_profile, dp_profile)
        .map(|(t_opt, dp_opt)| {
            t_opt.and_then(|t| {
                dp_opt.and_then(|dp| {
                    metfor::rh(t, dp)
                        // Ignore errors, if not possible to calculate just use missing value.
                        .ok()
                        .into()
                })
            })
        })
        .collect()
}

/// Given a sounding, calculate a profile of the potential temperature.
pub fn potential_temperature(snd: &Sounding) -> Vec<Optioned<f64>> {
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
                        .into()
                })
            })
        })
        .collect()
}

/// Given a sounding, calculate a profile of the equivalent potential temperature.
pub fn equivalent_potential_temperature(snd: &Sounding) -> Vec<Optioned<f64>> {
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
                            .into()
                    })
                })
            })
        })
        .collect()
}

/// Get a profile of the lapse rate between layers in &deg;C / km.
pub fn temperature_lapse_rate(snd: &Sounding) -> Vec<Optioned<f64>> {
    let t_profile = snd.get_profile(Temperature).iter().cloned();
    lapse_rate(snd, t_profile)
}

/// Get a profile of the average lapse rate from the surface to *, or the level on the y axis.
pub fn sfc_to_level_temperature_lapse_rate(snd: &Sounding) -> Vec<Optioned<f64>> {
    let z_profile = snd.get_profile(GeopotentialHeight);
    let t_profile = snd.get_profile(Temperature);

    let (t_sfc, z_sfc): (f64, f64) = if let (Some(t_sfc), Some(z_sfc)) = (
        snd.get_surface_value(Surface::Temperature).into(),
        snd.get_station_info().elevation().into(),
    ) {
        (t_sfc, z_sfc)
    } else {
        return vec![];
    };

    izip!(z_profile, t_profile)
        .map(|(z_opt, t_opt)| {
            let z_opt: Option<f64> = (*z_opt).into();
            let t_opt: Option<f64> = (*t_opt).into();
            if let (Some(z), Some(t)) = (z_opt, t_opt) {
                if (z - z_sfc).abs() < ::std::f64::EPSILON {
                    none()
                } else {
                    some((t - t_sfc) / (z - z_sfc) * 1000.0)
                }
            } else {
                none()
            }
        })
        .collect()
}

/// Get the lapse rate of equivalent potential temperature in &deg;K / km.
pub fn theta_e_lapse_rate(snd: &Sounding) -> Vec<Optioned<f64>> {
    let theta_e = snd.get_profile(ThetaE).iter().cloned();
    lapse_rate(snd, theta_e)
}

fn lapse_rate<I: Iterator<Item = Optioned<f64>>>(
    snd: &Sounding,
    v_profile: I,
) -> Vec<Optioned<f64>> {
    let z_profile = snd.get_profile(GeopotentialHeight);

    izip!(z_profile, v_profile)
        .scan((None, None), |prev_pair, (&z, v)| {
            let &mut (ref mut prev_z, ref mut prev_v) = prev_pair;

            let z: Option<f64> = z.into();
            let v: Option<f64> = v.into();

            let lapse_rate = if let (Some(ref prev_z), Some(ref prev_v), Some(ref z), Some(ref v)) =
                (*prev_z, *prev_v, z, v)
            {
                some((v - prev_v) / (z - prev_z) * 1000.0)
            } else {
                none()
            };

            *prev_z = z;
            *prev_v = v;

            Some(lapse_rate)
        })
        .collect()
}

/// Get the hydrolapse in (kg/kg)/km
pub fn hydrolapse(snd: &Sounding) -> Vec<Optioned<f64>> {
    let z_profile = snd.get_profile(GeopotentialHeight);
    let dp_profile = snd.get_profile(DewPoint);
    let p_profile = snd.get_profile(Pressure);

    izip!(p_profile, z_profile, dp_profile)
        .scan((None, None), |prev_pair, (&p, &z, &dp)| {
            let &mut (ref mut prev_z, ref mut prev_mw) = prev_pair;

            let p: Option<f64> = p.into();
            let z: Option<f64> = z.into();
            let dp: Option<f64> = dp.into();

            let mw = if let (Some(p), Some(dp)) = (p, dp) {
                ::metfor::mixing_ratio(dp, p).ok()
            } else {
                None
            };

            let mw_lapse_rate =
                if let (Some(p_z), Some(p_mw), Some(z), Some(mw)) = (*prev_z, *prev_mw, z, mw) {
                    some((mw - p_mw) / (z - p_z) * 1000.0)
                } else {
                    none()
                };

            *prev_z = z;
            *prev_mw = mw;

            Some(mw_lapse_rate)
        })
        .collect()
}

#[cfg(test)]
mod test_sounding_profiles {
    use super::*;

    fn make_test_sounding() -> Sounding {
        Sounding::new()
            .set_profile(Temperature, vec![some(9.8), some(0.0), some(-5.0)])
            .set_profile(
                GeopotentialHeight,
                vec![some(1000.0), some(2000.0), some(3000.0)],
            )
    }

    #[test]
    fn test_temperature_lapse_rate() {
        let snd = make_test_sounding();

        let lapse_rate = temperature_lapse_rate(&snd);
        println!("{:#?}", lapse_rate);
        assert!(lapse_rate.contains(&some(-9.8)));
        assert!(lapse_rate.contains(&some(-5.0)));
    }
}
