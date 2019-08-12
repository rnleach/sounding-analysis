//! Create profiles.
//!
//! The output will be at the same levels as the sounding and these are suitable to be set as a
//! profile in the sounding. For example, calculating a wet bulb or relative humidity profile from a
//! sounding with temperature and dew point. If one of the profiles required for the analysis in the
//! sounding is missing, the result cannot be calculated and an empty vector is returned.
//!
use crate::sounding::Sounding;
use itertools::izip;
use metfor::{
    self, Celsius, CelsiusDiff, CelsiusPKm, HydrolapsePKm, Kelvin, KelvinPKm, Km, Meters, Quantity,
    Temperature,
};
use optional::{none, some, Optioned};
use std::ops::Sub;

/// Given a sounding, calculate a profile of wet bulb temperature.
pub fn wet_bulb(snd: &Sounding) -> Vec<Optioned<Celsius>> {
    let p_profile = snd.pressure_profile();
    let t_profile = snd.temperature_profile();
    let dp_profile = snd.dew_point_profile();

    if p_profile.len().min(t_profile.len()).min(dp_profile.len()) == 0 {
        return vec![];
    }

    izip!(p_profile, t_profile, dp_profile)
        .map(|(p_opt, t_opt, dp_opt)| {
            p_opt.and_then(|p| {
                t_opt.and_then(|t| {
                    dp_opt.and_then(|dp| Optioned::<Celsius>::from(metfor::wet_bulb(t, dp, p)))
                })
            })
        })
        .collect()
}

/// Given a sounding, calculate a profile of relative humidity.
pub fn relative_humidity(snd: &Sounding) -> Vec<Optioned<f64>> {
    let t_profile = snd.temperature_profile();
    let dp_profile = snd.dew_point_profile();

    if t_profile.len().min(dp_profile.len()) == 0 {
        return vec![];
    }

    izip!(t_profile, dp_profile)
        .map(|(t_opt, dp_opt)| {
            t_opt.and_then(|t| dp_opt.and_then(|dp| Optioned::from(metfor::rh(t, dp))))
        })
        .collect()
}

/// Given a sounding, calculate a profile of relative humidity with respect to ice.
pub fn relative_humidity_ice(snd: &Sounding) -> Vec<Optioned<f64>> {
    let t_profile = snd.temperature_profile();
    let dp_profile = snd.dew_point_profile();

    if t_profile.len().min(dp_profile.len()) == 0 {
        return vec![];
    }

    izip!(t_profile, dp_profile)
        .map(|(t_opt, dp_opt)| {
            t_opt.and_then(|t| {
                dp_opt.and_then(|dp| {
                    // Use dew point to get vapor pressure of water - sounding more likely to have
                    // dew point information.
                    let vp_water = metfor::vapor_pressure_liquid_water(dp);
                    // Get the saturation vapor pressure relative to ice
                    let vp_sat_ice = metfor::vapor_pressure_ice(t);

                    let rh = vp_water.and_then(|vpw| vp_sat_ice.map(|vpsat| vpw / vpsat));
                    Optioned::from(rh)
                })
            })
        })
        .collect()
}

/// Given a sounding, calculate a profile of the potential temperature.
pub fn potential_temperature(snd: &Sounding) -> Vec<Optioned<Kelvin>> {
    let p_profile = snd.pressure_profile();
    let t_profile = snd.temperature_profile();

    if p_profile.len().min(t_profile.len()) == 0 {
        return vec![];
    }

    izip!(p_profile, t_profile)
        .map(|(p_opt, t_opt)| {
            p_opt.and_then(|p| t_opt.and_then(|t| Optioned::<Kelvin>::from(metfor::theta(p, t))))
        })
        .collect()
}

/// Given a sounding, calculate a profile of the equivalent potential temperature.
pub fn equivalent_potential_temperature(snd: &Sounding) -> Vec<Optioned<Kelvin>> {
    let p_profile = snd.pressure_profile();
    let t_profile = snd.temperature_profile();
    let dp_profile = snd.dew_point_profile();

    if p_profile.len().min(t_profile.len()).min(dp_profile.len()) == 0 {
        return vec![];
    }

    izip!(p_profile, t_profile, dp_profile)
        .map(|(p_opt, t_opt, dp_opt)| {
            p_opt.and_then(|p| {
                t_opt.and_then(|t| {
                    dp_opt.and_then(|dp| Optioned::<Kelvin>::from(metfor::theta_e(t, dp, p)))
                })
            })
        })
        .collect()
}

/// Get a profile of the lapse rate between layers in &deg;C / km.
pub fn temperature_lapse_rate(snd: &Sounding) -> Vec<Optioned<CelsiusPKm>> {
    let t_profile = snd.temperature_profile().iter().cloned();
    lapse_rate(snd, t_profile)
}

/// Get a profile of the average lapse rate from the surface to *, or the level on the y axis.
pub fn sfc_to_level_temperature_lapse_rate(snd: &Sounding) -> Vec<Optioned<CelsiusPKm>> {
    let z_profile = snd.height_profile();
    let t_profile = snd.temperature_profile();

    let (t_sfc, z_sfc): (Celsius, Meters) = if let (Some(t_sfc), Some(z_sfc)) = (
        snd.sfc_temperature().into(),
        snd.station_info().elevation().into(),
    ) {
        (t_sfc, z_sfc)
    } else {
        return vec![];
    };

    izip!(z_profile, t_profile)
        .map(|(z_opt, t_opt)| {
            if let (Some(z), Some(t)) = (z_opt.into_option(), t_opt.into_option()) {
                if (z - z_sfc).unpack().abs() < ::std::f64::EPSILON {
                    none()
                } else {
                    let CelsiusDiff(dt) = t - t_sfc;
                    let Km(dz) = Km::from(z - z_sfc);
                    some(CelsiusPKm(dt / dz))
                }
            } else {
                none()
            }
        })
        .collect()
}

/// Get the lapse rate of equivalent potential temperature in &deg;K / km.
pub fn theta_e_lapse_rate(snd: &Sounding) -> Vec<Optioned<KelvinPKm>> {
    let theta_e = snd.theta_e_profile().iter().cloned();
    lapse_rate(snd, theta_e)
}

fn lapse_rate<I, T>(snd: &Sounding, v_profile: I) -> Vec<Optioned<CelsiusPKm>>
where
    I: Iterator<Item = Optioned<T>>,
    T: Temperature + Sub<T> + optional::Noned,
    CelsiusDiff: From<<T as Sub<T>>::Output>,
{
    let z_profile = snd.height_profile();

    izip!(z_profile, v_profile)
        .scan((None, None), |prev_pair, (&z, v)| {
            let &mut (ref mut prev_z, ref mut prev_v) = prev_pair;

            let z: Option<Meters> = z.into_option();
            let v: Option<T> = v.into_option();

            let lapse_rate = if let (Some(ref prev_z), Some(ref prev_v), Some(ref z), Some(ref v)) =
                (*prev_z, *prev_v, z, v)
            {
                let CelsiusDiff(dt) = CelsiusDiff::from(*v - *prev_v);
                let Km(dz) = Km::from(*z - *prev_z);
                some(CelsiusPKm(dt / dz))
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
pub fn hydrolapse(snd: &Sounding) -> Vec<Optioned<HydrolapsePKm>> {
    let z_profile = snd.height_profile();
    let dp_profile = snd.dew_point_profile();
    let p_profile = snd.pressure_profile();

    izip!(p_profile, z_profile, dp_profile)
        .scan((None, None), |prev_pair, (&p, &z, &dp)| {
            let &mut (ref mut prev_z, ref mut prev_mw) = prev_pair;

            let p: Option<_> = p.into_option();
            let z: Option<_> = z.into_option();
            let dp: Option<_> = dp.into_option();

            let mw = if let (Some(p), Some(dp)) = (p, dp) {
                metfor::mixing_ratio(dp, p)
            } else {
                None
            };

            let mw_lapse_rate =
                if let (Some(p_z), Some(p_mw), Some(z), Some(mw)) = (*prev_z, *prev_mw, z, mw) {
                    let dmw = mw - p_mw;
                    let Km(dz) = Km::from(z - p_z);
                    some(HydrolapsePKm(dmw / dz))
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
    use metfor::{Celsius, Meters};

    fn make_test_sounding() -> Sounding {
        Sounding::new()
            .with_temperature_profile(vec![
                some(Celsius(9.8)),
                some(Celsius(0.0)),
                some(Celsius(-5.0)),
            ])
            .with_height_profile(vec![
                some(Meters(1000.0)),
                some(Meters(2000.0)),
                some(Meters(3000.0)),
            ])
    }

    #[test]
    fn test_temperature_lapse_rate() {
        let snd = make_test_sounding();

        let lapse_rate = temperature_lapse_rate(&snd);
        println!("{:#?}", lapse_rate);
        assert!(lapse_rate.contains(&some(CelsiusPKm(-9.8))));
        assert!(lapse_rate.contains(&some(CelsiusPKm(-5.0))));
    }
}
