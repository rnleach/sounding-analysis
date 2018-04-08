//! Create profiles.
//!
//! There are three kinds of profiles:
//!   - Those created from a sounding, the output will be at the same levels as the sounding and
//!     these are suitable to be set as a profile in the sounding. For example, calculating a
//!     wet bulb or relative humidity profile from a sounding with temperature and dew point. If one
//!     of the profiles required for the analysis in the sounding is missing, the result cannot be
//!     calculated and an empty vector is returned.
//!   - Those created for parcel analysis on a sounding. These require an initial parcel and a
//!     sounding. They are useful doing parcel analysis for variables such CAPE.
//!   - Those created for reference or perhaps drawing. These require a top and bottom pressure,
//!     constant value (e.g. potential temperature for a dry adiabat), and a buffer (or slice) to
//!     fill.
//!
// TODO: EXAMPLES HERE.

use metfor;
use sounding_base::Sounding;
use sounding_base::Profile::*;

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


// TODO: Richardson Number
// TODO: lapse rate
// TODO: theta-e lapse rate

// TODO: Dry adiabat
// TODO: Moist adiabat
// TODO: Constant MW
