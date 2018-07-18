//! Create and analyze a profile from lifting or descending a parcel.

use metfor;
use smallvec::SmallVec;
use sounding_base::{Profile::*, Sounding};

use error::*;
use interpolation::{linear_interp, linear_interpolate_sounding};
use parcel::Parcel;

/// Hold profiles for a parcel and it's environment.
#[derive(Debug, Clone)]
pub struct ParcelProfile {
    /// Pressure profile
    pub pressure: Vec<f64>,
    /// Height profile
    pub height: Vec<f64>,
    /// Parcel temperature profile
    pub parcel_t: Vec<f64>,
    /// Parcel Mixing ratio
    pub parcel_mw: Vec<f64>,
    /// Environment temperature profile
    pub environment_t: Vec<f64>,
    /// Environment mixing ratio
    pub environment_mw: Vec<f64>,
    /// The original parcel
    pub parcel: Parcel,
}

impl ParcelProfile {}

/// Lift a parcel
pub fn lift_parcel(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    //
    // Find the LCL
    //
    let (lcl_pressure, lcl_temperature) = metfor::pressure_and_temperature_at_lcl(
        parcel.temperature,
        parcel.dew_point,
        parcel.pressure,
    )?;

    let lcl_temperature = metfor::kelvin_to_celsius(lcl_temperature)?;
    let lcl_env = linear_interpolate_sounding(snd, lcl_pressure)?;
    let lcl_height = lcl_env.height.ok_or(AnalysisError::InvalidInput)?;
    let lcl_env_temperature = lcl_env.temperature.ok_or(AnalysisError::InvalidInput)?;
    let lcl_env_mw = lcl_env
        .dew_point
        .ok_or(AnalysisError::InvalidInput)
        .and_then(|dp| Ok(metfor::mixing_ratio(dp, lcl_pressure)?))?;

    //
    // The starting level to lift the parcel from
    //
    let parcel_start_data = linear_interpolate_sounding(snd, parcel.pressure)?;

    //
    // How to calculate a parcel temperature for a given pressure level
    //
    let theta = parcel.theta()?;
    let theta_e = parcel.theta_e()?;
    let dry_mw = parcel.mixing_ratio()?;
    let calc_parcel_t = |tgt_pres| {
        if tgt_pres > lcl_pressure {
            // Dry adiabatic lifting
            metfor::temperature_c_from_theta(theta, tgt_pres)
        } else {
            // Moist adiabatic lifting
            metfor::temperature_c_from_theta_e_saturated_and_pressure(tgt_pres, theta_e)
        }
    };

    let calc_parcel_mw = |tgt_pres, pcl_t| {
        if tgt_pres > lcl_pressure {
            // Dry adiabatic lifting
            Ok(dry_mw)
        } else {
            // Moist adiabatic lifting
            metfor::mixing_ratio(pcl_t, tgt_pres)
        }
    };

    //
    // Get the environment data to iterate over. We want the parcel profile to have all the same
    // pressure levels as the environmental sounding, plus a few special ones.
    //
    let snd_pressure = snd.get_profile(Pressure);
    let hgt = snd.get_profile(GeopotentialHeight);
    let env_t = snd.get_profile(Temperature);
    let env_dp = snd.get_profile(DewPoint);

    //
    // Allocate some buffers to hold the return values.
    //
    let mut pressure = Vec::with_capacity(snd_pressure.len() + 5);
    let mut height = Vec::with_capacity(snd_pressure.len() + 5);
    let mut parcel_t = Vec::with_capacity(snd_pressure.len() + 5);
    let mut parcel_mw = Vec::with_capacity(snd_pressure.len() + 5);
    let mut environment_t = Vec::with_capacity(snd_pressure.len() + 5);
    let mut environment_mw = Vec::with_capacity(snd_pressure.len() + 5);

    // Nested scope to limit closure borrows
    {
        // Helper function to add row to parcel profile
        let mut add_row = |pp, hh, pcl_tt, pcl_mw, env_tt, env_mw| {
            pressure.push(pp);
            height.push(hh);
            parcel_t.push(pcl_tt);
            parcel_mw.push(pcl_mw);
            environment_t.push(env_tt);
            environment_mw.push(env_mw);
        };

        // Start by adding the parcel level
        let mut p0 = parcel.pressure;
        let h0 = parcel_start_data.height.ok_or(AnalysisError::InvalidInput)?;
        let mut pcl_t0 = parcel.temperature;
        let pcl_mw0 = parcel.virtual_temperature_c()?;
        let mut env_t0 = parcel_start_data
            .temperature
            .ok_or(AnalysisError::InvalidInput)?;
        let env_mw0 = parcel_start_data
            .dew_point
            .ok_or(AnalysisError::InvalidInput)
            .and_then(|dp| Ok(metfor::virtual_temperature_c(env_t0, dp, p0)?))?;

        add_row(p0, h0, pcl_t0, pcl_mw0, env_t0, env_mw0);

        //
        // Construct an iterator that calculates the parcel values
        //
        let iter = izip!(snd_pressure, hgt, env_t, env_dp)
            // Remove rows with missing data and unpack options
            .filter_map(|(p, h, env_t, env_dp)|{
                if p.is_some() && h.is_some() && env_t.is_some() && env_dp.is_some(){
                    Some((p.unpack(), h.unpack(), env_t.unpack(), env_dp.unpack()))
                } else {
                    None
                }
            })
            // Remove rows at or below the parcel level
            .filter(move |(p, _, _, _)| *p < p0)
            // Calculate the parcel temperature, skip this level if there is an error
            .filter_map(|(p, h, env_t, env_dp)| {
                match calc_parcel_t(p) {
                    Ok(pcl_t) => Some((p, h, env_t, env_dp, pcl_t)),
                    Err(_) => None,
                }
            })
            // Calculate the parcel mixing ratio
            .filter_map(|(p, h, env_t, env_dp, pcl_t)|{
                match calc_parcel_mw(p, pcl_t) {
                    Ok(pcl_mw) => Some((p, h, env_t, env_dp, pcl_t, pcl_mw)),
                    Err(_) => None,
                }
            })
            // Calculate the environment mixing ratio
            .filter_map(|(p, h, env_t, env_dp, pcl_t, pcl_mw)|{
                match metfor::mixing_ratio(env_dp, p) {
                    Ok(env_mw) => Some((p, h, env_t, env_mw, pcl_t, pcl_mw)),
                    Err(_) => None,
                }
            });

        //
        // Pack the resulting values into their vectors.
        //
        for (p, h, env_t, env_mw, pcl_t, pcl_mw) in iter {
            // Check to see if we are passing the lcl
            if p0 > lcl_pressure && p < lcl_pressure {
                add_row(
                    lcl_pressure,
                    lcl_height,
                    lcl_temperature,
                    dry_mw,
                    lcl_env_temperature,
                    lcl_env_mw,
                );
            }

            // Check to see if the parcel and environment soundings have crossed
            if (pcl_t0 < env_t0 && pcl_t > env_t) || (pcl_t0 > env_t0 && pcl_t < env_t) {
                let tgt_pres = linear_interp(0.0, pcl_t - env_t, pcl_t0 - env_t0, p, p0);
                let tgt_row = linear_interpolate_sounding(snd, tgt_pres)?;
                let h2 = tgt_row.height.ok_or(AnalysisError::InvalidInput)?;
                let env_t2 = tgt_row.temperature.ok_or(AnalysisError::InvalidInput)?;
                let env_mw2 = tgt_row
                    .dew_point
                    .ok_or(AnalysisError::InvalidInput)
                    .and_then(|dp| Ok(metfor::mixing_ratio(dp, tgt_pres)?))?;
                let pcl_mw2 = if tgt_pres > lcl_pressure {
                    dry_mw
                } else {
                    metfor::mixing_ratio(env_t2, tgt_pres)?
                };

                add_row(tgt_pres, h2, env_t2, pcl_mw2, env_t2, env_mw2);
            }

            // Add the new values to the array
            add_row(p, h, pcl_t, pcl_mw, env_t, env_mw);

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
        parcel_mw,
        environment_t,
        environment_mw,
        parcel,
    })
}

/// Get the lfcs and el levels for a parcel.
///
/// Returns a list of (bottom_pressure, top_pressure) pairs.
pub fn cape_layers(parcel_profile: &ParcelProfile) -> SmallVec<[(f64, f64); ::VEC_SIZE]> {
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

    if !(parcel_profile.parcel_t.is_empty()
        || parcel_profile.environment_t.is_empty()
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
                to_ret.push((bottom, top));
                bottom = ::std::f64::MIN;
                top = ::std::f64::MAX;
            }
        });

    to_ret
}

/// Get the pressure layers with CIN in this profile.
///
/// Returns a list of (bottom_pressure, top_pressure) pairs.
pub fn cin_layers(parcel_profile: &ParcelProfile) -> SmallVec<[(f64, f64); ::VEC_SIZE]> {
    let mut to_ret = SmallVec::new();

    let l0 = izip!(
        &parcel_profile.pressure,
        &parcel_profile.parcel_t,
        &parcel_profile.environment_t
    );
    let l1 = izip!(&parcel_profile.parcel_t, &parcel_profile.environment_t).skip(1);

    let mut top = ::std::f64::MAX;
    let mut bottom = {
        let (p0, pt0, et0) = (
            parcel_profile.pressure[0],
            parcel_profile.parcel_t[0],
            parcel_profile.environment_t[0],
        );
        let (pt1, et1) = (parcel_profile.parcel_t[1], parcel_profile.environment_t[1]);
        if pt0 <= et0 && pt1 < et1 {
            p0
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

/// Descend a parcel dry adiabatically.
pub fn descend_dry(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let theta = parcel.theta()?;
    descend_parcel(
        parcel,
        snd,
        theta,
        ::metfor::temperature_c_from_theta,
        false,
    )
}

/// Descend a parcel moist adiabatically
pub fn descend_moist(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let theta = parcel.theta_e()?;

    let theta_func = |theta_e, press| {
        ::metfor::temperature_c_from_theta_e_saturated_and_pressure(press, theta_e)
    };

    descend_parcel(parcel, snd, theta, theta_func, true)
}

#[inline]
fn descend_parcel<F>(
    parcel: Parcel,
    snd: &Sounding,
    theta: f64,
    theta_func: F,
    saturated_descent: bool,
) -> Result<ParcelProfile>
where
    F: Fn(f64, f64) -> ::std::result::Result<f64, ::metfor::MetForErr>,
{
    let dry_mw = parcel.mixing_ratio()?;

    let mut pressure = Vec::new();
    let mut height = Vec::new();
    let mut parcel_t = Vec::new();
    let mut parcel_mw = Vec::new();
    let mut environment_t = Vec::new();
    let mut environment_mw = Vec::new();

    // Actually start at the bottom and work up.
    let press = snd.get_profile(Pressure);
    let env_t = snd.get_profile(Temperature);
    let env_dp = snd.get_profile(DewPoint);
    let hght = snd.get_profile(GeopotentialHeight);

    // Nested scope to limit borrows
    {
        // Helper function to add row to parcel profile
        let mut add_row = |pp, hh, pcl_tt, pcl_mw, env_tt, env_mw| {
            pressure.push(pp);
            height.push(hh);
            parcel_t.push(pcl_tt);
            parcel_mw.push(pcl_mw);
            environment_t.push(env_tt);
            environment_mw.push(env_mw);
        };

        izip!(press, hght, env_t, env_dp)
            .take_while(|(p_opt, _, _, _)| {
                if p_opt.is_some() {
                    p_opt.unpack() >= parcel.pressure
                } else {
                    true
                }
            })
            // Remove levels with missing data
            .filter_map(|(p_opt, h_opt, e_t_opt, e_dp_opt)| {
                if p_opt.is_some() && h_opt.is_some() && e_t_opt.is_some() && e_dp_opt.is_some() {
                    Some((p_opt.unpack(), h_opt.unpack(), e_t_opt.unpack(), e_dp_opt.unpack()))
                } else {
                    None
                }
            })
            // Get the environmental mixing ratio
            .filter_map(|(p, h, et, edp)|{
                match metfor::mixing_ratio(edp, p) {
                    Ok(env_mw) => Some((p, h, et, env_mw)),
                    Err(_) => None,
                }
            })
            // Get the parcel temperature
            .filter_map(|(p, h, et, e_mw)|{
                match theta_func(theta, p) {
                    Ok(pcl_t) => Some((p, h, pcl_t, et, e_mw)),
                    Err(_) => None,
                }
            })
            // Get the parcel mixing ratio
            .filter_map(|(p, h, pcl_t, et, e_mw)|{
                if saturated_descent {
                    match metfor::mixing_ratio(pcl_t, p) {
                        Ok(pcl_mw) => Some((p, h, pcl_t, pcl_mw, et, e_mw)),
                        Err(_) => None,
                    }
                } else {
                    Some((p, h, pcl_t, dry_mw, et, e_mw))
                }
            })
            .for_each(|(p, h, pt, pmw, et, e_mw)| {
                add_row(p, h, pt, pmw, et, e_mw);
            });

        // Add the parcel layer also
        let parcel_level = linear_interpolate_sounding(snd, parcel.pressure)?;
        let parcel_height = parcel_level.height.ok_or(AnalysisError::MissingValue)?;
        let env_t = parcel_level.temperature.ok_or(AnalysisError::MissingValue)?;
        let env_dp = parcel_level.dew_point.ok_or(AnalysisError::MissingValue)?;
        let env_mw = metfor::mixing_ratio(env_dp, parcel.pressure)?;
        add_row(
            parcel.pressure,
            parcel_height,
            parcel.temperature,
            parcel.mixing_ratio()?,
            env_t,
            env_mw,
        );
    }

    Ok(ParcelProfile {
        pressure,
        height,
        parcel_t,
        parcel_mw,
        environment_t,
        environment_mw,
        parcel,
    })
}

/// Convective available potential energy of a parcel in J/kg
pub fn cape(profile: &ParcelProfile) -> Result<f64> {
    let cape = cape_layers(profile)
        .iter()
        .fold(0.0, |acc, (bottom_p, top_p)| {
            let pressure = &profile.pressure;
            let height = &profile.height;
            let parcel_t = &profile.parcel_t;
            let parcel_mw = &profile.parcel_mw;
            let env_t = &profile.environment_t;
            let env_mw = &profile.environment_mw;

            let virtual_t = |t_k, mw| t_k * (1.0 + mw / metfor::epsilon) / (1.0 + mw);

            let layer_cape = izip!(pressure, height, parcel_t, parcel_mw, env_t, env_mw)
                .skip_while(|(p, _h, _pt, _pmw, _et, _emw)| *p > bottom_p)
                .take_while(|(p, _h, _pt, _pmw, _et, _emw)| *p >= top_p)
                .fold(
                    (0.0, 1_000_000.0, 0.0, 0.0),
                    |acc, (_, h, pt, pmw, et, emw)| {
                        let (mut cape, prev_h, prev_pt, prev_et) = acc;

                        let dz = h - prev_h;
                        let pt = metfor::celsius_to_kelvin(*pt)
                            .and_then(|t_k| Ok(virtual_t(t_k, pmw)))
                            .unwrap_or(0.0);
                        let et = metfor::celsius_to_kelvin(*et)
                            .and_then(|t_k| Ok(virtual_t(t_k, emw)))
                            .unwrap_or(0.0);

                        if pt < 1.0 || et < 1.0 {
                            // Skip over this layer, it has invalid data.
                            (cape, prev_h, prev_pt, prev_et)
                        } else if dz <= 0.0 {
                            // Must be just starting out, save the previous layer and move on
                            (cape, *h, pt, et)
                        } else {
                            cape += ((pt - et) / et + (prev_pt - prev_et) / prev_et) * dz;
                            (cape, *h, pt, et)
                        }
                    },
                )
                .0;

            acc + layer_cape
        });

    Ok(cape / 2.0 * 9.81)
}

// TODO: cin, el, lfc, lcl, dcape, ncape, hail zone cape
