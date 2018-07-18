//! Create and analyze a profile from lifting or descending a parcel.

use metfor;
use sounding_base::{Profile::*, Sounding};

use error::*;
use interpolation::{linear_interp, linear_interpolate_sounding};
use keys::ParcelIndex;
use parcel::Parcel;

/// Hold profiles for a parcel and it's environment.
#[derive(Debug, Clone)]
pub struct ParcelProfile {
    /// Pressure profile
    pub pressure: Vec<f64>,
    /// Height profile
    pub height: Vec<f64>,
    /// Parcel virtual temperature profile
    pub parcel_t: Vec<f64>,
    /// Environment virtual temperature profile
    pub environment_t: Vec<f64>,
}

impl ParcelProfile {}

/// Parcel analysis, this is a way to package the analysis of a parcel.
#[derive(Debug, Clone)]
pub struct ParcelAnalysis {
    // The orginal parcel and profile
    parcel: Parcel,
    profile: ParcelProfile,

    // Indicies from analysis
    cape: Option<f64>,
    hail_cape: Option<f64>,
    ncape: Option<f64>,
    lcl_height_agl: Option<f64>,  // cloud base for aviation
    lcl_pressure: Option<f64>,    // plotting on skew-t
    lcl_temperature: Option<f64>, // ice or ice/water cloud?
    cin: Option<f64>,
    el_pressure: Option<f64>,    // plotting on skew-t
    el_height_asl: Option<f64>,  // Calculating convective cloud tops for aviation
    el_temperature: Option<f64>, // useful for comparing to satellite
    lfc_pressure: Option<f64>,   // plotting on skew-t
}

impl ParcelAnalysis {
    /// Method to retrieve value from analysis.
    #[inline]
    pub fn get_index(&self, var: ParcelIndex) -> Option<f64> {
        use self::ParcelIndex::*;

        match var {
            LCLPressure => self.lcl_pressure,
            LCLHeightAGL => self.lcl_height_agl,
            LCLTemperature => self.lcl_temperature,
            CAPE => self.cape,
            CAPEHail => self.hail_cape,
            NCAPE => self.ncape,
            CIN => self.cin,
            ELPressure => self.el_pressure,
            ELHeightASL => self.el_height_asl,
            ELTemperature => self.el_temperature,
            LFC => self.lfc_pressure,
        }
    }
}

/// Lift a parcel
pub fn lift_parcel(parcel: Parcel, snd: &Sounding) -> Result<ParcelAnalysis> {
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
    let lcl_height = lcl_env.height.ok_or(AnalysisError::InterpolationError)?;
    let lcl_env_temperature = lcl_env
        .temperature
        .ok_or(AnalysisError::InterpolationError)?;
    let lcl_env_dp = lcl_env.dew_point.ok_or(AnalysisError::InterpolationError)?;

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
            metfor::temperature_c_from_theta(theta, tgt_pres).and_then(|t_c| {
                metfor::virtual_temperature_c(
                    t_c,
                    metfor::dew_point_from_p_and_mw(tgt_pres, dry_mw)?,
                    tgt_pres,
                )
            })
        } else {
            // Moist adiabatic lifting
            metfor::temperature_c_from_theta_e_saturated_and_pressure(tgt_pres, theta_e)
                .and_then(|t_c| metfor::virtual_temperature_c(t_c, t_c, tgt_pres))
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

    let mut lfc_pressure: Option<f64> = None;
    let mut el_pressure: Option<f64> = None;

    //
    // Allocate some buffers to hold the return values.
    //
    let mut pressure = Vec::with_capacity(snd_pressure.len() + 5);
    let mut height = Vec::with_capacity(snd_pressure.len() + 5);
    let mut parcel_t = Vec::with_capacity(snd_pressure.len() + 5);
    let mut environment_t = Vec::with_capacity(snd_pressure.len() + 5);

    // Nested scope to limit closure borrows
    {
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
        let mut pcl_t0 = parcel.virtual_temperature_c()?;
        let mut env_t0 = parcel_start_data
            .dew_point
            .ok_or(AnalysisError::InvalidInput)
            .and_then(|dp| {
                Ok(metfor::virtual_temperature_c(
                    parcel_start_data
                        .temperature
                        .ok_or(AnalysisError::InterpolationError)?,
                    dp,
                    p0,
                )?)
            })?;

        add_row(p0, h0, pcl_t0, env_t0);

        if pcl_t0 < env_t0 {
            el_pressure = Some(p0);
        } else {
            lfc_pressure = Some(p0);
        }

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
            // Calculate the environment virtual temperature
            .filter_map(|(p, h, env_t, env_dp, pcl_t)|{
                match metfor::virtual_temperature_c(env_t, env_dp, p) {
                    Ok(env_vt) => Some((p, h, env_vt, pcl_t)),
                    Err(_) => None,
                }
            });

        //
        // Pack the resulting values into their vectors.
        //
        for (p, h, env_t, pcl_t) in iter {
            // Check to see if we are passing the lcl
            if p0 > lcl_pressure && p < lcl_pressure {
                add_row(
                    lcl_pressure,
                    lcl_height,
                    metfor::virtual_temperature_c(lcl_temperature, lcl_temperature, lcl_pressure)?,
                    metfor::virtual_temperature_c(lcl_env_temperature, lcl_env_dp, lcl_pressure)?,
                );
            }

            // Check to see if the parcel and environment soundings have crossed
            if (pcl_t0 < env_t0 && pcl_t > env_t) || (pcl_t0 > env_t0 && pcl_t < env_t) {
                let tgt_pres = linear_interp(0.0, pcl_t - env_t, pcl_t0 - env_t0, p, p0);
                let tgt_row = linear_interpolate_sounding(snd, tgt_pres)?;
                let h2 = tgt_row.height.ok_or(AnalysisError::InterpolationError)?;
                let env_t2 = tgt_row
                    .temperature
                    .ok_or(AnalysisError::InterpolationError)
                    .and_then(|t_c| {
                        let env_dp2 = tgt_row.dew_point.ok_or(AnalysisError::InterpolationError)?;
                        metfor::virtual_temperature_c(t_c, env_dp2, tgt_pres)
                            .map_err(|err| AnalysisError::MetForError(err))
                    })?;

                add_row(tgt_pres, h2, env_t2, env_t2);

                if pcl_t0 < env_t0 && pcl_t > env_t {
                    // LFC crossing into positive bouyancy
                    lfc_pressure = Some(tgt_pres);
                } else {
                    // EL crossing into negative bouyancy
                    el_pressure = Some(tgt_pres);
                }
            }

            // Add the new values to the array
            add_row(p, h, pcl_t, env_t);

            // Remember them for the next iteration
            p0 = p;
            pcl_t0 = pcl_t;
            env_t0 = env_t;
        }
    }

    let profile = ParcelProfile {
        pressure,
        height,
        parcel_t,
        environment_t,
    };

    let lcl_pressure = Some(lcl_pressure);
    let lcl_temperature = Some(lcl_temperature);
    let lcl_height_agl = snd.get_station_info()
        .elevation()
        .into_option()
        .map(|elev| lcl_height - elev);

    // Check lfc and el for consistency
    match (lfc_pressure, el_pressure) {
        (Some(lfcp), Some(elp)) => {
            if lfcp < elp {
                lfc_pressure = None;
                el_pressure = None;
            }
        }
        (None, Some(_)) | (Some(_), None) => {
            lfc_pressure = None;
            el_pressure = None;
        }
        _ => {}
    }

    let (el_height_asl, el_temperature) = if let Some(elp) = el_pressure {
        let level = linear_interpolate_sounding(snd, elp);
        (
            level.ok().and_then(|lvl| lvl.height.into()),
            level.ok().and_then(|lvl| lvl.temperature.into()),
        )
    } else {
        (None, None)
    };

    let lfc_height_asl = if let Some(lfc) = lfc_pressure {
        let level = linear_interpolate_sounding(snd, lfc);
        level.ok().and_then(|lvl| lvl.height.into())
    } else {
        None
    };

    let (cape, cin, hail_cape) = match cape_cin(&profile, lfc_pressure, el_pressure) {
        Ok((cape, cin, hail_cape)) => (Some(cape), Some(cin), Some(hail_cape)),
        Err(_) => (None, None, None),
    };

    let ncape = cape.and_then(|cape| {
        lfc_height_asl
            .and_then(|lfc_h: f64| el_height_asl.and_then(|el_h| Some(cape / (el_h - lfc_h))))
    });

    Ok(ParcelAnalysis {
        parcel,
        profile,
        cape,
        hail_cape,
        ncape,
        lcl_height_agl,
        lcl_pressure,
        lcl_temperature,
        cin,
        el_pressure,
        el_height_asl,
        el_temperature,
        lfc_pressure,
    })
}

/// Descend a parcel dry adiabatically.
pub fn descend_dry(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let theta = parcel.theta()?;
    descend_parcel(parcel, snd, theta, ::metfor::temperature_c_from_theta)
}

/// Descend a parcel moist adiabatically
pub fn descend_moist(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let theta = parcel.theta_e()?;

    let theta_func = |theta_e, press| {
        ::metfor::temperature_c_from_theta_e_saturated_and_pressure(press, theta_e)
    };

    descend_parcel(parcel, snd, theta, theta_func)
}

#[inline]
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
            // Remove levels with missing data
            .filter_map(|(p_opt, h_opt, e_t_opt)| {
                if p_opt.is_some() && h_opt.is_some() && e_t_opt.is_some() {
                    Some((p_opt.unpack(), h_opt.unpack(), e_t_opt.unpack()))
                } else {
                    None
                }
            })
            // Get the parcel temperature
            .filter_map(|(p, h, et)|{
                match theta_func(theta, p) {
                    Ok(pcl_t) => Some((p, h, pcl_t, et)),
                    Err(_) => None,
                }
            })
            .for_each(|(p, h, pt, et)| {
                add_row(p, h, pt, et);
            });

        // Add the parcel layer also
        let parcel_level = linear_interpolate_sounding(snd, parcel.pressure)?;
        let parcel_height = parcel_level.height.ok_or(AnalysisError::MissingValue)?;
        let env_t = parcel_level.temperature.ok_or(AnalysisError::MissingValue)?;
        add_row(parcel.pressure, parcel_height, parcel.temperature, env_t);
    }

    Ok(ParcelProfile {
        pressure,
        height,
        parcel_t,
        environment_t,
    })
}

/// Convective available potential energy of a parcel in J/kg
///
/// Assumes the profile has virtual temperatures in it.
fn cape_cin(profile: &ParcelProfile, lfc: Option<f64>, el: Option<f64>) -> Result<(f64, f64, f64)> {
    let (lfc, el) = if let (Some(lfc), Some(el)) = (lfc, el) {
        (lfc, el)
    } else {
        return Err(AnalysisError::MissingValue);
    };

    let pressure = &profile.pressure;
    let height = &profile.height;
    let parcel_t = &profile.parcel_t;
    let env_t = &profile.environment_t;

    let (cape, cin, hail_zone_cape) = izip!(pressure, height, parcel_t, env_t)
        .take_while(|(&p, _h, _pt, _et)| p >= el)
        .fold(
            ((0.0, 0.0, 0.0), 1_000_000.0, 0.0, 0.0),
            |acc, (&p, &h, &pt, &et)| {
                let ((mut cape, mut cin, mut hail_cape), prev_h, prev_pt, prev_et) = acc;

                let dz = h - prev_h;

                if dz <= 0.0 {
                    // Must be just starting out, save the previous layer and move on
                    ((cape, cin, hail_cape), h, pt, et)
                } else {
                    let bouyancy = ((pt - et) / et + (prev_pt - prev_et) / prev_et) * dz;
                    if pt >= et && prev_pt >= prev_et && p <= lfc {
                        cape += bouyancy;
                        if pt <= -10.0 && pt >= -30.0 {
                            hail_cape += bouyancy;
                        }
                    } else if pt <= et && prev_pt <= prev_et {
                        cin += bouyancy
                    }
                    ((cape, cin, hail_cape), h, pt, et)
                }
            },
        )
        .0;

    Ok((
        cape / 2.0 * 9.81,
        cin / 2.0 * 9.81,
        hail_zone_cape / 2.0 * 9.81,
    ))
}

// TODO: dcape
