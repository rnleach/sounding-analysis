//! Create and analyze a profile from lifting or descending a parcel.

use metfor;
use sounding_base::{DataRow, Profile::*, Sounding};

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
    lifted_index: Option<f64>,
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
            LI => self.lifted_index,
        }
    }

    /// Retrieve the parcel's profile
    #[inline]
    pub fn get_profile(&self) -> &ParcelProfile {
        &self.profile
    }

    /// Retrieve the original parcel.
    #[inline]
    pub fn get_parcel(&self) -> &Parcel {
        &self.parcel
    }

    /// Calculate the parcel velocity at the equilibrium level. Note that this is most likely an
    /// over estimate due to the effects of entrainment and water/ice loading.
    ///
    /// Units are meters per second.
    #[inline]
    pub fn calculate_cape_velocity(&self) -> Option<f64> {
        self.cape.map(|cape| f64::sqrt(2.0 * cape))
    }
}

/// Lift a parcel for a convective parcel analysis.
///
/// The resulting `ParcelProfile` and analysis are based off of virtual temperatures and the idea
/// that if there is no *moist* convection, or convective cloud, then there is no CAPE or CIN.
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
    let (parcel_start_data, parcel) = find_parcel_start_data(snd, &parcel)?;

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

    //
    // Initialize some special levels/values we'll want to find during lifting
    //
    let mut lfc_pressure: Option<f64> = None;
    let mut el_pressure: Option<f64> = None;
    let mut lifted_index: Option<f64> = None;

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
        let mut h0 = parcel_start_data
            .height
            .ok_or(AnalysisError::InvalidInput)?;
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
        // Construct an iterator that selects the environment values and calculates the
        // corresponding parcel values.
        //
        let iter = izip!(snd_pressure, hgt, env_t, env_dp)
            // Remove rows with missing data and unpack options
            .filter_map(|(p, h, env_t, env_dp)| {
                if p.is_some() && h.is_some() && env_t.is_some() && env_dp.is_some() {
                    Some((p.unpack(), h.unpack(), env_t.unpack(), env_dp.unpack()))
                } else {
                    None
                }
            })
            // Remove rows at or below the parcel level
            .filter(move |(p, _, _, _)| *p < p0)
            // Calculate the parcel temperature, skip this level if there is an error
            .filter_map(|(p, h, env_t, env_dp)| match calc_parcel_t(p) {
                Ok(pcl_t) => Some((p, h, env_t, env_dp, pcl_t)),
                Err(_) => None,
            })
            // Calculate the environment virtual temperature, skip levels with errors
            .filter_map(|(p, h, env_t, env_dp, pcl_t)| {
                match metfor::virtual_temperature_c(env_t, env_dp, p) {
                    Ok(env_vt) => Some((p, h, env_vt, pcl_t)),
                    Err(_) => None,
                }
            });

        //
        // Pack the resulting values into their vectors and handle special levels
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
                let h2 = linear_interp(tgt_pres, p0, p, h0, h);
                let env_t2 = linear_interp(tgt_pres, p0, p, env_t0, env_t);

                add_row(tgt_pres, h2, env_t2, env_t2);

                if pcl_t0 < env_t0 && pcl_t > env_t {
                    // LFC crossing into positive bouyancy
                    if el_pressure.is_none() || el_pressure.unwrap() > lcl_pressure {
                        lfc_pressure = Some(tgt_pres);
                        el_pressure = None;
                    }
                } else {
                    // EL crossing into negative bouyancy
                    if el_pressure.is_none(){
                        el_pressure = Some(tgt_pres);
                    }
                }
            }

            // Check to see if the parcel is passing 500 hPa for the LI calculation
            if p0 >= 500.0 && p <= 500.0 {
                let tgt_et = linear_interp(500.0, p0, p, env_t0, env_t);
                let tgt_pt = linear_interp(500.0, p0, p, pcl_t0, pcl_t);
                lifted_index = Some(tgt_et - tgt_pt);
            }

            // Add the new values to the array
            add_row(p, h, pcl_t, env_t);

            // Remember them for the next iteration
            p0 = p;
            h0 = h;
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
    let lcl_height_agl = snd
        .get_station_info()
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
        _ => unreachable!(),
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

    let (cape, cin, hail_cape) = match cape_cin(&profile, lcl_pressure, lfc_pressure, el_pressure) {
        Ok((cape, cin, hail_cape)) => (Some(cape), Some(cin), Some(hail_cape)),
        Err(_) => (None, None, None),
    };

    let ncape = cape.and_then(|cape| {
        lfc_height_asl.and_then(|lfc_h: f64| el_height_asl.map(|el_h| cape / (el_h - lfc_h)))
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
        lifted_index,
    })
}

// In order for parcel lifting to work and create a parallel environmental profile, we need to
// start at a level in the sounding with pressure, height, temperature, and dew point. Otherwise
// we end up with too much missing data in the sounding.
fn find_parcel_start_data(snd: &Sounding, parcel: &Parcel) -> Result<(DataRow, Parcel)> {
    let good_row = |row: &DataRow| -> bool {
        row.temperature.is_some()
            && row.dew_point.is_some()
            && row.pressure.is_some()
            && row.height.is_some()
    };

    let first_guess = linear_interpolate_sounding(snd, parcel.pressure)?;
    if good_row(&first_guess) {
        return Ok((first_guess, *parcel));
    }

    let second_guess = snd
        .bottom_up()
        .filter(good_row)
        .nth(0)
        .ok_or(AnalysisError::NotEnoughData)?;

    let pressure = second_guess.pressure.ok_or(AnalysisError::InvalidInput)?;
    let theta = parcel.theta()?;
    let temperature = metfor::temperature_c_from_theta(theta, pressure)?;
    let mw = parcel.mixing_ratio()?;
    let dew_point = metfor::dew_point_from_p_and_mw(pressure, mw)?;
    let new_parcel = Parcel {
        pressure,
        temperature,
        dew_point,
    };

    Ok((second_guess, new_parcel))
}

/// Descend a parcel dry adiabatically.
///
/// The resulting `ParcelProfile` has actual temperatures and not virtual temperatures. This is for
/// analyzing inversions and visualizing what a sounding would look like if deep, dry mixing were
/// to occur from surface heating alone.
pub fn mix_down(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let theta = parcel.theta()?;
    descend_parcel(
        parcel,
        snd,
        theta,
        ::metfor::temperature_c_from_theta,
        false,
        false,
    )
}

/// Descend a parcel moist adiabatically.
///
/// The resulting `ParcelProfile` has virtual temperatures and is intended for calculating
/// DCAPE.
fn descend_moist(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let theta = parcel.theta_e()?;

    let theta_func = |theta_e, press| {
        ::metfor::temperature_c_from_theta_e_saturated_and_pressure(press, theta_e)
    };

    descend_parcel(parcel, snd, theta, theta_func, true, true)
}

#[inline]
fn descend_parcel<F>(
    parcel: Parcel,
    snd: &Sounding,
    theta: f64,
    theta_func: F,
    saturated: bool,
    virtual_t: bool,
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
    let env_dp = snd.get_profile(DewPoint);
    let hght = snd.get_profile(GeopotentialHeight);

    let pcl_mw = parcel.mixing_ratio()?;

    // Nested scope to limit borrows
    {
        // Helper function to add row to parcel profile
        let mut add_row = |pp, hh, pcl_tt, env_tt| {
            pressure.push(pp);
            height.push(hh);
            parcel_t.push(pcl_tt);
            environment_t.push(env_tt);
        };

        izip!(press, hght, env_t, env_dp)
            .take_while(|(p_opt, _, _, _)| {
                if p_opt.is_some() {
                    p_opt.unpack() >= parcel.pressure
                } else {
                    true // Just skip over levels with missing data
                }
            })
            // Remove levels with missing data
            .filter_map(|(p_opt, h_opt, e_t_opt, e_dp_opt)| {
                if p_opt.is_some() && h_opt.is_some() && e_t_opt.is_some() && e_dp_opt.is_some() {
                    Some((
                        p_opt.unpack(),
                        h_opt.unpack(),
                        e_t_opt.unpack(),
                        e_dp_opt.unpack(),
                    ))
                } else {
                    None
                }
            })
            // Get the parcel temperature
            .filter_map(|(p, h, e_t, e_dp)| match theta_func(theta, p) {
                Ok(pcl_t) => Some((p, h, pcl_t, e_t, e_dp)),
                Err(_) => None,
            })
            // Get the parcel dew point
            .filter_map(|(p, h, pcl_t, e_t, e_dp)| {
                let p_dp = if saturated {
                    pcl_t
                } else {
                    match metfor::dew_point_from_p_and_mw(p, pcl_mw) {
                        Ok(dp) => dp,
                        Err(_) => return None,
                    }
                };

                Some((p, h, pcl_t, p_dp, e_t, e_dp))
            })
            // Convert to virtual temperature if needed.
            .filter_map(|(p, h, pcl_t, p_dp, e_t, e_dp)| {
                let pcl_t = if virtual_t {
                    match metfor::virtual_temperature_c(pcl_t, p_dp, p) {
                        Ok(vt) => vt,
                        Err(_) => return None,
                    }
                } else {
                    pcl_t
                };

                let e_t = if virtual_t {
                    match metfor::virtual_temperature_c(e_t, e_dp, p) {
                        Ok(vt) => vt,
                        Err(_) => return None,
                    }
                } else {
                    e_t
                };

                Some((p, h, pcl_t, e_t))
            }).for_each(|(p, h, pt, et)| {
                add_row(p, h, pt, et);
            });

        // Add the parcel layer also
        let parcel_level = linear_interpolate_sounding(snd, parcel.pressure)?;
        let parcel_height = parcel_level.height.ok_or(AnalysisError::MissingValue)?;
        let env_t = parcel_level
            .temperature
            .ok_or(AnalysisError::MissingValue)?;
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
/// Assumes the profile has virtual temperatures in it. Returns a tuple with the values
/// (CAPE, CIN, hail zone cape)
fn cape_cin(
    profile: &ParcelProfile,
    lcl: Option<f64>,
    lfc: Option<f64>,
    el: Option<f64>,
) -> Result<(f64, f64, f64)> {
    let (lfc, el) = if let (Some(lcl), Some(lfc), Some(el)) = (lcl, lfc, el) {
        // If no LCL, then no moist convection, then don't mention CAPE/CIN
        if el < lcl {
            (lfc, el)
        } else {
            // No cloud, no moist convection
            return Ok((0.0, 0.0, 0.0));
        }
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

                let (pt, et) = match (metfor::celsius_to_kelvin(pt), metfor::celsius_to_kelvin(et))
                {
                    (Ok(pt), Ok(et)) => (pt, et),
                    _ => return ((cape, cin, hail_cape), prev_h, prev_pt, prev_et),
                };

                let dz = h - prev_h;

                if dz <= 0.0 {
                    // Must be just starting out, save the previous layer and move on
                    ((cape, cin, hail_cape), h, pt, et)
                } else {
                    let bouyancy = ((pt - et) / et + (prev_pt - prev_et) / prev_et) * dz;
                    if bouyancy > 0.0 && p <= lfc {
                        cape += bouyancy;
                        if pt <= -10.0 && pt >= -30.0 {
                            hail_cape += bouyancy;
                        }
                    } else if bouyancy < 0.0 {
                        cin += bouyancy
                    }
                    ((cape, cin, hail_cape), h, pt, et)
                }
            },
        ).0;

    Ok((
        cape / 2.0 * 9.81,
        cin / 2.0 * 9.81,
        hail_zone_cape / 2.0 * 9.81,
    ))
}

/// Downdraft CAPE.
///
/// Defined as the net area between a parcel descended moist adiabatically from the level of the
/// lowest theta-e in the lowest 400 hPa of the sounding.
///
/// Returns the profile, downdraft cape, and downrush temperature in a tuple.
pub fn dcape(snd: &Sounding) -> Result<(ParcelProfile, f64, f64)> {
    let t = snd.get_profile(Temperature);
    let dp = snd.get_profile(DewPoint);
    let p = snd.get_profile(Pressure);

    // Find the lowest pressure, 400 mb above the surface (or starting level)
    let top_p = -400.0 + p
        .iter()
        .filter_map(|p| p.into_option())
        .nth(0)
        .ok_or(AnalysisError::NotEnoughData)?;

    // Find the starting parcel.
    let pcl = izip!(p, t, dp)
        .filter_map(|(p, t, dp)| {
            if p.is_some() && t.is_some() && dp.is_some() {
                Some((p.unpack(), t.unpack(), dp.unpack()))
            } else {
                None
            }
        }).take_while(|(p, _, _)| *p >= top_p)
        .fold(
            Err(AnalysisError::NotEnoughData),
            |acc: Result<Parcel>, (p, t, dp)| match metfor::theta_e_kelvin(t, dp, p) {
                Ok(theta_e) => match acc {
                    Ok(parcel) => {
                        if let Ok(old_theta) = parcel.theta_e() {
                            if old_theta < theta_e {
                                Ok(parcel)
                            } else {
                                Ok(Parcel {
                                    temperature: t,
                                    dew_point: dp,
                                    pressure: p,
                                })
                            }
                        } else {
                            Ok(Parcel {
                                temperature: t,
                                dew_point: dp,
                                pressure: p,
                            })
                        }
                    }
                    Err(_) => Ok(Parcel {
                        temperature: t,
                        dew_point: dp,
                        pressure: p,
                    }),
                },
                Err(_) => acc,
            },
        )?;

    let profile = descend_moist(pcl, snd)?;

    let mut dcape = 0.0;
    let mut h0 = 1_000_000_000.0; // Big number
    let mut pt0 = 0.0;
    let mut et0 = 0.0;
    for (&p, &h, &pt, &et) in izip!(
        &profile.pressure,
        &profile.height,
        &profile.parcel_t,
        &profile.environment_t
    ) {
        let pt = match metfor::theta_kelvin(p, pt) {
            Ok(pt) => pt,
            Err(_) => continue,
        };

        let et = match metfor::theta_kelvin(p, et) {
            Ok(et) => et,
            Err(_) => continue,
        };

        let dz = h - h0;
        // we must be starting out, becuase h0 starts as a large positive number
        if dz <= 0.0 {
            h0 = h;
            pt0 = pt;
            et0 = et;
            continue;
        }

        dcape += ((pt - et) / et + (pt0 - et0) / et0) * dz;

        h0 = h;
        pt0 = pt;
        et0 = et;
    }

    // - for integration direction should be top down, 9.81 for gravity, and 2.0 for trapezoid rule.
    dcape *= -9.81 / 2.0;

    let downrush_t = *profile.parcel_t.get(0).ok_or(AnalysisError::MissingValue)?;

    Ok((profile, dcape, downrush_t))
}

/// Partition the CAPE between dry and moist ascent contributions. EXPERIMENTAL.
///
/// This is an experimental function that calculates how much CAPE there would be with a "dry"
/// ascent only. Above the LCL it keeps the parcel saturated but keeps lifting it at the dry
/// adiabatic lapse rate, and then calculates the CAPE of this profile. The difference between this
/// value and the CAPE is the amount of CAPE added by latent heat release. It isn't perfect, but
/// when applied to a convective parcel (think CCL and convective temperature) it can be used
/// to partition the energy contributed by heating the column from the sun and the energy added by
/// latent heat release. This can be useful for analyzing convection initiated by wildfire and
/// estimating how much the convective column is being driven by the surface heating and how much it
/// is being driven by latent heat release.
///
/// Returns a tuple with `(dry_cape, wet_cape)`
pub fn partition_cape(pa: &ParcelAnalysis) -> Result<(f64, f64)> {
    let total_cape = pa.cape.ok_or(AnalysisError::MissingValue)?;
    let lcl = pa.lcl_pressure.ok_or(AnalysisError::MissingValue)?;

    let parcel_theta = pa.parcel.theta()?;

    let lower_dry_profile = izip!(
        &pa.profile.pressure,
        &pa.profile.height,
        &pa.profile.parcel_t,
        &pa.profile.environment_t
    ).take_while(|(p, _, _, _)| **p >= lcl)
    .map(|(_, h, pt, et)| (*h, *pt, *et));

    let upper_dry_profile = izip!(
        &pa.profile.pressure,
        &pa.profile.height,
        &pa.profile.environment_t
    ).skip_while(|(p, _, _)| **p >= lcl)
    .filter_map(|(p, h, et)| {
        metfor::temperature_c_from_theta(parcel_theta, *p)
            .and_then(|t| metfor::virtual_temperature_c(t, t, *p))
            .ok()
            .map(|pt| (*p, *h, pt, *et))
    }).take_while(|(_, _, pt, et)| pt >= et)
    .map(|(_, h, pt, et)| (h, pt, et));

    let dry_profile = lower_dry_profile.chain(upper_dry_profile);

    let dry_cape = dry_profile
        .fold((0.0, 1_000_000.0, 0.0, 0.0), |acc, (h, pt, et)| {
            let (mut cape, prev_h, prev_pt, prev_et) = acc;

            let (pt, et) = match (metfor::celsius_to_kelvin(pt), metfor::celsius_to_kelvin(et)) {
                (Ok(pt), Ok(et)) => (pt, et),
                _ => return (cape, prev_h, prev_pt, prev_et),
            };

            let dz = h - prev_h;

            if dz <= 0.0 {
                // Must be just starting out, save the previous layer and move on
                (cape, h, pt, et)
            } else {
                let bouyancy = ((pt - et) / et + (prev_pt - prev_et) / prev_et) * dz;
                if bouyancy > 0.0 {
                    cape += bouyancy;
                }
                (cape, h, pt, et)
            }
        }).0;

    let dry_cape = (dry_cape / 2.0 * 9.81).min(total_cape);
    let wet_cape = total_cape - dry_cape;

    Ok((dry_cape, wet_cape))
}
