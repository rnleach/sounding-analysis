//! Create and analyze a profile from lifting or descending a parcel.

use metfor::{self, Celsius, CelsiusDiff, HectoPascal, JpKg, Kelvin, Meters, MetersPSec, Quantity};
use optional::{none, some, Optioned};
use sounding_base::{DataRow, Sounding};

use crate::error::*;
use crate::interpolation::{linear_interp, linear_interpolate_sounding};
use crate::parcel::Parcel;

/// Hold profiles for a parcel and it's environment.
#[derive(Debug, Clone)]
pub struct ParcelProfile {
    /// Pressure profile
    pub pressure: Vec<HectoPascal>,
    /// Height profile
    pub height: Vec<Meters>,
    /// Parcel virtual temperature profile
    pub parcel_t: Vec<Celsius>,
    /// Environment virtual temperature profile
    pub environment_t: Vec<Celsius>,
}

/// Parcel analysis, this is a way to package the analysis of a parcel.
#[derive(Debug, Clone)]
pub struct ParcelAnalysis {
    // The orginal parcel and profile
    parcel: Parcel,
    profile: ParcelProfile,

    // Indicies from analysis
    cape: Optioned<JpKg>,
    hail_cape: Optioned<JpKg>,
    ncape: Optioned<f64>,
    lcl_height_agl: Optioned<Meters>,    // cloud base for aviation
    lcl_pressure: Optioned<HectoPascal>, // plotting on skew-t
    lcl_temperature: Optioned<Celsius>,  // ice or ice/water cloud?
    cin: Optioned<JpKg>,
    el_pressure: Optioned<HectoPascal>,      // plotting on skew-t
    el_height_asl: Optioned<Meters>,         // Calculating convective cloud tops for aviation
    el_temperature: Optioned<Celsius>,       // useful for comparing to satellite
    lfc_pressure: Optioned<HectoPascal>,     // plotting on skew-t
    lfc_virt_temperature: Optioned<Celsius>, // plotting on skew-t
    lifted_index: Optioned<CelsiusDiff>,
}

impl ParcelAnalysis {
    /// Get the CAPE.
    pub fn cape(&self) -> Optioned<JpKg> {
        self.cape
    }

    /// Get the CAPE in the hail growth zone.
    pub fn hail_cape(&self) -> Optioned<JpKg> {
        self.hail_cape
    }

    /// Get the normalized cape.
    pub fn ncape(&self) -> Optioned<f64> {
        self.ncape
    }

    /// Get the LCL height AGL.
    pub fn lcl_height_agl(&self) -> Optioned<Meters> {
        self.lcl_height_agl
    }

    /// Get the LCL pressrue level.
    pub fn lcl_pressure(&self) -> Optioned<HectoPascal> {
        self.lcl_pressure
    }
    /// Get the temperature at the LCL.
    pub fn lcl_temperature(&self) -> Optioned<Celsius> {
        self.lcl_temperature
    }
    /// Get the CIN.
    pub fn cin(&self) -> Optioned<JpKg> {
        self.cin
    }
    /// Get the pressure at the equilibrium level.
    pub fn el_pressure(&self) -> Optioned<HectoPascal> {
        self.el_pressure
    }
    /// Get the height ASL of the equilibrium level.
    pub fn el_height_asl(&self) -> Optioned<Meters> {
        self.el_height_asl
    }
    /// Get the temperature at the equilibrium level.
    pub fn el_temperature(&self) -> Optioned<Celsius> {
        self.el_temperature
    }
    /// Get the pressure at the LFC.
    pub fn lfc_pressure(&self) -> Optioned<HectoPascal> {
        self.lfc_pressure
    }
    /// Get the virtual temperature at the LFC.
    pub fn lfc_virt_temperature(&self) -> Optioned<Celsius> {
        self.lfc_virt_temperature
    }
    /// Get the lifted index.
    pub fn lifted_index(&self) -> Optioned<CelsiusDiff> {
        self.lifted_index
    }

    /// Retrieve the parcel's profile
    #[inline]
    pub fn profile(&self) -> &ParcelProfile {
        &self.profile
    }

    /// Retrieve the original parcel.
    #[inline]
    pub fn parcel(&self) -> &Parcel {
        &self.parcel
    }

    /// Calculate the parcel vertical speed at the equilibrium level. Note that this is most likely
    /// an over estimate of updraft speed due to the effects of entrainment and water/ice loading.
    #[inline]
    pub fn calculate_cape_speed(&self) -> Option<MetersPSec> {
        self.cape
            .map(|cape| MetersPSec::pack(f64::sqrt(2.0 * cape.unpack())))
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
    )
    .ok_or(AnalysisError::MetForError)?;

    let lcl_temperature = Celsius::from(lcl_temperature);
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
    let theta = parcel.theta();
    let theta_e = parcel.theta_e()?;
    let dry_mw = parcel.mixing_ratio()?;
    let calc_parcel_t = |tgt_pres| {
        if tgt_pres > lcl_pressure {
            // Dry adiabatic lifting
            let t_k = metfor::temperature_from_theta(theta, tgt_pres);
            metfor::virtual_temperature(
                t_k,
                metfor::dew_point_from_p_and_mw(tgt_pres, dry_mw)?,
                tgt_pres,
            )
            .map(Celsius::from)
        } else {
            // Moist adiabatic lifting
            metfor::temperature_from_theta_e_saturated_and_pressure(tgt_pres, theta_e)
                .and_then(|t_c| metfor::virtual_temperature(t_c, t_c, tgt_pres))
                .map(Celsius::from)
        }
    };

    //
    // Get the environment data to iterate over. We want the parcel profile to have all the same
    // pressure levels as the environmental sounding, plus a few special ones.
    //
    let snd_pressure = snd.pressure_profile();
    let hgt = snd.height_profile();
    let env_t = snd.temperature_profile();
    let env_dp = snd.dew_point_profile();

    //
    // Initialize some special levels/values we'll want to find during lifting
    //
    let mut lfc_pressure: Optioned<HectoPascal> = none();
    let mut lfc_virt_temperature: Optioned<Celsius> = none();
    let mut el_pressure: Optioned<HectoPascal> = none();
    let mut lifted_index: Optioned<CelsiusDiff> = none();

    //
    // Allocate some buffers to hold the return values.
    //
    let mut pressure = Vec::with_capacity(snd_pressure.len() + 5);
    let mut height = Vec::with_capacity(snd_pressure.len() + 5);
    let mut parcel_t: Vec<Celsius> = Vec::with_capacity(snd_pressure.len() + 5);
    let mut environment_t: Vec<Celsius> = Vec::with_capacity(snd_pressure.len() + 5);

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
        let mut pcl_t0 = parcel.virtual_temperature().map(Celsius::from)?;
        let mut env_t0 = parcel_start_data
            .dew_point
            .ok_or(AnalysisError::InvalidInput)
            .and_then(|dp| {
                Ok(metfor::virtual_temperature(
                    parcel_start_data
                        .temperature
                        .ok_or(AnalysisError::InterpolationError)?,
                    dp,
                    p0,
                )
                .map(Celsius::from)
                .ok_or(AnalysisError::MetForError)?)
            })?;

        add_row(p0, h0, pcl_t0, env_t0);

        if pcl_t0 < env_t0 {
            el_pressure = some(p0);
        } else {
            lfc_pressure = some(p0);
            lfc_virt_temperature = some(pcl_t0);
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
            .filter_map(|(p, h, env_t, env_dp)| {
                calc_parcel_t(p).map(|pcl_t| (p, h, env_t, env_dp, pcl_t))
            })
            // Calculate the environment virtual temperature, skip levels with errors
            .filter_map(|(p, h, env_t, env_dp, pcl_t)| {
                metfor::virtual_temperature(env_t, env_dp, p)
                    .map(|env_vt| (p, h, Celsius::from(env_vt), pcl_t))
            });

        //
        // Pack the resulting values into their vectors and handle special levels
        //
        for (p, h, env_t, pcl_t) in iter {
            // Check to see if we are passing the lcl
            let lcl_data = if p0 > lcl_pressure && p < lcl_pressure {
                Some((
                    lcl_pressure,
                    lcl_height,
                    metfor::virtual_temperature(lcl_temperature, lcl_temperature, lcl_pressure)
                        .map(Celsius::from)
                        .ok_or(AnalysisError::MetForError)?,
                    metfor::virtual_temperature(lcl_env_temperature, lcl_env_dp, lcl_pressure)
                        .map(Celsius::from)
                        .ok_or(AnalysisError::MetForError)?,
                ))
            } else {
                None
            };

            // Check to see if the parcel and environment soundings have crossed
            let prof_cross_data =
                if (pcl_t0 < env_t0 && pcl_t > env_t) || (pcl_t0 > env_t0 && pcl_t < env_t) {
                    let tgt_pres =
                        linear_interp(CelsiusDiff(0.0), pcl_t - env_t, pcl_t0 - env_t0, p, p0);
                    let h2 = linear_interp(tgt_pres, p0, p, h0, h);
                    let env_t2 = linear_interp(tgt_pres, p0, p, env_t0, env_t);

                    Some((tgt_pres, h2, env_t2, env_t2))
                } else {
                    None
                };

            if let (Some((lclp, lclh, lclpt, lclet)), Some((cp, ch, cpt, cet))) =
                (lcl_data, prof_cross_data)
            {
                // Handle both with proper ordering
                if lclp > cp {
                    add_row(lclp, lclh, lclpt, lclet);
                    add_row(cp, ch, cpt, cet);
                } else {
                    add_row(cp, ch, cpt, cet);
                    add_row(lclp, lclh, lclpt, lclet);
                }
            } else if let Some((p, h, pt, et)) = lcl_data {
                // Just handle the lcl level alone
                add_row(p, h, pt, et);
            } else if let Some((p, h, pt, et)) = prof_cross_data {
                // Just handle adding the crossing level alone
                add_row(p, h, pt, et);
            }

            // In any ordering, deal with LFC and EL for crossings here
            if let Some((tgt_pres, _, vt, _)) = prof_cross_data {
                if pcl_t0 < env_t0 && pcl_t > env_t {
                    // LFC crossing into positive bouyancy
                    if el_pressure.is_none() || el_pressure.unwrap() > lcl_pressure {
                        lfc_pressure = some(tgt_pres);
                        lfc_virt_temperature = some(vt);
                        el_pressure = none();
                    }
                } else if el_pressure.is_none() {
                    // EL crossing into negative bouyancy
                    el_pressure = some(tgt_pres);
                }
            }

            // Check to see if the parcel is passing 500 hPa for the LI calculation
            if p0 >= HectoPascal(500.0) && p <= HectoPascal(500.0) {
                let tgt_et = linear_interp(HectoPascal(500.0), p0, p, env_t0, env_t);
                let tgt_pt = linear_interp(HectoPascal(500.0), p0, p, pcl_t0, pcl_t);
                lifted_index = some(tgt_et - tgt_pt);
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

    let lcl_pressure = some(lcl_pressure);
    let lcl_temperature = some(lcl_temperature);
    let lcl_height_agl: Optioned<Meters> = snd
        .station_info()
        .elevation()
        .into_option()
        .map(|elev| lcl_height - elev)
        .into();

    // Check lfc and el for consistency
    if lfc_pressure.is_some() && el_pressure.is_some() {
        if lfc_pressure < el_pressure {
            lfc_pressure = none();
            lfc_virt_temperature = none();
            el_pressure = none();
        }
    } else {
        lfc_pressure = none();
        lfc_virt_temperature = none();
        el_pressure = none();
    }

    let (el_height_asl, el_temperature): (Optioned<Meters>, Optioned<Celsius>) =
        if let Some(elp) = el_pressure.into_option() {
            let level = linear_interpolate_sounding(snd, elp);
            (
                level.ok().and_then(|lvl| lvl.height.into()).into(),
                level.ok().and_then(|lvl| lvl.temperature.into()).into(),
            )
        } else {
            (none(), none())
        };

    let lfc_height_asl: Optioned<Meters> = if let Some(lfc) = lfc_pressure.into_option() {
        let level = linear_interpolate_sounding(snd, lfc);
        Optioned::from(level.ok().and_then(|lvl| lvl.height.into()))
    } else {
        none()
    };

    let (cape, cin, hail_cape) = match cape_cin(&profile, lcl_pressure, lfc_pressure, el_pressure) {
        Ok((cape, cin, hail_cape)) => (some(cape), some(cin), some(hail_cape)),
        Err(_) => (none(), none(), none()),
    };

    let ncape = if let (Some(cape), Some(lfc_h), Some(el_h)) = (
        cape.into_option(),
        lfc_height_asl.into_option(),
        el_height_asl.into_option(),
    ) {
        some(cape.unpack() / (el_h - lfc_h).unpack())
    } else {
        none()
    };

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
        lfc_virt_temperature,
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
    let theta = parcel.theta();
    let temperature = Celsius::from(metfor::temperature_from_theta(theta, pressure));
    let mw = parcel.mixing_ratio()?;
    let dew_point =
        metfor::dew_point_from_p_and_mw(pressure, mw).ok_or(AnalysisError::MetForError)?;
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
    let theta = parcel.theta();
    let theta_func = |theta_val, press| {
        Some(Celsius::from(metfor::temperature_from_theta(
            theta_val, press,
        )))
    };

    descend_parcel(parcel, snd, theta, theta_func, false, false)
}

/// Descend a parcel moist adiabatically.
///
/// The resulting `ParcelProfile` has virtual temperatures and is intended for calculating
/// DCAPE.
fn descend_moist(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let theta = parcel.theta_e()?;

    let theta_func =
        |theta_e, press| metfor::temperature_from_theta_e_saturated_and_pressure(press, theta_e);

    descend_parcel(parcel, snd, theta, theta_func, true, true)
}

#[inline]
fn descend_parcel<F>(
    parcel: Parcel,
    snd: &Sounding,
    theta: Kelvin,
    theta_func: F,
    saturated: bool,
    virtual_t: bool,
) -> Result<ParcelProfile>
where
    F: Fn(Kelvin, HectoPascal) -> Option<Celsius>,
{
    let mut pressure = Vec::new();
    let mut height = Vec::new();
    let mut parcel_t = Vec::new();
    let mut environment_t = Vec::new();

    // Actually start at the bottom and work up.
    let press = snd.pressure_profile();
    let env_t = snd.temperature_profile();
    let env_dp = snd.dew_point_profile();
    let hght = snd.height_profile();

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
            .filter_map(|(p, h, e_t, e_dp)| {
                theta_func(theta, p).map(|pcl_t| (p, h, pcl_t, e_t, e_dp))
            })
            // Get the parcel dew point
            .filter_map(|(p, h, pcl_t, e_t, e_dp)| {
                let p_dp = if saturated {
                    pcl_t
                } else {
                    metfor::dew_point_from_p_and_mw(p, pcl_mw)?
                };

                Some((p, h, pcl_t, p_dp, e_t, e_dp))
            })
            // Convert to virtual temperature if needed.
            .filter_map(|(p, h, pcl_t, p_dp, e_t, e_dp)| {
                let pcl_t = if virtual_t {
                    Celsius::from(metfor::virtual_temperature(pcl_t, p_dp, p)?)
                } else {
                    pcl_t
                };

                let e_t = if virtual_t {
                    Celsius::from(metfor::virtual_temperature(e_t, e_dp, p)?)
                } else {
                    e_t
                };

                Some((p, h, pcl_t, e_t))
            })
            .for_each(|(p, h, pt, et)| {
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
    lcl: Optioned<HectoPascal>,
    lfc: Optioned<HectoPascal>,
    el: Optioned<HectoPascal>,
) -> Result<(JpKg, JpKg, JpKg)> {
    let (lfc, el) = if let (Some(lcl), Some(lfc), Some(el)) =
        (lcl.into_option(), lfc.into_option(), el.into_option())
    {
        // If no LCL, then no moist convection, then don't mention CAPE/CIN
        if el < lcl {
            (lfc, el)
        } else {
            // No cloud, no moist convection
            return Ok((JpKg(0.0), JpKg(0.0), JpKg(0.0)));
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
            (
                (0.0, 0.0, 0.0),
                Meters(std::f64::MAX),
                Kelvin(0.0),
                Kelvin(0.0),
            ),
            |acc, (&p, &h, &pt, &et)| {
                let ((mut cape, mut cin, mut hail_cape), prev_h, prev_pt, prev_et) = acc;

                let (pt, et) = (Kelvin::from(pt), Kelvin::from(et));

                let dz = h - prev_h;

                if dz <= Meters(0.0) {
                    // Must be just starting out, save the previous layer and move on
                    ((cape, cin, hail_cape), h, pt, et)
                } else {
                    let bouyancy = ((pt - et).unpack() / et.unpack()
                        + (prev_pt - prev_et).unpack() / prev_et.unpack())
                        * dz.unpack();
                    if bouyancy > 0.0 && p <= lfc {
                        cape += bouyancy;
                        if pt <= Celsius(-10.0) && pt >= Celsius(-30.0) {
                            hail_cape += bouyancy;
                        }
                    } else if bouyancy < 0.0 {
                        cin += bouyancy
                    }
                    ((cape, cin, hail_cape), h, pt, et)
                }
            },
        )
        .0;

    Ok((
        JpKg(cape / 2.0 * -metfor::g),
        JpKg(cin / 2.0 * -metfor::g),
        JpKg(hail_zone_cape / 2.0 * -metfor::g),
    ))
}

/// Downdraft CAPE.
///
/// Defined as the net area between a parcel descended moist adiabatically from the level of the
/// lowest theta-e in the lowest 400 hPa of the sounding.
///
/// Returns the profile, downdraft cape, and downrush temperature in a tuple.
pub fn dcape(snd: &Sounding) -> Result<(ParcelProfile, JpKg, Celsius)> {
    let t = snd.temperature_profile();
    let dp = snd.dew_point_profile();
    let p = snd.pressure_profile();

    // Find the lowest pressure, 400 mb above the surface (or starting level)
    let top_p = p
        .iter()
        .filter_map(|p| p.into_option())
        .nth(0)
        .ok_or(AnalysisError::NotEnoughData)?
        - HectoPascal(400.0);

    // Check that pressure is positive. This is theoretically possible if the surface pressure is
    // less than 400 hPa. Maybe in the Himalayas? Moisture is insignficant at 100 hPa anyway. This
    // is a really gross error check.
    debug_assert!(
        top_p.into_option().is_some() && top_p > HectoPascal(100.0),
        "surface pressure is below 500 mb"
    );

    // Find the starting parcel.
    let pcl = izip!(p, t, dp)
        .filter_map(|(p, t, dp)| {
            if p.is_some() && t.is_some() && dp.is_some() {
                Some((p.unpack(), t.unpack(), dp.unpack()))
            } else {
                None
            }
        })
        .take_while(|(p, _, _)| *p >= top_p)
        .fold(
            Err(AnalysisError::NotEnoughData),
            |acc: Result<Parcel>, (p, t, dp)| match metfor::theta_e(t, dp, p) {
                Some(th_e) => match acc {
                    Ok(parcel) => {
                        if let Ok(old_theta) = parcel.theta_e() {
                            if old_theta < th_e {
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
                None => acc,
            },
        )?;

    let profile = descend_moist(pcl, snd)?;

    let mut dcape = 0.0;
    let mut h0 = Meters(std::f64::MAX); // Big number
    let mut pt0 = Kelvin(0.0);
    let mut et0 = Kelvin(0.0);
    for (&p, &h, &pt, &et) in izip!(
        &profile.pressure,
        &profile.height,
        &profile.parcel_t,
        &profile.environment_t
    ) {
        let pt = metfor::theta(p, pt);
        let et = metfor::theta(p, et);

        let dz = h - h0;
        // we must be starting out, becuase h0 starts as a large positive number
        if dz <= Meters(0.0) {
            h0 = h;
            pt0 = pt;
            et0 = et;
            continue;
        }

        dcape +=
            ((pt - et).unpack() / et.unpack() + (pt0 - et0).unpack() / et0.unpack()) * dz.unpack();

        h0 = h;
        pt0 = pt;
        et0 = et;
    }

    // - for integration direction should be top down, 9.81 for gravity, and 2.0 for trapezoid rule.
    dcape *= metfor::g / 2.0;

    let downrush_t = *profile.parcel_t.get(0).ok_or(AnalysisError::MissingValue)?;

    Ok((profile, JpKg(dcape), downrush_t))
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
pub fn partition_cape(pa: &ParcelAnalysis) -> Result<(JpKg, JpKg)> {
    let lcl = pa.lcl_pressure.ok_or(AnalysisError::MissingValue)?;
    let el = pa.el_pressure.ok_or(AnalysisError::MissingValue)?;

    let parcel_theta = pa.parcel.theta();

    let lower_dry_profile = izip!(
        &pa.profile.pressure,
        &pa.profile.height,
        &pa.profile.parcel_t,
        &pa.profile.environment_t
    )
    .take_while(|(p, _, _, _)| **p >= lcl)
    .map(|(_, h, pt, et)| (*h, Kelvin::from(*pt), Kelvin::from(*et)));

    let upper_dry_profile = izip!(
        &pa.profile.pressure,
        &pa.profile.height,
        &pa.profile.environment_t
    )
    .skip_while(|(p, _, _)| **p >= lcl)
    .filter_map(|(p, h, et)| {
        let t_k = metfor::temperature_from_theta(parcel_theta, *p);
        metfor::virtual_temperature(t_k, t_k, *p).map(|pt_k| (*p, *h, pt_k, *et))
    })
    .take_while(|(_, _, pt, et)| pt >= et)
    .map(|(_, h, pt, et)| (h, pt, Kelvin::from(et)));

    let dry_profile = lower_dry_profile.chain(upper_dry_profile);

    let full_profile = izip!(
        &pa.profile.pressure,
        &pa.profile.height,
        &pa.profile.parcel_t,
        &pa.profile.environment_t
    )
    .take_while(|(p, _, _, _)| **p >= el)
    .map(|(_, h, pt, et)| (*h, Kelvin::from(*pt), Kelvin::from(*et)));

    fn calc_cape<T: Iterator<Item = (Meters, Kelvin, Kelvin)>>(iter: T) -> f64 {
        let cape = iter
            .fold(
                (0.0, Meters(std::f64::MAX), Kelvin(0.0), Kelvin(0.0)),
                |acc, (h, pt, et)| {
                    let (mut cape, prev_h, prev_pt, prev_et) = acc;

                    let dz = h - prev_h;

                    if dz <= Meters(0.0) {
                        // Must be just starting out, save the previous layer and move on
                        (cape, h, pt, et)
                    } else {
                        let bouyancy = ((pt - et).unpack() / et.unpack()
                            + (prev_pt - prev_et).unpack() / prev_et.unpack())
                            * dz.unpack();
                        cape += bouyancy;

                        (cape, h, pt, et)
                    }
                },
            )
            .0;

        cape / 2.0 * -metfor::g
    }

    let total_cape = calc_cape(full_profile).max(0.0);
    let dry_cape = calc_cape(dry_profile).max(0.0).min(total_cape);

    let wet_cape = total_cape - dry_cape;

    Ok((JpKg(dry_cape), JpKg(wet_cape)))
}
