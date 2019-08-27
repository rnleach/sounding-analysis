use super::{find_parcel_start_data, ParcelAscentAnalysis, ParcelProfile};
use crate::{
    error::{AnalysisError, Result},
    interpolation::{linear_interp, linear_interpolate, linear_interpolate_sounding},
    parcel::Parcel,
    sounding::Sounding,
};
use itertools::{izip, Itertools};
use metfor::{self, Celsius, CelsiusDiff, HectoPascal, JpKg, Kelvin, Meters, Quantity};
use optional::{none, some, Optioned};
use std::cmp::Ordering;

pub fn lift_parcel_new(parcel: Parcel, snd: &Sounding) -> Result<ParcelAscentAnalysis> {
    // Find the LCL
    let (pcl_lcl, lcl_temperature) = parcel_lcl(&parcel, snd)?;

    // The starting level to lift the parcel from
    let (parcel_start_data, parcel) = find_parcel_start_data(snd, &parcel)?;

    // How to calculate a parcel temperature for a given pressure level
    let parcel_calc_t = create_parcel_calc_t(parcel, pcl_lcl)?;
    let level_type_mapping = create_level_type_mapping(pcl_lcl);

    // Get the environment data to iterate over. We want the parcel profile to have all the same
    // pressure levels as the environmental sounding, plus a few special ones.
    let snd_pressure = snd.pressure_profile();
    let hgt = snd.height_profile();
    let env_t = snd.temperature_profile();
    let env_dp = snd.dew_point_profile();

    // Allocate some buffers to hold the return values.
    let mut pressure: Vec<HectoPascal> = Vec::with_capacity(snd_pressure.len() + 5);
    let mut height: Vec<Meters> = Vec::with_capacity(snd_pressure.len() + 5);
    let mut parcel_t: Vec<Celsius> = Vec::with_capacity(snd_pressure.len() + 5);
    let mut environment_t: Vec<Celsius> = Vec::with_capacity(snd_pressure.len() + 5);

    // Start by adding the parcel level
    let p0 = parcel.pressure;
    let h0 = parcel_start_data
        .height
        .ok_or(AnalysisError::InvalidInput)?;
    let pcl_t0 = parcel.virtual_temperature().map(Celsius::from)?;
    let env_t0 = parcel_start_data
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

    pressure.push(p0);
    height.push(h0);
    parcel_t.push(pcl_t0);
    environment_t.push(env_t0);

    // Construct an iterator that selects the environment values and calculates the
    // corresponding parcel values.
    let (lfc, el): (Option<AnalLevel>, Option<AnalLevel>) = izip!(snd_pressure, hgt, env_t, env_dp)
        // Remove rows with missing data
        .filter(|(p, h, t, dp)| p.is_some() && h.is_some() && t.is_some() && dp.is_some())
        // Unpack from the `Optioned` type
        .map(|(p, h, t, dp)| (p.unpack(), h.unpack(), t.unpack(), dp.unpack()))
        // Remove rows at or below the parcel level
        .filter(move |(p, _, _, _)| *p < p0)
        // Calculate the parcel temperature, skip this level if there is an error
        .filter_map(|(p, h, env_t, env_dp)| {
            parcel_calc_t(p).map(|pcl_virt_t| (p, h, env_t, env_dp, pcl_virt_t))
        })
        // Calculate the environment virtual temperature, skip levels with errors
        .filter_map(|(p, h, env_t, env_dp, pcl_virt_t)| {
            metfor::virtual_temperature(env_t, env_dp, p)
                .map(|env_vt| (p, h, Celsius::from(env_vt), pcl_virt_t))
        })
        // Wrap in the AnalLevel type
        .map(|(pressure, height, env_virt_t, pcl_virt_t)| AnalLevel {
            pressure,
            height,
            pcl_virt_t,
            env_virt_t,
        })
        // Look at them two levels at a time to check for crossing any special levels
        .tuple_windows::<(_, _)>()
        // Find the level type and insert special levels if needed.
        .flat_map(|(lvl0, lvl1)| level_type_mapping(lvl0, lvl1))
        // Add every level to the vectors.
        .scan(true, |start_flag, anal_level_type| {
            use AnalLevelType::*;

            let level_data: &AnalLevel = match &anal_level_type {
                Normal(data) | LFC(data) | LCL(data) | EL(data) => data,
            };

            pressure.push(level_data.pressure);
            height.push(level_data.height);
            parcel_t.push(level_data.pcl_virt_t);
            environment_t.push(level_data.env_virt_t);

            if *start_flag && level_data.pcl_virt_t >= level_data.env_virt_t {
                *start_flag = false;
                Some(LFC(*level_data))
            } else {
                Some(anal_level_type)
            }
        })
        // Analyze the levels to find the LCL, EL, LFC
        .fold((None, None), |acc, anal_level_type| {
            use AnalLevelType::*;

            let (mut lfc, mut el) = acc;

            match anal_level_type {
                Normal(_) => {}
                LFC(level_data) => {
                    if el.is_some() {
                        el = None;
                    }
                    lfc = Some(level_data);
                }
                EL(level_data) => {
                    if lfc.is_some() {
                        el = Some(level_data);
                    }
                }
                LCL(_) => {}
            };

            (lfc, el)
        });

    // Wrap the vectors into the ParcelProfile
    let profile = ParcelProfile {
        pressure,
        height,
        parcel_t,
        environment_t,
    };

    // Finalize the LCL variables.
    let lcl_pressure = some(pcl_lcl.pressure);
    let lcl_temperature = some(lcl_temperature);
    let lcl_height_agl: Optioned<Meters> = snd
        .station_info()
        .elevation()
        .into_option()
        .map(|elev| pcl_lcl.height - elev)
        .into();

    // Finalize the LFC and EL levels.
    let (lfc_pressure, lfc_virt_temperature, lfc_height_asl) = lfc
        .map(|lfc_level| {
            (
                some(lfc_level.pressure),
                some(lfc_level.env_virt_t),
                some(lfc_level.height),
            )
        })
        .unwrap_or((none(), none(), none()));

    let (el_pressure, el_height_asl) = el
        .map(|el_level| (some(el_level.pressure), some(el_level.height)))
        .unwrap_or((none(), none()));

    let el_temperature: Optioned<Celsius> =
        el_pressure.and_then(|pres| linear_interpolate(snd_pressure, env_t, pres));

    // Get the cape/cin values.
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

    Ok(ParcelAscentAnalysis {
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
    })
}

// A level in the analysis
#[derive(Clone, Copy)]
struct AnalLevel {
    pressure: HectoPascal,
    height: Meters,
    pcl_virt_t: Celsius,
    env_virt_t: Celsius,
}

enum AnalLevelType {
    Normal(AnalLevel),
    LFC(AnalLevel),
    LCL(AnalLevel),
    EL(AnalLevel),
}

struct AnalLevelTypeIterator {
    vals: [Option<AnalLevelType>; 4],
    next: usize,
}

fn parcel_lcl(parcel: &Parcel, snd: &Sounding) -> Result<(AnalLevel, Celsius)> {
    let (pressure, temperature) = metfor::pressure_and_temperature_at_lcl(
        parcel.temperature,
        parcel.dew_point,
        parcel.pressure,
    )
    .ok_or(AnalysisError::MetForError)?;

    let temperature = Celsius::from(temperature);
    let lcl_env = linear_interpolate_sounding(snd, pressure)?;
    let height = lcl_env.height.ok_or(AnalysisError::InterpolationError)?;
    let lcl_env_temperature = lcl_env
        .temperature
        .ok_or(AnalysisError::InterpolationError)?;
    let lcl_env_dp = lcl_env.dew_point.ok_or(AnalysisError::InterpolationError)?;
    let env_virt_t = Celsius::from(
        metfor::virtual_temperature(lcl_env_temperature, lcl_env_dp, pressure)
            .ok_or(AnalysisError::MetForError)?,
    );
    let pcl_virt_t = Celsius::from(
        metfor::virtual_temperature(temperature, temperature, pressure)
            .ok_or(AnalysisError::MetForError)?,
    );

    Ok((
        AnalLevel {
            pressure,
            height,
            pcl_virt_t,
            env_virt_t,
        },
        temperature,
    ))
}

fn create_parcel_calc_t(
    parcel: Parcel,
    lcl: AnalLevel,
) -> Result<impl Fn(HectoPascal) -> Option<Celsius>> {
    let theta = parcel.theta();
    let theta_e = parcel.theta_e()?;
    let dry_mw = parcel.mixing_ratio()?;

    Ok(move |tgt_pres| {
        if tgt_pres > lcl.pressure {
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
    })
}

fn create_level_type_mapping(
    lcl_info: AnalLevel,
) -> impl Fn(AnalLevel, AnalLevel) -> AnalLevelTypeIterator {
    move |lvl0: AnalLevel, lvl1: AnalLevel| -> AnalLevelTypeIterator {
        let mut iter = AnalLevelTypeIterator::default();
        let mut next_idx = 0usize;

        iter.vals[next_idx] = Some(AnalLevelType::Normal(lvl0));
        next_idx += 1;

        let AnalLevel {
            pcl_virt_t: pt0,
            env_virt_t: et0,
            pressure: p0,
            height: h0,
        } = lvl0;

        let AnalLevel {
            pcl_virt_t: pt1,
            env_virt_t: et1,
            pressure: p1,
            height: h1,
        } = lvl1;

        // Check to see if the parcel profile crossed over the environmental profile. Note that
        // this demarks a change in stability, either to stable or unstable.
        if (pt0 <= et0 && pt1 >= et1) || (pt0 >= et0 && pt1 <= et1) {
            let tgt_p = linear_interp(CelsiusDiff(0.0), pt0 - et0, pt1 - et1, p0, p1);
            let tgt_t = linear_interp(CelsiusDiff(0.0), pt0 - et0, pt1 - et1, pt0, pt1);
            let tgt_h = linear_interp(CelsiusDiff(0.0), pt0 - et0, pt1 - et1, h0, h1);

            let tgt_level = AnalLevel {
                pressure: tgt_p,
                height: tgt_h,
                pcl_virt_t: tgt_t,
                env_virt_t: tgt_t,
            };

            let tgt_level_type = if pt0 <= et0 && pt1 >= et1 {
                AnalLevelType::LFC(tgt_level)
            } else {
                AnalLevelType::EL(tgt_level)
            };

            iter.vals[next_idx] = Some(tgt_level_type);
            next_idx += 1;
        }

        // Check for the LCL, add it
        let AnalLevel {
            pressure: lcl_p, ..
        } = lcl_info;
        if p0 <= lcl_p && p1 >= lcl_p {
            iter.vals[next_idx] = Some(AnalLevelType::LCL(lcl_info));
        }

        // Sort the vals array in decreasing order by pressure
        iter.vals.sort_by(|a, b| {
            use AnalLevelType::*;

            let pa = match a {
                Some(Normal(p)) | Some(LCL(p)) | Some(LFC(p)) | Some(EL(p)) => p.pressure,
                None => HectoPascal(0.0),
            };

            let pb = match b {
                Some(Normal(p)) | Some(LCL(p)) | Some(LFC(p)) | Some(EL(p)) => p.pressure,
                None => HectoPascal(0.0),
            };

            // swap order of b and a to get decreasing sort.
            pb.partial_cmp(&pa).unwrap_or(Ordering::Equal)
        });

        iter
    }
}

impl Iterator for AnalLevelTypeIterator {
    type Item = AnalLevelType;

    fn next(&mut self) -> Option<Self::Item> {
        let item = if self.next >= 4 {
            None
        } else {
            self.vals[self.next].take()
        };

        self.next += 1;
        item
    }
}

impl Default for AnalLevelTypeIterator {
    fn default() -> Self {
        AnalLevelTypeIterator {
            vals: [None, None, None, None],
            next: 0,
        }
    }
}

#[deprecated]
#[allow(dead_code)]
// FIXME: delete when no longer needed.
pub fn lift_parcel_orig(parcel: Parcel, snd: &Sounding) -> Result<ParcelAscentAnalysis> {
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

    Ok(ParcelAscentAnalysis {
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
