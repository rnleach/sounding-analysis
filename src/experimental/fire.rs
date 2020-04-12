//! Experimental sounding analysis for fire weather.

use crate::{
    error::{AnalysisError, Result},
    interpolation::linear_interp,
    parcel::{mixed_layer_parcel, Parcel},
    parcel_profile::{
        find_parcel_start_data,
        lift::{
            create_level_type_mapping, create_parcel_calc_t, parcel_lcl, AnalLevel, AnalLevelType,
        },
        ParcelAscentAnalysis, ParcelProfile,
    },
    sounding::Sounding,
};
use itertools::{izip, Itertools};
use metfor::{self, Celsius, CelsiusDiff, JpKg, Kelvin, Meters, Quantity};

/// Result of lifting a parcel represntative of a plume core.
#[derive(Debug, Clone, Copy)]
pub struct PlumeAscentAnalysis {
    /// The parcel this analysis was completed for.
    pub parcel: Parcel,
    /// Net CAPE up to equilibrium level.
    pub net_cape: Option<JpKg>,
    /// The lifting condensation level of the parcel.
    pub lcl_height: Option<Meters>,
    /// The last EL reached before max_height
    pub el_height: Option<Meters>,
    /// The level where net CAPE becomes zero, the plume rises no more
    pub max_height: Option<Meters>,
}

/// This characterizes how much heating it would take to cause a plume to "blow up".
///
/// The blow up ΔT is the difference in temperature from the parcel that blows up to the parcel the
/// analysis started with (usually the mixed layer parcel). The blow up height is the difference in
/// the equilibrium level of the blow up ΔT plus 0.5C and the blow up ΔT minus 0.5C. This
/// difference was chosen because using the (numerical) derivative is extremely sensitive to the
/// stepsize used when calculating the derivative. This makes sense, since there is actually a step
/// discontinuity at the blow up temperature.
///
/// Two kinds of blow up are analyzed, the blow up of the plume top and the of the equilibrium
/// level. The plume top is defined as the level where the net CAPE becomes 0 again, implying that
/// that all the potential energy (CAPE) that was converted to kinetic energy (updraft) has been
/// used up due to negative bouyancy. At this level, parcels will start descending again.
///
/// The ΔT required to get the plume top over the LCL, and thus to create a cloud is also
/// calculated. This can be a useful for determining if a plume will have a cap cloud but will
/// not be unstable enough to blow up.
#[derive(Debug)]
pub struct BlowUpAnalysis {
    /// The original parcel we started with while searching for the blow up.
    pub starting_parcel: Parcel,
    /// The amount of warming required to cause a blow up of the equilibrium level.
    pub delta_t_el: CelsiusDiff,
    /// The amount of warming required to cause a blow up of the plume top.
    pub delta_t_top: CelsiusDiff,
    /// The amount of warming required to cause a cloud to form.
    pub delta_t_cloud: CelsiusDiff,
    /// The change in height from the blow up of the equilibrium level.
    pub delta_z_el: Meters,
    /// The change in height from the blow up of the plume top.
    pub delta_z_top: Meters,
}

/// Generate a series of `PlumeAscentAnalysis`s for a given sounding starting at the mixed layer
/// parcel and warming it by `increment` until it has gone `max_range` degrees.
///
/// Arguments:
/// * snd - the environmental sounding.
/// * max_range - the maximum amount of heating to apply.
/// * increment - the difference in temperature between parcels.
/// * moisture_ratio - a value of 10 means for every 10C of heating (dt), add 1 g/kg of moisture
///   to the parcel. If it is `None`, don't add any moisture.
///
pub fn calc_plumes(
    snd: &Sounding,
    increment: CelsiusDiff,
    max_range: CelsiusDiff,
    moisture_ratio: Option<f64>,
) -> Result<Vec<PlumeAscentAnalysis>> {
    let (_starting_parcel, parcels_iter) =
        plume_parcels(snd, max_range, increment, moisture_ratio)?;

    Ok(parcels_iter
        .filter_map(|(_, pcl)| analyze_plume_parcel(pcl, snd).ok())
        .skip_while(|anal| anal.el_height.is_none())
        .take_while(|anal| anal.el_height.is_some())
        .collect())
}

/// Find the parcel that causes the plume to blow up by finding the maximum derivative of the
/// equilibrium level vs parcel temperature.
///
/// Arguments:
/// * snd - the environmental sounding.
/// * moisture_ratio - a value of 10 means for every 10C of heating (dt), add 1 g/kg of moisture
///   to the parcel. If it is `None`, don't add any moisture.
///
pub fn blow_up(snd: &Sounding, moisture_ratio: Option<f64>) -> Result<BlowUpAnalysis> {
    const INCREMENT: CelsiusDiff = CelsiusDiff(0.1);
    const MAX_RANGE: CelsiusDiff = CelsiusDiff(20.0);

    let (starting_parcel, parcel_iter) = plume_parcels(snd, MAX_RANGE, INCREMENT, moisture_ratio)?;

    let (delta_t_cloud, delta_t_el, delta_t_top, _, _, delta_z_el, delta_z_top): (
        CelsiusDiff,
        CelsiusDiff,
        CelsiusDiff,
        f64,
        f64,
        Meters,
        Meters,
    ) = parcel_iter
        // Do the analysis, ignore errors.
        .filter_map(|(dt, pcl)| analyze_plume_parcel(pcl, snd).ok().map(|anal| (dt, anal)))
        // Skip levels with no useful information.
        .skip_while(|(_dt, anal)| {
            anal.el_height.is_none() && anal.max_height.is_none() && anal.lcl_height.is_none()
        })
        // Take while there is some useful information.
        .take_while(|(_dt, anal)| {
            anal.el_height.is_some() || anal.max_height.is_some() || anal.lcl_height.is_some()
        })
        // Filter out noise due to rounding errors etc. where warmer parcels don't rise as high.
        // This smooths it out so that deriviatives work well too.
        .scan(
            (Meters(0.0), Meters(0.0), Meters(0.0)),
            |(prev_lcl, prev_el, prev_max_z), (dt, anal)| {
                if let Some(lcl) = anal.lcl_height {
                    if lcl >= *prev_lcl {
                        *prev_lcl = lcl;
                    } else {
                        return Some(None);
                    }
                }

                if let Some(el) = anal.el_height {
                    if el >= *prev_el {
                        *prev_el = el;
                    } else {
                        return Some(None);
                    }
                }

                if let Some(max_z) = anal.max_height {
                    if max_z >= *prev_max_z {
                        *prev_max_z = max_z;
                    } else {
                        return Some(None);
                    }
                }

                Some(Some((dt, anal)))
            },
        )
        // Filter out "bad" layers.
        .filter_map(|opt| opt)
        // Pair up for simple derivative calculations.
        .tuple_windows::<(_, _)>()
        .fold(
            (
                CelsiusDiff(0.0),
                CelsiusDiff(0.0),
                CelsiusDiff(0.0),
                0.0,
                0.0,
                Meters(0.0),
                Meters(0.0),
            ),
            |acc, (lvl0, lvl1)| {
                let (
                    mut cloud_dt,
                    mut el_blow_up_dt,
                    mut max_z_blow_up_dt,
                    mut deriv_el,
                    mut deriv_max_z,
                    mut jump_el,
                    mut jump_max_z,
                ) = acc;

                let (dt0, anal0) = lvl0;
                let (dt1, anal1) = lvl1;

                debug_assert_ne!(dt0, dt1); // Required for division below
                let dx = (dt1 - dt0).unpack();

                if let (Some(el0), Some(el1)) = (anal0.el_height, anal1.el_height) {
                    let derivative = (el1 - el0).unpack() / dx;
                    if derivative > deriv_el {
                        el_blow_up_dt = CelsiusDiff((dt1 + dt0).unpack() / 2.0);
                        deriv_el = derivative;
                        jump_el = el1 - el0;
                    }
                }

                if let (Some(max_z_0), Some(max_z_1)) = (anal0.max_height, anal1.max_height) {
                    let derivative = (max_z_1 - max_z_0).unpack() / dx;
                    if derivative > deriv_max_z {
                        max_z_blow_up_dt = CelsiusDiff((dt1 + dt0).unpack() / 2.0);
                        deriv_max_z = derivative;
                        jump_max_z = max_z_1 - max_z_0;
                    }
                }

                if let (None, Some(_)) = (anal0.lcl_height, anal1.lcl_height) {
                    if cloud_dt == CelsiusDiff(0.0) {
                        cloud_dt = CelsiusDiff((dt1 + dt0).unpack() / 2.0);
                    }
                }

                (
                    cloud_dt,
                    el_blow_up_dt,
                    max_z_blow_up_dt,
                    deriv_el,
                    deriv_max_z,
                    jump_el,
                    jump_max_z,
                )
            },
        );

    Ok(BlowUpAnalysis {
        starting_parcel,
        delta_t_cloud,
        delta_t_el,
        delta_z_el,
        delta_t_top,
        delta_z_top,
    })
}

/// Lift a parcel until the net CAPE is zero.
pub fn analyze_plume_parcel(parcel: Parcel, snd: &Sounding) -> Result<PlumeAscentAnalysis> {
    // Get the starting parcel and the iterator to lift it.
    let (parcel, lift_iter) = lift_parcel(parcel, snd)?;

    let mut max_height: Option<Meters> = None;

    // Construct an iterator that selects the environment values and calculates the
    // corresponding parcel values.
    let (lcl_height, el_height, net_bouyancy) = lift_iter
        // Scan to get the max_height
        .scan(
            (0.0, Meters(0.0)),
            |(prev_bouyancy, prev_height), (int_bouyancy, anal_level_type)| {
                use crate::parcel_profile::lift::AnalLevelType::*;

                let height_val = match anal_level_type {
                    Normal(level) | LFC(level) | LCL(level) | EL(level) => level.height,
                };

                if *prev_bouyancy >= 0.0 && int_bouyancy < 0.0 {
                    let mx_height =
                        linear_interp(0.0, *prev_bouyancy, int_bouyancy, *prev_height, height_val);
                    max_height = Some(mx_height);
                }

                *prev_bouyancy = int_bouyancy;
                *prev_height = height_val;

                Some((int_bouyancy, anal_level_type))
            },
        )
        // Fold to get the EL Level, Max Height, and max integrated bouyancy
        .fold(
            (None, None, 0.0f64),
            |acc, (int_bouyancy, anal_level_type)| {
                use crate::parcel_profile::lift::AnalLevelType::*;

                let (mut lcl, mut el, mut max_bouyancy) = acc;

                max_bouyancy = max_bouyancy.max(int_bouyancy);
                match anal_level_type {
                    Normal(_) | LFC(_) => {}
                    LCL(level) => {
                        lcl = Some(level.height);
                    }
                    EL(level) => {
                        el = Some(level.height);
                    }
                }
                (lcl, el, max_bouyancy)
            },
        );

    let net_cape = if el_height.is_some() {
        Some(JpKg(net_bouyancy / 2.0 * -metfor::g))
    } else {
        None
    };

    // Make sure LCL is below max height
    let lcl_height = if let (Some(mxh), Some(lcl)) = (max_height, lcl_height) {
        if mxh > lcl {
            Some(lcl)
        } else {
            None
        }
    } else {
        lcl_height
    };

    Ok(PlumeAscentAnalysis {
        lcl_height,
        el_height,
        max_height,
        net_cape,
        parcel,
    })
}

/// Lift a parcel until the net CAPE is zero.
pub fn lift_plume_parcel(
    parcel: Parcel,
    snd: &Sounding,
) -> Result<(ParcelProfile, PlumeAscentAnalysis)> {
    // Get the starting parcel and the iterator to lift it.
    let (parcel, lift_iter) = lift_parcel(parcel, snd)?;

    let len = snd.pressure_profile().len();
    let mut pressure = Vec::with_capacity(len);
    let mut height = Vec::with_capacity(len);
    let mut parcel_t = Vec::with_capacity(len);
    let mut environment_t = Vec::with_capacity(len);

    let mut max_height: Option<Meters> = None;

    // Construct an iterator that selects the environment values and calculates the
    // corresponding parcel values.
    let (lcl_height, el_height, net_bouyancy) = lift_iter
        // Add the levels to the parcel profile.
        .scan((), |_dummy, (int_bouyancy, anal_level_type)| {
            use crate::parcel_profile::lift::AnalLevelType::*;
            match anal_level_type {
                Normal(lvl) | LFC(lvl) | LCL(lvl) | EL(lvl) => {
                    let AnalLevel {
                        pressure: p,
                        height: h,
                        pcl_virt_t,
                        env_virt_t,
                    } = lvl;
                    pressure.push(p);
                    height.push(h);
                    parcel_t.push(pcl_virt_t);
                    environment_t.push(env_virt_t);
                }
            }

            Some((int_bouyancy, anal_level_type))
        })
        // FIXME: STOP HERE, pass this iterator to analyze_plume_parcel to get the plume ascent
        // analysis. This will reduce the duplication of code!!!!
        // Scan to get the max_height
        .scan(
            (0.0, Meters(0.0)),
            |(prev_bouyancy, prev_height), (int_bouyancy, anal_level_type)| {
                use crate::parcel_profile::lift::AnalLevelType::*;

                let height_val = match anal_level_type {
                    Normal(level) | LFC(level) | LCL(level) | EL(level) => level.height,
                };

                if *prev_bouyancy >= 0.0 && int_bouyancy < 0.0 {
                    let mx_height =
                        linear_interp(0.0, *prev_bouyancy, int_bouyancy, *prev_height, height_val);
                    max_height = Some(mx_height);
                }

                *prev_bouyancy = int_bouyancy;
                *prev_height = height_val;

                Some((int_bouyancy, anal_level_type))
            },
        )
        // Fold to get the EL Level, Max Height, and max integrated bouyancy
        .fold(
            (None, None, 0.0f64),
            |acc, (int_bouyancy, anal_level_type)| {
                use crate::parcel_profile::lift::AnalLevelType::*;

                let (mut lcl, mut el, mut max_bouyancy) = acc;

                max_bouyancy = max_bouyancy.max(int_bouyancy);
                match anal_level_type {
                    Normal(_) | LFC(_) => {}
                    LCL(level) => {
                        lcl = Some(level.height);
                    }
                    EL(level) => {
                        el = Some(level.height);
                    }
                }
                (lcl, el, max_bouyancy)
            },
        );

    let net_cape = if el_height.is_some() {
        Some(JpKg(net_bouyancy / 2.0 * -metfor::g))
    } else {
        None
    };

    // Make sure LCL is below max height
    let lcl_height = if let (Some(mxh), Some(lcl)) = (max_height, lcl_height) {
        if mxh > lcl {
            Some(lcl)
        } else {
            None
        }
    } else {
        lcl_height
    };

    Ok((
        ParcelProfile {
            pressure,
            height,
            parcel_t,
            environment_t,
        },
        PlumeAscentAnalysis {
            lcl_height,
            el_height,
            max_height,
            net_cape,
            parcel,
        },
    ))
}

/// Get the starting parcel and build an iterator to lift it.
fn lift_parcel<'a>(
    parcel: Parcel,
    snd: &'a Sounding,
) -> Result<(Parcel, impl Iterator<Item = (f64, AnalLevelType)> + 'a)> {
    // Find the LCL
    let (pcl_lcl, _lcl_temperature) = parcel_lcl(&parcel, snd)?;

    // The starting level to lift the parcel from
    let (_parcel_start_data, parcel) = find_parcel_start_data(snd, &parcel)?;

    // How to calculate a parcel temperature for a given pressure level
    let parcel_calc_t = create_parcel_calc_t(parcel, pcl_lcl)?;
    let level_type_mapping = create_level_type_mapping(pcl_lcl);

    // Get the environment data to iterate over. We want the parcel profile to have all the same
    // pressure levels as the environmental sounding, plus a few special ones.
    let snd_pressure = snd.pressure_profile();
    let hgt = snd.height_profile();
    let env_t = snd.temperature_profile();
    let env_dp = snd.dew_point_profile();

    let p0 = parcel.pressure;

    // Construct an iterator that selects the environment values and calculates the
    // corresponding parcel values.
    let iter = izip!(snd_pressure, hgt, env_t, env_dp)
        // Remove rows with missing data
        .filter(|(p, h, t, dp)| p.is_some() && h.is_some() && t.is_some() && dp.is_some())
        // Unpack from the `Optioned` type
        .map(|(p, h, t, dp)| (p.unpack(), h.unpack(), t.unpack(), dp.unpack()))
        // Remove rows at or below the parcel level
        .filter(move |(p, _, _, _)| *p <= p0)
        // Calculate the parcel temperature, skip this level if there is an error
        .filter_map(move |(p, h, env_t, env_dp)| {
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
        .flat_map(move |(lvl0, lvl1)| level_type_mapping(lvl0, lvl1))
        // Pair the levels up to integrate the bouyancy.
        .tuple_windows::<(_, _)>()
        // Integrate the bouyancy.
        .scan(
            (0.0, 0.0, 0usize),
            |(prev_int_bouyancy, int_bouyancy, count), (anal_level_type0, anal_level_type1)| {
                use crate::parcel_profile::lift::AnalLevelType::*;

                let level_data0: &AnalLevel = match &anal_level_type0 {
                    Normal(data) | LFC(data) | LCL(data) | EL(data) => data,
                };
                let level_data1: &AnalLevel = match &anal_level_type1 {
                    Normal(data) | LFC(data) | LCL(data) | EL(data) => data,
                };

                let &AnalLevel {
                    height: h0,
                    pcl_virt_t: pcl0,
                    env_virt_t: env0,
                    ..
                } = level_data0;
                let &AnalLevel {
                    height: h1,
                    pcl_virt_t: pcl1,
                    env_virt_t: env1,
                    ..
                } = level_data1;

                let Meters(dz) = h1 - h0;
                debug_assert!(dz >= 0.0);

                let b0 = (pcl0 - env0).unpack() / Kelvin::from(env0).unpack();
                let b1 = (pcl1 - env1).unpack() / Kelvin::from(env1).unpack();
                let bouyancy = (b0 + b1) * dz;
                *prev_int_bouyancy = *int_bouyancy;
                *int_bouyancy += bouyancy;
                *count += 1;
                Some((
                    (*prev_int_bouyancy, *int_bouyancy, *count),
                    anal_level_type1,
                ))
            },
        )
        // Take until the bouyancy goes negative, then we're done, just run through the first five
        // to ensure we get through any goofy surface layers. This is also needed to get started if
        // the parcel temperature is the same as the sounding surface temperature. Also, for the
        // case of a fire plume, near the surface things are chaotic, so we should punch through
        // a shallow surface stable layer.
        //
        // Use the prev_int_bouyancy to at least one point past where the bouyancy becomes zero
        // so we have enough data to interpolate.
        .take_while(|((prev_int_bouyancy, _, count), _)| *prev_int_bouyancy >= 0.0 || *count < 5)
        .map(|((_, int_bouyancy, _), anal_level)| (int_bouyancy, anal_level));

    Ok((parcel, iter))
}

/// Given a sounding, return an iterator that creates parcels starting with the mixed layer parcel
/// and then incrementing the parcel temperature up to `plus_range` in increments of `increment`.
///
/// If the surface temperature is warmer than the mixed layer parcel temperature, then the starting
/// parcel will use the surface temperature with the mixed layer moisture.
///
/// Arguments:
/// * snd - the environmental sounding.
/// * plus_range - the maximum amount of heating to apply.
/// * increment - the difference in temperature between parcels.
/// * moisture_ratio - a value of 10 means for every 10C of heating (dt), add 1 g/kg of moisture
///   to the parcel. If it is `None`, don't add any moisture.
///
fn plume_parcels(
    snd: &Sounding,
    plus_range: CelsiusDiff,
    increment: CelsiusDiff,
    moisture_ratio: Option<f64>,
) -> Result<(Parcel, impl Iterator<Item = (CelsiusDiff, Parcel)>)> {
    let parcel = mixed_layer_parcel(snd)?;
    let (row, mut parcel) = find_parcel_start_data(snd, &parcel)?;
    if row.temperature.unwrap() > parcel.temperature {
        parcel.temperature = row.temperature.unwrap();
    }

    let CelsiusDiff(inc_val) = increment;
    let next_dt = CelsiusDiff(-inc_val);

    Ok((
        parcel,
        PlumeParcelIterator {
            starting_pcl: parcel,
            next_dt,
            max_dt: plus_range,
            increment,
            moisture_ratio,
        },
    ))
}

/// Iterator for `plume_parcels` function that generates increasingly warmer parcels with constant
/// moisture.
struct PlumeParcelIterator {
    starting_pcl: Parcel,
    next_dt: CelsiusDiff,
    max_dt: CelsiusDiff,
    increment: CelsiusDiff,
    moisture_ratio: Option<f64>,
}

impl Iterator for PlumeParcelIterator {
    type Item = (CelsiusDiff, Parcel);

    fn next(&mut self) -> Option<Self::Item> {
        self.next_dt += self.increment;
        if self.next_dt > self.max_dt {
            None
        } else {
            Some((
                self.next_dt,
                create_plume_parcel_from(self.starting_pcl, self.next_dt, self.moisture_ratio),
            ))
        }
    }
}

/// Create a new parcel assuming a starting parcel and a temperature increment.
///
/// Arguments:
/// * environment_parcel - the original starting parcel
/// * dt - the amount of heating for this parcel
/// * moisture_ratio - a value of 10 means for every 10C of heating (dt), add 1 g/kg of moisture
///   to the parcel. If it is `None`, don't add any moisture.
///
pub fn create_plume_parcel_from(
    environment_parcel: Parcel,
    dt: CelsiusDiff,
    moisture_ratio: Option<f64>,
) -> Parcel {
    let next_t = environment_parcel.temperature + dt;

    let next_td = if let Some(ratio) = moisture_ratio {
        let mw =
            metfor::specific_humidity(environment_parcel.dew_point, environment_parcel.pressure)
                .expect("error creating specific humidity.");

        let mw = mw + dt.unpack() / (1_000.0 * ratio);
        metfor::dew_point_from_p_and_mw(environment_parcel.pressure, mw)
            .expect("error creating dew point.")
    } else {
        environment_parcel.dew_point
    };

    Parcel {
        temperature: next_t,
        dew_point: next_td,
        ..environment_parcel
    }
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
pub fn partition_cape(pa: &ParcelAscentAnalysis) -> Result<(JpKg, JpKg)> {
    let lcl = pa.lcl_pressure().ok_or(AnalysisError::MissingValue)?;
    let el = pa.el_pressure().ok_or(AnalysisError::MissingValue)?;

    let parcel_theta = pa.parcel().theta();

    let profile = pa.profile();

    let lower_dry_profile = izip!(
        &profile.pressure,
        &profile.height,
        &profile.parcel_t,
        &profile.environment_t
    )
    .take_while(|(p, _, _, _)| **p >= lcl)
    .map(|(_, h, pt, et)| (*h, Kelvin::from(*pt), Kelvin::from(*et)));

    let upper_dry_profile = izip!(&profile.pressure, &profile.height, &profile.environment_t)
        .skip_while(|(p, _, _)| **p >= lcl)
        .filter_map(|(p, h, et)| {
            let t_k = metfor::temperature_from_theta(parcel_theta, *p);
            metfor::virtual_temperature(t_k, t_k, *p).map(|pt_k| (*p, *h, pt_k, *et))
        })
        .take_while(|(_, _, pt, et)| pt >= et)
        .map(|(_, h, pt, et)| (h, pt, Kelvin::from(et)));

    let dry_profile = lower_dry_profile.chain(upper_dry_profile);

    let full_profile = izip!(
        &profile.pressure,
        &profile.height,
        &profile.parcel_t,
        &profile.environment_t
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
