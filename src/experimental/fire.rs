//! Experimental sounding analysis for fire weather.

use crate::{
    error::{AnalysisError, Result},
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
use metfor::{Celsius, CelsiusDiff, JpKg, Kelvin, Meters, Quantity};

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

/// This characterizes how much heating it would take to cause a plume to "blow up" as defined by
/// the sudden change in the parcel equilibrium level with respect to parcel heating.
///
/// The blow up ΔT is the difference in temperature from the parcel that blows up to the parcel the
/// analysis started with (usually the mixed layer parcel). The blow up height is the difference in
/// the equilibrium level of the blow up ΔT plus 0.5C and the blow up ΔT minus 0.5C. This
/// difference was chosen because using the (numerical) derivative is extremely sensitive to the
/// stepsize used when calculating the derivative. This makes sense, since there is actually a step
/// discontinuity at the blow up temperature.
#[derive(Debug)]
pub struct BlowUpAnalysis {
    /// The original parcel we started with while searching for the blow up.
    pub starting_parcel: Parcel,
    /// The amount of warming required to cause a blow up.
    pub delta_t: CelsiusDiff,
    /// The change in height from the blow up.
    pub height: Meters,
}

/// Generate a series of `PlumeAscentAnalysis`s for a given sounding starting at the mixed layer
/// parcel and warming it by `increment` until it has gone `max_range` degrees.
pub fn calc_plumes(
    snd: &Sounding,
    increment: CelsiusDiff,
    max_range: CelsiusDiff,
) -> Result<Vec<PlumeAscentAnalysis>> {
    let (_starting_parcel, parcels_iter) = plume_parcels(snd, max_range, increment)?;

    Ok(parcels_iter
        .filter_map(|(_, pcl)| analyze_plume_parcel(pcl, snd).ok())
        .skip_while(|anal| anal.el_height.is_none())
        .take_while(|anal| anal.el_height.is_some())
        .collect())
}

/// Find the parcel that causes the plume to blow up by finding the maximum derivative of the
/// equilibrium level vs parcel temperature.
pub fn blow_up(snd: &Sounding) -> Result<BlowUpAnalysis> {
    const INCREMENT: CelsiusDiff = CelsiusDiff(0.1);
    const MAX_RANGE: CelsiusDiff = CelsiusDiff(20.0);

    let (starting_parcel, parcel_iter) = plume_parcels(snd, MAX_RANGE, INCREMENT)?;

    let (delta_t, _max_deriv, height): (CelsiusDiff, f64, Meters) = parcel_iter
        .filter_map(|(dt, pcl)| analyze_plume_parcel(pcl, snd).ok().map(|anal| (dt, anal)))
        .skip_while(|(_dt, anal)| anal.el_height.is_none())
        .take_while(|(_dt, anal)| anal.el_height.is_some())
        .filter_map(|(dt, anal)| anal.el_height.map(|h| (dt, h)))
        .scan(Meters(0.0), |prev, (dt, el)| {
            if el >= *prev {
                *prev = el;
                Some(Some((dt, el)))
            } else {
                Some(None)
            }
        })
        .filter_map(|opt| opt)
        .tuple_windows::<(_, _)>()
        .fold((CelsiusDiff(0.0), 0.0, Meters(0.0)), |acc, (lvl0, lvl1)| {
            let (blow_up_dt, deriv, jump) = acc;

            let (dt0, el0) = lvl0;
            let (dt1, el1) = lvl1;

            debug_assert_ne!(dt0, dt1); // Required for division below

            let derivative = (el1 - el0).unpack() / (dt1 - dt0).unpack();
            let diff = el1 - el0;
            if derivative > deriv {
                let actual_dt = CelsiusDiff((dt1 + dt0).unpack() / 2.0);
                (actual_dt, derivative, diff)
            } else {
                (blow_up_dt, deriv, jump)
            }
        });

    Ok(BlowUpAnalysis {
        starting_parcel,
        delta_t,
        height,
    })
}

/// Lift a parcel until the net CAPE is zero.
pub fn analyze_plume_parcel(parcel: Parcel, snd: &Sounding) -> Result<PlumeAscentAnalysis> {
    // Get the starting parcel and the iterator to lift it.
    let (parcel, lift_iter) = lift_parcel(parcel, snd)?;

    // Construct an iterator that selects the environment values and calculates the
    // corresponding parcel values.
    let (lcl_height, el_height, max_height, net_bouyancy) = lift_iter
        // Fold to get the EL Level, Max Height, and max integrated bouyancy
        .fold(
            (None, None, Meters(0.0), 0.0f64),
            |acc, (int_bouyancy, anal_level_type)| {
                use crate::parcel_profile::lift::AnalLevelType::*;

                let (lcl, el, _, max_bouyancy) = acc;

                match anal_level_type {
                    Normal(level) | LFC(level) => {
                        let max_hgt = level.height;
                        let max_bouyancy = max_bouyancy.max(int_bouyancy);
                        (lcl, el, max_hgt, max_bouyancy)
                    }
                    LCL(level) => {
                        let lcl = Some(level.height);
                        let max_hgt = level.height;
                        let max_bouyancy = max_bouyancy.max(int_bouyancy);
                        (lcl, el, max_hgt, max_bouyancy)
                    }
                    EL(level) => {
                        let el = Some(level.height);
                        let max_hgt = level.height;
                        let max_bouyancy = max_bouyancy.max(int_bouyancy);
                        (lcl, el, max_hgt, max_bouyancy)
                    }
                }
            },
        );

    let (net_cape, max_height) = if el_height.is_some() {
        (
            Some(JpKg(net_bouyancy / 2.0 * -metfor::g)),
            Some(max_height),
        )
    } else {
        (None, None)
    };

    Ok(PlumeAscentAnalysis {
        lcl_height,
        el_height,
        max_height,
        net_cape,
        parcel,
    })
}

/// Lift a parcel until the net CAPE is zero, don't analyze it just return the profile.
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

    // Construct an iterator that selects the environment values and calculates the
    // corresponding parcel values.
    let (lcl_height, el_height, max_height, net_bouyancy) = lift_iter
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
        // Fold to get the EL Level, Max Height, and max integrated bouyancy
        .fold(
            (None, None, Meters(0.0), 0.0f64),
            |acc, (int_bouyancy, anal_level_type)| {
                use crate::parcel_profile::lift::AnalLevelType::*;

                let (lcl, el, _, max_bouyancy) = acc;

                match anal_level_type {
                    Normal(level) | LFC(level) => {
                        let max_hgt = level.height;
                        let max_bouyancy = max_bouyancy.max(int_bouyancy);
                        (lcl, el, max_hgt, max_bouyancy)
                    }
                    LCL(level) => {
                        let lcl = Some(level.height);
                        let max_hgt = level.height;
                        let max_bouyancy = max_bouyancy.max(int_bouyancy);
                        (lcl, el, max_hgt, max_bouyancy)
                    }
                    EL(level) => {
                        let el = Some(level.height);
                        let max_hgt = level.height;
                        let max_bouyancy = max_bouyancy.max(int_bouyancy);
                        (lcl, el, max_hgt, max_bouyancy)
                    }
                }
            },
        );

    let (net_cape, max_height) = if el_height.is_some() {
        (
            Some(JpKg(net_bouyancy / 2.0 * -metfor::g)),
            Some(max_height),
        )
    } else {
        (None, None)
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
            (0.0, 0usize),
            |(int_bouyancy, count), (anal_level_type0, anal_level_type1)| {
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
                *int_bouyancy += bouyancy;
                *count += 1;
                Some(((*int_bouyancy, *count), anal_level_type1))
            },
        )
        // Take until the bouyancy goes negative, then we're done, just run through the first five
        // to ensure we get through any goofy surface layers. This is also needed to get started if
        // the parcel temperature is the same as the sounding surface temperature. Also, for the
        // case of a fire plume, near the surface things are chaotic, so we should punch through
        // a shallow surface stable layer.
        .take_while(|((int_bouyancy, count), _)| *int_bouyancy >= 0.0 || *count < 5)
        .map(|((int_bouyancy, _), anal_level)| (int_bouyancy, anal_level));

    Ok((parcel, iter))
}

/// Given a sounding, return an iterator that creates parcels starting with the mixed layer parcel
/// and then incrementing the parcel temperature up to `plus_range` in increments of `increment`.
///
/// If the surface temperature is warmer than the mixed layer parcel temperature, then the starting
/// parcel will use the surface temperature with the mixed layer moisture.
fn plume_parcels(
    snd: &Sounding,
    plus_range: CelsiusDiff,
    increment: CelsiusDiff,
) -> Result<(Parcel, impl Iterator<Item = (CelsiusDiff, Parcel)>)> {
    let parcel = mixed_layer_parcel(snd)?;
    let (row, mut parcel) = find_parcel_start_data(snd, &parcel)?;
    if row.temperature.unwrap() > parcel.temperature {
        parcel.temperature = row.temperature.unwrap();
    }

    let CelsiusDiff(inc_val) = increment;
    let next_dt = CelsiusDiff(-inc_val);

    let next_p = Parcel {
        temperature: parcel.temperature + next_dt,
        ..parcel
    };
    let max_t = parcel.temperature + plus_range;
    Ok((
        parcel,
        PlumeParcelIterator {
            next_p,
            next_dt,
            max_t,
            increment,
        },
    ))
}

/// Iterator for `plume_parcels` function that generates increasingly warmer parcels with constant
/// moisture.
struct PlumeParcelIterator {
    next_p: Parcel,
    next_dt: CelsiusDiff,
    max_t: Celsius,
    increment: CelsiusDiff,
}

impl Iterator for PlumeParcelIterator {
    type Item = (CelsiusDiff, Parcel);

    fn next(&mut self) -> Option<Self::Item> {
        let next_t = self.next_p.temperature + self.increment;
        self.next_dt += self.increment;
        if next_t > self.max_t {
            None
        } else {
            self.next_p = Parcel {
                temperature: next_t,
                ..self.next_p
            };
            Some((self.next_dt, self.next_p))
        }
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
