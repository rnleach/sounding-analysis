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
use metfor::{self, Celsius, CelsiusDiff, HectoPascal, JpKg, Kelvin, Meters, Quantity};

/// Result of lifting a parcel represntative of a fire plume core.
#[derive(Debug, Clone, Copy)]
pub struct PlumeAscentAnalysis {
    /// The parcel this analysis was completed for.
    pub parcel: Parcel,
    /// Maximum integrated buoyancy.
    pub max_int_buoyancy: Option<JpKg>,
    /// Level of maximum integrated buoyancy.
    pub level_max_int_buoyancy: Option<Meters>,
    /// Maximum integrated buoyancy without latent heat.
    pub max_dry_int_buoyancy: Option<JpKg>,
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
/// the level of max integrated buoyancy of the blow up ΔT plus 0.05C and the blow up ΔT minus
/// 0.05C. The level of max integrated buoyancy is also an equilibrium level (EL). Usually, but not
/// always, there is only a single EL which corresponds to the level of max integrated buoyancy. So
/// the level of max integrated buoyancy will be referred to as the EL, and in cases where there is
/// multiple ELs, it will be the one that corresponds to the max integrated buoyancy.
///
/// The ΔT required to get the plume top over the LCL, and thus to create a cloud is also
/// calculated. This can be a useful for determining if a plume will have a cap cloud but will
/// not be unstable enough to blow up.
///
/// After blow up of the equilibrium level, the maximum integrated buoyancy and the percentage of
/// it that is due to latent heat release are also recorded.
#[derive(Debug)]
pub struct BlowUpAnalysis {
    /// The original parcel we started with while searching for the blow up.
    pub starting_parcel: Parcel,
    /// The amount of warming required to cause a blow up of the level of maximum integrated
    /// buoyancy, or the equilibrium level.
    pub delta_t_el: CelsiusDiff,
    /// The amount of warming required to cause a cloud to form.
    pub delta_t_cloud: CelsiusDiff,
    /// The change in height from the blow up of the equilbrium level.
    pub delta_z_el: Meters,
    /// The maximum integrated buoyancy after blow up.
    pub mib: JpKg,
    /// The percentage of the `mib` that is due to latent heat release at `delta_t_el` + 1.0C
    pub pct_wet: f64,
}

impl BlowUpAnalysis {
    /// There's a lot of information packed into a `BlowUpAnalysis`, which is a simplification
    /// itself. This is an attempt to condense it even more so it could more easily used in things
    /// like climatology.
    ///
    /// The things we are looking for in an index include:
    ///  - High values for low ΔT values. Values above 10 should lead to a zero index because it is
    ///  very unlikely a blow up of any kind could occur. (In my experience so far).
    ///  - High values for high ΔZ values. Values above 10km should saturate the index because that
    ///  is a huge blow up.
    ///  - If there is no latent heat contribution, the index should be zero. It should maximize
    ///  where the contribution from latent heat balances that from sensible heat, and a blow up
    ///  that is dominated by latent heat should still show up in the index, but not as large as
    ///  one that is balanced between latent heat and sensible heat.
    pub fn as_index(&self) -> f64 {
        let dt_contribution = (10.0 - self.delta_t_el.unpack()) / 10.0;
        if dt_contribution < 0.0 {
            return 0.0;
        }

        let mut dz_contribution = self.delta_z_el.unpack() / 10_000.0;
        if dz_contribution > 1.0 {
            dz_contribution = 1.0;
        }

        debug_assert!(self.pct_wet <= 1.0 && self.pct_wet >= 0.0);
        let pct_wet_contribution = if self.pct_wet <= 0.5 {
            2.0 * self.pct_wet
        } else {
            1.5 - self.pct_wet
        };

        // Use the geometric mean of each contributing factor.
        (dt_contribution * dz_contribution * pct_wet_contribution).cbrt()
    }
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
/// level of maximum integrated buoyancy vs parcel temperature.
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

    let (delta_t_cloud, delta_t_el, _deriv_el, delta_z_el, mib, mut pct_wet): (
        CelsiusDiff,
        CelsiusDiff,
        f64,
        Meters,
        JpKg,
        f64,
    ) = parcel_iter
        // Do the analysis, ignore errors.
        .filter_map(|(dt, pcl)| analyze_plume_parcel(pcl, snd).ok().map(|anal| (dt, anal)))
        // Skip levels with no useful information.
        .skip_while(|(_dt, anal)| {
            anal.level_max_int_buoyancy.is_none()
                && anal.max_height.is_none()
                && anal.lcl_height.is_none()
        })
        // Take while there is some useful information.
        .take_while(|(_dt, anal)| {
            anal.level_max_int_buoyancy.is_some()
                || anal.max_height.is_some()
                || anal.lcl_height.is_some()
        })
        // Filter out noise due to rounding errors etc. where warmer parcels don't rise as high.
        // This smooths it out so that deriviatives work well too.
        .scan(
            (Meters(0.0), Meters(0.0), Meters(0.0)),
            |(prev_lcl, prev_lmib, prev_max_z), (dt, anal)| {
                if let Some(lcl) = anal.lcl_height {
                    if lcl >= *prev_lcl {
                        *prev_lcl = lcl;
                    } else {
                        return Some(None);
                    }
                }

                if let Some(lmib) = anal.level_max_int_buoyancy {
                    if lmib >= *prev_lmib {
                        *prev_lmib = lmib;
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
                0.0,
                Meters(0.0),
                JpKg(0.0),
                0.0,
            ),
            |acc, (lvl0, lvl1)| {
                let (
                    mut cloud_dt,
                    mut lmib_blow_up_dt,
                    mut deriv_lmib,
                    mut jump_lmib,
                    mut mib,
                    mut pct_wet,
                ) = acc;

                let (dt0, anal0) = lvl0;
                let (dt1, anal1) = lvl1;

                debug_assert_ne!(dt0, dt1); // Required for division below
                let dx = (dt1 - dt0).unpack();

                if let (Some(lmib0), Some(lmib1), Some(mib1), Some(dry_buoyancy1)) = (
                    anal0.level_max_int_buoyancy,
                    anal1.level_max_int_buoyancy,
                    anal1.max_int_buoyancy,
                    anal1.max_dry_int_buoyancy,
                ) {
                    let derivative = (lmib1 - lmib0).unpack() / dx;
                    let dt = CelsiusDiff((dt1 + dt0).unpack() / 2.0);
                    if derivative > deriv_lmib {
                        lmib_blow_up_dt = dt;
                        deriv_lmib = derivative;
                        jump_lmib = lmib1 - lmib0;
                        mib = mib1;
                    }
                    if dt.unpack() - 1.0 - lmib_blow_up_dt.unpack() <= 1.0e-4 {
                        pct_wet = (mib1 - dry_buoyancy1) / mib1;
                    }
                }

                if let (None, Some(_)) = (anal0.lcl_height, anal1.lcl_height) {
                    if cloud_dt == CelsiusDiff(0.0) {
                        cloud_dt = CelsiusDiff((dt1 + dt0).unpack() / 2.0);
                    }
                }

                (
                    cloud_dt,
                    lmib_blow_up_dt,
                    deriv_lmib,
                    jump_lmib,
                    mib,
                    pct_wet,
                )
            },
        );

    if pct_wet < 0.0 {
        pct_wet = 0.0
    } else if pct_wet > 1.0 {
        pct_wet = 1.0
    }

    Ok(BlowUpAnalysis {
        starting_parcel,
        delta_t_cloud,
        delta_t_el,
        delta_z_el,
        mib,
        pct_wet,
    })
}

/// Lift a parcel until the net CAPE is zero.
pub fn analyze_plume_parcel(parcel: Parcel, snd: &Sounding) -> Result<PlumeAscentAnalysis> {
    // Get the starting parcel and the iterator to lift it.
    let (parcel, lift_iter) = lift_parcel(parcel, snd)?;

    Ok(analyze_plume_parcel_iter(parcel, lift_iter))
}

/// Lift a parcel until the net CAPE is zero. Lift it at least 100 hPa above the surface.
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

    // An iterator that adds the values to the profile vectors as it goes.
    let lift_iter = lift_iter
        // Add the levels to the parcel profile.
        .scan(
            (),
            |_dummy, (int_buoyancy, dry_int_buoyancy, anal_level_type)| {
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

                Some((int_buoyancy, dry_int_buoyancy, anal_level_type))
            },
        );

    let plume_ascent_anal = analyze_plume_parcel_iter(parcel, lift_iter);

    Ok((
        ParcelProfile {
            pressure,
            height,
            parcel_t,
            environment_t,
        },
        plume_ascent_anal,
    ))
}

fn analyze_plume_parcel_iter(
    parcel: Parcel,
    iter: impl Iterator<Item = (f64, f64, AnalLevelType)>,
) -> PlumeAscentAnalysis {
    let mut max_height: Option<Meters> = None;

    // Construct an iterator that selects the environment values and calculates the
    // corresponding parcel values.
    let (lcl_height, el_height, max_int_buoyancy, level_max_int_buoyancy, dry_net_buoyancy) = iter
        // Scan to get the max_height
        .scan(
            (0.0, Meters(0.0)),
            |(prev_buoyancy, prev_height), (int_buoyancy, dry_int_buoyancy, anal_level_type)| {
                use crate::parcel_profile::lift::AnalLevelType::*;

                let height_val = match anal_level_type {
                    Normal(level) | LFC(level) | LCL(level) | EL(level) => level.height,
                };

                if *prev_buoyancy > 0.0 && int_buoyancy < 0.0 {
                    let mx_height =
                        linear_interp(0.0, *prev_buoyancy, int_buoyancy, *prev_height, height_val);
                    max_height = Some(mx_height);
                }

                *prev_buoyancy = int_buoyancy;
                *prev_height = height_val;

                Some((int_buoyancy, dry_int_buoyancy, anal_level_type))
            },
        )
        // Fold to get the EL Level, Max Height, and max integrated buoyancy
        .fold(
            (None, None, 0.0f64, None, 0.0f64),
            |acc, (int_buoyancy, dry_int_buoyancy, anal_level_type)| {
                use crate::parcel_profile::lift::AnalLevelType::*;

                let (
                    mut lcl,
                    mut el,
                    mut max_buoyancy,
                    mut level_max_buoyancy,
                    mut dry_max_buoyancy,
                ) = acc;

                let height = match anal_level_type {
                    Normal(level) | LFC(level) => level.height,
                    LCL(level) => {
                        lcl = Some(level.height);
                        level.height
                    }
                    EL(level) => {
                        el = Some(level.height);
                        level.height
                    }
                };

                dry_max_buoyancy = dry_max_buoyancy.max(dry_int_buoyancy);
                if int_buoyancy >= max_buoyancy {
                    max_buoyancy = int_buoyancy;
                    level_max_buoyancy = Some(height);
                }

                (lcl, el, max_buoyancy, level_max_buoyancy, dry_max_buoyancy)
            },
        );

    let (max_int_buoyancy, max_dry_int_buoyancy) = if el_height.is_some() {
        (
            Some(JpKg(max_int_buoyancy / 2.0 * -metfor::g)),
            Some(JpKg(dry_net_buoyancy / 2.0 * -metfor::g)),
        )
    } else {
        (None, None)
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

    PlumeAscentAnalysis {
        lcl_height,
        el_height,
        max_height,
        max_int_buoyancy,
        level_max_int_buoyancy,
        max_dry_int_buoyancy,
        parcel,
    }
}

/// Get the starting parcel and build an iterator to lift it.
fn lift_parcel<'a>(
    parcel: Parcel,
    snd: &'a Sounding,
) -> Result<(Parcel, impl Iterator<Item = (f64, f64, AnalLevelType)> + 'a)> {
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
    let theta0 = parcel.theta();

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
        // Pair the levels up to integrate the buoyancy.
        .tuple_windows::<(_, _)>()
        // Integrate the buoyancy.
        .scan(
            (0.0, 0.0, 0.0),
            move |(prev_int_buoyancy, int_buoyancy, dry_int_buoyancy),
                  (anal_level_type0, anal_level_type1)| {
                use crate::parcel_profile::lift::AnalLevelType::*;

                let level_data0: &AnalLevel = match &anal_level_type0 {
                    Normal(data) | LFC(data) | LCL(data) | EL(data) => data,
                };
                let level_data1: &AnalLevel = match &anal_level_type1 {
                    Normal(data) | LFC(data) | LCL(data) | EL(data) => data,
                };

                let &AnalLevel {
                    pressure: bottom_pres,
                    height: h0,
                    pcl_virt_t: pcl0,
                    env_virt_t: env0,
                } = level_data0;
                let &AnalLevel {
                    height: h1,
                    pcl_virt_t: pcl1,
                    env_virt_t: env1,
                    pressure: top_pres,
                } = level_data1;

                let dry_pcl0 = if bottom_pres > pcl_lcl.pressure {
                    Some(pcl0)
                } else {
                    let pcl_dry_t = metfor::temperature_from_theta(theta0, bottom_pres);
                    metfor::virtual_temperature(pcl_dry_t, pcl_dry_t, bottom_pres)
                        .map(Celsius::from)
                };

                let dry_pcl1 = if top_pres > pcl_lcl.pressure {
                    Some(pcl1)
                } else {
                    let pcl_dry_t = metfor::temperature_from_theta(theta0, top_pres);
                    metfor::virtual_temperature(pcl_dry_t, pcl_dry_t, bottom_pres)
                        .map(Celsius::from)
                };

                let Meters(dz) = h1 - h0;
                debug_assert!(dz >= 0.0);

                let b0 = (pcl0 - env0).unpack() / Kelvin::from(env0).unpack();
                let b1 = (pcl1 - env1).unpack() / Kelvin::from(env1).unpack();
                let buoyancy = (b0 + b1) * dz;

                if let (Some(dry_pcl0), Some(dry_pcl1)) = (dry_pcl0, dry_pcl1) {
                    let db0 = (dry_pcl0 - env0).unpack() / Kelvin::from(env0).unpack();
                    let db1 = (dry_pcl1 - env1).unpack() / Kelvin::from(env1).unpack();
                    let dry_buoyancy = (db0 + db1) * dz;
                    *dry_int_buoyancy += dry_buoyancy;
                }

                *prev_int_buoyancy = *int_buoyancy;
                *int_buoyancy += buoyancy;

                *dry_int_buoyancy = dry_int_buoyancy.min(*int_buoyancy);

                Some((
                    (
                        *prev_int_buoyancy,
                        *int_buoyancy,
                        *dry_int_buoyancy,
                        bottom_pres,
                    ),
                    anal_level_type1,
                ))
            },
        )
        // Take until the buoyancy goes negative, then we're done, just run through the first five
        // to ensure we get through any goofy surface layers. This is also needed to get started if
        // the parcel temperature is the same as the sounding surface temperature. Also, for the
        // case of a fire plume, near the surface things are chaotic, so we should punch through
        // a shallow surface stable layer.
        //
        // Use the prev_int_buoyancy to at least one point past where the buoyancy becomes zero
        // so we have enough data to interpolate.
        .take_while(move |((prev_int_buoyancy, _, _, pres), _)| {
            *prev_int_buoyancy >= 0.0 || *pres > p0 - HectoPascal(100.0)
        })
        .map(|((_, int_buoyancy, dry_int_buoyancy, _), anal_level)| {
            (int_buoyancy, dry_int_buoyancy, anal_level)
        });

    Ok((parcel, iter))
}

/// Given a sounding, return an iterator that creates parcels starting with the mixed layer parcel
/// and then incrementing the parcel temperature up to `plus_range` in increments of `increment`.
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
    let (_row, parcel) = find_parcel_start_data(snd, &parcel)?;

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
                        let buoyancy = ((pt - et).unpack() / et.unpack()
                            + (prev_pt - prev_et).unpack() / prev_et.unpack())
                            * dz.unpack();
                        cape += buoyancy;

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
