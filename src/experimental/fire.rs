//! Experimental sounding analysis for fire weather.

use crate::{
    error::Result,
    interpolation::linear_interp,
    parcel::{mixed_layer_parcel, Parcel},
    parcel_profile::{
        find_parcel_start_data,
        lift::{
            create_level_type_mapping, create_parcel_calc_t, parcel_lcl, AnalLevel, AnalLevelType,
        },
        ParcelProfile,
    },
    sounding::Sounding,
};
use itertools::{izip, Itertools};
use metfor::{self, Celsius, CelsiusDiff, HectoPascal, JpKg, Kelvin, Meters, Quantity};
use optional::{none, some, Optioned};

/// Result of lifting a parcel represntative of a fire plume core.
#[derive(Debug, Clone, Copy)]
pub struct PlumeAscentAnalysis {
    /// The parcel this analysis was completed for.
    pub parcel: Parcel,
    /// Maximum integrated buoyancy.
    pub max_int_buoyancy: Optioned<JpKg>,
    /// Maximum integrated buoyancy without latent heat.
    pub max_dry_int_buoyancy: Optioned<JpKg>,
    /// The lifting condensation level of the parcel.
    pub lcl_height: Optioned<Meters>,
    /// The equilibrium level. If there are multiple equilibrium levels, this is the one that
    /// corresponds to the maximum integrated buoyancy.
    pub el_height: Optioned<Meters>,
    /// The level where net CAPE becomes zero, the plume rises no more
    pub max_height: Optioned<Meters>,
    /// The elevation of the surface.
    pub sfc_height: Meters,
}

/// Various analysis results of lifting plumes parcels vs the heating supplied.
#[derive(Debug, Clone)]
pub struct PlumeHeatingAnalysis {
    /// The parcel this analysis was completed for.
    pub starting_parcel: Parcel,
    /// Ordinate
    pub dts: Vec<CelsiusDiff>,
    /// Buoyancies
    pub max_int_buoyancies: Vec<Optioned<JpKg>>,
    /// Wet Ratio 0-1
    pub wet_ratio: Vec<Optioned<f64>>,
    /// LCL
    pub lcl_heights: Vec<Optioned<Meters>>,
    /// EL
    pub el_heights: Vec<Optioned<Meters>>,
    /// Max Plume Height
    pub max_heights: Vec<Optioned<Meters>>,
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

impl BlowUpAnalysis {}

/// Generate a series of `PlumeAscentAnalysis`s for a given sounding starting at the mixed layer
/// parcel and warming it by `increment` until it has gone `max_range` degrees.
///
/// Arguments:
/// * snd - the environmental sounding.
/// * min_dt - the minimum amount of heating to apply, may be negative.
/// * max_dt - the maximum amount of heating to apply.
/// * increment - the difference in temperature between parcels.
/// * moisture_ratio - a value of 10 means for every 10C of heating (dt), add 1 g/kg of moisture
///   to the parcel. If it is `None`, don't add any moisture.
///
pub fn calc_plumes(
    snd: &Sounding,
    increment: CelsiusDiff,
    min_dt: CelsiusDiff,
    max_dt: CelsiusDiff,
    moisture_ratio: Option<f64>,
) -> Result<Vec<PlumeAscentAnalysis>> {
    let (_starting_parcel, parcels_iter) =
        plume_parcels(snd, min_dt, max_dt, increment, moisture_ratio)?;

    Ok(parcels_iter
        .filter_map(|(_, pcl)| analyze_plume_parcel(pcl, snd).ok())
        .skip_while(|anal| anal.el_height.is_none())
        .take_while(|anal| anal.el_height.is_some())
        .collect())
}

fn plumes_heating_iter(
    snd: &Sounding,
    moisture_ratio: Option<f64>,
) -> Result<(
    Parcel,
    impl Iterator<Item = (CelsiusDiff, PlumeAscentAnalysis)> + '_,
)> {
    const INCREMENT: CelsiusDiff = CelsiusDiff(0.1);
    const MIN_DT: CelsiusDiff = CelsiusDiff(-1.0);
    const MAX_DT: CelsiusDiff = CelsiusDiff(20.0);

    let (starting_parcel, parcel_iter) =
        plume_parcels(snd, MIN_DT, MAX_DT, INCREMENT, moisture_ratio)?;

    let anal_iter = parcel_iter
        // Do the analysis, ignore errors.
        .filter_map(move |(dt, pcl)| analyze_plume_parcel(pcl, snd).ok().map(|anal| (dt, anal)))
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
            |(prev_lcl, prev_lmib, prev_max_z), (dt, anal)| {
                if let (Some(lcl), Some(lmib), Some(max_z)) = (
                    anal.lcl_height.into_option(),
                    anal.el_height.into_option(),
                    anal.max_height.into_option(),
                ) {
                    if lcl >= *prev_lcl && lmib >= *prev_lmib && max_z >= *prev_max_z {
                        *prev_lcl = lcl;
                        *prev_lmib = lmib;
                        *prev_max_z = max_z;
                    } else {
                        return Some(None);
                    }
                }

                Some(Some((dt, anal)))
            },
        )
        // Filter out "bad" layers.
        .filter_map(|opt| opt);

    Ok((starting_parcel, anal_iter))
}

/// Find the parcel that causes the plume to blow up by finding the maximum derivative of the
/// equilibrium level vs parcel heating.
///
/// Arguments:
/// * snd - the environmental sounding.
/// * moisture_ratio - a value of 10 means for every 10C of heating (dt), add 1 g/kg of moisture
///   to the parcel. If it is `None`, don't add any moisture.
///
pub fn blow_up(snd: &Sounding, moisture_ratio: Option<f64>) -> Result<BlowUpAnalysis> {
    let (starting_parcel, anal_iter) = plumes_heating_iter(snd, moisture_ratio)?;

    let (anal_iter0, anal_iter) = anal_iter.tee();

    if let Some(immediate_blow_up) = check_for_immediate_blow_up(starting_parcel, anal_iter0) {
        return Ok(immediate_blow_up);
    }

    let (mut dt_el, mut delta_z_el, mut mib) = (CelsiusDiff(0.0), Meters(0.0), JpKg(0.0));
    let mut pct_wet = 0.0;
    let mut dt_cloud = CelsiusDiff(0.0);

    anal_iter
        // Extract necessary parts and filter out points missing critical data.
        .filter_map(|(dt, anal)| {
            if let (Some(el), Some(mib), Some(dry_mib)) = (
                anal.el_height.into_option(),
                anal.max_int_buoyancy.into_option(),
                anal.max_dry_int_buoyancy.into_option(),
            ) {
                Some((dt, el, mib, dry_mib, anal.lcl_height))
            } else {
                None
            }
        })
        // Pair up for simple derivative calculations.
        .tuple_windows::<(_, _)>()
        // Scan derivative for dt_el and delta_z_el
        .scan(0.0f64, |deriv, (lvl0, lvl1)| {
            let (dt0, el0, _, _, _) = lvl0;
            let (dt1, el1, mib1, _, _) = lvl1;

            debug_assert_ne!(dt0, dt1); // Required for division below
            let dx = dt1 - dt0;

            let derivative = (el1 - el0).unpack() / dx.unpack();
            let dt = (dt1 + dt0) / 2.0;
            if derivative > *deriv {
                dt_el = dt;
                *deriv = derivative;
                delta_z_el = el1 - el0;
                mib = mib1;
            }

            Some((dt, dt_el, lvl0, lvl1))
        })
        // Inspect to find pct_wet at dt_el + CelsiusDiff(1.0)
        .map(|(dt, curr_dt_el, lvl0, lvl1)| {
            let (_, _, mib0, dry0, lcl0) = lvl0;
            let (_, _, mib1, dry1, lcl1) = lvl1;

            if (dt - CelsiusDiff(1.0) - curr_dt_el).abs() <= CelsiusDiff(1.0e-4) {
                let avg_mib = (mib0 + mib1) / 2.0;
                let avg_dryb = (dry0 + dry1) / 2.0;
                pct_wet = (avg_mib - avg_dryb) / avg_mib;
            }

            (dt, lcl0, lcl1)
        })
        // Find dt_cloud
        .for_each(|(dt, lcl0, lcl1)| {
            if lcl0.is_none() && lcl1.is_some() && dt_cloud == CelsiusDiff(0.0) {
                dt_cloud = dt;
            }
        });

    Ok(BlowUpAnalysis {
        starting_parcel,
        delta_t_cloud: dt_cloud,
        delta_t_el: dt_el,
        delta_z_el,
        mib,
        pct_wet,
    })
}

/// This checks if the plume blows up at the very first dt step. This causes the derivative method
/// to fail, so it is a special case.
fn check_for_immediate_blow_up(
    starting_parcel: Parcel,
    mut iter: impl Iterator<Item = (CelsiusDiff, PlumeAscentAnalysis)>,
) -> Option<BlowUpAnalysis> {
    if let Some((dt0, el0, lcl_opt, mib0, dry_mib0, sfc_height)) = iter.next().and_then(|anal| {
        anal.1.el_height.into_option().map(|el| {
            (
                anal.0,
                el,
                anal.1.lcl_height,
                anal.1.max_int_buoyancy,
                anal.1.max_dry_int_buoyancy,
                anal.1.sfc_height,
            )
        })
    }) {
        if el0 > Meters(5000.0) {
            let dt_cloud;
            if lcl_opt.is_some() {
                dt_cloud = dt0;
            } else {
                dt_cloud = CelsiusDiff(0.0);
            }

            let dt_el = dt0;
            let delta_z_el = el0 - sfc_height;

            let mib;
            let pct_wet;
            if let Some(mib0) = mib0.into_option() {
                mib = mib0;

                if let Some(dry_mib0) = dry_mib0.into_option() {
                    pct_wet = (mib0 - dry_mib0) / mib0;
                } else {
                    pct_wet = 0.0;
                }
            } else {
                mib = JpKg(0.0);
                pct_wet = 0.0;
            }

            Some(BlowUpAnalysis {
                starting_parcel,
                delta_t_cloud: dt_cloud,
                delta_t_el: dt_el,
                delta_z_el,
                mib,
                pct_wet,
            })
        } else {
            None
        }
    } else {
        // Technically this is probably an error, but not the one we're looking for, so let the caller
        // detect and handle it.
        None
    }
}

/// Do a PlumeHeatingAnalysis.
pub fn plume_heating_analysis(
    snd: &Sounding,
    moisture_ratio: Option<f64>,
) -> Result<PlumeHeatingAnalysis> {
    let (starting_parcel, anal_iter) = plumes_heating_iter(snd, moisture_ratio)?;

    let mut dts: Vec<CelsiusDiff> = vec![];
    let mut max_int_buoyancies: Vec<Optioned<JpKg>> = vec![];
    let mut wet_ratio: Vec<Optioned<f64>> = vec![];
    let mut lcl_heights: Vec<Optioned<Meters>> = vec![];
    let mut el_heights: Vec<Optioned<Meters>> = vec![];
    let mut max_heights: Vec<Optioned<Meters>> = vec![];

    anal_iter
        .for_each(|(dt, anal)| {
            dts.push(dt);
            max_int_buoyancies.push(anal.max_int_buoyancy);

            let a_wet_ratio = anal.max_dry_int_buoyancy.and_then(|dry| {
                anal.max_int_buoyancy.map_t(|total| {
                    if total > JpKg(0.0) {
                        (total - dry) / total
                    } else {
                        0.0
                    }
                })
            });
            wet_ratio.push(a_wet_ratio);
            lcl_heights.push(anal.lcl_height);
            el_heights.push(anal.el_height);
            max_heights.push(anal.max_height);
        });

    Ok(PlumeHeatingAnalysis {
        starting_parcel,
        dts,
        max_int_buoyancies,
        wet_ratio,
        lcl_heights,
        el_heights,
        max_heights,
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
    let mut max_height: Optioned<Meters> = none();
    let mut sfc_height = Meters(std::f64::MAX);

    // Construct an iterator that selects the environment values and calculates the
    // corresponding parcel values.
    let (lcl_height, el_height, max_int_buoyancy, dry_net_buoyancy) = iter
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
                    max_height = some(mx_height);
                }

                // Find the lowest level and call it the surface
                if sfc_height > height_val {
                    sfc_height = height_val;
                }

                *prev_buoyancy = int_buoyancy;
                *prev_height = height_val;

                Some((int_buoyancy, dry_int_buoyancy, anal_level_type))
            },
        )
        // Fold to get the EL Level, LCL Height, and max integrated buoyancy
        .fold(
            (none(), none(), 0.0f64, 0.0f64),
            |acc, (int_buoyancy, dry_int_buoyancy, anal_level_type)| {
                use crate::parcel_profile::lift::AnalLevelType::*;

                let (mut lcl, mut el, mut max_buoyancy, mut dry_max_buoyancy) = acc;

                let height = match anal_level_type {
                    Normal(level) | LFC(level) | EL(level) => level.height,
                    LCL(level) => {
                        lcl = some(level.height);
                        level.height
                    }
                };

                dry_max_buoyancy = dry_max_buoyancy.max(dry_int_buoyancy);
                if int_buoyancy >= max_buoyancy {
                    max_buoyancy = int_buoyancy;
                    el = some(height);
                }

                (lcl, el, max_buoyancy, dry_max_buoyancy)
            },
        );

    let (max_int_buoyancy, max_dry_int_buoyancy) = if el_height.is_some() {
        (
            some(JpKg(max_int_buoyancy / 2.0 * -metfor::g)),
            some(JpKg(dry_net_buoyancy / 2.0 * -metfor::g)),
        )
    } else {
        (none(), none())
    };

    // Make sure LCL is below max height
    let lcl_height =
        if let (Some(mxh), Some(lcl)) = (max_height.into_option(), lcl_height.into_option()) {
            if mxh > lcl {
                some(lcl)
            } else {
                none()
            }
        } else {
            lcl_height
        };

    if sfc_height > Meters(std::f64::MAX / 2.0) {
        debug_assert!(
            sfc_height > Meters(std::f64::MAX / 2.0),
            "sfc_heigh never assigned."
        );
        sfc_height = Meters(0.0);
    }

    PlumeAscentAnalysis {
        lcl_height,
        el_height,
        max_height,
        max_int_buoyancy,
        max_dry_int_buoyancy,
        parcel,
        sfc_height,
    }
}

/// Get the starting parcel and build an iterator to lift it.
fn lift_parcel(
    parcel: Parcel,
    snd: &Sounding,
) -> Result<(Parcel, impl Iterator<Item = (f64, f64, AnalLevelType)> + '_)> {
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
                    let pcl_dry_t = metfor::temperature_from_pot_temp(theta0, bottom_pres);
                    metfor::virtual_temperature(pcl_dry_t, pcl_dry_t, bottom_pres)
                        .map(Celsius::from)
                };

                let dry_pcl1 = if top_pres > pcl_lcl.pressure {
                    Some(pcl1)
                } else {
                    let pcl_dry_t = metfor::temperature_from_pot_temp(theta0, top_pres);
                    metfor::virtual_temperature(pcl_dry_t, pcl_dry_t, bottom_pres)
                        .map(Celsius::from)
                };

                let Meters(dz) = h1 - h0;
                debug_assert!(dz >= 0.0);

                let b0 = (pcl0 - env0) / Kelvin::from(env0);
                let b1 = (pcl1 - env1) / Kelvin::from(env1);
                let buoyancy = (b0 + b1) * dz;

                if let (Some(dry_pcl0), Some(dry_pcl1)) = (dry_pcl0, dry_pcl1) {
                    let db0 = (dry_pcl0 - env0) / Kelvin::from(env0);
                    let db1 = (dry_pcl1 - env1) / Kelvin::from(env1);
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
/// * min_dt - the minimum amount of heating to apply, which may be negative.
/// * max_dt - the maximum amount of heating to apply.
/// * increment - the difference in temperature between parcels.
/// * moisture_ratio - a value of 10 means for every 10C of heating (dt), add 1 g/kg of moisture
///   to the parcel. If it is `None`, don't add any moisture.
///
fn plume_parcels(
    snd: &Sounding,
    min_dt: CelsiusDiff,
    max_dt: CelsiusDiff,
    increment: CelsiusDiff,
    moisture_ratio: Option<f64>,
) -> Result<(Parcel, impl Iterator<Item = (CelsiusDiff, Parcel)>)> {
    let parcel = mixed_layer_parcel(snd)?;
    let (_row, parcel) = find_parcel_start_data(snd, &parcel)?;

    let next_dt = min_dt - increment;

    Ok((
        parcel,
        PlumeParcelIterator {
            starting_pcl: parcel,
            next_dt,
            max_dt,
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
