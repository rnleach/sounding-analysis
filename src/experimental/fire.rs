//! Experimental sounding analysis for fire weather.

use crate::{
    error::{AnalysisError, Result},
    parcel::{convective_parcel, lowest_level_parcel, mixed_layer_parcel, Parcel},
    parcel_profile::{find_parcel_start_data, ParcelAscentAnalysis,
    lift::{parcel_lcl, create_parcel_calc_t},
    },
    sounding::Sounding,
};
use itertools::{izip, Itertools};
use metfor::{Celsius, CelsiusDiff, JpKg, Kelvin, Meters, Quantity};

/// Result of lifting a parcel represntative of a plume core.
pub struct PlumeAscentAnalysis {
    /// The parcel this analysis was completed for.
    pub parcel: Parcel,
    /// Net CAPE up to equilibrium level.
    pub net_cape: JpKg,
    /// The level where net CAPE becomes zero, the plume rises no more
    pub max_height: Meters,
    /// The last EL reached before max_height
    pub el_height: Meters,
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
pub struct BlowUpAnalysis {
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
    unimplemented!()
}

/// Find the parcel that causes the plume to blow up by finding the maximum derivative of the
/// equilibrium level vs parcel temperature.
pub fn blow_up(snd: &Sounding) -> Result<BlowUpAnalysis> {
    unimplemented!()
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
) -> Result<impl Iterator<Item = Parcel>> {
    let parcel = mixed_layer_parcel(snd)?;
    let (row, mut parcel) = find_parcel_start_data(snd, &parcel)?;
    if row.temperature.unwrap() > parcel.temperature {
        parcel.temperature = row.temperature.unwrap();
    }

    let max_t = parcel.temperature + plus_range;
    Ok(PlumeParcelIterator {
        next_p: parcel,
        max_t,
        increment,
    })
}

/// Iterator for `plume_parcels` function that generates increasingly warmer parcels with constant
/// moisture.
struct PlumeParcelIterator {
    next_p: Parcel,
    max_t: Celsius,
    increment: CelsiusDiff,
}

impl Iterator for PlumeParcelIterator {
    type Item = Parcel;
    fn next(&mut self) -> Option<Self::Item> {
        let next_t = self.next_p.temperature + self.increment;
        if next_t > self.max_t {
            None
        } else {
            self.next_p = Parcel {
                temperature: next_t,
                ..self.next_p
            };
            Some(self.next_p)
        }
    }
}

/// Lift a parcel until the net CAPE is zero.
pub fn lift_plume_parcel(parcel: Parcel, snd: &Sounding) -> Result<PlumeAscentAnalysis> {
    // Find the LCL
    let (pcl_lcl, lcl_temperature) = parcel_lcl(&parcel, snd)?;

    // The starting level to lift the parcel from
    let (parcel_start_data, parcel) = find_parcel_start_data(snd, &parcel)?;

    // How to calculate a parcel temperature for a given pressure level
    let parcel_calc_t = create_parcel_calc_t(parcel, pcl_lcl)?;


    unimplemented!()
}

// -----------------------------------------------------------------------------------------------
/// An analysis of the potential energy of a convective plume vs. a representative starting
/// temperature.
#[derive(Debug)]
pub struct PlumePotentialAnal {
    // Representative initial plume temperature and CAPE associated with that parcel.
    //
    // In reality it is much hotter at the surface, but there is a lot of entrainment to come still.
    plume_vals: Vec<(Celsius, JpKg)>,
}

// Create another lift algorithm that calculates CAPE no matter what.
// Once the LCL_P>EL_P, change indexes last_no_cloud and first_with_cloud

impl PlumePotentialAnal {
    /// Analyze a sounding for the plume potential
    pub fn analyze(snd: &Sounding) -> Result<PlumePotentialAnal> {
        const MAX_HEATING: CelsiusDiff = CelsiusDiff(30.0);
        const DT: CelsiusDiff = CelsiusDiff(0.1);

        let mut plume_vals = Vec::with_capacity(MAX_HEATING.unpack() as usize);

        for pcl in plume_parcels(snd, MAX_HEATING, DT)? {
            if let Ok((cape, _)) = lift_parcel(pcl, snd) {
                plume_vals.push((pcl.temperature, cape));
            } else {
                break;
            }
        }

        Ok(PlumePotentialAnal { plume_vals })
    }

    /// Get the values of The plume surface temperature and associated CAPE
    pub fn analysis_values(&self) -> &[(Celsius, JpKg)] {
        &self.plume_vals
    }
}

/// Analyze the sounding to get the values of (t˳, Δt˳, E˳, ΔE).
pub fn convective_parcel_initiation_energetics(
    snd: &Sounding,
) -> Result<(Celsius, CelsiusDiff, JpKg, JpKg)> {
    let starting_parcel = convective_parcel(snd)?;

    let mut no_cloud_pcl = starting_parcel;
    let mut no_cloud_pcl_data = lift_parcel(starting_parcel, snd)?;
    let mut cloud_pcl = starting_parcel;
    let mut cloud_pcl_data = no_cloud_pcl_data;

    // bracket the cloud/no cloud boundary
    let tgt_cloud_val = if no_cloud_pcl_data.1 {
        false
    } else if !cloud_pcl_data.1 {
        true
    } else {
        unreachable!()
    };

    if tgt_cloud_val {
        // Cloud parcel doesn't have a cloud! Increment until it does
        cloud_pcl.temperature += CelsiusDiff(1.0);
        cloud_pcl_data = lift_parcel(cloud_pcl, snd)?;
    } else {
        // No cloud parcel has a cloud! Decrement until it doesn't
        no_cloud_pcl.temperature += CelsiusDiff(-1.0);
        no_cloud_pcl_data = lift_parcel(no_cloud_pcl, snd)?;
    }
    while no_cloud_pcl_data.1 || !cloud_pcl_data.1 {
        if tgt_cloud_val {
            // Cloud parcel doesn't have a cloud! Increment until it does
            cloud_pcl.temperature += CelsiusDiff(1.0);
            cloud_pcl_data = lift_parcel(cloud_pcl, snd)?;
        } else {
            // No cloud parcel has a cloud! Decrement until it doesn't
            no_cloud_pcl.temperature += CelsiusDiff(-1.0);
            no_cloud_pcl_data = lift_parcel(no_cloud_pcl, snd)?;
        }
    }

    // use bisection to narrow the gap between no-cloud and cloud to 0.1°C
    let mut dt = cloud_pcl.temperature - no_cloud_pcl.temperature;
    while dt > CelsiusDiff(0.1) {
        let mid_t = no_cloud_pcl.temperature + CelsiusDiff(dt.unpack() / 2.0);
        let test_pcl = Parcel {
            temperature: mid_t,
            ..no_cloud_pcl
        };
        let test_pcl_data = lift_parcel(test_pcl, snd)?;

        if test_pcl_data.1 {
            // In Cloud!
            cloud_pcl = test_pcl;
            cloud_pcl_data = test_pcl_data;
        } else {
            // Not in a cloud
            no_cloud_pcl = test_pcl;
            no_cloud_pcl_data = test_pcl_data;
        }

        dt = cloud_pcl.temperature - no_cloud_pcl.temperature;
    }

    let t0 = no_cloud_pcl.temperature;
    let dt0 = t0 - lowest_level_parcel(snd)?.temperature;
    let e0 = no_cloud_pcl_data.0;
    let d_e = cloud_pcl_data.0 - e0;

    // return (t˳, E˳, ΔE)
    Ok((t0, dt0, e0, d_e))
}

fn lift_parcel(parcel: Parcel, snd: &Sounding) -> Result<(JpKg, bool)> {
    //
    // Find the LCL
    //
    let (lcl_pressure, _lcl_temperature) = metfor::pressure_and_temperature_at_lcl(
        parcel.temperature,
        parcel.dew_point,
        parcel.pressure,
    )
    .ok_or(AnalysisError::MetForError)?;

    //
    // The starting level to lift the parcel from
    //
    let (_, parcel) = find_parcel_start_data(snd, &parcel)?;

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
            .map(Kelvin::from)
        } else {
            // Moist adiabatic lifting
            metfor::temperature_from_theta_e_saturated_and_pressure(tgt_pres, theta_e)
                .and_then(|t_c| metfor::virtual_temperature(t_c, t_c, tgt_pres))
                .map(Kelvin::from)
        }
    };

    //
    // Get the environment data to iterate over.
    //
    let snd_pressure = snd.pressure_profile();
    let hgt = snd.height_profile();
    let env_t = snd.temperature_profile();
    let env_dp = snd.dew_point_profile();

    //
    // Construct an iterator that selects the environment values and calculates the
    // corresponding parcel values.
    //
    let (integrated_cape, cloud) = izip!(snd_pressure, hgt, env_t, env_dp)
        // Remove rows with missing data and unpack `optional::Optioned` types
        .filter_map(|(p, h, env_t, env_dp)| {
            if p.is_some() && h.is_some() && env_t.is_some() && env_dp.is_some() {
                Some((p.unpack(), h.unpack(), env_t.unpack(), env_dp.unpack()))
            } else {
                None
            }
        })
        // Remove rows at or below the parcel level
        .filter(move |(p, _, _, _)| *p < parcel.pressure)
        // Calculate the parcel temperature, skip this level if there is an error
        .filter_map(|(p, h, env_t, env_dp)| {
            calc_parcel_t(p).map(|pcl_t| (p, h, env_t, env_dp, pcl_t))
        })
        // Calculate the environment virtual temperature, skip levels with errors
        .filter_map(|(p, h, env_t, env_dp, pcl_t)| {
            metfor::virtual_temperature(env_t, env_dp, p).map(|env_vt| (p, h, env_vt, pcl_t))
        })
        // Take while parcel is warmer than the environment. We don't need to worry about
        // interpolating to the exact equilibrium level, because at the equilibrium level the
        // the contribution to the integraged cape is 0.0. So just don't get any negative by going
        // further!
        .take_while(|(_p, _h, env_t, pcl_t)| pcl_t > env_t)
        // Integrate witht the trapezoid rule
        .fold(
            (
                (0.0, false), // CAPE integral, cloud initiated
                Meters(std::f64::MAX),
                Kelvin(0.0),
                Kelvin(0.0),
            ),
            |acc, (p, h, env_t, pcl_t)| {
                let ((mut cape, mut cloud), prev_h, prev_env_t, prev_pcl_t) = acc;

                let dz = h - prev_h;

                if dz <= Meters(0.0) {
                    // Must be the first iteration, pass on the "old" values
                    ((cape, cloud), h, env_t, pcl_t)
                } else {
                    cape += ((pcl_t - env_t).unpack() / env_t.unpack()
                        + (prev_pcl_t - prev_env_t).unpack() / prev_env_t.unpack())
                        * dz.unpack();

                    cloud = p < lcl_pressure;

                    ((cape, cloud), h, env_t, pcl_t)
                }
            },
        )
        .0;

    Ok((JpKg(integrated_cape / 2.0 * -metfor::g), cloud))
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
