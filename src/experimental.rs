//! Functions for analysis functions that I'm experimenting with.

use crate::{
    error::{AnalysisError, Result},
    parcel::{mixed_layer_parcel, Parcel},
    parcel_profile::find_parcel_start_data,
};
use itertools::izip;
use metfor::{Celsius, CelsiusDiff, JpKg, Kelvin, Meters, Quantity};
use sounding_base::Sounding;

/// An analysis of the potential energy of a convective plume vs. a representative starting
/// temperature.
#[derive(Debug)]
pub struct PlumePotentialAnal {
    // Representative initial plume temperature and CAPE associated with that parcel.
    //
    // In reality it is much hotter at the surface, but there is a lot of entrainment to come still.
    plume_vals: Vec<(Celsius, JpKg)>,
    // Index of last of `plume_vals` which do not create a cloud.
    last_no_cloud: usize,
    // Index of first of `plum_vals` which does create a cloud.
    first_with_cloud: Option<usize>,
}

// Create another lift algorithm that calculates CAPE no matter what.
// Once the LCL_P>EL_P, change indexes last_no_cloud and first_with_cloud

impl PlumePotentialAnal {
    /// Analyze a sounding for the plume potential
    pub fn analyze(snd: &Sounding) -> Result<PlumePotentialAnal> {
        let mut last_no_cloud = 0;
        let mut first_with_cloud = None;
        let mut plume_vals = Vec::with_capacity(30);

        for (i, pcl) in plume_parcels(snd, CelsiusDiff(30.0), CelsiusDiff(1.0))?.enumerate() {
            let (cape, cloud) = lift_parcel(pcl, snd)?;
            plume_vals.push((pcl.temperature, cape));

            if !cloud && first_with_cloud.is_none() {
                last_no_cloud = i;
            } else if cloud && first_with_cloud.is_none() {
                first_with_cloud = Some(i);
            }
        }

        Ok(PlumePotentialAnal {
            plume_vals,
            last_no_cloud,
            first_with_cloud,
        })
    }

    /// Get the values of The plume surface temperature and associated CAPE
    pub fn analysis_values(&self) -> &[(Celsius, JpKg)] {
        &self.plume_vals
    }

    /// Get the temperature and energy associated with the last plume that did NOT create a cloud.
    pub fn t0_e0(&self) -> (Celsius, JpKg) {
        self.plume_vals[self.last_no_cloud]
    }

    /// Get the energy jump associated with the cloud creation. If no cloud was created ever, this
    /// will be `None`.
    pub fn delta_e(&self) -> Option<JpKg> {
        self.first_with_cloud.map(|first_with_cloud| {
            self.plume_vals[first_with_cloud].1 - self.plume_vals[self.last_no_cloud].1
        })
    }
}

// Given a sounding, return an iterator that creates parcels from the surface temperature to +30C
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
            metfor::virtual_temperature(env_t, env_dp, p)
                .map(|env_vt| (p, h, Kelvin::from(env_vt), pcl_t))
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
