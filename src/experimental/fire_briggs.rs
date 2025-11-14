//! Experimental sounding analysis for fire plumes utilizing the Briggs plume model for the
//! subcloud plume.

use crate::{
    error::{AnalysisError, Result},
    fire::{
        entrained_layer_mean_wind_speed, entrained_mixed_layer,
        potential_t_and_specific_humidity_to_pressure_and_temperature,
    },
    sounding::Sounding,
};
use itertools::{izip, Itertools};
use metfor::{self, Celsius, GigaWatts, HectoPascal, JpKg, Kelvin, Meters, Quantity};
use optional::{none, some, Optioned};

/// Various analysis results of lifting plumes parcels vs the heating supplied. The subcloud
/// portion of the plume is modeled with the Briggs plume model.
#[derive(Debug, Clone)]
pub struct BriggsPlumeHeatingAnalysis {
    /// How much heating does the fire need to produce for these results.
    pub fire_power: Vec<GigaWatts>,
    /// The distance along the SP-curve moved to generate each plume
    pub betas: Vec<f64>,
    /// Basically CAPE, but it takes the heating of the fire into account as well.
    pub max_int_buoyancies: Vec<Optioned<JpKg>>,
    /// The ratio of the CAPE that was from latent heat release.
    pub wet_ratio: Vec<Optioned<f64>>,
    /// LCL - cloud base.
    pub lcl_heights: Vec<Optioned<Meters>>,
    /// Equiplibrium level (leve of maximum integrated buoyancy.
    pub el_heights: Vec<Optioned<Meters>>,
    /// If we keep integrating above the EL, how far up before the buoyancy gets to zero. The max
    /// height minus the EL would be the depth of an overshooting top.
    pub max_heights: Vec<Optioned<Meters>>,
    /// The potential temperature of the entrained mixed layer and the start of the SP-curve
    pub starting_theta: Kelvin,
    /// The specific humidity of the entrained mixed layer and the start of the SP-curve
    pub starting_sh: f64,
    /// The moisture ratio used to calculate this plume.
    pub moisture_ratio: Option<f64>,
    /// The surface height used in the calculation 
    pub sfc_height: Meters,
    /// The surface pressure.
    pub p_sfc: HectoPascal,
}

/// Do a PlumeHeatingAnalysis.
pub fn briggs_plume_heating_analysis(
    snd: &Sounding,
    moisture_ratio: Option<f64>,
) -> Result<BriggsPlumeHeatingAnalysis> {
    let (sfc_height, p_sfc, starting_theta, starting_sh, anal_iter) = plumes_heating_iter(snd, moisture_ratio)?;

    let mut fire_power: Vec<GigaWatts> = vec![];
    let mut betas: Vec<f64> = vec![];
    let mut max_int_buoyancies: Vec<Optioned<JpKg>> = vec![];
    let mut wet_ratio: Vec<Optioned<f64>> = vec![];
    let mut lcl_heights: Vec<Optioned<Meters>> = vec![];
    let mut el_heights: Vec<Optioned<Meters>> = vec![];
    let mut max_heights: Vec<Optioned<Meters>> = vec![];

    anal_iter.for_each(|anal| {
        fire_power.push(anal.fire_power);
        betas.push(anal.beta);
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

    Ok(BriggsPlumeHeatingAnalysis {
        fire_power,
        betas,
        max_int_buoyancies,
        wet_ratio,
        lcl_heights,
        el_heights,
        max_heights,
        starting_theta,
        starting_sh,
        moisture_ratio, 
        sfc_height,
        p_sfc,
    })
}

/// Result of lifting a parcel represntative of a fire plume core.
#[derive(Debug, Clone )]
pub struct PlumeAscentAnalysis {
    /// The fire power required to make this plume
    pub fire_power: GigaWatts,
    /// The beta value used along the SP curve to calculate this plume.
    pub beta: f64,
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
    /// The pressure portion of the pluem profile.
    pub p_profile: Vec<HectoPascal>,
    /// The temperature portion of the plume profile.
    pub t_profile: Vec<Celsius>,
}

struct PlumeIterator<'a> {
    next_beta: f64,
    max_beta: f64,
    increment: f64,
    moisture_ratio: Option<f64>,
    starting_theta: Kelvin,
    starting_sh: f64,
    sfc_height: Meters,
    p_sfc: HectoPascal,
    snd: &'a Sounding,
}

impl Iterator for PlumeIterator<'_> {
    type Item = PlumeAscentAnalysis;

    fn next(&mut self) -> Option<Self::Item> {
        // Sounding profiles I'll use
        let heights = self.snd.height_profile();
        let pressures = self.snd.pressure_profile();
        let temperatures = self.snd.temperature_profile();

        loop {
            // Calcuate the next step along the SP curve
            self.next_beta += self.increment;
            if self.next_beta > self.max_beta {
                break;
            }

            let theta_lfc = Kelvin((1.0 + self.next_beta) * self.starting_theta.unpack());
            let sh_lfc = self.starting_sh
                + self.next_beta / self.moisture_ratio.unwrap_or(1.0) / 1000.0
                    * self.starting_theta.unpack();

            let (p_lfc, t_lfc) = match potential_t_and_specific_humidity_to_pressure_and_temperature(
                theta_lfc, sh_lfc,
            ) {
                Ok((p, t)) => (p, t),
                Err(_) => continue,
            };

            let theta_e_lfc = match metfor::equiv_pot_temperature(t_lfc, t_lfc, p_lfc) {
                Some(theta_e) => theta_e,
                None => continue,
            };

            let dtheta = theta_lfc - self.starting_theta;

            let height_asl_lfc =
                match crate::interpolation::linear_interpolate(pressures, heights, p_lfc)
                    .into_option()
                {
                    Some(h) => h,
                    None => continue,
                };

            let z_fc = height_asl_lfc - self.sfc_height;

            let u = match entrained_layer_mean_wind_speed(z_fc, self.snd) {
                Ok(u) => u,
                Err(_) => continue,
            };

            let fp = metfor::pft(z_fc, p_lfc, u, dtheta, theta_lfc, self.p_sfc);

            // Vectors for recording the plume profile
            let mut pcl_t: Vec<Celsius> = vec![];
            let mut env_t: Vec<Celsius> = vec![];
            let mut pcl_p: Vec<HectoPascal> = vec![];
            let mut pcl_h: Vec<Meters> = vec![];

            // Go up dry adiabatically until we reach the lfc using the lfc potential temperature
            // calculated above.
            for ((e_p, e_h, e_t, p_t), (e_p1, e_h1, e_t1, p_t1)) in
                izip!(pressures, heights, temperatures)
                    // Filter out levels with missing values
                    .filter(|(p, h, t)| p.is_some() && t.is_some() && h.is_some())
                    // Unpack from optioned
                    .map(|(p, h, t)| (p.unpack(), h.unpack(), t.unpack()))
                    // Calculate the parcel temperature
                    .map(|(p, h, t)| {
                        (
                            p,
                            h,
                            t,
                            Celsius::from(metfor::temperature_from_pot_temp(theta_lfc, p)),
                        )
                    })
                    // Pair them up to check for crossings
                    .tuple_windows::<(_, _)>()
                    // Take while below the lfc
                    .take_while(|((p, _h, _e_t, _p_t), (_p1, _h1, _e_t1, _p_t1))| *p > p_lfc)
            {
                pcl_t.push(p_t);
                env_t.push(e_t);
                pcl_p.push(e_p);
                pcl_h.push(e_h);

                // Check for a crossing! Do linear interpolation to add the crossing position.
                if p_t > e_t && p_t1 < e_t1 {
                    let alpha = (p_t - e_t).unpack()
                        / (p_t.unpack() - p_t1.unpack() + e_t1.unpack() - e_t.unpack());
                    let p_inc = HectoPascal(e_p.unpack() + alpha * (e_p1.unpack() - e_p.unpack()));
                    let t_inc = Celsius(e_t.unpack() + alpha * (e_t1.unpack() - e_t.unpack()));
                    let h_inc = Meters(e_h.unpack() + alpha * (e_h1.unpack() - e_h.unpack()));

                    pcl_t.push(t_inc);
                    env_t.push(t_inc);
                    pcl_p.push(p_inc);
                    pcl_h.push(h_inc);
                }
            }

            // Add the lfc values.
            match crate::interpolation::linear_interpolate(pressures, temperatures, p_lfc)
                .into_option()
            {
                Some(t) => env_t.push(t),
                None => continue,
            };
            pcl_t.push(t_lfc);
            pcl_p.push(p_lfc);
            pcl_h.push(height_asl_lfc);

            // Now, keep going up moist adiabatically.
            for ((e_p, e_h, e_t, p_t), (e_p1, e_h1, e_t1, p_t1)) in
                izip!(pressures, heights, temperatures)
                    // Filter out levels with missing values
                    .filter(|(p, h, t)| p.is_some() && h.is_some() && t.is_some())
                    // Unpack from optioned
                    .map(|(p, h, t)| (p.unpack(), h.unpack(), t.unpack()))
                    // Take while below the lfc
                    .skip_while(|(p, _h, _t)| *p > p_lfc)
                    // Calculate the equivalent potential temperature of the parcel.
                    .filter_map(|(p, h, t)| {
                        metfor::temperature_from_equiv_pot_temp_saturated_and_pressure(
                            p,
                            theta_e_lfc,
                        )
                        .map(|te| (p, h, t, te))
                    })
                    // Pair them up so we can detect when the parcel profile crosses the environmental
                    // profile
                    .tuple_windows::<(_, _)>()
            {
                pcl_t.push(p_t);
                env_t.push(e_t);
                pcl_p.push(e_p);
                pcl_h.push(e_h);

                // Check for a crossing! Do linear interpolation to add the crossing position.
                if p_t > e_t && p_t1 < e_t1 {
                    let alpha = (p_t - e_t).unpack()
                        / (p_t.unpack() - p_t1.unpack() + e_t1.unpack() - e_t.unpack());
                    let p_inc = HectoPascal(e_p.unpack() + alpha * (e_p1.unpack() - e_p.unpack()));
                    let t_inc = Celsius(e_t.unpack() + alpha * (e_t1.unpack() - e_t.unpack()));
                    let h_inc = Meters(e_h.unpack() + alpha * (e_h1.unpack() - e_h.unpack()));

                    pcl_t.push(t_inc);
                    env_t.push(t_inc);
                    pcl_p.push(p_inc);
                    pcl_h.push(h_inc);
                }
            }

            // Scan the resulting arrays to find critical levels
            let (max_ib, max_dry_ib, el, max_height) = izip!(pcl_p.iter(), pcl_h.iter(), pcl_t.iter(), env_t.iter())
                // Look at two levels at a time.
                .tuple_windows::<(_, _)>()
                // Scan through and build a running integral of the buoyancy.
                .scan(
                    (0.0, 0.0, 0.0),
                    |(prev_int_buoyancy, int_buoyancy, dry_int_buoyancy): &mut (f64, f64, f64),
                     ((p0, h0, pt0, et0), (p1, h1, pt1, et1))| {
                        let dry_pcl_t0 = if *p0 > p_lfc {
                            *pt0
                        } else {
                            Celsius::from(metfor::temperature_from_pot_temp(theta_lfc, *p0))
                        };

                        let dry_pcl_t1 = if *p1 > p_lfc {
                            *pt1
                        } else {
                            Celsius::from(metfor::temperature_from_pot_temp(theta_lfc, *p1))
                        };

                        let Meters(dz) = *h1 - *h0;
                        debug_assert!(dz >= 0.0);

                        let b0 = (*pt0 - *et0) / Kelvin::from(*et0);
                        let b1 = (*pt1 - *et1) / Kelvin::from(*et1);
                        let buoyancy = (b0 + b1) * dz;

                        let db0 = (dry_pcl_t0 - *et0) / Kelvin::from(*et0);
                        let db1 = (dry_pcl_t1 - *et1) / Kelvin::from(*et1);
                        let dry_buoyancy = (db0 + db1) * dz;

                        // Update internal state, the running integrals
                        *prev_int_buoyancy = *int_buoyancy;
                        *int_buoyancy += buoyancy;
                        *dry_int_buoyancy += dry_buoyancy;
                        *dry_int_buoyancy = dry_int_buoyancy.min(*int_buoyancy);

                        // Forward the results so far.
                        Some((*int_buoyancy, *prev_int_buoyancy, *dry_int_buoyancy, h1))
                    },
                )
                // Take until we cross into negative buoyancy
                .take_while(|(_ib, p_ib, _idb, _h)| *p_ib >= 0.0)
                // Find the max_int_buoyancy, its height (the el), the max height,
                .fold(
                    (0.0f64, 0.0f64, Meters(0.0), Meters(0.0)),
                    |acc, (ib, p_ib, idb, &h)| {
                        let (mut max_ib, mut max_idb, mut el, mut max_height) = acc;

                        if ib > max_ib {
                            max_ib = ib;
                            el = h;
                        }

                        max_idb = max_idb.max(idb);

                        if ib <= 0.0 && p_ib > 0.0 {
                            max_height = h;
                        }

                        (max_ib, max_idb, el, max_height)
                    },
                );

            let el_opt: Optioned<Meters> = if el > Meters(0.0) { some(el) } else { none() };

            let max_height_opt: Optioned<Meters> = if max_height > Meters(0.0) {
                some(max_height)
            } else {
                none()
            };

            // If a max height was found, make sure LCL is below it!
            let lcl_height_op = if let Some(mxh) = max_height_opt.into_option() {
                if mxh >= height_asl_lfc {
                    some(height_asl_lfc)
                } else {
                    none()
                }
            } else {
                some(height_asl_lfc)
            };

            let (max_int_buoyancy, max_dry_int_buoyancy) = if el_opt.is_some() {
                (
                    some(JpKg(max_ib / 2.0 * -metfor::g)),
                    some(JpKg(max_dry_ib / 2.0 * -metfor::g)),
                )
            } else {
                (none(), none())
            };

            return Some(PlumeAscentAnalysis {
                fire_power: fp,
                beta: self.next_beta,
                max_int_buoyancy,
                max_dry_int_buoyancy,
                lcl_height: lcl_height_op,
                el_height: el_opt,
                max_height: max_height_opt,
                p_profile: pcl_p,
                t_profile: pcl_t,
            });
        }

        None
    }
}

#[allow(missing_docs)]
pub fn plume_ascent_analysis(
    starting_theta: Kelvin,
    starting_sh: f64,
    beta: f64,
    moisture_ratio: Option<f64>,
    sfc_height: Meters,
    p_sfc: HectoPascal,
    snd: &Sounding,
) -> Option<PlumeAscentAnalysis> {

    // Sounding profiles I'll use
    let heights = snd.height_profile();
    let pressures = snd.pressure_profile();
    let temperatures = snd.temperature_profile();

    let theta_lfc = Kelvin((1.0 + beta) * starting_theta.unpack());
    let sh_lfc = starting_sh + beta / moisture_ratio.unwrap_or(1.0) / 1000.0 * starting_theta.unpack();

    let (p_lfc, t_lfc) = potential_t_and_specific_humidity_to_pressure_and_temperature(theta_lfc, sh_lfc).ok()?; 

    let theta_e_lfc = metfor::equiv_pot_temperature(t_lfc, t_lfc, p_lfc)?; 

    let dtheta = theta_lfc - starting_theta;

    let height_asl_lfc = crate::interpolation::linear_interpolate(pressures, heights, p_lfc).into_option()?; 

    let z_fc = height_asl_lfc - sfc_height;

    let u = entrained_layer_mean_wind_speed(z_fc, snd).ok()?; 

    let fp = metfor::pft(z_fc, p_lfc, u, dtheta, theta_lfc, p_sfc);

    // Vectors for recording the plume profile
    let mut pcl_t: Vec<Celsius> = vec![];
    let mut env_t: Vec<Celsius> = vec![];
    let mut pcl_p: Vec<HectoPascal> = vec![];
    let mut pcl_h: Vec<Meters> = vec![];

    // Go up dry adiabatically until we reach the lfc using the lfc potential temperature
    // calculated above.
    for ((e_p, e_h, e_t, p_t), (e_p1, e_h1, e_t1, p_t1)) in izip!(pressures, heights, temperatures)
        // Filter out levels with missing values
        .filter(|(p, h, t)| p.is_some() && t.is_some() && h.is_some())
        // Unpack from optioned
        .map(|(p, h, t)| (p.unpack(), h.unpack(), t.unpack()))
        // Calculate the parcel temperature
        .map(|(p, h, t)| {
            (
                p,
                h,
                t,
                Celsius::from(metfor::temperature_from_pot_temp(theta_lfc, p)),
            )
        })
        // Pair them up to check for crossings
        .tuple_windows::<(_, _)>()
        // Take while below the lfc
        .take_while(|((p, _h, _e_t, _p_t), (_p1, _h1, _e_t1, _p_t1))| *p > p_lfc)
    {
        pcl_t.push(p_t);
        env_t.push(e_t);
        pcl_p.push(e_p);
        pcl_h.push(e_h);

        // Check for a crossing! Do linear interpolation to add the crossing position.
        if p_t > e_t && p_t1 < e_t1 {
            let alpha = (p_t - e_t).unpack()
                / (p_t.unpack() - p_t1.unpack() + e_t1.unpack() - e_t.unpack());
            let p_inc = HectoPascal(e_p.unpack() + alpha * (e_p1.unpack() - e_p.unpack()));
            let t_inc = Celsius(e_t.unpack() + alpha * (e_t1.unpack() - e_t.unpack()));
            let h_inc = Meters(e_h.unpack() + alpha * (e_h1.unpack() - e_h.unpack()));

            pcl_t.push(t_inc);
            env_t.push(t_inc);
            pcl_p.push(p_inc);
            pcl_h.push(h_inc);
        }
    }

    // Add the lfc values.
    match crate::interpolation::linear_interpolate(pressures, temperatures, p_lfc).into_option() {
        Some(t) => env_t.push(t),
        None => return None,
    };

    pcl_t.push(t_lfc);
    pcl_p.push(p_lfc);
    pcl_h.push(height_asl_lfc);

    // Now, keep going up moist adiabatically.
    for ((e_p, e_h, e_t, p_t), (e_p1, e_h1, e_t1, p_t1)) in izip!(pressures, heights, temperatures)
        // Filter out levels with missing values
        .filter(|(p, h, t)| p.is_some() && h.is_some() && t.is_some())
        // Unpack from optioned
        .map(|(p, h, t)| (p.unpack(), h.unpack(), t.unpack()))
        // Take while below the lfc
        .skip_while(|(p, _h, _t)| *p > p_lfc)
        // Calculate the equivalent potential temperature of the parcel.
        .filter_map(|(p, h, t)| {
            metfor::temperature_from_equiv_pot_temp_saturated_and_pressure(p, theta_e_lfc)
                .map(|te| (p, h, t, te))
        })
        // Pair them up so we can detect when the parcel profile crosses the environmental
        // profile
        .tuple_windows::<(_, _)>()
    {
        pcl_t.push(p_t);
        env_t.push(e_t);
        pcl_p.push(e_p);
        pcl_h.push(e_h);

        // Check for a crossing! Do linear interpolation to add the crossing position.
        if p_t > e_t && p_t1 < e_t1 {
            let alpha = (p_t - e_t).unpack()
                / (p_t.unpack() - p_t1.unpack() + e_t1.unpack() - e_t.unpack());
            let p_inc = HectoPascal(e_p.unpack() + alpha * (e_p1.unpack() - e_p.unpack()));
            let t_inc = Celsius(e_t.unpack() + alpha * (e_t1.unpack() - e_t.unpack()));
            let h_inc = Meters(e_h.unpack() + alpha * (e_h1.unpack() - e_h.unpack()));

            pcl_t.push(t_inc);
            env_t.push(t_inc);
            pcl_p.push(p_inc);
            pcl_h.push(h_inc);
        }
    }

    // Scan the resulting arrays to find critical levels
    let (max_ib, max_dry_ib, el, max_height) = izip!(pcl_p.iter(), pcl_h.iter(), pcl_t.iter(), env_t.iter())
        // Look at two levels at a time.
        .tuple_windows::<(_, _)>()
        // Scan through and build a running integral of the buoyancy.
        .scan(
            (0.0, 0.0, 0.0),
            |(prev_int_buoyancy, int_buoyancy, dry_int_buoyancy): &mut (f64, f64, f64),
             ((&p0, &h0, &pt0, &et0), (&p1, &h1, &pt1, &et1)) | {
                let dry_pcl_t0 = if p0 > p_lfc {
                    pt0
                } else {
                    Celsius::from(metfor::temperature_from_pot_temp(theta_lfc, p0))
                };

                let dry_pcl_t1 = if p1 > p_lfc {
                    pt1
                } else {
                    Celsius::from(metfor::temperature_from_pot_temp(theta_lfc, p1))
                };

                let Meters(dz) = h1 - h0;
                debug_assert!(dz >= 0.0);

                let b0 = (pt0 - et0) / Kelvin::from(et0);
                let b1 = (pt1 - et1) / Kelvin::from(et1);
                let buoyancy = (b0 + b1) * dz;

                let db0 = (dry_pcl_t0 - et0) / Kelvin::from(et0);
                let db1 = (dry_pcl_t1 - et1) / Kelvin::from(et1);
                let dry_buoyancy = (db0 + db1) * dz;

                // Update internal state, the running integrals
                *prev_int_buoyancy = *int_buoyancy;
                *int_buoyancy += buoyancy;
                *dry_int_buoyancy += dry_buoyancy;
                *dry_int_buoyancy = dry_int_buoyancy.min(*int_buoyancy);

                // Forward the results so far.
                Some((*int_buoyancy, *prev_int_buoyancy, *dry_int_buoyancy, h1))
            },
        )
        // Take until we cross into negative buoyancy
        .take_while(|(_ib, p_ib, _idb, _h)| *p_ib >= 0.0)
        // Find the max_int_buoyancy, its height (the el), the max height,
        .fold(
            (0.0f64, 0.0f64, Meters(0.0), Meters(0.0)),
            |acc, (ib, p_ib, idb, h)| {
                let (mut max_ib, mut max_idb, mut el, mut max_height) = acc;

                if ib > max_ib {
                    max_ib = ib;
                    el = h;
                }

                max_idb = max_idb.max(idb);

                if ib <= 0.0 && p_ib > 0.0 {
                    max_height = h;
                }

                (max_ib, max_idb, el, max_height)
            },
        );

    let el_opt: Optioned<Meters> = if el > Meters(0.0) { some(el) } else { none() };

    let max_height_opt: Optioned<Meters> = if max_height > Meters(0.0) {
        some(max_height)
    } else {
        none()
    };

    // If a max height was found, make sure LCL is below it!
    let lcl_height_op = if let Some(mxh) = max_height_opt.into_option() {
        if mxh >= height_asl_lfc {
            some(height_asl_lfc)
        } else {
            none()
        }
    } else {
        some(height_asl_lfc)
    };

    let (max_int_buoyancy, max_dry_int_buoyancy) = if el_opt.is_some() {
        (
            some(JpKg(max_ib / 2.0 * -metfor::g)),
            some(JpKg(max_dry_ib / 2.0 * -metfor::g)),
        )
    } else {
        (none(), none())
    };

    Some(PlumeAscentAnalysis {
        fire_power: fp,
        beta,
        max_int_buoyancy,
        max_dry_int_buoyancy,
        lcl_height: lcl_height_op,
        el_height: el_opt,
        max_height: max_height_opt,
        p_profile: pcl_p,
        t_profile: pcl_t,
    })
}

fn plumes_heating_iter(
    snd: &Sounding,
    moisture_ratio: Option<f64>,
) -> Result<(
    Meters,
    HectoPascal,
    Kelvin,
    f64,
    impl Iterator<Item = PlumeAscentAnalysis> + use<'_>,
)> {
    const INCREMENT: f64 = 0.000_1;
    const MAX_BETA: f64 = 0.20;
    const MAX_FP: GigaWatts = GigaWatts(1_500.0);

    let (starting_theta, starting_sh, _z_ml, _bottom_p, _p_ml) = entrained_mixed_layer(snd)?;

    let (sfc_height, p_sfc) = izip!(snd.height_profile(), snd.pressure_profile())
        .filter(|(h, p)| h.is_some() && p.is_some())
        .map(|(h, p)| (h.unpack(), p.unpack()))
        .next()
        .ok_or(AnalysisError::NotEnoughData)?;

    let anal_iter = PlumeIterator {
        next_beta: 0.0 - INCREMENT,
        max_beta: MAX_BETA,
        increment: INCREMENT,
        moisture_ratio: moisture_ratio,
        starting_theta: starting_theta,
        starting_sh: starting_sh,
        sfc_height,
        p_sfc,
        snd,
    };

    // Short fuse it for high power situations.
    let anal_iter = anal_iter.take_while(|anal| anal.fire_power <= MAX_FP);

    Ok((sfc_height, p_sfc, starting_theta, starting_sh, anal_iter))
}
