//! Experimental sounding analysis for fire plumes utilizing the Briggs plume model for the
//! subcloud plume.

use crate::{
    error::{AnalysisError, Result},
    fire::{
        entrained_layer_mean_wind_speed, entrained_mixed_layer,
        potential_t_and_specific_humidity_to_pressure_and_temperature,
    },
    parcel::Parcel,
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
    /// The potential temperature of dry plume parcel
    pub thetas: Vec<Kelvin>,
    /// The specific humidity of the dry plume parcel
    pub specific_humidities: Vec<f64>,
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
    let (sfc_height, p_sfc, anal_iter) = plumes_heating_iter(snd, moisture_ratio)?;

    let mut fire_power: Vec<GigaWatts> = vec![];
    let mut betas: Vec<f64> = vec![];
    let mut max_int_buoyancies: Vec<Optioned<JpKg>> = vec![];
    let mut wet_ratio: Vec<Optioned<f64>> = vec![];
    let mut lcl_heights: Vec<Optioned<Meters>> = vec![];
    let mut el_heights: Vec<Optioned<Meters>> = vec![];
    let mut thetas: Vec<Kelvin> = vec![];
    let mut specific_humidities: Vec<f64> = vec![];

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
        thetas.push(anal.theta);
        specific_humidities.push(anal.specific_humidity);
    });

    Ok(BriggsPlumeHeatingAnalysis {
        fire_power,
        betas,
        max_int_buoyancies,
        wet_ratio,
        lcl_heights,
        el_heights,
        thetas,
        specific_humidities,
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
    /// The pressure portion of the pluem profile.
    pub p_profile: Vec<HectoPascal>,
    /// The temperature portion of the plume profile.
    pub t_profile: Vec<Celsius>,
    /// The starting potential Temperature
    pub theta: Kelvin,
    /// The starting specific humidity
    pub specific_humidity: f64,
}

struct PlumeIterator<'a> {
    next_beta: f64,
    max_beta: f64,
    beta_increment: f64,
    next_el: Meters,
    el_inc: Meters,
    moisture_ratio: Option<f64>,
    reached_sp_curve: bool,
    last_fire_power: GigaWatts,
    sp_theta: Kelvin,
    sp_sh: f64,
    sfc_height: Meters,
    p_sfc: HectoPascal,
    snd: &'a Sounding,
}

impl Iterator for PlumeIterator<'_> {
    type Item = PlumeAscentAnalysis;

    fn next(&mut self) -> Option<Self::Item> {

        if !self.reached_sp_curve {
            loop {
                self.next_el += self.el_inc;

                if let Ok(res_data) = entrained_mixed_layer_by_height_asl(self.next_el, self.snd) {
                    let (avg_theta, avg_sh, _h_agl, h_asl, _p_bottom, _p_top, theta_top) = res_data;

                    if h_asl < self.next_el {
                        // Reached the condensation level!
                        self.reached_sp_curve = true;
                        break;
                    }

                    if theta_top < avg_theta { continue; } // Skip this level if in absolutely unstable layer

                    let beta = theta_top.unpack() / avg_theta.unpack() - 1.0;
                    self.next_beta = beta;
                    
                    let paa = plume_ascent_analysis(
                        avg_theta,
                        avg_sh,
                        beta,
                        self.moisture_ratio,
                        self.sfc_height,
                        self.p_sfc,
                        self.snd);

                    match paa {
                        Some(p) => {
                            if p.fire_power <= self.last_fire_power {
                                continue;
                            } else {
                                self.last_fire_power = p.fire_power;
                                return Some(p);
                            }
                        },
                        None => continue,
                    };
                } else {
                    continue;
                }
            }
        }

        if self.reached_sp_curve {

            loop {
                // Calcuate the next step along the SP curve
                self.next_beta += self.beta_increment;
                if self.next_beta > self.max_beta {
                    break;
                }

                let paa = plume_ascent_analysis(
                    self.sp_theta,
                    self.sp_sh,
                    self.next_beta,
                    self.moisture_ratio,
                    self.sfc_height,
                    self.p_sfc,
                    self.snd);

                match paa {
                    Some(p) => {
                        if p.fire_power <= self.last_fire_power {
                            continue;
                        } else {
                            self.last_fire_power = p.fire_power;
                            return Some(p);
                        }
                    },
                    None => continue,
                };
            }
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

    let theta_lcl = Kelvin((1.0 + beta) * starting_theta.unpack());
    let sh_lcl = starting_sh + beta / moisture_ratio.unwrap_or(1.0) / 1000.0 * starting_theta.unpack();

    let (p_lcl, t_lcl) = potential_t_and_specific_humidity_to_pressure_and_temperature(theta_lcl, sh_lcl).ok()?; 

    let theta_e_lcl = metfor::equiv_pot_temperature(t_lcl, t_lcl, p_lcl)?; 

    let dtheta = theta_lcl - starting_theta;

    let height_asl_lcl = crate::interpolation::linear_interpolate(pressures, heights, p_lcl).into_option()?; 

    let z_fc = height_asl_lcl - sfc_height;

    let u = entrained_layer_mean_wind_speed(z_fc, snd).ok()?; 

    let fp = metfor::pft(z_fc, p_lcl, u, dtheta, theta_lcl, p_sfc);

    // Vectors for recording the plume profile
    let mut pcl_t: Vec<Celsius> = vec![];
    let mut env_t: Vec<Celsius> = vec![];
    let mut pcl_p: Vec<HectoPascal> = vec![];
    let mut pcl_h: Vec<Meters> = vec![];

    // Go up dry adiabatically until we reach the lcl using the lcl potential temperature
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
                Celsius::from(metfor::temperature_from_pot_temp(theta_lcl, p)),
            )
        })
        // Pair them up to check for crossings
        .tuple_windows::<(_, _)>()
        // Take while below the lcl
        .take_while(|((p, _h, _e_t, _p_t), (_p1, _h1, _e_t1, _p_t1))| *p > p_lcl)
    {
        assert!(e_h1 > e_h);
        assert!(height_asl_lcl > e_h);

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

            assert!(h_inc > e_h);
            assert!(h_inc < e_h1);

            if h_inc < height_asl_lcl {

                pcl_t.push(t_inc);
                env_t.push(t_inc);
                pcl_p.push(p_inc);
                pcl_h.push(h_inc);
            }
        }
    }

    // Add the lcl values.
    match crate::interpolation::linear_interpolate(pressures, temperatures, p_lcl).into_option() {
        Some(t) => env_t.push(t),
        None => return None,
    };

    pcl_t.push(t_lcl);
    pcl_p.push(p_lcl);
    pcl_h.push(height_asl_lcl);

    // Now, keep going up moist adiabatically.
    for ((e_p, e_h, e_t, p_t), (e_p1, e_h1, e_t1, p_t1)) in izip!(pressures, heights, temperatures)
        // Filter out levels with missing values
        .filter(|(p, h, t)| p.is_some() && h.is_some() && t.is_some())
        // Unpack from optioned
        .map(|(p, h, t)| (p.unpack(), h.unpack(), t.unpack()))
        // Take while below the lcl
        .skip_while(|(p, _h, _t)| *p > p_lcl)
        // Calculate the equivalent potential temperature of the parcel.
        .filter_map(|(p, h, t)| {
            metfor::temperature_from_equiv_pot_temp_saturated_and_pressure(p, theta_e_lcl)
                .map(|te| (p, h, t, te))
        })
        // Pair them up so we can detect when the parcel profile crosses the environmental
        // profile
        .tuple_windows::<(_, _)>()
    {
        assert!(e_h > height_asl_lcl);
        assert!(e_h1 > e_h);

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

            assert!(h_inc > e_h);
            assert!(h_inc < e_h1);
            assert!(h_inc > height_asl_lcl);

            pcl_t.push(t_inc);
            env_t.push(t_inc);
            pcl_p.push(p_inc);
            pcl_h.push(h_inc);
        }
    }

    // Scan the resulting arrays to find critical levels
    let (max_ib, max_dry_ib, el, _) = izip!(pcl_p.iter(), pcl_h.iter(), pcl_t.iter(), env_t.iter())
        // Look at two levels at a time.
        .tuple_windows::<(_, _)>()
        // Scan through and build a running integral of the buoyancy.
        .scan(
            (0.0, 0.0, 0.0),
            |(prev_int_buoyancy, int_buoyancy, dry_int_buoyancy): &mut (f64, f64, f64),
             ((&p0, &h0, &pt0, &et0), (&p1, &h1, &pt1, &et1)) | {
                let dry_pcl_t0 = if p0 > p_lcl {
                    pt0
                } else {
                    Celsius::from(metfor::temperature_from_pot_temp(theta_lcl, p0))
                };

                let dry_pcl_t1 = if p1 > p_lcl {
                    pt1
                } else {
                    Celsius::from(metfor::temperature_from_pot_temp(theta_lcl, p1))
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
            (0.0f64, 0.0f64, Meters(0.0), false),
            |acc, (ib, _p_ib, idb, h)| {
                let (mut max_ib, mut max_idb, mut el, mut found_first_el) = acc;

                if ib < max_ib && el > Meters(0.0) {
                    found_first_el = true;
                }

                if ib > max_ib {
                    max_ib = ib;

                    // If there is multiple EL's, only capture the first one.
                    if !found_first_el {
                        el = h;
                    }

                }

                max_idb = max_idb.max(idb);
                (max_ib, max_idb, el, found_first_el)
            },
        );

    let el_opt: Optioned<Meters> = if el > Meters(0.0) { some(el) } else { none() };
    let lcl_height_op = some(height_asl_lcl);

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
        theta: theta_lcl ,
        specific_humidity: sh_lcl,
        beta,
        max_int_buoyancy,
        max_dry_int_buoyancy,
        lcl_height: lcl_height_op,
        el_height: el_opt,
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
    impl Iterator<Item = PlumeAscentAnalysis> + use<'_>,
)> {
    const SP_INCREMENT: f64 = 0.000_1;
    const MAX_BETA: f64 = 0.20;
    const MAX_FP: GigaWatts = GigaWatts(1_500.0);
    const Z_INC: Meters = Meters(100.0);

    let (starting_theta, starting_sh, _z_ml, _bottom_p, _p_ml) = entrained_mixed_layer(snd)?;

    let (sfc_height, p_sfc) = izip!(snd.height_profile(), snd.pressure_profile())
        .filter(|(h, p)| h.is_some() && p.is_some())
        .map(|(h, p)| (h.unpack(), p.unpack()))
        .next()
        .ok_or(AnalysisError::NotEnoughData)?;

    let anal_iter = PlumeIterator {
        next_beta: 0.0 - SP_INCREMENT,
        max_beta: MAX_BETA,
        beta_increment: SP_INCREMENT,
        next_el: sfc_height,
        el_inc: Z_INC,
        moisture_ratio: moisture_ratio,
        reached_sp_curve: false,
        last_fire_power: GigaWatts(0.0),
        sp_theta: starting_theta,
        sp_sh: starting_sh,
        sfc_height,
        p_sfc,
        snd,
    };

    // Short fuse it for high power situations.
    let anal_iter = anal_iter.take_while(|anal| anal.fire_power <= MAX_FP);

    Ok((sfc_height, p_sfc, anal_iter))
}
/// This is like  step 1 from page 11 of Tory & Kepert, 2021. But instead of mixing up to the LCL,
/// we mix to a specific height or the LCL, whichever we reach first.
///
/// We also find the pressure level and height as they may be useful in later steps of the
/// algorithm.
///
/// # Returns
///
/// A tuple with the linear-height-weighted average potential temperature, linear-height-weighted
/// average specific humidity, height AGL of the top, height ASL of the top, and pressure of the
/// bottom and top of the mixing layer. The last member is the temperature at the top of the layer.
fn entrained_mixed_layer_by_height_asl(
    hgt_asl: Meters,
    snd: &Sounding,
) -> Result<(Kelvin, f64, Meters, Meters, HectoPascal, HectoPascal, Kelvin)> {

    const  MIN_DEPTH: Meters = Meters(100.0);

    let elevation: Meters = snd
        .station_info()
        .elevation()
        .ok_or(AnalysisError::NotEnoughData)?;

    let hgt_agl = hgt_asl - elevation;
    if hgt_agl < MIN_DEPTH {
        return Err(AnalysisError::NotEnoughData);
    }

    let mut ps = snd.pressure_profile();
    let mut zs = snd.height_profile();
    let mut ts = snd.temperature_profile();
    let mut dps = snd.dew_point_profile();

    // Check if the first two levels are the same level, and skip one if they are.
    if zs.len() >= 2 && zs[0] == zs[1] {
        ps = &ps[1..];
        zs = &zs[1..];
        ts = &ts[1..];
        dps = &dps[1..];
    }

    let bottom_p = ps
        .iter()
        .filter_map(|&p| p.into_option())
        .next()
        .ok_or(AnalysisError::NotEnoughData)?;

    match itertools::izip!(ps, zs, ts, dps)
        // Filter out levels with missing data
        .filter(|(p, z, t, dp)| p.is_some() && z.is_some() && t.is_some() && dp.is_some())
        // Unpack from the optional::Optioned type
        .map(|(p, z, t, dp)| (p.unpack(), z.unpack(), t.unpack(), dp.unpack()))
        // Map to height above ground.
        .map(|(p, z, t, dp)| (p, z - elevation, t, dp))
        // Add the potential temperature
        .map(|(p, h, t, dp)| (p, h, t, dp, metfor::potential_temperature(p, t)))
        // Transform to specific humidity
        .filter_map(|(p, h, t, dp, theta)| metfor::specific_humidity(dp, p).map(|q| (p, h, t, theta, q)))
        // Pair them up for integration with the trapezoid rule
        .tuple_windows::<(_, _)>()
        // Make a running integral
        .scan(
            (0.0, 0.0),
            |state, ((_p0, h0, _t0, theta0, q0), (p1, h1, t1, theta1, q1))| {
                let (sum_theta, sum_q): (&mut f64, &mut f64) = (&mut state.0, &mut state.1);

                let dh = (h1 - h0).unpack(); // meters
                debug_assert!(dh > 0.0);

                *sum_theta += (theta0.unpack() * h0.unpack() + theta1.unpack() * h1.unpack()) * dh;
                *sum_q += (q0 * h0.unpack() + q1 * h1.unpack()) * dh;

                // Divide by 2 for trapezoid rule and divide by 1/2 Z^2 for height weighting function
                let h_sq = h1.unpack() * h1.unpack();
                let avg_theta = Kelvin(*sum_theta / h_sq);
                let avg_q = *sum_q / h_sq;

                Some((p1, h1, t1, avg_theta, avg_q))
            },
        )
        // Convert into a parcel so we can lift it..
        .filter_map(|(p1, h1, t1, avg_theta, avg_q)| {
            let temperature = Celsius::from(metfor::temperature_from_pot_temp(avg_theta, bottom_p));
            let dew_point = metfor::dew_point_from_p_and_specific_humidity(bottom_p, avg_q)?;

            let pcl = Parcel {
                temperature,
                dew_point,
                pressure: bottom_p,
            };

            Some((p1, h1, t1, avg_theta, avg_q, pcl))
        })
        // Get the LCL pressure
        .filter_map(|(p1, h1, t1, avg_theta, avg_q, pcl)| {
            metfor::pressure_at_lcl(pcl.temperature, pcl.dew_point, pcl.pressure)
                .map(|lcl_pres| (p1, h1, t1, avg_theta, avg_q, lcl_pres))
        })
        // Pair up to compare layers
        .tuple_windows::<(_, _)>()
        // find where the mixed layer depth is just below the LCL or we've reached the desired height
        .find(|((p0, h0, _t0, _avg_theta0, _avg_q0, lcl_pres0), (p1, h1, _t1, _avg_theta1, _avg_q1, _))| {
            (p0 >= lcl_pres0 && p1 < lcl_pres0) || (h0 < &hgt_agl && h1 >= &hgt_agl)
        }) {
            Some(((p0, h0, t0, avg_theta0, avg_q0, lcl_pres0), (p1, h1, t1, avg_theta1, avg_q1, _))) => {
                if p0 >= lcl_pres0 && p1 < lcl_pres0 && h0 + elevation < hgt_asl {
                    // Hit the LCL
                    let theta_top = metfor::potential_temperature(p0, t0);
                    Ok((avg_theta0, avg_q0, h0, h0 + elevation, bottom_p, p0, theta_top))
                } else {
                    assert!(h0 <= hgt_agl && h1 >= hgt_agl);

                    let dh = (hgt_agl - h0).unpack();
                    let run = (h1 - h0).unpack();
                    let rise_theta = (avg_theta1 - avg_theta0).unpack();
                    let rise_q = avg_q1 - avg_q0;
                    let rise_p = (p1 - p0).unpack();
                    let rise_t = (t1 - t0).unpack();
                    let theta = Kelvin(avg_theta0.unpack() + rise_theta / run * dh);
                    

                    let q = avg_q0 + rise_q / run * dh;
                    let p = HectoPascal(p0.unpack() + rise_p / run * dh);
                    let t = Celsius(t0.unpack() + rise_t / run * dh);

                    let theta_top = metfor::potential_temperature(p, t);
                    Ok((theta, q, hgt_agl, hgt_asl, bottom_p, p, theta_top))
                }
            },
            None => Err(AnalysisError::FailedPrerequisite),
        }
}

