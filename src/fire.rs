//! Module for analysis related to fire weather and wildfire plumes.

use crate::{
    error::{AnalysisError, Result},
    parcel::Parcel,
    sounding::Sounding,
};
use itertools::{izip, Itertools};
use metfor::{
    Celsius, CelsiusDiff, GigaWatts, HectoPascal, Kelvin, KelvinDiff, Meters, MetersPSec, Quantity,
    WindSpdDir,
};

/// The Hot-Dry-Windy index
///
/// # References
///
/// Srock AF, Charney JJ, Potter BE, Goodrick SL. The Hot-Dry-Windy Index: A New Fire Weather Index.
///     Atmosphere. 2018; 9(7):279. https://doi.org/10.3390/atmos9070279
///
/// McDonald JM, Srock AF, Charney JJ. Development and Application of a Hot-Dry-Windy Index (HDW)
/// Climatology. Atmosphere. 2018; 9: 285. https://doi.org/10.3390/atmos9070285
#[inline]
pub fn hot_dry_windy(snd: &Sounding) -> Result<f64> {
    let elevation = if let Some(sfc_h) = snd.station_info().elevation().into_option() {
        sfc_h
    } else if let Some(lowest_h) = snd
        .height_profile()
        .iter()
        .filter_map(|optd| optd.into_option())
        .next()
    {
        lowest_h
    } else {
        return Err(AnalysisError::NotEnoughData);
    };

    let sfc_pressure: HectoPascal = snd
        .pressure_profile()
        .iter()
        .filter_map(|p| p.into_option())
        .next()
        .ok_or(AnalysisError::NotEnoughData)?;

    let p_profile = snd.pressure_profile();
    let h_profile = snd.height_profile();
    let t_profile = snd.temperature_profile();
    let dp_profile = snd.dew_point_profile();
    let ws_profile = snd.wind_profile();

    let (vpd, ws) = izip!(p_profile, h_profile, t_profile, dp_profile, ws_profile)
        // Remove rows with missing data
        .filter(|(p, h, t, dp, ws)| {
            p.is_some() && h.is_some() && t.is_some() && dp.is_some() && ws.is_some()
        })
        // Unpack from the Optioned type
        .map(|(p, h, t, dp, ws)| {
            (
                p.unpack(),
                h.unpack(),
                t.unpack(),
                dp.unpack(),
                ws.unpack().speed,
            )
        })
        // Only look up to 500 m above AGL
        .take_while(|(_, h, _, _, _)| *h <= elevation + Meters(500.0))
        // Calculate the surface adjusted temperature - for the surface adjusted VPD
        .map(|(p, _, t, dp, ws)| {
            (
                p,
                metfor::temperature_from_pot_temp(
                    metfor::potential_temperature(p, t),
                    sfc_pressure,
                ),
                dp,
                ws,
            )
        })
        // Calculate the surface adjusted dew point - for the surface adjusted VPD
        .filter_map(|(p, t, dp, ws)| {
            metfor::specific_humidity(dp, p)
                .and_then(|q| metfor::dew_point_from_p_and_specific_humidity(sfc_pressure, q))
                .map(|dp| (t, dp, ws))
        })
        // Convert t and dp to VPD
        .filter_map(|(t, dp, ws)| {
            metfor::vapor_pressure_water(t)
                .and_then(|sat_vap| metfor::vapor_pressure_water(dp).map(|vap| (sat_vap - vap, ws)))
        })
        // Convert knots to m/s and unpack all values from their Quantity types
        .map(|(vpd, ws)| (vpd.unpack(), MetersPSec::from(ws).unpack()))
        // Choose the max.
        .fold((0.0, 0.0), |(vpd_max, ws_max), (vpd, ws)| {
            (vpd.max(vpd_max), ws.max(ws_max))
        });

    Ok(vpd * ws)
}

/// A collection of parameters associated with a Pyrocumulonimbus Firepower Threshold (PFT)
/// analysis. See [pft] or [pft_analysis] for details.
#[derive(Clone, Debug)]
pub struct PFTAnalysis {
    /// The fire power required to cause a pyrocumulonimbus in Gigawatts. This is not the power per
    /// unit area, but the total power.
    pub pft: GigaWatts,
    /// The height weighted average potential temperature in the mixed layer.
    pub theta_ml: Kelvin,
    /// The height weighted average specific humidity in the mixed layer.
    pub q_ml: f64,
    /// The pressure at the top of the mixed layer.
    pub p_top_ml: HectoPascal,
    /// The pressure at the bottom of the pyroCb, where free convection starts.
    pub p_fc: HectoPascal,
    /// Coordinates along the SP-curve suitable for plotting on a skew-t log-p chart.
    pub sp_curve: Vec<(HectoPascal, Celsius)>,
    /// The minimum equivalent potential temperature of the plume element required to initiate a
    /// pyrocumulonimbus cloud.
    pub theta_e_fc: Kelvin,
}

/// Calculate the Pyrocumulonimbus Firepower Threshold (PFT).
///
/// The first reference below (Tory & Kepert, 2021) has most of the details about how to calculate
/// the PFT, the other paper (Tory et. al, 2018) outlines the model in general.
///
/// # References
///
/// Tory, K. J., & Kepert, J. D. (2021). **Pyrocumulonimbus Firepower Threshold: Assessing the
///     Atmospheric Potential for pyroCb**, Weather and Forecasting, 36(2), 439-456. Retrieved Jun 2,
///     2021, from https://journals.ametsoc.org/view/journals/wefo/36/2/WAF-D-20-0027.1.xml
///
/// Tory, K. J., Thurston, W., & Kepert, J. D. (2018). **Thermodynamics of Pyrocumulus: A Conceptual
///     Study**, Monthly Weather Review, 146(8), 2579-2598. Retrieved Jun 2, 2021, from
///     https://journals.ametsoc.org/view/journals/mwre/146/8/mwr-d-17-0377.1.xml
///
/// # Arguments
///  - snd is the sounding on which to do this calculation.
///  - moisture_ratio is the amount of heating required to add 1 g/kg to the specific humidity.
///
pub fn pft_analysis(snd: &Sounding, moisture_ratio: f64) -> Result<PFTAnalysis> {
    let (theta_ml, q_ml, _z_ml, p_sfc, p_top_ml) = entrained_mixed_layer(snd)?;

    let (z_fc, p_fc, theta_fc, d_theta_fc, theta_e_fc, sp_curve) =
        free_convection_level(snd, moisture_ratio, theta_ml, q_ml)?;

    let u_ml: MetersPSec = entrained_layer_mean_wind_speed(z_fc, snd)?;

    let pft = metfor::pft(z_fc, p_fc, u_ml, d_theta_fc, theta_fc, p_sfc);

    Ok(PFTAnalysis {
        pft,
        theta_ml,
        q_ml,
        p_top_ml,
        p_fc,
        sp_curve,
        theta_e_fc,
    })
}

/// Calculate the Pyrocumulonimbus Firepower Threshold (PFT).
///
/// The first reference below (Tory & Kepert, 2021) has most of the details about how to calculate
/// the PFT, the other paper (Tory et. al, 2018) outlines the model in general.
///
/// # References
///
/// Tory, K. J., & Kepert, J. D. (2021). **Pyrocumulonimbus Firepower Threshold: Assessing the
///     Atmospheric Potential for pyroCb**, Weather and Forecasting, 36(2), 439-456. Retrieved Jun 2,
///     2021, from https://journals.ametsoc.org/view/journals/wefo/36/2/WAF-D-20-0027.1.xml
///
/// Tory, K. J., Thurston, W., & Kepert, J. D. (2018). **Thermodynamics of Pyrocumulus: A Conceptual
///     Study**, Monthly Weather Review, 146(8), 2579-2598. Retrieved Jun 2, 2021, from
///     https://journals.ametsoc.org/view/journals/mwre/146/8/mwr-d-17-0377.1.xml
///
/// # Arguments
///  - snd is the sounding on which to do this calculation.
///  - moisture_ratio is the amount of heating required to add 1 g/kg to the specific humidity.
///
pub fn pft(snd: &Sounding, moisture_ratio: f64) -> Result<GigaWatts> {
    let (theta_ml, q_ml, _z_ml, p_sfc, _p_top_ml) = entrained_mixed_layer(snd)?;

    let (z_fc, p_fc, theta_fc, d_theta_fc, _theta_e, _sp_curve) =
        free_convection_level(snd, moisture_ratio, theta_ml, q_ml)?;

    let u_ml: MetersPSec = entrained_layer_mean_wind_speed(z_fc, snd)?;

    Ok(metfor::pft(z_fc, p_fc, u_ml, d_theta_fc, theta_fc, p_sfc))
}

/// This is step 1 from page 11 of Tory & Kepert, 2021.
///
/// We also find the pressure level and height as they may be useful in later steps of the
/// algorithm.
///
/// # Returns
///
/// A tuple with the linear-height-weighted average potential temperature, linear-height-weighted
/// average specific humidity, height AGL of the top, and pressure of the top of the mixing layer.
fn entrained_mixed_layer(
    snd: &Sounding,
) -> Result<(Kelvin, f64, Meters, HectoPascal, HectoPascal)> {
    let elevation: Meters = snd
        .station_info()
        .elevation()
        .ok_or(AnalysisError::NotEnoughData)?;

    let ps = snd.pressure_profile();
    let zs = snd.height_profile();
    let ts = snd.temperature_profile();
    let dps = snd.dew_point_profile();

    let bottom_p = ps
        .iter()
        .filter_map(|&p| p.into_option())
        .next()
        .ok_or(AnalysisError::NotEnoughData)?;

    let max_p = bottom_p + HectoPascal(50.0);

    match itertools::izip!(ps, zs, ts, dps)
        // Filter out levels with missing data
        .filter(|(p, z, t, dp)| p.is_some() && z.is_some() && t.is_some() && dp.is_some())
        // Unpack from the optional::Optioned type
        .map(|(p, z, t, dp)| (p.unpack(), z.unpack(), t.unpack(), dp.unpack()))
        // Map to height above ground.
        .map(|(p, z, t, dp)| (p, z - elevation, t, dp))
        // Transform to potential temperature
        .map(|(p, h, t, dp)| (p, h, metfor::potential_temperature(p, t), dp))
        // Transform to specific humidity
        .filter_map(|(p, h, theta, dp)| metfor::specific_humidity(dp, p).map(|q| (p, h, theta, q)))
        // Pair them up for integration with the trapezoid rule
        .tuple_windows::<(_, _)>()
        // Make a running integral
        .scan(
            (0.0, 0.0),
            |state, ((_p0, h0, theta0, q0), (p1, h1, theta1, q1))| {
                let (sum_theta, sum_q): (&mut f64, &mut f64) = (&mut state.0, &mut state.1);

                let dh = (h1 - h0).unpack(); // meters
                *sum_theta += (theta0.unpack() * h0.unpack() + theta1.unpack() * h1.unpack()) * dh;
                *sum_q += (q0 * h0.unpack() + q1 * h1.unpack()) * dh;

                // Divide by 2 for trapezoid rule and divide by 1/2 Z^2 for height weighting function
                let h_sq = h1.unpack() * h1.unpack();
                let avg_theta = Kelvin(*sum_theta / h_sq);
                let avg_q = *sum_q / h_sq;

                Some((p1, h1, avg_theta, avg_q))
            },
        )
        // Don't even start looking until we have at least a minimum thickness mixed layer.
        .filter(|(p1, _h1, _avg_theta, _avg_q)| *p1 <= max_p)
        // Convert into a parcel so we can lift it..
        .filter_map(|(p1, h1, avg_theta, avg_q)| {
            let temperature = Celsius::from(metfor::temperature_from_pot_temp(avg_theta, bottom_p));
            let dew_point = metfor::dew_point_from_p_and_specific_humidity(bottom_p, avg_q)?;

            let pcl = Parcel {
                temperature,
                dew_point,
                pressure: bottom_p,
            };

            Some((p1, h1, avg_theta, avg_q, pcl))
        })
        // Get the LCL pressure
        .filter_map(|(p1, h1, avg_theta, avg_q, pcl)| {
            metfor::pressure_at_lcl(pcl.temperature, pcl.dew_point, pcl.pressure)
                .map(|lcl_pres| (p1, h1, avg_theta, avg_q, lcl_pres))
        })
        // Pair up to compare layers
        .tuple_windows::<(_, _)>()
        // find where the mixed layer depth is just below the LCL
        .find(
            |(
                (p0, _h0, _avg_theta0, _avg_q0, lcl_pres0),
                (p1, _h1, _avg_theta1, _avg_q1, lcl_pres1),
            )| p0 >= lcl_pres0 && p1 < lcl_pres1,
        ) {
        Some(((p0, h0, avg_theta0, avg_q0, _), _)) => Ok((avg_theta0, avg_q0, h0, bottom_p, p0)),
        None => Err(AnalysisError::FailedPrerequisite),
    }
}

/// This is steps 3 and 4 from page 11 of Tory & Kepert, 2021.
///
/// We skipped step 2 to get here. Basically, use the max virtual temperature above z_ml
/// (or p_top_ml) to find the theta-e that is above that temperature by the threshold amount AND
/// tops out on the sounding at at least -20C or colder. Once we know what theta-e value we need,
/// we can find the level where that intersects the S-P curve.
///
/// Once we have our point on the S-P curve, zfc and d_theta_fc are easy to calculate.
fn free_convection_level(
    snd: &Sounding,
    moisture_ratio: f64,
    theta_ml: Kelvin,
    q_ml: f64,
) -> Result<(
    Meters,
    HectoPascal,
    Kelvin,
    KelvinDiff,
    Kelvin,
    Vec<(HectoPascal, Celsius)>,
)> {
    const SP_BETA_MAX: f64 = 0.25;
    const DELTA_BETA: f64 = 0.001;

    let apply_beta = move |beta| {
        let theta_sp = Kelvin((1.0 + beta) * theta_ml.unpack());
        let q_sp = q_ml + beta / moisture_ratio / 1000.0 * theta_ml.unpack();

        (theta_sp, q_sp)
    };

    let mut sp_curve: Vec<(HectoPascal, Celsius)> = Vec::with_capacity(100);

    let mut low_beta: f64 = f64::NAN;
    let mut high_beta: f64 = f64::NAN;
    let mut beta = 0.0;
    let mut i = 0;
    while beta <= SP_BETA_MAX {
        let (theta_sp, q_sp) = apply_beta(beta);

        let skew_t_coords =
            potential_t_and_specific_humidity_to_pressure_and_temperature(theta_sp, q_sp)?;

        sp_curve.push(skew_t_coords);

        let (p, t) = skew_t_coords;

        let sp_theta_e = match metfor::equiv_pot_temperature(t, t, p) {
            Some(val) => val,
            None => continue,
        };

        let meets_theta_e_requirements = is_free_convecting(snd, p, sp_theta_e)?;

        if !meets_theta_e_requirements && high_beta.is_nan() {
            low_beta = beta;
        } else if meets_theta_e_requirements && high_beta.is_nan() {
            high_beta = beta;
        }

        if !high_beta.is_nan() && i >= 100 {
            break;
        }

        beta += DELTA_BETA;
        i += 1;
    }

    if high_beta.is_nan() || low_beta.is_nan() {
        return Err(AnalysisError::FailedPrerequisite);
    }

    let beta = metfor::find_root(
        &|b| {
            let (theta_sp, q_sp) = apply_beta(b);
            let (p, t) =
                potential_t_and_specific_humidity_to_pressure_and_temperature(theta_sp, q_sp)
                    .ok()?;
            let sp_theta_e = metfor::equiv_pot_temperature(t, t, p)?;

            Some(
                (min_temperature_diff_to_max_cloud_top_temperature(snd, p, sp_theta_e).ok()?
                    - TEMPERATURE_BUFFER)
                    .unpack(),
            )
        },
        low_beta,
        high_beta,
    )
    .ok_or(AnalysisError::FailedPrerequisite)?;

    let (theta_sp, q_sp) = apply_beta(beta);
    let (p_lfc, t_lfc) =
        potential_t_and_specific_humidity_to_pressure_and_temperature(theta_sp, q_sp)?;
    let theta_e_lfc = metfor::equiv_pot_temperature(t_lfc, t_lfc, p_lfc)
        .ok_or(AnalysisError::FailedPrerequisite)?;

    let dtheta = theta_sp - theta_ml;

    let heights = snd.height_profile();
    let pressures = snd.pressure_profile();
    let height_asl_lfc = crate::interpolation::linear_interpolate(pressures, heights, p_lfc)
        .into_option()
        .ok_or(AnalysisError::InterpolationError)?;
    let sfc_height = heights
        .iter()
        .filter_map(|h| h.into_option())
        .next()
        .ok_or(AnalysisError::NotEnoughData)?;

    Ok((
        height_asl_lfc - sfc_height,
        p_lfc,
        theta_sp,
        dtheta,
        theta_e_lfc,
        sp_curve,
    ))
}

fn potential_t_and_specific_humidity_to_pressure_and_temperature(
    theta: Kelvin,
    specific_humidity: f64,
) -> Result<(HectoPascal, Celsius)> {
    let diff = |p_hpa: f64| -> Option<f64> {
        let t = metfor::temperature_from_pot_temp(theta, HectoPascal(p_hpa));
        let dp =
            metfor::dew_point_from_p_and_specific_humidity(HectoPascal(p_hpa), specific_humidity)?;

        Some((t - dp).unpack())
    };

    let p = metfor::find_root(
        &diff,
        HectoPascal(1080.0).unpack(),
        HectoPascal(100.0).unpack(),
    )
    .map(HectoPascal)
    .ok_or(AnalysisError::FailedPrerequisite)?;

    let t = Celsius::from(metfor::temperature_from_pot_temp(theta, p));

    Ok((p, t))
}

/*
/// This implements equations 14-17 from Tory et al, 2018.
fn potential_t_and_specific_humidity_to_pressure_and_temperature(
    theta: Kelvin,
    specific_humidity: f64,
) -> (HectoPascal, Celsius) {
    let theta_pl = Kelvin::from(theta).unpack();
    let t_pl = theta_pl;
    let p0 = HectoPascal(1000.0).unpack();
    let ps = p0;
    let mw_pl = metfor::mixing_ratio_from_specific_humidity(specific_humidity);

    let e_pl = ps / ((1.0 - metfor::epsilon) + metfor::epsilon / specific_humidity);
    let t_sp = 2840.0 / (3.5 * t_pl.ln() - e_pl.ln() - 4.805) + 55.0;

    let k = metfor::cpd / metfor::Rd * (1.0 + mw_pl * 1870.0 / 1005.7)
        / (1.0 + mw_pl / metfor::epsilon);

    let p_sp = p0 * (t_sp / theta_pl).powf(k);

    let t_sp = Celsius::from(Kelvin(t_sp));
    let p_sp = HectoPascal(p_sp);

    (p_sp, t_sp)
}
*/

/// Get the minimum difference between this theta-e and all those up to -20C.
fn min_temperature_diff_to_max_cloud_top_temperature(
    snd: &Sounding,
    starting_pressure: HectoPascal,
    starting_theta_e: Kelvin,
) -> Result<CelsiusDiff> {
    const MAX_PLUME_TOP_T: Celsius = Celsius(-20.0);

    let pp = snd.pressure_profile();
    let tp = snd.temperature_profile();

    let min_diff: CelsiusDiff = izip!(pp, tp)
        // Skip levels with missing data.
        .filter(|(p, t)| p.is_some() && t.is_some())
        //Unpack optioned values.
        .map(|(p, t)| (p.unpack(), t.unpack()))
        // Skip below the top of the mixed layer.
        .skip_while(|(p, _t)| *p > starting_pressure)
        // Get the parcel temperature
        .filter_map(|(p, t)| {
            metfor::temperature_from_equiv_pot_temp_saturated_and_pressure(p, starting_theta_e)
                .map(|pcl_t| (t, pcl_t))
        })
        // Ensure we get a least 2 levels
        .enumerate()
        .take_while(|(i, (_t, pcl_t))| *i < 2 || *pcl_t >= MAX_PLUME_TOP_T)
        // Find the minimum difference
        .fold(CelsiusDiff(500.0), |min_diff, (_i, (t, pcl_t))| {
            (pcl_t - t).min(min_diff)
        });

    if min_diff == CelsiusDiff(500.0) {
        Err(AnalysisError::FailedPrerequisite)
    } else {
        Ok(min_diff)
    }
}

/// Find the maximum equivalent_potential_temperature above the mixed layer that interesects the
/// sounding at or above the -20C level.
const TEMPERATURE_BUFFER: CelsiusDiff = CelsiusDiff(0.5);
fn is_free_convecting(
    snd: &Sounding,
    starting_pressure: HectoPascal,
    starting_theta_e: Kelvin,
) -> Result<bool> {
    Ok(
        min_temperature_diff_to_max_cloud_top_temperature(
            snd,
            starting_pressure,
            starting_theta_e,
        )? >= TEMPERATURE_BUFFER,
    )
}

/// Find the mean wind in the entrainment layer.
///
/// Remember, zfc is the level AGL.
fn entrained_layer_mean_wind_speed(zfc: Meters, snd: &Sounding) -> Result<MetersPSec> {
    let layer = crate::layer_agl(snd, zfc)?;

    let mean_wind = crate::wind::mean_wind(&layer, snd)?;
    let mean_wind = WindSpdDir::from(mean_wind).speed;

    Ok(mean_wind)
}
