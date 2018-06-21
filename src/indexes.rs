//! Indexes that are specific to a sounding, but not a particular parcel analysis of that sounding.

use sounding_base::{Profile, Sounding};

use error::*;
use interpolation::{linear_interp, linear_interpolate_sounding};
use parcel;
use parcel::ParcelProfile;

/// The showalter index, which is like the Lifted Index except for the 850 hPa parcel.
// TODO: This needs validated, it is very different than bufkit files.
#[inline]
pub fn showalter_index(snd: &Sounding) -> Result<f64> {
    let parcel = parcel::pressure_parcel(snd, 850.0)?;
    let profile = parcel::lift_parcel(parcel, snd)?;
    parcel_lifted_index(&profile)
}

#[inline]
pub(crate) fn parcel_lifted_index(parcel_profile: &ParcelProfile) -> Result<f64> {
    izip!(
        &parcel_profile.pressure,
        &parcel_profile.parcel_t,
        &parcel_profile.environment_t
    ).take_while(|&(p, _, _)| *p > 500.0)
        .last()
        .and_then(|(bottom_p, bottom_p_t, bottom_e_t)| {
            izip!(
                &parcel_profile.pressure,
                &parcel_profile.parcel_t,
                &parcel_profile.environment_t
            ).skip_while(|&(p, _, _)| *p > 500.0)
                .nth(0)
                .and_then(|(top_p, top_p_t, top_e_t)| {
                    let e_t = linear_interp(500.0, *bottom_p, *top_p, *bottom_e_t, *top_e_t);
                    let p_t = linear_interp(500.0, *bottom_p, *top_p, *bottom_p_t, *top_p_t);
                    Some(e_t - p_t)
                })
        })
        .ok_or(AnalysisError::MissingValue)
}

/// The Total Totals index
pub fn total_totals(snd: &Sounding) -> Result<f64> {
    let h5 = linear_interpolate_sounding(snd, 500.0)?;
    let h85 = linear_interpolate_sounding(snd, 850.0)?;
    let t_500 = h5.temperature.ok_or(AnalysisError::MissingValue)?;
    let t_850 = h85.temperature.ok_or(AnalysisError::MissingValue)?;
    let td_850 = h85.dew_point.ok_or(AnalysisError::MissingValue)?;

    let cross_totals = td_850 - t_500;
    let vertical_totals = t_850 - t_500;

    Ok(cross_totals + vertical_totals)
}

/// The SWeT (Severe Weather Threat) index
#[inline]
pub fn swet(snd: &Sounding) -> Result<f64> {
    let h5 = linear_interpolate_sounding(snd, 500.0)?;
    let h85 = linear_interpolate_sounding(snd, 850.0)?;

    let mut td_850 = h85.dew_point.ok_or(AnalysisError::MissingValue)?;
    if td_850 < 0.0 {
        td_850 = 0.0;
    }

    let v_850 = h85.speed.ok_or(AnalysisError::MissingValue)?;
    let d_850 = h85.direction.ok_or(AnalysisError::MissingValue)?;

    let v_500 = h5.speed.ok_or(AnalysisError::MissingValue)?;
    let d_500 = h5.direction.ok_or(AnalysisError::MissingValue)?;

    let mut total_totals = total_totals(snd)?;
    if total_totals < 49.0 {
        total_totals = 49.0;
    }

    let mut dir_component = (d_500 - d_850).to_radians().sin();
    if dir_component < 0.0 || (d_850 >= 130.0 && d_850 <= 250.0)
        || (d_500 >= 210.0 && d_500 <= 310.0) || (d_500 - d_850) >= 0.0
        || (v_850 >= 15.0 && v_500 >= 15.0)
    {
        dir_component = 0.0;
    } else {
        dir_component = 125.0 * (dir_component + 0.2);
    }

    Ok(12.0 * td_850 + 20.0 * (total_totals - 49.0) + 2.0 * v_850 + v_500 + dir_component)
}

/// The K-index
#[inline]
pub fn kindex(snd: &Sounding) -> Result<f64> {
    let h5 = linear_interpolate_sounding(snd, 500.0)?;
    let h7 = linear_interpolate_sounding(snd, 700.0)?;
    let h85 = linear_interpolate_sounding(snd, 850.0)?;

    let t5 = h5.temperature.ok_or(AnalysisError::MissingValue)?;
    let t7 = h7.temperature.ok_or(AnalysisError::MissingValue)?;
    let t85 = h85.temperature.ok_or(AnalysisError::MissingValue)?;

    let td7 = h7.dew_point.ok_or(AnalysisError::MissingValue)?;
    let td85 = h85.dew_point.ok_or(AnalysisError::MissingValue)?;

    Ok(t85 - t5 + td85 - t7 + td7)
}

/// Precipitable water (mm)
#[inline]
pub fn precipitable_water(snd: &Sounding) -> Result<f64> {
    let p_profile = snd.get_profile(Profile::Pressure);
    let dp_profile = snd.get_profile(Profile::DewPoint);

    let (integrated_mw, _, _) = izip!(p_profile, dp_profile)
        .filter_map(|pair| {
            let (p, dp) = pair;
            if p.is_some() && dp.is_some() {
                Some((p.unpack(), dp.unpack()))
            } else {
                None
            }
        })
        .filter_map(|(p, dp)| {
            ::metfor::mixing_ratio(dp, p)
                .ok()
                .and_then(|mw| Some((p, mw)))
        })
        .fold((0.0, 0.0, 0.0), |(mut acc_mw, prev_p, prev_mw), (p, mw)| {
            let dp = prev_p - p;
            if dp > 0.0 {
                acc_mw += (mw + prev_mw) * dp;
            }

            (acc_mw, p, mw)
        });

    Ok(integrated_mw / 9.81 / 997.0 * 100_000.0 / 2.0)
}

/// The Haines index for fire weather
#[inline]
pub fn haines(snd: &Sounding) -> Result<f64> {
    snd.get_station_info()
        .elevation()
        .into_option()
        .and_then(|elev_meters| {
            let (stability_term, stability_cutoffs, moisture_term, moisture_cutoffs) =
                if elev_meters <= 304.8 {
                    // Low elevation Haines
                    let level1 = linear_interpolate_sounding(snd, 950.0).ok()?;
                    let level2 = linear_interpolate_sounding(snd, 850.0).ok()?;

                    let t_low = level1.temperature.ok_or(AnalysisError::MissingValue).ok()?;
                    let t_hi = level2.temperature.ok_or(AnalysisError::MissingValue).ok()?;
                    let dp_hi = level2.dew_point.ok_or(AnalysisError::MissingValue).ok()?;

                    let stability_term = t_low - t_hi;
                    const STABILITY_CUTOFFS: (f64, f64) = (4.0, 7.0);

                    let moisture_term = t_hi - dp_hi;
                    const MOISTURE_CUTOFFS: (f64, f64) = (6.0, 9.0);

                    (
                        stability_term,
                        STABILITY_CUTOFFS,
                        moisture_term,
                        MOISTURE_CUTOFFS,
                    )
                } else if elev_meters <= 914.4 {
                    // Mid elevation Haines
                    let level1 = linear_interpolate_sounding(snd, 850.0).ok()?;
                    let level2 = linear_interpolate_sounding(snd, 700.0).ok()?;

                    let t_low = level1.temperature.ok_or(AnalysisError::MissingValue).ok()?;
                    let t_hi = level2.temperature.ok_or(AnalysisError::MissingValue).ok()?;
                    let dp_low = level2.dew_point.ok_or(AnalysisError::MissingValue).ok()?;

                    let stability_term = t_low - t_hi;
                    const STABILITY_CUTOFFS: (f64, f64) = (6.0, 10.0);

                    let moisture_term = t_low - dp_low;
                    const MOISTURE_CUTOFFS: (f64, f64) = (6.0, 12.0);

                    (
                        stability_term,
                        STABILITY_CUTOFFS,
                        moisture_term,
                        MOISTURE_CUTOFFS,
                    )
                } else {
                    // High elevation Haines
                    let level1 = linear_interpolate_sounding(snd, 700.0).ok()?;
                    let level2 = linear_interpolate_sounding(snd, 500.0).ok()?;

                    let t_low = level1.temperature.ok_or(AnalysisError::MissingValue).ok()?;
                    let t_hi = level2.temperature.ok_or(AnalysisError::MissingValue).ok()?;
                    let dp_low = level2.dew_point.ok_or(AnalysisError::MissingValue).ok()?;

                    let stability_term = t_low - t_hi;
                    const STABILITY_CUTOFFS: (f64, f64) = (18.0, 21.0);

                    let moisture_term = t_low - dp_low;
                    const MOISTURE_CUTOFFS: (f64, f64) = (15.0, 20.0);

                    (
                        stability_term,
                        STABILITY_CUTOFFS,
                        moisture_term,
                        MOISTURE_CUTOFFS,
                    )
                };

            let stability_part = if stability_term < stability_cutoffs.0 {
                1.0
            } else if stability_term < stability_cutoffs.1 {
                2.0
            } else {
                3.0
            };

            let moisture_part = if moisture_term < moisture_cutoffs.0 {
                1.0
            } else if moisture_term < moisture_cutoffs.1 {
                2.0
            } else {
                3.0
            };

            let haines_index = stability_part + moisture_part;

            debug_assert!(haines_index <= 6.0);
            Some(haines_index)
        })
        .ok_or(AnalysisError::MissingValue)
}
