//! Indexes that are specific to a sounding, but not a particular parcel analysis of that sounding.

use sounding_base::{Profile, Sounding};

use error::*;
use interpolation::{linear_interp, linear_interpolate_sounding};
use parcel;
use parcel::Parcel;

/// The showalter index, which is like the Lifted Index except for the 850 hPa parcel.
#[inline]
pub fn showalter_index(snd: &Sounding) -> Result<f64> {
    parcel_lifted_index(snd, &parcel::pressure_parcel(snd, 850.0)?)
}

#[inline]
pub(crate) fn parcel_lifted_index(snd: &Sounding, parcel: &Parcel) -> Result<f64> {
    let parcel_profile = parcel::lift_parcel(*parcel, snd)?;

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
        total_totals = 0.0;
    }

    let mut dir_component = (d_500 - d_850).to_radians().sin();
    if dir_component < 0.0 {
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

    let (integrated_mw, _) = izip!(p_profile, dp_profile)
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
        .fold((0.0, 0.0), |(mut acc_mw, prev_p), (p, mw)| {
            let dp = prev_p - p;
            if dp > 0.0 {
                acc_mw += mw * dp;
            }

            (acc_mw, p)
        });

    Ok(integrated_mw / 9.81 / 997.0 * 100_000.0)
}
