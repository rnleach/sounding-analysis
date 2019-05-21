//! Indexes that are specific to a sounding, but not a particular parcel analysis of that sounding.

use crate::{
    error::{AnalysisError, Result},
    interpolation::linear_interpolate_sounding,
};
use itertools::izip;
use metfor::{
    vapor_pressure_liquid_water, Celsius, CelsiusDiff, HectoPascal, Knots, Meters, MetersPSec, Mm,
    Quantity, WindSpdDir,
};
use sounding_base::Sounding;

/// The Total Totals index
#[inline]
pub fn total_totals(snd: &Sounding) -> Result<f64> {
    let h5 = linear_interpolate_sounding(snd, HectoPascal(500.0))?;
    let h85 = linear_interpolate_sounding(snd, HectoPascal(850.0))?;
    let t_500 = h5.temperature.ok_or(AnalysisError::MissingValue)?;
    let t_850 = h85.temperature.ok_or(AnalysisError::MissingValue)?;
    let td_850 = h85.dew_point.ok_or(AnalysisError::MissingValue)?;

    let cross_totals = td_850 - t_500;
    let vertical_totals = t_850 - t_500;

    let CelsiusDiff(tt) = cross_totals + vertical_totals;

    Ok(tt)
}

/// The SWeT (Severe Weather Threat) index
#[inline]
pub fn swet(snd: &Sounding) -> Result<f64> {
    let h5 = linear_interpolate_sounding(snd, HectoPascal(500.0))?;
    let h85 = linear_interpolate_sounding(snd, HectoPascal(850.0))?;

    let td_850 = h85
        .dew_point
        .map_t(|Celsius(dp)| dp)
        .map_t(|dp| if dp < 0.0 { 0.0 } else { dp })
        .ok_or(AnalysisError::MissingValue)?;

    let WindSpdDir {
        speed: v_850,
        direction: d_850,
    } = h85.wind.ok_or(AnalysisError::MissingValue)?;
    let WindSpdDir {
        speed: v_500,
        direction: d_500,
    } = h5.wind.ok_or(AnalysisError::MissingValue)?;

    let mut total_totals = total_totals(snd)?;
    if total_totals < 49.0 {
        total_totals = 49.0;
    }

    let mut dir_component = (d_500 - d_850).to_radians().sin();
    if dir_component < 0.0
        || (d_850 >= 130.0 && d_850 <= 250.0)
        || (d_500 >= 210.0 && d_500 <= 310.0)
        || (d_500 - d_850) >= 0.0
        || (v_850 >= Knots(15.0) && v_500 >= Knots(15.0))
    {
        dir_component = 0.0;
    } else {
        dir_component = 125.0 * (dir_component + 0.2);
    }

    Ok(12.0 * td_850
        + 20.0 * (total_totals - 49.0)
        + 2.0 * v_850.unpack()
        + v_500.unpack()
        + dir_component)
}

/// The K-index
#[inline]
pub fn kindex(snd: &Sounding) -> Result<Celsius> {
    let h5 = linear_interpolate_sounding(snd, HectoPascal(500.0))?;
    let h7 = linear_interpolate_sounding(snd, HectoPascal(700.0))?;
    let h85 = linear_interpolate_sounding(snd, HectoPascal(850.0))?;

    let t5 = h5.temperature.ok_or(AnalysisError::MissingValue)?;
    let t7 = h7.temperature.ok_or(AnalysisError::MissingValue)?;
    let t85 = h85.temperature.ok_or(AnalysisError::MissingValue)?;

    let td7 = h7.dew_point.ok_or(AnalysisError::MissingValue)?;
    let td85 = h85.dew_point.ok_or(AnalysisError::MissingValue)?;

    Ok(td7 + (t85 - t5) + (td85 - t7))
}

/// Precipitable water (mm)
#[inline]
pub fn precipitable_water(snd: &Sounding) -> Result<Mm> {
    let p_profile = snd.pressure_profile();
    let dp_profile = snd.dew_point_profile();

    let (integrated_mw, _, _) = izip!(p_profile, dp_profile)
        .filter_map(|pair| {
            let (p, dp) = pair;
            if p.is_some() && dp.is_some() {
                Some((p.unpack(), dp.unpack()))
            } else {
                None
            }
        })
        .filter_map(|(p, dp)| metfor::mixing_ratio(dp, p).and_then(|mw| Some((p, mw))))
        .fold(
            (0.0, HectoPascal(0.0), 0.0),
            |(mut acc_mw, prev_p, prev_mw), (p, mw)| {
                let dp = prev_p - p;
                if dp > HectoPascal(0.0) {
                    acc_mw += (mw + prev_mw) * dp.unpack();
                }

                (acc_mw, p, mw)
            },
        );

    Ok(Mm(integrated_mw / 9.81 / 997.0 * 100_000.0 / 2.0))
}

/// The Haines index for fire weather.
#[inline]
pub fn haines(snd: &Sounding) -> Result<u8> {
    snd.station_info()
        .elevation()
        .into_option()
        .ok_or(AnalysisError::MissingValue)
        .and_then(|elev| {
            if elev <= Meters(304.8) {
                haines_low(snd)
            } else if elev <= Meters(914.4) {
                haines_mid(snd)
            } else {
                haines_high(snd)
            }
        })
}

/// The low level version of the Haines index for fire weather.
#[inline]
pub fn haines_low(snd: &Sounding) -> Result<u8> {
    let level1 = linear_interpolate_sounding(snd, HectoPascal(950.0))
        .map_err(|_| AnalysisError::MissingValue)?;
    let level2 = linear_interpolate_sounding(snd, HectoPascal(850.0))
        .map_err(|_| AnalysisError::MissingValue)?;

    let Celsius(t_low) = level1.temperature.ok_or(AnalysisError::MissingValue)?;
    let Celsius(t_hi) = level2.temperature.ok_or(AnalysisError::MissingValue)?;
    let Celsius(dp_hi) = level2.dew_point.ok_or(AnalysisError::MissingValue)?;

    let stability_term = (t_low - t_hi).round();
    let stability_term = if stability_term >= 8.0 {
        3
    } else if stability_term > 3.0 {
        2
    } else {
        1
    };

    let moisture_term = (t_hi - dp_hi).round();
    let moisture_term = if moisture_term >= 10.0 {
        3
    } else if moisture_term > 5.0 {
        2
    } else {
        1
    };

    Ok(stability_term + moisture_term)
}

/// The mid level version of the Haines index for fire weather.
#[inline]
pub fn haines_mid(snd: &Sounding) -> Result<u8> {
    let level1 = linear_interpolate_sounding(snd, HectoPascal(850.0))
        .map_err(|_| AnalysisError::MissingValue)?;
    let level2 = linear_interpolate_sounding(snd, HectoPascal(700.0))
        .map_err(|_| AnalysisError::MissingValue)?;

    let Celsius(t_low) = level1.temperature.ok_or(AnalysisError::MissingValue)?;
    let Celsius(t_hi) = level2.temperature.ok_or(AnalysisError::MissingValue)?;
    let Celsius(dp_low) = level1.dew_point.ok_or(AnalysisError::MissingValue)?;

    let stability_term = (t_low - t_hi).round();
    let stability_term = if stability_term >= 11.0 {
        3
    } else if stability_term > 5.0 {
        2
    } else {
        1
    };

    let moisture_term = (t_low - dp_low).round();
    let moisture_term = if moisture_term >= 13.0 {
        3
    } else if moisture_term > 5.0 {
        2
    } else {
        1
    };

    Ok(stability_term + moisture_term)
}

/// The high level version of the Haines index for fire weather.
#[inline]
pub fn haines_high(snd: &Sounding) -> Result<u8> {
    let level1 = linear_interpolate_sounding(snd, HectoPascal(700.0))
        .map_err(|_| AnalysisError::MissingValue)?;
    let level2 = linear_interpolate_sounding(snd, HectoPascal(500.0))
        .map_err(|_| AnalysisError::MissingValue)?;

    let Celsius(t_low) = level1.temperature.ok_or(AnalysisError::MissingValue)?;
    let Celsius(t_hi) = level2.temperature.ok_or(AnalysisError::MissingValue)?;
    let Celsius(dp_low) = level1.dew_point.ok_or(AnalysisError::MissingValue)?;

    let stability_term = (t_low - t_hi).round();
    let stability_term = if stability_term >= 22.0 {
        3
    } else if stability_term > 17.0 {
        2
    } else {
        1
    };

    let moisture_term = (t_low - dp_low).round();
    let moisture_term = if moisture_term >= 21.0 {
        3
    } else if moisture_term > 14.0 {
        2
    } else {
        1
    };

    Ok(stability_term + moisture_term)
}

/// The Hot-Dry-Windy index
#[inline]
pub fn hot_dry_windy(snd: &Sounding) -> Result<f64> {
    let elevation = if let Some(sfc_h) = snd.station_info().elevation().into_option() {
        sfc_h
    } else if let Some(lowest_h) = snd
        .height_profile()
        .iter()
        .filter_map(|optd| optd.into_option())
        .nth(0)
    {
        lowest_h
    } else {
        return Err(AnalysisError::NotEnoughData);
    };

    let h_profile = snd.height_profile();
    let t_profile = snd.temperature_profile();
    let dp_profile = snd.dew_point_profile();
    let ws_profile = snd.wind_profile();

    let (vpd, ws) = izip!(h_profile, t_profile, dp_profile, ws_profile)
        // Remove rows with missing data
        .filter_map(|(h, t, dp, ws)| {
            if h.is_some() && t.is_some() && dp.is_some() && ws.is_some() {
                let WindSpdDir { speed: ws, .. } = ws.unpack();
                Some((h.unpack(), t.unpack(), dp.unpack(), ws)) // unpack optional::Optioned
            } else {
                None
            }
        })
        // Only look up to 500 m above AGL
        .take_while(|(h, _, _, _)| *h <= elevation + Meters(500.0))
        // Convert t and dp to VPD
        .filter_map(|(_, t, dp, ws)| {
            vapor_pressure_liquid_water(t).and_then(|sat_vap| {
                vapor_pressure_liquid_water(dp).and_then(|vap| Some((sat_vap - vap, ws)))
            })
        })
        // Convert knots to m/s and unpack all values.
        .map(|(vpd, ws)| (vpd.unpack(), MetersPSec::from(ws).unpack()))
        // Choose the max.
        .fold((0.0, 0.0), |(vpd_max, ws_max), (vpd, ws)| {
            (vpd.max(vpd_max), ws.max(ws_max))
        });

    Ok(vpd * ws)
}
