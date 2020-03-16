//! Indexes that are specific to a sounding, but not a particular parcel analysis of that sounding.

use crate::{
    error::{AnalysisError, Result},
    interpolation::linear_interpolate_sounding,
    sounding::Sounding,
};
use itertools::{izip, Itertools};
use metfor::{
    mixing_ratio, vapor_pressure_liquid_water, Celsius, HectoPascal, Meters, MetersPSec, Mm,
    Quantity,
};

/// Precipitable water (mm)
#[inline]
pub fn precipitable_water(snd: &Sounding) -> Result<Mm> {
    let p_profile = snd.pressure_profile();
    let dp_profile = snd.dew_point_profile();

    let integrated_mw = izip!(p_profile, dp_profile)
        // Remove levels with missing data
        .filter(|(p, dp)| p.is_some() && dp.is_some())
        // Unpack from the Optioned type
        .map(|(p, dp)| (p.unpack(), dp.unpack()))
        // Converte dew point to mixing ratio, removing failed levels.
        .filter_map(|(p, dp)| mixing_ratio(dp, p).map(|mw| (p, mw)))
        // View them as pairs for integration with the trapezoid method
        .tuple_windows::<(_, _)>()
        // Do the sum for integrating
        .fold(0.0, |mut acc_mw, ((p0, mw0), (p1, mw1))| {
            let dp = p0 - p1;
            acc_mw += (mw0 + mw1) * dp.unpack();

            acc_mw
        });

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
        .next()
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
        .filter(|(h, t, dp, ws)| h.is_some() && t.is_some() && dp.is_some() && ws.is_some())
        // Unpack from the Optioned type
        .map(|(h, t, dp, ws)| (h.unpack(), t.unpack(), dp.unpack(), ws.unpack().speed))
        // Only look up to 500 m above AGL
        .take_while(|(h, _, _, _)| *h <= elevation + Meters(500.0))
        // Convert t and dp to VPD, and remove any levels that error on calculating vapor pressure
        .filter_map(|(_, t, dp, ws)| {
            vapor_pressure_liquid_water(t)
                .and_then(|sat_vap| vapor_pressure_liquid_water(dp).map(|vap| (sat_vap - vap, ws)))
        })
        // Convert knots to m/s and unpack all values from their Quantity types
        .map(|(vpd, ws)| (vpd.unpack(), MetersPSec::from(ws).unpack()))
        // Choose the max.
        .fold((0.0, 0.0), |(vpd_max, ws_max), (vpd, ws)| {
            (vpd.max(vpd_max), ws.max(ws_max))
        });

    Ok(vpd * ws)
}
