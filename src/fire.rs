//! Module for analysis related to fire weather and wildfire plumes.

use crate::{
    error::{AnalysisError, Result},
    sounding::Sounding,
};
use itertools::izip;
use metfor::{HectoPascal, Meters, MetersPSec, Quantity};

/// The Hot-Dry-Windy index
///
/// Srock AF, Charney JJ, Potter BE, Goodrick SL. The Hot-Dry-Windy Index: A New Fire Weather Index.
///     Atmosphere. 2018; 9(7):279. https://doi.org/10.3390/atmos9070279
///
/// The implementation is described in the above cited paper. The paper below is an additional
/// reference.
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
