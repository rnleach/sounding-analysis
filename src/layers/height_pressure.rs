use super::Layer;
use crate::{
    error::{AnalysisError::MissingValue, Result},
    interpolation::linear_interpolate_sounding,
    levels::height_level,
    sounding::{DataRow, Sounding},
};
use metfor::{HectoPascal, Meters};

/// Get a layer that has a certain thickness, like 3km or 6km.
#[inline]
pub fn layer_agl(snd: &Sounding, meters_agl: Meters) -> Result<Layer> {
    let tgt_elev = snd
        .station_info()
        .elevation()
        .map(|e| e + meters_agl)
        .ok_or(MissingValue)?;

    // First row is surface data if present.
    let mut bottom = snd.data_row(0).unwrap_or_else(DataRow::default);

    if bottom.pressure.is_none()
        || bottom.temperature.is_none()
        || bottom.dew_point.is_none()
        || bottom.wind.is_none()
    {
        bottom = snd.data_row(1).unwrap_or_else(DataRow::default);
    }

    let top = height_level(tgt_elev, snd)?;
    Ok(Layer { bottom, top })
}

/// Get a layer defined by two pressure levels. `bottom_p` > `top_p`
#[inline]
pub fn pressure_layer(snd: &Sounding, bottom_p: HectoPascal, top_p: HectoPascal) -> Result<Layer> {
    debug_assert!(bottom_p > top_p);

    let bottom = linear_interpolate_sounding(snd, bottom_p)?;
    let top = linear_interpolate_sounding(snd, top_p)?;

    Ok(Layer { bottom, top })
}
