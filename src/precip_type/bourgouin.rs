//! This module implements the Bourgouin Precip Type algorithm.

use crate::{
    error::Result,
    layers::{
        cold_surface_temperature_layer, melting_freezing_energy_area,
        warm_surface_temperature_layer, warm_temperature_layer_aloft, Layer,
    },
    precip_type::{is_drizzler, PrecipType},
    sounding::Sounding,
};
use metfor::JpKg;

/// Analyze a sounding using the Bourgouin technique for precipitation type.
///
/// If the RH is below 80% in the dendritic layer, then it is assumed no cloud ice forms and an
/// appropriate drizzle type is selected. This does not apply for type D soundings (from the
/// original paper) since warm rain processes can not be ruled out.
///
/// Since the sounding generally doesn't come with information convective or stratiform and the
/// amount of precipitation, this function can only determine the phase (ice or liquid) and whether
/// it is likely to freeze once it reaches the surface (eg freezing rain). So all returned weather
/// codes will assume stratiform conditions and a light intensity.
pub fn bourgouin_precip_type(snd: &Sounding) -> Result<PrecipType> {
    let b_type = analyze_bourgouin_type(snd)?;
    let mut p_type = match b_type {
        BourgouinType::A {
            ref cold_surface,
            ref warm_layer_aloft,
        } => bourgouin_type_a_precip_type(snd, cold_surface, warm_layer_aloft)?,
        BourgouinType::B {
            ref warm_sfc_layer,
            ref warm_layer_aloft,
        } => bourgouin_type_b_precip_type(snd, warm_sfc_layer, warm_layer_aloft)?,
        BourgouinType::C { ref warm_sfc_layer } => {
            bourgouin_type_c_precip_type(snd, warm_sfc_layer)?
        }
        BourgouinType::D => bourgouin_type_d_precip_type(snd)?,
    };

    if is_drizzler(snd) {
        match b_type {
            BourgouinType::A { .. } | BourgouinType::D => p_type = PrecipType::LightFreezingDrizzle,
            BourgouinType::B { .. } => match p_type {
                PrecipType::LightIcePellets | PrecipType::LightSnow => {
                    p_type = PrecipType::LightFreezingDrizzle
                }
                PrecipType::LightRain | PrecipType::LightRainAndSnow => {
                    p_type = PrecipType::LightDrizzle
                }
                _ => {}
            },
            // Can't tell the difference between drizzle and warm rain processes.
            BourgouinType::C { .. } => {}
        }
    }

    Ok(p_type)
}

#[derive(Debug)]
#[allow(clippy::large_enum_variant)]
enum BourgouinType {
    A {
        cold_surface: Layer,
        warm_layer_aloft: Layer,
    },
    B {
        warm_sfc_layer: Layer,
        warm_layer_aloft: Layer,
    },
    C {
        warm_sfc_layer: Layer,
    },
    D,
}

fn analyze_bourgouin_type(snd: &Sounding) -> Result<BourgouinType> {
    let warm_layers = warm_temperature_layer_aloft(snd)?;
    let warm_layer_aloft: Option<Layer> = warm_layers.get(0).cloned();

    let b_type = match warm_layer_aloft {
        Some(warm_layer_aloft) => {
            if let Some(cold_surface) = cold_surface_temperature_layer(snd, &warm_layers)? {
                BourgouinType::A {
                    warm_layer_aloft,
                    cold_surface,
                }
            } else if let Some(warm_sfc_layer) = warm_surface_temperature_layer(snd)? {
                BourgouinType::B {
                    warm_layer_aloft,
                    warm_sfc_layer,
                }
            } else {
                unreachable!()
            }
        }
        None => {
            if let Some(warm_sfc_layer) = warm_surface_temperature_layer(snd)? {
                BourgouinType::C { warm_sfc_layer }
            } else {
                BourgouinType::D
            }
        }
    };

    Ok(b_type)
}

fn bourgouin_type_a_precip_type(
    snd: &Sounding,
    cold_surface: &Layer,
    warm_aloft: &Layer,
) -> Result<PrecipType> {
    let positive_area = melting_freezing_energy_area(snd, warm_aloft)?;
    debug_assert!(positive_area >= JpKg(0.0));

    if positive_area < JpKg(2.0) {
        return Ok(PrecipType::LightSnow);
    }

    let negative_area = melting_freezing_energy_area(snd, cold_surface)?;
    debug_assert!(negative_area <= JpKg(0.0));

    // 46 and 66 should be the values with anything in the middle being a mixed IP FZRA,
    // but I don't yet have a code for that.
    let threshold = JpKg(56.0) + positive_area * 0.66;

    let p_type = if negative_area > threshold {
        PrecipType::LightIcePellets
    } else {
        PrecipType::LightFreezingRain
    };

    Ok(p_type)
}

fn bourgouin_type_b_precip_type(
    snd: &Sounding,
    warm_surface: &Layer,
    warm_aloft: &Layer,
) -> Result<PrecipType> {
    let cold_layer_aloft = &Layer {
        bottom: warm_surface.top,
        top: warm_aloft.bottom,
    };

    let top_p_type = bourgouin_type_a_precip_type(snd, cold_layer_aloft, warm_aloft)?;

    let p_type = match top_p_type {
        PrecipType::LightSnow | PrecipType::LightIcePellets => {
            bourgouin_type_c_precip_type(snd, warm_surface)?
        }
        _ => PrecipType::LightRain,
    };

    Ok(p_type)
}

fn bourgouin_type_c_precip_type(snd: &Sounding, warm_surface: &Layer) -> Result<PrecipType> {
    let melting_energy = melting_freezing_energy_area(snd, warm_surface)?;
    debug_assert!(melting_energy >= JpKg(0.0));

    let p_type = if melting_energy < JpKg(5.6) {
        PrecipType::LightSnow
    } else if melting_energy > JpKg(13.2) {
        PrecipType::LightRain
    } else {
        PrecipType::LightRainAndSnow
    };

    Ok(p_type)
}

fn bourgouin_type_d_precip_type(_snd: &Sounding) -> Result<PrecipType> {
    Ok(PrecipType::LightSnow)
}
