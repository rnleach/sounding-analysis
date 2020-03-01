//! This module contains types and functions related to assessing the precipitation type from a
//! sounding. The types are based on the WMO Manual On Codes Vol I.1, Part A Alphanumeric Codes.
//!
//! The codes are drawn from table 4677 for manned weather stations. Note that not all observable
//! weather types are diagnosable from a sounding (eg blowing sand, dust, etc), and only
//! precipitation types are considered for this crate, so no fog or mist will be included. Most
//! algorithms for diagnosing precipitation type do not classify all possible precipitation types.

use crate::{
    error::Result,
    layers::{
        cold_surface_temperature_layer, melting_freezing_energy_area,
        warm_surface_temperature_layer, warm_temperature_layer_aloft, Layer,
    },
    sounding::Sounding,
};
use metfor::{JpKg, Mm, Quantity};
use std::convert::From;
use strum_macros::EnumIter;

/// Precipitation type enum. Values are meant to correspond to the code values from table 4677 in
/// the WMO Manual On Codes Vol I.1 Part A, Alphanumeric Codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, EnumIter)]
#[repr(u8)]
#[allow(missing_docs)]
pub enum PrecipType {
    // No weather
    None = 0,

    // Drizzle category
    IntermittentLightDrizzle = 50,
    ContinuousLightDrizzle = 51,
    IntermittentModerateDrizzle = 52,
    ContinuousModerateDrizzle = 53,
    IntermittentHeavyDrizzle = 54,
    ContinuousHeavyDrizzle = 55,
    LightFreezingDrizzle = 56,
    ModerateFreezingDrizzle = 57,

    // Rain category
    LightRainShowers = 60,
    LightRain = 61,
    ModerateRainShowers = 62,
    ModerateRain = 63,
    HeavyRainShowers = 64,
    HeavyRain = 65,
    LightFreezingRain = 66,
    ModerateFreezingRain = 67,
    LightRainAndSnow = 68,
    ModerateRainAndSnow = 69,

    // Snow/Ice category
    LightSnowShowers = 70,
    LightSnow = 71,
    ModerateSnowShowers = 72,
    ModerateSnow = 73,
    HeavySnowShowers = 74,
    HeavySnow = 75,
    DiamondDust = 76,
    SnowGrains = 77,
    IsolatedStarLikeSnowCrystals = 78,
    IcePellets = 79,

    // Catch all
    Unknown = 100,
}

impl From<u8> for PrecipType {
    fn from(val: u8) -> Self {
        use PrecipType::*;

        match val {
            0 => None,

            // Drizzle category
            50 => IntermittentLightDrizzle,
            51 => ContinuousLightDrizzle,
            52 => IntermittentModerateDrizzle,
            53 => ContinuousModerateDrizzle,
            54 => IntermittentHeavyDrizzle,
            55 => ContinuousHeavyDrizzle,
            56 => LightFreezingDrizzle,
            57 => ModerateFreezingDrizzle,

            // Rain category
            60 => LightRainShowers,
            61 => LightRain,
            62 => ModerateRainShowers,
            63 => ModerateRain,
            64 => HeavyRainShowers,
            65 => HeavyRain,
            66 => LightFreezingRain,
            67 => ModerateFreezingRain,
            68 => LightRainAndSnow,
            69 => ModerateRainAndSnow,

            // Snow/Ice category
            70 => LightSnowShowers,
            71 => LightSnow,
            72 => ModerateSnowShowers,
            73 => ModerateSnow,
            74 => HeavySnowShowers,
            75 => HeavySnow,
            76 => DiamondDust,
            77 => SnowGrains,
            78 => IsolatedStarLikeSnowCrystals,
            79 => IcePellets,

            _ => Unknown,
        }
    }
}

/// Analyze a sounding using the Bourgouin technique for precipitation type.
///
/// Since the sounding generally doesn't come with information convective or stratiform and the
/// amount of precipitation, this function can only determine the phase (ice or liquid) and whether
/// it is likely to freeze once it reaches the surface (eg freezing rain). So all returned weather
/// codes will assume convective or intermittent conditions (eg showers) and a light intensity.
/// This is consistent with the observed behavior of Bufkit.
pub fn bourgouin_precip_type(snd: &Sounding) -> Result<PrecipType> {
    let p_type = match analyze_bourgouin_type(snd)? {
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

    Ok(p_type)
}

#[derive(Debug)]
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
    let warm_layer_aloft: Option<Layer> = warm_layers.get(0).map(|lyr| lyr.clone());

    let b_type = match warm_layer_aloft {
        Some(warm_layer_aloft) => {
            if let Some(cold_surface) = cold_surface_temperature_layer(snd, &warm_layers)? {
                BourgouinType::A {
                    warm_layer_aloft,
                    cold_surface,
                }
            } else {
                if let Some(warm_sfc_layer) = warm_surface_temperature_layer(snd)? {
                    BourgouinType::B {
                        warm_layer_aloft,
                        warm_sfc_layer,
                    }
                } else {
                    unreachable!()
                }
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

    // FIXME: 46 and 66 should be the values with anything in the middle being a mixed IP FZRA,
    // but I don't yet have a code for that.
    let threshold = JpKg(56.0) + JpKg(0.66 * positive_area.unpack());

    let p_type = if negative_area > threshold {
        PrecipType::IcePellets
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
        PrecipType::LightSnow | PrecipType::IcePellets => {
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

/// Given a `PrecipType` and an hourly precipitation amount, adjust the intensity.
///
/// The adjustments are based on the definition of light, moderate, and heavy under the rain entry
/// in the American Meteorological Society Glossary of Meteorology.
///
/// This function is not exhuastive, but it is a best effort attempt to handle the most common
/// precipitation types: rain, snow, and freezing rain. It only works with the light intensity
/// level for these types since that is the default return type from the algorithms provided here
/// which have no way of knowing intensity. If a weather code provided by some other provider has
/// some more information, it should be used rather than adjusted by this routine.
///
/// Intensity for snow is actually determined by visibility, but that information is often not
/// available, so it is treated the same as rain which will probably cause a low bias in the
/// intensity of snow.
pub fn adjust_precip_type_intensity<L: Into<Mm>>(
    precip_type: PrecipType,
    hourly_precipitation: L,
) -> PrecipType {
    use PrecipType::*;

    let precipitation: Mm = hourly_precipitation.into();

    if precipitation <= Mm(0.0) {
        return PrecipType::None;
    }

    // Bail out if it isn't one of the cases which are light.
    match precip_type {
        LightRainShowers | LightRain | LightFreezingRain | LightSnowShowers | LightSnow => {}
        _ => return precip_type,
    }

    enum Intensity {
        Light,
        Moderate,
        Heavy,
    }

    let intensity = if precipitation <= Mm(2.5) {
        Intensity::Light
    } else if precipitation <= Mm(7.6) {
        Intensity::Moderate
    } else {
        Intensity::Heavy
    };

    match precip_type {
        LightRainShowers => match intensity {
            Intensity::Light => LightRainShowers,
            Intensity::Moderate => ModerateRainShowers,
            Intensity::Heavy => HeavyRainShowers,
        },
        LightRain => match intensity {
            Intensity::Light => LightRain,
            Intensity::Moderate => ModerateRain,
            Intensity::Heavy => HeavyRain,
        },
        LightFreezingRain => match intensity {
            Intensity::Light => LightFreezingRain,
            Intensity::Moderate => ModerateFreezingRain,
            Intensity::Heavy => ModerateFreezingRain,
        },
        LightSnowShowers => match intensity {
            Intensity::Light => LightSnowShowers,
            Intensity::Moderate => ModerateSnowShowers,
            Intensity::Heavy => HeavySnowShowers,
        },
        LightSnow => match intensity {
            Intensity::Light => LightSnow,
            Intensity::Moderate => ModerateSnow,
            Intensity::Heavy => HeavySnow,
        },
        _ => unreachable!(),
    }
}

#[cfg(test)]
mod test {
    use super::PrecipType::*;
    use super::*;
    use strum::IntoEnumIterator;

    #[test]
    #[rustfmt::skip]
    fn test_adjust_precip_type_intensity() {
        assert_eq!(adjust_precip_type_intensity(LightRainShowers, Mm(1.0)), LightRainShowers);
        assert_eq!(adjust_precip_type_intensity(LightRainShowers, Mm(5.0)), ModerateRainShowers);
        assert_eq!(adjust_precip_type_intensity(LightRainShowers, Mm(8.0)), HeavyRainShowers);

        assert_eq!(adjust_precip_type_intensity(LightRain, Mm(1.0)), LightRain);
        assert_eq!(adjust_precip_type_intensity(LightRain, Mm(5.0)), ModerateRain);
        assert_eq!(adjust_precip_type_intensity(LightRain, Mm(8.0)), HeavyRain);

        assert_eq!(adjust_precip_type_intensity(LightFreezingRain, Mm(1.0)), LightFreezingRain);
        assert_eq!(adjust_precip_type_intensity(LightFreezingRain, Mm(5.0)), ModerateFreezingRain);
        assert_eq!(adjust_precip_type_intensity(LightFreezingRain, Mm(8.0)), ModerateFreezingRain);

        assert_eq!(adjust_precip_type_intensity(LightSnowShowers, Mm(1.0)), LightSnowShowers);
        assert_eq!(adjust_precip_type_intensity(LightSnowShowers, Mm(5.0)), ModerateSnowShowers);
        assert_eq!(adjust_precip_type_intensity(LightSnowShowers, Mm(8.0)), HeavySnowShowers);

        assert_eq!(adjust_precip_type_intensity(LightSnow, Mm(1.0)), LightSnow);
        assert_eq!(adjust_precip_type_intensity(LightSnow, Mm(5.0)), ModerateSnow);
        assert_eq!(adjust_precip_type_intensity(LightSnow, Mm(8.0)), HeavySnow);

        for variant in PrecipType::iter() {
            // Skip special cases handled above
            match variant {
                LightRainShowers | LightRain | LightSnowShowers | LightSnow | LightFreezingRain => {
                    continue
                }
                _ => {}
            }

            for &qpf in &[Mm(1.0), Mm(5.0), Mm(8.0)] {
                assert_eq!(adjust_precip_type_intensity(variant, qpf), variant);
            }
        }
    }

    #[test]
    fn test_from_u8_and_back() {
        for variant in PrecipType::iter() {
            assert_eq!(variant, PrecipType::from(variant as u8));
        }
    }
}
