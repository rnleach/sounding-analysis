//! This module contains types and functions related to assessing the precipitation type from a
//! sounding. The types are based on the WMO Manual On Codes Vol I.1, Part A Alphanumeric Codes.
//!
//! The codes are drawn from table 4680 for automated weather stations. Note that not all observable
//! weather types are diagnosable from a sounding (eg blowing sand, dust, etc), and only
//! precipitation types are considered for this crate, so no fog or mist will be included. Most
//! algorithms for diagnosing precipitation type do not classify all possible precipitation types.
//!
//! Thunderstorms are also not diagnosed. Though it is possible that someday I may try to implement
//! algorithms to diagnose fog, mist, and thunder from soundings.

use crate::sounding::Sounding;
use itertools::izip;
use metfor::{Celsius, Mm, StatuteMiles};
use std::{convert::From, fmt::Display};
use strum_macros::EnumIter;

mod bourgouin;
pub use bourgouin::bourgouin_precip_type;
mod nssl;
pub use nssl::nssl_precip_type;

/// Precipitation type enum. Values are meant to correspond to the code values from table 4680 in
/// the WMO Manual On Codes Vol I.1 Part A, Alphanumeric Codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, EnumIter, Hash, PartialOrd, Ord)]
#[repr(u8)]
#[allow(missing_docs)]
pub enum PrecipType {
    // No weather
    None = 0,

    // Drizzle category
    Drizzle = 50,
    LightDrizzle = 51,
    ModerateDrizzle = 52,
    HeavyDrizzle = 53,
    LightFreezingDrizzle = 54,
    ModerateFreezingDrizzle = 55,
    HeavyFreezingDrizzle = 56,
    LightDrizzleAndRain = 57,
    ModerateDrizzleAndRain = 58,
    // DrizzleReserved = 59,

    // Rain category
    Rain = 60,
    LightRain = 61,
    ModerateRain = 62,
    HeavyRain = 63,
    LightFreezingRain = 64,
    ModerateFreezingRain = 65,
    HeavyFreezingRain = 66,
    LightRainAndSnow = 67,
    ModerateRainAndSnow = 68,
    // RainReserved = 69,

    // Snow/Ice category
    Snow = 70,
    LightSnow = 71,
    ModerateSnow = 72,
    HeavySnow = 73,
    LightIcePellets = 74,
    ModerateIcePellets = 75,
    HeavyIcePellets = 76,
    SnowGrains = 77,
    IceCrystals = 78,
    // SnowIceReserved = 79

    // Showers
    Showers = 80,
    LightRainShowers = 81,
    ModerateRainShowers = 82,
    HeavyRainShowers = 83,
    ViolentRainShowers = 84,
    LightSnowShowers = 85,
    ModerateSnowShowers = 86,
    HeavySnowShowers = 87,
    // ShowersReserved = 88,
    Hail = 89,

    // 90 series is for thunderstorms

    // Catch all
    Unknown = 100,
}

impl From<u8> for PrecipType {
    fn from(val: u8) -> Self {
        use PrecipType::*;

        match val {
            0 => None,

            // Drizzle category
            50 => Drizzle,
            51 => LightDrizzle,
            52 => ModerateDrizzle,
            53 => HeavyDrizzle,
            54 => LightFreezingDrizzle,
            55 => ModerateFreezingDrizzle,
            56 => HeavyFreezingDrizzle,
            57 => LightDrizzleAndRain,
            58 => ModerateDrizzleAndRain,

            // Rain category
            60 => Rain,
            61 => LightRain,
            62 => ModerateRain,
            63 => HeavyRain,
            64 => LightFreezingRain,
            65 => ModerateFreezingRain,
            66 => HeavyFreezingRain,
            67 => LightRainAndSnow,
            68 => ModerateRainAndSnow,

            // Snow/Ice category
            70 => Snow,
            71 => LightSnow,
            72 => ModerateSnow,
            73 => HeavySnow,
            74 => LightIcePellets,
            75 => ModerateIcePellets,
            76 => HeavyIcePellets,
            77 => SnowGrains,
            78 => IceCrystals,

            // Showers
            80 => Showers,
            81 => LightRainShowers,
            82 => ModerateRainShowers,
            83 => HeavyRainShowers,
            84 => ViolentRainShowers,
            85 => LightSnowShowers,
            86 => ModerateSnowShowers,
            87 => HeavySnowShowers,

            89 => Hail,

            _ => Unknown,
        }
    }
}

impl Display for PrecipType {
    fn fmt(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(formatter, "{} => {:?}", *self as u8, self)
    }
}

fn is_drizzler(snd: &Sounding) -> bool {
    let ts = snd.temperature_profile();
    let tds = snd.dew_point_profile();

    !izip!(ts, tds)
        // Remove levels with missing data
        .filter(|(t, td)| t.is_some() && td.is_some())
        // Unpack `Optioned` values
        .map(|(t, td)| (t.unpack(), td.unpack()))
        // Filter out non-dendritic zone temperatures
        .filter(|&(t, _)| t <= Celsius(-12.0) && t >= Celsius(-18.0))
        // Map to RH values
        .filter_map(|(t, td)| metfor::rh(t, td))
        // Find if any are over 90% rh
        .any(|rh_val| rh_val >= 0.80)
}

/// Given a `PrecipType` and potentially hourly precipitation, hourly convective precipitation, and
/// visibility, make a best effor to correct the weather type code, or PrecipType.
///
/// The adjustments are based on the definition of light, moderate, and heavy under the rain entry
/// in the American Meteorological Society Glossary of Meteorology, and if the original PrecipType
/// is a snow type and visibility information is available it is based off the visibility.
///
/// This function is not exhuastive, but it is a best effort attempt to handle the most common
/// precipitation types: rain, snow, freezing rain, ice pellets, and mixtures of these. It only
/// works with the light intensity level for these types since that is the default return type from
/// the algorithms provided here which have no way of knowing intensity. If a weather code provided
/// by some other provider has some more information, it should be used rather than adjusted by
/// this routine.
///
pub fn check_precip_type_intensity<S1, S2, V>(
    precip_type: PrecipType,
    hourly_precipitation: Option<S1>,
    convective_precipitation: Option<S2>,
    visibility: Option<V>,
) -> PrecipType
where
    S1: Into<Mm>,
    S2: Into<Mm>,
    V: Into<StatuteMiles>,
{
    use PrecipType::*;

    // Force into units we like.
    let hourly_precipitation = hourly_precipitation.map(Into::<Mm>::into);
    let convective_precipitation = convective_precipitation.map(Into::<Mm>::into);
    let visibility = visibility.map(Into::<StatuteMiles>::into);

    // Check whether we have enough information to make an adjustment.
    match check_preconditions_for_wx_type_adjustment(precip_type, hourly_precipitation, visibility)
    {
        PreconditionsCheckResult::Pass => {}
        PreconditionsCheckResult::Fail => return precip_type,
        PreconditionsCheckResult::NoQPF => return PrecipType::None,
    }

    let mode = get_precip_type_mode(hourly_precipitation, convective_precipitation);

    let intensity = get_precip_type_intensity(precip_type, hourly_precipitation, visibility);

    match precip_type {
        LightDrizzle => match intensity {
            Intensity::Light => LightDrizzle,
            Intensity::Moderate => ModerateDrizzle,
            Intensity::Heavy => HeavyDrizzle,
        },
        LightFreezingDrizzle => match intensity {
            Intensity::Light => LightFreezingDrizzle,
            Intensity::Moderate => ModerateFreezingDrizzle,
            Intensity::Heavy => HeavyFreezingDrizzle,
        },
        LightDrizzleAndRain => match intensity {
            Intensity::Light => LightDrizzleAndRain,
            Intensity::Moderate | Intensity::Heavy => ModerateDrizzleAndRain,
        },
        LightRain => match intensity {
            Intensity::Light => match mode {
                Mode::Convective => LightRainShowers,
                Mode::Stratiform => LightRain,
            },
            Intensity::Moderate => match mode {
                Mode::Convective => ModerateRainShowers,
                Mode::Stratiform => ModerateRain,
            },
            Intensity::Heavy => match mode {
                Mode::Convective => HeavyRainShowers,
                Mode::Stratiform => HeavyRain,
            },
        },
        LightFreezingRain => match intensity {
            Intensity::Light => LightFreezingRain,
            Intensity::Moderate => ModerateFreezingRain,
            Intensity::Heavy => HeavyFreezingRain,
        },
        LightRainAndSnow => match intensity {
            Intensity::Light => LightRainAndSnow,
            Intensity::Moderate | Intensity::Heavy => ModerateRainAndSnow,
        },
        LightSnow => match intensity {
            Intensity::Light => match mode {
                Mode::Convective => LightSnowShowers,
                Mode::Stratiform => LightSnow,
            },
            Intensity::Moderate => match mode {
                Mode::Convective => ModerateSnowShowers,
                Mode::Stratiform => ModerateSnow,
            },
            Intensity::Heavy => match mode {
                Mode::Convective => HeavySnowShowers,
                Mode::Stratiform => HeavySnow,
            },
        },
        LightIcePellets => match intensity {
            Intensity::Light => LightIcePellets,
            Intensity::Moderate => ModerateIcePellets,
            Intensity::Heavy => HeavyIcePellets,
        },
        // Anything else should be unreachable due to the preconditions checks.
        _ => unreachable!(),
    }
}

enum Mode {
    Stratiform,
    Convective,
}

enum Intensity {
    Light,
    Moderate,
    Heavy,
}

enum PreconditionsCheckResult {
    Pass,
    Fail,
    NoQPF,
}

/// Check all the preconditions for adjusting the wx type.
fn check_preconditions_for_wx_type_adjustment(
    precip_type: PrecipType,
    hourly_qpf: Option<Mm>,
    visibility: Option<StatuteMiles>,
) -> PreconditionsCheckResult {
    use PrecipType::*;
    use PreconditionsCheckResult::*;

    // Bail out if it isn't one of the cases which are light.
    match precip_type {
        LightDrizzle | LightFreezingDrizzle | LightDrizzleAndRain | LightRain
        | LightFreezingRain | LightRainAndSnow | LightSnow | LightIcePellets => {}
        _ => return Fail,
    }

    // Bail if there is not QPF
    match hourly_qpf {
        Some(pcp) => {
            if pcp <= Mm(0.0) {
                return NoQPF;
            }
        }
        std::option::Option::None => return Fail,
    }

    // Bail if we don't have enough information to make an adjustment
    match precip_type {
        LightDrizzle | LightFreezingDrizzle | LightDrizzleAndRain | LightRain
        | LightFreezingRain | LightRainAndSnow | LightIcePellets => {
            if hourly_qpf.is_none() {
                return Fail;
            }
        }
        LightSnow => {
            if hourly_qpf.is_none() && visibility.is_none() {
                return Fail;
            }
        }
        _ => unreachable!(),
    }
    Pass
}

fn get_precip_type_mode(hourly_qpf: Option<Mm>, convective_qpf: Option<Mm>) -> Mode {
    hourly_qpf
        .and_then(|total_pcp| {
            convective_qpf.map(|conv_pcp| {
                if conv_pcp > total_pcp * 0.5 {
                    Mode::Convective
                } else {
                    Mode::Stratiform
                }
            })
        })
        .unwrap_or(Mode::Stratiform)
}

fn get_precip_type_intensity(
    precip_type: PrecipType,
    hourly_qpf: Option<Mm>,
    visibility: Option<StatuteMiles>,
) -> Intensity {
    use PrecipType::*;

    match precip_type {
        LightDrizzleAndRain | LightRain | LightFreezingRain | LightRainAndSnow
        | LightIcePellets => {
            if let Some(pcp) = hourly_qpf {
                qpf_to_intensity_for_rain(pcp)
            } else {
                Intensity::Light
            }
        }
        LightDrizzle | LightFreezingDrizzle => {
            if let Some(pcp) = hourly_qpf {
                qpf_to_intensity_for_drizzle(pcp)
            } else {
                Intensity::Light
            }
        }
        LightSnow => {
            if let Some(vsby) = visibility {
                visibility_to_intensity(vsby)
            } else if let Some(pcp) = hourly_qpf {
                qpf_to_intensity_for_rain(pcp)
            } else {
                Intensity::Light
            }
        }
        // Should be unreachable because of the function precondition checks done earlier.
        _ => unreachable!(),
    }
}

// Definition based on the AMS Glossary of Meteorology definition for rain.
fn qpf_to_intensity_for_rain(qpf: Mm) -> Intensity {
    if qpf <= Mm(2.5) {
        Intensity::Light
    } else if qpf <= Mm(7.6) {
        Intensity::Moderate
    } else {
        Intensity::Heavy
    }
}

// Definition based on the AMS Glossary of Meteorology definition for drizzle.
fn qpf_to_intensity_for_drizzle(qpf: Mm) -> Intensity {
    if qpf <= Mm(0.3) {
        Intensity::Light
    } else if qpf <= Mm(0.5) {
        Intensity::Moderate
    } else {
        Intensity::Heavy
    }
}

// Definition base on the AMS Glossary of Meteorology definition for snow.
fn visibility_to_intensity(vsby: StatuteMiles) -> Intensity {
    if vsby < StatuteMiles(5.0 / 16.0) {
        Intensity::Heavy
    } else if vsby < StatuteMiles(5.0 / 8.0) {
        Intensity::Moderate
    } else {
        Intensity::Light
    }
}

#[cfg(test)]
mod test {
    use super::PrecipType::*;
    use super::*;
    use strum::IntoEnumIterator;

    #[test]
    #[rustfmt::skip]
    fn test_adjust_check_precip_type_intensity() {

        assert_eq!(check_precip_type_intensity(LightDrizzle, Some(Mm(0.01)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), LightDrizzle);
        assert_eq!(check_precip_type_intensity(LightDrizzle, Some(Mm(0.4)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), ModerateDrizzle);
        assert_eq!(check_precip_type_intensity(LightDrizzle, Some(Mm(0.6)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), HeavyDrizzle);

        assert_eq!(check_precip_type_intensity(LightFreezingDrizzle, Some(Mm(0.1)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), LightFreezingDrizzle);
        assert_eq!(check_precip_type_intensity(LightFreezingDrizzle, Some(Mm(0.4)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), ModerateFreezingDrizzle);
        assert_eq!(check_precip_type_intensity(LightFreezingDrizzle, Some(Mm(0.6)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), HeavyFreezingDrizzle);

        assert_eq!(check_precip_type_intensity(LightRain, Some(Mm(1.0)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), LightRain);
        assert_eq!(check_precip_type_intensity(LightRain, Some(Mm(5.0)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), ModerateRain);
        assert_eq!(check_precip_type_intensity(LightRain, Some(Mm(8.0)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), HeavyRain);

        assert_eq!(check_precip_type_intensity(LightFreezingRain, Some(Mm(1.0)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), LightFreezingRain);
        assert_eq!(check_precip_type_intensity(LightFreezingRain, Some(Mm(5.0)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), ModerateFreezingRain);
        assert_eq!(check_precip_type_intensity(LightFreezingRain, Some(Mm(8.0)), Some(Mm(0.0)), Some(StatuteMiles(10.0))), HeavyFreezingRain);

        assert_eq!(check_precip_type_intensity(LightSnow, Some(Mm(1.0)), Some(Mm(0.0)), std::option::Option::<StatuteMiles>::None), LightSnow);
        assert_eq!(check_precip_type_intensity(LightSnow, Some(Mm(5.0)), Some(Mm(0.0)), std::option::Option::<StatuteMiles>::None), ModerateSnow);
        assert_eq!(check_precip_type_intensity(LightSnow, Some(Mm(8.0)), Some(Mm(0.0)), std::option::Option::<StatuteMiles>::None), HeavySnow);

        for variant in PrecipType::iter() {
            // Skip special cases handled above
            match variant {
                LightIcePellets | LightRainAndSnow | LightDrizzleAndRain | LightDrizzle | 
                    LightFreezingDrizzle | LightRain | LightSnow | LightFreezingRain => {
                    continue
                }
                _ => {}
            }

            for &qpf in &[Mm(1.0), Mm(5.0), Mm(8.0)] {
                assert_eq!(check_precip_type_intensity(variant, Some(qpf), Some(Mm(0.0)), Some(StatuteMiles(10.0))), variant);
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
