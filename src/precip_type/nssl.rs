//! This module implements the NSSL Precip Type algorithm.

use crate::{
    error::{AnalysisError, Result},
    precip_type::{is_drizzler, PrecipType},
    sounding::Sounding,
};
use itertools::izip;
use metfor::{Celsius, Meters};

/// Analyze a sounding using the NSSL algorithm for precipitation type.
///
/// This algorithm also analyzes wet snow vs snow, but the precip codes we have in this module
/// do not allow that distinction, so wet snow will be analyzed as snow. It also may analyze mixed
/// freezing rain and ice pellets, for which we have no code. In that case we will analyze freezing
/// rain since that is generally considered the more severe weather type.
///
/// In addition to the NSSL algorithm if the type of the sounding is type 1 - 4 (i.e. the surface
/// wet bulb temperature is below 3C) and the RH is below 80% in the dendritic layer, then an
/// appropriate drizzle type is selected. If the surface wet bulb is warm enough, there is nothing
/// to prevent warm rain processes from generating rain so it is impossible to tell the difference
/// between drizzle and rain.
///
/// Since the sounding generally doesn't come with information convective or stratiform and the
/// amount of precipitation, this function can only determine the phase (ice or liquid) and whether
/// it is likely to freeze once it reaches the surface (eg freezing rain). So all returned weather
/// codes will assume stratiform conditions and a light intensity.
pub fn nssl_precip_type(snd: &Sounding) -> Result<PrecipType> {
    let n_type = analyze_nssl_type(snd)?;
    let mut p_type = match n_type {
        NsslType::WarmSfcRain => PrecipType::LightRain,
        NsslType::One => PrecipType::LightSnow,
        NsslType::Two { h0 } => analyze_nssl_type_2(h0),
        NsslType::Three { tw_max, tw_min } => analyze_nssl_type_3(tw_max, tw_min),
        NsslType::Four { tw_max, tw_min } => analyze_nssl_type_4(tw_max, tw_min),
    };

    if is_drizzler(snd) {
        match n_type {
            NsslType::WarmSfcRain => {}
            NsslType::One => p_type = PrecipType::LightFreezingDrizzle,
            NsslType::Two { .. } => {}
            NsslType::Three { .. } => {
                if p_type == PrecipType::LightRain {
                    p_type = PrecipType::LightDrizzle;
                } else {
                    p_type = PrecipType::LightFreezingDrizzle;
                }
            }
            NsslType::Four { .. } => p_type = PrecipType::LightFreezingDrizzle,
        }
    }

    Ok(p_type)
}

#[derive(Debug)]
enum NsslType {
    WarmSfcRain,
    One,
    Two { h0: Meters },
    Three { tw_max: Celsius, tw_min: Celsius },
    Four { tw_max: Celsius, tw_min: Celsius },
}

fn analyze_nssl_type(snd: &Sounding) -> Result<NsslType> {
    let sfc_wet_bulb = snd
        .wet_bulb_profile()
        .iter()
        .cloned()
        .filter_map(|opt_val| opt_val.into_option())
        .next()
        .ok_or(AnalysisError::MissingValue)?;

    if sfc_wet_bulb >= Celsius(3.0) {
        return Ok(NsslType::WarmSfcRain);
    }

    let wbzs = crate::wet_bulb_zero_levels(snd)?;

    let n_type = if sfc_wet_bulb > Celsius(0.0) {
        // Handle type 2 and type 3

        // There should only be 1 or 3 WBZ levels, but weird stuff happens, so lets handle an above
        // freezing sfc wet bulb and 2 crossings the same as 1 crossing.
        if wbzs.len() <= 2 {
            // Type 2
            let h0: Meters = wbzs[0]
                .height
                .into_option()
                .ok_or(AnalysisError::MissingValue)?;
            let h_surface = snd
                .station_info()
                .elevation()
                .into_option()
                .ok_or(AnalysisError::MissingValue)?;
            let h0 = h0 - h_surface;
            NsslType::Two { h0 }
        } else {
            assert!(wbzs.len() >= 3);
            // Type 3
            let h2: Meters = wbzs[0]
                .height
                .into_option()
                .ok_or(AnalysisError::MissingValue)?;
            let h1: Meters = wbzs[1]
                .height
                .into_option()
                .ok_or(AnalysisError::MissingValue)?;
            let h0: Meters = wbzs[2]
                .height
                .into_option()
                .ok_or(AnalysisError::MissingValue)?;

            let (tw_max, tw_min) = izip!(snd.wet_bulb_profile(), snd.height_profile())
                .filter(|(wb, h)| wb.is_some() && h.is_some())
                .map(|(wb, h)| (wb.unpack(), h.unpack()))
                .skip_while(|(_, h)| *h <= h2)
                .take_while(|(_, h)| *h <= h0)
                .fold((None, None), |acc, (wb, h)| {
                    let (mut tw_max, mut tw_min) = acc;

                    if h < h1 {
                        if let Some(min_val) = tw_min {
                            if wb < min_val {
                                tw_min = Some(wb);
                            }
                        } else {
                            tw_min = Some(wb);
                        }
                    } else {
                        if let Some(max_val) = tw_max {
                            if wb > max_val {
                                tw_max = Some(wb);
                            }
                        } else {
                            tw_max = Some(wb);
                        }
                    }

                    (tw_max, tw_min)
                });

            let tw_max = tw_max.ok_or(AnalysisError::MissingValue)?;
            let tw_min = tw_min.ok_or(AnalysisError::MissingValue)?;

            NsslType::Three { tw_max, tw_min }
        }
    } else {
        // Handle type 1 and type 4
        if wbzs.len() <= 1 {
            // I suppose the WB could start below freezing and jump above it with the
            // wet_bulb_zero_levels function above, but that would mean it was above freezing at
            // 500 hPa and below freezing somewhere below it. Which would be REALLY weird, but
            // that's why I look at the only 1 WBZ level case.
            NsslType::One
        } else {
            assert!(wbzs.len() >= 2);
            let h1: Meters = wbzs[0]
                .height
                .into_option()
                .ok_or(AnalysisError::MissingValue)?;
            let h0: Meters = wbzs[1]
                .height
                .into_option()
                .ok_or(AnalysisError::MissingValue)?;

            let (tw_max, tw_min) = izip!(snd.wet_bulb_profile(), snd.height_profile())
                .filter(|(wb, h)| wb.is_some() && h.is_some())
                .map(|(wb, h)| (wb.unpack(), h.unpack()))
                .take_while(|(_, h)| *h <= h0)
                .fold((None, None), |acc, (wb, h)| {
                    let (mut tw_max, mut tw_min) = acc;

                    if h < h1 {
                        if let Some(min_val) = tw_min {
                            if wb < min_val {
                                tw_min = Some(wb);
                            }
                        } else {
                            tw_min = Some(wb);
                        }
                    } else {
                        if let Some(max_val) = tw_max {
                            if wb > max_val {
                                tw_max = Some(wb);
                            }
                        } else {
                            tw_max = Some(wb);
                        }
                    }

                    (tw_max, tw_min)
                });

            let tw_max = tw_max.ok_or(AnalysisError::MissingValue)?;
            let tw_min = tw_min.ok_or(AnalysisError::MissingValue)?;

            NsslType::Four { tw_max, tw_min }
        }
    };

    Ok(n_type)
}

fn analyze_nssl_type_2(h0: Meters) -> PrecipType {
    if h0 < Meters(1000.0) {
        PrecipType::LightSnow
    } else {
        PrecipType::LightRain
    }
}

fn analyze_nssl_type_3(tw_max: Celsius, tw_min: Celsius) -> PrecipType {
    if tw_max >= Celsius(0.0) && tw_max < Celsius(2.0) && tw_min < Celsius(-5.0) {
        PrecipType::LightIcePellets
    } else {
        PrecipType::LightRain
    }
}

fn analyze_nssl_type_4(tw_max: Celsius, tw_min: Celsius) -> PrecipType {
    if tw_max > Celsius(2.0) && tw_min >= Celsius(-5.0) {
        PrecipType::LightFreezingRain
    } else if tw_max < Celsius(2.0) && tw_min < Celsius(-5.0) {
        PrecipType::LightIcePellets
    } else {
        // This should be FZRA mixed with IP, but we don't have that option.
        PrecipType::LightFreezingRain
    }
}
