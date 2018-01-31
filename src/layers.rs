//! This module finds significant layers such as the dendritic snow growth zone, the hail growth
//! zone, and inversions.

use error::*;
use smallvec::SmallVec;

use sounding_base::{Profile, Sounding};
use sounding_base::Profile::*;

use Layer;

/// Find the dendtritic growth zones throughout the profile. It is unusual, but possible there is
/// more than one.
///
/// # Errors
/// If the sounding is missing a temperature or pressure profile, `error::ErrorKind::MissingProfile`
/// is returned in the result. Otherwise, if no dendritic layers are found, an empty vector is
/// returned in the `Result`
pub fn dendritic_snow_zone(snd: &Sounding) -> Result<SmallVec<[Layer; ::VEC_SIZE]>, AnalysisError> {
    const ANALYSIS_NAME: &str = "Dendritic Snow Growth Zone(s)";

    let mut to_return: SmallVec<[Layer; ::VEC_SIZE]> = SmallVec::new();

    // Dendritic snow growth zone temperature range in C
    const WARM_SIDE: f64 = -12.0;
    const COLD_SIDE: f64 = -18.0;

    let t_profile = snd.get_profile(Temperature);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() {
        return Err(AnalysisError::MissingProfile("Temperature", ANALYSIS_NAME));
    }
    if p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile("Pressure", ANALYSIS_NAME));
    }

    let mut profile = t_profile.iter().zip(p_profile);

    let mut bottom_press = ::std::f64::MAX; // Only init because compiler can't tell value not used
    let mut top_press: f64;

    // Initialize the bottom of the sounding
    let mut last_t: f64;
    let mut last_press: f64;
    loop {
        if let Some((t, press)) = profile.by_ref().next() {
            if let (Some(t), Some(press)) = (*t, *press) {
                last_t = t;
                last_press = press;
                break;
            }
        } else {
            return Err(AnalysisError::NoDataProfile(
                "Temperature and Pressure",
                ANALYSIS_NAME,
            ));
        }
    }

    // Check to see if we are already in the dendtritic zone
    if last_t <= WARM_SIDE && last_t >= COLD_SIDE {
        bottom_press = last_press;
    }

    fn push_layer(
        bottom_press: f64,
        top_press: f64,
        snd: &Sounding,
        target_vec: &mut SmallVec<[Layer; ::VEC_SIZE]>,
    ) {
        let bottom = ::interpolation::linear_interpolate(snd, bottom_press);
        let top = ::interpolation::linear_interpolate(snd, top_press);
        target_vec.push(Layer { bottom, top });
    }

    for (t, press) in profile {
        if let (Some(t), Some(press)) = (*t, *press) {
            // Do not use if-else or continue statements because a layer might be so thin that
            // you cross into and out of it between levels.

            // Crossed into zone from warm side
            if last_t > WARM_SIDE && t <= WARM_SIDE {
                bottom_press =
                    ::interpolation::linear_interp(WARM_SIDE, last_t, t, last_press, press);
            }
            // Crossed into zone from cold side
            if last_t < COLD_SIDE && t >= COLD_SIDE {
                bottom_press =
                    ::interpolation::linear_interp(COLD_SIDE, last_t, t, last_press, press);
            }
            // Crossed out of zone to warm side
            if last_t <= WARM_SIDE && t > WARM_SIDE {
                top_press = ::interpolation::linear_interp(WARM_SIDE, last_t, t, last_press, press);
                push_layer(bottom_press, top_press, snd, &mut to_return);
            }
            // Crossed out of zone to cold side
            if last_t >= COLD_SIDE && t < COLD_SIDE {
                top_press = ::interpolation::linear_interp(COLD_SIDE, last_t, t, last_press, press);
                push_layer(bottom_press, top_press, snd, &mut to_return);
            }
            last_t = t;
            last_press = press;
        }
    }

    // Check to see if we ended in a dendtritic zone
    if last_t <= WARM_SIDE && last_t >= COLD_SIDE {
        top_press = last_press;

        push_layer(bottom_press, top_press, snd, &mut to_return);
    }

    Ok(to_return)
}

/// Assuming it is below freezing at the surface, this will find the warm layers aloft using the
/// dry bulb temperature. Does not look above 500 hPa.
pub fn warm_temperature_layer_aloft(
    snd: &Sounding,
) -> Result<SmallVec<[Layer; ::VEC_SIZE]>, AnalysisError> {
    const ANALYSIS_NAME: &str = "Warm temperature layer aloft.";

    warm_layer_aloft(snd, Temperature, ANALYSIS_NAME, "Temperature")
}

/// Assuming the wet bulb temperature is below freezing at the surface, this will find the warm
/// layers aloft using the wet bulb temperature. Does not look above 500 hPa.
pub fn warm_wet_bulb_layer_aloft(
    snd: &Sounding,
) -> Result<SmallVec<[Layer; ::VEC_SIZE]>, AnalysisError> {
    const ANALYSIS_NAME: &str = "Warm wet bulb layer aloft.";

    warm_layer_aloft(snd, WetBulb, ANALYSIS_NAME, "Wet Bulb Temperature")
}

fn warm_layer_aloft(
    snd: &Sounding,
    var: Profile,
    analysis_name: &'static str,
    profile_name: &'static str,
) -> Result<SmallVec<[Layer; ::VEC_SIZE]>, AnalysisError> {
    assert!(var == Temperature || var == WetBulb);

    let mut to_return: SmallVec<[Layer; ::VEC_SIZE]> = SmallVec::new();

    const FREEZING: f64 = 0.0;

    let t_profile = snd.get_profile(var);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() {
        return Err(AnalysisError::MissingProfile(profile_name, analysis_name));
    }
    if p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile("Pressure", analysis_name));
    }

    let mut profile = t_profile.iter().zip(p_profile);

    let mut bottom_press = ::std::f64::MAX; // Only init because compiler can't tell value not used
    let mut top_press: f64;

    // Initialize the bottom of the sounding
    let mut last_t: f64;
    let mut last_press: f64;
    loop {
        if let Some((t, press)) = profile.by_ref().next() {
            if let (Some(t), Some(press)) = (*t, *press) {
                last_t = t;
                last_press = press;
                break;
            }
        } else {
            match var {
                Temperature => {
                    return Err(AnalysisError::NoDataProfile(
                        "Pressure and temperature",
                        analysis_name,
                    ))
                }
                WetBulb => {
                    return Err(AnalysisError::NoDataProfile(
                        "Pressure, and wet bulb",
                        analysis_name,
                    ))
                }
                _ => unreachable!(),
            }
        }
    }

    // Check to see if we are below freezing at the bottom
    if last_t > FREEZING {
        return Ok(to_return);
    }

    let mut in_warm_zone = false;

    for (t, press) in profile {
        if let (Some(t), Some(press)) = (*t, *press) {
            if press < 500.0 {
                break;
            }

            if last_t <= FREEZING && t > FREEZING {
                bottom_press =
                    ::interpolation::linear_interp(FREEZING, last_t, t, last_press, press);
                in_warm_zone = true;
            }
            // Crossed out of zone to warm side
            if last_t > FREEZING && t <= FREEZING {
                top_press = ::interpolation::linear_interp(FREEZING, last_t, t, last_press, press);

                let bottom = ::interpolation::linear_interpolate(snd, bottom_press);
                let top = ::interpolation::linear_interpolate(snd, top_press);
                to_return.push(Layer { bottom, top });
                in_warm_zone = false;
            }
            last_t = t;
            last_press = press;
        }
    }

    // Check to see if we ended in a warm layer aloft
    if last_t > FREEZING && in_warm_zone {
        top_press = last_press;
        let bottom = ::interpolation::linear_interpolate(snd, bottom_press);
        let top = ::interpolation::linear_interpolate(snd, top_press);
        to_return.push(Layer { bottom, top });
    }

    Ok(to_return)
}

/// Assuming a warm layer aloft given by warm_layers, measure the cold surface layer.
pub fn cold_surface_temperature_layer(snd: &Sounding, warm_layers: &[Layer]) -> Option<Layer> {
    cold_surface_layer(snd, Temperature, warm_layers)
}

fn cold_surface_layer(snd: &Sounding, var: Profile, warm_layers: &[Layer]) -> Option<Layer> {
    assert!(var == Temperature || var == WetBulb);

    const FREEZING: f64 = 0.0;

    if warm_layers.is_empty() {
        return None;
    }

    let t_profile = snd.get_profile(var);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() || p_profile.is_empty() {
        return None; // Should not happen since we already used these to get warm layer
    }

    let mut profile = t_profile.iter().zip(p_profile);

    let last_t: f64;
    let last_press: f64;
    loop {
        if let Some((t, press)) = profile.next() {
            if let (Some(t), Some(press)) = (*t, *press) {
                last_t = t;
                last_press = press;
                break;
            }
        } else {
            return None;
        }
    }

    // Check to see if we are below freezing at the bottom
    if last_t > FREEZING {
        return None;
    }

    let bottom = ::interpolation::linear_interpolate(snd, last_press);

    Some(Layer {
        bottom,
        top: warm_layers[0].bottom,
    })
}

/// Get a layer that has a certain thickness, like 3km or 6km.
pub fn layer_agl(snd: &Sounding, meters_agl: f64) -> Result<Layer, AnalysisError> {
    const ANALYSIS: &str = "Layer AGL";

    let tgt_elev = if let Some(elev) = snd.get_location().2 {
        elev + meters_agl
    } else {
        return Err(AnalysisError::MissingValue("station elevation", ANALYSIS));
    };

    let h_profile = snd.get_profile(GeopotentialHeight);
    let p_profile = snd.get_profile(Pressure);

    if h_profile.is_empty() {
        return Err(AnalysisError::MissingProfile(
            "Geopotential Height",
            ANALYSIS,
        ));
    }
    if p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile("Pressure", ANALYSIS));
    }

    let mut profile = h_profile.iter().zip(p_profile);

    // Initialize the bottom of the sounding
    let mut last_h: f64;
    let mut last_press: f64;
    loop {
        if let Some((h, press)) = profile.by_ref().next() {
            if let (Some(h), Some(press)) = (*h, *press) {
                last_h = h;
                last_press = press;
                break;
            }
        } else {
            return Err(AnalysisError::NoDataProfile(
                "Geopotential Height",
                ANALYSIS,
            ));
        }
    }

    // Check we aren't too high already.
    if last_h > tgt_elev {
        return Err(AnalysisError::NotEnoughData(ANALYSIS));
    }

    let bottom_press = last_press; // Lowest level we could find.
    let bottom = ::interpolation::linear_interpolate(snd, bottom_press);

    for (h, press) in profile {
        if let (Some(h), Some(press)) = (*h, *press) {
            if last_h <= tgt_elev && h > tgt_elev {
                let top_press =
                    ::interpolation::linear_interp(tgt_elev, last_h, h, last_press, press);

                let top = ::interpolation::linear_interpolate(snd, top_press);

                return Ok(Layer { bottom, top });
            }
            last_h = h;
            last_press = press;
        }
    }

    Err(AnalysisError::NotEnoughData(ANALYSIS))
}

// TODO: Inversions.
