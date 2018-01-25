//! This module finds significant layers such as the dendritic snow growth zone, the hail growth
//! zone, and inversions.

use error::*;
use smallvec::SmallVec;

use sounding_base::Sounding;
use sounding_base::Profile::*;

use Layer;

/// Find the dendtritic growth zones throughout the profile. It is unusual, but possible there is
/// more than 1.
///
/// # Errors
/// If the sounding is missing a temperature or pressure profile, `error::ErrorKind::MissingProfile`
/// is returned in the result. Otherwise, if no dendritic layers are found, an empty vector is
/// returned in the `Result`
pub fn dendritic_snow_zone(snd: &Sounding) -> Result<SmallVec<[Layer; ::VEC_SIZE]>> {
    const ANALYSIS_NAME: &str = "Dendritic Snow Growth Zone(s)";

    let mut to_return: SmallVec<[Layer; ::VEC_SIZE]> = SmallVec::new();

    // Dendritic snow growth zone temperature range in C
    const WARM_SIDE: f64 = -12.0;
    const COLD_SIDE: f64 = -18.0;

    let t_profile = snd.get_profile(Temperature);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() {
        bail!(ErrorKind::MissingProfile("Temperature", ANALYSIS_NAME));
    }
    if p_profile.is_empty() {
        bail!(ErrorKind::MissingProfile("Pressure", ANALYSIS_NAME));
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
        }
    }

    // Check to see if we are already in the dendtritic zone
    if last_t <= WARM_SIDE && last_t >= COLD_SIDE {
        bottom_press = last_press;
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
                to_return.push(Layer {
                    bottom_press,
                    top_press,
                });
            }
            // Crossed out of zone to cold side
            if last_t >= COLD_SIDE && t < COLD_SIDE {
                top_press = ::interpolation::linear_interp(COLD_SIDE, last_t, t, last_press, press);
                to_return.push(Layer {
                    bottom_press,
                    top_press,
                });
            }
            last_t = t;
            last_press = press;
        }
    }

    // Check to see if we ended in a dendtritic zone
    if last_t <= WARM_SIDE && last_t >= COLD_SIDE {
        top_press = last_press;
        to_return.push(Layer {
            bottom_press,
            top_press,
        });
    }

    Ok(to_return)
}

/// Assuming it is below freezing at the surface, this will find the warm layers aloft using the
/// dry bulb temperature. Does not look above 500 hPa.
pub fn warm_temperature_layer_aloft(snd: &Sounding) -> Result<SmallVec<[Layer; ::VEC_SIZE]>> {
    const ANALYSIS_NAME: &str = "Warm temperature layer aloft.";

    let mut to_return: SmallVec<[Layer; ::VEC_SIZE]> = SmallVec::new();

    // Dendritic snow growth zone temperature range in C
    const FREEZING: f64 = 0.0;

    let t_profile = snd.get_profile(Temperature);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() {
        bail!(ErrorKind::MissingProfile("Temperature", ANALYSIS_NAME));
    }
    if p_profile.is_empty() {
        bail!(ErrorKind::MissingProfile("Pressure", ANALYSIS_NAME));
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
                to_return.push(Layer {
                    bottom_press,
                    top_press,
                });
                in_warm_zone = false;
            }
            last_t = t;
            last_press = press;
        }
    }

    // Check to see if we ended in a warm layer aloft
    if last_t > FREEZING && in_warm_zone {
        top_press = last_press;
        to_return.push(Layer {
            bottom_press,
            top_press,
        });
    }

    Ok(to_return)
}


// TODO: Cold layer at surface.
// TODO: Warm wet bulb layer aloft.
