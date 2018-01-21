use sounding_base::{Profile, Sounding};
use sounding_base::Profile::*;

/// Find the dendtritic growth zones throughout the profile. It is unusual, but possible there is
/// more than 1!
pub fn dendritic_growth_zone(snd: &Sounding, v_coord: Profile) -> Vec<(f64, f64)> {

    // FIXME: what happens if there is no profile for T or v_coord? Should return empty Vec.

    let mut result = Vec::with_capacity(2);

    // Dendritic snow growth zone temperature range in C
    const WARM_SIDE: f64 = -12.0;
    const COLD_SIDE: f64 = -18.0;

    let mut profile = snd.get_profile(Temperature)
        .iter()
        .zip(snd.get_profile(v_coord).iter());

    let mut bottom = ::std::f64::MAX; // Only initialize because compiler can't tell value not used
    let mut top: f64;

    // Initialize the bottom of the sounding
    let mut last_t: f64;
    let mut last_coord: f64;
    loop {
        if let Some((t, coord)) = profile.by_ref().next() {
            if let (Some(t), Some(coord)) = (*t, *coord) {
                last_t = t;
                last_coord = coord;
                break;
            }
        }
    }

    // Check to see if we are already in the dendtritic zone
    if last_t <= WARM_SIDE && last_t >= COLD_SIDE {
        bottom = last_coord;
    }

    for (t, coord) in profile {
        if let (Some(t), Some(coord)) = (*t, *coord) {
            // Crossed into zone from warm side
            if last_t > WARM_SIDE && t <= WARM_SIDE {
                bottom = ::interpolation::linear_interp(WARM_SIDE, last_t, t, last_coord, coord);
            }
            // Crossed into zone from cold side
            if last_t < COLD_SIDE && t >= COLD_SIDE {
                bottom = ::interpolation::linear_interp(COLD_SIDE, last_t, t, last_coord, coord);
            }
            // Crossed out of zone to warm side
            if last_t <= WARM_SIDE && t > WARM_SIDE {
                top = ::interpolation::linear_interp(WARM_SIDE, last_t, t, last_coord, coord);
                result.push((bottom, top));
            }
            // Crossed out of zone to cold side
            if last_t >= COLD_SIDE && t < COLD_SIDE {
                top = ::interpolation::linear_interp(COLD_SIDE, last_t, t, last_coord, coord);
                result.push((bottom, top));
            }
            last_t = t;
            last_coord = coord;
        }
    }

    // Check to see if we ended in a dendtritic zone
    if last_t <= WARM_SIDE && last_t >= COLD_SIDE {
        top = last_coord;
        result.push((bottom, top));
    }

    result
}

// TODO: Wet bulb zero height given v_coord Return multiple if needed.
// TODO: Freezing level given v_coord. Return multiple if needed.
// TODO: Metlting layer
// TODO: Warm layer aloft.
// TODO: Cold layer at surface.
