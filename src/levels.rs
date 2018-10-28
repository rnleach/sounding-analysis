//! This module finds significant levels such as the freezing level and wet bulb zero level. It also
//! has functions for finding critical values at a single level, such as the maximum wet bulb
//! temperature aloft.  It does not include functions for finding levels related to parcel analysis
//! and convection, those are found in the `parcel` module.

use smallvec::SmallVec;

use sounding_base::Profile::*;
use sounding_base::{DataRow, Profile, Sounding};

use error::AnalysisError::*;
use error::*;

use layers::Layer;

const FREEZING: f64 = 0.0;

/// A level in the atmosphere is described by a `DataRow` from a sounding.
pub type Level = DataRow;

/// A list of levels.
pub type Levels = SmallVec<[Level; ::VEC_SIZE]>;

/// Find the freezing/melting levels below 500 hPa.
pub fn freezing_levels(snd: &Sounding) -> Result<Levels> {
    find_temperature_levels(snd, Temperature, FREEZING)
}

/// Find the wet bulb zero levels
pub fn wet_bulb_zero_levels(snd: &Sounding) -> Result<Levels> {
    find_temperature_levels(snd, WetBulb, FREEZING)
}

fn find_temperature_levels(snd: &Sounding, var: Profile, target_t: f64) -> Result<Levels> {
    use interpolation::{linear_interp, linear_interpolate_sounding};

    debug_assert!(var == Temperature || var == WetBulb);

    let mut to_return: Levels = Levels::new();

    const TOP_PRESSURE: f64 = 500.0; // don't look above here.

    let p_profile = snd.get_profile(Pressure);
    let t_profile = snd.get_profile(var);

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    let mut iter = izip!(p_profile, t_profile).filter_map(|pair| {
        if pair.0.is_some() && pair.1.is_some() {
            let (p, t) = (pair.0.unpack(), pair.1.unpack());
            Some((p, t))
        } else {
            None
        }
    });

    let (bottom_p, bottom_t) = iter.by_ref().next().ok_or(NoDataProfile)?;

    iter
        // Don't bother looking above a certain level
        .take_while(|&(p, _)| p >= TOP_PRESSURE)
        // Reduce to get the temperature levels
        .fold(
            Ok((bottom_p, bottom_t)),
            |acc: Result<(f64, f64)>, (p, t)| {
                if let Ok((last_p, last_t)) = acc {
                    if last_t <= target_t && t > target_t || last_t > target_t && t <= target_t {
                        let target_p = linear_interp(target_t, last_t, t, last_p, p);
                        to_return.push(linear_interpolate_sounding(snd, target_p)?);
                    }
                    Ok((p, t))
                } else {
                    // Pass the error through
                    acc
                }
            },
        )
        // Swap my vector into the result
        .map(|_| to_return)
}

/// Maximum wet bulb temperature aloft.
pub fn max_wet_bulb_in_profile(snd: &Sounding) -> Result<Level> {
    max_t_aloft(snd, Profile::WetBulb)
}

/// Maximum temperature aloft.
pub fn max_temperature_in_profile(snd: &Sounding) -> Result<Level> {
    max_t_aloft(snd, Profile::Temperature)
}

// Only searches up to 500 hPa
fn max_t_aloft(snd: &Sounding, var: Profile) -> Result<Level> {
    use sounding_base::Profile::*;

    debug_assert!(var == Temperature || var == WetBulb);

    const TOP_PRESSURE: f64 = 500.0; // don't look above here.

    let p_profile = snd.get_profile(Pressure);
    let t_profile = snd.get_profile(var);

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    izip!(0usize.., p_profile, t_profile)
        .filter_map(|triple| {
            if triple.1.is_some() && triple.2.is_some() {
                let (i, p, t) = (triple.0, triple.1.unpack(), triple.2.unpack());
                if p >= TOP_PRESSURE {
                    Some((i, t))
                } else {
                    None
                }
            } else {
                None
            }
        }).fold(Ok((0, ::std::f64::MIN)), |acc: Result<_>, (i, t)| {
            if let Ok((_, mx_t)) = acc {
                if t > mx_t {
                    Ok((i, t))
                } else {
                    // Propagate most recent result through
                    acc
                }
            } else {
                // Propagate errors
                acc
            }
        })
        // Retrive the row
        .and_then(|(idx, _)| snd.get_data_row(idx).ok_or(InvalidInput))
}

/// Maximum temperature in a layer.
pub fn max_temperature_in_layer(snd: &Sounding, lyr: &Layer) -> Result<Level> {
    max_t_in_layer(snd, Profile::Temperature, lyr)
}

/// Maximum wet bulb temperature in a layer.
pub fn max_wet_bulb_in_layer(snd: &Sounding, lyr: &Layer) -> Result<Level> {
    max_t_in_layer(snd, Profile::WetBulb, lyr)
}

fn max_t_in_layer(snd: &Sounding, var: Profile, lyr: &Layer) -> Result<Level> {
    use sounding_base::Profile::*;

    debug_assert!(var == Temperature || var == WetBulb);

    let (bottom_p, top_p) = if lyr.bottom.pressure.is_some() && lyr.top.pressure.is_some() {
        (lyr.bottom.pressure.unpack(), lyr.top.pressure.unpack())
    } else {
        return Err(AnalysisError::InvalidInput);
    };

    let p_profile = snd.get_profile(Pressure);
    let t_profile = snd.get_profile(var);

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    izip!(0usize.., p_profile, t_profile)
        .filter_map(|triple| {
            if triple.1.is_some() && triple.2.is_some() {
                let (i, p, t) = (triple.0, triple.1.unpack(), triple.2.unpack());
                if p >= top_p && p <= bottom_p {
                    Some((i, t))
                } else {
                    None
                }
            } else {
                None
            }
        }).fold(Ok((0, ::std::f64::MIN)), |acc: Result<_>, (i, t)| {
            if let Ok((_, mx_t)) = acc {
                if t > mx_t {
                    Ok((i, t))
                } else {
                    // Propagate most recent result through
                    acc
                }
            } else {
                // Propagate errors
                acc
            }
        })
        // Retrive the row
        .and_then(|(idx, _)| snd.get_data_row(idx).ok_or(InvalidInput))
}
