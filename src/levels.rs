//! This module finds significant levels such as the freezing level and wet bulb zero level. It also
//! has functions for finding critical values at a single level, such as the maximum wet bulb
//! temperature aloft.  It does not include functions for finding levels related to parcel analysis
//! and convection, those are found in the `parcel` module.
use crate::{
    error::{
        AnalysisError::{InvalidInput, MissingProfile, NoDataProfile, NotEnoughData},
        {AnalysisError, Result},
    },
    layers::Layer,
    sounding::{DataRow, Sounding},
};
use itertools::izip;
use metfor::{Celsius, HectoPascal, Meters, FREEZING};
use optional::Optioned;
use smallvec::SmallVec;

/// A level in the atmosphere is described by a `DataRow` from a sounding.
pub type Level = DataRow;

/// A list of levels.
pub type Levels = SmallVec<[Level; crate::VEC_SIZE]>;

/// Find the freezing/melting levels below 500 hPa.
pub fn freezing_levels(snd: &Sounding) -> Result<Levels> {
    find_temperature_levels(
        FREEZING,
        snd.pressure_profile(),
        snd.temperature_profile(),
        snd,
    )
}

/// Find the wet bulb zero levels
pub fn wet_bulb_zero_levels(snd: &Sounding) -> Result<Levels> {
    find_temperature_levels(
        FREEZING,
        snd.pressure_profile(),
        snd.wet_bulb_profile(),
        snd,
    )
}

fn find_temperature_levels(
    target_t: Celsius,
    p_profile: &[Optioned<HectoPascal>],
    t_profile: &[Optioned<Celsius>],
    snd: &Sounding,
) -> Result<Levels> {
    use crate::interpolation::{linear_interp, linear_interpolate_sounding};

    let mut to_return: Levels = Levels::new();

    const TOP_PRESSURE: HectoPascal = HectoPascal(500.0); // don't look above here.

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
            |acc: Result<(HectoPascal, Celsius)>, (p, t)| {
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
    max_t_aloft(snd, snd.wet_bulb_profile())
}

/// Maximum temperature aloft.
pub fn max_temperature_in_profile(snd: &Sounding) -> Result<Level> {
    max_t_aloft(snd, snd.temperature_profile())
}

// Only searches up to 500 hPa
fn max_t_aloft(snd: &Sounding, t_profile: &[Optioned<Celsius>]) -> Result<Level> {
    const TOP_PRESSURE: HectoPascal = HectoPascal(500.0); // don't look above here.

    let p_profile = snd.pressure_profile();

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
        })
        .fold(
            Err(AnalysisError::NotEnoughData),
            |acc: Result<_>, (i, t)| {
                if let Ok((_, mx_t)) = acc {
                    if t > mx_t {
                        Ok((i, t))
                    } else {
                        // Propagate most recent result through
                        acc
                    }
                } else {
                    // just starting, so initialize result
                    Ok((i, t))
                }
            },
        )
        // Retrive the row
        .and_then(|(idx, _)| snd.data_row(idx).ok_or(InvalidInput))
}

/// Maximum temperature in a layer.
pub fn max_temperature_in_layer(snd: &Sounding, lyr: &Layer) -> Result<Level> {
    max_t_in_layer(snd, snd.temperature_profile(), lyr)
}

/// Maximum wet bulb temperature in a layer.
pub fn max_wet_bulb_in_layer(snd: &Sounding, lyr: &Layer) -> Result<Level> {
    max_t_in_layer(snd, snd.wet_bulb_profile(), lyr)
}

fn max_t_in_layer(snd: &Sounding, t_profile: &[Optioned<Celsius>], lyr: &Layer) -> Result<Level> {
    let (bottom_p, top_p) = if lyr.bottom.pressure.is_some() && lyr.top.pressure.is_some() {
        (lyr.bottom.pressure.unpack(), lyr.top.pressure.unpack())
    } else {
        return Err(AnalysisError::InvalidInput);
    };

    let p_profile = snd.pressure_profile();

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
        })
        .fold(
            Err(AnalysisError::NotEnoughData),
            |acc: Result<_>, (i, t)| {
                if let Ok((_, mx_t)) = acc {
                    if t > mx_t {
                        Ok((i, t))
                    } else {
                        // Propagate most recent result through
                        acc
                    }
                } else {
                    // We're just starting, so populate the result
                    Ok((i, t))
                }
            },
        )
        // Retrive the row
        .and_then(|(idx, _)| snd.data_row(idx).ok_or(InvalidInput))
}

/// Find a level at a specific geopotential height.
pub(crate) fn height_level(tgt_height: Meters, snd: &Sounding) -> Result<Level> {
    let h_profile = snd.height_profile();
    let p_profile = snd.pressure_profile();

    if h_profile.is_empty() || p_profile.is_empty() {
        return Err(MissingProfile);
    }

    izip!(p_profile, h_profile)
        // filter out levels with missing data
        .filter_map(|pair| {
            if pair.0.is_some() && pair.1.is_some() {
                let (p, h) = (pair.0.unpack(), pair.1.unpack());
                Some((p, h))
            } else {
                None
            }
        })
        // find the pressure at the target geopotential height, to be used later for interpolation.
        .fold(
            Ok((HectoPascal(std::f64::MAX), Meters(0.0f64), None)),
            |acc: Result<(_, _, Option<_>)>, (p, h)| {
                match acc {
                    // We have not yet found the target pressure to interpolate everything to, so
                    // check the current values.
                    Ok((last_p, last_h, None)) => {
                        if h > tgt_height {
                            // If we finally jumped above our target, we have it bracketed, interpolate
                            // and find target pressure.
                            let tgt_p = crate::interpolation::linear_interp(
                                tgt_height, last_h, h, last_p, p,
                            );
                            Ok((
                                HectoPascal(std::f64::MAX),
                                Meters(std::f64::MAX),
                                Some(tgt_p),
                            ))
                        } else {
                            // Keep climbing up the profile.
                            Ok((p, h, None))
                        }
                    }
                    // We have found the target pressure on the last iteration, pass it through
                    ok @ Ok((_, _, Some(_))) => ok,
                    // There was an error, keep passing it through.
                    e @ Err(_) => e,
                }
            },
        )
        // Extract the target pressure
        .and_then(|(_, _, opt)| opt.ok_or(NotEnoughData))
        // Do the interpolation.
        .and_then(|target_p| crate::interpolation::linear_interpolate_sounding(snd, target_p))
}
