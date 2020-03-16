//! This module finds significant levels such as the freezing level and wet bulb zero level. It also
//! has functions for finding critical values at a single level, such as the maximum wet bulb
//! temperature aloft.  It does not include functions for finding levels related to parcel analysis
//! and convection, those are found in the `parcel` module.
use crate::{
    error::{
        AnalysisError::{InvalidInput, MissingProfile},
        {AnalysisError, Result},
    },
    interpolation::{linear_interp, linear_interpolate_sounding},
    layers::Layer,
    sounding::{DataRow, Sounding},
};
use itertools::{izip, Itertools};
use metfor::{Celsius, HectoPascal, Meters, FREEZING};
use optional::Optioned;

/// A level in the atmosphere is described by a `DataRow` from a sounding.
pub type Level = DataRow;

/// A list of levels.
pub type Levels = Vec<Level>;

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
    const TOP_PRESSURE: HectoPascal = HectoPascal(500.0); // don't look above here

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    #[inline]
    fn crosses_target(t0: Celsius, t1: Celsius, target_t: Celsius) -> bool {
        (t0 <= target_t && t1 >= target_t) || (t0 >= target_t && t1 <= target_t)
    };

    izip!(p_profile, t_profile)
        // Remove levels with missing data
        .filter(|(p, t)| p.is_some() && t.is_some())
        // Unwrap from the Optioned type
        .map(|(p, t)| (p.unpack(), t.unpack()))
        // Don't bother looking above a certain level
        .take_while(|&(p, _)| p >= TOP_PRESSURE)
        // Look at them in pairs to find crossing the threshold
        .tuple_windows::<(_, _)>()
        // Filter out any pair that is not crossing the desired temperature level
        .filter(|&((_p0, t0), (_p1, t1))| crosses_target(t0, t1, target_t))
        // Interpolate crossings to get pressures values at the levels we want.
        .map(|((p0, t0), (p1, t1))| linear_interp(target_t, t0, t1, p0, p1))
        // Perform the interpolation on the full sounding to get the desired DataRow
        .map(|p_level| linear_interpolate_sounding(snd, p_level))
        .collect()
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
        // Filter out levels with missing values
        .filter(|(_, p, t)| p.is_some() && t.is_some())
        // unwrap the Optioned type
        .map(|(i, p, t)| (i, p.unpack(), t.unpack()))
        // Only look up to a certain level
        .take_while(|&(_, p, _)| p >= TOP_PRESSURE)
        // fold to get the maximum value
        .fold(
            Err(AnalysisError::NotEnoughData),
            |acc: Result<_>, (i, _, t)| {
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

#[inline]
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
        // Remove levels with missing values
        .filter(|(_, p, t)| p.is_some() && t.is_some())
        // Unwrap the Optioned type
        .map(|(i, p, t)| (i, p.unpack(), t.unpack()))
        // Skip values below the bottom of the layer
        .skip_while(|(_, p, _)| *p > bottom_p)
        // Only look up to the top of the layer
        .take_while(|(_, p, _)| *p >= top_p)
        // fold to find the max value
        .fold(
            Err(AnalysisError::NotEnoughData),
            |acc: Result<_>, (i, _, t)| {
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

    // Assumes h0 < h1, ie iterate from bottom to top of sounding
    #[inline]
    fn is_cross_over(h0: Meters, h1: Meters, tgt_height: Meters) -> bool {
        debug_assert!(h0 <= h1);
        h0 <= tgt_height && h1 >= tgt_height
    }

    izip!(p_profile, h_profile)
        // Filter out levels with missing data
        .filter(|(p, h)| p.is_some() && h.is_some())
        // Unpack from the Optioned type
        .map(|(p, h)| (p.unpack(), h.unpack()))
        // look at levels in pairs to find when the data crosses the desired height
        .tuple_windows::<(_, _)>()
        // Filter out all pairs that don't have a crossover
        .find(|&((_p0, h0), (_p1, h1))| is_cross_over(h0, h1, tgt_height))
        // Map the bracket into the pressure at the target height via interpolation
        .map(|((p0, h0), (p1, h1))| linear_interp(tgt_height, h0, h1, p0, p1))
        // map the option into a result
        .ok_or(AnalysisError::InvalidInput)
        // Interpolate the result
        .and_then(|level_p| linear_interpolate_sounding(snd, level_p))
}
