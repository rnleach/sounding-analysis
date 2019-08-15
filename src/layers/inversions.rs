use super::{Layer, Layers};
use crate::{
    error::{AnalysisError, Result},
    sounding::Sounding,
};
use itertools::izip;
use metfor::{Celsius, HectoPascal};

/// Get all inversion layers up to a specified pressure.
pub fn inversions(snd: &Sounding, top_p: HectoPascal) -> Result<Layers> {
    let mut to_return: Layers = Layers::new();

    let t_profile = snd.temperature_profile();
    let p_profile = snd.pressure_profile();

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    izip!(0usize.., p_profile, t_profile)
        // Filter out rows without both temperature and pressure.
        .filter_map(|triple| {
            if triple.1.is_some() && triple.2.is_some() {
                let (i, p, t) = (triple.0, triple.1.unpack(), triple.2.unpack());
                Some((i, p, t))
            } else {
                None
            }
        })
        // Filter out rows above the top pressure
        .filter_map(|(i, p, t)| if p < top_p { None } else { Some((i, t)) })
        // Capture the inversion layers
        .fold(
            (0, Celsius(std::f64::MAX), None),
            |(last_i, last_t, mut bottom_opt), (i, t)| {
                if bottom_opt.is_none() && last_t < t {
                    // Coming into an inversion
                    bottom_opt = snd.data_row(last_i);
                } else if bottom_opt.is_some() && last_t > t {
                    // Leaving an inversion
                    if let Some(layer) = bottom_opt.and_then(|bottom| {
                        snd.data_row(last_i)
                            .and_then(|top| Some(Layer { bottom, top }))
                    }) {
                        to_return.push(layer);
                        bottom_opt = None;
                    }
                }
                (i, t, bottom_opt)
            },
        );

    Ok(to_return)
}

/// Get a surface based inversion.
pub fn sfc_based_inversion(snd: &Sounding) -> Result<Option<Layer>> {
    let t_profile = snd.temperature_profile();
    let p_profile = snd.pressure_profile();
    let h_profile = snd.height_profile();

    if t_profile.is_empty() || p_profile.is_empty() || h_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    izip!(0usize.., p_profile, h_profile, t_profile)
        // Remove levels with missing data
        .filter_map(|tuple| {
            if tuple.1.is_some() && tuple.2.is_some() && tuple.2.is_some() {
                Some(tuple.0)
            } else {
                None
            }
        })
        // Get the first one.
        .nth(0)
        // Map Option to Result
        .ok_or(AnalysisError::NotEnoughData)
        // Map the result into a data row
        .and_then(|index| snd.data_row(index).ok_or(AnalysisError::MissingValue))
        // Now find the top
        .and_then(|bottom_row| {
            let sfc_t = bottom_row.temperature;
            if sfc_t.is_some() {
                let sfc_t = sfc_t.unpack();
                let val = izip!(0usize.., p_profile, t_profile, h_profile)
                    // Remove levels with missing data
                    .filter_map(|tuple| {
                        if tuple.1.is_some() && tuple.2.is_some() && tuple.3.is_some() {
                            let (i, p, t) = (tuple.0, tuple.1.unpack(), tuple.2.unpack());
                            Some((i, p, t))
                        } else {
                            None
                        }
                    })
                    // This is the first one!
                    .skip(1)
                    // Only look up to about 700 hPa
                    .take_while(|(_, p, _)| *p > HectoPascal(690.0))
                    // Remove those cooler than the surface
                    .filter(|(_, _, t)| *t > sfc_t)
                    .fold(None, |max_t_info, (i, _, t)| {
                        if let Some((max_t, _max_t_idx)) = max_t_info {
                            if t > max_t {
                                Some((t, i))
                            } else {
                                max_t_info
                            }
                        } else {
                            Some((t, i))
                        }
                    });

                match val {
                    Some((_, idx)) => snd
                        .data_row(idx)
                        .ok_or(AnalysisError::MissingValue)
                        .and_then(|top_row| {
                            Ok(Some(Layer {
                                bottom: bottom_row,
                                top: top_row,
                            }))
                        }),
                    None => Ok(None),
                }
            } else {
                Err(AnalysisError::MissingValue)
            }
        })
}
