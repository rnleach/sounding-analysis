use super::{Layer, Layers};
use crate::{
    error::{
        AnalysisError::{InvalidInput, MissingProfile, NoDataProfile},
        {AnalysisError, Result},
    },
    sounding::Sounding,
};
use itertools::izip;
use metfor::{Celsius, HectoPascal, FREEZING};
use optional::Optioned;

/// Find the dendtritic growth zones throughout the profile. It is unusual, but possible there is
/// more than one.
///
/// If there are none, then an empty vector is returned.
pub fn dendritic_snow_zone(snd: &Sounding) -> Result<Layers> {
    temperature_layer(snd, Celsius(-12.0), Celsius(-18.0), HectoPascal(300.0))
}

/// Find the hail growth zones throughout the profile. It is very unusual, but possible there is
/// more than one.
///
/// If there are none, then an empty vector is returned.
pub fn hail_growth_zone(snd: &Sounding) -> Result<Layers> {
    temperature_layer(snd, Celsius(-10.0), Celsius(-30.0), HectoPascal(1.0))
}

fn temperature_layer(
    snd: &Sounding,
    warm_side: Celsius,
    cold_side: Celsius,
    top_pressure: HectoPascal,
) -> Result<Layers> {
    use crate::interpolation::{linear_interp, linear_interpolate_sounding};
    let mut to_return: Layers = Layers::new();

    let t_profile = snd.temperature_profile();
    let p_profile = snd.pressure_profile();

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    izip!(p_profile, t_profile)
        // remove levels with missing values
        .filter_map(|pair| {
            if pair.0.is_some() && pair.1.is_some() {
                let (p, t) = (pair.0.unpack(), pair.1.unpack());
                Some((p, t))
            } else {
                None
            }
        })
        // Stop above a certain level
        .take_while(|&(p, _)| p > top_pressure)
        // find temperature layers
        .fold(
            Ok((None, None, None)),
            |acc: Result<(Option<HectoPascal>, Option<Celsius>, Option<_>)>, (p, t)| {
                match acc {
                    // We're not in a target layer currently
                    Ok((Some(last_p), Some(last_t), None)) => {
                        if last_t < cold_side && t >= cold_side && t <= warm_side {
                            // We crossed into a target layer from the cold side
                            let target_p = linear_interp(cold_side, last_t, t, last_p, p);
                            let bottom = linear_interpolate_sounding(snd, target_p)?;
                            Ok((Some(p), Some(t), Some(bottom)))
                        } else if last_t > warm_side && t >= cold_side && t <= warm_side {
                            // We crossed into a target layer from the warm side
                            let target_p = linear_interp(warm_side, last_t, t, last_p, p);
                            let bottom = linear_interpolate_sounding(snd, target_p)?;
                            Ok((Some(p), Some(t), Some(bottom)))
                        } else if (last_t < cold_side && t >= warm_side)
                            || (last_t > warm_side && t <= cold_side)
                        {
                            // We crossed completely through a target layer
                            let warm_p = linear_interp(warm_side, last_t, t, last_p, p);
                            let cold_p = linear_interp(cold_side, last_t, t, last_p, p);
                            let bottom = linear_interpolate_sounding(snd, warm_p.max(cold_p))?;
                            let top = linear_interpolate_sounding(snd, warm_p.min(cold_p))?;
                            to_return.push(Layer { bottom, top });
                            Ok((Some(p), Some(t), None))
                        } else {
                            // We weren't in a target layer
                            Ok((Some(p), Some(t), None))
                        }
                    }

                    // We're in a target layer, let's see if we passed out
                    Ok((Some(last_p), Some(last_t), Some(bottom))) => {
                        if t < cold_side {
                            // We crossed out of a target layer on the cold side
                            let target_p = linear_interp(cold_side, last_t, t, last_p, p);
                            let top = linear_interpolate_sounding(snd, target_p)?;
                            to_return.push(Layer { bottom, top });
                            Ok((Some(p), Some(t), None))
                        } else if t > warm_side {
                            // We crossed out of a target layer on the warm side
                            let target_p = linear_interp(warm_side, last_t, t, last_p, p);
                            let top = linear_interpolate_sounding(snd, target_p)?;
                            to_return.push(Layer { bottom, top });
                            Ok((Some(p), Some(t), None))
                        } else {
                            // We're still in a target layer
                            Ok((Some(p), Some(t), Some(bottom)))
                        }
                    }

                    // Propagate errors
                    e @ Err(_) => e,

                    // First row, lets get started
                    Ok((None, None, None)) => {
                        if t <= warm_side && t >= cold_side {
                            // Starting out in a target layer
                            let dr = linear_interpolate_sounding(snd, p)?;
                            Ok((Some(p), Some(t), Some(dr)))
                        } else {
                            // Not starting out in a target layer
                            Ok((Some(p), Some(t), None))
                        }
                    }

                    // No other combinations are possible
                    _ => unreachable!(),
                }
            },
        )
        // Swap my list into the result.
        .and_then(|_| Ok(to_return))
}

/// Assuming it is below freezing at the surface, this will find the warm layers aloft using the
/// dry bulb temperature. Does not look above 500 hPa.
pub fn warm_temperature_layer_aloft(snd: &Sounding) -> Result<Layers> {
    warm_layer_aloft(snd, snd.temperature_profile())
}

/// Assuming the wet bulb temperature is below freezing at the surface, this will find the warm
/// layers aloft using the wet bulb temperature. Does not look above 500 hPa.
pub fn warm_wet_bulb_layer_aloft(snd: &Sounding) -> Result<Layers> {
    warm_layer_aloft(snd, snd.wet_bulb_profile())
}

fn warm_layer_aloft(snd: &Sounding, t_profile: &[Optioned<Celsius>]) -> Result<Layers> {
    let mut to_return: Layers = Layers::new();

    let p_profile = snd.pressure_profile();

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    izip!(p_profile, t_profile)
        // Remove levels without pressure AND temperature data
        .filter_map(|pair| {
            if pair.0.is_some() && pair.1.is_some() {
                let (p, t) = (pair.0.unpack(), pair.1.unpack());
                Some((p, t))
            } else {
                None
            }
        })
        // Ignore anything above 500 hPa, extremely unlikely for a warm layer up there.
        .take_while(|&(p, _)| p > HectoPascal(500.0))
        // Find the warm layers!
        .fold(
            Ok((HectoPascal(std::f64::MAX), Celsius(std::f64::MAX), None)),
            |last_iter_res: Result<(_, _, _)>, (p, t)| {
                let (last_p, last_t, mut bottom) = last_iter_res?;
                if last_t <= FREEZING && t > FREEZING && bottom.is_none() {
                    // Entering a warm layer.
                    let bottom_p =
                        crate::interpolation::linear_interp(FREEZING, last_t, t, last_p, p);
                    bottom = Some(crate::interpolation::linear_interpolate_sounding(
                        snd, bottom_p,
                    )?);
                }
                if bottom.is_some() && last_t > FREEZING && t <= FREEZING {
                    // Crossed out of a warm layer
                    let top_p = crate::interpolation::linear_interp(FREEZING, last_t, t, last_p, p);
                    let top = crate::interpolation::linear_interpolate_sounding(snd, top_p)?;
                    {
                        let bottom = bottom.unwrap();
                        to_return.push(Layer { bottom, top });
                    }
                    bottom = None;
                }

                Ok((p, t, bottom))
            },
        )?;

    Ok(to_return)
}

/// Assuming a warm layer aloft given by `warm_layers`, measure the cold surface layer.
pub fn cold_surface_temperature_layer(snd: &Sounding, warm_layers: &[Layer]) -> Result<Layer> {
    cold_surface_layer(snd, snd.temperature_profile(), warm_layers)
}

fn cold_surface_layer(
    snd: &Sounding,
    t_profile: &[Optioned<Celsius>],
    warm_layers: &[Layer],
) -> Result<Layer> {
    if warm_layers.is_empty() {
        return Err(InvalidInput);
    }

    let p_profile = snd.pressure_profile();

    if t_profile.is_empty() || p_profile.is_empty() {
        // Should not happen since we SHOULD HAVE already used these to get the warm layers
        return Err(MissingProfile);
    }

    izip!(0usize.., p_profile, t_profile)
        // Remove levels with missing data
        .filter_map(|triplet| {
            if triplet.1.is_some() && triplet.2.is_some() {
                let (i, t) = (triplet.0, triplet.2.unpack());
                Some((i, t))
            } else {
                None
            }
        })
        // Map it to an error if the temperature is above freezing.
        .map(|(i, t)| {
            if t > FREEZING {
                Err(InvalidInput)
            } else {
                Ok(i)
            }
        })
        // Only take the first one, we want the surface layer, or the lowest layer available
        .next()
        // If there is nothing to get, there was no valid data in the profile.
        .unwrap_or(Err(NoDataProfile))
        // Map the result into a data row!
        .and_then(|index| snd.data_row(index).ok_or(InvalidInput))
        // Package it up in a layer
        .map(|bottom| Layer {
            bottom,
            top: warm_layers[0].bottom,
        })
}
