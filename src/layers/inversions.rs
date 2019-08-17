use super::{Layer, Layers};
use crate::{
    error::{AnalysisError, Result},
    sounding::Sounding,
};
use itertools::{izip, Itertools};
use metfor::HectoPascal;

/// Get all inversion layers up to a specified pressure.
pub fn inversions(snd: &Sounding, top_p: HectoPascal) -> Result<Layers> {
    let t_profile = snd.temperature_profile();
    let p_profile = snd.pressure_profile();

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    izip!(0usize.., p_profile, t_profile)
        // Filter out levels with missing data
        .filter(|&(_, p, t)| p.is_some() && t.is_some())
        // Unpack from the `Optioned` type
        .map(|(i, p, t)| (i, p.unpack(), t.unpack()))
        // Stop once we reach the top level pressure
        .take_while(|&(_, p, _)| p >= top_p)
        // Now all we need is the row index and the temperature
        .map(|(i, _, t)| (i, t))
        // Inspect the levels in pairs
        .tuple_windows::<(_, _)>()
        // Map into inversion types
        .map(|((i0, t0), (i1, t1))| {
            if t1 > t0 {
                LayerType::Inversion(i0, i1)
            } else {
                LayerType::NonInversion
            }
        })
        // Combine adjacent inversion layers into a single, larger inversion
        .scan(
            (None, None),
            |layers: &mut (Option<usize>, Option<usize>), layer: LayerType| {
                let bottom: &mut Option<usize> = &mut layers.0;
                let top: &mut Option<usize> = &mut layers.1;

                match layer {
                    LayerType::Inversion(i0, i1) => {
                        if bottom.is_none() {
                            *bottom = Some(i0);
                        }
                        *top = Some(i1);
                        Some(None)
                    }
                    LayerType::NonInversion => {
                        if let (Some(i0), Some(i1)) = (bottom.take(), top.take()) {
                            Some(Some((i0, i1)))
                        } else {
                            Some(None)
                        }
                    }
                }
            },
        )
        // Remove iterations that did not capture an inverison
        .filter_map(|opt| opt)
        // Map indexes into layers
        .map(|(i0, i1)| {
            snd.data_row(i0)
                .and_then(|bottom| snd.data_row(i1).map(|top| Layer { bottom, top }))
        })
        // Transform from `Option` to `Result`
        .map(|opt| opt.ok_or(AnalysisError::InvalidInput))
        .collect()
}

/// Get a surface based inversion.
pub fn sfc_based_inversion(snd: &Sounding) -> Result<Option<Layer>> {
    let t_profile = snd.temperature_profile();
    let p_profile = snd.pressure_profile();
    let h_profile = snd.height_profile();

    if t_profile.is_empty() || p_profile.is_empty() || h_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    let sfc_layer = izip!(0usize.., p_profile, h_profile, t_profile)
        // Remove levels with missing data
        .filter(|(_, p, h, t)| p.is_some() && h.is_some() && t.is_some())
        // Unpack from `Optioned` type and discard the height, we no longer need it.
        .map(|(i, p, _, t)| (i, p.unpack(), t.unpack()))
        // Only look up to about 700 hPa
        .take_while(|(_, p, _)| *p > HectoPascal(690.0))
        // Pair them up to look at them as a layer
        .tuple_windows::<(_, _)>()
        // Map into inversion types
        .map(|((i0, _, t0), (i1, _, t1))| {
            if t1 > t0 {
                LayerType::Inversion(i0, i1)
            } else {
                LayerType::NonInversion
            }
        })
        // Only take inversion layers to be combined later, if the first one isn't an inversion,
        // then this isn't a surface based inversion!
        .take_while(|lyr| {
            std::mem::discriminant(lyr) == std::mem::discriminant(&LayerType::Inversion(0, 0))
        })
        // Unwrap from the inversion type
        .map(|lyr| match lyr {
            LayerType::Inversion(i0, i1) => (i0, i1),
            LayerType::NonInversion => unreachable!(),
        })
        // Now combine the layers into one
        .fold(None, |combined_layers, (i0, i1)| match combined_layers {
            None => Some((i0, i1)),
            Some((bottom, _top)) => Some((bottom, i1)),
        }) // Option<(usize, usize)
        .and_then(|(bottom_i, top_i)| {
            let bottom = snd.data_row(bottom_i)?;
            let top = snd.data_row(top_i)?;
            Some(Layer { bottom, top })
        });

    Ok(sfc_layer)
}

enum LayerType {
    Inversion(usize, usize),
    NonInversion,
}
