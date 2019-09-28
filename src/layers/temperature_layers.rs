use super::{Layer, Layers};
use crate::{
    error::{
        AnalysisError::{InvalidInput, MissingProfile},
        {AnalysisError, Result},
    },
    interpolation::{linear_interp, linear_interpolate_sounding},
    sounding::Sounding,
};
use itertools::{izip, Itertools};
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

#[inline]
fn temperature_layer(
    snd: &Sounding,
    warm_side: Celsius,
    cold_side: Celsius,
    top_pressure: HectoPascal,
) -> Result<Layers> {
    let t_profile = snd.temperature_profile();
    let p_profile = snd.pressure_profile();

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    // Assuming were are iterating over the profile from the bottom (ground level) up, this type
    // records levels where the sounding crosses the targeted temperature layer.
    enum Crossing {
        // The level where the sounding crosses into a layer.
        Enter(HectoPascal),
        // The level where the sounding crosses out of a layer.
        Exit(HectoPascal),
        // The first element must be greater than the second! These levels completely jumped over
        // a temperature layer.
        Over(HectoPascal, HectoPascal),
        // This level is completely inside a temperature layer.
        In(HectoPascal),
    }

    // Take two levels and detect if traversing from the 0th to the 1st level crosses the
    // target temperature layer.
    fn to_crossing_type(
        pnt0: (HectoPascal, Celsius),
        pnt1: (HectoPascal, Celsius),
        warm_side: Celsius,
        cold_side: Celsius,
    ) -> Option<Crossing> {
        let (p0, t0) = pnt0;
        let (p1, t1) = pnt1;

        debug_assert!(p0 >= p1, "Pressure must decrease with height.");

        if t0 < cold_side && t1 >= cold_side && t1 <= warm_side {
            // Crossing into a layer from the cold side
            let cold_p = linear_interp(cold_side, t0, t1, p0, p1);
            Some(Crossing::Enter(cold_p))
        } else if t0 > warm_side && t1 <= warm_side && t1 >= cold_side {
            // Crossing into layer from the warm side
            let warm_p = linear_interp(warm_side, t0, t1, p0, p1);
            Some(Crossing::Enter(warm_p))
        } else if (t0 < cold_side && t1 > warm_side) || (t0 > warm_side && t1 < cold_side) {
            // Crossed over a layer
            let warm_p = linear_interp(warm_side, t0, t1, p0, p1);
            let cold_p = linear_interp(cold_side, t0, t1, p0, p1);
            Some(Crossing::Over(warm_p.max(cold_p), warm_p.min(cold_p)))
        } else if t0 > cold_side && t0 < warm_side && t1 > warm_side {
            // Crossed out of a layer into the warm side
            let warm_p = linear_interp(warm_side, t0, t1, p0, p1);
            Some(Crossing::Exit(warm_p))
        } else if t0 > cold_side && t0 < warm_side && t1 < cold_side {
            // Crossed out of a layer into the cold side
            let cold_p = linear_interp(cold_side, t0, t1, p0, p1);
            Some(Crossing::Exit(cold_p))
        } else if t0 >= cold_side && t0 <= warm_side && t1 >= cold_side && t1 <= warm_side {
            // We're in the midst of a layer
            Some(Crossing::In(p0))
        } else {
            None
        }
    }

    izip!(p_profile, t_profile)
        // Remove levels with missing values
        .filter(|(p, t)| p.is_some() && t.is_some())
        // Unwrap from the `Optioned` type
        .map(|(p, t)| (p.unpack(), t.unpack()))
        // Only take values up to a certain level
        .take_while(move |&(p, _)| p > top_pressure)
        // Make adjacent points into pairs so we can look at them two at a time
        .tuple_windows::<(_, _)>()
        // Map to a crossing type and filter out those that we don't need
        .filter_map(|(pnt0, pnt1)| to_crossing_type(pnt0, pnt1, warm_side, cold_side))
        // Scan the iterator and coalesce crossings into levels
        .scan(None, |bottom_p: &mut Option<_>, crossing_type: Crossing| {
            match crossing_type {
                Crossing::In(p) => {
                    if bottom_p.is_none() {
                        // to get here we started out in the layer and never had to CROSS into it
                        *bottom_p = Some(p);
                    }
                    // Yield this to indicate we're not done iterating, but we don't yet have a
                    // pair of values (top and bottom) to create a layer from.
                    Some(None)
                }
                Crossing::Enter(p) => {
                    *bottom_p = Some(p);
                    Some(None)
                }
                Crossing::Exit(top_p) => {
                    debug_assert!(bottom_p.is_some());
                    Some(Some((bottom_p.take().unwrap(), top_p)))
                }
                Crossing::Over(p0, p1) => {
                    debug_assert!(bottom_p.is_none());
                    Some(Some((p0, p1)))
                }
            }
        })
        // Filter out and unwrap steps that didn't yield a complete layer
        .filter_map(|opt| opt)
        // Interpolate the sounding and create a layer.
        .map(|(bottom_p, top_p)| {
            linear_interpolate_sounding(snd, bottom_p).and_then(|bottom| {
                linear_interpolate_sounding(snd, top_p).map(|top| Layer { bottom, top })
            })
        }) // Result<Layer>
        .collect()
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

#[inline]
fn warm_layer_aloft(snd: &Sounding, t_profile: &[Optioned<Celsius>]) -> Result<Layers> {
    let p_profile = snd.pressure_profile();

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    enum Crossing {
        IntoWarmLayer(HectoPascal),
        OutOfWarmLayer(HectoPascal),
    }

    fn crossing_type(
        pnt0: (HectoPascal, Celsius),
        pnt1: (HectoPascal, Celsius),
    ) -> Option<Crossing> {
        let (p0, t0) = pnt0;
        let (p1, t1) = pnt1;

        if t0 < FREEZING && t1 >= FREEZING {
            let crossing_p = linear_interp(FREEZING, t0, t1, p0, p1);
            Some(Crossing::IntoWarmLayer(crossing_p))
        } else if t0 >= FREEZING && t1 < FREEZING {
            let crossing_p = linear_interp(FREEZING, t0, t1, p0, p1);
            Some(Crossing::OutOfWarmLayer(crossing_p))
        } else {
            None
        }
    }

    izip!(p_profile, t_profile)
        // Remove levels with any missing temperature or pressure data
        .filter(|(p, t)| p.is_some() && t.is_some())
        // Unpack from the `Optioned` type
        .map(|(p, t)| (p.unpack(), t.unpack()))
        // Ignore anything above 500 hPa, extremely unlikely for a warm layer up there.
        .take_while(|&(p, _)| p > HectoPascal(500.0))
        // Skip any levels that start out in a surface based warm layer
        .skip_while(|&(_, t)| t > FREEZING)
        // Pair them up to look at two adjacent points at a time
        .tuple_windows::<(_, _)>()
        // Map into the crossing type, and filter out any levels that aren't a crossing
        .filter_map(|(pnt0, pnt1)| crossing_type(pnt0, pnt1))
        // Scan the crossings to create layers
        .scan(
            None,
            |bottom: &mut Option<_>, crossing: Crossing| match crossing {
                Crossing::IntoWarmLayer(bottom_p) => {
                    debug_assert!(bottom.is_none());
                    *bottom = Some(bottom_p);
                    Some(None)
                }
                Crossing::OutOfWarmLayer(top_p) => {
                    // If bottom is None, this is a surface based warm layer, not a warm layer
                    // aloft.
                    Some(bottom.take().map(|bottom_p| (bottom_p, top_p)))
                }
            },
        )
        // Filter out and unwrap steps that didn't yield a complete layer
        .filter_map(|opt| opt)
        // Interpolate the sounding and create a layer.
        .map(|(bottom_p, top_p)| {
            linear_interpolate_sounding(snd, bottom_p).and_then(|bottom| {
                linear_interpolate_sounding(snd, top_p).map(|top| Layer { bottom, top })
            })
        }) // Result<Layer>
        .collect()
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
        // Remove layers with missing data
        .filter(|(_, p, t)| p.is_some() && t.is_some())
        // Unpack from the `Optioned` type and discard the pressure, we don't need it anymore
        .map(|(i, _, t)| (i, t.unpack()))
        // For a surface based cold layer, we must start out below zero,
        .take_while(|(_, t)| *t <= FREEZING)
        // We only want the first one, the bottom of the layer, if it exists at all
        .nth(0) // Option<(i, t)>
        // translate it to a full data row in the sounding
        .and_then(|(i, _)| snd.data_row(i)) // Option<DataRow>
        // Add the top level, if it exists
        .and_then(|bottom| warm_layers.get(0).map(|lyr| (bottom, lyr.bottom)))
        // Package it up in a layer
        .map(|(bottom, top)| Layer { bottom, top }) // Option<Layer>
        // Transform from option into result
        .ok_or(InvalidInput)
}
