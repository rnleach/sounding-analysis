//! This module finds significant layers such as the dendritic snow growth zone, the hail growth
//! zone, and inversions.

use metfor::{
    Celsius, CelsiusPKm, HectoPascal, JpKg, Km, Meters, MetersPSec, Quantity, WindUV, FREEZING,
};
use optional::Optioned;
use smallvec::SmallVec;
use sounding_base::{DataRow, Sounding};

use crate::error::AnalysisError::*;
use crate::error::*;
use crate::levels::height_level;
use crate::parcel;
use crate::parcel_profile;

/// A layer in the atmosphere described by the values at the top and bottom.
#[derive(Debug, Clone, Copy)]
pub struct Layer {
    /// Pressure at the bottom of the layer.
    pub bottom: DataRow,
    /// Pressure at the top of the layer.
    pub top: DataRow,
}

/// A list of layers.
pub type Layers = SmallVec<[Layer; crate::VEC_SIZE]>;

impl Layer {
    /// Get the average lapse rate in C/km
    pub fn lapse_rate(&self) -> Option<CelsiusPKm> {
        let top_t = self.top.temperature.into_option()?;
        let bottom_t = self.bottom.temperature.into_option()?;

        let dt = (top_t - bottom_t).unpack();
        let dz = Km::from(self.height_thickness()?).unpack();

        Some(CelsiusPKm(dt / dz))
    }

    /// Get the height thickness in meters
    pub fn height_thickness(&self) -> Option<Meters> {
        let top = self.top.height.into_option()?;
        let bottom = self.bottom.height.into_option()?;
        if top == bottom {
            None
        } else {
            Some(top - bottom)
        }
    }

    /// Get the pressure thickness.
    pub fn pressure_thickness(&self) -> Option<HectoPascal> {
        let bottom_p = self.bottom.pressure.into_option()?;
        let top_p = self.top.pressure.into_option()?;
        if bottom_p == top_p {
            None
        } else {
            Some(bottom_p - top_p)
        }
    }

    /// Get the bulk wind shear (spd kts, direction degrees)
    pub fn wind_shear(&self) -> Option<WindUV<MetersPSec>> {
        let top = WindUV::from(self.top.wind.into_option()?);
        let bottom = WindUV::from(self.bottom.wind.into_option()?);

        Some(top - bottom)
    }
}

#[cfg(test)]
mod layer_tests {
    use super::*;
    use metfor::*;
    use optional::some;
    use sounding_base::DataRow;

    fn make_test_layer() -> Layer {
        let mut bottom = DataRow::default();
        bottom.pressure = some(HectoPascal(1000.0));
        bottom.temperature = some(Celsius(20.0));
        bottom.height = some(Meters(5.0));
        bottom.wind = some(WindSpdDir::<Knots> {
            speed: Knots(1.0),
            direction: 180.0,
        });

        let mut top = DataRow::default();
        top.pressure = some(HectoPascal(700.0));
        top.temperature = some(Celsius(-2.0));
        top.height = some(Meters(3012.0));
        top.wind = some(WindSpdDir::<Knots> {
            speed: Knots(1.0),
            direction: 90.0,
        });

        Layer { bottom, top }
    }

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn test_height_thickness() {
        let lyr = make_test_layer();
        println!("{:#?}", lyr);
        assert!(lyr
            .height_thickness()
            .unwrap()
            .approx_eq(Meters(3007.0), Meters(std::f64::EPSILON)));
    }

    #[test]
    fn test_pressure_thickness() {
        let lyr = make_test_layer();
        println!("{:#?}", lyr);
        assert!(lyr
            .pressure_thickness()
            .unwrap()
            .approx_eq(HectoPascal(300.0), HectoPascal(std::f64::EPSILON)));
    }

    #[test]
    fn test_lapse_rate() {
        let lyr = make_test_layer();
        println!(
            "{:#?}\n\n -- \n\n {:#?} \n\n --",
            lyr,
            lyr.lapse_rate().unwrap()
        );
        assert!(lyr
            .lapse_rate()
            .unwrap()
            .approx_eq(CelsiusPKm(-7.31626), CelsiusPKm(1.0e-5)));
    }

    #[test]
    fn test_wind_shear() {
        let lyr = make_test_layer();
        println!(
            "{:#?}\n\n -- \n\n {:#?} \n\n --",
            lyr,
            lyr.wind_shear().unwrap()
        );
        let shear = WindSpdDir::<Knots>::from(lyr.wind_shear().unwrap());
        let speed_shear = shear.abs();
        let WindSpdDir {
            direction: direction_shear,
            ..
        } = shear;

        assert!(speed_shear.approx_eq(Knots(::std::f64::consts::SQRT_2), Knots(1.0e-5)));
        assert!(approx_eq(direction_shear, 45.0, 1.0e-5));
    }
}

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

/// Get a layer that has a certain thickness, like 3km or 6km.
pub fn layer_agl(snd: &Sounding, meters_agl: Meters) -> Result<Layer> {
    let tgt_elev = {
        let elev = snd.station_info().elevation();
        if elev.is_some() {
            elev.unpack() + meters_agl
        } else {
            return Err(MissingValue);
        }
    };

    let bottom = snd.surface_as_data_row().unwrap_or_else(DataRow::default);
    let top = height_level(tgt_elev, snd)?;
    Ok(Layer { bottom, top })
}

/// Get a layer defined by two pressure levels. `bottom_p` > `top_p`
pub fn pressure_layer(snd: &Sounding, bottom_p: HectoPascal, top_p: HectoPascal) -> Result<Layer> {
    let sfc_pressure = snd
        .surface_as_data_row()
        .and_then(|row| row.pressure.into_option());

    if sfc_pressure.is_some() && sfc_pressure.unwrap() < bottom_p {
        return Err(InvalidInput);
    }

    let bottom = crate::interpolation::linear_interpolate_sounding(snd, bottom_p)?;
    let top = crate::interpolation::linear_interpolate_sounding(snd, top_p)?;

    Ok(Layer { bottom, top })
}

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

/// Get the effective inflow layer.
pub fn effective_inflow_layer(snd: &Sounding) -> Option<Layer> {
    let mut vals_iter = snd
        .bottom_up()
        // Convert rows to parcels
        .filter_map(|row| parcel::Parcel::from_datarow(row).map(|pcl| (row, pcl)))
        // Lift the parcel, skip if there is an error
        .filter_map(|(row, pcl)| {
            parcel_profile::lift_parcel(pcl, snd)
                .ok()
                .map(|pcl_anal| (row, pcl_anal))
        })
        // Get the CAPE and CIN
        .filter_map(|(row, pcl_anal)| {
            pcl_anal
                .cape()
                .into_option()
                .and_then(|cape| pcl_anal.cin().map(|cin| (row, cape, cin)))
        })
        // Skip levels until we get one that meets criteria
        .skip_while(|(_, cape, cin)| *cape < JpKg(100.0) || *cin < JpKg(-250.0))
        // Take levels as long as the meet criteria
        .take_while(|(_, cape, cin)| *cape >= JpKg(100.0) && *cin >= JpKg(-250.0))
        // Discard the cape and cin values, we only need the rows
        .map(|(row, _, _)| row);

    let bottom = vals_iter.nth(0);
    let top = vals_iter.last();

    if let (Some(bottom), Some(top)) = (bottom, top) {
        Some(Layer { bottom, top })
    } else {
        None
    }
}
