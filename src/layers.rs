//! This module finds significant layers such as the dendritic snow growth zone, the hail growth
//! zone, and inversions.

use error::*;
use error::AnalysisError::*;
use smallvec::SmallVec;

use sounding_base::{DataRow, Profile, Sounding};
use sounding_base::Profile::*;

const FREEZING: f64 = 0.0;

/// A layer in the atmosphere described by the values at the top and bottom.
#[derive(Debug, Clone, Copy)]
pub struct Layer {
    /// Pressure at the bottom of the layer.
    pub bottom: DataRow,
    /// Pressure at the top of the layer.
    pub top: DataRow,
}

/// A list of layers.
pub type Layers = SmallVec<[Layer; ::VEC_SIZE]>;

impl Layer {
    /// Get the average lapse rate in C/km
    pub fn lapse_rate(&self) -> Result<f64> {
        let top_t = self.top.temperature.ok_or(MissingValue)?;
        let bottom_t = self.bottom.temperature.ok_or(MissingValue)?;

        let dt = top_t - bottom_t;
        let dz = self.height_thickness()?;

        Ok(dt / dz * 1000.0)
    }

    /// Get the height thickness in meters
    #[cfg_attr(feature = "cargo-clippy", allow(float_cmp))]
    pub fn height_thickness(&self) -> Result<f64> {
        let top = self.top.height.ok_or(MissingValue)?;
        let bottom = self.bottom.height.ok_or(MissingValue)?;
        if top == bottom {
            Err(InvalidInput)
        } else {
            Ok(top - bottom)
        }
    }

    /// Get the pressure thickness.
    #[cfg_attr(feature = "cargo-clippy", allow(float_cmp))]
    pub fn pressure_thickness(&self) -> Result<f64> {
        let bottom_p = self.bottom.pressure.ok_or(MissingValue)?;
        let top_p = self.top.pressure.ok_or(MissingValue)?;
        if bottom_p == top_p {
            Err(InvalidInput)
        } else {
            Ok(bottom_p - top_p)
        }
    }

    /// Get the bulk wind shear (spd kts, direction degrees)
    pub fn wind_shear(&self) -> Result<(f64, f64)> {
        let top_spd = self.top.speed.ok_or(MissingValue)?;
        let top_dir = -(90.0 + self.top.direction.ok_or(MissingValue)?);
        let bottom_spd = self.bottom.speed.ok_or(MissingValue)?;
        let bottom_dir = -(90.0 + self.bottom.direction.ok_or(MissingValue)?);

        let top_u = top_dir.to_radians().cos() * top_spd;
        let top_v = top_dir.to_radians().sin() * top_spd;
        let bottom_u = bottom_dir.to_radians().cos() * bottom_spd;
        let bottom_v = bottom_dir.to_radians().sin() * bottom_spd;

        let du = top_u - bottom_u;
        let dv = top_v - bottom_v;

        let shear_spd = du.hypot(dv);
        let mut shear_dir = -(90.0 + dv.atan2(du).to_degrees());

        while shear_dir < 0.0 {
            shear_dir += 360.0;
        }
        while shear_dir > 360.0 {
            shear_dir -= 360.0;
        }

        Ok((shear_spd, shear_dir))
    }
}

#[cfg(test)]
mod layer_tests {
    use super::*;
    use sounding_base::DataRow;

    fn make_test_layer() -> Layer {
        let mut bottom = DataRow::default();
        bottom.pressure = Some(1000.0);
        bottom.temperature = Some(20.0);
        bottom.height = Some(5.0);
        bottom.speed = Some(1.0);
        bottom.direction = Some(180.0);

        let mut top = DataRow::default();
        top.pressure = Some(700.0);
        top.temperature = Some(-2.0);
        top.height = Some(3012.0);
        top.speed = Some(1.0);
        top.direction = Some(90.0);

        Layer { bottom, top }
    }

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn test_height_thickness() {
        let lyr = make_test_layer();
        println!("{:#?}", lyr);
        assert!(approx_eq(
            lyr.height_thickness().unwrap(),
            3007.0,
            ::std::f64::EPSILON
        ));
    }

    #[test]
    fn test_pressure_thickness() {
        let lyr = make_test_layer();
        println!("{:#?}", lyr);
        assert!(approx_eq(
            lyr.pressure_thickness().unwrap(),
            300.0,
            ::std::f64::EPSILON
        ));
    }

    #[test]
    fn test_lapse_rate() {
        let lyr = make_test_layer();
        println!(
            "{:#?}\n\n -- \n\n {} \n\n --",
            lyr,
            lyr.lapse_rate().unwrap()
        );
        assert!(approx_eq(lyr.lapse_rate().unwrap(), -7.31626, 1.0e-5));
    }

    #[test]
    fn test_wind_shear() {
        let lyr = make_test_layer();
        println!(
            "{:#?}\n\n -- \n\n {:#?} \n\n --",
            lyr,
            lyr.wind_shear().unwrap()
        );
        let (speed_shear, direction_shear) = lyr.wind_shear().unwrap();
        assert!(approx_eq(speed_shear, ::std::f64::consts::SQRT_2, 1.0e-5));
        assert!(approx_eq(direction_shear, 45.0, 1.0e-5));
    }
}

/// Find the dendtritic growth zones throughout the profile. It is unusual, but possible there is
/// more than one.
///
/// If there are none, then an empty vector is returned.
pub fn dendritic_snow_zone(snd: &Sounding) -> Result<Layers> {
    let mut to_return: SmallVec<[Layer; ::VEC_SIZE]> = SmallVec::new();

    // Dendritic snow growth zone temperature range in C
    const WARM_SIDE: f64 = -12.0;
    const COLD_SIDE: f64 = -18.0;
    const TOP_PRESSURE: f64 = 300.0; // don't look above here.

    let t_profile = snd.get_profile(Temperature);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    let mut profile = t_profile.iter().zip(p_profile);

    let mut bottom_press = ::std::f64::MAX;
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
        } else {
            return Err(AnalysisError::NoDataProfile);
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

            if press < TOP_PRESSURE {
                break;
            }

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
                push_layer(bottom_press, top_press, snd, &mut to_return)?;
            }
            // Crossed out of zone to cold side
            if last_t >= COLD_SIDE && t < COLD_SIDE {
                top_press = ::interpolation::linear_interp(COLD_SIDE, last_t, t, last_press, press);
                push_layer(bottom_press, top_press, snd, &mut to_return)?;
            }
            last_t = t;
            last_press = press;
        }
    }

    // Check to see if we ended in a dendtritic zone
    if last_t <= WARM_SIDE && last_t >= COLD_SIDE {
        top_press = last_press;

        push_layer(bottom_press, top_press, snd, &mut to_return)?;
    }

    Ok(to_return)
}

/// Assuming it is below freezing at the surface, this will find the warm layers aloft using the
/// dry bulb temperature. Does not look above 500 hPa.
pub fn warm_temperature_layer_aloft(snd: &Sounding) -> Result<Layers> {
    warm_layer_aloft(snd, Temperature)
}

/// Assuming the wet bulb temperature is below freezing at the surface, this will find the warm
/// layers aloft using the wet bulb temperature. Does not look above 500 hPa.
pub fn warm_wet_bulb_layer_aloft(snd: &Sounding) -> Result<Layers> {
    warm_layer_aloft(snd, WetBulb)
}

fn warm_layer_aloft(snd: &Sounding, var: Profile) -> Result<Layers> {
    use std::f64::MAX;

    assert!(var == Temperature || var == WetBulb);

    let mut to_return: Layers = Layers::new();

    let t_profile = snd.get_profile(var);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    p_profile
        .iter()
        .zip(t_profile)
        // Remove levels without pressure AND temperature data
        .filter_map(|pair|{
            if let (&Some(p), &Some(t)) = pair {
                Some((p,t))
            } else {
                None
            }
        })
        // Ignore anything above 500 hPa, extremely unlikely for a warm layer up there.
        .take_while(|&(p,_)| p > 500.0 )
        // Find the warm layers!
        .fold(Ok((MAX, MAX, None)), |last_iter_res: Result<(f64,f64, Option<_>)>, (p,t)|{
            let (last_p, last_t, mut bottom) = last_iter_res?;
            if last_t <= FREEZING && t > FREEZING && bottom.is_none() {
                // Entering a warm layer.
                let bottom_p = ::interpolation::linear_interp(FREEZING, last_t, t, last_p, p);
                bottom = Some(::interpolation::linear_interpolate(snd, bottom_p)?);
            }
            if bottom.is_some() && last_t > FREEZING && t <= FREEZING {
                // Crossed out of a warm layer
                let top_p = ::interpolation::linear_interp(FREEZING, last_t, t, last_p, p);
                let top = ::interpolation::linear_interpolate(snd, top_p)?;
                {
                let bottom = bottom.unwrap();
                to_return.push(Layer{bottom, top});}
                bottom = None;
            }

            Ok((p,t,bottom))
        })?;

    Ok(to_return)
}

/// Assuming a warm layer aloft given by `warm_layers`, measure the cold surface layer.
pub fn cold_surface_temperature_layer(snd: &Sounding, warm_layers: &[Layer]) -> Result<Layer> {
    cold_surface_layer(snd, Temperature, warm_layers)
}

fn cold_surface_layer(snd: &Sounding, var: Profile, warm_layers: &[Layer]) -> Result<Layer> {
    debug_assert!(var == Temperature || var == WetBulb);

    if warm_layers.is_empty() {
        return Err(InvalidInput);
    }

    let t_profile = snd.get_profile(var);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(MissingProfile); // Should not happen since we already used these to get warm layer
    }

    let mut profile = t_profile.iter().zip(p_profile);

    let last_t: f64;
    let last_press: f64;
    loop {
        if let Some((t, press)) = profile.next() {
            if let (Some(t), Some(press)) = (*t, *press) {
                last_t = t;
                last_press = press;
                break;
            }
        } else {
            return Err(NoDataProfile);
        }
    }

    // Check to see if we are below freezing at the bottom
    if last_t > FREEZING {
        return Err(InvalidInput);
    }

    let bottom = ::interpolation::linear_interpolate(snd, last_press)?;

    Ok(Layer {
        bottom,
        top: warm_layers[0].bottom,
    })
}

/// Get a layer that has a certain thickness, like 3km or 6km.
pub fn layer_agl(snd: &Sounding, meters_agl: f64) -> Result<Layer> {
    let tgt_elev = if let Some(elev) = snd.get_station_info().elevation() {
        elev + meters_agl
    } else {
        return Err(MissingValue);
    };

    let h_profile = snd.get_profile(GeopotentialHeight);
    let p_profile = snd.get_profile(Pressure);

    if h_profile.is_empty() || p_profile.is_empty() {
        return Err(MissingProfile);
    }

    let mut profile = h_profile.iter().zip(p_profile).filter_map(|pair| {
        if let (&Some(h), &Some(p)) = pair {
            Some((h, p))
        } else {
            None
        }
    });

    let bottom = snd.surface_as_data_row();

    // Initialize the bottom of the sounding
    let mut last_h: f64;
    let mut last_press: f64;
    // Try surface data
    if let (Some(h), Some(p)) = (bottom.height, bottom.pressure) {
        last_h = h;
        last_press = p;
    } else if let Some((h, press)) = profile.by_ref().next() {
        // Find lowest level in sounding
        last_h = h;
        last_press = press;
    } else {
        return Err(AnalysisError::NoDataProfile);
    }

    if last_h > tgt_elev {
        return Err(AnalysisError::NotEnoughData);
    }

    for (h, press) in profile {
        if last_h <= tgt_elev && h > tgt_elev {
            let top_press = ::interpolation::linear_interp(tgt_elev, last_h, h, last_press, press);

            let top = ::interpolation::linear_interpolate(snd, top_press)?;

            return Ok(Layer { bottom, top });
        }
        last_h = h;
        last_press = press;
    }

    Err(AnalysisError::NotEnoughData)
}

/// Get a layer defined by two pressure levels. `bottom_p` > `top_p`
pub fn pressure_layer(snd: &Sounding, bottom_p: f64, top_p: f64) -> Result<Layer> {
    let sfc_pressure = snd.surface_as_data_row().pressure;

    if sfc_pressure.is_some() && sfc_pressure.unwrap() < bottom_p {
        return Err(InvalidInput);
    }

    let bottom = ::interpolation::linear_interpolate(snd, bottom_p)?;
    let top = ::interpolation::linear_interpolate(snd, top_p)?;

    Ok(Layer { bottom, top })
}

/// Get all inversion layers up to 500 mb.
pub fn inversions(snd: &Sounding) -> Result<Layers> {
    let mut to_return: Layers = Layers::new();

    let t_profile = snd.get_profile(Temperature);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() || p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    p_profile
        .iter()
        .zip(t_profile)
        .enumerate()
        // Filter out rows without both temperature and pressure.
        .filter_map(|triple| {
            if let (i, (&Some(_), &Some(t))) = triple {
                Some((i,t))
            } else {
                None
            }
        })
        // Capture the inversion layers
        .fold((0,::std::f64::MAX, None),
            |(last_i, last_t, mut lyr), (i,t)|{

                if lyr.is_none() && last_t < t {
                    // Coming into an inversion
                    if let Some(dr) = snd.get_data_row(last_i){
                        lyr = Some(Layer{bottom: dr, top: DataRow::default()});
                    } else {
                        unreachable!();
                    }
                } else if lyr.is_some() && last_t > t {
                    // Leaving an inversion
                    if let Some(dr) = snd.get_data_row(last_i){
                        if let Some(Layer{bottom: btm, ..}) = lyr {
                            to_return.push(Layer{bottom: btm, top: dr});
                            lyr = None;
                        } else {
                            unreachable!();
                        }
                    } else {
                        unreachable!();
                    }
                }
                (i, t, lyr)
            });

    Ok(to_return)
}

fn push_layer(
    bottom_press: f64,
    top_press: f64,
    snd: &Sounding,
    target_vec: &mut SmallVec<[Layer; ::VEC_SIZE]>,
) -> Result<()> {
    let bottom = ::interpolation::linear_interpolate(snd, bottom_press)?;
    let top = ::interpolation::linear_interpolate(snd, top_press)?;
    target_vec.push(Layer { bottom, top });
    Ok(())
}
