//! This module finds significant layers such as the dendritic snow growth zone, the hail growth
//! zone, and inversions.

use error::*;
use smallvec::SmallVec;

use sounding_base::{DataRow, Profile, Sounding};
use sounding_base::Profile::*;

/// A layer in the atmosphere described by the values at the top and bottom.
#[derive(Debug, Clone, Copy)]
pub struct Layer {
    /// Pressure at the bottom of the layer.
    pub bottom: DataRow,
    /// Pressure at the top of the layer.
    pub top: DataRow,
}

impl Layer {
    /// Get the average lapse rate in C/km
    pub fn lapse_rate(&self) -> Option<f64> {
        let dt = self.top.temperature? - self.bottom.temperature?;
        let dz = self.height_thickness()?;
        Some(dt / dz * 1000.0)
    }

    /// Get the height thickness in meters
    pub fn height_thickness(&self) -> Option<f64> {
        Some(self.top.height? - self.bottom.height?)
    }

    /// Get the pressure thickness.
    pub fn pressure_thickness(&self) -> Option<f64> {
        Some(self.bottom.pressure? - self.top.pressure?)
    }

    /// Get the bulk wind shear (spd kts, direction degrees)
    pub fn wind_shear(&self) -> Option<(f64, f64)> {
        let top_spd = self.top.speed?;
        let top_dir = self.top.direction?;
        let bottom_spd = self.bottom.speed?;
        let bottom_dir = self.bottom.direction?;

        let top_u = top_dir.to_radians().cos() * top_spd;
        let top_v = top_dir.to_radians().sin() * top_spd;
        let bottom_u = bottom_dir.to_radians().cos() * bottom_spd;
        let bottom_v = bottom_dir.to_radians().sin() * bottom_spd;

        let du = top_u - bottom_u;
        let dv = top_v - bottom_v;

        let shear_spd = du.hypot(dv);
        let mut shear_dir = dv.atan2(du).to_degrees();

        while shear_dir < 0.0 {
            shear_dir += 360.0;
        }
        while shear_dir > 360.0 {
            shear_dir -= 360.0;
        }

        Some((shear_spd, shear_dir))
    }
}

/// Find the dendtritic growth zones throughout the profile. It is unusual, but possible there is
/// more than one.
///
/// # Errors
/// If the sounding is missing a temperature or pressure profile, `error::ErrorKind::MissingProfile`
/// is returned in the result. Otherwise, if no dendritic layers are found, an empty vector is
/// returned in the `Result`
pub fn dendritic_snow_zone(snd: &Sounding) -> Result<SmallVec<[Layer; ::VEC_SIZE]>> {
    let mut to_return: SmallVec<[Layer; ::VEC_SIZE]> = SmallVec::new();

    // Dendritic snow growth zone temperature range in C
    const WARM_SIDE: f64 = -12.0;
    const COLD_SIDE: f64 = -18.0;
    const TOP_PRESSURE: f64 = 300.0; // don't look above here.

    let t_profile = snd.get_profile(Temperature);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }
    if p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    let mut profile = t_profile.iter().zip(p_profile);

    let mut bottom_press = ::std::f64::MAX; // Only init because compiler can't tell value not used
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

    fn push_layer(
        bottom_press: f64,
        top_press: f64,
        snd: &Sounding,
        target_vec: &mut SmallVec<[Layer; ::VEC_SIZE]>,
    ) {
        let bottom = ::interpolation::linear_interpolate(snd, bottom_press);
        let top = ::interpolation::linear_interpolate(snd, top_press);
        target_vec.push(Layer { bottom, top });
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
                push_layer(bottom_press, top_press, snd, &mut to_return);
            }
            // Crossed out of zone to cold side
            if last_t >= COLD_SIDE && t < COLD_SIDE {
                top_press = ::interpolation::linear_interp(COLD_SIDE, last_t, t, last_press, press);
                push_layer(bottom_press, top_press, snd, &mut to_return);
            }
            last_t = t;
            last_press = press;
        }
    }

    // Check to see if we ended in a dendtritic zone
    if last_t <= WARM_SIDE && last_t >= COLD_SIDE {
        top_press = last_press;

        push_layer(bottom_press, top_press, snd, &mut to_return);
    }

    Ok(to_return)
}

/// Assuming it is below freezing at the surface, this will find the warm layers aloft using the
/// dry bulb temperature. Does not look above 500 hPa.
pub fn warm_temperature_layer_aloft(snd: &Sounding) -> Result<SmallVec<[Layer; ::VEC_SIZE]>> {
    warm_layer_aloft(snd, Temperature)
}

/// Assuming the wet bulb temperature is below freezing at the surface, this will find the warm
/// layers aloft using the wet bulb temperature. Does not look above 500 hPa.
pub fn warm_wet_bulb_layer_aloft(snd: &Sounding) -> Result<SmallVec<[Layer; ::VEC_SIZE]>> {
    warm_layer_aloft(snd, WetBulb)
}

fn warm_layer_aloft(snd: &Sounding, var: Profile) -> Result<SmallVec<[Layer; ::VEC_SIZE]>> {
    assert!(var == Temperature || var == WetBulb);

    let mut to_return: SmallVec<[Layer; ::VEC_SIZE]> = SmallVec::new();

    const FREEZING: f64 = 0.0;

    let t_profile = snd.get_profile(var);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }
    if p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    let mut profile = t_profile.iter().zip(p_profile);

    let mut bottom_press = ::std::f64::MAX; // Only init because compiler can't tell value not used
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
            match var {
                Temperature => return Err(AnalysisError::NoDataProfile),
                WetBulb => return Err(AnalysisError::NoDataProfile),
                _ => unreachable!(),
            }
        }
    }

    // Check to see if we are below freezing at the bottom
    if last_t > FREEZING {
        return Ok(to_return);
    }

    let mut in_warm_zone = false;

    for (t, press) in profile {
        if let (Some(t), Some(press)) = (*t, *press) {
            if press < 500.0 {
                break;
            }

            if last_t <= FREEZING && t > FREEZING {
                bottom_press =
                    ::interpolation::linear_interp(FREEZING, last_t, t, last_press, press);
                in_warm_zone = true;
            }
            // Crossed out of zone to warm side
            if last_t > FREEZING && t <= FREEZING {
                top_press = ::interpolation::linear_interp(FREEZING, last_t, t, last_press, press);

                let bottom = ::interpolation::linear_interpolate(snd, bottom_press);
                let top = ::interpolation::linear_interpolate(snd, top_press);
                to_return.push(Layer { bottom, top });
                in_warm_zone = false;
            }
            last_t = t;
            last_press = press;
        }
    }

    // Check to see if we ended in a warm layer aloft
    if last_t > FREEZING && in_warm_zone {
        top_press = last_press;
        let bottom = ::interpolation::linear_interpolate(snd, bottom_press);
        let top = ::interpolation::linear_interpolate(snd, top_press);
        to_return.push(Layer { bottom, top });
    }

    Ok(to_return)
}

/// Assuming a warm layer aloft given by warm_layers, measure the cold surface layer.
pub fn cold_surface_temperature_layer(snd: &Sounding, warm_layers: &[Layer]) -> Option<Layer> {
    cold_surface_layer(snd, Temperature, warm_layers)
}

fn cold_surface_layer(snd: &Sounding, var: Profile, warm_layers: &[Layer]) -> Option<Layer> {
    assert!(var == Temperature || var == WetBulb);

    const FREEZING: f64 = 0.0;

    if warm_layers.is_empty() {
        return None;
    }

    let t_profile = snd.get_profile(var);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() || p_profile.is_empty() {
        return None; // Should not happen since we already used these to get warm layer
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
            return None;
        }
    }

    // Check to see if we are below freezing at the bottom
    if last_t > FREEZING {
        return None;
    }

    let bottom = ::interpolation::linear_interpolate(snd, last_press);

    Some(Layer {
        bottom,
        top: warm_layers[0].bottom,
    })
}

/// Get a layer that has a certain thickness, like 3km or 6km.
pub fn layer_agl(snd: &Sounding, meters_agl: f64) -> Result<Layer> {
    let tgt_elev = if let Some(elev) = snd.get_station_info().elevation() {
        elev + meters_agl
    } else {
        return Err(AnalysisError::MissingValue);
    };

    let h_profile = snd.get_profile(GeopotentialHeight);
    let p_profile = snd.get_profile(Pressure);

    if h_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }
    if p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
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
    } else {
        // Find lowest level in sounding
        if let Some((h, press)) = profile.by_ref().next() {
            last_h = h;
            last_press = press;
        } else {
            return Err(AnalysisError::NoDataProfile);
        }
    }

    if last_h > tgt_elev {
        return Err(AnalysisError::NotEnoughData);
    }

    for (h, press) in profile {
        if last_h <= tgt_elev && h > tgt_elev {
            let top_press = ::interpolation::linear_interp(tgt_elev, last_h, h, last_press, press);

            let top = ::interpolation::linear_interpolate(snd, top_press);

            return Ok(Layer { bottom, top });
        }
        last_h = h;
        last_press = press;
    }

    Err(AnalysisError::NotEnoughData)
}

/// Get a layer defined by two pressure levels. `bottom_p` > `top_p`
pub fn pressure_layer(snd: &Sounding, bottom_p: f64, top_p: f64) -> Option<Layer> {
    let sfc = snd.surface_as_data_row();

    if sfc.pressure.is_some() && sfc.pressure.unwrap() < bottom_p {
        return None;
    }

    let bottom = ::interpolation::linear_interpolate(snd, bottom_p);
    let top = ::interpolation::linear_interpolate(snd, top_p);

    Some(Layer { bottom, top })
}

/// Get all inversion layers up 500 mb.
pub fn inversions(snd: &Sounding) -> Result<SmallVec<[Layer; ::VEC_SIZE]>> {
    let mut to_return: SmallVec<[Layer; ::VEC_SIZE]> = SmallVec::new();

    let t_profile = snd.get_profile(Temperature);
    let p_profile = snd.get_profile(Pressure);

    if t_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }
    if p_profile.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    let profile = t_profile
        .iter()
        .zip(p_profile)
        .enumerate()
        .filter_map(|triplet| {
            if let (i, (&Some(t), &Some(_))) = triplet {
                Some((i, t))
            } else {
                None
            }
        });

    let mut window = if let Some(window) = Window::new_with_iterator((0usize, 0.0f64), profile) {
        window
    } else {
        return Err(AnalysisError::NotEnoughData);
    };

    let mut bottom_idx = 0usize; // Only init because compiler can't tell value not used
    let mut top_idx: usize;
    let mut in_inversion = false;

    while window.slide() {
        let data = window.view();
        if !in_inversion {
            let mut all_increasing = true;
            let mut last_t = -::std::f64::MAX;
            for &(_, t) in data {
                all_increasing = all_increasing && t > last_t;
                last_t = t;
            }

            if all_increasing {
                bottom_idx = data[0].0;
                in_inversion = true;
            }
        } else {
            let mut all_decreasing = true;
            let mut last_t = ::std::f64::MAX;
            for &(_, t) in data {
                all_decreasing = all_decreasing && t < last_t;
                last_t = t;
            }

            if all_decreasing {
                top_idx = data[0].0;
                in_inversion = false;

                if let Some(bottom) = snd.get_data_row(bottom_idx) {
                    if let Some(top) = snd.get_data_row(top_idx) {
                        to_return.push(Layer { bottom, top });
                    }
                }
            }
        }
    }

    Ok(to_return)
}

const WINDOW_SIZE: usize = 3;
struct Window<T, I> {
    window: [T; WINDOW_SIZE],
    iter: I,
}

impl<T, I> Window<T, I>
where
    T: Copy,
    I: Iterator<Item = T>,
{
    fn new_with_iterator(seed: T, mut iter: I) -> Option<Self> {
        let mut window = [seed; WINDOW_SIZE];
        let mut count = 0;

        while let Some(val) = iter.by_ref().next() {
            window[count] = val;
            count += 1;
            if count == WINDOW_SIZE - 1 {
                break;
            }
        }

        if count == WINDOW_SIZE - 1 {
            Some(Window { window, iter })
        } else {
            None
        }
    }

    fn slide(&mut self) -> bool
    where
        I: Iterator<Item = T>,
    {
        for i in 0..(WINDOW_SIZE - 1) {
            self.window[i] = self.window[i + 1];
        }

        if let Some(val) = self.iter.next() {
            self.window[WINDOW_SIZE - 1] = val;
            true
        } else {
            false
        }
    }

    fn view(&self) -> &[T] {
        &self.window
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use test_data;

    #[test]
    fn simple_dendritic_layer() {
        let (snd, tgt_float_vals, tgt_int_vals) = test_data::load_test_file("standard.csv");

        // Check for number of dendritic growth zones.
        if let Some(num_dendritic_layers) = tgt_int_vals.get("num dendritic zones") {
            let num_dendritic_layers = *num_dendritic_layers as usize;
            let analyzed_num = dendritic_snow_zone(&snd).unwrap().len();
            assert!(num_dendritic_layers == analyzed_num);
        } else {
            panic!("No dendritic zones in test file.")
        }

        // Check the pressure levels of those growth zones.
        if let Some(dendritic_zone_pressures) = tgt_float_vals.get("dendritic zone pressures") {
            let dendritic_zone_pressures = dendritic_zone_pressures.chunks(2);
            let analyzed_layers = dendritic_snow_zone(&snd).unwrap();
            let analyzed_layers = analyzed_layers.iter();
            for (lyr, it) in analyzed_layers.zip(dendritic_zone_pressures) {
                println!(
                    "\nbottom {:?}  ---  {:?}",
                    lyr.bottom.pressure.unwrap(),
                    it[0]
                );
                assert!(test_data::approx_equal(
                    lyr.bottom.pressure.unwrap(),
                    it[0],
                    0.1
                ));
                println!("top {:?}  ---  {:?}", lyr.top.pressure.unwrap(), it[1]);
                assert!(test_data::approx_equal(
                    lyr.top.pressure.unwrap(),
                    it[1],
                    0.1
                ));
            }
        }
    }

    // TODO: Test complex dendritic layer.

    #[test]
    fn test_warm_layer_aloft() {
        let (snd, _tgt_float_vals, tgt_int_vals) = test_data::load_test_file("standard.csv");

        // Check for number of warm layers aloft.
        if let Some(num_warm_layers) = tgt_int_vals.get("num warm dry bulb aloft") {
            let num_warm_layers = *num_warm_layers as usize;
            let analyzed_num = warm_temperature_layer_aloft(&snd).unwrap().len();
            assert!(num_warm_layers == analyzed_num);
        } else {
            panic!("No warm dry bulb layer info in test file.")
        }

        // Check for number of warm layers aloft.
        if let Some(num_warm_layers) = tgt_int_vals.get("num warm wet bulb aloft") {
            let num_warm_layers = *num_warm_layers as usize;
            let analyzed_num = warm_wet_bulb_layer_aloft(&snd).unwrap().len();
            assert!(num_warm_layers == analyzed_num);
        } else {
            panic!("No warm wet bulb layer info in test file.")
        }
    }

    // TODO: test with actual frozen surface and warm layer.
}
