//! This module finds significant layers such as the dendritic snow growth zone, the hail growth
//! zone, and inversions.
use crate::sounding::DataRow;
use metfor::{CelsiusPKm, HectoPascal, Km, Meters, MetersPSec, Quantity, WindUV};

/// A layer in the atmosphere described by the values at the top and bottom.
#[derive(Debug, Clone, Copy)]
pub struct Layer {
    /// Sounding values at the bottom of the layer.
    pub bottom: DataRow,
    /// Sounding values at the top of the layer.
    pub top: DataRow,
}

/// A list of layers.
pub type Layers = Vec<Layer>;

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
    use crate::sounding::DataRow;
    use metfor::*;
    use optional::some;

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

mod temperature_layers;
pub use temperature_layers::{
    cold_surface_temperature_layer, dendritic_snow_zone, hail_growth_zone,
    warm_temperature_layer_aloft, warm_wet_bulb_layer_aloft,
};

mod height_pressure;
pub use height_pressure::{layer_agl, pressure_layer};

mod inversions;
pub use inversions::{inversions, sfc_based_inversion};

mod convective;
pub use convective::effective_inflow_layer;
