use metfor::{Celsius, HectoPascal, Kelvin, Knots, Meters, PaPS, WindSpdDir};
use optional::Optioned;

/// A copy of a row of the sounding data.
#[derive(Clone, Default, Copy, Debug, PartialEq)]
pub struct DataRow {
    /// Pressure in hPa
    pub pressure: Optioned<HectoPascal>,
    /// Temperature in C
    pub temperature: Optioned<Celsius>,
    /// Wet bulb temperature in C
    pub wet_bulb: Optioned<Celsius>,
    /// Dew point in C
    pub dew_point: Optioned<Celsius>,
    /// Equivalent potential temperature in Kelvin
    pub theta_e: Optioned<Kelvin>,
    /// Wind
    pub wind: Optioned<WindSpdDir<Knots>>,
    /// Pressure vertical velocity in Pa/sec
    pub pvv: Optioned<PaPS>,
    /// Geopotential Height in meters
    pub height: Optioned<Meters>,
    /// Cloud fraction in percent
    pub cloud_fraction: Optioned<f64>,
}
