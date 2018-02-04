//! Formulas for meteorological calculations.

/// Potential temperature in Kelvin
pub fn theta_kelvin(pressure_hpa: f64, temperature_c: f64) -> f64 {
    use std::f64;

    (temperature_c + 273.15) * f64::powf(1000.0 / pressure_hpa, 0.286)
}

/// Temperture in C
pub fn temperature_c_from_theta(theta_kelvin: f64, pressure_hpa: f64) -> f64 {
    use std::f64;

    theta_kelvin * f64::powf(pressure_hpa / 1000.0, 0.286) - 273.15
}

/// Get the vapor pressure of water as a function of temperature in hPa
pub fn vapor_pressure_water(temperature_c: f64) -> f64 {
    use std::f64;

    6.11 * f64::powf(10.0, 7.5 * temperature_c / (237.3 + temperature_c))
}

/// Given the temperature and dew point in celsius, calculate the relative humidity.
pub fn rh(temperature_c: f64, dew_point_c: f64) -> f64 {
    let e = vapor_pressure_water(dew_point_c);
    let es = vapor_pressure_water(temperature_c);
    e / es
}

/// Get the mixing ratio in g/kg.
pub fn mixing_ratio(temperature_c: f64, pressure_hpa: f64) -> f64 {
    let vp = vapor_pressure_water(temperature_c);
    621.97 * (vp / (pressure_hpa - vp))
}

/// Given a mixing ratio and pressure, calculate the temperature. The p is in hPa and the mw is in
/// g/kg. Assume 100% rh.
pub fn temperature_from_p_and_saturated_mw(p: f64, mw: f64) -> f64 {
    use std::f64;

    let z = mw * p / 6.11 / 621.97 / (1.0 + mw / 621.97);
    237.5 * f64::log10(z) / (7.5 - f64::log10(z))
}

/// Assuming a saturated parcel, calculate the equivalent potential temperature in Kelvin.
pub fn theta_e_saturated_kelvin(pressure_hpa: f64, temperature_c: f64) -> f64 {
    use std::f64;

    let theta = theta_kelvin(pressure_hpa, temperature_c);
    let mw = mixing_ratio(temperature_c, pressure_hpa) / 1000.0; // divide by 1000 to get kg/kg

    theta * f64::exp(2.6897e6 * mw / 1005.7 / (temperature_c + 273.15))
}

/// Convert celsius to fahrenheit.
pub fn celsius_to_f(temperature: f64) -> f64 {
    1.8 * temperature + 32.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use utility::test_tools::*;

    #[test]
    fn test_theta_kelvin() {
        assert!(approx_equal(
            32.0 + 273.15,
            theta_kelvin(1000.0, 32.0),
            1.0e-2
        ));
    }

    // TODO: Test all functions.
}
