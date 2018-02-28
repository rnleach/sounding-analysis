//! Functions for doing parcel analysis on a sounding, specifically related to convection.
//!

/// Variables defining a parcel as used in parcel analysis.
#[derive(Debug, Clone, Copy)]
pub struct Parcel {
    /// Temperature in C
    pub temperature: f64,
    /// Pressure in hPa
    pub pressure: f64,
    /// Dew point in C
    pub dew_point: f64,
}

// TODO: cape, cin, el, lfc, lcl, dcape, ncape, hail zone cape
