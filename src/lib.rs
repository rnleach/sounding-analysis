#![warn(missing_docs)]
//! Functions and data types for analyzing soundings from the
//! [sounding-base](https://github.com/rnleach/sounding-base.git) crate.

//
// API
//
pub use interpolation::linear_interpolate;

pub mod layers;
pub mod levels;
pub mod met_formulas;
pub mod parcel;
pub mod profile;

/// Utility functions and types. This will eventually be removed as it doesn't really belong in
/// this crate's API.
pub mod utility;

/// A layer in the atmosphere described by the top and bottom pressures.
#[derive(Debug, Clone, Copy)]
pub struct Layer {
    /// Pressure at the bottom of the layer.
    pub bottom: DataRow,
    /// Pressure at the top of the layer.
    pub top: DataRow,
}

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

//
// Internal use only
//

// 3rd party libs
extern crate failure;
#[macro_use]
extern crate failure_derive;
extern crate smallvec;

// framework libs
extern crate sounding_base;

// dev only libs
#[cfg(test)]
extern crate sounding_validate;

// Modules
mod error;
mod interpolation;
#[cfg(test)]
mod test_data;

// Internal use only
pub(crate) const VEC_SIZE: usize = sounding_base::SMALL_VEC_SIZE;

use sounding_base::DataRow;
