#![warn(missing_docs)]
//! Functions for analyzing soundings from the
//! [sounding-base](https://github.com/rnleach/sounding-base.git) crate.

extern crate failure;
#[macro_use] extern crate failure_derive;
extern crate smallvec;
extern crate sounding_base;

mod error;
mod interpolation;
pub use interpolation::linear_interpolate;

pub mod layers;
pub mod levels;
pub mod met_formulas;
pub mod parcel;
pub mod profile;

/// Utility functions and types. This will eventually be deprecated as it doesn't really belong in
/// this crate's API.
pub mod utility;

use sounding_base::DataRow;

/// A layer in the atmosphere described by the top and bottom pressures.
#[derive(Debug, Clone, Copy)]
pub struct Layer {
    /// Pressure at the bottom of the layer.
    pub bottom: DataRow,
    /// Pressure at the top of the layer.
    pub top: DataRow,
}

/// A special level, such as the freezing level or wet bulb zero level.
#[derive(Debug, Clone, Copy)]
pub struct Level {
    /// The pressure value at the level.
    pub pressure: DataRow,
}

pub(crate) const VEC_SIZE: usize = sounding_base::SMALL_VEC_SIZE;

/*
Types:
  Parcel: temperature, pressure, dew-point
  Layer: pressure_bottom, pressure_top
*/
