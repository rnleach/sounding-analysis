#![warn(missing_docs)]
//! Functions for analyzing soundings from the
//! [sounding-base](https://github.com/rnleach/sounding-base.git) crate.

#[macro_use]
extern crate error_chain;
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

/// A layer in the atmosphere described by the top and bottom pressures.
#[derive(Debug, Clone, Copy)]
pub struct Layer {
    /// Pressure at the bottom of the layer.
    pub bottom_press: f64,
    /// Pressure at the top of the layer.
    pub top_press: f64,
}

/// A special level, such as the freezing level or wet bulb zero level.
#[derive(Debug, Clone, Copy)]
pub struct Level {
    /// The pressure value at the level.
    pub pressure: f64,
}

pub(crate) const VEC_SIZE: usize = sounding_base::SMALL_VEC_SIZE;

/*
Types:
  Parcel: temperature, pressure, dew-point
  Layer: pressure_bottom, pressure_top
*/
