#![warn(missing_docs)]
//! Functions for analyzing soundings from the
//! [sounding-base](https://github.com/rnleach/sounding-base.git) crate.

extern crate sounding_base;

mod interpolation;
pub use interpolation::linear_interpolate;

mod snow;
pub use snow::dendritic_growth_zone;
