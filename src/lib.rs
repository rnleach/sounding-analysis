#![warn(missing_docs)]
/*!
Functions and data types for analyzing soundings from the 
[sounding-base](https://github.com/rnleach/sounding-base.git) crate.

## Purpose
Provides analysis capabilities for the [sounding-base](https://github.com/rnleach/sounding-base.git) 
crate.

*/

//
// API
//
pub use analysis::{Analysis, Index};

pub mod error;
pub use error::*;

pub use interpolation::linear_interpolate;

pub mod layers;
pub use layers::Layer;

pub mod levels;

pub mod parcel;
pub use parcel::Parcel;

pub mod profile;

//
// Internal use only
//

// 3rd party libs
extern crate failure;
#[macro_use]
extern crate failure_derive;
extern crate smallvec;

// framework libs
extern crate metfor;
extern crate sounding_base;

// dev only libs
#[cfg(test)]
extern crate sounding_validate;

// Modules
mod interpolation;
mod analysis;
#[cfg(test)]
mod test_data;

pub(crate) const VEC_SIZE: usize = 2;
