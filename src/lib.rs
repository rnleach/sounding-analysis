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
pub use analysis::{Analysis, ParcelAnalysis, ParcelIndex, ProfileIndex};

pub mod error;
pub use error::*;

pub mod indexes;

pub use interpolation::{linear_interpolate, linear_interpolate_sounding};

pub mod layers;
pub use layers::{Layer, Layers};

pub mod levels;
pub use levels::{Level, Levels};

pub mod parcel;
pub use parcel::{Parcel, ParcelProfile};

pub mod profile;

//
// Internal use only
//

// 3rd party libs
extern crate failure;
#[macro_use]
extern crate failure_derive;
#[macro_use]
extern crate itertools;
extern crate smallvec;

// framework libs
extern crate metfor;
extern crate sounding_base;

// dev only libs
#[cfg(test)]
extern crate sounding_validate;

// Modules
mod analysis;
mod interpolation;

pub(crate) const VEC_SIZE: usize = 2;
