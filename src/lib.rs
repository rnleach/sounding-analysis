#![warn(missing_docs)]
//! Functions for analyzing soundings from the
//! [sounding-base](https://github.com/rnleach/sounding-base.git) crate.

extern crate sounding_base;

mod interpolation;
pub use interpolation::linear_interpolate;

pub mod met_formulas;

mod snow;
pub use snow::dendritic_growth_zone;

pub mod utility;

/*

Types:
  Parcel: v_coord, temperature
  IndexesExt: Container for indexes not already in the `Sounding` data type

Take a `&mut Sounding` and fill in any missing information. Check element by element in the sounding
for missing fields, and try to fill them in. 

Not all analysis information will be in the `Sounding` structure. Or should it?

Indexes data structure - to return.

Create a parcel profile type.

Analyze cape-cin.
*/
