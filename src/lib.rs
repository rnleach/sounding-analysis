/*!
Functions and data types for analyzing soundings from the 
[sounding-base](https://github.com/rnleach/sounding-base.git) crate.

## Purpose
Provides analysis capabilities for the [sounding-base](https://github.com/rnleach/sounding-base.git) 
crate.

*/
#![warn(missing_docs)]

//
// API
//
pub use analysis::Analysis;
pub use error::{AnalysisError, Result};
pub use indexes::{haines, kindex, precipitable_water, swet, total_totals, hot_dry_windy};
pub use interpolation::{linear_interpolate, linear_interpolate_sounding};
pub use keys::{ParcelIndex, ProfileIndex};
pub use layers::{
    cold_surface_temperature_layer, dendritic_snow_zone, hail_growth_zone, inversions, layer_agl,
    pressure_layer, sfc_based_inversion, warm_temperature_layer_aloft, warm_wet_bulb_layer_aloft,
    Layer, Layers,
};
pub use levels::{
    freezing_levels, max_temperature_in_layer, max_temperature_in_profile, max_wet_bulb_in_layer,
    max_wet_bulb_in_profile, wet_bulb_zero_levels, Level, Levels,
};
pub use parcel::{
    mixed_layer_parcel, most_unstable_parcel, pressure_parcel, surface_parcel, Parcel,
};
pub use parcel_profile::{dcape, lift_parcel, mix_down, ParcelAnalysis, ParcelProfile};
pub use profile::{
    equivalent_potential_temperature, hydrolapse, potential_temperature, relative_humidity,
    sfc_to_level_temperature_lapse_rate, temperature_lapse_rate, theta_e_lapse_rate, wet_bulb,
};

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
extern crate optional;
extern crate sounding_base;

// dev only libs
#[cfg(test)]
extern crate sounding_validate;

// Modules
mod analysis;
mod error;
mod indexes;
mod interpolation;
mod keys;
mod layers;
mod levels;
mod parcel;
mod parcel_profile;
mod profile;

pub(crate) const VEC_SIZE: usize = 2;
