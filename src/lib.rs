/*!
Functions and data types for analyzing soundings from the
[sounding-base](https://github.com/rnleach/sounding-base.git) crate.

## Purpose
Provides analysis capabilities for the [sounding-base](https://github.com/rnleach/sounding-base.git)
crate.

*/
#![doc(test(attr(deny(warnings))))]
#![warn(missing_docs)]

//
// API
//
pub use crate::{
    analysis::Analysis,
    error::{AnalysisError, Result},
    experimental::PlumePotentialAnal,
    indexes::{
        haines, haines_high, haines_low, haines_mid, hot_dry_windy, kindex, precipitable_water,
        swet, total_totals,
    },
    interpolation::{linear_interpolate, linear_interpolate_sounding},
    layers::{
        cold_surface_temperature_layer, dendritic_snow_zone, effective_inflow_layer,
        hail_growth_zone, inversions, layer_agl, pressure_layer, sfc_based_inversion,
        warm_temperature_layer_aloft, warm_wet_bulb_layer_aloft, Layer, Layers,
    },
    levels::{
        freezing_levels, max_temperature_in_layer, max_temperature_in_profile,
        max_wet_bulb_in_layer, max_wet_bulb_in_profile, wet_bulb_zero_levels, Level, Levels,
    },
    parcel::{
        convective_parcel, lowest_level_parcel, mixed_layer_parcel, most_unstable_parcel,
        pressure_parcel, surface_parcel, Parcel,
    },
    parcel_profile::{
        dcape, lift_parcel, mix_down, partition_cape, robust_convective_parcel, ParcelAnalysis,
        ParcelProfile,
    },
    profile::{
        equivalent_potential_temperature, hydrolapse, potential_temperature, relative_humidity,
        relative_humidity_ice, sfc_to_level_temperature_lapse_rate, temperature_lapse_rate,
        theta_e_lapse_rate, wet_bulb,
    },
    wind::{bunkers_storm_motion, mean_wind, sr_helicity},
};

//
// Internal use only
//

// Modules
mod analysis;
mod error;
mod experimental;
mod indexes;
mod interpolation;
mod layers;
mod levels;
mod parcel;
mod parcel_profile;
mod profile;
mod wind;

pub(crate) const VEC_SIZE: usize = 2;
