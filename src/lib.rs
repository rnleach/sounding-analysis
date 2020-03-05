// FIXME: Add some examples of using analysis functions.
/*!
Functions and data types for analyzing soundings from radiosondes or models.

Library to represent an atmospheric sounding with pressure as the vertical coordinate.

# Examples
```
use optional::{Optioned, some};
use metfor::{HectoPascal, Celsius, Feet};

use sounding_analysis::{Sounding, StationInfo};

fn main() {

    // Create  pressure profile
    let pressure_profile: Vec<Optioned<HectoPascal>> =
        vec![1000.0, 925.0, 850.0, 700.0, 500.0, 300.0, 250.0, 100.0]
            .into_iter()
            .map(HectoPascal)
            .map(some)
            .collect();

    // Create a temperature profile
    let temperature_profile: Vec<Optioned<Celsius>> =
        vec![13.0, 7.0, 5.0, -4.5, -20.6, -44.0, -52.0, -56.5]
            .into_iter()
            .map(Celsius)
            .map(some)
            .collect();

    // Create some station info
    let stn = StationInfo::new_with_values(None, (45.6789, -115.6789), Feet(992.0));

    // Create a valid time. This uses a `chrono::NaiveDateTime`, and you should always assume
    // that valid times are in UTC.
    let vt = chrono::NaiveDate::from_ymd(2018,3,8).and_hms(12,0,0);

    // Use the builder pattern to construct a sounding.
    let snd = Sounding::new()
        .with_station_info(stn)
        .with_valid_time(vt)
        .with_lead_time(24)  // Lead time in hours for forecast soundings.
        .with_pressure_profile(pressure_profile)
        .with_temperature_profile(temperature_profile)
        .with_station_pressure(some(HectoPascal(1013.25)))
        .with_sfc_temperature(some(Celsius(15.0)));

    // Top down and bottom up iterators are provided. If surface data is available, it is
    // inserted into the profile.
    let mut iter = snd.top_down();

    let mut data_row = iter.next().unwrap();
    assert_eq!(data_row.pressure, some(HectoPascal(100.0)));
    assert_eq!(data_row.temperature, some(Celsius(-56.5)));

    data_row = iter.next().unwrap();
    assert_eq!(data_row.pressure, some(HectoPascal(250.0)));
    assert_eq!(data_row.temperature, some(Celsius(-52.0)));

    data_row = iter.next().unwrap();
    assert_eq!(data_row.pressure, some(HectoPascal(300.0)));
    assert_eq!(data_row.temperature, some(Celsius(-44.0)));

    data_row = iter.next().unwrap();
    assert_eq!(data_row.pressure, some(HectoPascal(500.0)));
    assert_eq!(data_row.temperature, some(Celsius(-20.6)));

    data_row = iter.next().unwrap();
    assert_eq!(data_row.pressure, some(HectoPascal(700.0)));
    assert_eq!(data_row.temperature, some(Celsius(-4.5)));

    data_row = iter.next().unwrap();
    assert_eq!(data_row.pressure, some(HectoPascal(850.0)));
    assert_eq!(data_row.temperature, some(Celsius(5.0)));

    data_row = iter.next().unwrap();
    assert_eq!(data_row.pressure, some(HectoPascal(925.0)));
    assert_eq!(data_row.temperature, some(Celsius(7.0)));

    data_row = iter.next().unwrap();
    assert_eq!(data_row.pressure, some(HectoPascal(1000.0)));
    assert_eq!(data_row.temperature, some(Celsius(13.0)));

    // THIS ONE IS THE SURFACE DATA!
    data_row = iter.next().unwrap();
    assert_eq!(data_row.pressure, some(HectoPascal(1013.25)));
    assert_eq!(data_row.temperature, some(Celsius(15.0)));

    assert_eq!(iter.next(), None);

    // Profiles and surface values can also be accessed via getter methods. Read the docs!
}
```

You probably noticed a lot of `optional::Optioned`s in the example. Basically, anything can be
missing, and missing values are common in upper air soundings. For example, at high altitude the
dew point or humidity are often missing (if not totally inaccurate).

*/
#![doc(test(attr(deny(warnings))))]
#![deny(missing_docs)]

//
// API
//
#[allow(deprecated)]
pub use crate::{
    error::{AnalysisError, Result},
    indexes::{haines, haines_high, haines_low, haines_mid, hot_dry_windy, precipitable_water},
    interpolation::{linear_interpolate, linear_interpolate_sounding},
    layers::{
        cold_surface_temperature_layer, dendritic_snow_zone, effective_inflow_layer,
        hail_growth_zone, inversions, layer_agl, melting_freezing_energy_area, pressure_layer,
        sfc_based_inversion, warm_surface_temperature_layer, warm_temperature_layer_aloft,
        warm_wet_bulb_layer_aloft, Layer, Layers,
    },
    levels::{
        freezing_levels, max_temperature_in_layer, max_temperature_in_profile,
        max_wet_bulb_in_layer, max_wet_bulb_in_profile, wet_bulb_zero_levels, Level, Levels,
    },
    parcel::{
        average_parcel, convective_parcel, effective_layer_parcel, lowest_level_parcel,
        mixed_layer_parcel, most_unstable_parcel, pressure_parcel, surface_parcel, Parcel,
    },
    parcel_profile::{
        dcape, lift_parcel, mix_down, robust_convective_parcel_ascent, ParcelAscentAnalysis,
        ParcelProfile,
    },
    precip_type::{
        adjust_precip_type_intensity, bourgouin_precip_type, check_precip_type_intensity,
        PrecipType,
    },
    profile::{
        equivalent_potential_temperature, hydrolapse, potential_temperature, relative_humidity,
        relative_humidity_ice, sfc_to_level_temperature_lapse_rate, temperature_lapse_rate,
        theta_e_lapse_rate, wet_bulb,
    },
    sounding::{DataRow, Sounding, StationInfo},
    wind::{bunkers_storm_motion, mean_wind, sr_helicity},
};

pub mod experimental;

#[doc(hidden)]
pub use crate::sounding::doctest;

//
// Internal use only
//

// Modules
mod error;
mod indexes;
mod interpolation;
mod layers;
mod levels;
mod parcel;
mod parcel_profile;
mod precip_type;
mod profile;
mod sounding;
mod wind;
