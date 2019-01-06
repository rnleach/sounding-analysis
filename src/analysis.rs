//! Data type and methods for building and describing an analysis.
//!
//! Not every possible analysis is in this data.
use std::collections::HashMap;

use sounding_base::Sounding;

use crate::indexes::{
    haines, haines_high, haines_low, haines_mid, hot_dry_windy, kindex, precipitable_water, swet,
    total_totals,
};
use crate::parcel::{convective_parcel, mixed_layer_parcel, most_unstable_parcel, surface_parcel};
use crate::parcel_profile::{dcape, lift_parcel, partition_cape, ParcelAnalysis, ParcelProfile};
use metfor::{Celsius, CelsiusDiff, JpKg, Length, Mm};
use optional::{none, some, Noned, Optioned};

/// Convenient package for commonly requested analysis values.
///
/// All parcel related values are assumed to be for the 100hPa mixed layer at the surface.
#[derive(Debug, Clone)]
pub struct Analysis {
    // Sounding used to make the analysis
    sounding: Sounding,

    // Profile specific indicies
    swet: Optioned<f64>,
    k_index: Optioned<Celsius>,
    total_totals: Optioned<f64>,
    precipitable_water: Optioned<Mm>,
    convective_t: Optioned<Celsius>,

    // Fire weather indicies
    haines: Optioned<u8>,
    haines_low: Optioned<u8>,
    haines_mid: Optioned<u8>,
    haines_high: Optioned<u8>,
    hdw: Optioned<f64>,
    convective_deficit: Optioned<CelsiusDiff>,
    cape_ratio: Optioned<f64>,

    // Downburst
    dcape: Optioned<JpKg>,
    downrush_t: Optioned<Celsius>,
    downburst_profile: Option<ParcelProfile>,

    // Parcel analysis
    mixed_layer: Option<ParcelAnalysis>,
    surface: Option<ParcelAnalysis>,
    most_unstable: Option<ParcelAnalysis>,
    convective: Option<ParcelAnalysis>,

    // Provider analysis
    provider_analysis: HashMap<&'static str, f64>,
}

impl Analysis {
    /// Create a new `Analysis`.
    pub fn new(snd: Sounding) -> Self {
        Analysis {
            sounding: snd,
            swet: none(),
            k_index: none(),
            total_totals: none(),
            precipitable_water: none(),
            convective_t: none(),

            haines: none(),
            haines_low: none(),
            haines_mid: none(),
            haines_high: none(),
            hdw: none(),
            convective_deficit: none(),
            cape_ratio: none(),

            dcape: none(),
            downrush_t: none(),
            downburst_profile: None,

            mixed_layer: None,
            surface: None,
            most_unstable: None,
            convective: None,

            provider_analysis: HashMap::new(),
        }
    }

    /// Builder method to set the SWeT
    pub fn with_swet<T>(self, value: T) -> Self
    where
        Optioned<f64>: From<T>,
    {
        Self {
            swet: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the K index
    pub fn with_k_index<T>(self, value: T) -> Self
    where
        Optioned<Celsius>: From<T>,
    {
        Self {
            k_index: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the Total Totals
    pub fn with_total_totals<T>(self, value: T) -> Self
    where
        Optioned<f64>: From<T>,
    {
        Self {
            total_totals: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the precipitable water
    pub fn with_pwat<T, U>(self, value: T) -> Self
    where
        Optioned<U>: From<T>,
        U: Length + Noned,
        Mm: From<U>,
    {
        let u: Optioned<U> = Optioned::from(value);
        let pw: Optioned<Mm> = u.map_t(Mm::from);
        Self {
            precipitable_water: pw,
            ..self
        }
    }

    /// Builder method to set the convective temperature
    pub fn with_convective_t<T>(self, value: T) -> Self
    where
        Optioned<Celsius>: From<T>,
    {
        Self {
            convective_t: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the haines index. High, mid, or low Haines should be based on the
    /// station elevation of the sounding.
    pub fn with_haines<T>(self, value: T) -> Self
    where
        Optioned<u8>: From<T>,
    {
        Self {
            haines: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the low level haines index.
    pub fn with_haines_low<T>(self, value: T) -> Self
    where
        Optioned<u8>: From<T>,
    {
        Self {
            haines_low: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the mid level haines index.
    pub fn with_haines_mid<T>(self, value: T) -> Self
    where
        Optioned<u8>: From<T>,
    {
        Self {
            haines_mid: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the high level haines index.
    pub fn with_haines_high<T>(self, value: T) -> Self
    where
        Optioned<u8>: From<T>,
    {
        Self {
            haines_high: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the hot-dry-windy index.
    pub fn with_hdw<T>(self, value: T) -> Self
    where
        Optioned<f64>: From<T>,
    {
        Self {
            hdw: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the wet/dry cape ratio. EXPERIMENTAL
    pub fn with_cape_ratio<T>(self, value: T) -> Self
    where
        Optioned<f64>: From<T>,
    {
        Self {
            cape_ratio: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the convective temperature deficit. EXPERIMENTAL
    pub fn with_convective_deficit<T>(self, value: T) -> Self
    where
        Optioned<CelsiusDiff>: From<T>,
    {
        Self {
            convective_deficit: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the DCAPE.
    pub fn with_dcape<T>(self, value: T) -> Self
    where
        Optioned<JpKg>: From<T>,
    {
        Self {
            dcape: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the downrush temperature of a wet micro-burst.
    pub fn with_downrush_t<T>(self, value: T) -> Self
    where
        Optioned<Celsius>: From<T>,
    {
        Self {
            downrush_t: Optioned::from(value),
            ..self
        }
    }

    /// Get the Swet
    pub fn swet(&self) -> Optioned<f64> {
        self.swet
    }

    /// Get the K index
    pub fn k_index(&self) -> Optioned<Celsius> {
        self.k_index
    }

    /// Get the Total Totals index
    pub fn total_totals(&self) -> Optioned<f64> {
        self.total_totals
    }

    /// Get the precipitable water.
    pub fn pwat(&self) -> Optioned<Mm> {
        self.precipitable_water
    }

    /// Get the convective temperature.
    pub fn convective_t(&self) -> Optioned<Celsius> {
        self.convective_t
    }

    /// Get the downrush temperature from a microburst.
    pub fn downrush_t(&self) -> Optioned<Celsius> {
        self.downrush_t
    }

    /// Get the Haines Index.
    pub fn haines(&self) -> Optioned<u8> {
        self.haines
    }

    /// Get the low level Haines Index.
    pub fn haines_low(&self) -> Optioned<u8> {
        self.haines_low
    }

    /// Get the mid level Haines Index.
    pub fn haines_mid(&self) -> Optioned<u8> {
        self.haines_mid
    }

    /// Get the high level Haines Index.
    pub fn haines_high(&self) -> Optioned<u8> {
        self.haines_high
    }

    /// Get the hot-dry-windy index.
    pub fn hdw(&self) -> Optioned<f64> {
        self.hdw
    }

    /// Get the wet/dry CAPE ratio. EXPERIMENTAL.
    pub fn cape_ratio(&self) -> Optioned<f64> {
        self.hdw
    }

    /// Get the convective temperature deficit.
    pub fn convective_deficit(&self) -> Optioned<CelsiusDiff> {
        self.convective_deficit
    }

    /// Get the DCAPE.
    pub fn dcape(&self) -> Optioned<JpKg> {
        self.dcape
    }

    /// Set the mixed layer parcel analysis.
    pub fn with_mixed_layer_parcel_analysis<T>(self, anal: T) -> Self
    where
        Option<ParcelAnalysis>: From<T>,
    {
        let mixed_layer = Option::from(anal);
        Analysis {
            mixed_layer,
            ..self
        }
    }

    /// Get the mixed layer parcel analysis
    pub fn mixed_layer_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
        self.mixed_layer.as_ref()
    }

    /// Set the surface parcel analysis.
    pub fn with_surface_parcel_analysis<T>(self, anal: T) -> Self
    where
        Option<ParcelAnalysis>: From<T>,
    {
        let surface = Option::from(anal);
        Analysis { surface, ..self }
    }

    /// Get the surface parcel analysis
    pub fn surface_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
        self.surface.as_ref()
    }

    /// Set the most unstable parcel analysis.
    pub fn with_most_unstable_parcel_analysis<T>(self, anal: T) -> Self
    where
        Option<ParcelAnalysis>: From<T>,
    {
        let most_unstable = Option::from(anal);
        Analysis {
            most_unstable,
            ..self
        }
    }

    /// Get the most unstable parcel analysis
    pub fn most_unstable_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
        self.most_unstable.as_ref()
    }

    /// Set the convective parcel analysis.
    pub fn with_convective_parcel_analysis<T>(self, anal: T) -> Self
    where
        Option<ParcelAnalysis>: From<T>,
    {
        let convective = Option::from(anal);
        Analysis { convective, ..self }
    }

    /// Get the convective parcel analysis
    pub fn convective_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
        self.convective.as_ref()
    }

    /// Set the downburst profile
    pub fn with_downburst_profile<T>(self, parcel_profile: T) -> Self
    where
        Option<ParcelProfile>: From<T>,
    {
        let downburst_profile = Option::from(parcel_profile);
        Analysis {
            downburst_profile,
            ..self
        }
    }

    /// Get the downburst profile
    pub fn downburst_profile(&self) -> Option<&ParcelProfile> {
        self.downburst_profile.as_ref()
    }

    /// Set the provider analysis.
    ///
    /// This is just a table of what ever values you want to store, it may be empty.
    pub fn with_provider_analysis(self, provider_analysis: HashMap<&'static str, f64>) -> Self {
        Analysis {
            provider_analysis,
            ..self
        }
    }

    /// Get a reference to the provider analysis so you can query it.
    pub fn provider_analysis(&self) -> &HashMap<&'static str, f64> {
        &self.provider_analysis
    }

    /// Get a mutable reference to the provider analysis so you can modify it.
    pub fn provider_analysis_mut(&mut self) -> &mut HashMap<&'static str, f64> {
        &mut self.provider_analysis
    }

    /// Get a reference to the sounding.
    pub fn sounding(&self) -> &Sounding {
        &self.sounding
    }

    /// Analyze the sounding to get as much information as you can.
    pub fn fill_in_missing_analysis(mut self) -> Self {
        self.swet = self
            .swet
            .or_else(|| Optioned::from(swet(&self.sounding).ok()));
        self.total_totals = self
            .total_totals
            .or_else(|| Optioned::from(total_totals(&self.sounding).ok()));
        self.k_index = self
            .k_index
            .or_else(|| Optioned::from(kindex(&self.sounding).ok()));
        self.precipitable_water = self
            .precipitable_water
            .or_else(|| Optioned::from(precipitable_water(&self.sounding).ok()));

        self.haines = self
            .haines
            .or_else(|| Optioned::from(haines(&self.sounding).ok()));
        self.haines_low = self
            .haines_low
            .or_else(|| Optioned::from(haines_low(&self.sounding).ok()));
        self.haines_mid = self
            .haines_mid
            .or_else(|| Optioned::from(haines_mid(&self.sounding).ok()));
        self.haines_high = self
            .haines_high
            .or_else(|| Optioned::from(haines_high(&self.sounding).ok()));
        self.hdw = self
            .hdw
            .or_else(|| Optioned::from(hot_dry_windy(&self.sounding).ok()));

        if self.dcape.is_none() || self.downrush_t.is_none() || self.downburst_profile.is_none() {
            let result = dcape(&self.sounding);

            if let Ok((pp, dcape, down_t)) = result {
                self.dcape = some(dcape);
                self.downrush_t = some(down_t);
                self.downburst_profile = Some(pp);
            }
        }

        if self.mixed_layer.is_none() {
            self.mixed_layer = match mixed_layer_parcel(&self.sounding) {
                Ok(parcel) => lift_parcel(parcel, &self.sounding).ok(),
                Err(_) => None,
            };
        }
        if self.most_unstable.is_none() {
            self.most_unstable = match most_unstable_parcel(&self.sounding) {
                Ok(parcel) => lift_parcel(parcel, &self.sounding).ok(),
                Err(_) => None,
            };
        }
        if self.surface.is_none() {
            self.surface = match surface_parcel(&self.sounding) {
                Ok(parcel) => lift_parcel(parcel, &self.sounding).ok(),
                Err(_) => None,
            };
        }
        if self.convective.is_none() {
            self.convective = match convective_parcel(&self.sounding) {
                Ok(parcel) => lift_parcel(parcel, &self.sounding).ok(),
                Err(_) => None,
            };
        }

        // Convective T
        if self.convective_t.is_none() {
            self.convective_t = self
                .convective
                .as_ref()
                .map(|parcel_anal| parcel_anal.parcel().temperature)
                .into();
        }

        // Convective deficit
        if self.convective_deficit.is_none() {
            self.convective_deficit = self.convective_t.and_then(|ct| {
                self.mixed_layer
                    .as_ref()
                    .map(|parcel_anal| ct - parcel_anal.parcel().temperature)
                    .into()
            });
        }

        // Cape ratio
        if self.cape_ratio.is_none() {
            self.cape_ratio = self
                .convective
                .as_ref()
                .and_then(|parcel_anal| partition_cape(parcel_anal).ok())
                .map(|(dry, wet)| wet / dry)
                .into();
        }

        self
    }
}
