//! Data type and methods for building and describing an analysis.
//!
//! Not every possible analysis is in this data.
use std::collections::HashMap;

use sounding_base::Sounding;

use indexes::{haines, kindex, precipitable_water, swet, total_totals};
use keys::ProfileIndex;
use parcel::{convective_parcel, mixed_layer_parcel, most_unstable_parcel, surface_parcel};
use parcel_profile::{dcape, lift_parcel, ParcelAnalysis, ParcelProfile};

/// Convenient package for commonly requested analysis values.
///
/// All parcel related values are assumed to be for the 100hPa mixed layer at the surface.
#[derive(Debug, Clone)]
pub struct Analysis {
    // Sounding used to make the analysis
    sounding: Sounding,

    // Profile specific indicies
    swet: Option<f64>,
    k_index: Option<f64>,
    precipitable_water: Option<f64>,
    total_totals: Option<f64>,
    haines: Option<f64>,

    //
    // TODO: Add all three of the haines indexes?
    //

    // Downburst
    dcape: Option<f64>,
    downrush_t: Option<f64>,
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
            swet: None,
            k_index: None,
            precipitable_water: None,
            total_totals: None,
            haines: None,

            dcape: None,
            downrush_t: None,
            downburst_profile: None,

            mixed_layer: None,
            surface: None,
            most_unstable: None,
            convective: None,

            provider_analysis: HashMap::new(),
        }
    }

    /// Set a value in the analysis
    pub fn with_profile_index<T>(self, var: ProfileIndex, value: T) -> Self
    where
        Option<f64>: From<T>,
    {
        use self::ProfileIndex::*;

        let opt = Option::from(value);

        match var {
            SWeT => Analysis { swet: opt, ..self },
            K => Analysis {
                k_index: opt,
                ..self
            },
            PWAT => Analysis {
                precipitable_water: opt,
                ..self
            },
            TotalTotals => Analysis {
                total_totals: opt,
                ..self
            },
            DCAPE => Analysis { dcape: opt, ..self },
            DownrushT => Analysis {
                downrush_t: opt,
                ..self
            },
            Haines => Analysis {
                haines: opt,
                ..self
            },
        }
    }

    /// Method to retrieve value from analysis.
    pub fn get_profile_index(&self, var: ProfileIndex) -> Option<f64> {
        use self::ProfileIndex::*;

        match var {
            SWeT => self.swet,
            K => self.k_index,
            PWAT => self.precipitable_water,
            TotalTotals => self.total_totals,
            DCAPE => self.dcape,
            DownrushT => self.downrush_t,
            Haines => self.haines,
        }
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
    pub fn get_mixed_layer_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
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
    pub fn get_surface_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
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
    pub fn get_most_unstable_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
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
    pub fn get_convective_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
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
    pub fn get_downburst_profile(&self) -> Option<&ParcelProfile> {
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
        self.swet = self.swet.or_else(|| swet(&self.sounding).ok());
        self.total_totals = self
            .total_totals
            .or_else(|| total_totals(&self.sounding).ok());
        self.k_index = self.k_index.or_else(|| kindex(&self.sounding).ok());
        self.precipitable_water = self
            .precipitable_water
            .or_else(|| precipitable_water(&self.sounding).ok());
        self.haines = self.haines.or_else(|| haines(&self.sounding).ok());
        if self.dcape.is_none() || self.downrush_t.is_none() || self.downburst_profile.is_none() {
            let result = dcape(&self.sounding);

            if let Ok((pp, dcape, down_t)) = result {
                self.dcape = Some(dcape);
                self.downrush_t = Some(down_t);
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

        self
    }
}
