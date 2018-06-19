//! Data type and methods for building and describing an analysis.
//!
//! Not every possible analysis is in this data.
use std::collections::HashMap;

use sounding_base::Sounding;

use error::*;
use indexes::{haines, kindex, parcel_lifted_index, precipitable_water, showalter_index, swet,
              total_totals};
use parcel::{Parcel, ParcelProfile, lift_parcel, mixed_layer_parcel, surface_parcel, most_unstable_parcel};

/// Sounding indexes calculated from the sounding and not any particular profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProfileIndex {
    /// Showalter index
    Showalter,
    /// Severe Weather Threat Index
    SWeT,
    /// K-index
    K,
    /// Precipitable Water (mm)
    PWAT,
    /// Total-Totals
    TotalTotals,
    /// Bulk Richardson Number
    BulkRichardsonNumber,
    /// Haines index
    Haines,
    // TODO: DCAPE
}

/// Indexes from a parcel analysis of a sounding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParcelIndex {
    /// Lifted index
    LI,
    /// Lifting Condensation Level, or LCL (hPa), pressure vertical coordinate.
    LCL,
    /// Convective Available Potential Energy, or CAPE. (J/kg)
    CAPE,
    /// Temperature at LCL (K)
    LCLTemperature,
    /// Convective Inhibitive Energy, or CIN (J/kg)
    CIN,
    /// Equilibrium Level (hPa), pressure vertical coordinate
    EquilibriumLevel,
    /// Level of Free Convection (hPa), pressure vertical coordinate
    LFC,
    // TODO: NCAPE, hail zone cape
}

/// Convenient package for commonly requested analysis values.
///
/// All parcel related values are assumed to be for the 100hPa mixed layer at the surface.
#[derive(Debug, Clone)]
pub struct Analysis {
    // Sounding used to make the analysis
    sounding: Sounding,

    // Profile specific indicies
    showalter: Option<f64>,
    swet: Option<f64>,
    k_index: Option<f64>,
    precipitable_water: Option<f64>,
    total_totals: Option<f64>,
    bulk_richardson_number: Option<f64>,
    haines: Option<f64>,

    // Parcel analysis
    mixed_layer: Option<ParcelAnalysis>,
    surface: Option<ParcelAnalysis>,
    most_unstable: Option<ParcelAnalysis>,

    // Provider analysis
    provider_analysis: HashMap<&'static str, f64>,
}

impl Analysis {
    /// Create a new `Analysis`.
    pub fn new(snd: Sounding) -> Self {
        Analysis {
            sounding: snd,
            showalter: None,
            swet: None,
            k_index: None,
            precipitable_water: None,
            total_totals: None,
            bulk_richardson_number: None,
            haines: None,

            mixed_layer: None,
            surface: None,
            most_unstable: None,

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
            Showalter => Analysis {
                showalter: opt,
                ..self
            },
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
            BulkRichardsonNumber => Analysis {
                bulk_richardson_number: opt,
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
            Showalter => self.showalter,
            SWeT => self.swet,
            K => self.k_index,
            PWAT => self.precipitable_water,
            TotalTotals => self.total_totals,
            BulkRichardsonNumber => self.bulk_richardson_number,
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

    /// Get the surface parcel analysis
    pub fn get_most_unstable_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
        self.most_unstable.as_ref()
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
        self.showalter = self.showalter
            .or_else(|| showalter_index(&self.sounding).ok());
        self.swet = self.swet.or_else(|| swet(&self.sounding).ok());
        self.total_totals = self.total_totals
            .or_else(|| total_totals(&self.sounding).ok());
        self.k_index = self.k_index.or_else(|| kindex(&self.sounding).ok());
        self.precipitable_water = self.precipitable_water
            .or_else(|| precipitable_water(&self.sounding).ok());
        self.haines = self.haines.or_else(|| haines(&self.sounding).ok());
        // TODO: bulk richardson number

        if self.mixed_layer.is_none(){
            self.mixed_layer = match mixed_layer_parcel(&self.sounding){
                Ok(parcel) => ParcelAnalysis::create(parcel, &self.sounding).ok(),
                Err(_) => None,
            };
        }
        if self.most_unstable.is_none(){
            self.most_unstable = match most_unstable_parcel(&self.sounding){
                Ok(parcel) => ParcelAnalysis::create(parcel, &self.sounding).ok(),
                Err(_) => None,
            };
        }
        if self.surface.is_none(){
            self.surface = match surface_parcel(&self.sounding){
                Ok(parcel) => ParcelAnalysis::create(parcel, &self.sounding).ok(),
                Err(_) => None,
            };
        }

        self
    }
}

/// Parcel analysis, this is a way to package the analysis of a parcel.
#[derive(Debug, Clone)]
pub struct ParcelAnalysis {
    // The orginal parcel and profile
    parcel: Parcel,
    profile: ParcelProfile,

    // Indicies from analysis
    li: Option<f64>,
    cape: Option<f64>,
    lcl: Option<f64>,
    lcl_temperature: Option<f64>,
    cin: Option<f64>,
    el: Option<f64>,
    lfc: Option<f64>,
}

impl ParcelAnalysis {
    /// Create a new empty `Analysis`.
    pub fn new(parcel: Parcel, profile: ParcelProfile) -> Self {
        ParcelAnalysis {
            parcel,
            profile,
            li: None,
            cape: None,
            lcl: None,
            lcl_temperature: None,
            cin: None,
            el: None,
            lfc: None,
        }
    }

    /// Create a new analysis and fill it by doing all the parcel analysis'
    #[inline]
    pub fn create(parcel: Parcel, snd: &Sounding) -> Result<Self> {
        let profile = lift_parcel(parcel, snd)?;
        // FIXME: finish building these!
        let li = parcel_lifted_index(&profile).ok();
        let cape = None;
        let lcl = None;
        let lcl_temperature = None;
        let cin = None;
        let el = None;
        let lfc = None;

        Ok(ParcelAnalysis{
            parcel, profile, li, cape, lcl, lcl_temperature, cin, el, lfc,
        })
    }

    /// Set a value in the analysis
    #[inline]
    pub fn set_index<T>(self, var: ParcelIndex, value: T) -> Self
    where
        Option<f64>: From<T>,
    {
        use self::ParcelIndex::*;

        let opt = Option::from(value);

        match var {
            LI => ParcelAnalysis { li: opt, ..self },
            LCL => ParcelAnalysis { lcl: opt, ..self },
            LCLTemperature => ParcelAnalysis {
                lcl_temperature: opt,
                ..self
            },
            CAPE => ParcelAnalysis { cape: opt, ..self },
            CIN => ParcelAnalysis { cin: opt, ..self },
            EquilibriumLevel => ParcelAnalysis { el: opt, ..self },
            LFC => ParcelAnalysis { lfc: opt, ..self },
        }
    }

    /// Method to retrieve value from analysis.
    #[inline]
    pub fn get_index(&self, var: ParcelIndex) -> Option<f64> {
        use self::ParcelIndex::*;

        match var {
            LI => self.li,
            LCL => self.lcl,
            LCLTemperature => self.lcl_temperature,
            CAPE => self.cape,
            CIN => self.cin,
            EquilibriumLevel => self.el,
            LFC => self.lfc,
        }
    }
}
