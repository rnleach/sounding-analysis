//! Data type and methods for building and describing an analysis.
//!
//! Not every possible analysis is in this data.

/// Sounding indexes.
///
/// and silently fail in release mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Index {
    /// Showalter index
    Showalter,
    /// Lifted index
    LI,
    /// Severe Weather Threat Index
    SWeT,
    /// K-index
    K,
    /// Lifting Condensation Level, or LCL (hPa), pressure vertical coordinate.
    LCL,
    /// Precipitable Water (mm)
    PWAT,
    /// Total-Totals
    TotalTotals,
    /// Convective Available Potential Energy, or CAPE. (J/kg)
    CAPE,
    /// Temperature at LCL (K)
    LCLTemperature,
    /// Convective Inhibitive Energy, or CIN (J/kg)
    CIN,
    /// Equilibrium Level (hPa), pressure vertical coordinate
    EquilibrimLevel,
    /// Level of Free Convection (hPa), pressure vertical coordinate
    LFC,
    /// Bulk Richardson Number
    BulkRichardsonNumber,
    /// Haines index
    Haines,
}

/// Convenient package for commonly requested analysis values.
///
/// All parcel related values are assumed to be for the 100hPa mixed layer at the surface.
#[derive(Debug, Clone, Default)]
pub struct Analysis {
    showalter: Option<f64>,
    mixed_layer_lifted_index: Option<f64>,
    swet: Option<f64>,
    k_index: Option<f64>,
    precipitable_water: Option<f64>,
    total_totals: Option<f64>,
    mixed_layer_cape: Option<f64>,
    mixed_layer_lcl: Option<f64>,
    mixed_layer_lcl_temperature: Option<f64>,
    mixed_layer_cin: Option<f64>,
    mixed_layer_equilibrium_level: Option<f64>,
    mixed_layer_lfc: Option<f64>,
    bulk_richardson_number: Option<f64>,
    haines: Option<f64>,
}

impl Analysis {
    /// Create a new empty `Analysis`.
    pub fn new() -> Analysis {
        Analysis::default()
    }

    /// Set a value in the analysis
    pub fn set<T>(self, var: Index, value: T) -> Self
    where
        Option<f64>: From<T>,
    {
        use Index::*;

        let opt = Option::from(value);

        match var {
            Showalter => Analysis {
                showalter: opt,
                ..self
            },
            LI => Analysis {
                mixed_layer_lifted_index: opt,
                ..self
            },
            SWeT => Analysis { swet: opt, ..self },
            K => Analysis {
                k_index: opt,
                ..self
            },
            LCL => Analysis {
                mixed_layer_lcl: opt,
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
            CAPE => Analysis {
                mixed_layer_cape: opt,
                ..self
            },
            LCLTemperature => Analysis {
                mixed_layer_lcl_temperature: opt,
                ..self
            },
            CIN => Analysis {
                mixed_layer_cin: opt,
                ..self
            },
            EquilibrimLevel => Analysis {
                mixed_layer_equilibrium_level: opt,
                ..self
            },
            LFC => Analysis {
                mixed_layer_lfc: opt,
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
    pub fn get(&self, var: Index) -> Option<f64> {
        use Index::*;

        match var {
            Showalter => self.showalter,
            LI => self.mixed_layer_lifted_index,
            SWeT => self.swet,
            K => self.k_index,
            LCL => self.mixed_layer_lcl,
            PWAT => self.precipitable_water,
            TotalTotals => self.total_totals,
            CAPE => self.mixed_layer_cape,
            LCLTemperature => self.mixed_layer_lcl_temperature,
            CIN => self.mixed_layer_cin,
            EquilibrimLevel => self.mixed_layer_equilibrium_level,
            LFC => self.mixed_layer_lfc,
            BulkRichardsonNumber => self.bulk_richardson_number,
            Haines => self.haines,
        }
    }
}
