//! Enums used as keys for setting options in functions.

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
