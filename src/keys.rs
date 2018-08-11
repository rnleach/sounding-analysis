//! Enums used as keys for setting options in functions.

/// Sounding indexes calculated from the sounding and not any particular profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProfileIndex {
    /// Severe Weather Threat Index
    SWeT,
    /// K-index
    K,
    /// Precipitable Water (mm)
    PWAT,
    /// Total-Totals
    TotalTotals,
    /// Haines index
    Haines,
    /// Downward CAPE
    DCAPE,
    /// Downrush temperature. The temperature of a saturated downburst from parcel theory.
    DownrushT,
}

/// Indexes from a parcel analysis of a sounding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParcelIndex {
    /// Lifting Condensation Level meters AGL
    LCLHeightAGL,
    /// Lifting Condensation Level, or LCL (hPa), pressure vertical coordinate.
    LCLPressure,
    /// Convective Available Potential Energy, or CAPE. (J/kg)
    CAPE,
    /// CAPE in the hail growth zone
    CAPEHail,
    /// Temperature at LCL (C)
    LCLTemperature,
    /// Convective Inhibitive Energy, or CIN (J/kg)
    CIN,
    /// Equilibrium Level (hPa), pressure vertical coordinate
    ELPressure,
    /// Eqilibrium level height (meters ASL)
    ELHeightASL,
    /// Equilibrium level temperature (degrees C)
    ELTemperature,
    /// Level of Free Convection (hPa), pressure vertical coordinate
    LFC,
    /// Normalized CAPE
    NCAPE,
    /// Lifted Index
    LI,
}
