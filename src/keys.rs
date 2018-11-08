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
    /// Haines index, whichever version is auto-calculated from the sounding elevation.
    Haines,
    /// Haines index, low version
    HainesLow,
    /// Haines index, mid level version
    HainesMid,
    /// Haines index, high version
    HainesHigh,
    /// Hot-dry-windy index
    Hdw,
    /// Downward CAPE
    DCAPE,
    /// Downrush temperature. The temperature of a saturated downburst from parcel theory.
    DownrushT,
    /// Mixed layer temperatures for the convective parcel.
    ConvectiveT,
    /// Differenence in `ConvectiveT` and the temperature of the mixed layer parcel.
    ConvectiveDeficit,
    /// Ratio of wet/dry cape from the convective parcel analysis.
    CapeRatio,
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
    /// Level of free convection (hPa), pressure vertical coordinate
    LFC,
    /// Level of free convection in the temperature coordinate.
    // Useful for plotting on skew-t
    LFCVirtualTemperature,
    /// Normalized CAPE
    NCAPE,
    /// Lifted Index
    LI,
}
