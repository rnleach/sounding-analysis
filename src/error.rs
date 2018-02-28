//! Error types for the sounding-analysis crate.

/// Error type for the crate.
#[derive(Debug, Fail)]
pub enum AnalysisError {
    /// A profile that is required for this analysis is missing.
    #[fail(display = "Missing profile ({}) required for the analysis ({}).", _0, _1)]
    MissingProfile(&'static str, &'static str),
    /// A value (surface value, index, location, etc) that is required is not available.
    #[fail(display = "Missing value ({}) required for analysis ({}).", _0, _1)]
    MissingValue(&'static str, &'static str),
    /// Not enough data available for anlaysis
    #[fail(display = "Not enough data available for analysis ({}).", _0)]
    NotEnoughData(&'static str),
    /// There is no data available that meets the requirements.
    #[fail(display = "Profile for {} is full of missing values, cannot do {} analysis.", _0, _1)]
    NoDataProfile(&'static str, &'static str),
    /// Bad or invalid input.
    #[fail(display = "Invalid input to {}: {}", _0, _1)]
    InvalidInput(&'static str, &'static str),
}

/// Shorthand for results.
pub type Result<T> = ::std::result::Result<T, AnalysisError>;
