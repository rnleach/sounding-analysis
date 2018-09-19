//! Error types for the sounding-analysis crate.
use metfor;

/// Error type for the crate.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Fail)]
pub enum AnalysisError {
    /// A profile that is required for this analysis is missing.
    #[fail(display = "Missing profile required for the analysis.")]
    MissingProfile,
    /// A value (surface value, index, location, etc) that is required is not available.
    #[fail(display = "Missing value required for analysis.")]
    MissingValue,
    /// Not enough data available for anlaysis
    #[fail(display = "Not enough data available for analysis.")]
    NotEnoughData,
    /// There is no data available that meets the requirements.
    #[fail(display = "Profile is full of missing values, cannot do analysis.")]
    NoDataProfile,
    /// Bad or invalid input.
    #[fail(display = "Invalid input.")]
    InvalidInput,
    /// Missing data during interpolation, or it would have been extrapolation
    #[fail(display = "None value encountered during interpolation.")]
    InterpolationError,

    /// Forward an error from the metfor crate
    #[fail(display = "Error bubbled up from metfor crate: {}", _0)]
    MetForError(metfor::MetForErr),
}

/// Shorthand for results.
pub type Result<T> = ::std::result::Result<T, AnalysisError>;

impl From<metfor::MetForErr> for AnalysisError {
    fn from(err: metfor::MetForErr) -> Self {
        AnalysisError::MetForError(err)
    }
}
