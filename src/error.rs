//! Error types for the sounding-analysis crate.
use metfor;
use std::fmt::Display;

/// Shorthand for results.
pub type Result<T> = ::std::result::Result<T, AnalysisError>;

/// Error type for the crate.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum AnalysisError {
    /// A profile that is required for this analysis is missing.
    MissingProfile,
    /// A value (surface value, index, location, etc) that is required is not available.
    MissingValue,
    /// Not enough data available for anlaysis
    NotEnoughData,
    /// There is no data available that meets the requirements.
    NoDataProfile,
    /// Bad or invalid input.
    InvalidInput,
    /// Missing data during interpolation, or it would have been extrapolation
    InterpolationError,
    /// Forward an error from the metfor crate
    MetForError(metfor::MetForErr),
}

impl Display for AnalysisError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        use crate::AnalysisError::*;
        match self {
            MissingProfile => write!(f, "missing profile required for the analysis"),
            MissingValue => write!(f, "missing value required for analysis"),
            NotEnoughData => write!(f, "not enough data available for analysis"),
            NoDataProfile => write!(f, "profile is full of missing values, cannot do analysis"),
            InvalidInput => write!(f, "invalid input"),
            InterpolationError => write!(f, "none value encountered during interpolation"),
            MetForError(err) => write!(f, "error bubbled up from metfor crate: {}", err),
        }
    }
}

impl From<metfor::MetForErr> for AnalysisError {
    fn from(err: metfor::MetForErr) -> Self {
        AnalysisError::MetForError(err)
    }
}
