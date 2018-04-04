//! Error types for the sounding-analysis crate.

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
}

/// Shorthand for results.
pub type Result<T> = ::std::result::Result<T, AnalysisError>;

impl AnalysisError {
    #[allow(missing_docs)]
    pub fn tag(self, file: &str, line: u32) -> Self {
        if cfg!(debug_assertions) {
            println!(
                "\n******\n\n\nError initiated in {} on line {}: {}\n\n\n******",
                file, line, self
            );
        }
        self
    }
}
