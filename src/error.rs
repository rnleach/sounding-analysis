//! Error types for the sounding-analysis crate.

#[derive(Debug, Fail)]
pub enum AnalysisError {
    /// A profile that is required for this analysis is missing.
    #[fail(display="Missing profile ({}) required for the analysis ({}).", _0, _1)]
    MissingProfile(&'static str, &'static str),
    /// There is no data available that meets the requirements.
    
    /// Bad or invalid input.
}
error_chain!{
    errors {
        /// A profile that is required for this analysis is missing.
        MissingProfile(profile: &'static str, analysis: &'static str) {
            description("Missing required profile."),
            display("Profile for {} is missing, cannot do {} analysis.", profile, analysis),
        }

        /// There is no data available that meets the requirements.
        NoDataProfile(profile: &'static str, analysis: &'static str) {
            description("Profile has all missing values."),
            display("Profile for {} is full of missing values, cannot do {} analysis.", profile, analysis),
        }

        /// Bad or invalid input.
        InvalidInput(msg: &'static str, analysis: &'static str) {
            description("Invalid input to algorithm."),
            display("Invalid input to {}: {}", analysis, msg),
        }
    }
}
