#![allow(missing_docs, unused_doc_comment)]
//! Error types for the sounding-analysis crate.

error_chain!{
    errors {
        /// A profile that is required for this analysis is missing.
        MissingProfile(profile: &'static str, analysis: &'static str) {
            description("Missing required profile."),
            display("Profile for {} is missing, cannot do {} analysis.", profile, analysis),
        }

        /// A profile that is required for this analysis is missing.
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
