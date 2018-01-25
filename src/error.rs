#![allow(missing_docs, unused_doc_comment)]
//! Error types for the sounding-analysis crate.

error_chain!{
    errors {
        /// A profile that is required for this analysis is missing.
        MissingProfile(profile: &'static str, analysis: &'static str) {
            description("Missing required profile."),
            display("Profile for {} is missing, cannot do {} analysis.", profile, analysis),
        }
    }
}
