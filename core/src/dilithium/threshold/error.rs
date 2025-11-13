use std::error::Error as StdError;
use std::fmt;

/// Errors specific to the threshold Dilithium module.
//Add thiserror variants as needed.
#[non_exhaustive]
#[derive(Debug)]
pub enum ThresholdError {
    PartialSignatureChallengeMismatch,
}

impl fmt::Display for ThresholdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ThresholdError::PartialSignatureChallengeMismatch => {
                write!(f, "All partial signatures must use the same challenge")
            }
        }
    }
}

impl StdError for ThresholdError {}
