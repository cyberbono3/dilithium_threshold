use thiserror::Error;

/// Errors specific to the threshold Dilithium module.
#[derive(Debug, Error)]
pub enum ThresholdError {
    #[error("All partial signatures must use the same challenge")]
    PartialSignatureChallengeMismatch,
}
