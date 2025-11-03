use thiserror::Error;

/// Result type specialized for threshold operations.
pub type ThresholdResult<T> = std::result::Result<T, ThresholdError>;

/// Errors that can arise while executing the threshold Dilithium protocol.
#[non_exhaustive]
#[derive(Debug, Error)]
pub enum ThresholdError {
    #[error(
        "Invalid threshold configuration: threshold {0} exceeds participant count {1}"
    )]
    InvalidThreshold(usize, usize),
    #[error("Insufficient shares: need {0}, got {1}")]
    InsufficientShares(usize, usize),
    #[error("Invalid participant ID: {0}")]
    InvalidParticipantId(usize),
    #[error("Inconsistent share lengths")]
    InconsistentShareLengths,
    #[error("Modular inverse does not exist")]
    ModularInverseError,
    #[error("Signature generation failed after maximum attempts")]
    SignatureGenerationFailed,
    #[error("All partial signatures must use the same challenge")]
    PartialSignatureChallengeMismatch,
    #[error("Invalid signature bounds")]
    InvalidSignatureBounds,
    #[error("Invalid polynomial index: {0} >= {1}")]
    InvalidIndex(usize, usize),
}
