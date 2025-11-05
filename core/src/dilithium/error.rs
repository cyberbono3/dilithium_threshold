use crate::dilithium::shamir::error::ShamirError;
use crate::dilithium::threshold::error::ThresholdError;
use math::error::MatrixError;
use thiserror::Error;

/// Result type specialized for Dilithium operations.
pub type DilithiumResult<T> = std::result::Result<T, DilithiumError>;

/// Errors that can arise while executing the threshold Dilithium protocol.
#[non_exhaustive]
#[derive(Debug, Error)]
pub enum DilithiumError {
    #[error(
        "Invalid threshold configuration: threshold {0} exceeds participant count {1}"
    )]
    InvalidThreshold(usize, usize),
    #[error("Insufficient shares: need {0}, got {1}")]
    InsufficientShares(usize, usize),
    #[error(transparent)]
    Shamir(#[from] ShamirError),
    #[error("Signature generation failed after maximum attempts")]
    SignatureGenerationFailed,
    #[error(transparent)]
    Threshold(#[from] ThresholdError),
    #[error(transparent)]
    Matrix(#[from] MatrixError),
}
