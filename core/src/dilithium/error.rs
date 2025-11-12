use crate::dilithium::shamir::error::ShamirError;
use crate::dilithium::threshold::error::ThresholdError;
use math::error::MatrixError;
use std::error::Error as StdError;
use std::fmt;

/// Result type specialized for Dilithium operations.
pub type DilithiumResult<T> = std::result::Result<T, DilithiumError>;

/// Errors that can arise while executing the threshold Dilithium protocol.
#[non_exhaustive]
#[derive(Debug)]
pub enum DilithiumError {
    InvalidThreshold(usize, usize),
    InsufficientShares(usize, usize),
    Shamir(ShamirError),
    SignatureGenerationFailed,
    Threshold(ThresholdError),
    Matrix(MatrixError),
}

impl fmt::Display for DilithiumError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DilithiumError::InvalidThreshold(threshold, participants) => {
                write!(
                    f,
                    "Invalid threshold configuration: threshold {} exceeds participant count {}",
                    threshold, participants
                )
            }
            DilithiumError::InsufficientShares(expected, actual) => {
                write!(f, "Insufficient shares: need {expected}, got {actual}")
            }
            DilithiumError::Shamir(err) => err.fmt(f),
            DilithiumError::SignatureGenerationFailed => {
                write!(f, "Signature generation failed after maximum attempts")
            }
            DilithiumError::Threshold(err) => err.fmt(f),
            DilithiumError::Matrix(err) => err.fmt(f),
        }
    }
}

impl StdError for DilithiumError {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        match self {
            DilithiumError::Shamir(err) => Some(err),
            DilithiumError::Threshold(err) => Some(err),
            DilithiumError::Matrix(err) => Some(err),
            _ => None,
        }
    }
}

impl From<ShamirError> for DilithiumError {
    fn from(value: ShamirError) -> Self {
        DilithiumError::Shamir(value)
    }
}

impl From<ThresholdError> for DilithiumError {
    fn from(value: ThresholdError) -> Self {
        DilithiumError::Threshold(value)
    }
}

impl From<MatrixError> for DilithiumError {
    fn from(value: MatrixError) -> Self {
        DilithiumError::Matrix(value)
    }
}
