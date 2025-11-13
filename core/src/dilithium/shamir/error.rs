use thiserror::Error;

/// Result type specialized for Shamir operations.
pub type ShamirResult<T> = Result<T, ShamirError>;

/// Errors originating from the Shamir secret sharing module.
#[derive(Debug, Error)]
pub enum ShamirError {
    #[error("Invalid participant ID: {0}")]
    InvalidParticipantId(usize),
    #[error("Invalid index: {0} >= {1}")]
    InvalidIndex(usize, usize),
    #[error("Inconsistent share lengths")]
    InconsistentShareLengths,
}
