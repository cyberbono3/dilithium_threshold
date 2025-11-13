use std::error::Error as StdError;
use std::fmt;

/// Result type specialized for Shamir operations.
pub type ShamirResult<T> = Result<T, ShamirError>;

/// Errors originating from the Shamir secret sharing module.
// Add thiserror variants as needed.
#[non_exhaustive]
#[derive(Debug)]
pub enum ShamirError {
    InvalidParticipantId(usize),
    InvalidIndex(usize, usize),
    InconsistentShareLengths,
}

impl fmt::Display for ShamirError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ShamirError::InvalidParticipantId(id) => {
                write!(f, "Invalid participant ID: {id}")
            }
            ShamirError::InvalidIndex(idx, len) => {
                write!(f, "Invalid index: {idx} >= {len}")
            }
            ShamirError::InconsistentShareLengths => {
                write!(f, "Inconsistent share lengths")
            }
        }
    }
}

impl StdError for ShamirError {}
