use std::fmt;

#[derive(Debug)]
pub enum ThresholdError {
    InvalidThreshold {
        threshold: usize,
        participant_number: usize,
    },
    InsufficientShares {
        required: usize,
        provided: usize,
    },
    InvalidParticipantId(usize),
    InconsistentShareLengths,
    ModularInverseError,
    SignatureGenerationFailed,
    PartialSignatureChallengeMismatch,
    InvalidSignatureBounds,
    InvalidIndex {
        index: usize,
        length: usize,
    },
}

impl fmt::Display for ThresholdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ThresholdError::InvalidThreshold {
                threshold,
                participant_number,
            } => {
                write!(f, "Invalid threshold configuration: threshold {} > participant_number {}", 
                    threshold, participant_number)
            }
            ThresholdError::InsufficientShares { required, provided } => {
                write!(
                    f,
                    "Insufficient shares: need {}, got {}",
                    required, provided
                )
            }
            ThresholdError::InvalidParticipantId(id) => {
                write!(f, "Invalid participant ID: {}", id)
            }
            ThresholdError::InconsistentShareLengths => {
                write!(f, "Inconsistent share lengths")
            }
            ThresholdError::ModularInverseError => {
                write!(f, "Modular inverse does not exist")
            }
            ThresholdError::SignatureGenerationFailed => {
                write!(f, "Signature generation failed after maximum attempts")
            }
            ThresholdError::PartialSignatureChallengeMismatch => {
                write!(f, "All partial signatures must use the same challenge")
            }
            ThresholdError::InvalidSignatureBounds => {
                write!(f, "Invalid signature bounds")
            }
            ThresholdError::InvalidIndex { index, length } => {
                write!(f, "Invalid polynomial index: {} >= {}", index, length)
            }
        }
    }
}

impl std::error::Error for ThresholdError {}

pub type Result<T> = std::result::Result<T, ThresholdError>;
