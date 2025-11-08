use std::str::FromStr;

use thiserror::Error;

use crate::field_element::FieldElement;

pub mod polynomial {
    use thiserror::Error;

    #[derive(Debug, Clone, PartialEq, Eq, Error)]
    pub enum Error {
        #[error("coefficient out of range")]
        CoefficientOutOfRange,
    }
}

pub mod matrix {
    use thiserror::Error;

    #[derive(Debug, Clone, PartialEq, Eq, Error)]
    #[non_exhaustive]
    pub enum Error {
        #[error("Matrix cannot be empty")]
        Empty,
        #[error("matrix is ragged: row {row} has {found} columns but expected {expected}")]
        Ragged {
            row: usize,
            expected: usize,
            found: usize,
        },
        #[error(
            "Matrix multiplication produced {found} entries but expected {expected}"
        )]
        OutputLengthMismatch { expected: usize, found: usize },
        #[error("Matrix columns ({matrix_cols}) must match vector length ({vector_len})")]
        VectorShapeMismatch {
            matrix_cols: usize,
            vector_len: usize,
        },
    }
}

pub use matrix::Error as MatrixError;
pub use polynomial::Error as PolynomialError;

/// Common result type used across this crate.
pub type Result<T, E = MathError> = core::result::Result<T, E>;

/// Top-level error type to keep error management simple for users.
#[derive(Debug, Clone, Eq, PartialEq, Error)]
#[non_exhaustive]
pub enum MathError {
    #[error(transparent)]
    ParseFieldElement(#[from] ParseFieldElementError),
    #[error(transparent)]
    Ntt(#[from] NttError),
    #[error("polynomial error: {0:?}")]
    Polynomial(PolynomialError),
    #[error(transparent)]
    Matrix(#[from] MatrixError),
}

pub type Error = MathError;

/// Errors returned by NTT/INTT helpers.
#[derive(Debug, Clone, Eq, PartialEq, Error)]
#[non_exhaustive]
pub enum NttError {
    #[error("length must be a power of two, got {0}")]
    NonPowerOfTwo(usize),
    #[error("length {0} exceeds u32::MAX")]
    TooLarge(usize),
    #[error("missing primitive root of unity of order {0}")]
    MissingPrimitiveRoot(u32),
}

#[derive(Debug, Clone, Eq, PartialEq, Error)]
#[non_exhaustive]
pub enum ParseFieldElementError {
    #[error("invalid `u64`")]
    ParseU64Error(#[source] <u64 as FromStr>::Err),
    #[error(
        "incorrect number of bytes: {0} != {bytes} == `FieldElement::BYTES`",
        bytes = FieldElement::BYTES
    )]
    InvalidNumBytes(usize),
    #[error("non-canonical {0} >= {p} == `FieldElement::P`", p = FieldElement::P)]
    NotCanonical(u64),
}
