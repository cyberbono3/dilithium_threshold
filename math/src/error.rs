use std::str::FromStr;

use thiserror::Error;

use crate::field_element::FieldElement;

/// Common result type used across this crate.
pub type Result<T, E = Error> = core::result::Result<T, E>;

/// Top-level error type to keep error management simple for users.
#[derive(Debug, Clone, Eq, PartialEq, Error)]
#[non_exhaustive]
pub enum Error {
    #[error(transparent)]
    ParseFieldElement(#[from] ParseFieldElementError),
    #[error(transparent)]
    Ntt(#[from] NttError),
    #[error(transparent)]
    Polynomial(PolynomialError),
    #[error(transparent)]
    Matrix(#[from] MatrixError),
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum PolynomialError {
    #[error("coefficient out of range")]
    CoefficientOutOfRange,
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum MatrixError {
    #[error("matrix cannot be empty during {operation}")]
    Empty { operation: &'static str },
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
    #[error(
        "matrix/vector shape mismatch during {operation}: matrix is {matrix_rows}x{matrix_cols} but vector length is {vector_len}"
    )]
    VectorShapeMismatch {
        operation: &'static str,
        matrix_rows: usize,
        matrix_cols: usize,
        vector_len: usize,
    },
}

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
