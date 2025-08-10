use std::str::FromStr;

use thiserror::Error;

use crate::field_element::FieldElement;

#[derive(Debug, Clone, PartialEq)]
pub enum PolynomialError {
    CoefficientOutOfRange,
}

#[derive(Debug, Clone, Eq, PartialEq, Error)]
#[non_exhaustive]
pub enum ParseFieldElementError {
    #[error("invalid `u64`")]
    ParseU64Error(#[source] <u64 as FromStr>::Err),

    #[error("non-canonical {0} >= {p} == `FieldElement::P`", p = FieldElement::P)]
    NotCanonical(u64),

    #[error(
        "incorrect number of bytes: {0} != {bytes} == `FieldElement::BYTES`",
        bytes = FieldElement::BYTES
    )]
    InvalidNumBytes(usize),
}
