use std::convert::TryFrom;


#[derive(Debug, Clone, PartialEq)]
pub enum PolynomialError {
    CoefficientOutOfRange,
}

