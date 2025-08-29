pub use crate::{fe, fe_array, fe_vec, poly, poly_vec};
pub use crate::{
    field_element::FieldElement,
    ntt::{intt, ntt, try_intt, try_ntt, Transform},
    poly_vector::PolynomialVector,
    polynomial::Polynomial,
};

/// Dilithium prime modulus (alias to the single source of truth).
pub const Q: i32 = FieldElement::P as i32;
pub const N: usize = 256;
