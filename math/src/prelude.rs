pub use crate::{fe, fe_array, fe_vec, poly, poly_vec};
pub use crate::{
    field_element::FieldElement,
    ntt::{intt, ntt},
    poly_vector::PolynomialVector,
    polynomial::Polynomial,
};

pub const Q: i32 = 8380417; // Dilithium prime modulus
pub const N: usize = 256;
