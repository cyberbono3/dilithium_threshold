//! Helper builders shared by the math primitives to keep macros/lightweight constructors DRY.

use crate::{
    field_element::FieldElement, poly::Polynomial,
    poly_vector::PolynomialVector, traits::FiniteField,
};

/// Build a vector by repeating a value that converts into `FieldElement`.
pub fn repeat_field_element<T>(value: T, count: usize) -> Vec<FieldElement>
where
    T: Into<FieldElement> + Copy,
{
    (0..count).map(|_| value.into()).collect()
}

/// Build a polynomial vector by repeating an existing polynomial.
pub fn repeat_polynomial<FF>(
    poly: Polynomial<'static, FF>,
    count: usize,
) -> PolynomialVector<'static, FF>
where
    FF: FiniteField + Clone,
{
    PolynomialVector::new(vec![poly; count])
}

/// Build an empty polynomial vector with a given capacity, avoiding repeated allocations.
pub fn zeroed_vector<FF>(length: usize) -> PolynomialVector<'static, FF>
where
    FF: FiniteField,
{
    PolynomialVector::zero(length)
}

/// Compute the maximum infinity norm given an iterator of per-polynomial norms.
pub fn max_infinity_norm_from_values(
    norms: impl IntoIterator<Item = u32>,
) -> u32 {
    norms.into_iter().max().unwrap_or(0)
}
