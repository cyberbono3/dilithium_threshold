use std::ops::{Add, Deref, DerefMut, Index, IndexMut, Mul, Sub};

use crate::{error::Result, poly::Polynomial, poly_vec, traits::FiniteField};

use num_traits::Zero;

/// Represents a vector of polynomials in Rq.
/// Used for representing keys and intermediate values in Dilithium.
#[derive(Clone, Debug, PartialEq)]
pub struct PolynomialVector<'coeffs, FF: FiniteField> {
    polys: Vec<Polynomial<'coeffs, FF>>,
}

impl<'coeffs, FF: FiniteField> PolynomialVector<'coeffs, FF> {
    /// Construct a polynomial vector from owned polynomials.
    pub fn from_vec(polynomials: Vec<Polynomial<'coeffs, FF>>) -> Self {
        Self { polys: polynomials }
    }

    /// Borrow the underlying polynomials.
    pub fn as_slice(&self) -> &[Polynomial<'coeffs, FF>] {
        &self.polys
    }

    /// Borrow a polynomial at the given index.
    pub fn get(&self, index: usize) -> Option<&Polynomial<'coeffs, FF>> {
        self.polys.get(index)
    }

    /// Mutably borrow a polynomial at the given index.
    pub fn get_mut(
        &mut self,
        index: usize,
    ) -> Option<&mut Polynomial<'coeffs, FF>> {
        self.polys.get_mut(index)
    }

    /// Replace the polynomial at `index`.
    pub fn set(
        &mut self,
        index: usize,
        poly: Polynomial<'coeffs, FF>,
    ) -> Result<(), &'static str> {
        if index >= self.polys.len() {
            return Err("Index out of bounds");
        }
        self.polys[index] = poly;
        Ok(())
    }

    /// Length of the vector.
    pub fn len(&self) -> usize {
        self.polys.len()
    }

    /// Whether the vector contains no polynomials.
    pub fn is_empty(&self) -> bool {
        self.polys.is_empty()
    }

    /// Maximum infinity norm over the contained polynomials.
    pub fn norm_infinity(&self) -> u32 {
        crate::utils::max_infinity_norm_from_values(
            self.polys.iter().map(|p| p.norm_infinity()),
        )
    }
}

impl<FF: FiniteField> PolynomialVector<'static, FF> {
    /// Initialize a polynomial vector with owned polynomials.
    pub fn new(polynomials: Vec<Polynomial<'static, FF>>) -> Self {
        Self::from_vec(polynomials)
    }

    /// Create a zero vector of length `length`.
    pub fn zero(length: usize) -> Self {
        let zero_poly = Polynomial::<FF>::zero();
        Self {
            polys: vec![zero_poly; length],
        }
    }
}

impl<'coeffs, FF: FiniteField> IntoIterator for PolynomialVector<'coeffs, FF> {
    type Item = Polynomial<'coeffs, FF>;
    type IntoIter = std::vec::IntoIter<Polynomial<'coeffs, FF>>;

    fn into_iter(self) -> Self::IntoIter {
        self.polys.into_iter()
    }
}

impl<'a, 'coeffs, FF: FiniteField> IntoIterator
    for &'a PolynomialVector<'coeffs, FF>
{
    type Item = &'a Polynomial<'coeffs, FF>;
    type IntoIter = std::slice::Iter<'a, Polynomial<'coeffs, FF>>;

    fn into_iter(self) -> Self::IntoIter {
        self.polys.iter()
    }
}

impl<'coeffs, FF: FiniteField> From<Vec<Polynomial<'coeffs, FF>>>
    for PolynomialVector<'coeffs, FF>
{
    fn from(polys: Vec<Polynomial<'coeffs, FF>>) -> Self {
        Self::from_vec(polys)
    }
}

impl<'coeffs, FF: FiniteField> FromIterator<Polynomial<'coeffs, FF>>
    for PolynomialVector<'coeffs, FF>
{
    fn from_iter<T: IntoIterator<Item = Polynomial<'coeffs, FF>>>(
        iter: T,
    ) -> Self {
        Self::from_vec(iter.into_iter().collect())
    }
}

impl<'coeffs, FF: FiniteField> Extend<Polynomial<'coeffs, FF>>
    for PolynomialVector<'coeffs, FF>
{
    fn extend<T: IntoIterator<Item = Polynomial<'coeffs, FF>>>(
        &mut self,
        iter: T,
    ) {
        self.polys.extend(iter);
    }
}

impl<'coeffs, FF: FiniteField> Deref for PolynomialVector<'coeffs, FF> {
    type Target = [Polynomial<'coeffs, FF>];
    fn deref(&self) -> &Self::Target {
        &self.polys
    }
}

impl<'coeffs, FF: FiniteField> DerefMut for PolynomialVector<'coeffs, FF> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.polys
    }
}

impl<'coeffs, FF: FiniteField> Index<usize> for PolynomialVector<'coeffs, FF> {
    type Output = Polynomial<'coeffs, FF>;
    fn index(&self, i: usize) -> &Self::Output {
        &self.polys[i]
    }
}

impl<'coeffs, FF: FiniteField> IndexMut<usize>
    for PolynomialVector<'coeffs, FF>
{
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.polys[i]
    }
}

impl<FF: FiniteField> Add for PolynomialVector<'static, FF> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        let polys = combine_owned_vectors(
            self.polys,
            other.polys,
            "addition",
            |a, b| a + b,
        );
        Self { polys }
    }
}

impl<FF: FiniteField> Add<&PolynomialVector<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        let polys = combine_owned_with_slice(
            self.polys,
            &other.polys,
            "addition",
            |a, b| a + b.clone(),
        );

        Self { polys }
    }
}

impl<FF: FiniteField> Sub for PolynomialVector<'static, FF> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        let polys = combine_owned_vectors(
            self.polys,
            other.polys,
            "subtraction",
            |a, b| a - b,
        );

        Self { polys }
    }
}

impl<FF: FiniteField> Sub<&PolynomialVector<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        let polys = combine_owned_with_slice(
            self.polys,
            &other.polys,
            "subtraction",
            |a, b| a - b.clone(),
        );

        Self { polys }
    }
}
// TODO add proper testing
impl<FF: FiniteField> Mul for PolynomialVector<'static, FF> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        let polys = combine_owned_vectors(
            self.polys,
            other.polys,
            "multiplication",
            |a, b| a * b,
        );

        Self { polys }
    }
}

fn expect_same_len(lhs: usize, rhs: usize, context: &str) {
    assert!(
        lhs == rhs,
        "Polynomial vector lengths must match for {context}: {lhs} vs {rhs}"
    );
}

fn combine_owned_vectors<FF, F>(
    lhs: Vec<Polynomial<'static, FF>>,
    rhs: Vec<Polynomial<'static, FF>>,
    context: &str,
    mut op: F,
) -> Vec<Polynomial<'static, FF>>
where
    FF: FiniteField,
    F: FnMut(
        Polynomial<'static, FF>,
        Polynomial<'static, FF>,
    ) -> Polynomial<'static, FF>,
{
    expect_same_len(lhs.len(), rhs.len(), context);
    lhs.into_iter().zip(rhs).map(|(a, b)| op(a, b)).collect()
}

fn combine_owned_with_slice<FF, F>(
    lhs: Vec<Polynomial<'static, FF>>,
    rhs: &[Polynomial<'static, FF>],
    context: &str,
    mut op: F,
) -> Vec<Polynomial<'static, FF>>
where
    FF: FiniteField,
    F: FnMut(
        Polynomial<'static, FF>,
        &Polynomial<'static, FF>,
    ) -> Polynomial<'static, FF>,
{
    expect_same_len(lhs.len(), rhs.len(), context);
    lhs.into_iter()
        .zip(rhs.iter())
        .map(|(a, b)| op(a, b))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    fn poly_from_u32s(coeffs: &[u32]) -> Polynomial<'static, FieldElement> {
        let elems = coeffs
            .iter()
            .copied()
            .map(FieldElement::from)
            .collect::<Vec<_>>();
        Polynomial::from(elems)
    }

    #[test]
    fn collect_into_polynomial_vector() {
        let polys = vec![poly_from_u32s(&[1]), poly_from_u32s(&[2])];
        let vec: PolynomialVector<'static, FieldElement> =
            polys.clone().into_iter().collect();
        assert_eq!(vec.len(), polys.len());
        assert_eq!(vec.as_slice(), polys.as_slice());
    }

    #[test]
    fn extend_appends_polynomials() {
        let mut vec = PolynomialVector::from_vec(vec![
            poly_from_u32s(&[1]),
            poly_from_u32s(&[2]),
        ]);
        vec.extend(vec![poly_from_u32s(&[3]), poly_from_u32s(&[4])]);
        assert_eq!(vec.len(), 4);
        assert_eq!(vec[0], poly_from_u32s(&[1]));
        assert_eq!(vec[3], poly_from_u32s(&[4]));
    }

    #[test]
    fn into_iter_consumes_vector() {
        let vec = PolynomialVector::from_vec(vec![
            poly_from_u32s(&[5]),
            poly_from_u32s(&[6]),
        ]);
        let collected: Vec<_> = vec.into_iter().collect();
        assert_eq!(collected, vec![poly_from_u32s(&[5]), poly_from_u32s(&[6])]);
    }

    #[test]
    fn norm_infinity_empty_is_zero() {
        let vec: PolynomialVector<'static, FieldElement> =
            PolynomialVector::new(Vec::new());
        assert_eq!(0, vec.norm_infinity());
    }

    #[test]
    fn norm_infinity_reports_max_coeff_norm() {
        let polys = vec![
            poly_from_u32s(&[1, 2, 3]),    // max |coeff| = 3
            poly_from_u32s(&[10, 20]),     // max = 20
            poly_from_u32s(&[4, 5, 6, 7]), // max = 7
        ];
        let vec = PolynomialVector::from_vec(polys);
        assert_eq!(20, vec.norm_infinity());
    }
}

impl<FF: FiniteField> Mul<u64> for PolynomialVector<'static, FF> {
    type Output = PolynomialVector<'static, FF>;

    fn mul(self, scalar: u64) -> Self::Output {
        let scalar_ff = FF::from(scalar);

        let polys = self
            .polys
            .iter()
            .map(|p| {
                let scaled_coeffs: Vec<FF> =
                    p.coefficients().iter().map(|&c| c * scalar_ff).collect();
                Polynomial::new(scaled_coeffs)
            })
            .collect();

        PolynomialVector { polys }
    }
}

// impl<FF: FiniteField> Mul<PolynomialVector<'static, FF>>
//     for &Vec<Vec<Polynomial<'static, FF>>>
// {
//     type Output = PolynomialVector<'static, FF>;

//     fn mul(self, vector: PolynomialVector<'static, FF>) -> Self::Output {
//         matrix_vector_multiply(self, &vector)
//     }
// }

// impl<FF: FiniteField> Mul<&PolynomialVector<'static, FF>>
//     for &Vec<Vec<Polynomial<'static, FF>>>
// {
//     type Output = PolynomialVector<'static, FF>;

//     fn mul(self, vector: &PolynomialVector<'static, FF>) -> Self::Output {
//         matrix_vector_multiply(self, vector)
//     }
// }

pub fn matrix_vector_multiply<FF: FiniteField>(
    m: &[Vec<Polynomial<'static, FF>>],
    v: &PolynomialVector<'static, FF>,
) -> PolynomialVector<'static, FF> {
    // Check matrix is not empty
    assert!(!m.is_empty(), "Matrix cannot be empty");

    let cols = m[0].len();

    // Check matrix columns match vector length
    assert_eq!(
        cols,
        v.len(),
        "Matrix columns ({}) must match vector length ({})",
        cols,
        v.len()
    );

    // Check all rows have same length (rectangular matrix)
    assert!(
        m.iter().all(|row| row.len() == cols),
        "All matrix rows must have the same length"
    );

    let result = m
        .iter()
        .map(|row| {
            row.iter().zip(v.as_slice()).fold(
                Polynomial::<FF>::zero(),
                |mut acc, (a_ij, v_j)| {
                    // TODO make &Polynomial be multiplied &POlynomial, remove .clone()
                    acc += a_ij.clone() * v_j.clone();
                    acc
                },
            )
        })
        .collect();

    poly_vec!(result)
}

impl<FF: FiniteField> Mul<Polynomial<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn mul(self, poly: Polynomial<'static, FF>) -> Self {
        Self {
            polys: self.polys.into_iter().map(|p| p * poly.clone()).collect(),
        }
    }
}

impl<FF: FiniteField> Mul<&Polynomial<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn mul(self, poly: &Polynomial<'static, FF>) -> Self {
        Self {
            polys: self.polys.into_iter().map(|p| p * poly.clone()).collect(),
        }
    }
}

//TODO add test coeverage
