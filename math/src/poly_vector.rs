use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Sub, SubAssign};

use crate::{
    error::{PolynomialVectorError, Result},
    poly::Polynomial,
    traits::FiniteField,
};

use num_traits::Zero;

// Macro for convenient polynomial vector creation
///
/// Creates a `PolynomialVector` with the given polynomials. This macro provides
/// several convenient ways to construct polynomial vectors for use in Dilithium
/// and threshold signature operations.
///
/// # Examples
///
/// ## Basic Usage
///
/// Creating a vector from existing polynomials:
/// ```
/// use math::prelude::*;
///
/// let p1: Polynomial<'static, FieldElement> = poly![1, 2, 3];
/// let p2: Polynomial<'static, FieldElement>  = poly![4, 5, 6];
/// let vec: PolynomialVector<'static, FieldElement> = poly_vec!(p1.clone(), p2.clone());
///
/// assert_eq!(vec.len(), 2);
/// assert_eq!(vec.get(0), Some(&p1));
/// assert_eq!(vec.get(1), Some(&p2));
/// ```
///
/// ## Creating Zero Vector
///
/// Creating a zero polynomial vector of specified length:
/// ```
/// use math::prelude::*;
/// use num_traits::identities::Zero;
///
/// // Create zero vector with 3 polynomials
/// let zero_vec: PolynomialVector<'static, FieldElement> = poly_vec![0; 3];
/// assert_eq!(zero_vec.len(), 3);
/// assert!(zero_vec.get(0).unwrap().is_zero());
/// assert!(zero_vec.get(1).unwrap().is_zero());
/// assert!(zero_vec.get(2).unwrap().is_zero());
///
/// // Empty vector
/// let empty: PolynomialVector<'_, FieldElement> = poly_vec![];
/// assert!(empty.is_empty());
/// ```
///
/// ## From Polynomial Expressions
///
/// Creating vectors with polynomial expressions:
/// ```
/// use math::prelude::*;
/// use num_traits::identities::Zero;
///
/// // Using poly! macro inline
/// let vec1: PolynomialVector<'static, FieldElement> = poly_vec![
///     poly![1, 2, 3],
///     poly![4, 5, 6],
///     poly![7, 8, 9]
/// ];
/// assert_eq!(vec1.len(), 3);
///
/// // With trailing comma
/// let vec2: PolynomialVector<'static, FieldElement> = poly_vec![
///     poly![10, 20],
///     poly![30, 40],
/// ];
/// assert_eq!(vec2.len(), 2);
///
/// // Mixed polynomial types
/// let p1: Polynomial<'static, FieldElement> = poly![1; 5];  // 5 ones
/// let p2: Polynomial<'static, FieldElement> = poly![];      // zero polynomial
/// let vec3: PolynomialVector<'static, FieldElement> = poly_vec![p1, p2, poly![1, 2, 3]];
/// assert_eq!(vec3.len(), 3);
/// ```
///
/// # Examples
///
/// ```
/// use math::prelude::*;
///
/// let pv: PolynomialVector<'static, FieldElement> =
///     poly_vec![poly![1, 2], poly![3, 4]];
/// assert_eq!(pv.len(), 2);
/// assert_eq!(pv[0], poly![1, 2]);
///
/// let zeros: PolynomialVector<'static, FieldElement> = poly_vec![0; 3];
/// assert!(zeros.iter().all(|p| p.coefficients().is_empty()));
/// ```
#[macro_export]
macro_rules! poly_vec {
    // Empty case - empty polynomial vector
    () => {
        $crate::poly_vector::PolynomialVector::from_vec(vec![])
    };

    //Single expression that evaluates to Vec<Polynomial>
    ($vec:expr) => {
        $crate::poly_vector::PolynomialVector::from_vec($vec)
    };

    // Repeated zero polynomial case: poly_vec![0; count]
    (0; $count:expr) => {
        $crate::poly_vector::PolynomialVector::zero($count)
    };

    // Repeated polynomial case: poly_vec![poly; count]
    ($poly:expr; $count:expr) => {{
        let value = $poly;
        let mut polys = Vec::with_capacity($count);
        polys.resize_with($count, || value.clone());
        $crate::poly_vector::PolynomialVector::from_vec(polys)
    }};

    // List of polynomials
    ($($poly:expr),+ $(,)?) => {
        $crate::poly_vector::PolynomialVector::from_vec(vec![$($poly),+])
    };
}

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
    ) -> Result<(), PolynomialVectorError> {
        if index >= self.polys.len() {
            return Err(PolynomialVectorError::IndexOutOfBounds {
                index,
                len: self.polys.len(),
            });
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
        self.polys
            .iter()
            .map(|p| p.norm_infinity())
            .max()
            .unwrap_or(0)
    }
}

impl<FF: FiniteField> PolynomialVector<'static, FF> {
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
        let mut lhs = self;
        lhs += &other;
        lhs
    }
}

impl<FF: FiniteField> Add<&PolynomialVector<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        let mut lhs = self;
        lhs += other;
        lhs
    }
}

impl<FF: FiniteField> AddAssign<&PolynomialVector<'static, FF>>
    for PolynomialVector<'static, FF>
{
    fn add_assign(&mut self, rhs: &Self) {
        expect_same_len(self.polys.len(), rhs.polys.len(), "addition");
        for (lhs_poly, rhs_poly) in self.polys.iter_mut().zip(rhs.polys.iter())
        {
            *lhs_poly += rhs_poly;
        }
    }
}

impl<FF: FiniteField> Sub for PolynomialVector<'static, FF> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        let mut lhs = self;
        lhs -= &other;
        lhs
    }
}

impl<FF: FiniteField> Sub<&PolynomialVector<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        let mut lhs = self;
        lhs -= other;
        lhs
    }
}

impl<FF: FiniteField> SubAssign<&PolynomialVector<'static, FF>>
    for PolynomialVector<'static, FF>
{
    fn sub_assign(&mut self, rhs: &Self) {
        expect_same_len(self.polys.len(), rhs.polys.len(), "subtraction");
        for (lhs_poly, rhs_poly) in self.polys.iter_mut().zip(rhs.polys.iter())
        {
            *lhs_poly -= rhs_poly;
        }
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
            PolynomialVector::from_vec(Vec::new());
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

    #[test]
    fn add_assign_with_borrow_updates_in_place() {
        let mut lhs = PolynomialVector::from_vec(vec![
            poly_from_u32s(&[1, 2]),
            poly_from_u32s(&[3, 4]),
        ]);
        let rhs = PolynomialVector::from_vec(vec![
            poly_from_u32s(&[5, 6]),
            poly_from_u32s(&[7, 8]),
        ]);

        lhs += &rhs;

        assert_eq!(lhs[0], poly_from_u32s(&[6, 8]));
        assert_eq!(lhs[1], poly_from_u32s(&[10, 12]));
    }

    #[test]
    fn sub_assign_with_borrow_updates_in_place() {
        let mut lhs = PolynomialVector::from_vec(vec![
            poly_from_u32s(&[9, 9]),
            poly_from_u32s(&[5, 5]),
        ]);
        let rhs = PolynomialVector::from_vec(vec![
            poly_from_u32s(&[4, 1]),
            poly_from_u32s(&[2, 2]),
        ]);

        lhs -= &rhs;

        assert_eq!(lhs[0], poly_from_u32s(&[5, 8]));
        assert_eq!(lhs[1], poly_from_u32s(&[3, 3]));
    }

    #[test]
    fn add_and_sub_behave_like_componentwise_ops() {
        let a = PolynomialVector::from_vec(vec![
            poly_from_u32s(&[1, 2, 3]),
            poly_from_u32s(&[4, 5]),
        ]);
        let b = PolynomialVector::from_vec(vec![
            poly_from_u32s(&[6, 7, 8]),
            poly_from_u32s(&[9, 10]),
        ]);

        let sum = a.clone() + b.clone();
        assert_eq!(sum[0], poly_from_u32s(&[7, 9, 11]));
        assert_eq!(sum[1], poly_from_u32s(&[13, 15]));

        let diff = b - a;
        assert_eq!(diff[0], poly_from_u32s(&[5, 5, 5]));
        assert_eq!(diff[1], poly_from_u32s(&[5, 5]));
    }

    #[test]
    #[should_panic(expected = "Polynomial vector lengths must match for addition")]
    fn add_assign_panics_on_length_mismatch() {
        let mut lhs = PolynomialVector::from_vec(vec![poly_from_u32s(&[1])]);
        let rhs =
            PolynomialVector::from_vec(vec![poly_from_u32s(&[2]), poly_from_u32s(&[3])]);
        lhs += &rhs;
    }

    #[test]
    fn elementwise_multiplication_matches_manual() {
        let a = PolynomialVector::from_vec(vec![
            poly_from_u32s(&[1, 0]),
            poly_from_u32s(&[2]),
        ]);
        let b = PolynomialVector::from_vec(vec![
            poly_from_u32s(&[3]),
            poly_from_u32s(&[4, 5]),
        ]);

        let product = a * b;
        assert_eq!(product[0], poly_from_u32s(&[3, 0]));
        assert_eq!(product[1], poly_from_u32s(&[8, 10]));
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

impl<FF: FiniteField> Mul<Polynomial<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn mul(self, poly: Polynomial<'static, FF>) -> Self {
        self * &poly
    }
}

impl<FF: FiniteField> Mul<&Polynomial<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn mul(self, poly: &Polynomial<'static, FF>) -> Self {
        Self {
            polys: self.polys.into_iter().map(|p| p * poly).collect(),
        }
    }
}

//TODO add test coeverage
