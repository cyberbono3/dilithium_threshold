
use std::ops::{Add, Sub, Mul, Index, IndexMut, Deref, DerefMut};

use crate::{
    polynomial::Polynomial,
    poly_vector::PolynomialVector,
    traits::FiniteField,
    error::{Result, Error},
};

use num_traits::Zero;

/// A simple, rectangular matrix of polynomials over a finite field.
///
/// This mirrors the style and ergonomics of `PolynomialVector`:
/// - 'static lifetime for owned coefficient storage in mutating/arith ops
/// - fallible (`try_*`) and panicking variants for shape-checked operations
/// - standard operator trait impls where ergonomic
#[derive(Clone, Debug, PartialEq)]
pub struct Matrix<'coeffs, FF: FiniteField> {
    rows: Vec<Vec<Polynomial<'coeffs, FF>>>,
}

impl<FF: FiniteField> Matrix<'static, FF> {
    /// Construct a new matrix from rows. Panics if rows have differing lengths.
    pub fn new(rows: Vec<Vec<Polynomial<'static, FF>>>) -> Self {
        if !rows.is_empty() {
            let cols = rows[0].len();
            assert!(
                rows.iter().all(|r| r.len() == cols),
                "All matrix rows must have the same length"
            );
        }
        Self { rows }
    }

    /// Borrow the underlying rows.
    pub fn as_slice(&self) -> &[Vec<Polynomial<'static, FF>>] { &self.rows }

    /// Number of rows.
    pub fn rows(&self) -> usize { self.rows.len() }

    /// Number of columns (0 if empty).
    pub fn cols(&self) -> usize { if self.rows.is_empty() { 0 } else { self.rows[0].len() } }

    /// (rows, cols)
    pub fn shape(&self) -> (usize, usize) { (self.rows(), self.cols()) }

    /// Create a zero matrix of size (rows x cols).
    pub fn zeros(rows: usize, cols: usize) -> Self {
        let mut data = Vec::with_capacity(rows);
        for _ in 0..rows {
            let mut r = Vec::with_capacity(cols);
            for _ in 0..cols {
                r.push(Polynomial::<FF>::zero());
            }
            data.push(r);
        }
        Self { rows: data }
    }

    /// Infinity norm over entries, using the polynomial infinity norm.
    pub fn norm_infinity(&self) -> u32 {
        self.rows
            .iter()
            .flat_map(|r| r.iter())
            .map(|p| p.norm_infinity())
            .max()
            .unwrap_or(0)
    }

    /// Fallible matrix–vector multiplication with shape checks.
    ///
    /// Returns an error on any shape mismatch. For consistency with the existing
    /// `PolynomialVector` utilities, uses the same error variant currently used
    /// in `poly_vector.rs` for shape issues.
    pub fn try_mul_vector(
        &self,
        v: &PolynomialVector<'static, FF>,
    ) -> Result<PolynomialVector<'static, FF>> {
        // mirror the style from `poly_vector.rs`
        if self.rows.is_empty() {
            return Err(Error::Polynomial(super::error::PolynomialError::CoefficientOutOfRange));
        }
        let cols = self.cols();
        if cols != v.len() {
            return Err(Error::Polynomial(super::error::PolynomialError::CoefficientOutOfRange));
        }
        if !self.rows.iter().all(|row| row.len() == cols) {
            return Err(Error::Polynomial(super::error::PolynomialError::CoefficientOutOfRange));
        }
        Ok(self.mul_vector(v))
    }

    /// Panicking matrix–vector multiplication (asserts on shape mismatch).
    pub fn mul_vector(&self, v: &PolynomialVector<'static, FF>) -> PolynomialVector<'static, FF> {
        // Check matrix is not empty
        assert!(!self.rows.is_empty(), "Matrix cannot be empty");

        let cols = self.cols();
        // Check matrix columns match vector length
        assert_eq!(
            cols, v.len(),
            "Matrix columns ({}) must match vector length ({})",
            cols, v.len()
        );

        // Check all rows have same length (rectangular matrix)
        assert!(
            self.rows.iter().all(|row| row.len() == cols),
            "All matrix rows must have the same length"
        );

        let polys = self
            .rows
            .iter()
            .map(|row| {
                row.iter().zip(v.as_slice()).fold(
                    Polynomial::<FF>::zero(),
                    |mut acc, (a_ij, v_j)| {
                        // TODO: enable &Polynomial * &Polynomial to avoid clones
                        acc += a_ij.clone() * v_j.clone();
                        acc
                    },
                )
            })
            .collect();

        PolynomialVector::new(polys)
    }
}

/// Immutable indexing by row.
impl<FF: FiniteField> Index<usize> for Matrix<'static, FF> {
    type Output = Vec<Polynomial<'static, FF>>;
    fn index(&self, i: usize) -> &Self::Output { &self.rows[i] }
}

/// Mutable indexing by row.
impl<FF: FiniteField> IndexMut<usize> for Matrix<'static, FF> {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output { &mut self.rows[i] }
}

/// Deref to a slice of rows.
impl<FF: FiniteField> Deref for Matrix<'static, FF> {
    type Target = [Vec<Polynomial<'static, FF>>];
    fn deref(&self) -> &Self::Target { &self.rows }
}

/// DerefMut to a slice of rows.
impl<FF: FiniteField> DerefMut for Matrix<'static, FF> {
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.rows }
}

/// Element-wise matrix addition. Panics if shapes are not equal.
impl<FF: FiniteField> Add for Matrix<'static, FF> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let (r1, c1) = self.shape();
        let (r2, c2) = rhs.shape();
        assert!(
            r1 == r2 && c1 == c2,
            "Matrix shapes must match for addition: ({}, {}) vs ({}, {})",
            r1, c1, r2, c2
        );

        let rows = self
            .rows
            .into_iter()
            .zip(rhs.rows)
            .map(|(row_a, row_b)| {
                row_a.into_iter()
                     .zip(row_b.into_iter())
                     .map(|(a, b)| a + b)
                     .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        Self { rows }
    }
}

/// Element-wise matrix subtraction. Panics if shapes are not equal.
impl<FF: FiniteField> Sub for Matrix<'static, FF> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let (r1, c1) = self.shape();
        let (r2, c2) = rhs.shape();
        assert!(
            r1 == r2 && c1 == c2,
            "Matrix shapes  must match for subtraction: ({}, {}) vs ({}, {})",
            r1, c1, r2, c2
        );

        let rows = self
            .rows
            .into_iter()
            .zip(rhs.rows)
            .map(|(row_a, row_b)| {
                row_a.into_iter()
                     .zip(row_b.into_iter())
                     .map(|(a, b)| a - b)
                     .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        Self { rows }
    }
}

/// Scalar (polynomial) multiplication: each entry multiplied by the same polynomial.
impl<FF: FiniteField> Mul<&Polynomial<'static, FF>> for Matrix<'static, FF> {
    type Output = Self;

    fn mul(self, scalar: &Polynomial<'static, FF>) -> Self::Output {
        let rows = self.rows
            .into_iter()
            .map(|row| {
                row.into_iter()
                   .map(|p| p * scalar.clone())
                   .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
            Self { rows }
    }
}

/// Matrix–vector multiplication via `*` operator.
impl<FF: FiniteField> Mul<PolynomialVector<'static, FF>> for Matrix<'static, FF> {
    type Output = PolynomialVector<'static, FF>;

    fn mul(self, v: PolynomialVector<'static, FF>) -> Self::Output {
        self.mul_vector(&v)
    }
}

/// Also support `&Matrix * &PolynomialVector` ergonomically.
impl<FF: FiniteField> Mul<&PolynomialVector<'static, FF>> for &Matrix<'static, FF> {
    type Output = PolynomialVector<'static, FF>;

    fn mul(self, v: &PolynomialVector<'static, FF>) -> Self::Output {
        self.mul_vector(v)
    }
}

/// Convenience free functions with the same naming pattern used in `poly_vector.rs`.
pub fn try_matrix_vector_multiply<FF: FiniteField>(
    m: &Matrix<'static, FF>,
    v: &PolynomialVector<'static, FF>,
) -> Result<PolynomialVector<'static, FF>> {
    m.try_mul_vector(v)
}

pub fn matrix_vector_multiply<FF: FiniteField>(
    m: &Matrix<'static, FF>,
    v: &PolynomialVector<'static, FF>,
) -> PolynomialVector<'static, FF> {
    m.mul_vector(v)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        error::{Error as MathError, PolynomialError},
        field_element::FieldElement,
        poly_vector::PolynomialVector,
        polynomial::Polynomial,
    };
    use num_traits::Zero;

    type FE   = FieldElement;
    type Poly = Polynomial<'static, FE>;
    type Mat  = Matrix<'static, FE>;
    type PV   = PolynomialVector<'static, FE>;

    #[inline]
    fn c(v: u32) -> Poly { Polynomial::from(vec![v]) }

    #[inline]
    fn pv_from_consts(vals: &[u32]) -> PV {
        PolynomialVector::new(vals.iter().cloned().map(c).collect())
    }

    // --------------------------
    // shape / construction
    // --------------------------

    #[test]
    fn zeros_shape_and_entries_are_zero() {
        let m: Mat = Mat::zeros(2, 3);
        assert_eq!(m.rows(), 2);
        assert_eq!(m.cols(), 3);
        assert_eq!(m.shape(), (2, 3));
        for i in 0..2 {
            for j in 0..3 {
                assert!(m[i][j].is_zero(), "entry ({},{}) not zero", i, j);
            }
        }
    }

    #[test]
    #[should_panic(expected = "All matrix rows must have the same length")]
    fn new_ragged_panics() {
        // First row has 2 entries, second row has 1 -> should panic.
        let _ = Mat::new(vec![
            vec![c(1), c(2)],
            vec![c(3)],
        ]);
    }

    #[test]
    fn shape_rows_cols_work() {
        let m = Mat::new(vec![
            vec![c(1), Polynomial::from(vec![2u32, 3u32])],
            vec![Polynomial::from(vec![4u32, 5u32, 6u32]), c(7)],
        ]);
        assert_eq!(m.rows(), 2);
        assert_eq!(m.cols(), 2);
        assert_eq!(m.shape(), (2, 2));

        // as_slice is a thin borrow of underlying rows
        let rows_slice: &[Vec<Poly>] = m.as_slice();
        assert_eq!(rows_slice.len(), m.rows());
    }

    // --------------------------
    // indexing / mutation
    // --------------------------

    #[test]
    fn indexing_and_mutation() {
        let mut m = Mat::zeros(2, 2);
        assert!(m[1][0].is_zero());
        m[1][0] = Polynomial::from(vec![5u32, 6u32]);
        assert_eq!(m[1][0], Polynomial::from(vec![5u32, 6u32]));
    }

    // --------------------------
    // norms
    // --------------------------

    #[test]
    fn norm_infinity_empty_is_zero() {
        let m: Mat = Mat::new(vec![]);
        assert_eq!(m.norm_infinity(), 0);
    }

    #[test]
    fn norm_infinity_nonempty() {
        // Entries with max centered abs coeffs: [1,2,3] -> 3, [1] -> 1, [4,10] -> 10, [3] -> 3
        let m = Mat::new(vec![
            vec![Polynomial::from(vec![1u32, 2u32, 3u32]), c(1)],
            vec![Polynomial::from(vec![4u32, 10u32]), c(3)],
        ]);
        assert_eq!(m.norm_infinity(), 10);
    }

    // --------------------------
    // element-wise + and -
    // --------------------------

    #[test]
    fn elementwise_add() {
        let a = Mat::new(vec![
            vec![c(1), c(2)],
            vec![c(3), c(4)],
        ]);
        let b = Mat::new(vec![
            vec![c(10), c(20)],
            vec![c(30), c(40)],
        ]);
        let sum = a + b;
        assert_eq!(sum[0][0], c(11));
        assert_eq!(sum[0][1], c(22));
        assert_eq!(sum[1][0], c(33));
        assert_eq!(sum[1][1], c(44));
    }

    #[test]
    #[should_panic(expected = "Matrix shapes must match for addition")]
    fn elementwise_add_mismatched_shapes_panics() {
        let a = Mat::new(vec![vec![c(1), c(2)]]);
        let b = Mat::new(vec![vec![c(3)], vec![c(4)]]);
        let _ = a + b;
    }

    #[test]
    fn elementwise_sub() {
        let a = Mat::new(vec![
            vec![Polynomial::from(vec![10u32]), Polynomial::from(vec![20u32, 1u32])],
            vec![c(30), c(40)],
        ]);
        let b = Mat::new(vec![
            vec![c(7), c(5)],
            vec![c(12), c(40)],
        ]);
        let d = a - b;
        assert_eq!(d[0][0], c(3));
        //  (20 + x) - 5 = (15 + x)
        assert_eq!(d[0][1], Polynomial::from(vec![15u32, 1u32]));
        assert_eq!(d[1][0], c(18));
        assert_eq!(d[1][1], c(0));
    }

    #[test]
    #[should_panic(expected = "Matrix shapes  must match for subtraction")]
    fn elementwise_sub_mismatched_shapes_panics() {
        let a = Mat::new(vec![vec![c(1)], vec![c(2)]]);
        let b = Mat::new(vec![vec![c(1), c(2)]]);
        let _ = a - b;
    }

    // --------------------------
    // scalar (polynomial) multiplication
    // --------------------------

    #[test]
    fn scalar_multiply_by_constant_poly() {
        let m = Mat::new(vec![
            vec![Polynomial::from(vec![1u32, 1u32]), c(3)],
            vec![c(4), c(5)],
        ]);
        let s = c(2);
        let out = m.clone() * &s;
        assert_eq!(out[0][0], Polynomial::from(vec![2u32, 2u32]));
        assert_eq!(out[0][1], c(6));
        assert_eq!(out[1][0], c(8));
        assert_eq!(out[1][1], c(10));

        // property: out[i][j] == m[i][j] * s
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(out[i][j], m[i][j].clone() * s.clone());
            }
        }
    }

    // --------------------------
    // matrix–vector multiplication
    // --------------------------

    #[test]
    fn mul_vector_with_constants() {
        // 2x3 * 3 -> 2
        let m = Mat::new(vec![
            vec![c(1), c(2), c(3)],
            vec![c(4), c(5), c(6)],
        ]);
        let v = pv_from_consts(&[7, 8, 9]);

        let out = m.mul_vector(&v);
        assert_eq!(out.len(), 2);
        // Row 0: 1*7 + 2*8 + 3*9 = 50
        assert_eq!(out[0], c(50));
        // Row 1: 4*7 + 5*8 + 6*9 = 122
        assert_eq!(out[1], c(122));

        // operator form should agree
        let out2 = m.clone() * v.clone();
        assert_eq!(out2, out);
    }

    #[test]
    fn mul_vector_with_nonconstant_polynomials() {
        // Let x = [0,1], then:
        // Row 0: [x, 1, 0] • [1, 2, 3] = x*1 + 1*2 + 0*3 = x + 2
        // Row 1: [0, x, 3] • [1, 2, 3] = 0*1 + x*2 + 3*3 = 2x + 9
        let x   = Polynomial::x_to_the(1);
        let one = c(1);
        let two = c(2);
        let three = c(3);

        let m = Mat::new(vec![
            vec![x.clone(), one.clone(), Polynomial::zero()],
            vec![Polynomial::zero(), x.clone(), three.clone()],
        ]);
        let v = PV::new(vec![one.clone(), two.clone(), three.clone()]);

        let out = m.mul_vector(&v);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0], x.clone() + two.clone());              // x + 2
        assert_eq!(out[1], (x * two.clone()) + c(9));             // 2x + 9
    }

    #[test]
    #[should_panic(expected = "Matrix cannot be empty")]
    fn mul_vector_panics_on_empty_matrix() {
        let m: Mat = Mat::new(vec![]);
        let v = PV::new(vec![]);
        let _ = m.mul_vector(&v);
    }

    #[test]
    #[should_panic(expected = "Matrix columns (3) must match vector length (2)")]
    fn mul_vector_panics_on_mismatched_shapes() {
        let m = Mat::new(vec![
            vec![c(1), c(2), c(3)],
            vec![c(4), c(5), c(6)],
        ]);
        let v = pv_from_consts(&[7, 8]); // length 2, but matrix has 3 columns
        let _ = m.mul_vector(&v);
    }

    // --------------------------
    // try_mul_vector (fallible)
    // --------------------------

    #[test]
    fn try_mul_vector_errors_on_empty_or_mismatch() {
        // empty matrix
        let m_empty: Mat = Mat::new(vec![]);
        let v_empty = PV::new(vec![]);
        let err = m_empty.try_mul_vector(&v_empty).unwrap_err();
        assert_eq!(err, MathError::Polynomial(PolynomialError::CoefficientOutOfRange));

        // column/length mismatch
        let m = Mat::new(vec![
            vec![c(1), c(2), c(3)],
            vec![c(4), c(5), c(6)],
        ]);
        let v_bad = pv_from_consts(&[7, 8]);
        let err2 = m.try_mul_vector(&v_bad).unwrap_err();
        assert_eq!(err2, MathError::Polynomial(PolynomialError::CoefficientOutOfRange));
    }

    // --------------------------
    // wrapper free functions
    // --------------------------

    #[test]
    fn wrapper_free_functions_match_methods() {
        let m = Mat::new(vec![
            vec![c(1), c(2)],
            vec![c(3), c(4)],
        ]);
        let v = pv_from_consts(&[5, 6]);

        let out_method = m.mul_vector(&v);
        let out_free = super::matrix_vector_multiply(&m, &v);
        assert_eq!(out_method, out_free);

        // fallible wrapper
        let out_ok = super::try_matrix_vector_multiply(&m, &v).unwrap();
        assert_eq!(out_ok, out_method);
    }
}

