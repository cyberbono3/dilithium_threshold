use crate::hash::shake128;
use crate::params::{K, L, N};
use crate::utils::zero_polyvec;
use math::{poly::Polynomial, traits::FiniteField};
use num_traits::Zero;
use std::ops::Mul;

#[derive(Clone, Debug, PartialEq)]
pub struct MatrixA<'a, FF: FiniteField> {
    pub rows: Vec<Vec<Polynomial<'a, FF>>>, // K x L
}

impl<FF: FiniteField> MatrixA<'static, FF> {
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
    pub fn as_slice(&self) -> &[Vec<Polynomial<'static, FF>>] {
        &self.rows
    }

    /// Number of rows.
    pub fn rows(&self) -> usize {
        self.rows.len()
    }

    /// Number of columns (0 if empty).
    pub fn cols(&self) -> usize {
        if self.rows.is_empty() {
            0
        } else {
            self.rows[0].len()
        }
    }

    /// (rows, cols)
    pub fn shape(&self) -> (usize, usize) {
        (self.rows(), self.cols())
    }

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
    /// // TODO add error handling
    // pub fn try_mul_vector(
    //     &self,
    //     v: &[Polynomial<'static, FF>],
    // ) -> Vec<Polynomial<'static, FF>> {
    //     // mirror the style from `poly_vector.rs`
    //     // if self.rows.is_empty() {
    //     //     return Err(Error::Polynomial(
    //     //         PolynomialError::CoefficientOutOfRange,
    //     //     ));
    //     // }
    //     assert!(!self.rows.is_empty(), "No rows in matrix");
    //     let cols = self.cols();
    //     assert_eq!(cols, v.len(), "Number of columns must match the length of vector");
    //     // let cols = self.cols();
    //     // if cols != v.len() {
    //     //     return Err(Error::Polynomial(
    //     //         PolynomialError::CoefficientOutOfRange,
    //     //     ));
    //     // }
    //     assert!(!self.rows.iter().all(|row| row.len() == cols), "Length of every row  must match number of columns");
    //     // if !self.rows.iter().all(|row| row.len() == cols) {
    //     //     return Err(Error::Polynomial(
    //     //         PolynomialError::CoefficientOutOfRange,
    //     //     ));
    //     // }
    //     //Ok(self.mul_vector(v))
    //     self.mul_vector(v)
    // }

    // TODO add proper error handling and testing
    /// Panicking matrix–vector multiplication (asserts on shape mismatch).
    pub fn mul_vector(
        &self,
        v: &[Polynomial<'static, FF>],
    ) -> Vec<Polynomial<'static, FF>> {
        // Check matrix is not empty
        assert!(!self.rows.is_empty(), "Matrix cannot be empty");

        let cols = self.cols();
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
            self.rows.iter().all(|row| row.len() == cols),
            "All matrix rows must have the same length"
        );

        self.rows.iter().map(|row| row_mul(row, v)).collect()
    }
}

impl<'a, 'b, FF: FiniteField + 'static> Mul<&[Polynomial<'b, FF>; L]>
    for &MatrixA<'a, FF>
{
    type Output = [Polynomial<'static, FF>; K];

    fn mul(self, rhs: &[Polynomial<'b, FF>; L]) -> Self::Output {
        let mut result = zero_polyvec::<K, FF>();
        for (slot, row) in result.iter_mut().zip(&self.rows) {
            *slot = row_mul(row, rhs);
        }
        result
    }
}

fn row_mul<'a, 'b, FF: FiniteField>(
    row: &[Polynomial<'a, FF>],
    vec: &[Polynomial<'b, FF>],
) -> Polynomial<'static, FF> {
    row.iter().zip(vec.iter()).fold(
        Polynomial::<FF>::zero(),
        |mut acc, (a, b)| {
            acc += a.clone() * b.clone();
            acc
        },
    )
}

/// Expand A from rho using SHAKE128 as XOF (educational: uses modulo reduction).
pub fn expand_a_from_rho<FF: FiniteField + std::convert::From<i64>>(
    rho: [u8; 32],
) -> MatrixA<'static, FF> {
    let mut rows = Vec::with_capacity(K);
    for i in 0..K {
        let mut row = Vec::with_capacity(L);
        for j in 0..L {
            // Domain sep: rho || i || j
            let mut seed = Vec::with_capacity(32 + 4);
            seed.extend_from_slice(&rho);
            seed.extend_from_slice(&(i as u16).to_le_bytes());
            seed.extend_from_slice(&(j as u16).to_le_bytes());
            let stream = shake128(4 * N, &seed); // 4*N bytes -> N u32s
            let mut coeffs = [0u32; N];
            for t in 0..N {
                let b = &stream[4 * t..4 * t + 4];
                let v = u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
                coeffs[t] = v;
            }
            row.push(Polynomial::from(coeffs));
        }
        rows.push(row);
    }
    MatrixA { rows }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{K, L, Q};
    use crate::utils::zero_polyvec;
    use math::field_element::FieldElement;
    use num_traits::Zero;

    fn poly(coeffs: &[i64]) -> Polynomial<'static, FieldElement> {
        Polynomial::from(coeffs.to_vec())
    }

    #[test]
    fn expand_a_is_deterministic_and_modq() {
        let rho = [7u8; 32];
        let a1 = expand_a_from_rho::<FieldElement>(rho);
        let a2 = expand_a_from_rho::<FieldElement>(rho);
        assert_eq!(a1.rows.len(), K);
        assert_eq!(a2.rows.len(), K);
        for i in 0..K {
            assert_eq!(a1.rows[i].len(), L);
            assert_eq!(a2.rows[i].len(), L);
            for j in 0..L {
                // determinism
                assert_eq!(a1.rows[i][j], a2.rows[i][j]);
                // coeffs in [0, Q)
                for &c in a1.rows[i][j].coefficients() {
                    let num: i64 = c.into();
                    assert!((0..Q).contains(&num), "coefficient out of range");
                }
            }
        }
    }

    #[test]
    fn matrix_constructors_and_accessors() {
        let row = vec![Polynomial::<FieldElement>::zero(); L];
        let rows = vec![row.clone(); K];
        let matrix = MatrixA::new(rows.clone());

        assert_eq!(matrix.shape(), (K, L));
        assert_eq!(matrix.rows(), K);
        assert_eq!(matrix.cols(), L);
        assert_eq!(matrix.as_slice(), &rows);
    }

    #[test]
    fn matrix_zeros_produces_zero_matrix() {
        let matrix = MatrixA::<FieldElement>::zeros(3, 2);
        assert_eq!(matrix.shape(), (3, 2));
        assert!(
            matrix
                .as_slice()
                .iter()
                .all(|row| row.iter().all(|poly| poly.is_zero()))
        );
    }

    #[test]
    fn norm_infinity_returns_max_entry_norm() {
        let mut matrix = MatrixA::<FieldElement>::zeros(2, 2);
        matrix.rows[0][0] = Polynomial::from([FieldElement::from(5i64)]);
        matrix.rows[1][1] = Polynomial::from([FieldElement::from(10i64)]);

        assert_eq!(matrix.norm_infinity(), 10);
    }

    #[test]
    fn mul_vector_matches_row_mul() {
        let a = expand_a_from_rho::<FieldElement>([11u8; 32]);
        let y = zero_polyvec::<L, FieldElement>();

        let via_method = a.mul_vector(&y);
        let via_operator = (&a) * &y;

        assert_eq!(via_method.as_slice(), via_operator.as_slice());
    }

    #[test]
    #[should_panic(
        expected = "Matrix columns (2) must match vector length (3)"
    )]
    fn mul_vector_panics_on_shape_mismatch() {
        let matrix = MatrixA::<FieldElement>::zeros(2, 2);
        let vector =
            [Polynomial::zero(), Polynomial::zero(), Polynomial::zero()];
        let _ = matrix.mul_vector(&vector);
    }

    #[test]
    fn row_mul_matches_manual_computation() {
        let row = vec![
            Polynomial::from([FieldElement::from(2i64)]),
            Polynomial::from([FieldElement::from(3i64)]),
        ];
        let vec = vec![
            Polynomial::from([FieldElement::from(5i64)]),
            Polynomial::from([FieldElement::from(7i64)]),
        ];

        let result: Polynomial<'static, FieldElement> =
            super::row_mul(&row, &vec);

        let mut expected = Polynomial::from([FieldElement::from(10i64)]);
        expected += Polynomial::from([FieldElement::from(21i64)]);

        assert_eq!(result, expected);
    }

    #[test]
    fn row_mul_with_zero_vector_is_zero() {
        let row = vec![
            Polynomial::from([FieldElement::from(2i64)]),
            Polynomial::from([FieldElement::from(3i64)]),
        ];
        let zero_vec = vec![Polynomial::zero(), Polynomial::zero()];

        let result: Polynomial<'static, FieldElement> =
            super::row_mul(&row, &zero_vec);

        assert!(result.is_zero());
    }

    #[test]
    fn expand_a_changes_with_rho() {
        let a1 = expand_a_from_rho::<FieldElement>([1u8; 32]);
        let a2 = expand_a_from_rho::<FieldElement>([2u8; 32]);
        // It's overwhelmingly likely that at least one position differs.
        let mut any_diff = false;
        'outer: for i in 0..K {
            for j in 0..L {
                if a1.rows[i][j] != a2.rows[i][j] {
                    any_diff = true;
                    break 'outer;
                }
            }
        }
        assert!(
            any_diff,
            "Matrices for different rho should differ (with overwhelming probability)"
        );
    }

    #[test]
    fn mul_operator_zero_is_zero() {
        let a = expand_a_from_rho::<FieldElement>([42u8; 32]);
        let y = zero_polyvec::<L, FieldElement>();
        let w = (&a) * &y;
        for i in 0..K {
            assert_eq!(w[i], Polynomial::zero());
        }
    }

    #[test]
    fn mul_operator_produces_expected_polynomials() {
        let rows = vec![
            vec![poly(&[1]), poly(&[2]), poly(&[3]), poly(&[0])],
            vec![poly(&[0, 1]), poly(&[2]), poly(&[0]), poly(&[0])],
            vec![poly(&[1]), poly(&[1]), poly(&[1]), poly(&[0])],
            vec![poly(&[0]), poly(&[0]), poly(&[0]), poly(&[5])],
        ];
        let matrix = MatrixA::new(rows.clone());

        let rhs = [poly(&[1, 1]), poly(&[2]), poly(&[0, 1]), poly(&[1])];

        let product = (&matrix) * &rhs;

        let mut expected_row0 = rows[0][0].clone() * rhs[0].clone();
        expected_row0 += rows[0][1].clone() * rhs[1].clone();
        expected_row0 += rows[0][2].clone() * rhs[2].clone();
        expected_row0 += rows[0][3].clone() * rhs[3].clone();

        let mut expected_row1 = rows[1][0].clone() * rhs[0].clone();
        expected_row1 += rows[1][1].clone() * rhs[1].clone();
        expected_row1 += rows[1][2].clone() * rhs[2].clone();
        expected_row1 += rows[1][3].clone() * rhs[3].clone();

        let mut expected_row2 = rows[2][0].clone() * rhs[0].clone();
        expected_row2 += rows[2][1].clone() * rhs[1].clone();
        expected_row2 += rows[2][2].clone() * rhs[2].clone();
        expected_row2 += rows[2][3].clone() * rhs[3].clone();

        let mut expected_row3 = rows[3][0].clone() * rhs[0].clone();
        expected_row3 += rows[3][1].clone() * rhs[1].clone();
        expected_row3 += rows[3][2].clone() * rhs[2].clone();
        expected_row3 += rows[3][3].clone() * rhs[3].clone();

        let expected =
            [expected_row0, expected_row1, expected_row2, expected_row3];

        assert_eq!(product, expected);
    }

    #[test]
    fn mul_operator_does_not_mutate_inputs() {
        let rows = vec![
            vec![poly(&[1]), poly(&[2]), poly(&[3]), poly(&[0])],
            vec![poly(&[0, 1]), poly(&[2]), poly(&[0]), poly(&[0])],
            vec![poly(&[1]), poly(&[1]), poly(&[1]), poly(&[0])],
            vec![poly(&[0]), poly(&[0]), poly(&[0]), poly(&[5])],
        ];
        let matrix = MatrixA::new(rows);
        let matrix_snapshot = matrix.clone();

        let rhs = [poly(&[1, 1]), poly(&[2]), poly(&[0, 1]), poly(&[1])];
        let rhs_snapshot = rhs.clone();

        let _ = (&matrix) * &rhs;

        assert_eq!(matrix, matrix_snapshot);
        assert_eq!(rhs, rhs_snapshot);
    }
}
