use crate::hash::shake128;
use crate::params::{K, L, N};
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

        self.rows
            .iter()
            .map(|row| {
                row.iter().zip(v).fold(
                    Polynomial::<FF>::zero(),
                    |mut acc, (a_ij, v_j)| {
                        // TODO: enable &Polynomial * &Polynomial to avoid clones
                        acc += a_ij.clone() * v_j.clone();
                        acc
                    },
                )
            })
            .collect()
    }
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

// TODO declare the trait
pub fn mat_vec_mul<FF: FiniteField>(
    a: &MatrixA<FF>,
    y: &[Polynomial<FF>; L],
) -> [Polynomial<'static, FF>; K] {
    std::array::from_fn(|i| {
        a.rows[i].iter().zip(y.iter()).fold(
            Polynomial::zero(),
            |mut acc, (aij, yj)| {
                acc += aij.clone() * yj.clone();
                acc
            },
        )
    })
}

// impl<FF: FiniteField> Mul<&Polynomial<'static, FF>> for Vec<Polynomial<'static, FF>> {
//     type Output = Self;

//     fn mul(self, poly: Self) -> Self::Output {
//         assert!(
//             self.len() == poly.len(),
//             "Vector lengths must match for element-wise multiplication"
//         );

//         self
//             .into_iter()
//             .map(|p1| p1 * poly)
//             .collect()

//     }
// }

// TODO add more tests
#[cfg(test)]
mod tests {
    use crate::matrix::{expand_a_from_rho, mat_vec_mul};
    use crate::params::{K, L, N, Q};
    //use crate::poly::{Poly, mod_q};
    use math::poly::Polynomial;
    use num_traits::Zero;

    use math::field_element::FieldElement;
    use math::traits::FiniteField;

    fn zero_y<FF: FiniteField>() -> [Polynomial<'static, FF>; L] {
        [
            Polynomial::zero(),
            Polynomial::zero(),
            Polynomial::zero(),
            Polynomial::zero(),
        ]
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
    fn mat_vec_mul_zero_is_zero() {
        let a = expand_a_from_rho::<FieldElement>([42u8; 32]);
        let y = zero_y();
        let w = mat_vec_mul(&a, &y);
        for i in 0..K {
            assert_eq!(w[i], Polynomial::zero());
        }
    }
}
