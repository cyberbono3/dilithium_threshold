use std::convert::TryInto;

use math::{error::MatrixError, poly::Polynomial, traits::FiniteField, Matrix};
use num_traits::Zero;

use crate::dilithium::params::{K, L, N, Q};
use crate::matrix::hash::shake128;

/// Type alias to the shared matrix implementation in the math crate.
pub type MatrixA<'a, FF> = Matrix<'a, FF>;

/// Types that can be produced from matrix-vector multiplication results.
pub trait MatrixMulOutput<FF: FiniteField>: Sized {
    fn from_vec(
        vec: Vec<Polynomial<'static, FF>>,
    ) -> Result<Self, MatrixError>;
}

impl<FF: FiniteField> MatrixMulOutput<FF> for Vec<Polynomial<'static, FF>> {
    fn from_vec(
        vec: Vec<Polynomial<'static, FF>>,
    ) -> Result<Self, MatrixError> {
        Ok(vec)
    }
}

impl<FF: FiniteField> MatrixMulOutput<FF> for [Polynomial<'static, FF>; K] {
    fn from_vec(
        vec: Vec<Polynomial<'static, FF>>,
    ) -> Result<Self, MatrixError> {
        vec.try_into().map_err(|vec: Vec<Polynomial<'static, FF>>| {
            MatrixError::OutputLengthMismatch {
                expected: K,
                found: vec.len(),
            }
        })
    }
}

/// Additional helpers used by the Dilithium core on top of the generic matrix.
pub trait MatrixAExt<FF: FiniteField> {
    /// Multiply the matrix by a collection of polynomials and return the requested container.
    fn mul_polynomials<'poly, Input, Output>(
        &self,
        vec: Input,
    ) -> Result<Output, MatrixError>
    where
        Input: AsRef<[Polynomial<'poly, FF>]>,
        Output: MatrixMulOutput<FF>,
        FF: 'poly;
}

impl<FF: FiniteField + 'static> MatrixAExt<FF> for Matrix<'static, FF> {
    fn mul_polynomials<'poly, Input, Output>(
        &self,
        vec: Input,
    ) -> Result<Output, MatrixError>
    where
        Input: AsRef<[Polynomial<'poly, FF>]>,
        Output: MatrixMulOutput<FF>,
        FF: 'poly,
    {
        let owned: Vec<Polynomial<'static, FF>> = vec
            .as_ref()
            .iter()
            .map(|poly| Polynomial::from(poly.coefficients().to_vec()))
            .collect();
        let cols = self.cols();
        if cols != owned.len() {
            return Err(MatrixError::VectorShapeMismatch {
                matrix_cols: cols,
                vector_len: owned.len(),
            });
        }

        let result = self
            .as_slice()
            .iter()
            .map(|row| {
                row.iter()
                    .zip(owned.iter())
                    .fold(Polynomial::<FF>::zero(), |mut acc, (a_ij, v_j)| {
                        acc += a_ij.clone() * v_j.clone();
                        acc
                    })
            })
            .collect();

        Output::from_vec(result)
    }
}

/// Expand the public matrix `A` from the seed `rho` using SHAKE128 (educational variant).
pub fn expand_a_from_rho<FF: FiniteField + From<i64>>(
    rho: [u8; 32],
) -> MatrixA<'static, FF> {
    let modulus = Q as u32;
    let mut rows = Vec::with_capacity(K);
    for i in 0..K {
        let mut row = Vec::with_capacity(L);
        for j in 0..L {
            // Domain separation: rho || i || j
            let mut seed = Vec::with_capacity(32 + 4);
            seed.extend_from_slice(&rho);
            seed.extend_from_slice(&(i as u16).to_le_bytes());
            seed.extend_from_slice(&(j as u16).to_le_bytes());
            let stream = shake128(4 * N, &seed); // 4 bytes per coefficient
            let coeffs = (0..N)
                .map(|index| {
                    let start = 4 * index;
                    let word = u32::from_le_bytes([
                        stream[start],
                        stream[start + 1],
                        stream[start + 2],
                        stream[start + 3],
                    ]);
                    (word % modulus) as i64
                })
                .map(FF::from)
                .collect::<Vec<_>>();
            row.push(Polynomial::from(coeffs));
        }
        rows.push(row);
    }
    Matrix::new(rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dilithium::params::{K, L};
    use math::field_element::FieldElement;
    use math::poly::Polynomial;
    use num_traits::Zero;

    type FE = FieldElement;

    #[test]
    fn mul_polynomials_matches_manual_computation() {
        let matrix = expand_a_from_rho::<FE>([7u8; 32]);
        let input: Vec<Polynomial<'static, FE>> =
            vec![Polynomial::zero(); L];

        let slice_result: Vec<_> =
            matrix.mul_polynomials(input.as_slice()).unwrap();
        let manual: Vec<_> = matrix
            .as_slice()
            .iter()
            .map(|row| {
                row.iter()
                    .zip(input.iter())
                    .fold(Polynomial::zero(), |mut acc, (a, b)| {
                        acc += a.clone() * b.clone();
                        acc
                    })
            })
            .collect();

        assert_eq!(slice_result, manual);
    }

    #[test]
    fn mul_polynomials_computes_expected_values() {
        let poly = |value: u32| -> Polynomial<'static, FE> {
            Polynomial::from(vec![FE::from(value)])
        };

        let a11 = poly(1);
        let a12 = poly(2);
        let a21 = poly(3);
        let a22 = poly(4);
        let x0 = poly(5);
        let x1 = poly(6);

        let matrix = Matrix::new(vec![
            vec![a11.clone(), a12.clone()],
            vec![a21.clone(), a22.clone()],
        ]);
        let input = vec![x0.clone(), x1.clone()];

        let result: Vec<Polynomial<'static, FE>> =
            matrix.mul_polynomials(input.as_slice()).unwrap();

        let expected = vec![
            a11.clone() * x0.clone() + a12.clone() * x1.clone(),
            a21.clone() * x0 + a22.clone() * x1,
        ];

        assert_eq!(result, expected);
    }

    #[test]
    fn mul_polynomials_returns_fixed_length_array() {
        let matrix = expand_a_from_rho::<FE>([13u8; 32]);
        let vec: [Polynomial<'static, FE>; L] =
            std::array::from_fn(|_| Polynomial::zero());

        let array: [Polynomial<'static, FE>; K] =
            matrix.mul_polynomials(&vec).unwrap();
        assert_eq!(array.len(), K);
    }

    #[test]
    fn mul_polynomials_errors_on_shape_mismatch() {
        let matrix: Matrix<'static, FieldElement> = Matrix::zeros(2, 2);
        let vec = vec![Polynomial::zero(); 3];
        let result: Result<
            Vec<Polynomial<'static, FieldElement>>,
            MatrixError,
        > = matrix.mul_polynomials(vec.as_slice());
        let err = result.unwrap_err();
        assert!(matches!(
            err,
            MatrixError::VectorShapeMismatch {
                matrix_cols: 2,
                vector_len: 3
            }
        ));
    }

    #[test]
    fn mul_polynomials_array_output_mismatch_errors() {
        let matrix: Matrix<'static, FieldElement> = Matrix::zeros(K + 1, L);
        let vec: [Polynomial<'static, FieldElement>; L] =
            std::array::from_fn(|_| Polynomial::zero());

        let result: Result<[Polynomial<'static, FieldElement>; K], MatrixError> =
            matrix.mul_polynomials(&vec);

        let err = result.unwrap_err();
        assert!(matches!(
            err,
            MatrixError::OutputLengthMismatch {
                expected: K,
                found
            } if found == K + 1
        ));
    }

    #[test]
    fn expand_a_changes_with_rho() {
        let a1 = expand_a_from_rho::<FE>([1u8; 32]);
        let a2 = expand_a_from_rho::<FE>([2u8; 32]);
        assert_ne!(a1, a2);
    }
}
