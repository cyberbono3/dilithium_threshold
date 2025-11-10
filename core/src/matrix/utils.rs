use math::{Matrix, poly::Polynomial, traits::FiniteField};
use num_traits::Zero;

use crate::dilithium::params::{K, L, N, Q};
use crate::matrix::hash::shake128;

/// Convert a slice of polynomials into owned `'static` polynomials.
pub fn to_owned_polynomials<'poly, FF: FiniteField>(
    slice: &[Polynomial<'poly, FF>],
) -> Vec<Polynomial<'static, FF>> {
    slice
        .iter()
        .map(|poly| Polynomial::from(poly.coefficients().to_vec()))
        .collect()
}

const MATRIX_POLY_OP: &str = "matrix polynomial multiplication";

/// Ensure the matrix column count matches the vector length.
pub fn ensure_cols_match(
    matrix_rows: usize,
    matrix_cols: usize,
    vec_len: usize,
) -> Result<(), math::error::MatrixError> {
    if matrix_cols != vec_len {
        return Err(math::error::MatrixError::VectorShapeMismatch {
            operation: MATRIX_POLY_OP,
            matrix_rows,
            matrix_cols,
            vector_len: vec_len,
        });
    }
    Ok(())
}

/// Multiply rows by the provided vector, returning the resulting polynomials.
pub fn multiply_rows<FF: FiniteField>(
    rows: &[Vec<Polynomial<'static, FF>>],
    vec: &[Polynomial<'static, FF>],
) -> Vec<Polynomial<'static, FF>> {
    rows.iter()
        .map(|row| {
            row.iter().zip(vec.iter()).fold(
                Polynomial::<FF>::zero(),
                |mut acc, (a, b)| {
                    acc += a.clone() * b.clone();
                    acc
                },
            )
        })
        .collect()
}

/// Expand the public matrix `A` from the seed `rho` using SHAKE128 (educational variant).
pub fn expand_a_from_rho<FF: FiniteField + From<i64>>(
    rho: [u8; 32],
) -> Matrix<'static, FF> {
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
    use math::field_element::FieldElement;

    type FE = FieldElement;

    #[test]
    fn multiply_rows_matches_manual_computation() {
        let poly = |value: u32| -> Polynomial<'static, FE> {
            Polynomial::from(vec![FE::from(value)])
        };

        let a11 = poly(1);
        let a12 = poly(2);
        let a21 = poly(3);
        let a22 = poly(4);
        let x0 = poly(5);
        let x1 = poly(6);

        let rows = vec![
            vec![a11.clone(), a12.clone()],
            vec![a21.clone(), a22.clone()],
        ];
        let input = vec![x0.clone(), x1.clone()];

        let result = multiply_rows(&rows, &input);
        let expected = vec![
            a11.clone() * x0.clone() + a12.clone() * x1.clone(),
            a21.clone() * x0 + a22.clone() * x1,
        ];

        assert_eq!(result, expected);
    }

    #[test]
    fn expand_a_changes_with_rho() {
        let a1 = expand_a_from_rho::<FE>([1u8; 32]);
        let a2 = expand_a_from_rho::<FE>([2u8; 32]);
        assert_ne!(a1, a2);
    }
}
