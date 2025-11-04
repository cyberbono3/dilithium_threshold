use std::convert::TryInto;

use math::{
    poly::Polynomial,
    poly_vector::PolynomialVector,
    traits::FiniteField,
    Matrix,
};

use crate::dilithium::params::{K, L, N, Q};
use crate::matrix::hash::shake128;

/// Type alias to the shared matrix implementation in the math crate.
pub type MatrixA<'a, FF> = Matrix<'a, FF>;

/// Additional helpers used by the Dilithium core on top of the generic matrix.
pub trait MatrixAExt<FF: FiniteField> {
    /// Multiply the matrix by a slice of polynomials, returning the resulting vector.
    fn mul_poly_slice(
        &self,
        vec: &[Polynomial<'_, FF>],
    ) -> Vec<Polynomial<'static, FF>>;

    /// Multiply the matrix by a fixed-length polynomial array, returning a fixed-length result.
    fn mul_poly_array(
        &self,
        vec: &[Polynomial<'_, FF>; L],
    ) -> [Polynomial<'static, FF>; K];
}

impl<FF: FiniteField + 'static> MatrixAExt<FF> for Matrix<'static, FF> {
    fn mul_poly_slice(
        &self,
        vec: &[Polynomial<'_, FF>],
    ) -> Vec<Polynomial<'static, FF>> {
        let owned: Vec<Polynomial<'static, FF>> = vec
            .iter()
            .map(|poly| Polynomial::from(poly.coefficients().to_vec()))
            .collect();
        let vector = PolynomialVector::new(owned);
        self.mul_vector(&vector).into_iter().collect()
    }

    fn mul_poly_array(
        &self,
        vec: &[Polynomial<'_, FF>; L],
    ) -> [Polynomial<'static, FF>; K] {
        self.mul_poly_slice(vec)
            .try_into()
            .expect("matrix multiplication should yield K entries")
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
    use math::poly_vector::PolynomialVector;
    use num_traits::Zero;

    type FE = FieldElement;

    #[test]
    fn mul_poly_slice_matches_polynomial_vector() {
        let matrix = expand_a_from_rho::<FE>([7u8; 32]);
        let vec = PolynomialVector::new(vec![Polynomial::zero(); L]);

        let slice_result = matrix.mul_poly_slice(vec.as_slice());
        let vector_result: Vec<_> = matrix.mul_vector(&vec).into_iter().collect();

        assert_eq!(slice_result, vector_result);
    }

    #[test]
    fn mul_poly_array_returns_fixed_length() {
        let matrix = expand_a_from_rho::<FE>([13u8; 32]);
        let vec = std::array::from_fn(|_| Polynomial::zero());

        let array = matrix.mul_poly_array(&vec);
        assert_eq!(array.len(), K);
    }

    #[test]
    #[should_panic(expected = "Matrix columns")]
    fn mul_poly_slice_panics_on_shape_mismatch() {
        let matrix: Matrix<'static, FieldElement> = Matrix::zeros(2, 2);
        let vec = vec![Polynomial::zero(); 3];
        let _ = matrix.mul_poly_slice(vec.as_slice());
    }

    #[test]
    fn expand_a_changes_with_rho() {
        let a1 = expand_a_from_rho::<FE>([1u8; 32]);
        let a2 = expand_a_from_rho::<FE>([2u8; 32]);
        assert_ne!(a1, a2);
    }
}
