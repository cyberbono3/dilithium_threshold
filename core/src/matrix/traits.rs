use math::{
    error::MatrixError,
    poly::Polynomial,
    poly_vector::PolynomialVector,
    traits::FiniteField,
    Matrix,
};

use super::utils::{ensure_cols_match, multiply_rows, to_owned_polynomials};

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

impl<FF: FiniteField, const LEN: usize> MatrixMulOutput<FF>
    for [Polynomial<'static, FF>; LEN]
{
    fn from_vec(
        vec: Vec<Polynomial<'static, FF>>,
    ) -> Result<Self, MatrixError> {
        vec.try_into().map_err(|vec: Vec<Polynomial<'static, FF>>| {
            MatrixError::OutputLengthMismatch {
                expected: LEN,
                found: vec.len(),
            }
        })
    }
}

impl<FF: FiniteField> MatrixMulOutput<FF> for PolynomialVector<'static, FF> {
    fn from_vec(
        vec: Vec<Polynomial<'static, FF>>,
    ) -> Result<Self, MatrixError> {
        Ok(PolynomialVector::new(vec))
    }
}

/// Additional helpers used by the Dilithium core on top of the generic matrix.
pub trait MatrixMulExt<FF: FiniteField> {
    /// Multiply the matrix by a collection of polynomials and return the requested container.
    fn mul_polynomials<'poly, Input, Output>(
        &self,
        vec: Input,
    ) -> Result<Output, MatrixError>
    where
        Input: AsRef<[Polynomial<'poly, FF>]>,
        Output: MatrixMulOutput<FF>,
        FF: 'poly;

    /// Convenience helper returning `None` when the multiplication fails.
    fn matrix_mul_output<'poly, Input, Output>(
        &self,
        vec: Input,
    ) -> Option<Output>
    where
        Input: AsRef<[Polynomial<'poly, FF>]>,
        Output: MatrixMulOutput<FF>,
        FF: 'poly;
}

impl<FF: FiniteField + 'static> MatrixMulExt<FF> for Matrix<'static, FF> {
    fn mul_polynomials<'poly, Input, Output>(
        &self,
        vec: Input,
    ) -> Result<Output, MatrixError>
    where
        Input: AsRef<[Polynomial<'poly, FF>]>,
        Output: MatrixMulOutput<FF>,
        FF: 'poly,
    {
        let owned = to_owned_polynomials(vec.as_ref());
        ensure_cols_match(self.rows(), self.cols(), owned.len())?;
        Output::from_vec(multiply_rows(self.as_slice(), &owned))
    }

    fn matrix_mul_output<'poly, Input, Output>(
        &self,
        vec: Input,
    ) -> Option<Output>
    where
        Input: AsRef<[Polynomial<'poly, FF>]>,
        Output: MatrixMulOutput<FF>,
        FF: 'poly,
    {
        self.mul_polynomials(vec).ok()
    }
}
