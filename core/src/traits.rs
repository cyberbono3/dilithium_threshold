use math::{poly::Polynomial, traits::FiniteField};

/// Abstract source for polynomial evaluation points used during reconstruction.
pub trait PointSource<FF: FiniteField> {
    /// The x-coordinate corresponding to this point.
    fn x(&self) -> FF;

    /// Return the polynomial at a given vector index (if present).
    fn poly_at(&self, index: usize) -> Option<&Polynomial<'static, FF>>;

    /// Total number of polynomials available in this source.
    fn poly_count(&self) -> usize;
}
