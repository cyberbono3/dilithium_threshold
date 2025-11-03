use math::{prelude::*, traits::FiniteField};

pub trait PointSource<FF: FiniteField> {
    /// The x-coordinate of this participant/point (e.g., participant id as field element).
    fn x(&self) -> FF;

    /// Return the polynomial at a given vector index (if present).
    fn poly_at(&self, index: usize) -> Option<&Polynomial<'static, FF>>;

    /// Total number of polynomials in the vector (for bounds/error reporting).
    fn poly_count(&self) -> usize;
}
