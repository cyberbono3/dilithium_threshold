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

/// Shared abstraction for types that expose a participant id and associated polynomial vector.
pub trait PolyVectorSource<FF: FiniteField> {
    fn participant_id(&self) -> usize;
    fn polynomials(&self) -> &[Polynomial<'static, FF>];
}

impl<FF, T> PointSource<FF> for T
where
    FF: FiniteField + 'static,
    T: PolyVectorSource<FF>,
{
    fn x(&self) -> FF {
        let id = self.participant_id();
        (id as u64).into()
    }

    fn poly_at(&self, index: usize) -> Option<&Polynomial<'static, FF>> {
        self.polynomials().get(index)
    }

    fn poly_count(&self) -> usize {
        self.polynomials().len()
    }
}
