use math::traits::FiniteField;

use math::prelude::*;

/// Dilithium signature containing vectors z and h, and challenge c.
//TODO replace z,c with Signature
#[derive(Clone, Debug, PartialEq)]
pub struct DilithiumSignature<'a, FF: FiniteField> {
    pub z: Vec<Polynomial<'a, FF>>,
    pub h: Vec<Polynomial<'a, FF>>,
    pub c: Polynomial<'a, FF>,
}

impl<'a, FF: FiniteField> DilithiumSignature<'a, FF> {
    /// Construct a signature from its response, hint, and challenge components.
    pub fn new(
        z: Vec<Polynomial<'a, FF>>, // TODO consider to use arrays of length L instead
        h: Vec<Polynomial<'a, FF>>, // TODO consider to use arrays of length K instead
        c: Polynomial<'a, FF>,
    ) -> Self {
        DilithiumSignature { z, h, c }
    }
}
