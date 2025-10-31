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
    /// Initialize signature.
    pub fn new(
        z: Vec<Polynomial<'a, FF>>, // TODO consider to use arrays instead of length L
        h: Vec<Polynomial<'a, FF>>, // TODO consider to use arrays instead of length K
        c: Polynomial<'a, FF>,
    ) -> Self {
        DilithiumSignature { z, h, c }
    }
}
