use sha3::{
    Shake256,
    digest::{ExtendableOutput, Update, XofReader},
};

use crate::{
    params::{GAMMA1, K, L},
    utils::get_hash_reader,
};
use math::{matrix::Matrix, traits::FiniteField};

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

/// Sample polynomial with coefficients in [-gamma1, gamma1].
// TODO move it to utils section
// TODO add test cases
pub fn sample_gamma1<FF: FiniteField>(seed: &[u8]) -> Polynomial<'static, FF> {
    let mut reader = get_hash_reader(seed);

    let mut bytes = vec![0u8; N * 4];
    reader.read(&mut bytes);

    //let mut coeffs = vec![0i32; N];
    // TODO speed up
    let coeffs: Vec<FF> = (0..N)
        .map(|i| {
            let idx = i * 4;
            let val = u32::from_le_bytes([
                bytes[idx],
                bytes[idx + 1],
                bytes[idx + 2],
                bytes[idx + 3],
            ]);
            let sample = val % (2 * GAMMA1 as u32 + 1);
            let coeff = sample as i32 - GAMMA1 as i32;
            <i32 as Into<FF>>::into(coeff)
        })
        .collect();

    poly!(coeffs)
}
