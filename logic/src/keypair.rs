use crate::hash::shake256;
use crate::matrix::{MatrixA, expand_a_from_rho};
use crate::params::{K, L, N};
use crate::utils::{random_bytes, zero_polyvec};
use math::{poly::Polynomial, traits::FiniteField};
use std::ops::Mul;

#[derive(Clone, Debug, PartialEq)]
pub struct PublicKey<'a, FF: FiniteField> {
    pub a: MatrixA<'a, FF>, // uncompressed: include A directly
    pub t: [Polynomial<'a, FF>; K], // t = A*s1 + s2
    pub rho: [u8; 32],      // seed used for A (kept for provenance)
}

impl<'a, FF: FiniteField> PublicKey<'a, FF> {
    /// Construct a public key from its parts.
    /// Note: this does not verify that `t == A*s1 + s2`—it just packages fields.
    pub fn new(
        a: MatrixA<'a, FF>,
        t: [Polynomial<'a, FF>; K],
        rho: [u8; 32],
    ) -> Self {
        Self { a, t, rho }
    }
}

#[derive(Clone, Debug)]
pub struct PrivateKey<'a, FF: FiniteField> {
    pub a: MatrixA<'a, FF>, // include A here for convenience
    pub s1: [Polynomial<'a, FF>; L],
    pub s2: [Polynomial<'a, FF>; K],
}

impl<'a, FF: FiniteField> PrivateKey<'a, FF> {
    /// Construct a secret key from its parts.
    pub fn new(
        a: MatrixA<'a, FF>,
        s1: [Polynomial<'a, FF>; L],
        s2: [Polynomial<'a, FF>; K],
    ) -> Self {
        Self { a, s1, s2 }
    }
}

/// CBD for η=2 from an XOF stream.
fn cbd_eta2<FF: FiniteField + From<i64>>(
    stream: &[u8],
) -> Polynomial<'static, FF> {
    // Each coefficient uses 4 bits: (b0 + b1) - (b2 + b3)
    let mut out = [0i64; N];
    let mut bitpos = 0usize;
    for out_i in out.iter_mut().take(N) {
        let byte_idx = bitpos / 8;
        let two =
            u16::from_le_bytes([stream[byte_idx], stream[byte_idx + 1]]) as u32;
        let bits = two;
        let a0 = (bits & 1) + ((bits >> 1) & 1);
        let a1 = ((bits >> 2) & 1) + ((bits >> 3) & 1);
        *out_i = (a0 as i64) - (a1 as i64); // in [-2,2]
        bitpos += 4;
    }
    out.into()
}

pub fn keygen<FF: FiniteField + From<i64>>()
-> (PublicKey<'static, FF>, PrivateKey<'static, FF>) {
    let rho = random_bytes();
    let s1_seed = random_bytes();
    let s2_seed = random_bytes();

    keygen_with_seeds::<FF>(rho, s1_seed, s2_seed)
}

#[inline]
fn expand_secret_array<const LEN: usize, FF>(
    seed: impl AsRef<[u8]>,
) -> [Polynomial<'static, FF>; LEN]
where
    FF: math::traits::FiniteField + From<i64>,
{
    let seed = seed.as_ref();
    let mut input = Vec::with_capacity(seed.len() + 2);
    input.extend_from_slice(seed);
    input.extend_from_slice(&[0u8; 2]);

    let mut secrets = zero_polyvec::<LEN, FF>();
    for (j, slot) in secrets.iter_mut().enumerate() {
        let idx_bytes = (j as u16).to_le_bytes();
        let offset = input.len() - 2;
        input[offset..].copy_from_slice(&idx_bytes);
        let bytes = shake256(2 * N, &input);
        *slot = cbd_eta2::<FF>(&bytes);
    }

    secrets
}

pub fn keygen_with_seeds<FF: FiniteField + From<i64>>(
    rho: [u8; 32],
    s1_seed: [u8; 32],
    s2_seed: [u8; 32],
) -> (PublicKey<'static, FF>, PrivateKey<'static, FF>) {
    let a = expand_a_from_rho(rho);
    let s1 = expand_secret_array::<L, FF>(s1_seed);
    let s2 = expand_secret_array::<K, FF>(s2_seed);

    let mut t = a.mul(&s1);
    t.iter_mut()
        .zip(s2.iter())
        .for_each(|(dest, addend)| *dest += addend.clone());

    let public_key = PublicKey::new(a.clone(), t, rho);
    let private_key = PrivateKey::new(a, s1, s2);

    (public_key, private_key)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::expand_a_from_rho;
    use crate::params::{K, L};
    use crate::utils::zero_polyvec;
    use math::field_element::FieldElement;

    #[test]
    fn keygen_shapes_and_secret_bounds() {
        let (pk, sk) = keygen::<FieldElement>();

        // Shapes
        assert_eq!(pk.a.rows.len(), K);
        for i in 0..K {
            assert_eq!(pk.a.rows[i].len(), L);
        }
        assert_eq!(pk.t.len(), K);
        assert_eq!(sk.s1.len(), L);
        assert_eq!(sk.s2.len(), K);
    }

    #[test]
    fn public_matrix_matches_rho_and_secret_matrix() {
        let (pk, sk) = keygen::<FieldElement>();

        // pk.a should equal expand_a_from_rho(pk.rho)
        let a_from_rho = expand_a_from_rho(pk.rho);
        for i in 0..K {
            for j in 0..L {
                assert_eq!(
                    pk.a.rows[i][j], a_from_rho.rows[i][j],
                    "A mismatch at {},{}",
                    i, j
                );
            }
        }

        // pk.a and sk.a should match
        for i in 0..K {
            for j in 0..L {
                assert_eq!(
                    pk.a.rows[i][j], sk.a.rows[i][j],
                    "A(pk) != A(sk) at {},{}",
                    i, j
                );
            }
        }
    }

    #[test]
    fn t_equals_a_times_s1_plus_s2() {
        let (pk, sk) = keygen::<FieldElement>();

        let as1 = sk.a.mul(&sk.s1);
        let mut expected_t = zero_polyvec::<K, FieldElement>();
        for i in 0..K {
            let mut sum = as1[i].clone();
            sum += sk.s2[i].clone();
            expected_t[i] = sum;
        }

        for i in 0..K {
            assert_eq!(pk.t[i], expected_t[i], "t mismatch at row {}", i);
        }
    }

    #[test]
    fn keygen_is_deterministic_given_seeds() {
        let rho = [3u8; 32];
        let s1 = [5u8; 32];
        let s2 = [7u8; 32];
        let (pk1, sk1) = keygen_with_seeds::<FieldElement>(rho, s1, s2);
        let (pk2, sk2) = keygen_with_seeds::<FieldElement>(rho, s1, s2);

        // Matrices and secrets identical
        assert_eq!(pk1.a.rows, pk2.a.rows);
        for i in 0..K {
            assert_eq!(sk1.s2[i], sk2.s2[i]);
        }
        for j in 0..L {
            assert_eq!(sk1.s1[j], sk2.s1[j]);
        }
        for i in 0..K {
            assert_eq!(pk1.t[i], pk2.t[i]);
        }
    }
}
