use crate::basic::sign::{DilithiumSignature, SigningField};
use crate::basic::utils::{
    all_infty_norm_below, derive_challenge, pack_w1_for_hash,
    polyvec_sub_scaled, use_hints,
};
use crate::dilithium::params::{BETA, GAMMA1, K, L};
use crate::matrix::MatrixMulExt;
use math::{Matrix, poly::Polynomial, traits::FiniteField};

#[derive(Clone, Debug, PartialEq)]
pub struct PublicKey<'a, FF: FiniteField> {
    pub a: Matrix<'a, FF>, // uncompressed: include A directly
    pub t: [Polynomial<'a, FF>; K], // t = A*s1 + s2
    pub rho: [u8; 32],     // seed used for A (kept for provenance)
}

impl<'a, FF: FiniteField> PublicKey<'a, FF> {
    /// Construct a public key from its parts.
    /// Note: this does not verify that `t == A*s1 + s2`—it just packages fields.
    pub fn new(
        a: Matrix<'a, FF>,
        t: [Polynomial<'a, FF>; K],
        rho: [u8; 32],
    ) -> Self {
        Self { a, t, rho }
    }
}

impl<FF> PublicKey<'static, FF>
where
    FF: SigningField,
    i64: From<FF>,
{
    /// Verify `sig` against `msg` using only the public key.
    pub fn verify(&self, msg: &[u8], sig: &DilithiumSignature<'_, FF>) -> bool {
        if !all_infty_norm_below::<FF, L>(&sig.z, GAMMA1 - BETA) {
            return false;
        }

        let Some(az) = self.a.matrix_mul_output(&sig.z) else {
            return false;
        };
        let az_minus_ct = polyvec_sub_scaled::<FF, K>(&az, &sig.c, &self.t);
        let w1_prime = use_hints(&sig.h, &az_minus_ct);

        let derived = derive_challenge(msg, &pack_w1_for_hash(&w1_prime));
        derived == sig.c
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::keys::test_helpers;
    use crate::basic::{KeyPair, keygen_with_seeds};
    use crate::dilithium::params::{K, L};
    use crate::matrix::expand_a_from_rho;
    use math::{field_element::FieldElement, poly::Polynomial};

    #[test]
    fn public_key_has_expected_dimensions() {
        let KeyPair { public: pk, .. } = test_helpers::fixture_keypair();

        assert_eq!(pk.a.rows(), K);
        for row in pk.a.as_slice() {
            assert_eq!(row.len(), L);
        }
        assert_eq!(pk.t.len(), K);
    }

    #[test]
    fn public_matrix_matches_rho_seed() {
        let KeyPair { public: pk, .. } = test_helpers::fixture_keypair();

        let from_rho = expand_a_from_rho(pk.rho);
        assert_eq!(pk.a.as_slice(), from_rho.as_slice());
    }

    #[test]
    fn t_equals_a_times_s1_plus_s2() {
        let KeyPair {
            public: pk,
            private: sk,
        } = test_helpers::fixture_keypair();

        let as1: [Polynomial<'static, FieldElement>; K] =
            sk.a.mul_polynomials(&sk.s1)
                .expect("matrix-vector multiplication should succeed");
        let expected: Vec<_> = as1
            .into_iter()
            .zip(sk.s2.iter())
            .map(|(row, s2)| row + s2)
            .collect();

        assert_eq!(pk.t.as_slice(), expected.as_slice());
    }

    #[test]
    fn public_key_is_deterministic_given_seeds() {
        let seeds = test_helpers::fixture_seeds();
        let KeyPair { public: pk1, .. } =
            keygen_with_seeds::<FieldElement>(seeds)
                .expect("key generation should succeed");
        let KeyPair { public: pk2, .. } =
            keygen_with_seeds::<FieldElement>(seeds)
                .expect("key generation should succeed");

        assert_eq!(pk1, pk2);
    }
}
