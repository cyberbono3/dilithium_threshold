use crate::basic::keypair::PrivateKey;
use crate::dilithium::params::{BETA, GAMMA1, GAMMA2, K, L};
use crate::matrix::{MatrixMulExt, hash::shake256};

use math::prelude::FieldElement;
use math::{poly::Polynomial, traits::FiniteField};
use std::array::from_fn;
/// Uncompressed Dilithium signature with response vector `z`, hint vector `h`, and challenge `c`.
#[derive(Clone, Debug, PartialEq)]
pub struct DilithiumSignature<'a, FF: FiniteField> {
    pub z: [Polynomial<'a, FF>; L],
    pub h: [Polynomial<'a, FF>; K],
    pub c: Polynomial<'a, FF>, // challenge polynomial with TAU non-zeros in {-1,0,1}
}

impl<'a, FF: FiniteField> DilithiumSignature<'a, FF> {
    /// Construct a signature from its response, hint, and challenge components.
    pub fn new(
        z: [Polynomial<'a, FF>; L],
        h: [Polynomial<'a, FF>; K],
        c: Polynomial<'a, FF>,
    ) -> Self {
        Self { z, h, c }
    }
}

pub trait SigningField:
    FiniteField + From<i64> + Into<[u8; FieldElement::BYTES]> + 'static
{
}

impl<T> SigningField for T where
    T: FiniteField + From<i64> + Into<[u8; FieldElement::BYTES]> + 'static
{
}

use crate::basic::utils::{
    all_infty_norm_below, derive_challenge, make_hints, pack_w1_for_hash,
    poly_high, poly_low, polyvec_add_scaled, polyvec_sub_scaled, sample_y,
};

/// Internal helper driving the rejection-sampling signing loop.
pub(crate) struct SigningEngine<'pk, 'msg, FF>
where
    FF: SigningField,
    i64: From<FF>,
{
    priv_key: &'pk PrivateKey<'static, FF>,
    msg: &'msg [u8],
    y_seed: Vec<u8>,
}

impl<'pk, 'msg, FF> SigningEngine<'pk, 'msg, FF>
where
    FF: SigningField,
    i64: From<FF>,
{
    pub(crate) const REJECTION_LIMIT: u32 = 10000;

    #[inline]
    pub(crate) fn rejection_attempts() -> std::ops::Range<u32> {
        0..Self::REJECTION_LIMIT
    }

    /// Construct a new engine for the supplied private key and message.
    pub(crate) fn new(
        priv_key: &'pk PrivateKey<'static, FF>,
        msg: &'msg [u8],
    ) -> Self {
        Self {
            priv_key,
            msg,
            y_seed: shake256(32, msg),
        }
    }

    /// Attempt one signing iteration; returns `None` if bounds are violated.
    pub(crate) fn try_with_counter(
        &self,
        ctr: u32,
    ) -> Option<DilithiumSignature<'static, FF>> {
        let y = sample_y::<FF>(&self.y_seed, ctr);
        let w: [Polynomial<'static, FF>; K] =
            self.priv_key.a.matrix_mul_output(&y)?;
        let w1 = from_fn(|idx| poly_high(&w[idx]));
        let challenge = derive_challenge(self.msg, &pack_w1_for_hash(&w1));

        let z = polyvec_add_scaled::<FF, L>(&y, &challenge, &self.priv_key.s1);
        if !all_infty_norm_below::<FF, L>(&z, GAMMA1 - BETA) {
            return None;
        }

        let ay_minus_cs2 =
            polyvec_sub_scaled::<FF, K>(&w, &challenge, &self.priv_key.s2);
        let w0 = from_fn(|idx| poly_low(&ay_minus_cs2[idx]));
        if !all_infty_norm_below::<FF, K>(&w0, GAMMA2 - BETA) {
            return None;
        }

        let hints = make_hints::<FF>(&w1, &ay_minus_cs2);

        Some(DilithiumSignature {
            z,
            h: hints,
            c: challenge,
        })
    }
}

pub(crate) type SigningEngineConfig<FF> = SigningEngine<'static, 'static, FF>;

#[cfg(test)]
mod tests {
    use super::*; // brings sign/verify + private helpers into scope
    use crate::basic::keypair::*;
    use crate::basic::utils;
    use crate::dilithium::params::{BETA, GAMMA1, L, N};
    use crate::dilithium::utils::zero_polyvec;
    use math::field_element::FieldElement;

    /// Deterministic fixture producing a reproducible keypair.
    fn keypair_fixture_1<FF: FiniteField + From<i64>>() -> KeyPair<'static, FF>
    {
        let rho = [0x42u8; 32];
        let s1 = [0x24u8; 32];
        let s2 = [0x18u8; 32];
        keygen_with_seeds::<FF>(KeypairSeeds::new(rho, s1, s2))
            .expect("key generation should succeed")
    }

    /// Alternate deterministic keypair fixture used across tests.
    fn keypair_fixture_2<FF: FiniteField + From<i64>>() -> KeyPair<'static, FF>
    {
        let rho = [0xA5u8; 32];
        let s1 = [0x5Au8; 32];
        let s2 = [0x33u8; 32];
        keygen_with_seeds::<FF>(KeypairSeeds::new(rho, s1, s2))
            .expect("key generation should succeed")
    }

    mod hint_tests {
        use super::*;
        use crate::dilithium::params::ALPHA;

        #[test]
        fn hints_restore_reference_high_bits() {
            let mut reference = zero_polyvec::<K, FieldElement>();
            let mut base = zero_polyvec::<K, FieldElement>();
            // Choose coefficients that differ by a single ALPHA window.
            reference[0] = vec![FieldElement::from(ALPHA + 5)].into();
            base[0] = vec![FieldElement::from(5)].into();

            let target_high =
                std::array::from_fn(|idx| poly_high(&reference[idx]));
            let hints = utils::make_hints(&target_high, &base);
            let recovered = utils::use_hints(&hints, &base);

            assert_eq!(recovered[0], target_high[0]);
        }
    }

    mod signing_engine_tests {
        use super::*;

        /// Ensure the signing engine produces a signature before hitting the rejection limit.
        #[test]
        fn engine_produces_signatures_within_limit() {
            let message = b"engine-test";
            let KeyPair {
                public: pub_key,
                private: priv_key,
            } = keypair_fixture_1::<FieldElement>();
            let engine = SigningEngine::new(&priv_key, message);

            let mut produced = None;
            for ctr in
                super::SigningEngineConfig::<FieldElement>::rejection_attempts()
            {
                if let Some(sig) = engine.try_with_counter(ctr) {
                    produced = Some(sig);
                    break;
                }
            }

            let signature =
                produced.expect("expected a signature before limit");
            assert!(
                utils::all_infty_norm_below::<FieldElement, L>(
                    &signature.z,
                    GAMMA1 - BETA
                ),
                "z exceeds bound"
            );
            assert!(
                pub_key.verify(message, &signature),
                "engine output should verify against the public key"
            );
        }

        /// Check that the engine yields identical signatures for identical inputs.
        #[test]
        fn engine_is_deterministic_for_same_inputs() {
            let message = b"deterministic-engine";
            let seeds =
                KeypairSeeds::new([0x42u8; 32], [0x24u8; 32], [0x18u8; 32]);

            let priv_key_a = keygen_with_seeds::<FieldElement>(seeds)
                .expect("key generation should succeed")
                .private;
            let priv_key_b = keygen_with_seeds::<FieldElement>(seeds)
                .expect("key generation should succeed")
                .private;

            let engine_a = SigningEngine::new(&priv_key_a, message);
            let engine_b = SigningEngine::new(&priv_key_b, message);

            let ctr = 5;
            assert_eq!(
                engine_a.try_with_counter(ctr),
                engine_b.try_with_counter(ctr),
                "same inputs must yield identical signatures"
            );
        }

        /// Ensure the engine rejects samples when masks violate bounds.
        #[test]
        fn engine_rejects_when_mask_norm_exceeds_bounds() {
            let message = b"rejection-test";
            let KeyPair {
                public: _,
                private: mut priv_key,
            } = keypair_fixture_1::<FieldElement>();

            for poly in priv_key.s1.iter_mut() {
                let coeffs = vec![FieldElement::from(GAMMA1 + BETA + 1); N];
                *poly = coeffs.into();
            }

            let engine = SigningEngine::new(&priv_key, message);

            for ctr in 0..10 {
                assert!(
                    engine.try_with_counter(ctr).is_none(),
                    "expected rejection for counter {ctr}"
                );
            }
        }
    }

    mod sample_y_tests {
        use super::*;

        /// Convert a field element into a centered signed integer.
        fn centered_value(fe: FieldElement) -> i64 {
            let mut v = i64::from(fe);
            let p = FieldElement::P as i64;
            if v > p / 2 {
                v -= p;
            }
            v
        }

        /// Return true if all coefficients stay within the ±GAMMA1 bound.
        fn coeffs_within_bounds(
            polys: &[Polynomial<'_, FieldElement>],
        ) -> bool {
            polys
                .iter()
                .flat_map(|poly| poly.coefficients())
                .map(|&c| centered_value(c).abs())
                .all(|abs| abs <= GAMMA1)
        }

        /// Sampling with identical seed and counter should be deterministic.
        #[test]
        fn deterministic_per_seed_and_counter() {
            let seed = [0x11u8; 32];
            let first = utils::sample_y::<FieldElement>(&seed, 7);
            let second = utils::sample_y::<FieldElement>(&seed, 7);
            assert_eq!(first, second);
        }

        /// Different counters must yield distinct sampled vectors.
        #[test]
        fn different_counters_produce_distinct_vectors() {
            let seed = [0x77u8; 32];
            let first = utils::sample_y::<FieldElement>(&seed, 1);
            let second = utils::sample_y::<FieldElement>(&seed, 2);
            assert_ne!(first, second);
        }

        /// All coefficients should remain centered within the permitted bound.
        #[test]
        fn coefficients_are_centered() {
            let seed = [0xC3u8; 32];
            let polys = utils::sample_y::<FieldElement>(&seed, 5);
            assert!(coeffs_within_bounds(&polys));
        }
    }

    mod pack_w1_tests {
        use super::*;

        /// Helper to create a polynomial from integer coefficients.
        fn simple_poly(coeffs: &[i64]) -> Polynomial<'static, FieldElement> {
            Polynomial::from(
                coeffs
                    .iter()
                    .map(|&c| FieldElement::from(c))
                    .collect::<Vec<_>>(),
            )
        }

        /// Packing should be deterministic when the input polynomials are fixed.
        #[test]
        fn deterministic_for_fixed_polys() {
            let polys = [
                simple_poly(&[1, 2, 3]),
                simple_poly(&[4, 5, 6]),
                simple_poly(&[7, 8, 9]),
                simple_poly(&[10, 11, 12]),
            ];

            let first = utils::pack_w1_for_hash(&polys);
            let second = utils::pack_w1_for_hash(&polys);
            assert_eq!(first, second);
        }

        /// Different polynomials must lead to different packed byte arrays.
        #[test]
        fn different_inputs_produce_different_outputs() {
            let polys_a = [
                simple_poly(&[1, 2, 3]),
                simple_poly(&[4, 5, 6]),
                simple_poly(&[7, 8, 9]),
                simple_poly(&[10, 11, 12]),
            ];

            let polys_b = [
                simple_poly(&[1, 2, 3]),
                simple_poly(&[4, 5, 6]),
                simple_poly(&[7, 8, 9]),
                simple_poly(&[10, 11, 13]),
            ];

            let hash_a = utils::pack_w1_for_hash(&polys_a);
            let hash_b = utils::pack_w1_for_hash(&polys_b);
            assert_ne!(hash_a, hash_b);
        }
    }

    mod derive_challenge_tests {
        use super::*;
        use crate::dilithium::params::TAU;
        use num_traits::ConstZero;

        /// Derive challenge should be deterministic for identical inputs.
        #[test]
        fn deterministic_for_same_inputs() {
            let msg = b"challenge";
            let hash = vec![0x42; 32];
            let c1 = utils::derive_challenge::<FieldElement>(msg, &hash);
            let c2 = utils::derive_challenge::<FieldElement>(msg, &hash);
            assert_eq!(c1, c2);
        }

        /// Different messages should produce different challenge polynomials.
        #[test]
        fn changing_message_changes_challenge() {
            let hash = vec![0x77; 64];
            let c1 = utils::derive_challenge::<FieldElement>(b"m1", &hash);
            let c2 = utils::derive_challenge::<FieldElement>(b"m2", &hash);
            assert_ne!(c1, c2);
        }

        /// The derived challenge must contain exactly TAU non-zero coefficients.
        #[test]
        fn challenge_has_tau_non_zero_entries() {
            let msg = b"nonzero-count";
            let hash = vec![0xAB; 128];
            let challenge = utils::derive_challenge::<FieldElement>(msg, &hash);
            let non_zero = challenge
                .coefficients()
                .iter()
                .filter(|&&c| !FieldElement::ZERO.eq(&c))
                .count();
            assert_eq!(non_zero, TAU);
        }
    }

    mod infty_norm_tests {
        use super::*;

        /// Helper to construct polynomials for infinity norm tests.
        fn make_poly(coeffs: &[i64]) -> Polynomial<'static, FieldElement> {
            Polynomial::from(
                coeffs
                    .iter()
                    .map(|&c| FieldElement::from(c))
                    .collect::<Vec<_>>(),
            )
        }

        #[test]
        fn polyvec_add_scaled_adds_scaled_component() {
            let base = [
                Polynomial::from(vec![FieldElement::from(1i32)]),
                Polynomial::from(vec![FieldElement::from(2i32)]),
            ];
            let mult = [
                Polynomial::from(vec![FieldElement::from(3i32)]),
                Polynomial::from(vec![FieldElement::from(-1i32)]),
            ];
            let scale = Polynomial::from(vec![FieldElement::from(2i32)]);

            let result = utils::polyvec_add_scaled::<FieldElement, 2>(
                &base, &scale, &mult,
            );

            let expected = [
                Polynomial::from(vec![FieldElement::from(7i32)]),
                Polynomial::from(vec![FieldElement::from(0i32)]),
            ];

            assert_eq!(result, expected);
        }

        #[test]
        fn polyvec_sub_scaled_subtracts_scaled_component() {
            let base = [
                Polynomial::from(vec![FieldElement::from(5i32)]),
                Polynomial::from(vec![FieldElement::from(-3i32)]),
            ];
            let mult = [
                Polynomial::from(vec![FieldElement::from(2i32)]),
                Polynomial::from(vec![FieldElement::from(4i32)]),
            ];
            let scale = Polynomial::from(vec![FieldElement::from(3i32)]);

            let result = utils::polyvec_sub_scaled::<FieldElement, 2>(
                &base, &scale, &mult,
            );

            let expected = [
                Polynomial::from(vec![FieldElement::from(-1i32)]),
                Polynomial::from(vec![FieldElement::from(-15i32)]),
            ];

            assert_eq!(result, expected);
        }

        /// Detect cases where all polynomials stay within the supplied bound.
        #[test]
        fn detects_within_bound() {
            let polys = [make_poly(&[1, -2, 3]), make_poly(&[0, 0, 4])];
            assert!(utils::all_infty_norm_below(&polys, 5));
        }

        /// Detect cases where the infinity norm exceeds the bound.
        #[test]
        fn detects_violation() {
            let polys = [make_poly(&[1, -2, 3]), make_poly(&[0, 0, 6])];
            assert!(!utils::all_infty_norm_below(&polys, 5));
        }
    }

    /// Signing then verifying with the same key should succeed.
    #[test]
    fn sign_and_verify_roundtrip() {
        let KeyPair {
            public: pub_key,
            private: priv_key,
        } = keypair_fixture_1::<FieldElement>();
        let msg = b"hello, sign+verify!";
        let sig = priv_key.sign(msg).expect("signing should succeed");
        assert!(pub_key.verify(msg, &sig));
    }

    /// Verification must fail if the message changes after signing.
    #[test]
    fn verify_rejects_modified_message() {
        let KeyPair {
            public: pub_key,
            private: priv_key,
        } = keypair_fixture_1::<FieldElement>();
        let msg = b"immutable message";
        let sig = priv_key.sign(msg).expect("signing should succeed");

        // Verify against a different message
        let other = b"immutable message (edited)";
        assert!(
            !pub_key.verify(other, &sig),
            "verification must fail if the message changes"
        );
    }

    /// Tampering with z should cause signature verification to fail.
    #[test]
    fn verify_rejects_tampered_z() {
        let KeyPair {
            public: pub_key,
            private: priv_key,
        } = keypair_fixture_1::<FieldElement>();
        let msg = b"tamper z";
        let mut sig = priv_key.sign(msg).expect("signing should succeed");

        // Flip one coefficient in z[0]
        let mut plus_one = [0i64; N];
        plus_one[0] = 1;
        let add1: Polynomial<'static, FieldElement> = plus_one.into();
        sig.z[0] += &add1;

        assert!(
            !pub_key.verify(msg, &sig),
            "verification must fail when z is altered"
        );
    }

    /// Tampering with the challenge polynomial must invalidate the signature.
    #[test]
    fn verify_rejects_tampered_c() {
        let KeyPair {
            public: pub_key,
            private: priv_key,
        } = keypair_fixture_1::<FieldElement>();
        let msg = b"tamper c";
        let mut sig = priv_key.sign(msg).expect("signing should succeed");

        // Replace the challenge with a sparse poly that's *not* the derived one
        let mut c = [0i64; N];
        c[0] = 1;
        sig.c = c.into();

        assert!(
            !pub_key.verify(msg, &sig),
            "verification must fail when c is altered"
        );
    }

    /// Signatures must not validate under a different public key.
    #[test]
    fn verify_fails_with_wrong_public_key() {
        let KeyPair {
            public: pub_key1,
            private: priv_key1,
        } = keypair_fixture_1::<FieldElement>();
        let KeyPair {
            public: pub_key2,
            private: _,
        } = keypair_fixture_2::<FieldElement>();

        let msg = b"use the right key, please";
        let sig = priv_key1.sign(msg).expect("signing should succeed");

        assert!(pub_key1.verify(msg, &sig));
        assert!(
            !pub_key2.verify(msg, &sig),
            "verification must fail under a different public key"
        );
    }

    /// Verification should fail when checked against a different message.
    #[test]
    fn verify_rejects_wrong_message() {
        let KeyPair {
            public: pub_key,
            private: priv_key,
        } = keypair_fixture_1::<FieldElement>();
        let message = b"hello world";
        let signature = priv_key.sign(message).expect("signing should succeed");

        let wrong = b"hello wurld";
        assert!(
            !pub_key.verify(wrong, &signature),
            "verification must fail for a different message"
        );
    }

    /// Signing the same message twice should yield identical signatures.
    #[test]
    fn signatures_are_deterministic_for_same_message() {
        // With the current design (y derived from msg via SHAKE256), signing is deterministic.
        let KeyPair {
            public: _,
            private: priv_key1,
        } = keypair_fixture_1::<FieldElement>();
        let msg = b"deterministic";

        let sig1 = priv_key1.sign(msg).expect("signing should succeed");
        let sig2 = priv_key1.sign(msg).expect("signing should succeed");

        // Compare c and each z[i] explicitly to highlight determinism.
        assert_eq!(sig1.c, sig2.c);
        for i in 0..L {
            assert_eq!(sig1.z[i], sig2.z[i], "z[{}] differs", i);
        }
    }

    /// Signing should work for both empty and long messages.
    #[test]
    fn handles_empty_and_long_messages() {
        let KeyPair {
            public: pub_key1,
            private: priv_key1,
        } = keypair_fixture_1::<FieldElement>();

        // Empty message
        let empty = b"";
        let sig_empty = priv_key1.sign(empty).expect("signing should succeed");
        assert!(pub_key1.verify(empty, &sig_empty));

        // Long message
        let long = vec![0xABu8; 8192];
        let sig_long = priv_key1.sign(&long).expect("signing should succeed");
        assert!(pub_key1.verify(&long, &sig_long));
    }

    /// The infinity norm of z should always remain within the specified bound.
    #[test]
    fn z_infinity_norm_is_within_bound() {
        let KeyPair {
            public: _,
            private: priv_key1,
        } = keypair_fixture_1::<FieldElement>();
        let msg = b"bounds check";
        let sig = priv_key1.sign(msg).expect("signing should succeed");

        for (i, poly) in sig.z.iter().enumerate() {
            let norm = poly.norm_infinity() as i64;
            assert!(
                norm < GAMMA1 - BETA,
                "z[{}] ||.||_inf = {} exceeds GAMMA1 - BETA",
                i,
                norm
            );
        }
    }
}
