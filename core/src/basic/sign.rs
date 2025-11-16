use crate::basic::keys::PrivateKey;
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
    use crate::basic::utils;
    use crate::basic::{KeyPair, KeypairSeeds, keygen_with_seeds};
    use crate::dilithium::params::{BETA, GAMMA1, L, N};
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
