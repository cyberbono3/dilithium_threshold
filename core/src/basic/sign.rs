use crate::basic::keypair::PrivateKey;
use crate::dilithium::params::{ALPHA, BETA, GAMMA1, GAMMA2, K, L, N, Q};
use crate::matrix::{MatrixMulExt, hash::shake256};

use math::prelude::FieldElement;
use math::{poly::Polynomial, traits::FiniteField};
use std::array::from_fn;
const HINT_MODULUS: i64 = (Q - 1) / ALPHA;

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

pub(crate) mod utils {
    use super::*;
    use crate::dilithium::utils::{
        derive_challenge_polynomial, shake256_squeezed, zero_polyvec,
    };

    /// Split a polynomial into its high and low components.
    pub(crate) fn poly_high_low<FF: FiniteField + From<i64>>(
        p: &Polynomial<'_, FF>,
    ) -> (Polynomial<'static, FF>, Polynomial<'static, FF>)
    where
        i64: From<FF>,
    {
        let (hi, lo): (Vec<_>, Vec<_>) = padded_coefficients(p)
            .map(|coeff| high_low_bits(i64::from(coeff)))
            .unzip();
        (hi.into(), lo.into())
    }

    /// Extract only the high component of a polynomial used for hashing.
    pub(crate) fn poly_high<FF: FiniteField + From<i64>>(
        p: &Polynomial<'_, FF>,
    ) -> Polynomial<'static, FF>
    where
        i64: From<FF>,
    {
        let (high, _) = poly_high_low(p);
        high
    }

    /// Extract only the low component of a polynomial.
    pub(crate) fn poly_low<FF: FiniteField + From<i64>>(
        p: &Polynomial<'_, FF>,
    ) -> Polynomial<'static, FF>
    where
        i64: From<FF>,
    {
        let (_, low) = poly_high_low(p);
        low
    }

    /// Serialize the w1 vector into bytes that feed the challenge hash.
    pub(crate) fn pack_w1_for_hash<
        FF: FiniteField + Into<[u8; FieldElement::BYTES]>,
    >(
        w1: &[Polynomial<'_, FF>; K],
    ) -> Vec<u8> {
        let mut out = Vec::with_capacity(K * N * 2);
        for poly in w1 {
            for v in poly.coefficients() {
                let arr: [u8; FieldElement::BYTES] = (*v).into();
                out.extend_from_slice(&arr);
            }
        }
        out
    }

    /// Sample the y polynomial vector for a given message seed and counter.
    pub(crate) fn sample_y<FF: FiniteField + From<i64>>(
        seed: &[u8],
        ctr: u32,
    ) -> [Polynomial<'static, FF>; L] {
        let mut out = zero_polyvec::<L, FF>();

        for (j, slot) in out.iter_mut().enumerate() {
            let idx_bytes = (j as u16).to_le_bytes();
            let ctr_bytes = ctr.to_le_bytes();
            let bytes =
                shake256_squeezed(seed, &[&idx_bytes, &ctr_bytes], 3 * N);
            let mut coeffs = [0i64; N];
            for (i, chunk) in bytes.chunks_exact(3).take(N).enumerate() {
                let v = (chunk[0] as i64)
                    | ((chunk[1] as i64) << 8)
                    | ((chunk[2] as i64) << 16);
                let modulus = 2 * GAMMA1 + 1;
                let centered = (v % modulus + modulus) % modulus - GAMMA1;
                coeffs[i] = centered;
            }
            *slot = coeffs.into();
        }
        out
    }

    /// Derive the sparse challenge polynomial from message and w1 bytes.
    pub(crate) fn derive_challenge<FF: FiniteField + From<i64>>(
        message: &[u8],
        w1_pack: &[u8],
    ) -> Polynomial<'static, FF> {
        let mut seed = Vec::with_capacity(message.len() + w1_pack.len());
        seed.extend_from_slice(message);
        seed.extend_from_slice(w1_pack);
        derive_challenge_polynomial::<FF>(&seed)
    }

    #[inline]
    fn polyvec_add_scaled_with_sign<
        FF: FiniteField + 'static,
        const LEN: usize,
    >(
        base: &[Polynomial<'static, FF>; LEN],
        scale: &Polynomial<'_, FF>,
        mult: &[Polynomial<'static, FF>; LEN],
        negate: bool,
    ) -> [Polynomial<'static, FF>; LEN] {
        let mut dest = zero_polyvec::<LEN, FF>();
        for ((slot, base_poly), mult_poly) in
            dest.iter_mut().zip(base).zip(mult)
        {
            slot.clone_from(base_poly);
            let mut scaled = mult_poly * scale;
            if negate {
                scaled = -scaled;
            }
            *slot += scaled;
        }
        dest
    }

    /// Compute element-wise addition of `base` and `scale * mult`.
    #[inline]
    pub(crate) fn polyvec_add_scaled<
        FF: FiniteField + 'static,
        const LEN: usize,
    >(
        base: &[Polynomial<'static, FF>; LEN],
        scale: &Polynomial<'_, FF>,
        mult: &[Polynomial<'static, FF>; LEN],
    ) -> [Polynomial<'static, FF>; LEN] {
        polyvec_add_scaled_with_sign(base, scale, mult, false)
    }

    /// Compute element-wise subtraction of `scale * mult` from `base`.
    #[inline]
    pub(crate) fn polyvec_sub_scaled<
        FF: FiniteField + 'static,
        const LEN: usize,
    >(
        base: &[Polynomial<'static, FF>; LEN],
        scale: &Polynomial<'_, FF>,
        mult: &[Polynomial<'static, FF>; LEN],
    ) -> [Polynomial<'static, FF>; LEN] {
        polyvec_add_scaled_with_sign(base, scale, mult, true)
    }

    /// Fetch the coefficient at `idx`, returning zero when out of bounds.
    fn coeff_at_or_zero<FF: FiniteField>(
        poly: &Polynomial<'_, FF>,
        idx: usize,
    ) -> FF {
        poly.coefficients().get(idx).copied().unwrap_or(FF::ZERO)
    }

    fn make_hint_poly<FF>(
        target_high: &Polynomial<'_, FF>,
        base: &Polynomial<'_, FF>,
    ) -> Polynomial<'static, FF>
    where
        FF: FiniteField + From<i64>,
        i64: From<FF>,
    {
        let mut coeffs = vec![FF::ZERO; N];
        for (idx, coeff) in coeffs.iter_mut().enumerate() {
            let target_hi = i64::from(coeff_at_or_zero(target_high, idx));
            let base_value = i64::from(coeff_at_or_zero(base, idx));
            let (base_hi, _) = high_low_bits(base_value);
            if target_hi != base_hi {
                *coeff = FF::ONE;
            }
        }
        coeffs.into()
    }

    /// Generate per-coefficient hints enabling reconstruction of the reference high bits.
    pub(crate) fn make_hints<'a, FF>(
        target_high: &[Polynomial<'a, FF>; K],
        base: &[Polynomial<'a, FF>; K],
    ) -> [Polynomial<'static, FF>; K]
    where
        FF: FiniteField + From<i64>,
        i64: From<FF>,
    {
        std::array::from_fn(|idx| make_hint_poly(&target_high[idx], &base[idx]))
    }

    fn use_hint_poly<FF>(
        hint: &Polynomial<'_, FF>,
        base: &Polynomial<'_, FF>,
    ) -> Polynomial<'static, FF>
    where
        FF: FiniteField + From<i64>,
        i64: From<FF>,
    {
        let mut coeffs = vec![FF::ZERO; N];
        for (idx, coeff) in coeffs.iter_mut().enumerate() {
            let hint_bit = !coeff_at_or_zero(hint, idx).is_zero();
            let base_value = i64::from(coeff_at_or_zero(base, idx));
            let (mut hi, lo) = high_low_bits(base_value);
            if hint_bit {
                if lo > 0 {
                    hi += 1;
                } else {
                    hi -= 1;
                }
            }
            let modulus = super::HINT_MODULUS;
            hi = ((hi % modulus) + modulus) % modulus;
            *coeff = FF::from(hi);
        }
        coeffs.into()
    }

    /// Apply stored hints to recover the original high-bit decomposition.
    pub(crate) fn use_hints<'hint, 'poly, FF>(
        hints: &[Polynomial<'hint, FF>; K],
        base: &[Polynomial<'poly, FF>; K],
    ) -> [Polynomial<'static, FF>; K]
    where
        FF: FiniteField + From<i64>,
        i64: From<FF>,
    {
        std::array::from_fn(|idx| use_hint_poly(&hints[idx], &base[idx]))
    }

    /// Check whether each polynomial's infinity norm stays below `bound`.
    #[inline]
    pub(crate) fn all_infty_norm_below<FF: FiniteField, const LEN: usize>(
        polys: &[Polynomial<'_, FF>; LEN],
        bound: i64,
    ) -> bool
    where
        i64: From<FF>,
    {
        polys.iter().all(|p| (p.norm_infinity() as i64) < bound)
    }

    /// Split a coefficient into high and low parts relative to ALPHA.
    fn high_low_bits(x: i64) -> (i64, i64) {
        let w1 = x / ALPHA;
        let mut w0 = x - w1 * ALPHA;
        if w0 > ALPHA / 2 {
            w0 -= ALPHA;
        }
        (w1, w0)
    }

    /// Iterate coefficients padded with zeros up to length N.
    fn padded_coefficients<'a, FF: FiniteField>(
        poly: &'a Polynomial<'a, FF>,
    ) -> impl Iterator<Item = FF> + 'a {
        poly.coefficients()
            .iter()
            .copied()
            .chain(std::iter::repeat(FF::ZERO))
            .take(N)
    }
}

use utils::{
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
            let first = super::utils::sample_y::<FieldElement>(&seed, 7);
            let second = super::utils::sample_y::<FieldElement>(&seed, 7);
            assert_eq!(first, second);
        }

        /// Different counters must yield distinct sampled vectors.
        #[test]
        fn different_counters_produce_distinct_vectors() {
            let seed = [0x77u8; 32];
            let first = super::utils::sample_y::<FieldElement>(&seed, 1);
            let second = super::utils::sample_y::<FieldElement>(&seed, 2);
            assert_ne!(first, second);
        }

        /// All coefficients should remain centered within the permitted bound.
        #[test]
        fn coefficients_are_centered() {
            let seed = [0xC3u8; 32];
            let polys = super::utils::sample_y::<FieldElement>(&seed, 5);
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

            let first = super::utils::pack_w1_for_hash(&polys);
            let second = super::utils::pack_w1_for_hash(&polys);
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

            let hash_a = super::utils::pack_w1_for_hash(&polys_a);
            let hash_b = super::utils::pack_w1_for_hash(&polys_b);
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
            let c1 = super::utils::derive_challenge::<FieldElement>(msg, &hash);
            let c2 = super::utils::derive_challenge::<FieldElement>(msg, &hash);
            assert_eq!(c1, c2);
        }

        /// Different messages should produce different challenge polynomials.
        #[test]
        fn changing_message_changes_challenge() {
            let hash = vec![0x77; 64];
            let c1 =
                super::utils::derive_challenge::<FieldElement>(b"m1", &hash);
            let c2 =
                super::utils::derive_challenge::<FieldElement>(b"m2", &hash);
            assert_ne!(c1, c2);
        }

        /// The derived challenge must contain exactly TAU non-zero coefficients.
        #[test]
        fn challenge_has_tau_non_zero_entries() {
            let msg = b"nonzero-count";
            let hash = vec![0xAB; 128];
            let challenge =
                super::utils::derive_challenge::<FieldElement>(msg, &hash);
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

            let result = super::utils::polyvec_add_scaled::<FieldElement, 2>(
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

            let result = super::utils::polyvec_sub_scaled::<FieldElement, 2>(
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
            assert!(super::utils::all_infty_norm_below(&polys, 5));
        }

        /// Detect cases where the infinity norm exceeds the bound.
        #[test]
        fn detects_violation() {
            let polys = [make_poly(&[1, -2, 3]), make_poly(&[0, 0, 6])];
            assert!(!super::utils::all_infty_norm_below(&polys, 5));
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
