use crate::dilithium::error::DilithiumError;
use crate::dilithium::params::{ALPHA, BETA, GAMMA1, GAMMA2, K, L, N};
use crate::dilithium::utils::zero_polyvec;
use crate::matrix::hash::shake256;
use crate::basic::keypair::{PrivateKey, PublicKey};

use math::prelude::FieldElement;
use math::{poly::Polynomial, traits::FiniteField};
use std::array::from_fn;
const REJECTION_LIMIT: u32 = 10000;

#[derive(Clone, Debug, PartialEq)]
pub struct Signature<'a, FF: FiniteField> {
    pub z: [Polynomial<'a, FF>; L],
    pub c: Polynomial<'a, FF>, // challenge polynomial with TAU non-zeros in {-1,0,1}
}

mod utils {
    use super::*;
    use crate::dilithium::utils::derive_challenge_polynomial;

    pub(super) fn poly_high_low<FF: FiniteField + From<i64>>(
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

    pub(super) fn poly_high<FF: FiniteField + From<i64>>(
        p: &Polynomial<'_, FF>,
    ) -> Polynomial<'static, FF>
    where
        i64: From<FF>,
    {
        let (high, _) = poly_high_low(p);
        high
    }

    pub(super) fn poly_low<FF: FiniteField + From<i64>>(
        p: &Polynomial<'_, FF>,
    ) -> Polynomial<'static, FF>
    where
        i64: From<FF>,
    {
        let (_, low) = poly_high_low(p);
        low
    }

    pub(super) fn pack_w1_for_hash<
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

    pub(super) fn sample_y<FF: FiniteField + From<i64>>(
        seed: &[u8],
        ctr: u32,
    ) -> [Polynomial<'static, FF>; L] {
        let mut out = zero_polyvec::<L, FF>();

        for (j, slot) in out.iter_mut().enumerate() {
            let mut inp = Vec::with_capacity(seed.len() + 6);
            inp.extend_from_slice(seed);
            inp.extend_from_slice(&(j as u16).to_le_bytes());
            inp.extend_from_slice(&ctr.to_le_bytes());
            let bytes = shake256(3 * N, &inp); // 3 bytes per coeff -> 24-bit
            let mut coeffs = [0i64; N];
            for (i, chunk) in bytes.chunks_exact(3).take(N).enumerate() {
                let v = (chunk[0] as i64)
                    | ((chunk[1] as i64) << 8)
                    | ((chunk[2] as i64) << 16);
                let modulus = (2 * GAMMA1 + 1) as i64;
                let centered = (v % modulus + modulus) % modulus - GAMMA1;
                coeffs[i] = centered;
            }
            *slot = coeffs.into();
        }
        out
    }

    pub(super) fn derive_challenge<FF: FiniteField + From<i64>>(
        message: &[u8],
        w1_pack: &[u8],
    ) -> Polynomial<'static, FF> {
        let mut seed = Vec::with_capacity(message.len() + w1_pack.len());
        seed.extend_from_slice(message);
        seed.extend_from_slice(w1_pack);
        derive_challenge_polynomial::<FF>(&seed)
    }

    #[inline]
    pub(super) fn polyvec_add_scaled<
        FF: FiniteField + 'static,
        const LEN: usize,
    >(
        base: &[Polynomial<'static, FF>; LEN],
        scale: &Polynomial<'_, FF>,
        mult: &[Polynomial<'static, FF>; LEN],
    ) -> [Polynomial<'static, FF>; LEN] {
        let mut dest = zero_polyvec::<LEN, FF>();
        for ((slot, base_poly), mult_poly) in
            dest.iter_mut().zip(base.iter()).zip(mult.iter())
        {
            slot.clone_from(base_poly);
            let scaled = mult_poly.clone() * scale.clone();
            *slot += scaled;
        }
        dest
    }

    #[inline]
    pub(super) fn polyvec_sub_scaled<
        FF: FiniteField + 'static,
        const LEN: usize,
    >(
        base: &[Polynomial<'static, FF>; LEN],
        scale: &Polynomial<'_, FF>,
        mult: &[Polynomial<'static, FF>; LEN],
    ) -> [Polynomial<'static, FF>; LEN] {
        let mut dest = zero_polyvec::<LEN, FF>();
        for ((slot, base_poly), mult_poly) in
            dest.iter_mut().zip(base.iter()).zip(mult.iter())
        {
            slot.clone_from(base_poly);
            let scaled = mult_poly.clone() * scale.clone();
            *slot += -scaled;
        }
        dest
    }

    #[inline]
    pub(super) fn all_infty_norm_below<FF: FiniteField, const LEN: usize>(
        polys: &[Polynomial<'_, FF>; LEN],
        bound: i64,
    ) -> bool
    where
        i64: From<FF>,
    {
        polys.iter().all(|p| (p.norm_infinity() as i64) < bound)
    }

    fn high_low_bits(x: i64) -> (i64, i64) {
        let w1 = x / ALPHA;
        let mut w0 = x - w1 * ALPHA;
        if w0 > ALPHA / 2 {
            w0 -= ALPHA;
        }
        (w1, w0)
    }

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
    all_infty_norm_below, derive_challenge, pack_w1_for_hash, poly_high,
    poly_low, polyvec_add_scaled, polyvec_sub_scaled, sample_y,
};

struct SigningEngine<'pk, 'msg, FF>
where
    FF: FiniteField + Into<[u8; FieldElement::BYTES]> + From<i64> + 'static,
    i64: From<FF>,
{
    priv_key: &'pk PrivateKey<'static, FF>,
    msg: &'msg [u8],
    y_seed: Vec<u8>,
}

impl<'pk, 'msg, FF> SigningEngine<'pk, 'msg, FF>
where
    FF: FiniteField + Into<[u8; FieldElement::BYTES]> + From<i64> + 'static,
    i64: From<FF>,
{
    fn new(priv_key: &'pk PrivateKey<'static, FF>, msg: &'msg [u8]) -> Self {
        Self {
            priv_key,
            msg,
            y_seed: shake256(32, msg),
        }
    }

    fn try_with_counter(&self, ctr: u32) -> Option<Signature<'static, FF>> {
        let y = sample_y::<FF>(&self.y_seed, ctr);
        let w = self.priv_key.a.mul(&y);
        let w1 = from_fn(|idx| poly_high(&w[idx]));
        let challenge = derive_challenge(self.msg, &pack_w1_for_hash(&w1));

        let z = polyvec_add_scaled::<FF, L>(&y, &challenge, &self.priv_key.s1);
        if !all_infty_norm_below::<FF, L>(&z, GAMMA1 - BETA) {
            return None;
        }

        let ay_minus_cs2 =
            polyvec_sub_scaled::<FF, K>(&w, &challenge, &self.priv_key.s2);
        let w0 = ay_minus_cs2.map(|p| poly_low(&p));
        if !all_infty_norm_below::<FF, K>(&w0, GAMMA2 - BETA) {
            return None;
        }

        Some(Signature { z, c: challenge })
    }
}

pub fn sign<FF: FiniteField + Into<[u8; FieldElement::BYTES]> + From<i64>>(
    priv_key: &PrivateKey<'static, FF>,
    _t: &[Polynomial<'_, FF>; K], // kept for API parity
    msg: &[u8],
) -> Result<Signature<'static, FF>, DilithiumError>
where
    i64: From<FF>,
{
    let engine = SigningEngine::new(priv_key, msg);

    (0u32..REJECTION_LIMIT)
        .find_map(|ctr| engine.try_with_counter(ctr))
        .ok_or(DilithiumError::SignatureGenerationFailed)
}

/// Verify (uncompressed template):
/// 1) compute w1' = HighBits(Az - c t, 2*GAMMA2)
/// 2) check ||z||_∞ < GAMMA1 - BETA
/// 3) check c == H(M || pack(w1'))
//

pub fn verify<
    FF: FiniteField + Into<[u8; FieldElement::BYTES]> + From<i64> + 'static,
>(
    pub_key: &PublicKey<'static, FF>,
    msg: &[u8],
    sig: &Signature<'_, FF>,
) -> bool
where
    i64: From<FF>,
{
    if !all_infty_norm_below::<FF, L>(&sig.z, GAMMA1 - BETA) {
        return false;
    }

    let az = pub_key.a.mul(&sig.z);
    let w1_prime = {
        let az_minus_ct = polyvec_sub_scaled::<FF, K>(&az, &sig.c, &pub_key.t);
        az_minus_ct.map(|poly| poly_high(&poly))
    };

    let derived = derive_challenge(msg, &pack_w1_for_hash(&w1_prime));
    derived == sig.c
}

#[cfg(test)]
mod tests {
    use super::*; // brings sign/verify + private helpers into scope
    use crate::dilithium::params::{BETA, GAMMA1, L, N};
    use crate::basic::keypair::*;
    use math::field_element::FieldElement;

    // Deterministic fixtures for reproducible tests
    fn keypair_fixture_1<FF: FiniteField + From<i64>>()
    -> (PublicKey<'static, FF>, PrivateKey<'static, FF>) {
        let rho = [0x42u8; 32];
        let s1 = [0x24u8; 32];
        let s2 = [0x18u8; 32];
        keygen_with_seeds::<FF>(rho, s1, s2)
    }

    fn keypair_fixture_2<FF: FiniteField + From<i64>>()
    -> (PublicKey<'static, FF>, PrivateKey<'static, FF>) {
        let rho = [0xA5u8; 32];
        let s1 = [0x5Au8; 32];
        let s2 = [0x33u8; 32];
        keygen_with_seeds::<FF>(rho, s1, s2)
    }

    mod signing_engine_tests {
        use super::*;

        #[test]
        fn engine_produces_signatures_within_limit() {
            let message = b"engine-test";
            let (pub_key, priv_key) = keypair_fixture_1::<FieldElement>();
            let engine = SigningEngine::new(&priv_key, message);

            let mut produced = None;
            for ctr in 0..REJECTION_LIMIT {
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
                verify::<FieldElement>(&pub_key, message, &signature),
                "engine output should verify against the public key"
            );
        }

        #[test]
        fn engine_is_deterministic_for_same_inputs() {
            let message = b"deterministic-engine";
            let seeds = ([0x42u8; 32], [0x24u8; 32], [0x18u8; 32]);

            let (_, priv_key_a) =
                keygen_with_seeds::<FieldElement>(seeds.0, seeds.1, seeds.2);
            let (_, priv_key_b) =
                keygen_with_seeds::<FieldElement>(seeds.0, seeds.1, seeds.2);

            let engine_a = SigningEngine::new(&priv_key_a, message);
            let engine_b = SigningEngine::new(&priv_key_b, message);

            let ctr = 5;
            assert_eq!(
                engine_a.try_with_counter(ctr),
                engine_b.try_with_counter(ctr),
                "same inputs must yield identical signatures"
            );
        }

        #[test]
        fn engine_rejects_when_mask_norm_exceeds_bounds() {
            let message = b"rejection-test";
            let (_, mut priv_key) = keypair_fixture_1::<FieldElement>();

            for poly in priv_key.s1.iter_mut() {
                let coeffs =
                    vec![FieldElement::from((GAMMA1 + BETA + 1) as i64); N];
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

        fn centered_value(fe: FieldElement) -> i64 {
            let mut v = i64::from(fe);
            let p = FieldElement::P as i64;
            if v > p / 2 {
                v -= p;
            }
            v
        }

        fn coeffs_within_bounds(
            polys: &[Polynomial<'_, FieldElement>],
        ) -> bool {
            polys
                .iter()
                .flat_map(|poly| poly.coefficients())
                .map(|&c| centered_value(c).abs())
                .all(|abs| abs <= GAMMA1)
        }

        #[test]
        fn deterministic_per_seed_and_counter() {
            let seed = [0x11u8; 32];
            let first = super::utils::sample_y::<FieldElement>(&seed, 7);
            let second = super::utils::sample_y::<FieldElement>(&seed, 7);
            assert_eq!(first, second);
        }

        #[test]
        fn different_counters_produce_distinct_vectors() {
            let seed = [0x77u8; 32];
            let first = super::utils::sample_y::<FieldElement>(&seed, 1);
            let second = super::utils::sample_y::<FieldElement>(&seed, 2);
            assert_ne!(first, second);
        }

        #[test]
        fn coefficients_are_centered() {
            let seed = [0xC3u8; 32];
            let polys = super::utils::sample_y::<FieldElement>(&seed, 5);
            assert!(coeffs_within_bounds(&polys));
        }
    }

    mod pack_w1_tests {
        use super::*;

        fn simple_poly(coeffs: &[i64]) -> Polynomial<'static, FieldElement> {
            Polynomial::from(
                coeffs
                    .iter()
                    .map(|&c| FieldElement::from(c))
                    .collect::<Vec<_>>(),
            )
        }

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

        #[test]
        fn deterministic_for_same_inputs() {
            let msg = b"challenge";
            let hash = vec![0x42; 32];
            let c1 = super::utils::derive_challenge::<FieldElement>(msg, &hash);
            let c2 = super::utils::derive_challenge::<FieldElement>(msg, &hash);
            assert_eq!(c1, c2);
        }

        #[test]
        fn changing_message_changes_challenge() {
            let hash = vec![0x77; 64];
            let c1 =
                super::utils::derive_challenge::<FieldElement>(b"m1", &hash);
            let c2 =
                super::utils::derive_challenge::<FieldElement>(b"m2", &hash);
            assert_ne!(c1, c2);
        }

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

        fn make_poly(coeffs: &[i64]) -> Polynomial<'static, FieldElement> {
            Polynomial::from(
                coeffs
                    .iter()
                    .map(|&c| FieldElement::from(c))
                    .collect::<Vec<_>>(),
            )
        }

        #[test]
        fn detects_within_bound() {
            let polys = [make_poly(&[1, -2, 3]), make_poly(&[0, 0, 4])];
            assert!(super::utils::all_infty_norm_below(&polys, 5));
        }

        #[test]
        fn detects_violation() {
            let polys = [make_poly(&[1, -2, 3]), make_poly(&[0, 0, 6])];
            assert!(!super::utils::all_infty_norm_below(&polys, 5));
        }
    }

    #[test]
    fn sign_and_verify_roundtrip() {
        let (pub_key, priv_key) = keypair_fixture_1();
        let msg = b"hello, sign+verify!";
        let sig = super::sign::<FieldElement>(&priv_key, &pub_key.t, msg)
            .expect("signing should succeed");
        assert!(super::verify::<FieldElement>(&pub_key, msg, &sig));
    }

    #[test]
    fn verify_rejects_modified_message() {
        let (pub_key, priv_key) = keypair_fixture_1();
        let msg = b"immutable message";
        let sig = super::sign::<FieldElement>(&priv_key, &pub_key.t, msg)
            .expect("signing should succeed");

        // Verify against a different message
        let other = b"immutable message (edited)";
        assert!(
            !super::verify::<FieldElement>(&pub_key, other, &sig),
            "verification must fail if the message changes"
        );
    }

    #[test]
    fn verify_rejects_tampered_z() {
        let (pub_key, priv_key) = keypair_fixture_1();
        let msg = b"tamper z";
        let mut sig = super::sign::<FieldElement>(&priv_key, &pub_key.t, msg)
            .expect("signing should succeed");

        // Flip one coefficient in z[0]
        let mut plus_one = [0i64; N];
        plus_one[0] = 1;
        let add1: Polynomial<'static, FieldElement> = plus_one.into();
        sig.z[0] = sig.z[0].clone() + add1;

        assert!(
            !super::verify::<FieldElement>(&pub_key, msg, &sig),
            "verification must fail when z is altered"
        );
    }

    #[test]
    fn verify_rejects_tampered_c() {
        let (pub_key, priv_key) = keypair_fixture_1();
        let msg = b"tamper c";
        let mut sig = super::sign::<FieldElement>(&priv_key, &pub_key.t, msg)
            .expect("signing should succeed");

        // Replace the challenge with a sparse poly that's *not* the derived one
        let mut c = [0i64; N];
        c[0] = 1;
        sig.c = c.into();

        assert!(
            !super::verify::<FieldElement>(&pub_key, msg, &sig),
            "verification must fail when c is altered"
        );
    }

    #[test]
    fn verify_fails_with_wrong_public_key() {
        let (pub_key1, priv_key1) = keypair_fixture_1();
        let (pub_key2, _) = keypair_fixture_2();

        let msg = b"use the right key, please";
        let sig = super::sign::<FieldElement>(&priv_key1, &pub_key1.t, msg)
            .expect("signing should succeed");

        assert!(super::verify::<FieldElement>(&pub_key1, msg, &sig));
        assert!(
            !super::verify::<FieldElement>(&pub_key2, msg, &sig),
            "verification must fail under a different public key"
        );
    }

    #[test]
    fn verify_rejects_wrong_message() {
        let (pub_key, priv_key) = keypair_fixture_1();
        let message = b"hello world";
        let signature =
            super::sign::<FieldElement>(&priv_key, &pub_key.t, message)
                .expect("signing should succeed");

        let wrong = b"hello wurld";
        assert!(
            !super::verify::<FieldElement>(&pub_key, wrong, &signature),
            "verification must fail for a different message"
        );
    }

    #[test]
    fn signatures_are_deterministic_for_same_message() {
        // With the current design (y derived from msg via SHAKE256), signing is deterministic.
        let (pub_key1, priv_key1) = keypair_fixture_1();
        let msg = b"deterministic";

        let sig1 = super::sign::<FieldElement>(&priv_key1, &pub_key1.t, msg)
            .expect("signing should succeed");
        let sig2 = super::sign::<FieldElement>(&priv_key1, &pub_key1.t, msg)
            .expect("signing should succeed");

        // Compare c and each z[i] (Signature doesn't implement PartialEq)
        assert_eq!(sig1.c, sig2.c);
        for i in 0..L {
            assert_eq!(sig1.z[i], sig2.z[i], "z[{}] differs", i);
        }
    }

    #[test]
    fn handles_empty_and_long_messages() {
        let (pub_key1, priv_key1) = keypair_fixture_1();

        // Empty message
        let empty = b"";
        let sig_empty =
            super::sign::<FieldElement>(&priv_key1, &pub_key1.t, empty)
                .expect("signing should succeed");
        assert!(super::verify::<FieldElement>(&pub_key1, empty, &sig_empty));

        // Long message
        let long = vec![0xABu8; 8192];
        let sig_long =
            super::sign::<FieldElement>(&priv_key1, &pub_key1.t, &long)
                .expect("signing should succeed");
        assert!(super::verify::<FieldElement>(&pub_key1, &long, &sig_long));
    }

    #[test]
    fn z_infinity_norm_is_within_bound() {
        let (pub_key1, priv_key1) = keypair_fixture_1();
        let msg = b"bounds check";
        let sig = super::sign::<FieldElement>(&priv_key1, &pub_key1.t, msg)
            .expect("signing should succeed");

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
