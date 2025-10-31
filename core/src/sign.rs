use crate::error::ThresholdError;
use crate::hash::shake256;
use crate::keypair::{PrivateKey, PublicKey};
use crate::params::{ALPHA, BETA, GAMMA1, GAMMA2, K, L, N, TAU};
use crate::utils::zero_polyvec;

use math::prelude::FieldElement;
use math::{poly::Polynomial, traits::FiniteField};
use std::ops::Mul;

const REJECTION_LIMIT: u32 = 10000;

#[derive(Clone, Debug)]
pub struct Signature<'a, FF: FiniteField> {
    pub z: [Polynomial<'a, FF>; L],
    pub c: Polynomial<'a, FF>, // challenge polynomial with TAU non-zeros in {-1,0,1}
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

fn poly_high_low<FF: FiniteField + From<i64>>(
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

fn poly_high<FF: FiniteField + From<i64>>(
    p: &Polynomial<'_, FF>,
) -> Polynomial<'static, FF>
where
    i64: From<FF>,
{
    let (high, _) = poly_high_low(p);
    high
}

fn poly_low<FF: FiniteField + From<i64>>(
    p: &Polynomial<'_, FF>,
) -> Polynomial<'static, FF>
where
    i64: From<FF>,
{
    let (_, low) = poly_high_low(p);
    low
}

/// Pack w1 (k polys of small integers) into bytes for hashing.
/// For ML-DSA-44, w1 coeffs are in 0..=43 (6 bits). We'll store as u16 for simplicity.
fn pack_w1_for_hash<FF: FiniteField + Into<[u8; FieldElement::BYTES]>>(
    w1: &[Polynomial<'_, FF>; K],
) -> Vec<u8> {
    let mut out = Vec::with_capacity(K * N * 2);
    for i in 0..K {
        for v in w1[i].coefficients() {
            //let u: u16 = (*v).try_into()?;
            let arr: [u8; FieldElement::BYTES] = (*v).into();
            out.extend_from_slice(&arr);
        }
    }
    out
}

/// Sample masking y in [-GAMMA1, GAMMA1] deterministically from a seed.
fn sample_y<FF: FiniteField + From<i64>>(
    seed: &[u8],
    ctr: u32,
) -> [Polynomial<'static, FF>; L] {
    let mut out = zero_polyvec::<L, FF>();

    for j in 0..L {
        let mut inp = Vec::new();
        inp.extend_from_slice(seed);
        inp.extend_from_slice(&(j as u16).to_le_bytes());
        inp.extend_from_slice(&ctr.to_le_bytes());
        let bytes = shake256(3 * N, &inp); // 3 bytes per coeff -> 24-bit
        let mut c = [0i64; N];
        for i in 0..N {
            let b = &bytes[3 * i..3 * i + 3];
            let v = ((b[0] as u32)
                | ((b[1] as u32) << 8)
                | ((b[2] as u32) << 16)) as i64;
            // Map to [-GAMMA1, GAMMA1)
            let m = (2 * GAMMA1 + 1) as i64;
            let r = (v % m + m) % m;
            c[i] = r - GAMMA1; // centered
        }
        out[j] = c.into();
    }
    out
}

/// Build a sparse ternary polynomial c with TAU nonzeros in {-1, +1},
/// derived from H(M || w1_pack). Uses SHAKE256 as a PRF.
fn derive_challenge<FF: FiniteField + From<i64>>(
    m: &[u8],
    w1_pack: &[u8],
) -> Polynomial<'static, FF> {
    let mut seed = Vec::with_capacity(m.len() + w1_pack.len());
    seed.extend_from_slice(m);
    seed.extend_from_slice(w1_pack);
    let mut stream = shake256(4 * TAU + 1024, &seed); // extra bytes for indices

    let mut used = [false; N];
    let mut c = [0i64; N];
    let mut filled = 0usize;
    let mut idx = 0usize;

    while filled < TAU {
        if idx + 3 >= stream.len() {
            let more = shake256(1024, &stream);
            stream.extend_from_slice(&more);
        }
        // 2 bytes for index, 1 byte for sign
        let j = u16::from_le_bytes([stream[idx], stream[idx + 1]]) as usize % N;
        let sbit = stream[idx + 2] & 1;
        idx += 3;

        if !used[j] {
            used[j] = true;
            c[j] = if sbit == 1 { 1 } else { -1 };
            filled += 1;
        }
    }
    c.into()
}

#[inline]
fn polyvec_add_scaled_in_place<
    'a,
    FF: FiniteField + 'static,
    const LEN: usize,
>(
    dest: &mut [Polynomial<'a, FF>; LEN],
    base: &[Polynomial<'a, FF>; LEN],
    scale: &Polynomial<'a, FF>,
    mult: &[Polynomial<'a, FF>; LEN],
) {
    dest.iter_mut().zip(base.iter()).zip(mult.iter()).for_each(
        |((slot, base_poly), mult_poly)| {
            // slot = base + scale * mult
            slot.clone_from(base_poly);
            *slot += scale.clone() * mult_poly.clone();
        },
    );
}

#[inline]
fn polyvec_sub_scaled_in_place<
    'a,
    FF: FiniteField + 'static,
    const LEN: usize,
>(
    dest: &mut [Polynomial<'a, FF>; LEN],
    base: &[Polynomial<'a, FF>; LEN],
    scale: &Polynomial<'a, FF>,
    mult: &[Polynomial<'a, FF>; LEN],
) {
    let neg_scale = -scale.clone();
    polyvec_add_scaled_in_place(dest, base, &neg_scale, mult);
}

#[inline]
fn all_infty_norm_below<FF: FiniteField, const LEN: usize>(
    polys: &[Polynomial<'_, FF>; LEN],
    bound: i64,
) -> bool
where
    i64: From<FF>,
{
    polys.iter().all(|p| (p.norm_infinity() as i64) < bound)
}

pub fn sign<FF: FiniteField + Into<[u8; FieldElement::BYTES]> + From<i64>>(
    priv_key: &PrivateKey<'static, FF>,
    _t: &[Polynomial<'_, FF>; K], // kept for API parity
    msg: &[u8],
) -> Result<Signature<'static, FF>, ThresholdError>
where
    i64: From<FF>,
{
    let y_seed = shake256(32, msg);

    (0u32..REJECTION_LIMIT)
        .find_map(|ctr| {
            let y = sample_y(&y_seed, ctr);
            let w = priv_key.a.mul(&y);

            let w1 = w.clone().map(|p| poly_high(&p));
            let w1_pack = pack_w1_for_hash(&w1);
            let c = derive_challenge(msg, &w1_pack);

            let mut z = zero_polyvec::<L, FF>();
            polyvec_add_scaled_in_place::<FF, L>(&mut z, &y, &c, &priv_key.s1);

            let ok1 = all_infty_norm_below::<FF, L>(&z, GAMMA1 - BETA);

            let mut ay_minus_cs2 = zero_polyvec::<L, FF>();
            polyvec_sub_scaled_in_place::<FF, K>(
                &mut ay_minus_cs2,
                &w,
                &c,
                &priv_key.s2,
            );
            let w0 = ay_minus_cs2.map(|p| poly_low(&p));
            let ok2 = all_infty_norm_below::<FF, K>(&w0, GAMMA2 - BETA);

            ok1.then_some((z, c)).filter(|_| ok2)
        })
        .map(|(z, c)| Signature { z, c })
        .ok_or(ThresholdError::SignatureGenerationFailed)
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
    all_infty_norm_below::<FF, L>(&sig.z, GAMMA1 - BETA)
        && derive_challenge(
            msg,
            &pack_w1_for_hash(&{
                let az = pub_key.a.mul(&sig.z);
                let mut az_minus_ct = zero_polyvec::<L, FF>();
                polyvec_sub_scaled_in_place::<FF, K>(
                    &mut az_minus_ct,
                    &az,
                    &sig.c,
                    &pub_key.t,
                );
                az_minus_ct.map(|p| poly_high(&p))
            }),
        ) == sig.c
}

#[cfg(test)]
mod tests {
    use super::*; // brings sign/verify + private helpers into scope
    use crate::keypair::*;
    use crate::params::{BETA, GAMMA1, L, N};
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
            let first = super::sample_y::<FieldElement>(&seed, 7);
            let second = super::sample_y::<FieldElement>(&seed, 7);
            assert_eq!(first, second);
        }

        #[test]
        fn different_counters_produce_distinct_vectors() {
            let seed = [0x77u8; 32];
            let first = super::sample_y::<FieldElement>(&seed, 1);
            let second = super::sample_y::<FieldElement>(&seed, 2);
            assert_ne!(first, second);
        }

        #[test]
        fn coefficients_are_centered() {
            let seed = [0xC3u8; 32];
            let polys = super::sample_y::<FieldElement>(&seed, 5);
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

            let first = super::pack_w1_for_hash(&polys);
            let second = super::pack_w1_for_hash(&polys);
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

            let hash_a = super::pack_w1_for_hash(&polys_a);
            let hash_b = super::pack_w1_for_hash(&polys_b);
            assert_ne!(hash_a, hash_b);
        }
    }

    mod derive_challenge_tests {
        use super::*;
        use num_traits::ConstZero;

        #[test]
        fn deterministic_for_same_inputs() {
            let msg = b"challenge";
            let hash = vec![0x42; 32];
            let c1 = super::derive_challenge::<FieldElement>(msg, &hash);
            let c2 = super::derive_challenge::<FieldElement>(msg, &hash);
            assert_eq!(c1, c2);
        }

        #[test]
        fn changing_message_changes_challenge() {
            let hash = vec![0x77; 64];
            let c1 = super::derive_challenge::<FieldElement>(b"m1", &hash);
            let c2 = super::derive_challenge::<FieldElement>(b"m2", &hash);
            assert_ne!(c1, c2);
        }

        #[test]
        fn challenge_has_tau_non_zero_entries() {
            let msg = b"nonzero-count";
            let hash = vec![0xAB; 128];
            let challenge = super::derive_challenge::<FieldElement>(msg, &hash);
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
            assert!(super::all_infty_norm_below(&polys, 5));
        }

        #[test]
        fn detects_violation() {
            let polys = [make_poly(&[1, -2, 3]), make_poly(&[0, 0, 6])];
            assert!(!super::all_infty_norm_below(&polys, 5));
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
