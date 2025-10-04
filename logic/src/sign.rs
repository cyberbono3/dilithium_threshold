use core::convert::Into;
use std::convert::From;

use crate::hash::shake256;
use crate::keypair::{PrivateKey, PublicKey};
use crate::matrix::{MatrixA, mat_vec_mul};
use crate::params::{ALPHA, BETA, GAMMA1, GAMMA2, K, L, N, TAU};
use crate::utils::zero_polyvec;

use math::prelude::FieldElement;
use math::{poly::Polynomial, traits::FiniteField};

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

fn poly_high_low<FF: FiniteField + From<FF> + From<i64>>(
    p: &Polynomial<'_, FF>,
) -> (Polynomial<'static, FF>, Polynomial<'static, FF>)
where
    i64: From<FF>,
{
    let mut hi = [0i64; N];
    let mut lo = [0i64; N];
    for i in 0..N {
        // TODO add error handling
        let val = p.coefficients()[i];
        let (h, l) = high_low_bits(val.into());
        hi[i] = h;
        lo[i] = l;
    }
    (hi.into(), lo.into())
}

fn poly_high<FF: FiniteField + From<FF> + From<i64>>(
    p: &Polynomial<'_, FF>,
) -> Polynomial<'static, FF>
where
    i64: From<FF>,
{
    poly_high_low(p).0
}

fn poly_low<FF: FiniteField + From<i64>>(
    p: &Polynomial<'_, FF>,
) -> Polynomial<'static, FF>
where
    i64: From<FF>,
{
    poly_high_low(p).1
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
fn polyvec_add_scaled_in_place<FF: FiniteField + 'static, const LEN: usize>(
    dest: &mut [Polynomial<'_, FF>; LEN],
    base: &[Polynomial<'static, FF>; LEN],
    scale: &Polynomial<'static, FF>,
    mult: &[Polynomial<'static, FF>; LEN],
) {
    for i in 0..LEN {
        let mut tmp = base[i].clone();
        tmp += scale.clone() * mult[i].clone();
        dest[i] = tmp;
    }
}

#[inline]
fn polyvec_sub_scaled_in_place<FF: FiniteField + 'static, const LEN: usize>(
    dest: &mut [Polynomial<'_, FF>; LEN],
    base: &[Polynomial<'_, FF>; LEN],
    scale: &Polynomial<'_, FF>,
    mult: &[Polynomial<'_, FF>; LEN],
) {
    for i in 0..LEN {
        dest[i] = base[i].clone() - scale.clone() * mult[i].clone();
    }
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
) -> Signature<'static, FF>
where
    i64: From<FF>,
{
    let y_seed = shake256(32, msg);

    for ctr in 0u32..REJECTION_LIMIT {
        let y = sample_y(&y_seed, ctr);
        let w = mat_vec_mul(&priv_key.a, &y);

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

        if ok1 && ok2 {
            return Signature { z, c };
        }
    }

    panic!("sign: rejection sampling failed to converge (demo limit)");
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

    let az = mat_vec_mul(&pub_key.a, &sig.z);
    let mut az_minus_ct = zero_polyvec::<L, FF>();
    polyvec_sub_scaled_in_place::<FF, K>(
        &mut az_minus_ct,
        &az,
        &sig.c,
        &pub_key.t,
    );

    let w1p = az_minus_ct.map(|p| poly_high(&p));
    let w1p_pack = pack_w1_for_hash(&w1p);

    let c2 = derive_challenge(msg, &w1p_pack);
    c2 == sig.c
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

    #[test]
    fn sign_and_verify_roundtrip() {
        let (pub_key, priv_key) = keypair_fixture_1();
        let msg = b"hello, sign+verify!";
        let sig = super::sign::<FieldElement>(&priv_key, &pub_key.t, msg);
        assert!(super::verify::<FieldElement>(&pub_key, msg, &sig));
    }

    #[test]
    fn verify_rejects_modified_message() {
        let (pub_key, priv_key) = keypair_fixture_1();
        let msg = b"immutable message";
        let sig = super::sign::<FieldElement>(&priv_key, &pub_key.t, msg);

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
        let mut sig = super::sign::<FieldElement>(&priv_key, &pub_key.t, msg);

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
        let mut sig = super::sign::<FieldElement>(&priv_key, &pub_key.t, msg);

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
        let sig = super::sign::<FieldElement>(&priv_key1, &pub_key1.t, msg);

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

        let sig1 = super::sign::<FieldElement>(&priv_key1, &pub_key1.t, msg);
        let sig2 = super::sign::<FieldElement>(&priv_key1, &pub_key1.t, msg);

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
            super::sign::<FieldElement>(&priv_key1, &pub_key1.t, empty);
        assert!(super::verify::<FieldElement>(&pub_key1, empty, &sig_empty));

        // Long message
        let long = vec![0xABu8; 8192];
        let sig_long =
            super::sign::<FieldElement>(&priv_key1, &pub_key1.t, &long);
        assert!(super::verify::<FieldElement>(&pub_key1, &long, &sig_long));
    }

    #[test]
    fn z_infinity_norm_is_within_bound() {
        let (pub_key1, priv_key1) = keypair_fixture_1();
        let msg = b"bounds check";
        let sig = super::sign::<FieldElement>(&priv_key1, &pub_key1.t, msg);

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
