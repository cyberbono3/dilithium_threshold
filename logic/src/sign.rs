use std::convert::{From, Into};


use crate::hash::shake256;
use crate::matrix::{MatrixA, mat_vec_mul};
use crate::params::{ALPHA, BETA, ETA, GAMMA1, GAMMA2, K, L, N, Q, TAU};
use crate::poly::mod_q;
use math::{
    poly::Polynomial,
    traits::FiniteField,
};



use rand::RngCore;
use num_traits::Zero;


#[derive(Clone, Debug)]
pub struct Signature<'a, FF: FiniteField> {
    pub z: [Polynomial<'a, FF>; L],
    pub c: Polynomial<'a, FF>, // challenge polynomial with TAU non-zeros in {-1,0,1}
}

fn high_low_bits(x: i64) -> (i64, i64) {
    // Input x is taken modulo q
    //TODO mod_q does not needed, address it 
    let x0 = mod_q(x);
    // w1 = floor(x / ALPHA), w0 = x - w1*ALPHA adjusted to centered (-ALPHA/2, ALPHA/2]
    let w1 = x0 / ALPHA;
    let mut w0 = x0 - w1 * ALPHA;
    if w0 > ALPHA / 2 {
        w0 -= ALPHA;
    }
    (w1, w0)
}

fn poly_high_low<FF: FiniteField + From<FF> + From<i64>  >(p: &Polynomial<'_, FF>) -> (Polynomial<'static, FF>, Polynomial<'static, FF>) where i64: From<FF>  {
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

fn poly_high<FF: FiniteField + From<FF> + From<i64>>(p: &Polynomial<'_, FF>) -> Polynomial<'static, FF> where i64: From<FF> {
    poly_high_low(p).0
}

fn poly_low<FF: FiniteField + From<i64>>(p: &Polynomial<'_, FF>) -> Polynomial<'static, FF> where i64: From<FF>  {
    poly_high_low(p).1
}

/// Pack w1 (k polys of small integers) into bytes for hashing.
/// For ML-DSA-44, w1 coeffs are in 0..=43 (6 bits). We'll store as u16 for simplicity.
fn pack_w1_for_hash<FF: FiniteField + Into<u16> >(w1: &[Polynomial<'_, FF>; K]) -> Vec<u8> {
    let mut out = Vec::with_capacity(K * N * 2);
    for i in 0..K {
        for v in w1[i].coefficients() {
            let u: u16 = (*v).into();
            out.extend_from_slice(&u.to_le_bytes());
        }
    }
    out
}

/// Sample masking y in [-GAMMA1, GAMMA1] deterministically from a seed.
fn sample_y<FF: FiniteField + From<i64>>(seed: &[u8], ctr: u32) -> [Polynomial<'static, FF>; L] {
    let mut out = [Polynomial::zero(), Polynomial::zero(), Polynomial::zero(), Polynomial::zero()];
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
fn derive_challenge<FF: FiniteField + From<i64> >(m: &[u8], w1_pack: &[u8]) -> Polynomial<'static, FF> {
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

/// Sign according to the (uncompressed) template in Fig. 1 of the spec/paper.
/// z = y + c*s1; checks: ||z||_∞ < GAMMA1 - BETA,  ||LowBits(Ay - c*s2)||_∞ < GAMMA2 - BETA
pub fn sign<FF: FiniteField + Into<u16> + From<i64> >(
    sk_a: &MatrixA<'_, FF>,
    s1: &[Polynomial<'_, FF>; L],
    s2: &[Polynomial<'_, FF>; K],
    t: &[Polynomial<'_, FF>; K],
    msg: &[u8],
) -> Signature<'static, FF> where i64: From<FF>  {
    // Hedged/deterministic seed (simple approach: H(msg) as seed for y)
    let y_seed = shake256(32, msg);

    // Rejection loop
    for ctr in 0u32..10_000 {
        // bound to avoid infinite loop in a demo
        let y = sample_y(&y_seed, ctr);

        // w = A*y
        let w = mat_vec_mul(sk_a, &y);

        // w1 = HighBits(w, 2*GAMMA2)
        let w1 = w.clone().map(|p| poly_high(&p));
        let w1_pack = pack_w1_for_hash(&w1);

        // challenge c from H(M || w1)
        let c = derive_challenge(msg, &w1_pack);

        // z = y + c*s1 (component-wise in Rq)
        let mut z = [Polynomial::zero(), Polynomial::zero(), Polynomial::zero(), Polynomial::zero()];
        for j in 0..L {
            // TODO remove remove by immplementing new traits for Polynomial
            let cj_s1 = c.clone() * s1[j].clone();
            let mut z_j = y[j].clone();
            z_j += cj_s1.clone();
            z[j] = z_j;
        }

        // First reject check: ||z||_∞ < GAMMA1 - BETA
        let ok1 = z.iter().all(|p|(p.norm_infinity() as i64) < (GAMMA1 - BETA));

        // Second reject check: ||LowBits(Ay - c*s2, 2*GAMMA2)||_∞ < GAMMA2 - BETA
        let mut cs2 = [Polynomial::zero(), Polynomial::zero(), Polynomial::zero(), Polynomial::zero()];
        for i in 0..K {
            cs2[i] = c.clone() * s2[i].clone();
        }
        let mut ay_minus_cs2 = w.clone();
        for i in 0..K {
            // TODO implemenet SubAssign for 
            ay_minus_cs2[i] = ay_minus_cs2[i].clone() - cs2[i].clone();
        }
        let w0 = ay_minus_cs2.map(|p| poly_low(&p));
        let ok2 = w0.iter().all(|p| (p.norm_infinity() as i64)  < (GAMMA2 - BETA));

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
pub fn verify<FF: FiniteField + From<i64> + Into<u16> + From<FF> + Into<i64> + 'static >(
    pk_a: &MatrixA<'_, FF>,
    t: &[Polynomial<'_, FF>; K],
    msg: &[u8],
    sig: &Signature<'_, FF>,
) -> bool where i64: From<FF>{
    // ||z||_∞ < GAMMA1 - BETA
    if !sig.z.iter().all(|p| (p.norm_infinity() as i64)  < (GAMMA1 - BETA)) {
        return false;
    }

    // Compute Az
    let az = mat_vec_mul(pk_a, &sig.z);

    // Compute c*t (k polys), subtract: Az - c*t
    let mut c_t = [Polynomial::zero(), Polynomial::zero(), Polynomial::zero(), Polynomial::zero()];
    for i in 0..K {
        // TODO make trait for mul &
        c_t[i] = sig.c.clone() * t[i].clone();
    }
    let mut az_minus_ct = az.clone();
    for i in 0..K {
        az_minus_ct[i] = az_minus_ct[i].clone() - c_t[i].clone();
    }

    // w1' high bits
    let w1p = az_minus_ct.map(|p| poly_high(&p));
    let w1p_pack = super::sign::pack_w1_for_hash(&w1p);

    // recompute challenge
    let c2 = derive_challenge(msg, &w1p_pack);

    c2 == sig.c
}


// #[cfg(test)]
// mod tests {
//     use crate::keypair::keygen;
//     use crate::params::{BETA, GAMMA1, K, L};
//     use crate::poly::mod_q;
//     use crate::sign::{sign, verify};

//     #[test]
//     fn sign_round_trip_ok() {
//         let (pk, sk) = keygen();
//         let msg = b"round trip OK";
//         let sig = sign(&sk.a, &sk.s1, &sk.s2, &pk.t, msg);
//         assert!(verify(&pk.a, &pk.t, msg, &sig));
//     }

//     #[test]
//     fn sign_is_deterministic_given_keys_and_message() {
//         let (pk, sk) = keygen();
//         let msg = b"deterministic signature";
//         let s1 = sign(&sk.a, &sk.s1, &sk.s2, &pk.t, msg);
//         let s2 = sign(&sk.a, &sk.s1, &sk.s2, &pk.t, msg);

//         // Compare field-by-field (Signature doesn't derive PartialEq)
//         assert_eq!(s1.c, s2.c);
//         for j in 0..L {
//             assert_eq!(s1.z[j], s2.z[j], "z differs at index {}", j);
//         }
//     }

//     #[test]
//     fn verify_fails_if_message_changes() {
//         let (pk, sk) = keygen();
//         let sig = sign(&sk.a, &sk.s1, &sk.s2, &pk.t, b"hello");
//         assert!(!verify(&pk.a, &pk.t, b"hullo", &sig));
//     }

//     #[test]
//     fn verify_fails_if_z_is_tampered() {
//         let (pk, sk) = keygen();
//         let msg = b"tamper z";
//         let mut sig = sign(&sk.a, &sk.s1, &sk.s2, &pk.t, msg);

//         // Flip one coefficient but keep it in [0,q)
//         sig.z[0].c[0] = mod_q(sig.z[0].c[0] + 1);

//         assert!(!verify(&pk.a, &pk.t, msg, &sig));
//     }

//     #[test]
//     fn verify_fails_if_c_is_tampered() {
//         let (pk, sk) = keygen();
//         let msg = b"tamper c";
//         let mut sig = sign(&sk.a, &sk.s1, &sk.s2, &pk.t, msg);

//         // Find a nonzero location in c and flip its sign
//         if let Some(pos) = sig.c.c.iter().position(|&v| v != 0) {
//             sig.c.c[pos] = -sig.c.c[pos];
//         } else {
//             panic!("challenge has no nonzeros, unexpected for TAU > 0");
//         }
//         assert!(!verify(&pk.a, &pk.t, msg, &sig));
//     }

//     #[test]
//     fn verify_fails_with_wrong_public_components() {
//         let (pk1, sk1) = keygen();
//         let (pk2, _sk2) = keygen();
//         let msg = b"wrong pk";
//         let sig = sign(&sk1.a, &sk1.s1, &sk1.s2, &pk1.t, msg);

//         // Wrong A and wrong t should both fail
//         assert!(!verify(&pk2.a, &pk2.t, msg, &sig));
//         // Mixed components should also fail
//         assert!(!verify(&pk1.a, &pk2.t, msg, &sig));
//         assert!(!verify(&pk2.a, &pk1.t, msg, &sig));
//     }

//     #[test]
//     fn verify_rejects_if_z_norm_bound_is_violated() {
//         let (pk, sk) = keygen();
//         let msg = b"bound check";
//         let mut sig = sign(&sk.a, &sk.s1, &sk.s2, &pk.t, msg);

//         // Force a coefficient just outside the allowed infinity norm bound.
//         sig.z[0].c[0] = (GAMMA1 - BETA + 1) as i64;

//         assert!(!verify(&pk.a, &pk.t, msg, &sig));
//     }

//     #[test]
//     fn sign_verify_handles_empty_and_long_messages() {
//         let (pk, sk) = keygen();

//         // Empty message
//         let empty = b"";
//         let sig_empty = sign(&sk.a, &sk.s1, &sk.s2, &pk.t, empty);
//         assert!(verify(&pk.a, &pk.t, empty, &sig_empty));

//         // Long message
//         let long = vec![0xABu8; 8192];
//         let sig_long = sign(&sk.a, &sk.s1, &sk.s2, &pk.t, &long);
//         assert!(verify(&pk.a, &pk.t, &long, &sig_long));
//     }
// }
