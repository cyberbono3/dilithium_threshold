/*
use crate::hash::shake256;
use crate::matrix::{MatrixA, mat_vec_mul};
use crate::params::{ALPHA, BETA, ETA, GAMMA1, GAMMA2, K, L, N, Q, TAU};
use crate::poly::{Poly, centered, mod_q};
use rand::RngCore;

#[derive(Clone, Debug)]
pub struct Signature {
    pub z: [Poly; L],
    pub c: Poly, // challenge polynomial with TAU non-zeros in {-1,0,1}
}

fn high_low_bits(x: i64) -> (i64, i64) {
    // Input x is taken modulo q
    let x0 = mod_q(x);
    // w1 = floor(x / ALPHA), w0 = x - w1*ALPHA adjusted to centered (-ALPHA/2, ALPHA/2]
    let w1 = x0 / ALPHA;
    let mut w0 = x0 - w1 * ALPHA;
    if w0 > ALPHA / 2 {
        w0 -= ALPHA;
    }
    (w1, w0)
}

fn poly_high_low(p: &Poly) -> (Poly, Poly) {
    let mut hi = [0i64; N];
    let mut lo = [0i64; N];
    for i in 0..N {
        let (h, l) = high_low_bits(p.c[i]);
        hi[i] = h;
        lo[i] = l;
    }
    (Poly { c: hi }, Poly { c: lo })
}

fn poly_high(p: &Poly) -> Poly {
    poly_high_low(p).0
}

fn poly_low(p: &Poly) -> Poly {
    poly_high_low(p).1
}

/// Pack w1 (k polys of small integers) into bytes for hashing.
/// For ML-DSA-44, w1 coeffs are in 0..=43 (6 bits). We'll store as u16 for simplicity.
fn pack_w1_for_hash(w1: &[Poly; K]) -> Vec<u8> {
    let mut out = Vec::with_capacity(K * N * 2);
    for i in 0..K {
        for &v in &w1[i].c {
            let u = v as u16;
            out.extend_from_slice(&u.to_le_bytes());
        }
    }
    out
}

/// Sample masking y in [-GAMMA1, GAMMA1] deterministically from a seed.
fn sample_y(seed: &[u8], ctr: u32) -> [Poly; L] {
    let mut out = [Poly::zero(), Poly::zero(), Poly::zero(), Poly::zero()];
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
        out[j] = Poly { c };
    }
    out
}

/// Build a sparse ternary polynomial c with TAU nonzeros in {-1, +1},
/// derived from H(M || w1_pack). Uses SHAKE256 as a PRF.
fn derive_challenge(m: &[u8], w1_pack: &[u8]) -> Poly {
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
    Poly { c }
}

/// Sign according to the (uncompressed) template in Fig. 1 of the spec/paper.
/// z = y + c*s1; checks: ||z||_∞ < GAMMA1 - BETA,  ||LowBits(Ay - c*s2)||_∞ < GAMMA2 - BETA
pub fn sign(
    sk_a: &MatrixA,
    s1: &[Poly; L],
    s2: &[Poly; K],
    t: &[Poly; K],
    msg: &[u8],
) -> Signature {
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
        let mut z = [Poly::zero(), Poly::zero(), Poly::zero(), Poly::zero()];
        for j in 0..L {
            let cj_s1 = c.mul(&s1[j]);
            let mut z_j = y[j].clone();
            z_j.add_assign(&cj_s1);
            z[j] = z_j;
        }

        // First reject check: ||z||_∞ < GAMMA1 - BETA
        let ok1 = z.iter().all(|p| p.norm_inf() < (GAMMA1 - BETA));

        // Second reject check: ||LowBits(Ay - c*s2, 2*GAMMA2)||_∞ < GAMMA2 - BETA
        let mut cs2 = [Poly::zero(), Poly::zero(), Poly::zero(), Poly::zero()];
        for i in 0..K {
            cs2[i] = c.mul(&s2[i]);
        }
        let mut ay_minus_cs2 = w.clone();
        for i in 0..K {
            ay_minus_cs2[i].sub_assign(&cs2[i]);
        }
        let w0 = ay_minus_cs2.map(|p| poly_low(&p));
        let ok2 = w0.iter().all(|p| p.norm_inf() < (GAMMA2 - BETA));

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
pub fn verify(
    pk_a: &MatrixA,
    t: &[Poly; K],
    msg: &[u8],
    sig: &Signature,
) -> bool {
    // ||z||_∞ < GAMMA1 - BETA
    if !sig.z.iter().all(|p| p.norm_inf() < (GAMMA1 - BETA)) {
        return false;
    }

    // Compute Az
    let az = mat_vec_mul(pk_a, &sig.z);

    // Compute c*t (k polys), subtract: Az - c*t
    let mut c_t = [Poly::zero(), Poly::zero(), Poly::zero(), Poly::zero()];
    for i in 0..K {
        c_t[i] = sig.c.mul(&t[i]);
    }
    let mut az_minus_ct = az.clone();
    for i in 0..K {
        az_minus_ct[i].sub_assign(&c_t[i]);
    }

    // w1' high bits
    let w1p = az_minus_ct.map(|p| poly_high(&p));
    let w1p_pack = super::sign::pack_w1_for_hash(&w1p);

    // recompute challenge
    let c2 = derive_challenge(msg, &w1p_pack);

    c2 == sig.c
}

*/

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
