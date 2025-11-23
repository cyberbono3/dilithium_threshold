//! Correctness tests based on Section 4.2 of the CRYSTALS–Dilithium paper https://repository.ubn.ru.nl/bitstream/handle/2066/191703/191703.pdf
//!
//!
//! 1. Basic correctness: for all (pk, sk) ← KeyGen and messages m,
//!    Verify(pk, m, Sign(sk, m)) = 1.
//! 2. Equation (1) from Section 4.2:
//!    w - c s2 = A z - c t
//!    where:
//!    w = A y, z = y + c s1, t = A s1 + s2.
//!
//! The second test reconstructs w and both sides of (1) from the
//! secret key and the final signature and checks they are identical.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha20Rng;

use dilithium_core::basic::keys::KeypairSeeds;
use dilithium_core::basic::sign::DilithiumSignature;
use dilithium_core::basic::{KeyPair, PublicKey, keygen_with_seeds};
use dilithium_core::dilithium::params::{K, L};
use dilithium_core::matrix::{MatrixMulExt, expand_a_from_rho};
use math::{Matrix, field_element::FieldElement, poly::Polynomial};

fn generate_keypair(rng: &mut ChaCha20Rng) -> KeyPair<'static, FieldElement> {
    let mut rho = [0u8; 32];
    let mut s1 = [0u8; 32];
    let mut s2 = [0u8; 32];
    rng.fill_bytes(&mut rho);
    rng.fill_bytes(&mut s1);
    rng.fill_bytes(&mut s2);

    keygen_with_seeds::<FieldElement>(KeypairSeeds::new(rho, s1, s2))
        .expect("key generation should succeed")
}

fn sign_message(
    rng: &mut ChaCha20Rng,
    keypair: &KeyPair<'static, FieldElement>,
    msg: &[u8],
) -> DilithiumSignature<'static, FieldElement> {
    let _ = rng; // signatures are deterministic; RNG only drives keygen seeds
    keypair
        .sign(msg)
        .expect("signature generation should succeed")
}

fn verify_signature(
    pk: &PublicKey<'static, FieldElement>,
    msg: &[u8],
    sig: &DilithiumSignature<'static, FieldElement>,
) -> bool {
    pk.verify(msg, sig)
}

// ──────────────────────────────────────────────────────────────
//  EXTRA HELPERS FOR EQUATION (1)
// ──────────────────────────────────────────────────────────────
//
// Equation (1) in Section 4.2: :contentReference[oaicite:2]{index=2}
//
//   w - c s2 = A y - c s2
//             = A (z - c s1) - c s2
//             = A z - c t
//
// We can reconstruct all the pieces for the *final* accepted
// signature (z, h, c) using only the secret key and the signature.
//
// Needed building blocks (adapt these to match your actual API):
//
//   • Extract (ρ, s1, s2, t) from SecretKey.
//   • Extract (z, c) from Signature.
//   • A := matrix_from_rho(ρ)
//   • y := z - c * s1
//   • w := A * y
//   • lhs := w - c * s2
//   • rhs := A * z - c * t
//   • assert_eq!(lhs, rhs)

/// Matrix type for A, and vector/polynomial vector types for
/// s1, s2, t, y, z, w.
type MatrixA = Matrix<'static, FieldElement>;
type PolyVecL<'a> = [Polynomial<'a, FieldElement>; L];
type PolyVecK<'a> = [Polynomial<'a, FieldElement>; K];
type Challenge<'a> = Polynomial<'a, FieldElement>;

/// Rebuild A from the public seed ρ stored in the secret key.
/// In the paper: A ∼ R_q^{k×ℓ} ← Sam(ρ). :contentReference[oaicite:3]{index=3}
fn matrix_from_rho(rho: &[u8; 32]) -> MatrixA {
    expand_a_from_rho::<FieldElement>(*rho)
}

/// Multiply matrix A by a length-ℓ vector y to get w = A y in R_q^k.
fn mat_vec_mul(a: &MatrixA, y: &PolyVecL<'_>) -> PolyVecK<'static> {
    a.matrix_mul_output(y)
        .expect("matrix-vector multiplication should succeed")
}

/// Multiply a challenge `c` (the sparse polynomial) with a length-ℓ
/// vector `s1`, returning c * s1 (component-wise polynomial mult).
fn mul_challenge_vec_l<'a>(
    c: &Challenge<'a>,
    v: &PolyVecL<'a>,
) -> PolyVecL<'a> {
    std::array::from_fn(|idx| v[idx].clone() * c.clone())
}

/// Same but for length-k vectors (s2, t).
fn mul_challenge_vec_k<'a>(
    c: &Challenge<'a>,
    v: &PolyVecK<'a>,
) -> PolyVecK<'a> {
    std::array::from_fn(|idx| v[idx].clone() * c.clone())
}

/// Subtract two length-ℓ vectors: a - b.
fn sub_vec_l<'a>(a: &PolyVecL<'a>, b: &PolyVecL<'a>) -> PolyVecL<'a> {
    std::array::from_fn(|idx| a[idx].clone() - b[idx].clone())
}

/// Subtract two length-k vectors: a - b.
fn sub_vec_k<'a>(a: &PolyVecK<'a>, b: &PolyVecK<'a>) -> PolyVecK<'a> {
    std::array::from_fn(|idx| a[idx].clone() - b[idx].clone())
}

// You may prefer to use your own `eq` / `is_zero` methods instead.
fn vec_k_eq<'a>(a: &PolyVecK<'a>, b: &PolyVecK<'a>) -> bool {
    a == b
}

// ──────────────────────────────────────────────────────────────
//  TEST 1: Sign-then-Verify Correctness (Section 4.2 summary)
// ──────────────────────────────────────────────────────────────

#[test]
fn dilithium_sign_verify_roundtrip_is_correct() {
    // Deterministic RNG so the test is reproducible.
    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

    const NUM_KEYS: usize = 8;
    const MSGS_PER_KEY: usize = 8;

    for key_index in 0..NUM_KEYS {
        let keypair = generate_keypair(&mut rng);

        for msg_index in 0..MSGS_PER_KEY {
            let mut msg = [0u8; 64];
            rng.fill_bytes(&mut msg);

            let sig = sign_message(&mut rng, &keypair, &msg);
            let ok = verify_signature(&keypair.public, &msg, &sig);

            assert!(
                ok,
                "Correctness failed for key #{key_index}, message #{msg_index}: \
                 Verify(pk, m, Sign(sk, m)) != 1"
            );
        }
    }
}

// ──────────────────────────────────────────────────────────────
//  TEST 2: Equation (1) from Section 4.2
//          w - c s2 = A z - c t
// ──────────────────────────────────────────────────────────────
//
// For each random key/message/signature, we reconstruct:
//   • z, c from the Signature,
//   • ρ, s1, s2, t from the SecretKey,
//   • A from ρ,
// then derive:
//   y  = z - c*s1
//   w  = A*y
//   LHS = w - c*s2
//   RHS = A*z - c*t
// and assert LHS == RHS in R_q^k.
//
// This directly checks the algebraic identity proved as Eq. (1)
// in Section 4.2. :contentReference[oaicite:4]{index=4}

#[test]
fn equation_1_holds_for_random_signatures() {
    let mut rng = ChaCha20Rng::from_seed([1u8; 32]);

    const NUM_TRIALS: usize = 32;

    for trial in 0..NUM_TRIALS {
        let keypair = generate_keypair(&mut rng);

        let mut msg = [0u8; 64];
        rng.fill_bytes(&mut msg);

        let sig = sign_message(&mut rng, &keypair, &msg);
        let sk = keypair.secret();

        // ── 1. Extract secret-key components: ρ, s1, s2, t ───────
        //
        // In the paper: sk = (ρ, s1, s2, t). :contentReference[oaicite:5]{index=5}
        let rho = &keypair.public.rho;
        let s1 = sk.s1;
        let s2 = sk.s2;
        let t = &keypair.public.t;

        // ── 2. Extract signature components: z, c ───────────────
        //
        // In the paper: σ = (z, h, c). We only need z, c here.
        let z = &sig.z;
        let c = &sig.c;

        // ── 3. Rebuild A from ρ ─────────────────────────────────
        let a = matrix_from_rho(rho);

        // ── 4. Recover y from z and c,s1: y = z - c*s1 ─────────
        let cs1 = mul_challenge_vec_l(c, s1);
        let y = sub_vec_l(z, &cs1);

        // ── 5. Compute w = A*y ─────────────────────────────────
        let w = mat_vec_mul(&a, &y);

        // ── 6. Compute LHS = w - c*s2 ──────────────────────────
        let cs2 = mul_challenge_vec_k(c, s2);
        let lhs = sub_vec_k(&w, &cs2);

        // ── 7. Compute RHS = A*z - c*t ─────────────────────────
        let az = mat_vec_mul(&a, z);
        let ct = mul_challenge_vec_k(c, t);
        let rhs = sub_vec_k(&az, &ct);

        assert!(
            vec_k_eq(&lhs, &rhs),
            "Equation (1) failed in trial {trial}: w - c*s2 != A*z - c*t"
        );

        // Optional extra sanity: the signature should still verify.
        assert!(
            verify_signature(&keypair.public, &msg, &sig),
            "Signature that satisfied Eq. (1) did not verify (this should not happen)"
        );
    }
}
