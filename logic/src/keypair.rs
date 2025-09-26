use std::convert::From;

use num_traits::Zero;
use rand::RngCore;

use crate::hash::shake256;
use crate::matrix::{MatrixA, expand_a_from_rho, mat_vec_mul};
use crate::params::{K, L, N};
use crate::utils::random_bytes;
use math::{poly::Polynomial, traits::FiniteField};

#[derive(Clone, Debug)]
pub struct PublicKey<'a, FF: FiniteField> {
    pub a: MatrixA<'a, FF>, // uncompressed: include A directly
    pub t: [Polynomial<'a, FF>; K], // t = A*s1 + s2
    pub rho: [u8; 32],      // seed used for A (kept for provenance)
}

impl<'a, FF: FiniteField> PublicKey<'a, FF> {
    /// Construct a public key from its parts.
    /// Note: this does not verify that `t == A*s1 + s2`—it just packages fields.
    pub fn new(
        a: MatrixA<'a, FF>,
        t: [Polynomial<'a, FF>; K],
        rho: [u8; 32],
    ) -> Self {
        Self { a, t, rho }
    }
}

#[derive(Clone, Debug)]
pub struct SecretKey<'a, FF: FiniteField> {
    pub a: MatrixA<'a, FF>, // include A here for convenience
    pub s1: [Polynomial<'a, FF>; L],
    pub s2: [Polynomial<'a, FF>; K],
}

impl<'a, FF: FiniteField> SecretKey<'a, FF> {
    /// Construct a secret key from its parts.
    pub fn new(
        a: MatrixA<'a, FF>,
        s1: [Polynomial<'a, FF>; L],
        s2: [Polynomial<'a, FF>; K],
    ) -> Self {
        Self { a, s1, s2 }
    }
}

/// CBD for η=2 from an XOF stream
fn cbd_eta2<FF: FiniteField + From<i64>>(
    stream: &[u8],
) -> Polynomial<'static, FF> {
    // Each coefficient uses 4 bits: (b0 + b1) - (b2 + b3)
    let mut out = [0i64; N];
    let mut bitpos = 0usize;
    for out_i in out.iter_mut().take(N) {
        let mut bits = 0u32;
        // pull 2 bytes (16 bits) for simplicity
        let byte_idx = bitpos / 8;
        let two =
            u16::from_le_bytes([stream[byte_idx], stream[byte_idx + 1]]) as u32;
        bits = two;
        let a0 = (bits & 1) + ((bits >> 1) & 1);
        let a1 = ((bits >> 2) & 1) + ((bits >> 3) & 1);
        *out_i = (a0 as i64) - (a1 as i64); // in [-2,2]
        bitpos += 4;
    }
    out.into()
}

// Put this near the top of the file (or in a utils module) so `keygen` can call it.

/// Fill `out` with polynomials sampled as:
///   poly_i = cbd_eta2( SHAKE256(2*N, seed || tag(i)) )
///
/// - `out`: destination slice of polynomials (e.g., s1 or s2)
/// - `seed`: 32‑byte seed (or any length; we only append one tag byte)
/// - `make_tag`: computes the last byte appended to the seed for index `i`
fn fill_secret_polys<F, FF: FiniteField + From<i64>>(
    out: &mut [Polynomial<'static, FF>],
    seed: &[u8],
    make_tag: F,
) where
    F: Fn(usize) -> u8,
{
    // Build "seed || tag" once; mutate the last byte per iteration.
    let mut inbuf = Vec::with_capacity(seed.len() + 1);
    inbuf.extend_from_slice(seed);
    inbuf.push(0u8); // placeholder tag; overwritten below

    for (i, dst) in out.iter_mut().enumerate() {
        // Set the per-index tag byte.
        *inbuf.last_mut().expect("tag byte exists") = make_tag(i);

        // Expand and sample.
        let stream = shake256(2 * N, &inbuf);
        *dst = cbd_eta2(&stream);
    }
}

pub fn keygen<FF: FiniteField + From<i64>>()
-> (PublicKey<'static, FF>, SecretKey<'static, FF>) {
    // Generate rho (seed for A), and seeds for s1, s2
    let mut rng = rand::thread_rng();
    let mut rho = [0u8; 32];
    rng.fill_bytes(&mut rho);
    let a = expand_a_from_rho(rho);

    // Expand s1, s2 with SHAKE256 streams (deterministic from seeds)
    let s1_seed = random_bytes();
    let s2_seed = random_bytes();

    let mut s1 = [
        Polynomial::zero(),
        Polynomial::zero(),
        Polynomial::zero(),
        Polynomial::zero(),
    ];
    let mut s2 = [
        Polynomial::zero(),
        Polynomial::zero(),
        Polynomial::zero(),
        Polynomial::zero(),
    ];

    // s1: tag = j as u8
    fill_secret_polys(&mut s1, &s1_seed, |j| j as u8);

    // s2: tag = (i as u8) ^ 0xA5
    fill_secret_polys(&mut s2, &s2_seed, |i| (i as u8) ^ 0xA5);

    // t = A*s1 + s2
    let y = s1.clone(); // reuse shape
    let t_vec = mat_vec_mul(&a, &y).map(|p| p); // A*s1
    let mut t = [
        Polynomial::zero(),
        Polynomial::zero(),
        Polynomial::zero(),
        Polynomial::zero(),
    ];

    for ((out, a), b) in t.iter_mut().zip(&t_vec).zip(&s2) {
        *out = std::iter::once(b.clone()).fold(a.clone(), |mut acc, x| {
            acc += x;
            acc
        })
    }

    (PublicKey::new(a.clone(), t, rho), SecretKey::new(a, s1, s2))
}



#[inline]
fn expand_secret_array<const LEN: usize, FF>(
    seed: impl AsRef<[u8]>,
) -> [math::poly::Polynomial<'static, FF>; LEN]
where
    FF: math::traits::FiniteField + From<i64>,
{
    use std::array::from_fn;
    let seed = seed.as_ref();

    from_fn(|j| {
        let mut inp = Vec::with_capacity(seed.len() + 2);
        inp.extend_from_slice(seed);
        // match keygen_with_seeds: index encoded as little-endian u16
        inp.extend_from_slice(&(j as u16).to_le_bytes());
        let bytes = crate::hash::shake256(2 * crate::params::N, &inp);
        cbd_eta2::<FF>(&bytes)
    })
}

pub fn keygen_with_seeds<FF: FiniteField + From<i64>>(
    rho: [u8; 32],
    s1_seed: [u8; 32],
    s2_seed: [u8; 32],
) -> (PublicKey<'static, FF>, SecretKey<'static, FF>) {
    let a = expand_a_from_rho(rho);

    let s1: [Polynomial<'static, FF>; L] =
        expand_secret_array::<L, FF>(&s1_seed);

    let s2: [Polynomial<'static, FF>; K] =
        expand_secret_array::<K, FF>(&s2_seed);

    let t_vec = mat_vec_mul(&a, &s1);
    let t: [Polynomial<'static, FF>; K] = std::array::from_fn(|i| {
        let mut sum = t_vec[i].clone();
        sum += s2[i].clone();
        sum
    });

    (
        PublicKey {
            a: a.clone(),
            t,
            rho,
        },
        SecretKey { a, s1, s2 },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::expand_a_from_rho;
    use crate::matrix::mat_vec_mul;
    use crate::params::{K, L};
    use math::field_element::FieldElement;
    use math::poly::Polynomial;

    #[test]
    fn keygen_shapes_and_secret_bounds() {
        let (pk, sk) = keygen::<FieldElement>();

        // Shapes
        assert_eq!(pk.a.a.len(), K);
        for i in 0..K {
            assert_eq!(pk.a.a[i].len(), L);
        }
        assert_eq!(pk.t.len(), K);
        assert_eq!(sk.s1.len(), L);
        assert_eq!(sk.s2.len(), K);
    }

    #[test]
    fn public_matrix_matches_rho_and_secret_matrix() {
        let (pk, sk) = keygen::<FieldElement>();

        // pk.a should equal expand_a_from_rho(pk.rho)
        let a_from_rho = expand_a_from_rho(pk.rho);
        for i in 0..K {
            for j in 0..L {
                assert_eq!(
                    pk.a.a[i][j], a_from_rho.a[i][j],
                    "A mismatch at {},{}",
                    i, j
                );
            }
        }

        // pk.a and sk.a should match
        for i in 0..K {
            for j in 0..L {
                assert_eq!(
                    pk.a.a[i][j], sk.a.a[i][j],
                    "A(pk) != A(sk) at {},{}",
                    i, j
                );
            }
        }
    }

    #[test]
    fn t_equals_a_times_s1_plus_s2() {
        let (pk, sk) = keygen::<FieldElement>();

        let as1 = mat_vec_mul(&sk.a, &sk.s1);
        let mut expected_t = [
            Polynomial::zero(),
            Polynomial::zero(),
            Polynomial::zero(),
            Polynomial::zero(),
        ];

        for i in 0..K {
            expected_t[i] = as1[i].clone();
            expected_t[i] += sk.s2[i].clone();
        }
        // expected_t[..K].clone_from_slice(&as1[..K]);
        // expected_t[..K].clone_from_slice(&sk.s2[..K]);
        for i in 0..K {
            assert_eq!(pk.t[i], expected_t[i], "t mismatch at row {}", i);
        }
    }

    #[test]
    fn keygen_is_deterministic_given_seeds() {
        let rho = [3u8; 32];
        let s1 = [5u8; 32];
        let s2 = [7u8; 32];
        let (pk1, sk1) = keygen_with_seeds::<FieldElement>(rho, s1, s2);
        let (pk2, sk2) = keygen_with_seeds::<FieldElement>(rho, s1, s2);

        // Matrices and secrets identical
        assert_eq!(pk1.a.a, pk2.a.a);
        for i in 0..K {
            assert_eq!(sk1.s2[i], sk2.s2[i]);
        }
        for j in 0..L {
            assert_eq!(sk1.s1[j], sk2.s1[j]);
        }
        for i in 0..K {
            assert_eq!(pk1.t[i], pk2.t[i]);
        }
    }
}
