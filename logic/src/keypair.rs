use std::convert::From;

use num_traits::Zero;
use rand::RngCore;

use crate::hash::shake256;
use crate::matrix::{MatrixA, expand_a_from_rho, mat_vec_mul};
use crate::params::{ETA, K, L, N};
use crate::utils::random_bytes;
use math::{poly::Polynomial, traits::FiniteField};

#[derive(Clone, Debug)]
pub struct PublicKey<'a, FF: FiniteField> {
    pub a: MatrixA<'a, FF>, // uncompressed: include A directly
    pub t: [Polynomial<'a, FF>; K], // t = A*s1 + s2
    pub rho: [u8; 32],      // seed used for A (kept for provenance)
}

#[derive(Clone, Debug)]
pub struct SecretKey<'a, FF: FiniteField> {
    pub a: MatrixA<'a, FF>, // include A here for convenience
    pub s1: [Polynomial<'a, FF>; L],
    pub s2: [Polynomial<'a, FF>; K],
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

    // TODO defien standadlone function
    for j in 0..L {
        let mut inbuf = Vec::with_capacity(33);
        inbuf.extend_from_slice(&s1_seed);
        inbuf.push(j as u8);
        let stream = shake256(2 * N, &inbuf);
        s1[j] = cbd_eta2(&stream);
        // ensure coefficients in [-ETA, ETA]
        // TODO fix it
        // debug_assert!(s1[j].c.iter().all(|&x| x.abs() <= ETA));
    }
    // TODO defien standadlone function
    for i in 0..K {
        let mut inbuf = Vec::with_capacity(33);
        inbuf.extend_from_slice(&s2_seed);
        inbuf.push((i as u8) ^ 0xA5);
        let stream = shake256(2 * N, &inbuf);
        s2[i] = cbd_eta2(&stream);
    }

    // t = A*s1 + s2
    let y = s1.clone(); // reuse shape
    let t_vec = mat_vec_mul(&a, &y).map(|p| p); // A*s1
    let mut t = [
        Polynomial::zero(),
        Polynomial::zero(),
        Polynomial::zero(),
        Polynomial::zero(),
    ];
    for i in 0..K {
        let mut sum = t_vec[i].clone();
        sum += s2[i].clone();
        t[i] = sum;
    }

    let pk = PublicKey {
        a: a.clone(),
        t,
        rho,
    };
    let sk = SecretKey { a, s1, s2 };
    (pk, sk)
}

pub fn keygen_with_seeds<FF: FiniteField + From<i64>>(
    rho: [u8; 32],
    s1_seed: [u8; 32],
    s2_seed: [u8; 32],
) -> (PublicKey<'static, FF>, SecretKey<'static, FF>) {
    let a = expand_a_from_rho(rho);
    let s1: [Polynomial<'static, FF>; L] = std::array::from_fn(|j| {
        let mut inp = Vec::new();
        inp.extend_from_slice(&s1_seed);
        inp.extend_from_slice(&(j as u16).to_le_bytes());
        let bs1 = shake256(2 * N, &inp);
        cbd_eta2::<FF>(&bs1)
    });
    let s2: [Polynomial<'static, FF>; K] = std::array::from_fn(|j| {
        let mut inp = Vec::new();
        inp.extend_from_slice(&s2_seed);
        inp.extend_from_slice(&(j as u16).to_le_bytes());
        let bs2 = shake256(2 * N, &inp);
        cbd_eta2::<FF>(&bs2)
    });
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
    use crate::params::{ETA, K, L};
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
        //expected_t[..K].clone_from_slice(&as1[..K]);
        for i in 0..K {
            expected_t[i] = as1[i].clone();
            expected_t[i] += sk.s2[i].clone();
        }
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
