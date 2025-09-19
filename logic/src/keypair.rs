use crate::hash::shake256;
use crate::matrix::{MatrixA, expand_a_from_rho, mat_vec_mul};
use crate::params::{ETA, K, L, N};
use crate::poly::Poly;
use rand::RngCore;

#[derive(Clone, Debug)]
pub struct PublicKey {
    pub a: MatrixA,    // uncompressed: include A directly
    pub t: [Poly; K],  // t = A*s1 + s2
    pub rho: [u8; 32], // seed used for A (kept for provenance)
}

#[derive(Clone, Debug)]
pub struct SecretKey {
    pub a: MatrixA, // include A here for convenience
    pub s1: [Poly; L],
    pub s2: [Poly; K],
}

/// CBD for η=2 from an XOF stream
fn cbd_eta2(stream: &[u8]) -> Poly {
    // Each coefficient uses 4 bits: (b0 + b1) - (b2 + b3)
    let mut out = [0i64; N];
    let mut bitpos = 0usize;
    for i in 0..N {
        let mut bits = 0u32;
        // pull 2 bytes (16 bits) for simplicity
        let byte_idx = bitpos / 8;
        let two =
            u16::from_le_bytes([stream[byte_idx], stream[byte_idx + 1]]) as u32;
        bits = two;
        let a0 = ((bits >> 0) & 1) + ((bits >> 1) & 1);
        let a1 = ((bits >> 2) & 1) + ((bits >> 3) & 1);
        out[i] = (a0 as i64) - (a1 as i64); // in [-2,2]
        bitpos += 4;
    }
    Poly { c: out }
}

pub fn keygen() -> (PublicKey, SecretKey) {
    // Generate rho (seed for A), and seeds for s1, s2
    let mut rng = rand::thread_rng();
    let mut rho = [0u8; 32];
    rng.fill_bytes(&mut rho);
    let a = expand_a_from_rho(rho);

    // Expand s1, s2 with SHAKE256 streams (deterministic from seeds)
    let s1_seed = {
        let mut tmp = [0u8; 32];
        rng.fill_bytes(&mut tmp);
        tmp
    };
    let s2_seed = {
        let mut tmp = [0u8; 32];
        rng.fill_bytes(&mut tmp);
        tmp
    };

    let mut s1 = [Poly::zero(), Poly::zero(), Poly::zero(), Poly::zero()];
    let mut s2 = [Poly::zero(), Poly::zero(), Poly::zero(), Poly::zero()];

    for j in 0..L {
        let mut inbuf = Vec::with_capacity(33);
        inbuf.extend_from_slice(&s1_seed);
        inbuf.push(j as u8);
        let stream = shake256(2 * N, &inbuf);
        s1[j] = cbd_eta2(&stream);
        // ensure coefficients in [-ETA, ETA]
        debug_assert!(s1[j].c.iter().all(|&x| x.abs() <= ETA));
    }
    for i in 0..K {
        let mut inbuf = Vec::with_capacity(33);
        inbuf.extend_from_slice(&s2_seed);
        inbuf.push((i as u8) ^ 0xA5);
        let stream = shake256(2 * N, &inbuf);
        s2[i] = cbd_eta2(&stream);
        debug_assert!(s2[i].c.iter().all(|&x| x.abs() <= ETA));
    }

    // t = A*s1 + s2
    let y = s1.clone(); // reuse shape
    let t_vec = mat_vec_mul(&a, &y).map(|p| p); // A*s1
    let mut t = [Poly::zero(), Poly::zero(), Poly::zero(), Poly::zero()];
    for i in 0..K {
        let mut sum = t_vec[i].clone();
        sum.add_assign(&s2[i]);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::expand_a_from_rho;
    use crate::matrix::mat_vec_mul;
    use crate::params::{ETA, K, L};
    use crate::poly::Poly;

    #[test]
    fn keygen_shapes_and_secret_bounds() {
        let (pk, sk) = keygen();

        // Shapes
        assert_eq!(pk.a.a.len(), K);
        for i in 0..K {
            assert_eq!(pk.a.a[i].len(), L);
        }
        assert_eq!(pk.t.len(), K);
        assert_eq!(sk.s1.len(), L);
        assert_eq!(sk.s2.len(), K);

        // Secret coeff bounds
        for j in 0..L {
            assert!(sk.s1[j].c.iter().all(|&x| x.abs() <= ETA));
        }
        for i in 0..K {
            assert!(sk.s2[i].c.iter().all(|&x| x.abs() <= ETA));
        }
    }

    #[test]
    fn public_matrix_matches_rho_and_secret_matrix() {
        let (pk, sk) = keygen();

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
        let (pk, sk) = keygen();

        let as1 = mat_vec_mul(&sk.a, &sk.s1);
        let mut expected_t =
            [Poly::zero(), Poly::zero(), Poly::zero(), Poly::zero()];
        for i in 0..K {
            expected_t[i] = as1[i].clone();
            expected_t[i].add_assign(&sk.s2[i]);
        }
        for i in 0..K {
            assert_eq!(pk.t[i], expected_t[i], "t mismatch at row {}", i);
        }
    }
}
