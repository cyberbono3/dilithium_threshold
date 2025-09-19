use crate::hash::shake128;
use crate::params::{K, L, N, Q};
use crate::poly::{Poly, mod_q};

#[derive(Clone, Debug)]
pub struct MatrixA {
    pub a: Vec<Vec<Poly>>, // K x L
}

/// Expand A from rho using SHAKE128 as XOF (educational: uses modulo reduction).
pub fn expand_a_from_rho(rho: [u8; 32]) -> MatrixA {
    let mut mat = Vec::with_capacity(K);
    for i in 0..K {
        let mut row = Vec::with_capacity(L);
        for j in 0..L {
            // Domain sep: rho || i || j
            let mut seed = Vec::with_capacity(32 + 4);
            seed.extend_from_slice(&rho);
            seed.extend_from_slice(&(i as u16).to_le_bytes());
            seed.extend_from_slice(&(j as u16).to_le_bytes());
            let stream = shake128(4 * N, &seed); // 4*N bytes -> N u32s
            let mut coeffs = [0i64; N];
            for t in 0..N {
                let b = &stream[4 * t..4 * t + 4];
                let v = u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as i64;
                coeffs[t] = mod_q(v % Q); // simple mod (note: slight bias; ok educationally)
            }
            row.push(Poly { c: coeffs });
        }
        mat.push(row);
    }
    MatrixA { a: mat }
}

pub fn mat_vec_mul(a: &MatrixA, y: &[Poly; L]) -> [Poly; K] {
    let mut out = [Poly::zero(), Poly::zero(), Poly::zero(), Poly::zero()];
    for i in 0..K {
        let mut acc = Poly::zero();
        for j in 0..L {
            let prod = a.a[i][j].mul(&y[j]);
            acc.add_assign(&prod);
        }
        out[i] = acc;
    }
    out
}

#[cfg(test)]
mod tests {
    use crate::matrix::{expand_a_from_rho, mat_vec_mul};
    use crate::params::{K, L, N, Q};
    use crate::poly::{Poly, mod_q};

    fn unit_poly_at(idx: usize, val: i64) -> Poly {
        let mut p = Poly::zero();
        p.c[idx] = mod_q(val);
        p
    }

    fn zero_y() -> [Poly; L] {
        [Poly::zero(), Poly::zero(), Poly::zero(), Poly::zero()]
    }

    #[test]
    fn expand_a_is_deterministic_and_modq() {
        let rho = [7u8; 32];
        let a1 = expand_a_from_rho(rho);
        let a2 = expand_a_from_rho(rho);
        assert_eq!(a1.a.len(), K);
        assert_eq!(a2.a.len(), K);
        for i in 0..K {
            assert_eq!(a1.a[i].len(), L);
            assert_eq!(a2.a[i].len(), L);
            for j in 0..L {
                // determinism
                assert_eq!(a1.a[i][j], a2.a[i][j]);
                // coeffs in [0, Q)
                for &c in &a1.a[i][j].c {
                    assert!(0 <= c && c < Q, "coefficient out of range");
                }
            }
        }
    }

    #[test]
    fn expand_a_changes_with_rho() {
        let a1 = expand_a_from_rho([1u8; 32]);
        let a2 = expand_a_from_rho([2u8; 32]);
        // It's overwhelmingly likely that at least one position differs.
        let mut any_diff = false;
        'outer: for i in 0..K {
            for j in 0..L {
                if a1.a[i][j] != a2.a[i][j] {
                    any_diff = true;
                    break 'outer;
                }
            }
        }
        assert!(
            any_diff,
            "Matrices for different rho should differ (with overwhelming probability)"
        );
    }

    #[test]
    fn mat_vec_mul_zero_is_zero() {
        let a = expand_a_from_rho([42u8; 32]);
        let y = zero_y();
        let w = mat_vec_mul(&a, &y);
        for i in 0..K {
            assert_eq!(w[i], Poly::zero());
        }
    }

    #[test]
    fn mat_vec_mul_basis_column_selects_column() {
        // y[j] = 1 (monomial), others zero => output column j of A
        let a = expand_a_from_rho([3u8; 32]);
        for j in 0..L {
            let mut y = zero_y();
            y[j] = unit_poly_at(0, 1); // 1 at X^0
            let w = mat_vec_mul(&a, &y);
            for i in 0..K {
                assert_eq!(
                    w[i], a.a[i][j],
                    "Mismatch at column {} row {}",
                    j, i
                );
            }
        }
    }

    #[test]
    fn mat_vec_mul_is_linear() {
        let a = expand_a_from_rho([99u8; 32]);

        // Build two deterministic y vectors with small support
        let mut y1 = zero_y();
        y1[0] = unit_poly_at(0, 1);
        y1[1] = unit_poly_at(5, 2);
        y1[2] = unit_poly_at(7, Q - 3); // -3 (mod q)
        y1[3] = unit_poly_at(1, 4);

        let mut y2 = zero_y();
        y2[0] = unit_poly_at(2, 9);
        y2[1] = unit_poly_at(0, 3);
        y2[2] = unit_poly_at(N - 1, 8);
        y2[3] = unit_poly_at(3, 1);

        let mut y_sum = y1.clone();
        for j in 0..L {
            y_sum[j].add_assign(&y2[j]);
        }

        let wy1 = mat_vec_mul(&a, &y1);
        let wy2 = mat_vec_mul(&a, &y2);
        let wysum = mat_vec_mul(&a, &y_sum);

        let mut wy1_plus_wy2 = wy1.clone();
        for i in 0..K {
            wy1_plus_wy2[i].add_assign(&wy2[i]);
        }

        assert_eq!(wysum, wy1_plus_wy2);
    }

    #[test]
    fn mat_vec_mul_outputs_are_mod_q() {
        let a = expand_a_from_rho([5u8; 32]);
        let mut y = zero_y();
        // Dense small values
        let mut p = Poly::zero();
        for i in 0..N {
            p.c[i] = mod_q((i as i64 % 11) - 5);
        }
        for j in 0..L {
            y[j] = p.clone();
        }

        let w = mat_vec_mul(&a, &y);
        for i in 0..K {
            for &c in &w[i].c {
                assert!(0 <= c && c < Q);
            }
        }
    }
}
