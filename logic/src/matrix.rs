use crate::hash::shake128;
use crate::params::{K, L, N};
use math::{poly::Polynomial, traits::FiniteField};
use num_traits::Zero;

#[derive(Clone, Debug)]
pub struct MatrixA<'a, FF: FiniteField> {
    pub a: Vec<Vec<Polynomial<'a, FF>>>, // K x L
}

/// Expand A from rho using SHAKE128 as XOF (educational: uses modulo reduction).
pub fn expand_a_from_rho<FF: FiniteField + std::convert::From<i64>>(
    rho: [u8; 32],
) -> MatrixA<'static, FF> {
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
            let mut coeffs = [0u32; N];
            for t in 0..N {
                let b = &stream[4 * t..4 * t + 4];
                let v = u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
                coeffs[t] = v;
            }
            row.push(Polynomial::from(coeffs));
        }
        mat.push(row);
    }
    MatrixA { a: mat }
}

pub fn mat_vec_mul<FF: FiniteField>(
    a: &MatrixA<FF>,
    y: &[Polynomial<FF>; L],
) -> [Polynomial<'static, FF>; K] {
    std::array::from_fn(|i| {
        a.a[i].iter().zip(y.iter()).fold(
            Polynomial::zero(),
            |mut acc, (aij, yj)| {
                acc += aij.clone() * yj.clone();
                acc
            },
        )
    })
}

// TODO add more tests
#[cfg(test)]
mod tests {
    use crate::matrix::{expand_a_from_rho, mat_vec_mul};
    use crate::params::{K, L, N, Q};
    //use crate::poly::{Poly, mod_q};
    use math::poly::Polynomial;
    use num_traits::Zero;

    use math::field_element::FieldElement;
    use math::traits::FiniteField;

    // TODO fix it
    // fn unit_poly_at<FF: FiniteField>(idx: usize, val: i64) -> Polynomial<'static, FF> {
    //     let mut p = Polynomial::zero();
    //     p.c[idx] = mod_q(val);
    //     p
    // }

    fn zero_y<FF: FiniteField>() -> [Polynomial<'static, FF>; L] {
        [
            Polynomial::zero(),
            Polynomial::zero(),
            Polynomial::zero(),
            Polynomial::zero(),
        ]
    }

    #[test]
    fn expand_a_is_deterministic_and_modq() {
        let rho = [7u8; 32];
        let a1 = expand_a_from_rho::<FieldElement>(rho);
        let a2 = expand_a_from_rho::<FieldElement>(rho);
        assert_eq!(a1.a.len(), K);
        assert_eq!(a2.a.len(), K);
        for i in 0..K {
            assert_eq!(a1.a[i].len(), L);
            assert_eq!(a2.a[i].len(), L);
            for j in 0..L {
                // determinism
                assert_eq!(a1.a[i][j], a2.a[i][j]);
                // coeffs in [0, Q)
                for &c in a1.a[i][j].coefficients() {
                    let num: i64 = c.into();
                    assert!((0..Q).contains(&num), "coefficient out of range");
                }
            }
        }
    }

    #[test]
    fn expand_a_changes_with_rho() {
        let a1 = expand_a_from_rho::<FieldElement>([1u8; 32]);
        let a2 = expand_a_from_rho::<FieldElement>([2u8; 32]);
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
        let a = expand_a_from_rho::<FieldElement>([42u8; 32]);
        let y = zero_y();
        let w = mat_vec_mul(&a, &y);
        for i in 0..K {
            assert_eq!(w[i], Polynomial::zero());
        }
    }
}
