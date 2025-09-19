use crate::params::{N, Q};
use core::ops::{Add, Sub};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly {
    pub c: [i64; N],
}

impl Poly {
    pub fn zero() -> Self {
        Self { c: [0; N] }
    }
    pub fn from_slice(sl: &[i64]) -> Self {
        assert_eq!(sl.len(), N);
        let mut p = Self::zero();
        p.c.copy_from_slice(sl);
        p
    }
    pub fn reduce_mod_q(&mut self) {
        for x in &mut self.c {
            *x = mod_q(*x);
        }
    }
    pub fn add_assign(&mut self, other: &Poly) {
        for i in 0..N {
            self.c[i] = mod_q(self.c[i] + other.c[i]);
        }
    }
    pub fn sub_assign(&mut self, other: &Poly) {
        for i in 0..N {
            self.c[i] = mod_q(self.c[i] - other.c[i]);
        }
    }
    pub fn mul(&self, rhs: &Poly) -> Poly {
        // Schoolbook negacyclic convolution in Z_q[X]/(X^N + 1)
        let mut acc = [0i128; N];
        for i in 0..N {
            let ai = self.c[i] as i128;
            for j in 0..N {
                let idx = i + j;
                let prod = ai * (rhs.c[j] as i128);
                if idx < N {
                    acc[idx] += prod;
                } else {
                    acc[idx - N] -= prod; // negacyclic: X^N == -1
                }
            }
        }
        let mut out = Poly::zero();
        for i in 0..N {
            // reduce to [0, Q)
            let mut v = acc[i] % (Q as i128);
            if v < 0 {
                v += Q as i128;
            }
            out.c[i] = v as i64;
        }
        out
    }
    pub fn mul_assign(&mut self, rhs: &Poly) {
        let tmp = self.mul(rhs);
        *self = tmp;
    }
    pub fn scale_add_assign(&mut self, rhs: &Poly, k: i64) {
        for i in 0..N {
            self.c[i] = mod_q(self.c[i] + k * rhs.c[i]);
        }
    }
    pub fn norm_inf(&self) -> i64 {
        self.c.iter().map(|&x| centered(x).abs()).max().unwrap_or(0)
    }
}

/// Modular reduction to [0, Q)
pub fn mod_q(x: i64) -> i64 {
    let mut r = x % Q;
    if r < 0 {
        r += Q;
    }
    r
}

/// Center to (-Q/2, Q/2]
pub fn centered(x: i64) -> i64 {
    let r = mod_q(x);
    if r > Q / 2 { r - Q } else { r }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{N, Q};

    // Helper: build a polynomial with given (index, value) pairs.
    // Values are reduced mod q so they’re safe for mul().
    fn p(vals: &[(usize, i64)]) -> Poly {
        let mut c = [0i64; N];
        for &(i, v) in vals {
            c[i] = mod_q(v);
        }
        Poly { c }
    }

    #[test]
    fn mod_q_basic() {
        assert_eq!(mod_q(0), 0);
        assert_eq!(mod_q(Q), 0);
        assert_eq!(mod_q(1), 1);
        assert_eq!(mod_q(Q - 1), Q - 1);
        assert_eq!(mod_q(-1), Q - 1);
        assert_eq!(mod_q(Q + 5), 5);
        assert_eq!(mod_q(-Q - 5), Q - 5);
        assert_eq!(mod_q(2 * Q + 7), 7);
        assert_eq!(mod_q(-2 * Q - 7), Q - 7);
    }

    #[test]
    fn centered_basic() {
        let q2 = Q / 2;
        assert_eq!(centered(0), 0);
        assert_eq!(centered(Q), 0); // wraps to 0 then centered
        assert_eq!(centered(1), 1);
        assert_eq!(centered(Q - 1), -1);
        assert_eq!(centered(q2), q2); // tie goes to positive side
        assert_eq!(centered(q2 + 1), -q2); // just over half rounds negative
    }

    #[test]
    fn reduce_mod_q_normalizes() {
        let mut a = Poly::zero();
        a.c[0] = -1;
        a.c[1] = Q;
        a.c[2] = Q + 2;
        a.c[3] = -Q - 3;
        a.c[4] = 0;
        a.c[5] = Q - 1;

        a.reduce_mod_q();

        assert_eq!(a.c[0], Q - 1);
        assert_eq!(a.c[1], 0);
        assert_eq!(a.c[2], 2);
        assert_eq!(a.c[3], Q - 3);
        assert_eq!(a.c[4], 0);
        assert_eq!(a.c[5], Q - 1);
    }

    #[test]
    fn add_sub_wrap() {
        let a = p(&[(0, Q - 1), (1, 2)]);
        let b = p(&[(0, 2), (1, Q - 1)]);

        let mut sum = a.clone();
        sum.add_assign(&b);
        assert_eq!(sum.c[0], 1); // (Q-1 + 2) mod Q = 1
        assert_eq!(sum.c[1], 1); // (2 + Q-1) mod Q = 1

        let mut diff = a.clone();
        diff.sub_assign(&b);
        assert_eq!(diff.c[0], Q - 3); // (Q-1 - 2) mod Q = Q-3
        assert_eq!(diff.c[1], 3); // (2 - (Q-1)) mod Q = 3
    }

    #[test]
    fn mul_zero_one_constant() {
        // a(x) = 1 + 2x + 3x^2 + 4x^{N-1}
        let a = p(&[(0, 1), (1, 2), (2, 3), (N - 1, 4)]);

        // Zero
        let z = Poly::zero();
        let prod0 = a.mul(&z);
        assert_eq!(prod0, Poly::zero());

        // One
        let one = p(&[(0, 1)]);
        let prod1 = a.mul(&one);
        assert_eq!(prod1, a);

        // Constant 5
        let five = p(&[(0, 5)]);
        let prod5 = a.mul(&five);
        let mut exp = [0i64; N];
        exp[0] = mod_q(5 * 1);
        exp[1] = mod_q(5 * 2);
        exp[2] = mod_q(5 * 3);
        exp[N - 1] = mod_q(5 * 4);
        assert_eq!(prod5.c, exp);
    }

    #[test]
    fn mul_assign_matches_mul() {
        let a = p(&[(1, 7), (3, 11), (N - 1, 5)]);
        let b = p(&[(0, 9), (1, 2)]);

        let x = a.mul(&b);
        let mut y = a.clone();
        y.mul_assign(&b);
        assert_eq!(x, y);
    }

    #[test]
    fn mul_by_x_is_negacyclic_rotation() {
        // a(x) = 1 + 2x + 3x^2 + 4x^{N-1}
        let a = p(&[(0, 1), (1, 2), (2, 3), (N - 1, 4)]);
        let x = p(&[(1, 1)]); // X

        let ax = a.mul(&x);

        // Expected: multiply by X rotates left by 1 with negacyclic wrap:
        // coeff[0] = -a[N-1]  (mod q) = Q-4
        // coeff[i] = a[i-1] for i >= 1
        let mut exp = [0i64; N];
        exp[0] = mod_q(-(4));
        exp[1] = mod_q(1);
        exp[2] = mod_q(2);
        exp[3] = mod_q(3);
        assert_eq!(ax.c[0], exp[0]);
        assert_eq!(ax.c[1], exp[1]);
        assert_eq!(ax.c[2], exp[2]);
        assert_eq!(ax.c[3], exp[3]);

        // All other coefficients should be 0.
        for i in 4..(N - 1) {
            assert_eq!(ax.c[i], 0, "unexpected nonzero at {}", i);
        }
        assert_eq!(ax.c[N - 1], 0);
    }

    #[test]
    fn x_times_xn_minus_1_is_minus_one() {
        let x = p(&[(1, 1)]); // X
        let xn1 = p(&[(N - 1, 1)]); // X^{N-1}
        let prod = x.mul(&xn1); // X^N ≡ -1 (mod X^N + 1)

        let mut exp = [0i64; N];
        exp[0] = mod_q(-1);
        let expected = Poly { c: exp };
        assert_eq!(prod, expected);
    }

    #[test]
    fn distributivity_holds() {
        // a(x), b(x), c(x) with a few nonzero terms
        let a = p(&[(0, 1), (1, 2), (5, 7), (N - 1, 9)]);
        let b = p(&[(0, 3), (2, 4), (10, 5)]);
        let c = p(&[(0, 6), (1, 1), (3, 2), (N - 2, 7)]);

        let mut a_plus_b = a.clone();
        a_plus_b.add_assign(&b);
        let left = a_plus_b.mul(&c);

        let mut right = a.mul(&c);
        let bc = b.mul(&c);
        right.add_assign(&bc);

        assert_eq!(left, right);
    }

    #[test]
    fn scale_add_assign_works() {
        let a = p(&[(0, 10), (3, 20), (7, Q - 5)]);
        let b = p(&[(0, 1), (3, 2), (7, 3)]);

        let mut s = a.clone();
        s.scale_add_assign(&b, 2); // a + 2*b

        assert_eq!(s.c[0], mod_q(10 + 2 * 1));
        assert_eq!(s.c[3], mod_q(20 + 2 * 2));
        assert_eq!(s.c[7], mod_q((Q - 5) + 2 * 3));
    }

    #[test]
    fn norm_inf_is_from_centered_coeffs() {
        let mut c = [0i64; N];
        c[0] = 0; // centered 0
        c[1] = Q / 2; // centered +Q/2
        c[2] = Q / 2 + 1; // centered -Q/2
        c[3] = Q - 1; // centered -1
        c[4] = 5; // centered +5
        let p = Poly { c };
        assert_eq!(p.norm_inf(), Q / 2);
    }

    #[test]
    fn from_slice_copies_exactly() {
        let mut arr = [0i64; N];
        arr[0] = 123;
        arr[10] = Q - 2;
        arr[N - 1] = 42;
        let p = Poly::from_slice(&arr);
        assert_eq!(p.c, arr);
    }
}
