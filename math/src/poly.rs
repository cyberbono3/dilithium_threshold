use std::ops::{Add, AddAssign, Mul, Neg, Sub};

use rand::Rng;

use super::ntt::{intt, montgomery_reduce, ntt};

/// Constants for the Dilithium algorithm
pub const Q: i32 = 8380417; // Dilithium's prime modulus
pub const N: usize = 256; // Polynomial degree bound

/// Represents a polynomial in Rq = Zq[X]/(X^256 + 1).
///
/// Coefficients are stored as an array of integers modulo Q.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Polynomial {
    // TODO make field elements insterad of i32
    coeffs: [i32; N],
}

impl Polynomial {
    /// Initialize polynomial with given coefficients.
    /// Handles reduction modulo X^N + 1 for coefficients beyond degree N.
    pub fn new(coeffs: Vec<i32>) -> Self {
        let mut result = [0i32; N];

        if coeffs.len() <= N {
            // Simple case: just copy and pad with zeros
            result[..coeffs.len()].copy_from_slice(&coeffs);
        } else {
            // Reduce modulo X^N + 1
            for (i, &coeff) in coeffs.iter().enumerate() {
                let pos = i % N;
                let quotient = i / N;

                if quotient % 2 == 0 {
                    // Even powers of X^N contribute positively
                    result[pos] =
                        Self::mod_reduce(result[pos] as i64 + coeff as i64);
                } else {
                    // Odd powers of X^N contribute negatively (since X^N = -1)
                    result[pos] =
                        Self::mod_reduce(result[pos] as i64 - coeff as i64);
                }
            }
        }

        // Ensure all coefficients are in [0, Q)
        for coeff in &mut result {
            *coeff = Self::mod_reduce(*coeff as i64);
        }

        Self { coeffs: result }
    }

    pub fn coeffs(&self) -> [i32; N] {
        self.coeffs
    }

    /// Simple modular reduction ensuring result is in [0, Q)
    #[inline(always)]
    fn mod_reduce(a: i64) -> i32 {
        let mut t = (a % (Q as i64)) as i32;
        t += (t >> 31) & Q; // Add Q if t is negative
        t
    }

    /// Create zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: [0; N] }
    }

    /// Create polynomial representing 1.
    pub fn one() -> Self {
        let mut coeffs = [0; N];
        coeffs[0] = 1;
        Self { coeffs }
    }

    /// Generate random polynomial with coefficients in [0, bound).
    pub fn random(bound: i32) -> Self {
        let mut rng = rand::thread_rng();
        let mut coeffs = [0i32; N];
        for coeff in &mut coeffs {
            *coeff = rng.gen_range(0..bound);
        }
        Self::from(coeffs)
    }

    /// NTT-based multiplication (fast polynomial multiplication)
    pub fn ntt_multiply(&self, other: &Self) -> Self {
        // Copy coefficients to avoid modifying the original polynomials
        let mut a_ntt = self.coeffs;
        let mut b_ntt = other.coeffs;

        // Forward NTT transforms both polynomials to NTT domain
        ntt(&mut a_ntt);
        ntt(&mut b_ntt);

        // Pointwise multiplication in NTT domain
        for i in 0..N {
            a_ntt[i] = montgomery_reduce(a_ntt[i] as i64 * b_ntt[i] as i64);
        }

        // Inverse NTT to get back to coefficient domain
        intt(&mut a_ntt);

        Self { coeffs: a_ntt }
    }

    /// Compute infinity norm of polynomial.
    /// Returns the maximum absolute value of coefficients when centered around 0.
    /// TODO test it properly
    pub fn norm_infinity(&self) -> i32 {
        let q_half = Q / 2;
        self.coeffs
            .iter()
            .map(|&c| {
                let signed = if c > q_half { c - Q } else { c };
                signed.abs()
            })
            .max()
            .unwrap_or(0)
    }

    /// Compute L2 norm of polynomial.
    pub fn norm_l2(&self) -> f64 {
        let q_half = Q / 2;
        let sum: i64 = self
            .coeffs
            .iter()
            .map(|&c| {
                let signed = if c > q_half { c - Q } else { c };
                (signed as i64) * (signed as i64)
            })
            .sum();

        (sum as f64).sqrt()
    }

    /// Check if polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&c| c == 0)
    }

    /// Get degree of polynomial.
    /// Returns -1 for the zero polynomial.
    pub fn degree(&self) -> i32 {
        for i in (0..N).rev() {
            if self.coeffs[i] != 0 {
                return i as i32;
            }
        }
        -1
    }
}

// Conversion implementations
impl From<[i32; N]> for Polynomial {
    fn from(coeffs: [i32; N]) -> Self {
        let mut result = coeffs;
        for coeff in &mut result {
            *coeff = Self::mod_reduce(*coeff as i64);
        }
        Self { coeffs: result }
    }
}

impl From<&[i32]> for Polynomial {
    fn from(coeffs: &[i32]) -> Self {
        let mut result = [0i32; N];
        let len = coeffs.len().min(N);
        result[..len].copy_from_slice(&coeffs[..len]);

        for coeff in &mut result {
            *coeff = Self::mod_reduce(*coeff as i64);
        }
        Self { coeffs: result }
    }
}

impl From<Vec<i32>> for Polynomial {
    fn from(coeffs: Vec<i32>) -> Self {
        Self::new(coeffs)
    }
}

impl From<&Vec<i32>> for Polynomial {
    fn from(coeffs: &Vec<i32>) -> Self {
        Self::from(coeffs.as_slice())
    }
}

// Arithmetic operations
impl Add for Polynomial {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut coeffs = [0i32; N];
        for i in 0..N {
            coeffs[i] = Self::mod_reduce(
                self.coeffs[i] as i64 + other.coeffs[i] as i64,
            );
        }
        Self { coeffs }
    }
}

impl AddAssign for Polynomial {
    fn add_assign(&mut self, other: Self) {
        for i in 0..N {
            self.coeffs[i] = Self::mod_reduce(
                self.coeffs[i] as i64 + other.coeffs[i] as i64,
            );
        }
    }
}

impl Sub for Polynomial {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let mut coeffs = [0i32; N];
        for i in 0..N {
            coeffs[i] = Self::mod_reduce(
                self.coeffs[i] as i64 - other.coeffs[i] as i64,
            );
        }
        Self { coeffs }
    }
}

impl Mul<i32> for Polynomial {
    type Output = Self;

    fn mul(self, scalar: i32) -> Self {
        let mut coeffs = [0i32; N];
        for i in 0..N {
            coeffs[i] = Self::mod_reduce(self.coeffs[i] as i64 * scalar as i64);
        }
        Self { coeffs }
    }
}

impl Mul<Polynomial> for Polynomial {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.ntt_multiply(&other)
    }
}

impl Mul<&Polynomial> for Polynomial {
    type Output = Self;

    fn mul(self, other: &Polynomial) -> Self {
        self.ntt_multiply(other)
    }
}

impl Neg for Polynomial {
    type Output = Self;

    fn neg(self) -> Self {
        let mut coeffs = [0i32; N];
        for i in 0..N {
            coeffs[i] = Self::mod_reduce(-(self.coeffs[i] as i64));
        }
        Self { coeffs }
    }
}

impl Default for Polynomial {
    fn default() -> Self {
        Self { coeffs: [0i32; N] }
    }
}

#[cfg(test)]
mod prop_tests {
    use super::*;
    use quickcheck::{Arbitrary, Gen, TestResult};
    use quickcheck_macros::quickcheck;

    #[test]
    fn test_add_assign_basic() {
        let mut p1 = Polynomial::from(vec![1, 2, 3]);
        let p2 = Polynomial::from(vec![4, 5, 6]);

        p1 += p2;

        assert_eq!(p1.coeffs[0], 5);
        assert_eq!(p1.coeffs[1], 7);
        assert_eq!(p1.coeffs[2], 9);
    }

    #[test]
    fn test_polynomial_new_exceeding_n() {
        // Test 1: Basic case with N+1 coefficients
        // X^N = -1 in the ring R[X]/(X^N + 1)
        // So coefficient at position N should be subtracted from coefficient at position 0
        {
            let mut coeffs = vec![0; N + 1];
            coeffs[0] = 5;
            coeffs[N] = 3; // This represents 3*X^N = -3

            let poly = Polynomial::new(coeffs);
            assert_eq!(poly.coeffs[0], 2); // 5 - 3 = 2
            assert_eq!(poly.coeffs[1], 0);
        }

        // Test 2: Coefficients at 2N (even multiple of N)
        // X^(2N) = (X^N)^2 = (-1)^2 = 1
        // So coefficient at position 2N should be added to coefficient at position 0
        {
            let mut coeffs = vec![0; 2 * N + 1];
            coeffs[0] = 10;
            coeffs[2 * N] = 7; // This represents 7*X^(2N) = 7

            let poly = Polynomial::new(coeffs);
            assert_eq!(poly.coeffs[0], 17); // 10 + 7 = 17
        }

        // Test 3: Coefficients at 3N (odd multiple of N)
        // X^(3N) = (X^N)^3 = (-1)^3 = -1
        // So coefficient at position 3N should be subtracted from coefficient at position 0
        {
            let mut coeffs = vec![0; 3 * N + 1];
            coeffs[0] = 20;
            coeffs[3 * N] = 8; // This represents 8*X^(3N) = -8

            let poly = Polynomial::new(coeffs);
            assert_eq!(poly.coeffs[0], 12); // 20 - 8 = 12
        }

        // Test 4: Multiple overlapping coefficients
        // Testing reduction at various positions
        {
            let mut coeffs = vec![0; 2 * N + 5];
            coeffs[0] = 100;
            coeffs[1] = 50;
            coeffs[N] = 30; // -30 at position 0
            coeffs[N + 1] = 40; // -40 at position 1
            coeffs[2 * N] = 20; // +20 at position 0
            coeffs[2 * N + 1] = 10; // +10 at position 1

            let poly = Polynomial::new(coeffs);
            assert_eq!(poly.coeffs[0], 90); // 100 - 30 + 20 = 90
            assert_eq!(poly.coeffs[1], 20); // 50 - 40 + 10 = 20
        }

        // Test 5: Negative coefficients with reduction
        {
            let mut coeffs = vec![0; N + 2];
            coeffs[0] = -50;
            coeffs[N] = -30; // Represents -30*X^N = 30

            let poly = Polynomial::new(coeffs);
            // -50 - (-30) = -50 + 30 = -20
            // Then mod_reduce ensures it's in [0, Q)
            assert_eq!(poly.coeffs[0], Q - 20);
        }

        // Test 6: Large coefficients requiring modular reduction
        {
            let mut coeffs = vec![0; N + 1];
            coeffs[0] = Q + 100;
            coeffs[N] = Q + 50; // Will be subtracted

            let poly = Polynomial::new(coeffs);
            // (Q + 100) - (Q + 50) = 50
            assert_eq!(poly.coeffs[0], 50);
        }

        // Test 7: Complex case with multiple positions and reductions
        {
            let mut coeffs = vec![0; 3 * N + 10];
            // Set up various coefficients that will map to the same positions
            for i in 0..10 {
                coeffs[i] = i as i32 * 10;
                coeffs[N + i] = i as i32 * 5; // Odd power - subtract
                coeffs[2 * N + i] = i as i32 * 3; // Even power - add
                coeffs[3 * N + i] = i as i32 * 2; // Odd power - subtract
            }

            let poly = Polynomial::new(coeffs);

            // Verify first few coefficients
            // Position 0: 0 - 0 + 0 - 0 = 0
            assert_eq!(poly.coeffs[0], 0);

            // Position 1: 10 - 5 + 3 - 2 = 6
            assert_eq!(poly.coeffs[1], 6);

            // Position 2: 20 - 10 + 6 - 4 = 12
            assert_eq!(poly.coeffs[2], 12);

            // Position 3: 30 - 15 + 9 - 6 = 18
            assert_eq!(poly.coeffs[3], 18);
        }

        // Test 8: Edge case with exactly N coefficients (no reduction needed)
        {
            let coeffs: Vec<i32> = (0..N as i32).collect();
            let poly = Polynomial::new(coeffs.clone());

            for i in 0..N {
                assert_eq!(poly.coeffs[i], i as i32);
            }
        }

        // Test 9: Very large coefficient vector
        {
            let mut coeffs = vec![1; 5 * N];
            let poly = Polynomial::new(coeffs);

            // Each position should have accumulated:
            // 1 (original) - 1 (from N) + 1 (from 2N) - 1 (from 3N) + 1 (from 4N) = 1
            for i in 0..N {
                assert_eq!(poly.coeffs[i], 1);
            }
        }

        // Test 10: Mixed positive and negative with wraparound
        {
            let mut coeffs = vec![0; 2 * N + 3];
            coeffs[N - 1] = 1000;
            coeffs[2 * N - 1] = 2000; // Position 2N-1 has quotient 1 (odd), so SUBTRACTS
            coeffs[N] = -500; // Odd power, subtracts from position 0
            coeffs[0] = 100;

            let poly = Polynomial::new(coeffs);

            // Position N-1: 1000 - 2000 = -1000
            // -1000 mod Q = Q - 1000 = 8379417
            assert_eq!(poly.coeffs[N - 1], Q - 1000);

            // For position 0: 100 - (-500) = 100 + 500 = 600
            assert_eq!(poly.coeffs[0], 600);
        }
    }

    #[test]
    fn test_polynomial_new_modular_arithmetic() {
        // Test that all coefficients are properly reduced to [0, Q)
        let mut coeffs = vec![0; 2 * N];

        // Set various coefficients that need reduction
        coeffs[0] = Q + 1;
        coeffs[1] = 2 * Q + 5;
        coeffs[2] = -10;
        coeffs[3] = -Q - 20;
        coeffs[N] = Q; // Will be subtracted from position 0

        let poly = Polynomial::new(coeffs);

        // Verify all coefficients are in valid range [0, Q)
        for &coeff in &poly.coeffs {
            assert!(
                coeff >= 0 && coeff < Q,
                "Coefficient {} is out of range [0, {})",
                coeff,
                Q
            );
        }

        // Check specific values after reduction
        // Position 0: (Q + 1) - Q = 1
        assert_eq!(poly.coeffs[0], 1);

        // Position 1: (2 * Q + 5) mod Q = 5
        assert_eq!(poly.coeffs[1], 5);

        // Position 2: -10 mod Q = Q - 10
        assert_eq!(poly.coeffs[2], Q - 10);

        // Position 3: (-Q - 20) mod Q = Q - 20
        assert_eq!(poly.coeffs[3], Q - 20);
    }

    // Implement Arbitrary for Polynomial to generate random test cases
    impl Arbitrary for Polynomial {
        fn arbitrary(g: &mut Gen) -> Self {
            // Generate polynomials with varying strategies
            let strategy = u8::arbitrary(g) % 4;

            match strategy {
                0 => {
                    // Zero polynomial
                    Polynomial::zero()
                }
                1 => {
                    // Sparse polynomial (few non-zero coefficients)
                    let mut coeffs = [0i32; N];
                    let num_nonzero = usize::arbitrary(g) % 10 + 1;
                    for _ in 0..num_nonzero {
                        let idx = usize::arbitrary(g) % N;
                        coeffs[idx] = (i32::arbitrary(g) % Q).abs();
                    }
                    Polynomial::from(coeffs)
                }
                2 => {
                    // Dense polynomial with small coefficients
                    let mut coeffs = [0i32; N];
                    for i in 0..N {
                        if bool::arbitrary(g) {
                            coeffs[i] = (i32::arbitrary(g) % 1000).abs();
                        }
                    }
                    Polynomial::from(coeffs)
                }
                _ => {
                    // Fully random polynomial
                    let mut coeffs = [0i32; N];
                    for i in 0..N {
                        coeffs[i] = (i32::arbitrary(g) % Q).abs();
                    }
                    Polynomial::from(coeffs)
                }
            }
        }

        fn shrink(&self) -> Box<dyn Iterator<Item = Self>> {
            // Shrink by reducing coefficient values and making polynomial sparser
            let coeffs = self.coeffs;
            let shrunk: Vec<Polynomial> = (0..N)
                .filter(|&i| coeffs[i] != 0)
                .flat_map(move |i| {
                    let mut new_coeffs = coeffs;
                    vec![
                        // Set coefficient to zero
                        {
                            new_coeffs[i] = 0;
                            Polynomial::from(new_coeffs)
                        },
                        // Halve the coefficient
                        {
                            new_coeffs[i] = coeffs[i] / 2;
                            Polynomial::from(new_coeffs)
                        },
                    ]
                })
                .collect();

            Box::new(shrunk.into_iter())
        }
    }

    // Property: Zero is the additive identity
    #[quickcheck]
    fn prop_addition_identity(p: Polynomial) -> bool {
        let zero = Polynomial::zero();
        let sum1 = p + zero;
        let sum2 = zero + p;
        sum1 == p && sum2 == p
    }

    // Property: Every polynomial has an additive inverse
    #[quickcheck]
    fn prop_addition_inverse(p: Polynomial) -> bool {
        let neg_p = -p;
        let sum = p + neg_p;
        sum.is_zero()
    }

    // Property: Subtraction is equivalent to adding the negative
    #[quickcheck]
    fn prop_subtraction_as_addition(p1: Polynomial, p2: Polynomial) -> bool {
        let diff = p1 - p2;
        let sum = p1 + (-p2);
        diff == sum
    }

    // Property: Multiplication is commutative
    #[quickcheck]
    fn prop_multiplication_commutative(
        p1: Polynomial,
        p2: Polynomial,
    ) -> TestResult {
        // Skip if polynomials have high degree to avoid slow tests
        if p1.degree() > 50 || p2.degree() > 50 {
            return TestResult::discard();
        }

        let prod1 = p1 * p2;
        let prod2 = p2 * p1;
        TestResult::from_bool(prod1 == prod2)
    }

    // Property: Scalar multiplication is distributive
    #[quickcheck]
    fn prop_scalar_multiplication_distributive(
        p1: Polynomial,
        p2: Polynomial,
        scalar: i32,
    ) -> bool {
        let scalar = scalar % Q;
        let left = (p1 + p2) * scalar;
        let right = (p1 * scalar) + (p2 * scalar);
        left == right
    }

    // Property: Degree of sum is at most max of degrees
    #[quickcheck]
    fn prop_degree_of_sum(p1: Polynomial, p2: Polynomial) -> bool {
        let sum = p1 + p2;
        let max_degree = std::cmp::max(p1.degree(), p2.degree());
        sum.degree() <= max_degree
    }

    // Property: Degree of product is at most sum of degrees (in quotient ring)
    #[quickcheck]
    fn prop_degree_of_product(p1: Polynomial, p2: Polynomial) -> TestResult {
        if p1.is_zero() || p2.is_zero() {
            return TestResult::discard();
        }

        let prod = p1 * p2;
        // In the quotient ring R[X]/(X^N + 1), degree is always < N
        TestResult::from_bool(prod.degree() < N as i32)
    }

    // Property: Norm properties
    #[quickcheck]
    fn prop_norm_infinity_triangle_inequality(
        p1: Polynomial,
        p2: Polynomial,
    ) -> bool {
        let sum = p1 + p2;
        sum.norm_infinity() <= p1.norm_infinity() + p2.norm_infinity()
    }

    // Property: Zero polynomial has zero norm
    #[quickcheck]
    fn prop_zero_norm() -> bool {
        let zero = Polynomial::zero();
        zero.norm_infinity() == 0 && zero.norm_l2() == 0.0
    }

    // Property: Negation preserves norm
    #[quickcheck]
    fn prop_negation_preserves_norm(p: Polynomial) -> bool {
        let neg_p = -p;
        p.norm_infinity() == neg_p.norm_infinity()
            && (p.norm_l2() - neg_p.norm_l2()).abs() < 1e-10
    }

    // Property: Conversion from different types preserves values
    #[quickcheck]
    fn prop_from_slice_conversion(coeffs: Vec<i32>) -> TestResult {
        if coeffs.len() > N {
            return TestResult::discard();
        }

        let p1 = Polynomial::from(coeffs.as_slice());
        let p2 = Polynomial::from(&coeffs);
        let p3 = Polynomial::from(coeffs.clone());

        TestResult::from_bool(p1 == p2 && p2 == p3)
    }

    // Property: NTT multiplication gives same result as definition
    #[quickcheck]
    fn prop_ntt_multiplication_correctness(
        p1: Polynomial,
        p2: Polynomial,
    ) -> TestResult {
        // Only test with small polynomials for efficiency
        if p1.degree() > 10 || p2.degree() > 10 {
            return TestResult::discard();
        }

        let ntt_result = p1.ntt_multiply(&p2);
        let op_result = p1 * p2;

        TestResult::from_bool(ntt_result == op_result)
    }
}
