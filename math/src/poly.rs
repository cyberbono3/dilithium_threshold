use rand::Rng;
use std::ops::{Add, Mul, Neg, Sub};

use super::ntt::{intt, montgomery_reduce, ntt};

/// Constants for the Dilithium algorithm
pub const Q: i32 = 8380417; // Dilithium's prime modulus
pub const N: usize = 256; // Polynomial degree bound
const ROOT_OF_UNITY: i32 = 1753; // primitive 256th root of unity mod Q

/// Simple modular reduction ensuring result is in [0, Q)
#[inline(always)]
fn mod_reduce(a: i64) -> i32 {
    let mut t = (a % (Q as i64)) as i32;
    t += (t >> 31) & Q; // Add Q if t is negative
    t
}

/// Compute a^b mod Q using fast exponentiation
fn pow_mod(a: i32, mut b: u32) -> i32 {
    let mut result = 1i64;
    let a_i64 = if a < 0 { a as i64 + Q as i64 } else { a as i64 };
    let mut base = a_i64;

    while b > 0 {
        if b & 1 == 1 {
            result = (result * base) % (Q as i64);
        }
        base = (base * base) % (Q as i64);
        b >>= 1;
    }
    result as i32
}

/// Compute modular inverse using Fermat's little theorem
fn mod_inverse(a: i32) -> i32 {
    pow_mod(a, (Q - 2) as u32)
}

/// Represents a polynomial in Rq = Zq[X]/(X^256 + 1).
///
/// Coefficients are stored as an array of integers modulo Q.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polynomial {
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
                    result[pos] = mod_reduce(result[pos] as i64 + coeff as i64);
                } else {
                    // Odd powers of X^N contribute negatively (since X^N = -1)
                    result[pos] = mod_reduce(result[pos] as i64 - coeff as i64);
                }
            }
        }

        // Ensure all coefficients are in [0, Q)
        for coeff in &mut result {
            *coeff = mod_reduce(*coeff as i64);
        }

        Self { coeffs: result }
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
            *coeff = mod_reduce(*coeff as i64);
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
            *coeff = mod_reduce(*coeff as i64);
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
            coeffs[i] = mod_reduce(self.coeffs[i] as i64 + other.coeffs[i] as i64);
        }
        Self { coeffs }
    }
}

impl Sub for Polynomial {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let mut coeffs = [0i32; N];
        for i in 0..N {
            coeffs[i] = mod_reduce(self.coeffs[i] as i64 - other.coeffs[i] as i64);
        }
        Self { coeffs }
    }
}

impl Mul<i32> for Polynomial {
    type Output = Self;

    fn mul(self, scalar: i32) -> Self {
        let mut coeffs = [0i32; N];
        for i in 0..N {
            coeffs[i] = mod_reduce(self.coeffs[i] as i64 * scalar as i64);
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

impl Neg for Polynomial {
    type Output = Self;

    fn neg(self) -> Self {
        let mut coeffs = [0i32; N];
        for i in 0..N {
            coeffs[i] = mod_reduce(-(self.coeffs[i] as i64));
        }
        Self { coeffs }
    }
}

impl Default for Polynomial {
    fn default() -> Self {
        Self{coeffs: [0i32;N]}
    }
}

/// Represents a vector of polynomials in Rq.
///
/// Used for representing keys and intermediate values in Dilithium.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolynomialVector {
    polys: Vec<Polynomial>,
}

impl PolynomialVector {
    /// Initialize polynomial vector.
    pub fn new(polynomials: Vec<Polynomial>) -> Self {
        Self { polys: polynomials }
    }

    /// Get polynomial at index.
    pub fn get(&self, index: usize) -> Option<&Polynomial> {
        self.polys.get(index)
    }

    /// Get mutable polynomial at index.
    pub fn get_mut(&mut self, index: usize) -> Option<&mut Polynomial> {
        self.polys.get_mut(index)
    }

    /// Set polynomial at index.
    pub fn set(&mut self, index: usize, poly: Polynomial) -> Result<(), &'static str> {
        if index >= self.polys.len() {
            return Err("Index out of bounds");
        }
        self.polys[index] = poly;
        Ok(())
    }

    /// Get length of vector.
    pub fn len(&self) -> usize {
        self.polys.len()
    }

    /// Check if vector is empty.
    pub fn is_empty(&self) -> bool {
        self.polys.is_empty()
    }

    /// Compute infinity norm of vector.
    pub fn norm_infinity(&self) -> i32 {
        self.polys
            .iter()
            .map(|p| p.norm_infinity())
            .max()
            .unwrap_or(0)
    }

    /// Create zero vector of given length.
    pub fn zero(length: usize) -> Self {
        Self {
            polys: vec![Polynomial::zero(); length],
        }
    }

    /// Generate random polynomial vector.
    pub fn random(length: usize, bound: i32) -> Self {
        Self {
            polys: (0..length).map(|_| Polynomial::random(bound)).collect(),
        }
    }
}

impl Add for PolynomialVector {
    type Output = Result<Self, &'static str>;

    fn add(self, other: Self) -> Self::Output {
        if self.len() != other.len() {
            return Err("Vector lengths must match");
        }

        let polys: Vec<Polynomial> = self
            .polys
            .into_iter()
            .zip(other.polys)
            .map(|(p1, p2)| p1 + p2)
            .collect();

        Ok(Self { polys })
    }
}

impl Sub for PolynomialVector {
    type Output = Result<Self, &'static str>;

    fn sub(self, other: Self) -> Self::Output {
        if self.len() != other.len() {
            return Err("Vector lengths must match");
        }

        let polys: Vec<Polynomial> = self
            .polys
            .into_iter()
            .zip(other.polys)
            .map(|(p1, p2)| p1 - p2)
            .collect();

        Ok(Self { polys })
    }
}

impl Mul<i32> for PolynomialVector {
    type Output = Self;

    fn mul(self, scalar: i32) -> Self {
        Self {
            polys: self.polys.into_iter().map(|p| p * scalar).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_creation() {
        let p = Polynomial::new(vec![1, 2, 3]);
        assert_eq!(p.coeffs.len(), N);
        assert_eq!(p.coeffs[0], 1);
        assert_eq!(p.coeffs[1], 2);
        assert_eq!(p.coeffs[2], 3);
    }

    #[test]
    fn test_polynomial_addition() {
        let p1 = Polynomial::new(vec![1, 2, 3]);
        let p2 = Polynomial::new(vec![4, 5, 6]);
        let p3 = p1 + p2;
        assert_eq!(p3.coeffs[0], 5);
        assert_eq!(p3.coeffs[1], 7);
        assert_eq!(p3.coeffs[2], 9);
    }

    #[test]
    fn test_polynomial_multiplication() {
        let p1 = Polynomial::new(vec![1, 2]);
        let p2 = Polynomial::new(vec![3, 4]);
        let p3 = p1 * p2;
        // (1 + 2x)(3 + 4x) = 3 + 4x + 6x + 8x^2 = 3 + 10x + 8x^2
        assert_eq!(p3.coeffs[0], 3);
        assert_eq!(p3.coeffs[1], 10);
        assert_eq!(p3.coeffs[2], 8);
    }

    #[test]
    fn test_ntt_multiplication() {
        let p1 = Polynomial::new(vec![1, 2]);
        let p2 = Polynomial::new(vec![3, 4]);
        let p3 = p1.ntt_multiply(&p2);
        // (1 + 2x)(3 + 4x) = 3 + 10x + 8x^2
        assert_eq!(p3.coeffs[0], 3);
        assert_eq!(p3.coeffs[1], 10);
        assert_eq!(p3.coeffs[2], 8);
    }

    #[test]
    fn test_vector_operations() {
        let v1 = PolynomialVector::random(3, 100);
        let v2 = PolynomialVector::random(3, 100);

        let v3 = (v1.clone() + v2.clone()).unwrap();
        assert_eq!(v3.len(), 3);

        let v4 = v1 * 5;
        assert_eq!(v4.len(), 3);
    }

    #[test]
    fn test_polynomial_degree() {
        let p1 = Polynomial::zero();
        assert_eq!(p1.degree(), -1);

        let p2 = Polynomial::new(vec![1, 2, 3, 0, 0]);
        assert_eq!(p2.degree(), 2);

        let p3 = Polynomial::one();
        assert_eq!(p3.degree(), 0);
    }

    #[test]
    fn test_polynomial_norms() {
        let p = Polynomial::new(vec![1, 2, 3]);
        assert!(p.norm_infinity() > 0);
        assert!(p.norm_l2() > 0.0);
    }

    #[test]
    fn test_polynomial_negation() {
        let p1 = Polynomial::new(vec![1, 2, 3]);
        let p2 = -p1.clone();
        let p3 = p1 + p2;
        assert!(p3.is_zero());
    }
}