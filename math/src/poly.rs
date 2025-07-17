/// Polynomial operations in the ring Rq = Zq[X]/(X^256 + 1).
///
/// This module implements polynomial arithmetic operations needed for the
/// Dilithium algorithm, including addition, multiplication, and reduction
/// operations in the quotient ring.

use std::ops::{Add, Sub, Mul, Neg};
use rand::Rng;

// Constants (these would typically be in a separate constants module)
const Q: i32 = 8380417; // Dilithium's prime modulus
const N: usize = 256;   // Polynomial degree bound

/// Represents a polynomial in Rq = Zq[X]/(X^256 + 1).
/// 
/// Coefficients are stored as a vector of integers modulo Q.
#[derive(Clone, Debug)]
pub struct Polynomial {
    coeffs: Vec<i32>,
}

impl Polynomial {
    /// Initialize polynomial with given coefficients.
    pub fn new(coeffs: Vec<i32>) -> Self {
        let mut poly = if coeffs.len() > N {
            // Reduce modulo X^N + 1
            Self {
                coeffs: Self::reduce_mod_xn_plus_1(&coeffs),
            }
        } else if coeffs.len() < N {
            // Pad with zeros
            let mut padded = vec![0; N];
            padded[..coeffs.len()].copy_from_slice(&coeffs);
            Self { coeffs: padded }
        } else {
            Self { coeffs }
        };
        
        // Reduce coefficients modulo Q
        for coeff in &mut poly.coeffs {
            *coeff = coeff.rem_euclid(Q);
        }
        
        poly
    }
    
    /// Reduce polynomial modulo X^N + 1.
    fn reduce_mod_xn_plus_1(coeffs: &[i32]) -> Vec<i32> {
        let mut result = vec![0; N];
        
        for (i, &coeff) in coeffs.iter().enumerate() {
            let pos = i % N;
            if i >= N && (i / N) % 2 == 1 {
                // X^N = -1, so X^(N+k) = -X^k
                result[pos] = (result[pos] - coeff).rem_euclid(Q);
            } else {
                result[pos] = (result[pos] + coeff).rem_euclid(Q);
            }
        }
        
        result
    }
    
    /// Multiply two polynomials in Rq.
    /// 
    /// This is a naive O(N^2) implementation. For better performance,
    /// use NTT-based multiplication.
    fn poly_multiply(&self, other: &Polynomial) -> Polynomial {
        let mut result_coeffs = vec![0i64; 2 * N - 1];
        
        for (i, &a) in self.coeffs.iter().enumerate() {
            for (j, &b) in other.coeffs.iter().enumerate() {
                result_coeffs[i + j] += (a as i64) * (b as i64);
            }
        }
        
        // Convert back to i32 after reducing modulo Q
        let reduced_coeffs: Vec<i32> = result_coeffs
            .iter()
            .map(|&c| (c % Q as i64) as i32)
            .collect();
        
        Polynomial::new(reduced_coeffs)
    }
    
    /// Compute infinity norm of polynomial.
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
        let sum: i64 = self.coeffs
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
    pub fn degree(&self) -> i32 {
        if self.is_zero() {
            return -1;
        }
        
        for i in (0..N).rev() {
            if self.coeffs[i] != 0 {
                return i as i32;
            }
        }
        -1
    }
    
    /// Create zero polynomial.
    pub fn zero() -> Self {
        Self {
            coeffs: vec![0; N],
        }
    }
    
    /// Create polynomial representing 1.
    pub fn one() -> Self {
        let mut coeffs = vec![0; N];
        coeffs[0] = 1;
        Self { coeffs }
    }
    
    /// Generate random polynomial with coefficients in [0, bound).
    pub fn random(bound: i32) -> Self {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<i32> = (0..N)
            .map(|_| rng.gen_range(0..bound))
            .collect();
        Self::new(coeffs)
    }
    
    /// Create a copy of the polynomial.
    pub fn copy(&self) -> Self {
        self.clone()
    }
}

impl PartialEq for Polynomial {
    fn eq(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }
}

impl Eq for Polynomial {}

impl Add for Polynomial {
    type Output = Self;
    
    fn add(self, other: Self) -> Self {
        let coeffs: Vec<i32> = self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| (a + b).rem_euclid(Q))
            .collect();
        
        Self { coeffs }
    }
}

impl Sub for Polynomial {
    type Output = Self;
    
    fn sub(self, other: Self) -> Self {
        let coeffs: Vec<i32> = self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| (a - b).rem_euclid(Q))
            .collect();
        
        Self { coeffs }
    }
}

impl Mul<i32> for Polynomial {
    type Output = Self;
    
    fn mul(self, scalar: i32) -> Self {
        let coeffs: Vec<i32> = self.coeffs
            .iter()
            .map(|&c| ((c as i64 * scalar as i64) % Q as i64) as i32)
            .collect();
        
        Self { coeffs }
    }
}

impl Mul<Polynomial> for Polynomial {
    type Output = Self;
    
    fn mul(self, other: Polynomial) -> Self {
        self.poly_multiply(&other)
    }
}

impl Neg for Polynomial {
    type Output = Self;
    
    fn neg(self) -> Self {
        let coeffs: Vec<i32> = self.coeffs
            .iter()
            .map(|&c| (-c).rem_euclid(Q))
            .collect();
        
        Self { coeffs }
    }
}

/// Represents a vector of polynomials in Rq.
/// 
/// Used for representing keys and intermediate values in Dilithium.
#[derive(Clone, Debug)]
pub struct PolynomialVector {
    polys: Vec<Polynomial>,
    length: usize,
}

impl PolynomialVector {
    /// Initialize polynomial vector.
    pub fn new(polynomials: Vec<Polynomial>) -> Self {
        let length = polynomials.len();
        Self {
            polys: polynomials,
            length,
        }
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
        if index >= self.length {
            return Err("Index out of bounds");
        }
        self.polys[index] = poly;
        Ok(())
    }
    
    /// Get length of vector.
    pub fn len(&self) -> usize {
        self.length
    }
    
    /// Check if vector is empty.
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }
    
    /// Compute infinity norm of vector.
    pub fn norm_infinity(&self) -> i32 {
        self.polys
            .iter()
            .map(|p| p.norm_infinity())
            .max()
            .unwrap_or(0)
    }
    
    /// Create a copy of the vector.
    pub fn copy(&self) -> Self {
        self.clone()
    }
    
    /// Create zero vector of given length.
    pub fn zero(length: usize) -> Self {
        let polys = vec![Polynomial::zero(); length];
        Self { polys, length }
    }
    
    /// Generate random polynomial vector.
    pub fn random(length: usize, bound: i32) -> Self {
        let polys: Vec<Polynomial> = (0..length)
            .map(|_| Polynomial::random(bound))
            .collect();
        Self { polys, length }
    }
}

impl PartialEq for PolynomialVector {
    fn eq(&self, other: &Self) -> bool {
        self.length == other.length && self.polys == other.polys
    }
}

impl Eq for PolynomialVector {}

impl Add for PolynomialVector {
    type Output = Result<Self, &'static str>;
    
    fn add(self, other: Self) -> Self::Output {
        if self.length != other.length {
            return Err("Vector lengths must match");
        }
        
        let polys: Vec<Polynomial> = self.polys
            .into_iter()
            .zip(other.polys.into_iter())
            .map(|(p1, p2)| p1 + p2)
            .collect();
        
        Ok(Self {
            polys,
            length: self.length,
        })
    }
}

impl Sub for PolynomialVector {
    type Output = Result<Self, &'static str>;
    
    fn sub(self, other: Self) -> Self::Output {
        if self.length != other.length {
            return Err("Vector lengths must match");
        }
        
        let polys: Vec<Polynomial> = self.polys
            .into_iter()
            .zip(other.polys.into_iter())
            .map(|(p1, p2)| p1 - p2)
            .collect();
        
        Ok(Self {
            polys,
            length: self.length,
        })
    }
}

impl Mul<i32> for PolynomialVector {
    type Output = Self;
    
    fn mul(self, scalar: i32) -> Self {
        let polys: Vec<Polynomial> = self.polys
            .into_iter()
            .map(|p| p * scalar)
            .collect();
        
        Self {
            polys,
            length: self.length,
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
    fn test_vector_operations() {
        let v1 = PolynomialVector::random(3, 100);
        let v2 = PolynomialVector::random(3, 100);
        
        let v3 = (v1.clone() + v2.clone()).unwrap();
        assert_eq!(v3.len(), 3);
        
        let v4 = v1 * 5;
        assert_eq!(v4.len(), 3);
    }
}