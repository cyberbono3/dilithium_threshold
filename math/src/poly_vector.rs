use std::ops::{Add, Mul, Sub};

use super::poly::Polynomial;

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
    pub fn set(
        &mut self,
        index: usize,
        poly: Polynomial,
    ) -> Result<(), &'static str> {
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
