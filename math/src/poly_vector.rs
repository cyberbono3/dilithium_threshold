use std::ops::{Add, Mul, Sub};

use super::polynomial::Polynomial;

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

    /// Get a slice of the underlying polynomial data
    pub fn as_slice(&self) -> &[Polynomial] {
        &self.polys
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
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        assert!(self.len() == other.len());

        let polys: Vec<Polynomial> = self
            .polys
            .into_iter()
            .zip(other.polys)
            .map(|(p1, p2)| p1 + p2)
            .collect();

        Self { polys }
    }
}

impl Add<&PolynomialVector> for PolynomialVector {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        assert!(
            self.len() == other.len(),
            "Vector lengths  must match for addition"
        );

        let polys: Vec<Polynomial> = self
            .polys
            .into_iter()
            .zip(other.polys.clone())
            .map(|(p1, p2)| p1 + p2)
            .collect();

        Self { polys }
    }
}

impl Sub for PolynomialVector {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        assert!(
            self.len() == other.len(),
            "Vector lengths  must match for subtraction"
        );

        let polys: Vec<Polynomial> = self
            .polys
            .into_iter()
            .zip(other.polys)
            .map(|(p1, p2)| p1 - p2)
            .collect();

        Self { polys }
    }
}

impl Sub<&PolynomialVector> for PolynomialVector {
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        assert!(
            self.len() == other.len(),
            "Vector lengths  must match for subtraction"
        );

        let polys: Vec<Polynomial> = self
            .polys
            .into_iter()
            .zip(other.polys.clone())
            .map(|(p1, p2)| p1 - p2)
            .collect();

        Self { polys }
    }
}
// TODO add proper testing
impl Mul for PolynomialVector {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        assert!(
            self.len() == other.len(),
            "Vector lengths must match for element-wise multiplication"
        );

        let polys: Vec<Polynomial> = self
            .polys
            .into_iter()
            .zip(other.polys)
            .map(|(p1, p2)| p1 * p2)
            .collect();

        Self { polys }
    }
}

impl Mul<PolynomialVector> for &Vec<Vec<Polynomial>> {
    type Output = PolynomialVector;

    fn mul(self, vector: PolynomialVector) -> Self::Output {
        matrix_vector_multiply(self, &vector)
    }
}

impl Mul<&PolynomialVector> for &Vec<Vec<Polynomial>> {
    type Output = PolynomialVector;

    fn mul(self, vector: &PolynomialVector) -> Self::Output {
        matrix_vector_multiply(self, vector)
    }
}

pub fn matrix_vector_multiply(
    m: &[Vec<Polynomial>],
    v: &PolynomialVector,
) -> PolynomialVector {
    // Check matrix is not empty
    assert!(!m.is_empty(), "Matrix cannot be empty");

    let cols = m[0].len();

    // Check matrix columns match vector length
    assert_eq!(
        cols,
        v.len(),
        "Matrix columns ({}) must match vector length ({})",
        cols,
        v.len()
    );

    // Check all rows have same length (rectangular matrix)
    assert!(
        m.iter().all(|row| row.len() == cols),
        "All matrix rows must have the same length"
    );

    let result = m
        .iter()
        .map(|row| {
            row.iter().zip(v.as_slice()).fold(
                Polynomial::zero(),
                |mut acc, (a_ij, v_j)| {
                    acc += *a_ij * v_j;
                    acc
                },
            )
        })
        .collect();

    PolynomialVector::new(result)
}

impl Mul<i32> for PolynomialVector {
    type Output = Self;

    fn mul(self, scalar: i32) -> Self {
        Self {
            polys: self.polys.into_iter().map(|p| p * scalar).collect(),
        }
    }
}

impl Mul<Polynomial> for PolynomialVector {
    type Output = Self;

    fn mul(self, poly: Polynomial) -> Self {
        Self {
            polys: self.polys.into_iter().map(|p| p * poly).collect(),
        }
    }
}

impl Mul<&Polynomial> for PolynomialVector {
    type Output = Self;

    fn mul(self, poly: &Polynomial) -> Self {
        Self {
            polys: self.polys.into_iter().map(|p| p * poly).collect(),
        }
    }
}
