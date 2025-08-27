use std::ops::{Add, Mul, Sub};

use super::{polynomial::Polynomial, traits::FiniteField};

use num_traits::Zero;

// Macro for convenient polynomial vector creation
///
/// Creates a `PolynomialVector` with the given polynomials. This macro provides
/// several convenient ways to construct polynomial vectors for use in Dilithium
/// and threshold signature operations.
///
/// # Examples
///
/// ## Basic Usage
///
/// Creating a vector from existing polynomials:
/// ```
/// use math::prelude::*;
///
/// let p1: Polynomial<'static, FieldElement> = poly![1, 2, 3];
/// let p2: Polynomial<'static, FieldElement>  = poly![4, 5, 6];
/// let vec: PolynomialVector<'static, Polynomial<'static, FieldElement>> = poly_vec!(p1, p2);
///
/// assert_eq!(vec.len(), 2);
/// assert_eq!(vec.get(0), Some(&p1));
/// assert_eq!(vec.get(1), Some(&p2));
/// ```
///
/// ## Creating Zero Vector
///
/// Creating a zero polynomial vector of specified length:
/// ```
/// use math::prelude::*;
///
/// // Create zero vector with 3 polynomials
/// let zero_vec = poly_vec![0; 3];
/// assert_eq!(zero_vec.len(), 3);
/// assert!(zero_vec.get(0).unwrap().is_zero());
/// assert!(zero_vec.get(1).unwrap().is_zero());
/// assert!(zero_vec.get(2).unwrap().is_zero());
///
/// // Empty vector
/// let empty = poly_vec![];
/// assert!(empty.is_empty());
/// ```
///
/// ## From Polynomial Expressions
///
/// Creating vectors with polynomial expressions:
/// ```
/// use math::prelude::*;
///
/// // Using poly! macro inline
/// let vec1 = poly_vec![
///     poly![1, 2, 3],
///     poly![4, 5, 6],
///     poly![7, 8, 9]
/// ];
/// assert_eq!(vec1.len(), 3);
///
/// // With trailing comma
/// let vec2 = poly_vec![
///     poly![10, 20],
///     poly![30, 40],
/// ];
/// assert_eq!(vec2.len(), 2);
///
/// // Mixed polynomial types
/// let p1 = poly![1; 5];  // 5 ones
/// let p2 = poly![];      // zero polynomial
/// let vec3 = poly_vec![p1, p2, poly![1, 2, 3]];
/// assert_eq!(vec3.len(), 3);
/// ```
///
/// ## Creating Repeated Polynomials
///
/// Creating a vector with repeated polynomials:
/// ```
/// use math::prelude::*;
///
/// // Repeat the same polynomial 4 times
/// let p = poly![42, 17, 99];
/// let vec = poly_vec![p; 4];
///
/// assert_eq!(vec.len(), 4);
/// for i in 0..4 {
///     assert_eq!(vec.get(i), Some(&p));
/// }
/// ```
///
/// ## From a Vec of Polynomials
///
/// Creating from an existing vector of polynomials:
/// ```
/// use math::prelude::*;
///
/// let polys = vec![
///     poly![1, 2],
///     poly![3, 4],
///     poly![5, 6]
/// ];
/// let vec = poly_vec![polys];
/// assert_eq!(vec.len(), 3);
/// ```
///
/// ## Complex Examples
///
/// More complex usage patterns:
/// ```
/// use math::prelude::*;
///
/// // Using expressions to create polynomials
/// let a = 5;
/// let vec1 = poly_vec![
///     poly![a, a*2, a*3],
///     poly![a+1, a+2, a+3],
///     poly![-a, -2*a, -3*a]  // Will be reduced mod Q
/// ];
/// assert_eq!(vec1.len(), 3);
///
/// // Creating from function calls
/// fn make_poly(start: i32) -> Polynomial {
///     poly![start, start+1, start+2]
/// }
///
/// let vec2 = poly_vec![
///     make_poly(10),
///     make_poly(20),
///     make_poly(30)
/// ];
/// assert_eq!(vec2.len(), 3);
/// assert_eq!(vec2.get(0).unwrap().coeffs()[0], 10);
/// assert_eq!(vec2.get(1).unwrap().coeffs()[0], 20);
/// assert_eq!(vec2.get(2).unwrap().coeffs()[0], 30);
///
/// // Creating uniform vector with closures
/// let vec3 = poly_vec![(0..5).map(|i| poly![i, i*i, i*i*i]).collect::<Vec<_>>()];
/// assert_eq!(vec3.len(), 5);
/// ```
///
/// ## Usage in Algorithms
///
///
/// // Creating commitment vectors
/// let commitments = poly_vec![0; 5];  // 5 zero polynomials
/// assert_eq!(commitments.len(), 5);
///
/// // Vector arithmetic
/// let v1 = poly_vec![poly![1, 2], poly![3, 4]];
/// let v2 = poly_vec![poly![5, 6], poly![7, 8]];
/// let sum = v1 + v2;
/// assert_eq!(sum.get(0).unwrap().coeffs()[0], 6);  // 1 + 5
/// assert_eq!(sum.get(1).unwrap().coeffs()[0], 10); // 3 + 7
/// ```
#[macro_export]
macro_rules! poly_vec {
    // Empty case - empty polynomial vector
    () => {
        $crate::poly_vector::PolynomialVector::new(vec![])
    };

    //Single expression that evaluates to Vec<Polynomial>
    ($vec:expr) => {
        $crate::poly_vector::PolynomialVector::new($vec)
    };

    // Repeated zero polynomial case: poly_vec![0; count]
    (0; $count:expr) => {
        $crate::poly_vector::PolynomialVector::zero($count)
    };

    // Repeated polynomial case: poly_vec![poly; count]
    ($poly:expr; $count:expr) => {{
        let p = $poly;
        $crate::poly_vector::PolynomialVector::new(vec![p; $count])
    }};

    // List of polynomials
    ($($poly:expr),+ $(,)?) => {
        $crate::poly_vector::PolynomialVector::new(vec![$($poly),+])
    };
}

/// Represents a vector of polynomials in Rq.

/// Used for representing keys and intermediate values in Dilithium.
#[derive(Clone, Debug, PartialEq)]
pub struct PolynomialVector<'coeffs, FF: FiniteField> {
    polys: Vec<Polynomial<'coeffs, FF>>,
}

impl<FF: FiniteField> PolynomialVector<'static, FF> {
    /// Initialize polynomial vector.
    pub fn new(polynomials: Vec<Polynomial<'static, FF>>) -> Self {
        Self { polys: polynomials }
    }

    /// Get a slice of the underlying polynomial data
    pub fn as_slice(&self) -> &[Polynomial<'static, FF>] {
        &self.polys
    }

    /// Get polynomial at index.
    pub fn get(&self, index: usize) -> Option<&Polynomial<'static, FF>> {
        self.polys.get(index)
    }

    /// Get mutable polynomial at index.
    pub fn get_mut(
        &mut self,
        index: usize,
    ) -> Option<&mut Polynomial<'static, FF>> {
        self.polys.get_mut(index)
    }

    /// Set polynomial at index.
    pub fn set(
        &mut self,
        index: usize,
        poly: Polynomial<'static, FF>,
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

    /// Create zero vector of given length.
    pub fn zero(length: usize) -> Self {
        let zero_poly = Polynomial::<FF>::zero();
        Self {
            polys: vec![zero_poly; length],
        }
    }

    // /// Generate random polynomial vector.
    // pub fn random(length: usize, bound: i32) -> Self {
    //     Self {
    //         polys: (0..length).map(|_| Polynomial::random(bound)).collect(),
    //     }
    // }

    pub fn norm_infinity(&self) -> u32 {
        self.polys
            .iter()
            .map(|p| p.norm_infinity())
            .max()
            .unwrap_or(0)
    }
}

impl<FF: FiniteField> Add for PolynomialVector<'static, FF> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        assert!(self.len() == other.len());

        let polys: Vec<Polynomial<'static, FF>> = self
            .polys
            .into_iter()
            .zip(other.polys)
            .map(|(p1, p2)| p1 + p2)
            .collect();

        Self { polys }
    }
}

impl<FF: FiniteField> Add<&PolynomialVector<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        assert!(
            self.len() == other.len(),
            "Vector lengths  must match for addition"
        );

        let polys: Vec<Polynomial<'static, FF>> = self
            .polys
            .into_iter()
            .zip(other.polys.clone())
            .map(|(p1, p2)| p1 + p2)
            .collect();

        Self { polys }
    }
}

impl<FF: FiniteField> Sub for PolynomialVector<'static, FF> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        assert!(
            self.len() == other.len(),
            "Vector lengths  must match for subtraction"
        );

        let polys: Vec<Polynomial<'static, FF>> = self
            .polys
            .into_iter()
            .zip(other.polys)
            .map(|(p1, p2)| p1 - p2)
            .collect();

        Self { polys }
    }
}

impl<FF: FiniteField> Sub<&PolynomialVector<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        assert!(
            self.len() == other.len(),
            "Vector lengths  must match for subtraction"
        );

        let polys: Vec<Polynomial<'_, FF>> = self
            .polys
            .into_iter()
            .zip(other.polys.clone())
            .map(|(p1, p2)| p1 - p2)
            .collect();

        Self { polys }
    }
}
// TODO add proper testing
impl<FF: FiniteField> Mul for PolynomialVector<'static, FF> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        assert!(
            self.len() == other.len(),
            "Vector lengths must match for element-wise multiplication"
        );

        let polys: Vec<Polynomial<'static, FF>> = self
            .polys
            .into_iter()
            .zip(other.polys)
            .map(|(p1, p2)| p1 * p2)
            .collect();

        Self { polys }
    }
}

impl<FF: FiniteField> Mul<u64> for PolynomialVector<'static, FF> {
    type Output = PolynomialVector<'static, FF>;

    fn mul(self, scalar: u64) -> Self::Output {
        let scalar_ff = FF::from(scalar);

        let polys = self
            .polys
            .iter()
            .map(|p| {
                let scaled_coeffs: Vec<FF> =
                    p.coefficients().iter().map(|&c| c * scalar_ff).collect();
                Polynomial::new(scaled_coeffs)
            })
            .collect();

        PolynomialVector { polys }
    }
}

impl<FF: FiniteField> Mul<PolynomialVector<'static, FF>>
    for &Vec<Vec<Polynomial<'static, FF>>>
{
    type Output = PolynomialVector<'static, FF>;

    fn mul(self, vector: PolynomialVector<'static, FF>) -> Self::Output {
        matrix_vector_multiply(self, &vector)
    }
}

impl<FF: FiniteField> Mul<&PolynomialVector<'static, FF>>
    for &Vec<Vec<Polynomial<'static, FF>>>
{
    type Output = PolynomialVector<'static, FF>;

    fn mul(self, vector: &PolynomialVector<'static, FF>) -> Self::Output {
        matrix_vector_multiply(self, vector)
    }
}

pub fn matrix_vector_multiply<FF: FiniteField>(
    m: &[Vec<Polynomial<'static, FF>>],
    v: &PolynomialVector<'static, FF>,
) -> PolynomialVector<'static, FF> {
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
                Polynomial::<FF>::zero(),
                |mut acc, (a_ij, v_j)| {
                    // TODO make &Polynomial be multiplied &POlynomial, remove .clone()
                    acc += a_ij.clone() * v_j.clone();
                    acc
                },
            )
        })
        .collect();

    poly_vec!(result)
}

impl<FF: FiniteField> Mul<Polynomial<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn mul(self, poly: Polynomial<'static, FF>) -> Self {
        Self {
            polys: self.polys.into_iter().map(|p| p * poly.clone()).collect(),
        }
    }
}

impl<FF: FiniteField> Mul<&Polynomial<'static, FF>>
    for PolynomialVector<'static, FF>
{
    type Output = Self;

    fn mul(self, poly: &Polynomial<'static, FF>) -> Self {
        Self {
            polys: self.polys.into_iter().map(|p| p * poly.clone()).collect(),
        }
    }
}

//TODO add test coeverage
