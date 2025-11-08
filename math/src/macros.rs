//! Shared macros for constructing core math primitives.
//!
//! These macros delegate to the types they create, which keeps the public API
//! concise and avoids duplicating builder logic across the crate.

/// Simplifies constructing [`FieldElement`](crate::field_element::FieldElement)s.
///
/// The type `FieldElement` must be in scope for this macro to work.
/// See [`FieldElement::from`](crate::field_element::FieldElement::from) for supported types.
///
/// ```
/// use math::prelude::*;
///
/// let a = fe!(42);
/// assert_eq!(a, FieldElement::from(42));
/// ```
#[macro_export]
macro_rules! fe {
    ($value:expr) => {
        $crate::field_element::FieldElement::from($value)
    };
}

/// Create a [`Vec`] of [`FieldElement`](crate::field_element::FieldElement)s.
///
/// ```
/// use math::prelude::*;
///
/// let repeated = fe_vec![7; 3];
/// assert_eq!(repeated, vec![FieldElement::from(7); 3]);
/// ```
#[macro_export]
macro_rules! fe_vec {
    ($b:expr; $n:expr) => {
        $crate::builders::repeat_field_element($b, $n)
    };
    ($($b:expr),* $(,)?) => {
        vec![$($crate::field_element::FieldElement::from($b)),*]
    };
}

/// Create an array of [`FieldElement`](crate::field_element::FieldElement)s.
///
/// ```
/// use math::prelude::*;
///
/// let arr = fe_array![1, 2, 3];
/// assert_eq!(arr, [FieldElement::from(1), FieldElement::from(2), FieldElement::from(3)]);
/// ```
#[macro_export]
macro_rules! fe_array {
    ($b:expr; $n:expr) => {
        [$crate::field_element::FieldElement::from($b); $n]
    };
    ($($b:expr),* $(,)?) => {
        [$($crate::field_element::FieldElement::from($b)),*]
    };
}

/// Construct a [`Polynomial`](crate::poly::Polynomial) from coefficients.
///
/// ```
/// use math::prelude::*;
///
/// let poly: Polynomial<'static, FieldElement> = poly![1, 2, 3];
/// assert_eq!(poly.coefficients(), &[fe!(1), fe!(2), fe!(3)]);
/// ```
#[macro_export]
macro_rules! poly {
    () => {{
        $crate::prelude::Polynomial::zero()
    }};
    ($arr:expr) => {{
        $crate::prelude::Polynomial::from($arr)
    }};
    ($val:expr; $count:expr) => {{
        let value = $val;
        let count = $count;
        $crate::prelude::Polynomial::from(vec![value; count])
    }};
    ($($coeff:expr),+ $(,)?) => {{
        $crate::prelude::Polynomial::from(vec![$($coeff),+])
    }};
}

/// Construct a [`PolynomialVector`](crate::poly_vector::PolynomialVector).
///
/// ```
/// use math::prelude::*;
///
/// let vec: PolynomialVector<'static, FieldElement> =
///     poly_vec![poly![1, 2], poly![3, 4]];
/// assert_eq!(vec.len(), 2);
/// ```
#[macro_export]
macro_rules! poly_vec {
    () => {
        $crate::poly_vector::PolynomialVector::new(vec![])
    };
    ($vec:expr) => {
        $crate::poly_vector::PolynomialVector::new($vec)
    };
    (0; $count:expr) => {
        $crate::builders::zeroed_vector($count)
    };
    ($poly:expr; $count:expr) => {{
        let p = $poly;
        $crate::builders::repeat_polynomial(p, $count)
    }};
    ($($poly:expr),+ $(,)?) => {
        $crate::poly_vector::PolynomialVector::new(vec![$($poly),+])
    };
}
