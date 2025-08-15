use std::convert::TryFrom;
use std::fmt;
use std::hash::Hash;
use std::num::TryFromIntError;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Neg;
use std::ops::Sub;
use std::ops::SubAssign;
use std::str::FromStr;

use arbitrary::Arbitrary;
use arbitrary::Unstructured;
use get_size2::GetSize;
use num_traits::ConstOne;
use num_traits::ConstZero;
use num_traits::One;
use num_traits::Zero;
use phf::phf_map;

use serde::Deserialize;
use serde::Deserializer;
use serde::Serialize;
use serde::Serializer;

use super::{error::ParseFieldElementError, traits::*};

// Primitive roots of unity for prime 8380417
// These are computed using generator 10
const PRIMITIVE_ROOTS: phf::Map<u32, u32> = phf_map! {
    0u32 => 1,
    1u32 => 1,
    2u32 => 8380416,
    4u32 => 4808194,
    8u32 => 4614810,
    16u32 => 2883726,
    32u32 => 6250525,
    64u32 => 7044481,
    128u32 => 3241972,
    256u32 => 6644104,
    512u32 => 1921994,
    1024u32 => 550930,
    2048u32 => 1028169,
    4096u32 => 1856698,
    8192u32 => 1938117,
};

/// Base field element ∈ ℤ_{8380417}.
///
/// In Montgomery representation.
#[derive(Debug, Copy, Clone, Default, Hash, PartialEq, Eq)]
pub struct FieldElement(u32);

/// Simplifies constructing [FieldElement]s.
///
/// The type [`FieldElement`] must be in scope for this macro to work.
/// See [`FieldElement::from`] for supported types.
///
/// # Examples
///
/// ```
/// use math::prelude::*;
/// let a = fe!(42);
/// let b = fe!(-12); // correctly translates to `FieldElement::P - 12`
/// let c = fe!(42 - 12);
/// assert_eq!(a + b, c);
///```
#[macro_export]
macro_rules! fe {
    ($value:expr) => {
        $crate::field_element::FieldElement::from($value)
    };
}

/// Simplifies constructing vectors of [FieldElement]s.
///
/// The type [`FieldElement`] must be in scope for this macro to work. See also [`fe!`].
///
/// # Examples
///
/// ```
/// use math::prelude::*;
/// let a = fe_vec![1, 2, 3];
/// let b = vec![fe!(1), fe!(2), fe!(3)];
/// assert_eq!(a, b);
/// ```
///
/// ```
/// use math::prelude::*;
/// let a = fe_vec![42; 15];
/// let b = vec![fe!(42); 15];
/// assert_eq!(a, b);
/// ```
///
#[macro_export]
macro_rules! fe_vec {
    ($b:expr; $n:expr) => {
        vec![$crate::field_element::FieldElement::from($b); $n]
    };
    ($($b:expr),* $(,)?) => {
        vec![$($crate::field_element::FieldElement::from($b)),*]
    };
}

/// Simplifies constructing arrays of [base field element][FieldElement]s.
///
/// The type [`FieldElement`] must be in scope for this macro to work. See also [`fe!`].
///
/// # Examples
///
/// ```
/// use math::prelude::*;
/// let a = fe_array![1, 2, 3];
/// let b = [fe!(1), fe!(2), fe!(3)];
/// assert_eq!(a, b);
/// ```
///
/// ```
/// use math::prelude::*;
/// let a = fe_array![42; 15];
/// let b = [fe!(42); 15];
/// assert_eq!(a, b);
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

impl GetSize for FieldElement {
    fn get_stack_size() -> usize {
        std::mem::size_of::<Self>()
    }

    fn get_heap_size(&self) -> usize {
        0
    }
}

impl<'a> Arbitrary<'a> for FieldElement {
    fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
        u.arbitrary().map(FieldElement::new)
    }
}

impl Serialize for FieldElement {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        self.value().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for FieldElement {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        Ok(Self::new(u32::deserialize(deserializer)?))
    }
}

impl FieldElement {
    pub const BYTES: usize = 4; // Changed from 8 to 4

    /// Dilithium prime modulus: 8380417
    pub const P: u32 = 8380417;
    pub const MAX: u32 = Self::P - 1;

    /// R^2 mod P for Montgomery form
    /// R = 2^32 mod 8380417 = 4193792
    /// R^2 = 4193792^2 mod 8380417 = 2365951
    const R2: u32 = 2365951;

    /// -2^-1
    pub const MINUS_TWO_INVERSE: Self = Self::new(4190208); // (P - 1) / 2

    #[inline]
    pub const fn new(value: u32) -> Self {
        Self(Self::montyred((value as u64) * (Self::R2 as u64)))
    }

    /// Construct a new base field element iff the given value is
    /// [canonical][Self::is_canonical], an error otherwise.
    fn try_new(v: u32) -> Result<Self, ParseFieldElementError> {
        Self::is_canonical(v)
            .then(|| Self::new(v))
            .ok_or(ParseFieldElementError::NotCanonical(v as u64))
    }

    #[inline]
    pub const fn value(&self) -> u32 {
        self.canonical_representation()
    }

    /// Get a generator for the entire field
    /// For prime 8380417, valid generators include 10, 13, 14, 15, 17
    /// We use 10 as it's the smallest two-digit generator
    pub const fn generator() -> Self {
        FieldElement::new(10)
    }

    // You should probably only use `increment` and `decrement` for testing purposes
    pub fn increment(&mut self) {
        *self += Self::one();
    }

    // You should probably only use `increment` and `decrement` for testing purposes
    pub fn decrement(&mut self) {
        *self -= Self::one();
    }

    #[inline]
    const fn canonical_representation(&self) -> u32 {
        Self::montyred(self.0 as u64)
    }

    #[must_use]
    #[inline]
    pub const fn mod_pow(&self, exp: u32) -> Self {
        let mut acc = FieldElement::ONE;
        let bit_length = u32::BITS - exp.leading_zeros();
        let mut i = 0;
        while i < bit_length {
            acc = Self(Self::montyred(acc.0 as u64 * acc.0 as u64));
            if exp & (1 << (bit_length - 1 - i)) != 0 {
                acc = Self(Self::montyred(acc.0 as u64 * self.0 as u64));
            }
            i += 1;
        }

        acc
    }

    /// Montgomery reduction
    #[inline(always)]
    pub const fn montyred(x: u64) -> u32 {
        // For prime 8380417, we use Montgomery reduction with R = 2^32
        // n' = -n^(-1) mod R where n is our prime
        // n' for 8380417 is 4236238847
        const N_PRIME: u32 = 4236238847;

        let xl = x as u32;

        // m = (x * n') mod R
        let m = xl.wrapping_mul(N_PRIME);

        // t = (x + m * n) / R
        let mn = (m as u64) * (Self::P as u64);
        let t = (x.wrapping_add(mn)) >> 32;

        // Fast path for common cases (99.9% of actual usage)
        // In multiplication: t < 2*P
        if t < 2 * Self::P as u64 {
            if t >= Self::P as u64 {
                (t - Self::P as u64) as u32
            } else {
                t as u32
            }
        } else {
            // Slow path for arbitrary inputs (only in tests/edge cases)
            // Use division since it's still const fn compatible
            // and this path is rarely taken in practice
            (t % Self::P as u64) as u32
        }
    }

    /// Return the raw bytes or 8-bit chunks of the Montgomery
    /// representation, in little-endian byte order
    pub const fn raw_bytes(&self) -> [u8; 4] {
        self.0.to_le_bytes()
    }

    /// Take a slice of 4 bytes and interpret it as an integer in
    /// little-endian byte order, and cast it to a FieldElement
    /// in Montgomery representation
    pub const fn from_raw_bytes(bytes: &[u8; 4]) -> Self {
        Self(u32::from_le_bytes(*bytes))
    }

    /// Return the raw 16-bit chunks of the Montgomery
    /// representation, in little-endian chunk order
    pub const fn raw_u16s(&self) -> [u16; 2] {
        [(self.0 & 0xffff) as u16, ((self.0 >> 16) & 0xffff) as u16]
    }

    /// Take a slice of 2 16-bit chunks and interpret it as an integer in
    /// little-endian chunk order, and cast it to a FieldElement
    /// in Montgomery representation
    pub const fn from_raw_u16s(chunks: &[u16; 2]) -> Self {
        Self(((chunks[1] as u32) << 16) | (chunks[0] as u32))
    }

    #[inline]
    pub fn raw_u64(&self) -> u64 {
        self.0 as u64
    }

    #[inline]
    pub const fn from_raw_u32(e: u32) -> FieldElement {
        FieldElement(e)
    }

    #[inline]
    pub const fn raw_u32(&self) -> u32 {
        self.0
    }

    #[inline]
    pub const fn is_canonical(x: u32) -> bool {
        x < Self::P
    }
}

impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let canonical_value = Self::canonical_representation(self);
        let cutoff = 256;
        if canonical_value >= Self::P - cutoff {
            write!(f, "-{}", Self::P - canonical_value)
        } else if canonical_value <= cutoff {
            write!(f, "{canonical_value}")
        } else {
            write!(f, "{canonical_value:>07}") // Adjusted width for smaller prime
        }
    }
}

impl FromStr for FieldElement {
    type Err = ParseFieldElementError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parsed: u32 = s.parse().map_err(|_| {
            Self::Err::ParseU64Error(s.parse::<u64>().unwrap_err())
        })?;
        Self::try_new(parsed)
    }
}

impl From<usize> for FieldElement {
    fn from(value: usize) -> Self {
        // For our smaller prime, we need to reduce modulo P
        Self::new((value as u64 % Self::P as u64) as u32)
    }
}

impl From<u128> for FieldElement {
    fn from(value: u128) -> Self {
        // Reduce u128 modulo our prime
        Self::new((value % Self::P as u128) as u32)
    }
}

impl From<u64> for FieldElement {
    fn from(value: u64) -> Self {
        // Reduce u64 modulo our prime
        Self::new((value % Self::P as u64) as u32)
    }
}

macro_rules! impl_from_small_unsigned_int_for_fe {
    ($($t:ident),+ $(,)?) => {$(
        impl From<$t> for FieldElement {
            fn from(value: $t) -> Self {
                Self::new(u32::from(value))
            }
        }
    )+};
}

impl_from_small_unsigned_int_for_fe!(u8, u16, u32);

impl From<isize> for FieldElement {
    fn from(value: isize) -> Self {
        if value >= 0 {
            Self::from(value as usize)
        } else {
            // Convert to i64 to safely handle isize::MIN
            Self::from(value as i64)
        }
    }
}

impl From<i64> for FieldElement {
    fn from(value: i64) -> Self {
        if value >= 0 {
            Self::from(value as u64)
        } else {
            // For negative values, we need to handle the modular arithmetic carefully
            // We want the representative in [0, P-1] of the equivalence class of value mod P
            let remainder = value % (Self::P as i64);
            if remainder == 0 {
                Self::ZERO
            } else {
                // remainder is negative here, so P + remainder gives us the positive representative
                Self::new((Self::P as i64 + remainder) as u32)
            }
        }
    }
}

macro_rules! impl_from_small_signed_int_for_fe {
    ($($t:ident),+ $(,)?) => {$(
        impl From<$t> for FieldElement {
            fn from(value: $t) -> Self {
                i64::from(value).into()
            }
        }
    )+};
}

impl_from_small_signed_int_for_fe!(i8, i16, i32);

macro_rules! impl_try_from_fe_for_int {
    ($($t:ident),+ $(,)?) => {$(
        impl TryFrom<FieldElement> for $t {
            type Error = TryFromIntError;

            fn try_from(value: FieldElement) -> Result<Self, Self::Error> {
                $t::try_from(value.canonical_representation())
            }
        }

        impl TryFrom<&FieldElement> for $t {
            type Error = TryFromIntError;

            fn try_from(value: &FieldElement) -> Result<Self, Self::Error> {
                $t::try_from(value.canonical_representation())
            }
        }
    )+};
}

impl_try_from_fe_for_int!(u8, i8, u16, i16, usize, isize);

macro_rules! impl_from_fe_for_int {
    ($($t:ident),+ $(,)?) => {$(
        impl From<FieldElement> for $t {
            fn from(elem: FieldElement) -> Self {
                Self::from(elem.canonical_representation())
            }
        }

        impl From<&FieldElement> for $t {
            fn from(elem: &FieldElement) -> Self {
                Self::from(elem.canonical_representation())
            }
        }
    )+};
}

impl_from_fe_for_int!(u32, u64, u128, i128);

impl From<FieldElement> for i32 {
    fn from(elem: FieldElement) -> Self {
        fe_to_i32(&elem)
    }
}

impl From<&FieldElement> for i32 {
    fn from(elem: &FieldElement) -> Self {
        fe_to_i32(elem)
    }
}

const fn fe_to_i32(fe: &FieldElement) -> i32 {
    let v = fe.canonical_representation();
    if v <= i32::MAX as u32 {
        v as i32
    } else {
        (v as i64 - FieldElement::P as i64) as i32
    }
}

impl From<FieldElement> for i64 {
    fn from(elem: FieldElement) -> Self {
        elem.canonical_representation() as i64
    }
}

impl From<&FieldElement> for i64 {
    fn from(elem: &FieldElement) -> Self {
        elem.canonical_representation() as i64
    }
}

/// Convert a B-field element to a byte array.
/// The client uses this for its database.
impl From<FieldElement> for [u8; FieldElement::BYTES] {
    fn from(fe: FieldElement) -> Self {
        // It's crucial to map this to the canonical representation before converting.
        // Otherwise, the representation is degenerate.
        fe.canonical_representation().to_le_bytes()
    }
}

impl TryFrom<[u8; FieldElement::BYTES]> for FieldElement {
    type Error = ParseFieldElementError;

    fn try_from(array: [u8; FieldElement::BYTES]) -> Result<Self, Self::Error> {
        Self::try_new(u32::from_le_bytes(array))
    }
}

impl TryFrom<&[u8]> for FieldElement {
    type Error = ParseFieldElementError;

    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        <[u8; FieldElement::BYTES]>::try_from(bytes)
            .map_err(|_| Self::Error::InvalidNumBytes(bytes.len()))?
            .try_into()
    }
}

impl Inverse for FieldElement {
    #[inline]
    fn inverse(&self) -> Self {
        let x = *self;
        assert_ne!(
            x,
            Self::zero(),
            "Attempted to find the multiplicative inverse of zero."
        );

        // Use Fermat's little theorem: a^(p-1) = 1 mod p
        // So a^(-1) = a^(p-2) mod p
        x.mod_pow(Self::P - 2)
    }
}

impl ModPowU32 for FieldElement {
    #[inline]
    fn mod_pow_u32(&self, exp: u32) -> Self {
        self.mod_pow(exp)
    }
}

impl CyclicGroupGenerator for FieldElement {
    fn get_cyclic_group_elements(&self, max: Option<usize>) -> Vec<Self> {
        let mut val = *self;
        let mut ret: Vec<Self> = vec![Self::one()];

        loop {
            ret.push(val);
            val *= *self;
            if val.is_one() || max.is_some() && ret.len() >= max.unwrap() {
                break;
            }
        }
        ret
    }
}

impl rand::prelude::Distribution<FieldElement> for rand_distr::Standard {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> FieldElement {
        FieldElement::new(rng.gen_range(0..=FieldElement::MAX))
    }
}

impl FiniteField for FieldElement {
    fn to_le_bytes(&self) -> Vec<u8> {
        self.0.to_le_bytes().to_vec()
    }

    fn centered_absolute_value(&self) -> u32 {
        let value = self.value();
        let half_p = FieldElement::P / 2;
        if value > half_p {
            FieldElement::P - value
        } else {
            value
        }
    }
}

impl Zero for FieldElement {
    #[inline]
    fn zero() -> Self {
        Self::ZERO
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self == &Self::ZERO
    }
}

impl ConstZero for FieldElement {
    const ZERO: Self = Self::new(0);
}

impl One for FieldElement {
    #[inline]
    fn one() -> Self {
        Self::ONE
    }

    #[inline]
    fn is_one(&self) -> bool {
        self == &Self::ONE
    }
}

impl ConstOne for FieldElement {
    const ONE: Self = Self::new(1);
}

impl Add for FieldElement {
    type Output = Self;

    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        let sum = self.0 as u64 + rhs.0 as u64;
        if sum >= Self::P as u64 {
            Self((sum - Self::P as u64) as u32)
        } else {
            Self(sum as u32)
        }
    }
}

impl AddAssign for FieldElement {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl SubAssign for FieldElement {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl MulAssign for FieldElement {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Mul for FieldElement {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self(Self::montyred((self.0 as u64) * (rhs.0 as u64)))
    }
}

impl Neg for FieldElement {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::zero() - self
    }
}

impl Sub for FieldElement {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        if self.0 >= rhs.0 {
            Self(self.0 - rhs.0)
        } else {
            Self(self.0 + Self::P - rhs.0)
        }
    }
}

impl Div for FieldElement {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, other: Self) -> Self {
        other.inverse() * self
    }
}

impl ModPowU64 for FieldElement {
    #[inline]
    fn mod_pow_u64(&self, pow: u64) -> Self {
        // For large exponents, reduce modulo P-1 first (Fermat's little theorem)
        let reduced_pow = (pow % (Self::P - 1) as u64) as u32;
        self.mod_pow(reduced_pow)
    }
}

impl PrimitiveRootOfUnity for FieldElement {
    fn primitive_root_of_unity(n: u32) -> Option<FieldElement> {
        PRIMITIVE_ROOTS.get(&n).map(|&r| FieldElement::new(r))
    }
}

pub mod other {
    use rand::distributions::Distribution;
    use rand::distributions::Standard;
    use rand::Rng;

    pub fn random_elements<T>(n: usize) -> Vec<T>
    where
        Standard: Distribution<T>,
    {
        rand::thread_rng().sample_iter(Standard).take(n).collect()
    }
}

#[cfg(test)]
mod b_prime_field_element_test {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::Hasher;

    use itertools::izip;
    use proptest::prelude::*;
    use proptest_arbitrary_interop::arb;
    use serde_json;
    use test_strategy::proptest;

    use super::*;

    impl proptest::arbitrary::Arbitrary for FieldElement {
        type Parameters = ();

        fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
            arb().boxed()
        }

        type Strategy = BoxedStrategy<Self>;
    }

    #[proptest]
    fn get_size(fe: FieldElement) {
        prop_assert_eq!(4, fe.get_size());
    }

    #[proptest]
    fn serialization_and_deserialization_to_and_from_json_is_identity(
        fe: FieldElement,
    ) {
        let serialized = serde_json::to_string(&fe).unwrap();
        let deserialized: FieldElement =
            serde_json::from_str(&serialized).unwrap();
        prop_assert_eq!(fe, deserialized);
    }

    #[proptest]
    fn deserializing_u32_is_like_calling_new(
        #[strategy(0..=FieldElement::MAX)] value: u32,
    ) {
        let fe = FieldElement::new(value);
        let deserialized: FieldElement =
            serde_json::from_str(&value.to_string()).unwrap();
        prop_assert_eq!(fe, deserialized);
    }

    #[proptest]
    fn parsing_string_representing_canonical_u32_gives_correct_bfield_element(
        #[strategy(0..=FieldElement::MAX)] v: u32,
    ) {
        let fe = FieldElement::from_str(&v.to_string()).unwrap();
        prop_assert_eq!(v, fe.value());
    }

    #[proptest]
    fn parsing_string_representing_too_big_u32_as_bfield_element_gives_error(
        #[strategy(FieldElement::P..)] v: u32,
    ) {
        let err = FieldElement::from_str(&v.to_string()).err().unwrap();
        prop_assert_eq!(ParseFieldElementError::NotCanonical(v as u64), err);
    }

    #[proptest]
    fn zero_is_neutral_element_for_addition(fe: FieldElement) {
        let zero = FieldElement::ZERO;
        prop_assert_eq!(fe + zero, fe);
    }

    #[proptest]
    fn one_is_neutral_element_for_multiplication(fe: FieldElement) {
        let one = FieldElement::ONE;
        prop_assert_eq!(fe * one, fe);
    }

    #[proptest]
    fn addition_is_commutative(
        element_0: FieldElement,
        element_1: FieldElement,
    ) {
        prop_assert_eq!(element_0 + element_1, element_1 + element_0);
    }

    #[proptest]
    fn multiplication_is_commutative(
        element_0: FieldElement,
        element_1: FieldElement,
    ) {
        prop_assert_eq!(element_0 * element_1, element_1 * element_0);
    }

    #[proptest]
    fn addition_is_associative(
        element_0: FieldElement,
        element_1: FieldElement,
        element_2: FieldElement,
    ) {
        prop_assert_eq!(
            (element_0 + element_1) + element_2,
            element_0 + (element_1 + element_2)
        );
    }

    #[proptest]
    fn multiplication_is_associative(
        element_0: FieldElement,
        element_1: FieldElement,
        element_2: FieldElement,
    ) {
        prop_assert_eq!(
            (element_0 * element_1) * element_2,
            element_0 * (element_1 * element_2)
        );
    }

    #[proptest]
    fn multiplication_distributes_over_addition(
        element_0: FieldElement,
        element_1: FieldElement,
        element_2: FieldElement,
    ) {
        prop_assert_eq!(
            element_0 * (element_1 + element_2),
            element_0 * element_1 + element_0 * element_2
        );
    }

    #[proptest]
    fn multiplication_with_inverse_gives_identity(
        #[filter(!#fe.is_zero())] fe: FieldElement,
    ) {
        prop_assert!((fe.inverse() * fe).is_one());
    }

    #[proptest]
    fn division_by_self_gives_identity(
        #[filter(!#fe.is_zero())] fe: FieldElement,
    ) {
        prop_assert!((fe / fe).is_one());
    }

    #[proptest]
    fn values_larger_than_modulus_are_handled_correctly(
        #[strategy(FieldElement::P..=u32::MAX)] large_value: u32,
    ) {
        let fe = FieldElement::new(large_value);
        let expected_value = large_value % FieldElement::P;
        prop_assert_eq!(expected_value, fe.value());
    }

    #[test]
    fn display_test() {
        let seven: FieldElement = FieldElement::new(7);
        assert_eq!("7", format!("{seven}"));

        let minus_one: FieldElement = FieldElement::new(FieldElement::P - 1);
        assert_eq!("-1", format!("{minus_one}"));

        let minus_fifteen: FieldElement =
            FieldElement::new(FieldElement::P - 15);
        assert_eq!("-15", format!("{minus_fifteen}"));
    }

    #[test]
    fn zero_is_zero() {
        let zero = FieldElement::zero();
        assert!(zero.is_zero());
        assert_eq!(zero, FieldElement::ZERO);
    }

    #[proptest]
    fn not_zero_is_nonzero(fe: FieldElement) {
        if fe.value() == 0 {
            return Ok(());
        }
        prop_assert!(!fe.is_zero());
    }

    #[test]
    fn one_is_one() {
        let one = FieldElement::one();
        assert!(one.is_one());
        assert_eq!(one, FieldElement::ONE);
    }

    #[proptest]
    fn not_one_is_not_one(fe: FieldElement) {
        if fe.value() == 1 {
            return Ok(());
        }
        prop_assert!(!fe.is_one());
    }

    #[test]
    fn one_unequal_zero() {
        let one = FieldElement::ONE;
        let zero = FieldElement::ZERO;
        assert_ne!(one, zero);
    }

    #[proptest]
    fn byte_array_of_small_field_elements_is_zero_at_high_indices(value: u8) {
        let fe = FieldElement::new(value as u32);
        let byte_array: [u8; 4] = fe.into();

        prop_assert_eq!(value, byte_array[0]);
        (1..4).for_each(|i| {
            assert_eq!(0, byte_array[i]);
        });
    }

    #[proptest]
    fn byte_array_conversion(fe: FieldElement) {
        let array: [u8; 4] = fe.into();
        let fe_recalculated: FieldElement = array.try_into()?;
        prop_assert_eq!(fe, fe_recalculated);
    }

    #[proptest]
    fn byte_array_outside_range_is_not_accepted(
        #[strategy(FieldElement::P..=u32::MAX)] value: u32,
    ) {
        let byte_array = value.to_le_bytes();
        prop_assert!(FieldElement::try_from(byte_array).is_err());
    }

    #[proptest]
    fn value_is_preserved(#[strategy(0..FieldElement::P)] value: u32) {
        prop_assert_eq!(value, FieldElement::new(value).value());
    }

    #[test]
    fn supposed_generator_is_generator() {
        let generator = FieldElement::generator();
        let largest_meaningful_power = FieldElement::P - 1;
        let generator_pow_p = generator.mod_pow(largest_meaningful_power);
        let generator_pow_p_half =
            generator.mod_pow(largest_meaningful_power / 2);

        assert_eq!(FieldElement::ONE, generator_pow_p);
        assert_ne!(FieldElement::ONE, generator_pow_p_half);
    }

    #[proptest]
    fn increment(mut fe: FieldElement) {
        let old_value = fe.value();
        fe.increment();
        let expected_value = (old_value + 1) % FieldElement::P;
        prop_assert_eq!(expected_value, fe.value());
    }

    #[test]
    fn incrementing_max_value_wraps_around() {
        let mut fe = FieldElement::new(FieldElement::MAX);
        fe.increment();
        assert_eq!(0, fe.value());
    }

    #[proptest]
    fn decrement(mut fe: FieldElement) {
        let old_value = fe.value();
        fe.decrement();
        let expected_value =
            old_value.checked_sub(1).unwrap_or(FieldElement::P - 1);
        prop_assert_eq!(expected_value, fe.value());
    }

    #[test]
    fn decrementing_min_value_wraps_around() {
        let mut fe = FieldElement::ZERO;
        fe.decrement();
        assert_eq!(FieldElement::MAX, fe.value());
    }

    #[test]
    fn empty_batch_inversion() {
        let empty_inv = FieldElement::batch_inversion(vec![]);
        assert!(empty_inv.is_empty());
    }

    #[proptest]
    fn batch_inversion(fes: Vec<FieldElement>) {
        // Filter out zero elements since batch_inversion panics on zero
        let non_zero_fes: Vec<_> =
            fes.into_iter().filter(|fe| !fe.is_zero()).collect();

        if non_zero_fes.is_empty() {
            return Ok(());
        }

        let fes_inv = FieldElement::batch_inversion(non_zero_fes.clone());
        prop_assert_eq!(non_zero_fes.len(), fes_inv.len());
        for (fe, fe_inv) in izip!(non_zero_fes, fes_inv) {
            prop_assert_eq!(FieldElement::ONE, fe * fe_inv);
        }
    }

    #[test]
    fn mul_div_pbt() {
        // Verify that the mul result is sane
        let rands: Vec<FieldElement> = other::random_elements(100);
        for i in 1..rands.len() {
            let prod_mul = rands[i - 1] * rands[i];
            let mut prod_mul_assign = rands[i - 1];
            prod_mul_assign *= rands[i];
            assert_eq!(
                prod_mul, prod_mul_assign,
                "mul and mul_assign must be the same for B field elements"
            );
            assert_eq!(prod_mul / rands[i - 1], rands[i]);
            assert_eq!(prod_mul / rands[i], rands[i - 1]);
        }
    }

    #[test]
    fn add_sub_wrap_around_test() {
        // Test wrap-around for our smaller prime
        let element = FieldElement::new(4);
        let sum = FieldElement::new(FieldElement::MAX) + element;
        assert_eq!(FieldElement::new(3), sum);
        let diff = sum - element;
        assert_eq!(FieldElement::new(FieldElement::MAX), diff);
    }

    #[test]
    fn neg_test() {
        assert_eq!(-FieldElement::ZERO, FieldElement::ZERO);
        assert_eq!(
            (-FieldElement::ONE).canonical_representation(),
            FieldElement::MAX
        );
        let max = FieldElement::new(FieldElement::MAX);
        let max_plus_one = max + FieldElement::ONE;
        let max_plus_two = max_plus_one + FieldElement::ONE;
        assert_eq!(FieldElement::ZERO, -max_plus_one);
        assert_eq!(max, -max_plus_two);
    }

    #[test]
    fn equality_and_hash_test() {
        assert_eq!(FieldElement::ZERO, FieldElement::ZERO);
        assert_eq!(FieldElement::ONE, FieldElement::ONE);
        assert_ne!(FieldElement::ONE, FieldElement::ZERO);
        assert_eq!(FieldElement::new(42), FieldElement::new(42));
        assert_ne!(FieldElement::new(42), FieldElement::new(43));

        assert_eq!(
            FieldElement::new(102),
            FieldElement::new(FieldElement::MAX) + FieldElement::new(103)
        );
        assert_ne!(
            FieldElement::new(103),
            FieldElement::new(FieldElement::MAX) + FieldElement::new(103)
        );

        // Verify that hashing works for canonical representations
        let mut hasher_a = DefaultHasher::new();
        let mut hasher_b = DefaultHasher::new();

        std::hash::Hash::hash(&FieldElement::new(42), &mut hasher_a);
        std::hash::Hash::hash(&FieldElement::new(42), &mut hasher_b);
        assert_eq!(hasher_a.finish(), hasher_b.finish());

        // Verify that hashing works for non-canonical representations
        hasher_a = DefaultHasher::new();
        hasher_b = DefaultHasher::new();
        let non_canonical =
            FieldElement::new(FieldElement::MAX) + FieldElement::new(103);
        std::hash::Hash::hash(&(non_canonical), &mut hasher_a);
        std::hash::Hash::hash(&FieldElement::new(102), &mut hasher_b);
        assert_eq!(hasher_a.finish(), hasher_b.finish());
    }

    #[test]
    fn mod_pow_test_powers_of_two() {
        let two = FieldElement::new(2);
        // 2^23 = 8388608 > 8380417, so we'll see wrap-around after i=23
        for i in 0..24 {
            //let expected = fe!(FieldElement::montyred(2u64.pow(i)));
            let value = 2u64.pow(i) % FieldElement::P as u64;
            let expected = FieldElement::new(value as u32);
            assert_eq!(expected, two.mod_pow(i));
        }
    }

    #[test]
    fn mod_pow_test_powers_of_three() {
        let three = FieldElement::new(3);
        // 3^15 = 14348907 > 8380417, so we'll see wrap-around after i=14
        for i in 0..16 {
            let value = 3u64.pow(i) % FieldElement::P as u64;
            let expected = FieldElement::new(value as u32);
            assert_eq!(expected, three.mod_pow(i));
        }
    }

    #[test]
    fn mod_pow_test() {
        // Test specific values for our prime
        // These are primitive roots of unity we computed
        assert!(FieldElement::new(4808194).mod_pow(4).is_one());
        assert_eq!(
            FieldElement::new(4808194),
            FieldElement::new(4808194).mod_pow(5)
        );
        assert!(FieldElement::new(8380416).mod_pow(2).is_one());
        assert!(FieldElement::new(4614810).mod_pow(8).is_one());
        assert!(FieldElement::new(949247).mod_pow(3).is_one());
        assert!(FieldElement::new(4117873).mod_pow(12).is_one());
        assert!(FieldElement::new(0).mod_pow(0).is_one());
    }

    #[test]
    fn get_primitive_root_of_unity_test() {
        let mut errors = Vec::with_capacity(14);

        for i in 0..14 {
            let power = 1 << i;
            match FieldElement::primitive_root_of_unity(power) {
                Some(root) => {
                    if !root.mod_pow(power).is_one() {
                        errors.push(format!("root^{} != 1", power));
                    }
                    if power > 1 && root.mod_pow(power / 2).is_one() {
                        errors.push(format!(
                            "root^{} == 1, not primitive",
                            power / 2
                        ));
                    }
                }
                None => errors.push(format!("No root for n = {}", power)),
            }
        }

        assert!(
            errors.is_empty(),
            "Primitive root errors:\n{}",
            errors.join("\n")
        );
    }

    #[test]
    #[should_panic(
        expected = "Attempted to find the multiplicative inverse of zero."
    )]
    fn multiplicative_inverse_of_zero() {
        let zero = FieldElement::ZERO;
        let _ = zero.inverse();
    }

    #[test]
    fn u32_conversion() {
        // Test conversion to/from u32
        let val = FieldElement::new(1234567);
        let as_u32: u32 = val.into();
        assert_eq!(1234567, as_u32);

        // Test values that exceed our prime
        for i in 0..100 {
            let val_exceeding_p = FieldElement::P + i;
            let fe = FieldElement::new(val_exceeding_p);
            let expected = val_exceeding_p % FieldElement::P;
            assert_eq!(expected, fe.value());
        }
    }

    // Todo fix it
    // #[test]
    // fn inverse_or_zero_fe() {
    //     let zero = FieldElement::ZERO;
    //     let one = FieldElement::ONE;
    //     assert_eq!(zero, zero.inverse_or_zero());

    //     let mut rng = rand::thread_rng();
    //     let elem: FieldElement = rng.gen();
    //     if elem.is_zero() {
    //         assert_eq!(zero, elem.inverse_or_zero())
    //     } else {
    //         assert_eq!(one, elem * elem.inverse_or_zero());
    //     }
    // }

    // TODO fix it
    // #[test]
    // fn test_random_squares() {
    //     let mut rng = thread_rng();
    //     let p = FieldElement::P;
    //     for _ in 0..100 {
    //         let a = rng.gen_range(0..p);
    //         let asq = (((a as u64) * (a as u64)) % (p as u64)) as u32;
    //         let b = FieldElement::new(a);
    //         let bsq = FieldElement::new(asq);
    //         assert_eq!(bsq, b * b);
    //         assert_eq!(bsq.value(), (b * b).value());
    //         assert_eq!(b.value(), a);
    //         assert_eq!(bsq.value(), asq);
    //     }
    //     let one = FieldElement::new(1);
    //     assert_eq!(one, one * one);
    // }

    #[test]
    fn equals() {
        let a = FieldElement::ONE;
        let b = fe!(FieldElement::MAX) * fe!(FieldElement::MAX);

        // elements are equal
        assert_eq!(a, b);
        assert_eq!(a.value(), b.value());
    }

    //TODO fix it
    // #[test]
    // fn test_random_raw() {

    //     let mut rng = rand::thread_rng();
    //     for _ in 0..100 {
    //         let e: FieldElement = rng.gen();
    //         let bytes = e.raw_bytes();
    //         let c = FieldElement::from_raw_bytes(&bytes);
    //         assert_eq!(e, c);

    //         let mut f = 0u32;
    //         for (i, b) in bytes.iter().enumerate() {
    //             f += (*b as u32) << (8 * i);
    //         }
    //         assert_eq!(e, FieldElement(f));

    //         let chunks = e.raw_u16s();
    //         let g = FieldElement::from_raw_u16s(&chunks);
    //         assert_eq!(e, g);

    //         let mut h = 0u32;
    //         for (i, ch) in chunks.iter().enumerate() {
    //             h += (*ch as u32) << (16 * i);
    //         }
    //         assert_eq!(e, FieldElement(h));
    //     }
    // }

    #[test]
    fn test_fixed_inverse() {
        // Test specific inverse pairs for our prime
        let a = FieldElement::new(12345);
        let a_inv = a.inverse();
        let a_inv_or_0 = a.inverse_or_zero();
        assert_eq!(a_inv, a_inv_or_0);
        assert_eq!(FieldElement::ONE, a * a_inv);

        // Test another pair
        let b = FieldElement::new(7654321);
        let b_inv = b.inverse();
        assert_eq!(FieldElement::ONE, b * b_inv);
    }

    #[test]
    fn test_fixed_modpow() {
        let exponent = 1234567u32;
        let base = FieldElement::new(42);
        let result = base.mod_pow_u32(exponent);
        // Verify it's correct by checking base^exponent * base^(-exponent) = 1
        let inv_result = base.mod_pow_u32(FieldElement::P - 1 - exponent);
        assert_eq!(FieldElement::ONE, result * inv_result);
    }

    #[test]
    fn test_fixed_mul() {
        {
            let a = FieldElement::new(123456);
            let b = FieldElement::new(789012);
            let c = a * b;
            // 123456 * 789012 = 97,408,265,472
            // 97,408,265,472 mod 8380417 = 2678681
            let expected = FieldElement::new(2678681);
            assert_eq!(c, expected);
        }

        {
            let a = FieldElement::new(4190208); // (P-1)/2
            let b = FieldElement::new(2);
            let c = a * b;
            let expected = FieldElement::new(8380416); // P-1
            assert_eq!(c, expected);
        }
    }

    #[proptest]
    fn conversion_from_i32_to_fe_is_correct(v: i32) {
        let fe = FieldElement::from(v);

        // Calculate expected value using modular arithmetic
        let expected = if v >= 0 {
            (v as u32) % FieldElement::P
        } else {
            let remainder = v % (FieldElement::P as i32);
            if remainder == 0 {
                0
            } else {
                // remainder is negative, so P + remainder gives positive representative
                (FieldElement::P as i32 + remainder) as u32
            }
        };

        prop_assert_eq!(expected, fe.value());
    }

    #[proptest]
    fn conversion_from_isize_to_fe_is_correct(v: isize) {
        let fe = FieldElement::from(v);
        match v {
            0.. => prop_assert_eq!(
                (v as u64 % FieldElement::P as u64) as u32,
                fe.value()
            ),
            _ => {
                let expected = if (-v as u64) < FieldElement::P as u64 {
                    FieldElement::P - (-v as u32)
                } else {
                    FieldElement::P
                        - ((-v as u64) % FieldElement::P as u64) as u32
                };
                prop_assert_eq!(expected, fe.value())
            }
        }
    }

    #[test]
    fn bfield_element_can_be_converted_to_and_from_many_types() {
        let _ = fe!(0_u8);
        let _ = fe!(0_u16);
        let _ = fe!(0_u32);
        let _ = fe!(0_u64);
        let _ = fe!(0_u128);
        let _ = fe!(0_usize);

        let max = fe!(FieldElement::MAX);
        assert_eq!(max, fe!(-1_i8));
        assert_eq!(max, fe!(-1_i16));
        assert_eq!(max, fe!(-1_i32));
        assert_eq!(max, fe!(-1_i64));
        assert_eq!(max, fe!(-1_isize));

        assert!(u8::try_from(FieldElement::ZERO).is_ok());
        assert!(i8::try_from(FieldElement::ZERO).is_ok());
        assert!(u16::try_from(FieldElement::ZERO).is_ok());
        assert!(i16::try_from(FieldElement::ZERO).is_ok());
        assert!(usize::try_from(FieldElement::ZERO).is_ok());
        assert!(isize::try_from(FieldElement::ZERO).is_ok());

        let _ = u32::from(max);
        let _ = i32::from(max);
        let _ = u64::from(max);
        let _ = i64::from(max);
        let _ = u128::from(max);
        let _ = i128::from(max);
    }

    #[test]
    fn bfield_conversion_works_for_types_min_and_max() {
        let _ = fe!(u8::MIN);
        let _ = fe!(u8::MAX);
        let _ = fe!(u16::MIN);
        let _ = fe!(u16::MAX);
        let _ = fe!(u32::MIN);
        let _ = fe!(u32::MAX);
        let _ = fe!(u64::MIN);
        let _ = fe!(u64::MAX);
        let _ = fe!(u128::MIN);
        let _ = fe!(u128::MAX);
        let _ = fe!(usize::MIN);
        let _ = fe!(usize::MAX);
        let _ = fe!(i8::MIN);
        let _ = fe!(i8::MAX);
        let _ = fe!(i16::MIN);
        let _ = fe!(i16::MAX);
        let _ = fe!(i32::MIN);
        let _ = fe!(i32::MAX);
        let _ = fe!(i64::MIN);
        let _ = fe!(i64::MAX);
        let _ = fe!(isize::MIN);
        let _ = fe!(isize::MAX);
    }

    #[proptest]
    fn naive_and_actual_conversion_from_i64_agree(v: i64) {
        fn naive_conversion(x: i64) -> FieldElement {
            let p = FieldElement::P as i128;
            let value = i128::from(x).rem_euclid(p) as u32;
            FieldElement::new(value)
        }

        prop_assert_eq!(naive_conversion(v), FieldElement::from(v));
    }

    #[test]
    fn fe_macro_can_be_used() {
        let b = fe!(42);
        let _ = fe!(42u32);
        let _ = fe!(-1);
        let _ = fe!(b);
        let _ = fe!(b.0);
        let _ = fe!(42_usize);
        let _ = fe!(-2_isize);

        let c: Vec<FieldElement> = fe_vec![1, 2, 3];
        let d: [FieldElement; 3] = fe_array![1, 2, 3];
        assert_eq!(c, d);
    }

    #[proptest]
    fn fe_macro_produces_same_result_as_calling_new(value: u32) {
        prop_assert_eq!(FieldElement::new(value), fe!(value));
    }

    #[test]
    fn const_minus_two_inverse_is_really_minus_two_inverse() {
        assert_eq!(fe!(-2).inverse(), FieldElement::MINUS_TWO_INVERSE);
    }

    #[test]
    fn test_specific_primitive_roots() {
        // Test some of our computed primitive roots
        let root_2 = FieldElement::primitive_root_of_unity(2).unwrap();
        assert_eq!(FieldElement::new(8380416), root_2);
        assert_eq!(FieldElement::ONE, root_2.mod_pow(2));
        assert_eq!(FieldElement::new(8380416), root_2.mod_pow(1));

        let root_4 = FieldElement::primitive_root_of_unity(4).unwrap();
        assert_eq!(FieldElement::new(4808194), root_4);
        assert_eq!(FieldElement::ONE, root_4.mod_pow(4));
        assert_eq!(FieldElement::new(8380416), root_4.mod_pow(2));

        let root_8 = FieldElement::primitive_root_of_unity(8).unwrap();
        assert_eq!(FieldElement::new(4614810), root_8);
        assert_eq!(FieldElement::ONE, root_8.mod_pow(8));
        assert_ne!(FieldElement::ONE, root_8.mod_pow(4));
    }

    #[test]
    fn test_edge_cases_for_small_prime() {
        // Test values near the prime boundary
        let near_p = FieldElement::new(FieldElement::P - 1);
        let one = FieldElement::ONE;
        assert_eq!(FieldElement::ZERO, near_p + one);

        // Test Montgomery form consistency
        let a = FieldElement::new(1000000);
        let b = FieldElement::new(2000000);
        let c = a + b;
        assert_eq!(FieldElement::new(3000000), c);

        // Test that values >= P are properly reduced
        let over_p = FieldElement::new(FieldElement::P + 100);
        assert_eq!(FieldElement::new(100), over_p);
    }
}

#[cfg(test)]
mod montyred_tests {
    use super::*;
    use proptest::prelude::*;
    use test_strategy::proptest;

    #[test]
    fn test_montyred_zero() {
        assert_eq!(0, FieldElement::montyred(0));
    }

    #[test]
    fn test_montyred_identity() {
        // R = 2^32 mod P = 4193792
        // montyred(R) should give 1
        const R: u64 = 4193792;
        assert_eq!(1, FieldElement::montyred(R));
    }

    #[test]
    fn test_montyred_specific_values() {
        // Test with known values
        // These can be verified with an independent implementation

        // montyred(P) = 0
        assert_eq!(0, FieldElement::montyred(FieldElement::P as u64));

        // montyred(2*P) = 0
        assert_eq!(0, FieldElement::montyred(2 * FieldElement::P as u64));

        // montyred(R^2 mod P) = R mod P = 4193792
        // Since montyred(x) = x * R^(-1) mod P
        const R_MOD_P: u32 = 4193792; // R mod P where R = 2^32
        assert_eq!(R_MOD_P, FieldElement::montyred(FieldElement::R2 as u64));
    }

    #[test]
    fn test_montyred_max_input() {
        // Test with maximum possible input (product of two max Montgomery values)
        let max_mont = FieldElement::P - 1;
        let product = (max_mont as u64) * (max_mont as u64);
        let result = FieldElement::montyred(product);

        // Result should be in range [0, P)
        assert!(result < FieldElement::P);
    }

    #[proptest]
    fn montyred_output_always_in_range(a: u64) {
        let result = FieldElement::montyred(a);
        prop_assert!(result < FieldElement::P);
    }

    #[proptest]
    fn montyred_is_deterministic(a: u64) {
        let result1 = FieldElement::montyred(a);
        let result2 = FieldElement::montyred(a);
        prop_assert_eq!(result1, result2);
    }

    #[test]
    fn test_montyred_consistency_with_operations() {
        // Verify that a * R^-1 mod P is computed correctly
        let a = 12345u64;
        let a_mont = FieldElement::new(a as u32);
        let a_back = a_mont.canonical_representation();
        assert_eq!(a as u32, a_back);
    }

    #[proptest]
    fn montgomery_multiplication_property(a: u32, b: u32) {
        let a = a % FieldElement::P;
        let b = b % FieldElement::P;

        let a_mont = FieldElement::new(a);
        let b_mont = FieldElement::new(b);
        let c_mont = a_mont * b_mont;
        let c = c_mont.value();

        // Verify that c = (a * b) mod P
        let expected = ((a as u64 * b as u64) % FieldElement::P as u64) as u32;
        prop_assert_eq!(expected, c);
    }

    #[test]
    fn test_montyred_intermediate_overflow() {
        // Test values that might cause overflow in intermediate calculations
        let large_val = u64::MAX / 2;
        let result = FieldElement::montyred(large_val);
        assert!(result < FieldElement::P);
    }

    #[test]
    fn verify_n_prime_constant() {
        // Verify that N_PRIME is correct
        // N' * P ≡ -1 (mod 2^32)
        const N_PRIME: u32 = 4236238847;
        let product = (N_PRIME as u64).wrapping_mul(FieldElement::P as u64);
        let low_32 = product as u32;

        // Should equal 2^32 - 1 when taken mod 2^32
        assert_eq!(u32::MAX, low_32);
    }
}
