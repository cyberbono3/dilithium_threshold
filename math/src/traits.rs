use std::fmt::Debug;
use std::fmt::Display;
use std::hash::Hash;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Neg;
use std::ops::Sub;
use std::ops::SubAssign;

use num_traits::ConstOne;
use num_traits::ConstZero;
use num_traits::One;
use num_traits::Zero;
use serde::de::DeserializeOwned;
use serde::Serialize;

use crate::field_element::FieldElement;

pub trait CyclicGroupGenerator
where
    Self: Sized,
{
    fn get_cyclic_group_elements(&self, max: Option<usize>) -> Vec<Self>;
}

pub trait PrimitiveRootOfUnity
where
    Self: Sized,
{
    // Changed from u64 to u32 since our prime only supports roots up to 2^13
    fn primitive_root_of_unity(n: u32) -> Option<Self>;
}

pub trait ModPowU64 {
    #[must_use]
    fn mod_pow_u64(&self, pow: u64) -> Self;
}

pub trait ModPowU32 {
    #[must_use]
    fn mod_pow_u32(&self, exp: u32) -> Self;
}

pub trait Inverse
where
    Self: Sized + Zero,
{
    fn inverse(&self) -> Self;

    fn inverse_or_zero(&self) -> Self {
        if self.is_zero() {
            Self::zero()
        } else {
            self.inverse()
        }
    }
}

pub trait AddGroup:
    Copy
    + ConstZero
    + Zero
    + Add<Output = Self>
    + Sub<Output = Self>
    + AddAssign
    + SubAssign
{
}

impl<T> AddGroup for T where
    T: Copy
        + ConstZero
        + Zero
        + Add<Output = T>
        + Sub<Output = T>
        + AddAssign
        + SubAssign
{
}

pub trait MulGroup:
    Copy + ConstOne + One + Mul<Output = Self> + Div<Output = Self> + MulAssign
{
}

impl<T> MulGroup for T where
    T: Copy + ConstOne + One + Mul<Output = T> + Div<Output = T> + MulAssign
{
}

pub trait CanonicalEncoding {
    fn to_le_bytes(&self) -> Vec<u8>;
    fn centered_absolute_value(&self) -> u32;
}

pub trait FieldCore:
    Debug + Display + Default + Eq + Serialize + DeserializeOwned + Hash
{
}

impl<T> FieldCore for T where
    T: Debug + Display + Default + Eq + Serialize + DeserializeOwned + Hash
{
}

pub trait FieldConversions:
    From<u64> + From<i32> + From<u32> + MulAssign<FieldElement>
{
}

impl<T> FieldConversions for T where
    T: From<u64> + From<i32> + From<u32> + MulAssign<FieldElement>
{
}

pub trait FiniteField:
    AddGroup
    + MulGroup
    + FieldCore
    + CanonicalEncoding
    + Neg<Output = Self>
    + Inverse
    + ModPowU32
    + PrimitiveRootOfUnity
    + CyclicGroupGenerator
    + FieldConversions
    + Send
    + Sync
{
    /// Attempt to compute the Montgomery batch inversion of the provided elements.
    ///
    /// Returns `None` when any of the inputs is zero. Adapted from
    /// <https://paulmillr.com/posts/noble-secp256k1-fast-ecc/#batch-inversion>.
    fn try_batch_inversion(mut input: Vec<Self>) -> Option<Vec<Self>> {
        let input_length = input.len();
        if input_length == 0 {
            return Some(Vec::new());
        }

        let one = Self::one();
        let mut scratch: Vec<Self> = Vec::with_capacity(input_length);
        let mut acc = one;

        for value in &input {
            if value.is_zero() {
                return None;
            }
            scratch.push(acc);
            acc *= *value;
        }

        acc = acc.inverse();

        for (value, prefix) in
            input.iter_mut().rev().zip(scratch.into_iter().rev())
        {
            let current = *value;
            *value = acc * prefix;
            acc *= current;
        }

        Some(input)
    }

    /// Montgomery batch inversion, panicking if any input is zero.
    fn batch_inversion(input: Vec<Self>) -> Vec<Self> {
        Self::try_batch_inversion(input)
            .expect("batch_inversion: cannot invert zero element")
    }

    #[inline(always)]
    fn square(self) -> Self {
        self * self
    }
}
