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


pub trait FiniteField:
    Copy
    + Debug
    + Display
    + Default
    + Eq
    + Serialize
    + DeserializeOwned
    + Hash
    + ConstZero
    + ConstOne
    + Add<Output = Self>
    + Mul<Output = Self>
    + Sub<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + MulAssign
    + MulAssign<FieldElement>
    + SubAssign
    + CyclicGroupGenerator
    + PrimitiveRootOfUnity
    + Inverse
    + ModPowU32
    + From<u64> 
    + From<i32>
    + From<u32>
    + Send
    + Sync
{
    /// Montgomery Batch Inversion
    // Adapted from https://paulmillr.com/posts/noble-secp256k1-fast-ecc/#batch-inversion
    fn batch_inversion(input: Vec<Self>) -> Vec<Self> {
        let input_length = input.len();
        if input_length == 0 {
            return Vec::<Self>::new();
        }

        let zero = Self::zero();
        let one = Self::one();
        let mut scratch: Vec<Self> = vec![zero; input_length];
        let mut acc = one;
        scratch[0] = input[0];

        for i in 0..input_length {
            assert!(!input[i].is_zero(), "Cannot do batch inversion on zero");
            scratch[i] = acc;
            acc *= input[i];
        }

        acc = acc.inverse();

        let mut res = input;
        for i in (0..input_length).rev() {
            let tmp = acc * res[i];
            res[i] = acc * scratch[i];
            acc = tmp;
        }

        res
    }

    #[inline(always)]
    fn square(self) -> Self {
        self * self
    }
    fn to_le_bytes(&self) -> Vec<u8>;
    fn centered_absolute_value(&self) -> u32;
}
