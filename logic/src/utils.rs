//! Small helpers to keep code DRY and straightforward.

use num_traits::Zero;
use rand::RngCore;


use math::{
    poly::Polynomial,
    traits::FiniteField,
};



/// Build an array by mapping a function over indices.
#[inline]
pub fn arr_from_fn<T, const M: usize>(mut f: impl FnMut(usize) -> T) -> [T; M] {
    std::array::from_fn(|i| f(i))
}

/// Map over a fixed-size array by reference, returning a new array.
#[inline]
pub fn arr_map_ref<T, U, const M: usize>(xs: &[T; M], mut f: impl FnMut(&T) -> U) -> [U; M] {
    std::array::from_fn(|i| f(&xs[i]))
}

/// Zip-map over two fixed-size arrays by reference, returning a new array.
#[inline]
pub fn arr_zip_map_ref<A, B, U, const M: usize>(
    a: &[A; M],
    b: &[B; M],
    mut f: impl FnMut(&A, &B) -> U,
) -> [U; M] {
    std::array::from_fn(|i| f(&a[i], &b[i]))
}

/// A common pattern in this codebase: create a `[Polynomial; M]` of zeros.
#[inline]
pub fn zero_polys<FF: FiniteField, const M: usize>() -> [Polynomial<'static, FF>; M] {
    arr_from_fn(|_| Polynomial::zero())
}

/// Expand (seed || idx) using SHAKE256 to `len` bytes.
#[inline]
pub fn shake256_expand_idx(seed: &[u8], idx: u16, len: usize) -> Vec<u8> {
    let mut inp = Vec::with_capacity(seed.len() + 2);
    inp.extend_from_slice(seed);
    inp.extend_from_slice(&idx.to_le_bytes());
    crate::hash::shake256(len, &inp)
}

/// Expand (seed || idx || ctr) using SHAKE256 to `len` bytes.
#[inline]
pub fn shake256_expand_idx_ctr(seed: &[u8], idx: u16, ctr: u32, len: usize) -> Vec<u8> {
    let mut inp = Vec::with_capacity(seed.len() + 2 + 4);
    inp.extend_from_slice(seed);
    inp.extend_from_slice(&idx.to_le_bytes());
    inp.extend_from_slice(&ctr.to_le_bytes());
    crate::hash::shake256(len, &inp)
}


pub fn random_bytes()-> [u8;32] {
    let mut rng = rand::thread_rng();
    let mut tmp = [0u8; 32];
    rng.fill_bytes(&mut tmp);
    tmp
}
