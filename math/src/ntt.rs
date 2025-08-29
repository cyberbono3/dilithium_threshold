use std::ops::MulAssign;

use num_traits::ConstOne;
use num_traits::ConstZero;

use crate::{
    error::{NttError, Result},
    field_element::FieldElement,
    traits::{FiniteField, Inverse, ModPowU32, PrimitiveRootOfUnity},
};

/// Direction for a number-theoretic transform.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum Transform {
    Forward,
    Inverse,
}

/// Fallible NTT (no panics).
pub fn try_ntt<FF>(x: &mut [FF]) -> Result<()>
where
    FF: FiniteField + MulAssign<FieldElement>,
{
    ntt_in_place(x, Transform::Forward)
}

/// Fallible INTT (no panics).
pub fn try_intt<FF>(x: &mut [FF]) -> Result<()>
where
    FF: FiniteField + MulAssign<FieldElement>,
{
    ntt_in_place(x, Transform::Inverse)
}

/// Backwards-compatible wrappers that panic on invalid input.
pub fn ntt<FF>(x: &mut [FF])
where
    FF: FiniteField + MulAssign<FieldElement>,
{
    try_ntt(x).expect("ntt: slice length must be a power of two <= u32::MAX and have a root of unity");
}

pub fn intt<FF>(x: &mut [FF])
where
    FF: FiniteField + MulAssign<FieldElement>,
{
    try_intt(x).expect("intt: slice length must be a power of two <= u32::MAX and have a root of unity");
}

/// Like [NTT][self::ntt], but with greater control over the root of unity that
/// is to be used.
///
/// Does _not_ check whether
/// - the passed-in root of unity is indeed a primitive root of unity of the
///   appropriate order, or whether
/// - the passed-in log₂ of the slice length matches.
///
/// Use [NTT][self:ntt] if you want a nicer interface.
#[expect(clippy::many_single_char_names)]
#[inline]
fn ntt_unchecked<FF>(x: &mut [FF], omega: FieldElement, log2_slice_len: u32)
where
    FF: FiniteField + MulAssign<FieldElement>,
{
    let slice_len = x.len() as u32;

    for k in 0..slice_len {
        let rk = bitreverse_u32(k, log2_slice_len);
        if k < rk {
            x.swap(rk as usize, k as usize);
        }
    }

    let mut m = 1;
    for _ in 0..log2_slice_len {
        let w_m = omega.mod_pow_u32(slice_len / (2 * m));
        let mut k = 0;
        while k < slice_len {
            let mut w = FieldElement::ONE;
            for j in 0..m {
                let u = x[(k + j) as usize];
                let mut v = x[(k + j + m) as usize];
                v *= w;
                x[(k + j) as usize] = u + v;
                x[(k + j + m) as usize] = u - v;
                w *= w_m;
            }

            k += 2 * m;
        }

        m *= 2;
    }
}

#[inline]
pub fn bitreverse_usize(mut n: usize, l: usize) -> usize {
    let mut r = 0;
    for _ in 0..l {
        r = (r << 1) | (n & 1);
        n >>= 1;
    }
    r
}

pub fn bitreverse_order<FF>(array: &mut [FF]) {
    let mut logn = 0;
    while (1 << logn) < array.len() {
        logn += 1;
    }

    for k in 0..array.len() {
        let rk = bitreverse_usize(k, logn);
        if k < rk {
            array.swap(rk, k);
        }
    }
}

/// Compute the [NTT][self::ntt], but leave the array in
/// [bitreversed order][self::bitreverse_order].
///
/// This method can be expected to outperform regular NTT when
///  - it is followed up by [INTT][self::intt] (e.g. for fast multiplication)
///  - the `powers_of_omega_bitreversed` can be precomputed (which is not the
///    case here).
///
/// In that case, be sure to use the matching [`intt_noswap`] and don't forget
/// to unscale by `n`, e.g. using [`unscale`].
pub fn ntt_noswap<FF>(x: &mut [FF])
where
    FF: FiniteField + MulAssign<FieldElement>,
{
    let n: usize = x.len();
    debug_assert!(n.is_power_of_two());

    let root_order = n.try_into().unwrap();
    let omega = FieldElement::primitive_root_of_unity(root_order).unwrap();

    let mut logn: usize = 0;
    while (1 << logn) < x.len() {
        logn += 1;
    }

    let mut powers_of_omega_bitreversed = vec![FieldElement::ZERO; n];
    let mut omegai = FieldElement::ONE;
    for i in 0..n / 2 {
        powers_of_omega_bitreversed[bitreverse_usize(i, logn - 1)] = omegai;
        omegai *= omega;
    }

    let mut m: usize = 1;
    let mut t: usize = n;
    while m < n {
        t >>= 1;

        for (i, zeta) in powers_of_omega_bitreversed.iter().enumerate().take(m)
        {
            let s = i * t * 2;
            for j in s..(s + t) {
                let u = x[j];
                let mut v = x[j + t];
                v *= *zeta;
                x[j] = u + v;
                x[j + t] = u - v;
            }
        }

        m *= 2;
    }
}

/// Compute the [inverse NTT][self::intt], assuming that the array is presented
/// in [bitreversed order][self::bitreverse_order]. Also, don't unscale by `n`
/// afterward.
///
/// See also [`ntt_noswap`].
pub fn intt_noswap<FF>(x: &mut [FF])
where
    FF: FiniteField + MulAssign<FieldElement>,
{
    let n = x.len();
    debug_assert!(n.is_power_of_two());

    let root_order = n.try_into().unwrap();
    let omega = FieldElement::primitive_root_of_unity(root_order).unwrap();
    let omega_inverse = omega.inverse();

    let mut logn: usize = 0;
    while (1 << logn) < x.len() {
        logn += 1;
    }

    let mut m = 1;
    for _ in 0..logn {
        let w_m = omega_inverse.mod_pow_u32((n / (2 * m)).try_into().unwrap());
        let mut k = 0;
        while k < n {
            let mut w = FieldElement::ONE;
            for j in 0..m {
                let u = x[k + j];
                let mut v = x[k + j + m];
                v *= w;
                x[k + j] = u + v;
                x[k + j + m] = u - v;
                w *= w_m;
            }

            k += 2 * m;
        }

        m *= 2;
    }
}

/// Unscale the array by multiplying every element by the
/// inverse of the array's length. Useful for following up intt.
pub fn unscale(array: &mut [FieldElement]) {
    let ninv = FieldElement::new(array.len() as u32).inverse();
    for a in array.iter_mut() {
        *a *= ninv;
    }
}

#[inline]
fn bitreverse_u32(mut n: u32, l: u32) -> u32 {
    let mut r = 0;
    for _ in 0..l {
        r = (r << 1) | (n & 1);
        n >>= 1;
    }
    r
}

#[inline]
fn ilog2_pow2_u32(n: u32) -> u32 {
    debug_assert!(n.is_power_of_two());
    n.ilog2()
}

fn ntt_in_place<FF>(x: &mut [FF], direction: Transform) -> Result<()>
where
    FF: FiniteField + MulAssign<FieldElement>,
{
    let n_usize = x.len();
    let n_u32 =
        u32::try_from(n_usize).map_err(|_| NttError::TooLarge(n_usize))?;
    if n_u32 != 0 && !n_u32.is_power_of_two() {
        return Err(NttError::NonPowerOfTwo(n_usize).into());
    }
    if n_u32 == 0 {
        return Ok(()); // nothing to do
    }

    let root_order = n_u32;
    let omega = FieldElement::primitive_root_of_unity(root_order)
        .ok_or(NttError::MissingPrimitiveRoot(root_order))?;
    let w = match direction {
        Transform::Forward => omega,
        Transform::Inverse => omega.inverse(),
    };
    let log2 = ilog2_pow2_u32(n_u32);

    // bit-reverse permutation
    for k in 0..n_u32 {
        let rk = bitreverse_u32(k, log2);
        if k < rk {
            x.swap(rk as usize, k as usize);
        }
    }

    // Cooley–Tukey butterflies
    let mut m = 1u32;
    for _stage in 0..log2 {
        let w_m = w.mod_pow_u32(n_u32 / (2 * m));
        let mut k = 0u32;
        while k < n_u32 {
            let mut twiddle = FieldElement::ONE;
            for j in 0..m {
                let u = x[(k + j) as usize];
                let mut v = x[(k + j + m) as usize];
                v *= twiddle;
                x[(k + j) as usize] = u + v;
                x[(k + j + m) as usize] = u - v;
                twiddle *= w_m;
            }
            k += 2 * m;
        }
        m *= 2;
    }

    if matches!(direction, Transform::Inverse) {
        unscale_ffi(x, n_usize);
    }
    Ok(())
}

#[inline]
fn unscale_ffi<FF>(x: &mut [FF], n: usize)
where
    FF: FiniteField + MulAssign<FieldElement>,
{
    let ninv = FieldElement::new(n as u32).inverse();
    for a in x.iter_mut() {
        *a *= ninv;
    }
}

#[cfg(test)]
mod fast_ntt_attempt_tests {
    use proptest::collection::vec;
    use proptest::prelude::*;
    use proptest_arbitrary_interop::arb;
    use test_strategy::proptest;

    use super::*;
    use crate::field_element::other::random_elements;
    use crate::prelude::*;
    use crate::traits::PrimitiveRootOfUnity;

    #[test]
    fn chu_ntt_b_field_prop_test() {
        for log_2_n in 1..10 {
            let n = 1 << log_2_n;
            for _ in 0..10 {
                let mut values = random_elements(n);
                let original_values = values.clone();
                ntt::<FieldElement>(&mut values);
                assert_ne!(original_values, values);
                intt::<FieldElement>(&mut values);
                assert_eq!(original_values, values);

                values[0] = fe!(FieldElement::MAX);
                let original_values_with_max_element = values.clone();
                ntt::<FieldElement>(&mut values);
                assert_ne!(original_values, values);
                intt::<FieldElement>(&mut values);
                assert_eq!(original_values_with_max_element, values);
            }
        }
    }

    #[test]
    fn field_basic_test_of_chu_ntt() {
        let mut input_output = fe_vec!(1, 4, 0, 0);

        let original_input = input_output.clone();

        // For the field with prime 8380417, we need to calculate the expected values
        // The NTT of [1, 4, 0, 0] with n=4 uses omega = primitive_root_of_unity(4) = 4808194
        //
        // NTT formula: X[k] = sum(x[j] * omega^(j*k)) for j=0 to n-1
        //
        // X[0] = 1*1 + 4*1 + 0*1 + 0*1 = 5
        // X[1] = 1*1 + 4*omega + 0*omega^2 + 0*omega^3 = 1 + 4*4808194 = 1 + 19232776 mod 8380417 = 2472943
        // X[2] = 1*1 + 4*omega^2 + 0*omega^4 + 0*omega^6 = 1 + 4*8380416 = 1 + 33521664 mod 8380417 = 8380413
        // X[3] = 1*1 + 4*omega^3 + 0*omega^6 + 0*omega^9 = 1 + 4*3572223 = 1 + 14288892 mod 8380417 = 5908476

        let expected = fe_vec!(5, 2471943, 8380414, 5908476);

        ntt::<FieldElement>(&mut input_output);
        assert_eq!(expected, input_output);

        // Verify that INTT(NTT(x)) = x
        intt::<FieldElement>(&mut input_output);
        assert_eq!(original_input, input_output);
    }

    #[test]
    fn bfield_max_value_test_of_chu_ntt() {
        let mut input_output = fe_vec!(FieldElement::MAX, 0, 0, 0);
        let original_input = input_output.clone();
        let expected = fe_vec!(
            FieldElement::MAX,
            FieldElement::MAX,
            FieldElement::MAX,
            FieldElement::MAX
        );

        ntt::<FieldElement>(&mut input_output);
        assert_eq!(expected, input_output);

        // Verify that INTT(NTT(x)) = x
        intt::<FieldElement>(&mut input_output);
        assert_eq!(original_input, input_output);
    }

    #[test]
    fn ntt_on_empty_input() {
        let mut input_output = vec![];
        let original_input = input_output.clone();

        ntt::<FieldElement>(&mut input_output);
        assert_eq!(0, input_output.len());

        // Verify that INTT(NTT(x)) = x
        intt::<FieldElement>(&mut input_output);
        assert_eq!(original_input, input_output);
    }

    #[proptest]
    fn ntt_on_input_of_length_one(fe: FieldElement) {
        let mut test_vector = vec![fe];
        ntt(&mut test_vector);
        assert_eq!(vec![fe], test_vector);
    }

    #[proptest(cases = 10)]
    fn ntt_then_intt_is_identity_operation(
        #[strategy((0_usize..=13).prop_map(|l| 1 << l))] _vector_length: usize,
        #[strategy(vec(arb(), #_vector_length))] mut input: Vec<FieldElement>,
    ) {
        let original_input = input.clone();
        ntt::<FieldElement>(&mut input);
        intt::<FieldElement>(&mut input);
        assert_eq!(original_input, input);
    }

    #[test]
    fn b_field_ntt_with_length_32() {
        let mut input_output = fe_vec![
            1, 4, 0, 0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0,
            0, 0, 1, 4, 0, 0, 0, 0, 0, 0,
        ];
        let original_input = input_output.clone();

        // Transform under test
        ntt::<FieldElement>(&mut input_output);
        println!("actual_output = {input_output:?}");

        // Reference result = naive DFT at powers of the SAME ω used by NTT
        let n: u32 = original_input.len() as u32; // 32
        let omega = FieldElement::primitive_root_of_unity(n).unwrap();

        let expected: Vec<FieldElement> = (0..n as usize)
            .map(|k| {
                let mut acc = FieldElement::ZERO;
                for (j, &xj) in original_input.iter().enumerate() {
                    let pow = ((k as u32) * (j as u32)) % n;
                    acc += xj * omega.mod_pow_u32(pow);
                }
                acc
            })
            .collect();

        assert_eq!(expected, input_output, "NTT output mismatch vs naive DFT");

        // Verify that INTT(NTT(x)) = x
        let mut back = input_output.clone();
        intt::<FieldElement>(&mut back);
        assert_eq!(original_input, back, "INTT(NTT(x)) != x");
    }

    #[test]
    fn test_compare_ntt_to_eval() {
        // cover sizes 2..512
        for log_size in 1..10 {
            let n = 1usize << log_size;

            // Random polynomial in coefficient form
            let coeffs =
                crate::field_element::other::random_elements::<FieldElement>(n);
            let mut got = coeffs.clone();
            ntt(&mut got);

            // ----- Recover S_k and z_k from columns for e0 and e1 -----

            // Column for e0 (x = [1,0,0,...]): y_k = S_k * T_0
            let mut e0 = vec![FieldElement::ZERO; n];
            e0[0] = FieldElement::ONE;
            ntt(&mut e0);
            let s: Vec<FieldElement> = e0; // we'll set T_0 := 1, so s_k := S_k

            // Column for e1 (x = [0,1,0,...]): y_k = S_k * T_1 * z_k
            // choose T_1 := 1 (any nonzero choice works consistently);
            // then z_k = col1_k / s_k
            let mut e1 = vec![FieldElement::ZERO; n];
            e1[1] = FieldElement::ONE;
            ntt(&mut e1);
            let z: Vec<FieldElement> =
                (0..n).map(|k| e1[k] * s[k].inverse()).collect();

            // ----- Recover T_j from the first row (k = 0) -----
            // col_j[0] = S_0 * T_j * z_0^j  =>  T_j = col_j[0] / (S_0 * z_0^j)
            let s0 = s[0];
            let z0 = z[0];
            let mut t: Vec<FieldElement> = vec![FieldElement::ZERO; n];
            t[0] = FieldElement::ONE; // by our convention T_0 := 1

            for j in 1..n {
                let mut ej = vec![FieldElement::ZERO; n];
                ej[j] = FieldElement::ONE;
                ntt(&mut ej);
                let denom = s0 * z0.mod_pow_u32(j as u32);
                t[j] = ej[0] * denom.inverse();
            }

            // ----- Evaluate with recovered (S, z, T) and compare -----
            let mut expected = vec![FieldElement::ZERO; n];
            for k in 0..n {
                let zk = z[k];
                let mut zk_pow = FieldElement::ONE;
                let mut acc = FieldElement::ZERO;
                for j in 0..n {
                    acc += coeffs[j] * t[j] * zk_pow;
                    zk_pow *= zk;
                }
                expected[k] = s[k] * acc;
            }

            assert_eq!(
                expected, got,
                "NTT must equal scaled evaluation at recovered points"
            );
        }
    }

    #[test]
    fn test_ntt_noswap() {
        for log_size in 1..8 {
            let size = 1 << log_size;
            println!("size: {size}");
            let a: Vec<FieldElement> = random_elements(size);
            let mut a1 = a.clone();
            ntt(&mut a1);
            let mut a2 = a.clone();
            ntt_noswap(&mut a2);
            bitreverse_order(&mut a2);
            assert_eq!(a1, a2);

            intt(&mut a1);
            bitreverse_order(&mut a2);
            intt_noswap(&mut a2);
            for a2e in a2.iter_mut() {
                *a2e *= FieldElement::new(size.try_into().unwrap()).inverse();
            }
            assert_eq!(a1, a2);
        }
    }
}

#[cfg(test)]
mod ntt_in_place_tests {
    use super::*; // brings private items (like ntt_in_place and Transform) into scope
    use crate::field_element::other::random_elements;
    use crate::field_element::FieldElement;
    use crate::prelude::*;
    use crate::traits::ModPowU32;

    // --- Edge cases ---

    #[test]
    fn length_zero_is_ok_and_noop_both_directions() {
        let mut v: Vec<FieldElement> = vec![];
        let before = v.clone();
        ntt_in_place(&mut v, Transform::Forward).unwrap();
        assert_eq!(before, v);
        ntt_in_place(&mut v, Transform::Inverse).unwrap();
        assert_eq!(before, v);
    }

    #[test]
    fn length_one_is_identity_both_directions() {
        let x = FieldElement::new(1234567);
        let mut v = vec![x];
        ntt_in_place(&mut v, Transform::Forward).unwrap();
        assert_eq!(vec![x], v);
        ntt_in_place(&mut v, Transform::Inverse).unwrap();
        assert_eq!(vec![x], v);
    }

    // --- Correctness on tiny sizes ---

    #[test]
    fn length_two_matches_simple_butterfly() {
        let a = FieldElement::new(5);
        let b = FieldElement::new(9);
        let mut v = vec![a, b];

        // NTT_2([a, b]) = [a+b, a-b]
        ntt_in_place(&mut v, Transform::Forward).unwrap();
        assert_eq!(vec![a + b, a - b], v);

        // INTT_2 should bring us back (internal scaling by 1/2 is handled)
        ntt_in_place(&mut v, Transform::Inverse).unwrap();
        assert_eq!(vec![a, b], v);
    }

    #[test]
    fn four_point_known_example_forward_and_inverse() {
        // Example for the commonly used field with p = 8380417.
        // primitive_root_of_unity(4) = 4808194 and the NTT of [1,4,0,0] is:
        // [5, 2471943, 8380414, 5908476]
        let mut v = fe_vec!(1, 4, 0, 0);
        let expected = fe_vec!(5, 2471943, 8380414, 5908476);

        ntt_in_place(&mut v, Transform::Forward).unwrap();
        assert_eq!(expected, v);

        ntt_in_place(&mut v, Transform::Inverse).unwrap();
        assert_eq!(fe_vec!(1, 4, 0, 0), v);
    }

    #[test]
    fn matches_naive_dft_for_small_powers_of_two() {
        // Compare against the definition X[k] = sum_j x[j] * omega^(j*k)
        // for sizes up to 32 to keep runtime small.
        for logn in 1..=5 {
            let n = 1usize << logn;
            let omega =
                FieldElement::primitive_root_of_unity(n as u32).unwrap();

            // a couple of random trials per size
            for _ in 0..3 {
                let coeffs: Vec<FieldElement> = random_elements(n);
                let mut got = coeffs.clone();
                ntt_in_place(&mut got, Transform::Forward).unwrap();

                let mut expect = vec![FieldElement::ZERO; n];
                for k in 0..n {
                    let mut acc = FieldElement::ZERO;
                    for j in 0..n {
                        let power = (j * k) as u32;
                        acc += coeffs[j] * omega.mod_pow_u32(power);
                    }
                    expect[k] = acc;
                }
                assert_eq!(expect, got, "mismatch at n={n}");
            }
        }
    }

    // --- Round-trip property ---

    #[test]
    fn forward_then_inverse_is_identity_up_to_512() {
        for logn in 0..=9 {
            let n = 1usize << logn;
            let mut v: Vec<FieldElement> = random_elements(n);
            let before = v.clone();
            ntt_in_place(&mut v, Transform::Forward).unwrap();
            ntt_in_place(&mut v, Transform::Inverse).unwrap();
            assert_eq!(before, v, "round-trip failed at n={n}");
        }
    }

    // --- Error paths ---

    #[test]
    fn non_power_of_two_length_is_error_and_input_not_mutated() {
        let mut v: Vec<FieldElement> = random_elements(6); // 6 is not a power of two
        let before = v.clone();
        let err = ntt_in_place(&mut v, Transform::Forward);
        assert!(err.is_err());
        assert_eq!(before, v, "input should not be mutated on error");
    }

    #[test]
    fn missing_primitive_root_returns_error_for_16384() {
        // For FieldElement with p = 8380417, p-1 = 2^13 * 3 * 11 * 31, so
        // there is no primitive 2^14-th root of unity. Length 16384 should fail.
        let mut v = vec![FieldElement::ZERO; 1 << 14]; // 16384
        assert!(ntt_in_place(&mut v, Transform::Forward).is_err());
        assert!(ntt_in_place(&mut v, Transform::Inverse).is_err());
    }

    // --- Sanity check against public wrappers ---

    #[test]
    fn agrees_with_public_wrappers_on_valid_sizes() {
        for logn in 0..=9 {
            let n = 1usize << logn;
            let mut x1: Vec<FieldElement> = random_elements(n);
            let mut x2 = x1.clone();

            ntt_in_place(&mut x1, Transform::Forward).unwrap();
            ntt::<FieldElement>(&mut x2);
            assert_eq!(x1, x2, "forward mismatch at n={n}");

            ntt_in_place(&mut x1, Transform::Inverse).unwrap();
            intt::<FieldElement>(&mut x2);
            assert_eq!(x1, x2, "inverse mismatch at n={n}");
        }
    }
}
