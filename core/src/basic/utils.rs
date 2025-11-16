use crate::dilithium::params::{ALPHA, GAMMA1, K, L, N, Q};
use crate::dilithium::utils::{
    derive_challenge_polynomial, shake256_squeezed, zero_polyvec,
};
use math::prelude::FieldElement;
use math::{poly::Polynomial, traits::FiniteField};

const HINT_MODULUS: i64 = (Q - 1) / ALPHA;

/// Split a polynomial into its high and low components.
pub(crate) fn poly_high_low<FF: FiniteField + From<i64>>(
    p: &Polynomial<'_, FF>,
) -> (Polynomial<'static, FF>, Polynomial<'static, FF>)
where
    i64: From<FF>,
{
    let (hi, lo): (Vec<_>, Vec<_>) = padded_coefficients(p)
        .map(|coeff| high_low_bits(i64::from(coeff)))
        .unzip();
    (hi.into(), lo.into())
}

/// Extract only the high component of a polynomial used for hashing.
pub(crate) fn poly_high<FF: FiniteField + From<i64>>(
    p: &Polynomial<'_, FF>,
) -> Polynomial<'static, FF>
where
    i64: From<FF>,
{
    let (high, _) = poly_high_low(p);
    high
}

/// Extract only the low component of a polynomial.
pub(crate) fn poly_low<FF: FiniteField + From<i64>>(
    p: &Polynomial<'_, FF>,
) -> Polynomial<'static, FF>
where
    i64: From<FF>,
{
    let (_, low) = poly_high_low(p);
    low
}

/// Serialize the w1 vector into bytes that feed the challenge hash.
pub(crate) fn pack_w1_for_hash<
    FF: FiniteField + Into<[u8; FieldElement::BYTES]>,
>(
    w1: &[Polynomial<'_, FF>; K],
) -> Vec<u8> {
    let mut out = Vec::with_capacity(K * N * 2);
    for poly in w1 {
        for v in poly.coefficients() {
            let arr: [u8; FieldElement::BYTES] = (*v).into();
            out.extend_from_slice(&arr);
        }
    }
    out
}

/// Sample the y polynomial vector for a given message seed and counter.
pub(crate) fn sample_y<FF: FiniteField + From<i64>>(
    seed: &[u8],
    ctr: u32,
) -> [Polynomial<'static, FF>; L] {
    let mut out = zero_polyvec::<L, FF>();

    for (j, slot) in out.iter_mut().enumerate() {
        let idx_bytes = (j as u16).to_le_bytes();
        let ctr_bytes = ctr.to_le_bytes();
        let bytes = shake256_squeezed(seed, &[&idx_bytes, &ctr_bytes], 3 * N);
        let mut coeffs = [0i64; N];
        for (i, chunk) in bytes.chunks_exact(3).take(N).enumerate() {
            let v = (chunk[0] as i64)
                | ((chunk[1] as i64) << 8)
                | ((chunk[2] as i64) << 16);
            let modulus = 2 * GAMMA1 + 1;
            let centered = (v % modulus + modulus) % modulus - GAMMA1;
            coeffs[i] = centered;
        }
        *slot = coeffs.into();
    }
    out
}

/// Expand a 32-byte seed into an η=2 vector of secret polynomials.
pub(crate) fn expand_secret_vector<const LEN: usize, FF>(
    seed: &[u8; 32],
) -> [Polynomial<'static, FF>; LEN]
where
    FF: FiniteField + From<i64>,
{
    let mut secrets = zero_polyvec::<LEN, FF>();
    for (index, slot) in secrets.iter_mut().enumerate() {
        let idx_bytes = (index as u16).to_le_bytes();
        let bytes = shake256_squeezed(seed, &[&idx_bytes], 2 * N);
        *slot = cbd_eta2::<FF>(&bytes);
    }
    secrets
}

/// CBD for η=2 from an XOF stream.
pub(crate) fn cbd_eta2<FF: FiniteField + From<i64>>(
    stream: &[u8],
) -> Polynomial<'static, FF> {
    let mut out = [0i64; N];
    let mut bitpos = 0usize;
    for out_i in out.iter_mut().take(N) {
        let byte_idx = bitpos / 8;
        let two =
            u16::from_le_bytes([stream[byte_idx], stream[byte_idx + 1]]) as u32;
        let bits = two;
        let a0 = (bits & 1) + ((bits >> 1) & 1);
        let a1 = ((bits >> 2) & 1) + ((bits >> 3) & 1);
        *out_i = (a0 as i64) - (a1 as i64);
        bitpos += 4;
    }
    out.into()
}

/// Derive the sparse challenge polynomial from message and w1 bytes.
pub(crate) fn derive_challenge<FF: FiniteField + From<i64>>(
    message: &[u8],
    w1_pack: &[u8],
) -> Polynomial<'static, FF> {
    let mut seed = Vec::with_capacity(message.len() + w1_pack.len());
    seed.extend_from_slice(message);
    seed.extend_from_slice(w1_pack);
    derive_challenge_polynomial::<FF>(&seed)
}

#[inline]
fn polyvec_add_scaled_with_sign<FF: FiniteField + 'static, const LEN: usize>(
    base: &[Polynomial<'static, FF>; LEN],
    scale: &Polynomial<'_, FF>,
    mult: &[Polynomial<'static, FF>; LEN],
    negate: bool,
) -> [Polynomial<'static, FF>; LEN] {
    let mut dest = zero_polyvec::<LEN, FF>();
    for ((slot, base_poly), mult_poly) in dest.iter_mut().zip(base).zip(mult) {
        slot.clone_from(base_poly);
        let mut scaled = mult_poly * scale;
        if negate {
            scaled = -scaled;
        }
        *slot += scaled;
    }
    dest
}

/// Compute element-wise addition of `base` and `scale * mult`.
#[inline]
pub(crate) fn polyvec_add_scaled<
    FF: FiniteField + 'static,
    const LEN: usize,
>(
    base: &[Polynomial<'static, FF>; LEN],
    scale: &Polynomial<'_, FF>,
    mult: &[Polynomial<'static, FF>; LEN],
) -> [Polynomial<'static, FF>; LEN] {
    polyvec_add_scaled_with_sign(base, scale, mult, false)
}

/// Compute element-wise subtraction of `scale * mult` from `base`.
#[inline]
pub(crate) fn polyvec_sub_scaled<
    FF: FiniteField + 'static,
    const LEN: usize,
>(
    base: &[Polynomial<'static, FF>; LEN],
    scale: &Polynomial<'_, FF>,
    mult: &[Polynomial<'static, FF>; LEN],
) -> [Polynomial<'static, FF>; LEN] {
    polyvec_add_scaled_with_sign(base, scale, mult, true)
}

/// Generate per-coefficient hints enabling reconstruction of the reference high bits.
pub(crate) fn make_hints<'a, FF>(
    target_high: &[Polynomial<'a, FF>; K],
    base: &[Polynomial<'a, FF>; K],
) -> [Polynomial<'static, FF>; K]
where
    FF: FiniteField + From<i64>,
    i64: From<FF>,
{
    std::array::from_fn(|idx| make_hint_poly(&target_high[idx], &base[idx]))
}

/// Apply stored hints to recover the original high-bit decomposition.
pub(crate) fn use_hints<'hint, 'poly, FF>(
    hints: &[Polynomial<'hint, FF>; K],
    base: &[Polynomial<'poly, FF>; K],
) -> [Polynomial<'static, FF>; K]
where
    FF: FiniteField + From<i64>,
    i64: From<FF>,
{
    std::array::from_fn(|idx| use_hint_poly(&hints[idx], &base[idx]))
}

/// Check whether each polynomial's infinity norm stays below `bound`.
#[inline]
pub(crate) fn all_infty_norm_below<FF: FiniteField, const LEN: usize>(
    polys: &[Polynomial<'_, FF>; LEN],
    bound: i64,
) -> bool
where
    i64: From<FF>,
{
    polys.iter().all(|p| (p.norm_infinity() as i64) < bound)
}

fn coeff_at_or_zero<FF: FiniteField>(
    poly: &Polynomial<'_, FF>,
    idx: usize,
) -> FF {
    poly.coefficients().get(idx).copied().unwrap_or(FF::ZERO)
}

fn make_hint_poly<FF>(
    target_high: &Polynomial<'_, FF>,
    base: &Polynomial<'_, FF>,
) -> Polynomial<'static, FF>
where
    FF: FiniteField + From<i64>,
    i64: From<FF>,
{
    let mut coeffs = vec![FF::ZERO; N];
    for (idx, coeff) in coeffs.iter_mut().enumerate() {
        let target_hi = i64::from(coeff_at_or_zero(target_high, idx));
        let base_value = i64::from(coeff_at_or_zero(base, idx));
        let (base_hi, _) = high_low_bits(base_value);
        if target_hi != base_hi {
            *coeff = FF::ONE;
        }
    }
    coeffs.into()
}

fn use_hint_poly<FF>(
    hint: &Polynomial<'_, FF>,
    base: &Polynomial<'_, FF>,
) -> Polynomial<'static, FF>
where
    FF: FiniteField + From<i64>,
    i64: From<FF>,
{
    let mut coeffs = vec![FF::ZERO; N];
    for (idx, coeff) in coeffs.iter_mut().enumerate() {
        let hint_bit = !coeff_at_or_zero(hint, idx).is_zero();
        let base_value = i64::from(coeff_at_or_zero(base, idx));
        let (mut hi, lo) = high_low_bits(base_value);
        if hint_bit {
            if lo > 0 {
                hi += 1;
            } else {
                hi -= 1;
            }
        }
        let modulus = HINT_MODULUS;
        hi = ((hi % modulus) + modulus) % modulus;
        *coeff = FF::from(hi);
    }
    coeffs.into()
}

/// Split a coefficient into high and low parts relative to ALPHA.
fn high_low_bits(x: i64) -> (i64, i64) {
    let w1 = x / ALPHA;
    let mut w0 = x - w1 * ALPHA;
    if w0 > ALPHA / 2 {
        w0 -= ALPHA;
    }
    (w1, w0)
}

/// Iterate coefficients padded with zeros up to length N.
fn padded_coefficients<'a, FF: FiniteField>(
    poly: &'a Polynomial<'a, FF>,
) -> impl Iterator<Item = FF> + 'a {
    poly.coefficients()
        .iter()
        .copied()
        .chain(std::iter::repeat(FF::ZERO))
        .take(N)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dilithium::params::{ALPHA, GAMMA1, K, TAU};
    use crate::dilithium::utils::zero_polyvec;
    use math::field_element::FieldElement;

    mod hint_tests {
        use super::*;

        #[test]
        fn hints_restore_reference_high_bits() {
            let mut reference = zero_polyvec::<K, FieldElement>();
            let mut base = zero_polyvec::<K, FieldElement>();
            // Choose coefficients that differ by a single ALPHA window.
            reference[0] = vec![FieldElement::from(ALPHA + 5)].into();
            base[0] = vec![FieldElement::from(5)].into();

            let target_high =
                std::array::from_fn(|idx| poly_high(&reference[idx]));
            let hints = make_hints(&target_high, &base);
            let recovered = use_hints(&hints, &base);

            assert_eq!(recovered[0], target_high[0]);
        }
    }

    mod sample_y_tests {
        use super::*;

        /// Convert a field element into a centered signed integer.
        fn centered_value(fe: FieldElement) -> i64 {
            let mut v = i64::from(fe);
            let p = FieldElement::P as i64;
            if v > p / 2 {
                v -= p;
            }
            v
        }

        /// Return true if all coefficients stay within the ±GAMMA1 bound.
        fn coeffs_within_bounds(
            polys: &[Polynomial<'_, FieldElement>],
        ) -> bool {
            polys
                .iter()
                .flat_map(|poly| poly.coefficients())
                .map(|&c| centered_value(c).abs())
                .all(|abs| abs <= GAMMA1)
        }

        /// Sampling with identical seed and counter should be deterministic.
        #[test]
        fn deterministic_per_seed_and_counter() {
            let seed = [0x11u8; 32];
            let first = sample_y::<FieldElement>(&seed, 7);
            let second = sample_y::<FieldElement>(&seed, 7);
            assert_eq!(first, second);
        }

        /// Different counters must yield distinct sampled vectors.
        #[test]
        fn different_counters_produce_distinct_vectors() {
            let seed = [0x77u8; 32];
            let first = sample_y::<FieldElement>(&seed, 1);
            let second = sample_y::<FieldElement>(&seed, 2);
            assert_ne!(first, second);
        }

        /// All coefficients should remain centered within the permitted bound.
        #[test]
        fn coefficients_are_centered() {
            let seed = [0xC3u8; 32];
            let polys = sample_y::<FieldElement>(&seed, 5);
            assert!(coeffs_within_bounds(&polys));
        }
    }

    mod pack_w1_tests {
        use super::*;

        /// Helper to create a polynomial from integer coefficients.
        fn simple_poly(coeffs: &[i64]) -> Polynomial<'static, FieldElement> {
            Polynomial::from(
                coeffs
                    .iter()
                    .map(|&c| FieldElement::from(c))
                    .collect::<Vec<_>>(),
            )
        }

        /// Packing should be deterministic when the input polynomials are fixed.
        #[test]
        fn deterministic_for_fixed_polys() {
            let polys = [
                simple_poly(&[1, 2, 3]),
                simple_poly(&[4, 5, 6]),
                simple_poly(&[7, 8, 9]),
                simple_poly(&[10, 11, 12]),
            ];

            let first = pack_w1_for_hash(&polys);
            let second = pack_w1_for_hash(&polys);
            assert_eq!(first, second);
        }

        /// Different polynomials must lead to different packed byte arrays.
        #[test]
        fn different_inputs_produce_different_outputs() {
            let polys_a = [
                simple_poly(&[1, 2, 3]),
                simple_poly(&[4, 5, 6]),
                simple_poly(&[7, 8, 9]),
                simple_poly(&[10, 11, 12]),
            ];

            let polys_b = [
                simple_poly(&[1, 2, 3]),
                simple_poly(&[4, 5, 6]),
                simple_poly(&[7, 8, 9]),
                simple_poly(&[10, 11, 13]),
            ];

            let hash_a = pack_w1_for_hash(&polys_a);
            let hash_b = pack_w1_for_hash(&polys_b);
            assert_ne!(hash_a, hash_b);
        }
    }

    mod derive_challenge_tests {
        use super::*;

        /// Derive challenge should be deterministic for identical inputs.
        #[test]
        fn deterministic_for_same_inputs() {
            let msg = b"challenge";
            let hash = vec![0x42; 32];
            let c1 = derive_challenge::<FieldElement>(msg, &hash);
            let c2 = derive_challenge::<FieldElement>(msg, &hash);
            assert_eq!(c1, c2);
        }

        /// Different messages should produce different challenge polynomials.
        #[test]
        fn changing_message_changes_challenge() {
            let hash = vec![0x77; 64];
            let c1 = derive_challenge::<FieldElement>(b"m1", &hash);
            let c2 = derive_challenge::<FieldElement>(b"m2", &hash);
            assert_ne!(c1, c2);
        }

        /// The derived challenge must contain exactly TAU non-zero coefficients.
        #[test]
        fn challenge_has_tau_non_zero_entries() {
            use num_traits::ConstZero;

            let msg = b"nonzero-count";
            let hash = vec![0xAB; 128];
            let challenge = derive_challenge::<FieldElement>(msg, &hash);
            let non_zero = challenge
                .coefficients()
                .iter()
                .filter(|&&c| !FieldElement::ZERO.eq(&c))
                .count();
            assert_eq!(non_zero, TAU);
        }
    }

    mod infty_norm_tests {
        use super::*;

        /// Helper to construct polynomials for infinity norm tests.
        fn make_poly(coeffs: &[i64]) -> Polynomial<'static, FieldElement> {
            Polynomial::from(
                coeffs
                    .iter()
                    .map(|&c| FieldElement::from(c))
                    .collect::<Vec<_>>(),
            )
        }

        #[test]
        fn polyvec_add_scaled_adds_scaled_component() {
            let base = [
                Polynomial::from(vec![FieldElement::from(1i32)]),
                Polynomial::from(vec![FieldElement::from(2i32)]),
            ];
            let mult = [
                Polynomial::from(vec![FieldElement::from(3i32)]),
                Polynomial::from(vec![FieldElement::from(-1i32)]),
            ];
            let scale = Polynomial::from(vec![FieldElement::from(2i32)]);

            let result =
                polyvec_add_scaled::<FieldElement, 2>(&base, &scale, &mult);

            let expected = [
                Polynomial::from(vec![FieldElement::from(7i32)]),
                Polynomial::from(vec![FieldElement::from(0i32)]),
            ];

            assert_eq!(result, expected);
        }

        #[test]
        fn polyvec_sub_scaled_subtracts_scaled_component() {
            let base = [
                Polynomial::from(vec![FieldElement::from(5i32)]),
                Polynomial::from(vec![FieldElement::from(-3i32)]),
            ];
            let mult = [
                Polynomial::from(vec![FieldElement::from(2i32)]),
                Polynomial::from(vec![FieldElement::from(4i32)]),
            ];
            let scale = Polynomial::from(vec![FieldElement::from(3i32)]);

            let result =
                polyvec_sub_scaled::<FieldElement, 2>(&base, &scale, &mult);

            let expected = [
                Polynomial::from(vec![FieldElement::from(-1i32)]),
                Polynomial::from(vec![FieldElement::from(-15i32)]),
            ];

            assert_eq!(result, expected);
        }

        /// Detect cases where all polynomials stay within the supplied bound.
        #[test]
        fn detects_within_bound() {
            let polys = [make_poly(&[1, -2, 3]), make_poly(&[0, 0, 4])];
            assert!(all_infty_norm_below(&polys, 5));
        }

        /// Detect cases where the infinity norm exceeds the bound.
        #[test]
        fn detects_violation() {
            let polys = [make_poly(&[1, -2, 3]), make_poly(&[0, 0, 6])];
            assert!(!all_infty_norm_below(&polys, 5));
        }
    }
}
