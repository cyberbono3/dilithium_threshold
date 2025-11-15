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
