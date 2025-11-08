//! Miscellaneous helpers shared across math modules.

/// Compute the maximum infinity norm given an iterator of per-polynomial norms.
pub fn max_infinity_norm_from_values(
    norms: impl IntoIterator<Item = u32>,
) -> u32 {
    norms.into_iter().max().unwrap_or(0)
}
