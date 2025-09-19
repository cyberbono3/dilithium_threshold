/// Round-3 / ML-DSA-44 (a.k.a. Dilithium-2) parameter set.
/// See Table 2 in the Round-3 specification. (Educational; uncompressed variant.)
pub const N: usize = 256;
pub const Q: i64 = 8_380_417; // 2^23 - 2^13 + 1
pub const D: u32 = 13; // dropped bits for t in compressed scheme (unused here but kept for reference)
pub const TAU: usize = 39; // challenge weight (±1s)
pub const ETA: i64 = 2; // SK coeff range via CBD(η=2)
pub const K: usize = 4; // A is K x L
pub const L: usize = 4;

pub const GAMMA1: i64 = 1 << 17; // 131072
pub const GAMMA2: i64 = (Q - 1) / 88; // 95_232
pub const ALPHA: i64 = 2 * GAMMA2; // 190_464
pub const BETA: i64 = (TAU as i64) * ETA; // 78
