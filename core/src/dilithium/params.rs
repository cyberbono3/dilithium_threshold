use std::fmt;

/// Default security level used throughout the logic layer (Dilithium-2).
pub const DEFAULT_SECURITY_LEVEL: usize = SecurityLevel::Level2.as_usize();

/// Supported Dilithium security levels as defined in the ML-DSA specification.
pub const SUPPORTED_SECURITY_LEVELS: [SecurityLevel; 3] = [
    SecurityLevel::Level2,
    SecurityLevel::Level3,
    SecurityLevel::Level5,
];

/// Enumerates the supported ML-DSA / Dilithium security levels.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(usize)]
pub enum SecurityLevel {
    Level2 = 2,
    Level3 = 3,
    Level5 = 5,
}

impl SecurityLevel {
    /// Return the numeric identifier used by the Dilithium specification.
    #[inline]
    pub const fn as_usize(self) -> usize {
        self as usize
    }

    /// Retrieve the parameter set associated with this security level.
    #[inline]
    pub const fn config(self) -> DilithiumConfig {
        match self {
            SecurityLevel::Level2 => DILITHIUM_L2_CONFIG,
            SecurityLevel::Level3 => DILITHIUM_L3_CONFIG,
            SecurityLevel::Level5 => DILITHIUM_L5_CONFIG,
        }
    }
}

impl TryFrom<usize> for SecurityLevel {
    type Error = InvalidSecurityLevel;

    fn try_from(value: usize) -> Result<Self, Self::Error> {
        match value {
            x if x == SecurityLevel::Level2 as usize => {
                Ok(SecurityLevel::Level2)
            }
            x if x == SecurityLevel::Level3 as usize => {
                Ok(SecurityLevel::Level3)
            }
            x if x == SecurityLevel::Level5 as usize => {
                Ok(SecurityLevel::Level5)
            }
            other => Err(InvalidSecurityLevel(other)),
        }
    }
}

/// Error raised when an unsupported security level is requested.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InvalidSecurityLevel(pub usize);

impl fmt::Display for InvalidSecurityLevel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "unsupported Dilithium security level {}", self.0)
    }
}

impl std::error::Error for InvalidSecurityLevel {}

/// Dilithium parameter collection associated with a concrete security level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DilithiumConfig {
    pub n: usize,
    pub q: i64,
    pub d: u32,
    pub k: usize,
    pub l: usize,
    pub eta: u32,
    pub tau: usize,
    pub beta: u32,
    pub gamma1: u32,
    pub gamma2: u32,
}

/// Helper container with literal parameters used to build [`DilithiumConfig`] instances.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DilithiumConfigValues {
    pub n: usize,
    pub q: i64,
    pub d: u32,
    pub k: usize,
    pub l: usize,
    pub eta: u32,
    pub tau: usize,
    pub beta: u32,
    pub gamma1: u32,
    pub gamma2: u32,
}

impl DilithiumConfig {
    /// Create a new parameter set from literal values (const-friendly).
    pub const fn from_values(values: DilithiumConfigValues) -> Self {
        Self {
            n: values.n,
            q: values.q,
            d: values.d,
            k: values.k,
            l: values.l,
            eta: values.eta,
            tau: values.tau,
            beta: values.beta,
            gamma1: values.gamma1,
            gamma2: values.gamma2,
        }
    }

    /// Return the parameter set for a given [`SecurityLevel`].
    #[inline]
    pub const fn for_level(level: SecurityLevel) -> Self {
        level.config()
    }

    /// Convenience helper that preserves the historical `usize` API.
    ///
    /// Returns an error when `security_level` does not correspond to a supported set.
    pub fn new(security_level: usize) -> Result<Self, InvalidSecurityLevel> {
        Self::try_from(security_level)
    }

    /// Enumerate the security levels known at compile time.
    pub const fn supported_levels() -> &'static [SecurityLevel; 3] {
        &SUPPORTED_SECURITY_LEVELS
    }
}

impl Default for DilithiumConfig {
    fn default() -> Self {
        SecurityLevel::Level2.config()
    }
}

impl TryFrom<usize> for DilithiumConfig {
    type Error = InvalidSecurityLevel;

    fn try_from(value: usize) -> Result<Self, Self::Error> {
        SecurityLevel::try_from(value).map(DilithiumConfig::for_level)
    }
}

/// Round-3 / ML-DSA-44 (Dilithium-2) parameter set (educational, uncompressed).
const DILITHIUM_L2_CONFIG: DilithiumConfig =
    DilithiumConfig::from_values(DilithiumConfigValues {
        n: 256,
        q: 8_380_417,
        d: 13,
        k: 4,
        l: 4,
        eta: 2,
        tau: 39,
        beta: 78,
        gamma1: 131_072,
        gamma2: 95_232,
    });

/// Round-3 / ML-DSA-65 (Dilithium-3) parameter set (educational, uncompressed).
const DILITHIUM_L3_CONFIG: DilithiumConfig =
    DilithiumConfig::from_values(DilithiumConfigValues {
        n: 256,
        q: 8_380_417,
        d: 13,
        k: 6,
        l: 5,
        eta: 4,
        tau: 49,
        beta: 196,
        gamma1: 524_288,
        gamma2: 261_888,
    });

/// Round-3 / ML-DSA-87 (Dilithium-5) parameter set (educational, uncompressed).
const DILITHIUM_L5_CONFIG: DilithiumConfig =
    DilithiumConfig::from_values(DilithiumConfigValues {
        n: 256,
        q: 8_380_417,
        d: 13,
        k: 8,
        l: 7,
        eta: 2,
        tau: 60,
        beta: 120,
        gamma1: 524_288,
        gamma2: 261_888,
    });

/// Default parameter set backing the module-level constants.
pub const DEFAULT_CONFIG: DilithiumConfig = DILITHIUM_L2_CONFIG;

pub const N: usize = DEFAULT_CONFIG.n;
pub const Q: i64 = DEFAULT_CONFIG.q;
pub const D: u32 = DEFAULT_CONFIG.d;
pub const K: usize = DEFAULT_CONFIG.k;
pub const L: usize = DEFAULT_CONFIG.l;
pub const ETA: i64 = DEFAULT_CONFIG.eta as i64;
pub const TAU: usize = DEFAULT_CONFIG.tau;
pub const BETA: i64 = DEFAULT_CONFIG.beta as i64;
pub const GAMMA1: i64 = DEFAULT_CONFIG.gamma1 as i64;
pub const GAMMA2: i64 = DEFAULT_CONFIG.gamma2 as i64;
pub const ALPHA: i64 = 2 * GAMMA2;

/// Validate the relation between threshold and participant counts.
pub fn validate_threshold_config(
    threshold: usize,
    participants: usize,
) -> bool {
    (2..=participants).contains(&threshold) && (2..N).contains(&participants)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EXPECTED_LEVEL_CONFIGS: [(SecurityLevel, DilithiumConfig); 3] = [
        (SecurityLevel::Level2, DILITHIUM_L2_CONFIG),
        (SecurityLevel::Level3, DILITHIUM_L3_CONFIG),
        (SecurityLevel::Level5, DILITHIUM_L5_CONFIG),
    ];

    #[test]
    fn default_config_matches_level2() {
        assert_eq!(DilithiumConfig::default(), DILITHIUM_L2_CONFIG);
        assert_eq!(DEFAULT_CONFIG, DILITHIUM_L2_CONFIG);
        assert_eq!(DEFAULT_SECURITY_LEVEL, SecurityLevel::Level2.as_usize());
    }

    #[test]
    fn security_level_iteration_is_complete() {
        let from_const: Vec<_> = SUPPORTED_SECURITY_LEVELS.to_vec();
        let from_helper: Vec<_> = DilithiumConfig::supported_levels().to_vec();
        assert_eq!(from_const, from_helper);
        assert_eq!(from_const.len(), EXPECTED_LEVEL_CONFIGS.len());
        for level in &from_const {
            assert!(
                EXPECTED_LEVEL_CONFIGS
                    .iter()
                    .any(|(expected_level, _)| expected_level == level),
                "level {:?} missing in expectations",
                level
            );
        }
    }

    #[test]
    fn config_lookup_per_security_level() {
        for (level, expected) in EXPECTED_LEVEL_CONFIGS {
            assert_eq!(level.config(), expected);
            assert_eq!(DilithiumConfig::for_level(level), expected);
            assert_eq!(
                DilithiumConfig::try_from(level.as_usize()).unwrap(),
                expected
            );
        }
    }

    #[test]
    fn invalid_security_levels_are_rejected() {
        for level in [0, 1, 4, 6, 100] {
            assert!(SecurityLevel::try_from(level).is_err());
            assert!(DilithiumConfig::try_from(level).is_err());
        }
    }

    #[test]
    fn new_rejects_invalid_levels() {
        for level in [0, 1, 4, 6, 100] {
            assert!(
                DilithiumConfig::new(level).is_err(),
                "expected level {level} to be rejected"
            );
        }
    }

    #[test]
    fn debug_includes_field_information() {
        let config = DilithiumConfig::default();
        let debug_str = format!("{config:?}");
        assert!(debug_str.contains("DilithiumConfig"));
        assert!(debug_str.contains("k: 4"));
        assert!(debug_str.contains("l: 4"));
    }

    #[test]
    fn default_security_level_constant_round_trip() {
        let config_from_const =
            DilithiumConfig::new(DEFAULT_SECURITY_LEVEL).unwrap();
        assert_eq!(config_from_const, DilithiumConfig::default());
    }

    #[test]
    fn validate_threshold_config_accepts_expected_inputs() {
        assert!(validate_threshold_config(2, 2));
        assert!(validate_threshold_config(3, 3));
        assert!(validate_threshold_config(2, 10));
        assert!(validate_threshold_config(5, 10));
        assert!(validate_threshold_config(10, 10));
        assert!(validate_threshold_config(2, 100));
        assert!(validate_threshold_config(50, 100));
        assert!(validate_threshold_config(100, 100));
    }

    #[test]
    fn validate_threshold_config_rejects_small_thresholds() {
        assert!(!validate_threshold_config(0, 5));
        assert!(!validate_threshold_config(1, 5));
        assert!(!validate_threshold_config(0, 2));
        assert!(!validate_threshold_config(1, 2));
    }

    #[test]
    fn validate_threshold_config_rejects_threshold_above_participants() {
        assert!(!validate_threshold_config(3, 2));
        assert!(!validate_threshold_config(11, 10));
        assert!(!validate_threshold_config(101, 100));
    }

    #[test]
    fn validate_threshold_config_rejects_small_participant_sets() {
        assert!(!validate_threshold_config(2, 0));
        assert!(!validate_threshold_config(2, 1));
        assert!(!validate_threshold_config(0, 0));
        assert!(!validate_threshold_config(1, 1));
    }

    #[test]
    fn validate_threshold_config_edge_cases() {
        assert!(validate_threshold_config(2, 2));
        assert!(!validate_threshold_config(2, 1));
        assert!(!validate_threshold_config(0, 0));
        assert!(validate_threshold_config(50, 50));
        assert!(validate_threshold_config(49, 50));
    }

    #[test]
    fn validate_threshold_config_typical_scenarios() {
        assert!(validate_threshold_config(2, 3));
        assert!(validate_threshold_config(3, 5));
        assert!(validate_threshold_config(6, 11));
        assert!(validate_threshold_config(51, 100));
        assert!(validate_threshold_config(7, 10));
        assert!(validate_threshold_config(67, 100));
    }
}
