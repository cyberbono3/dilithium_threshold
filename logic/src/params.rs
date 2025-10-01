pub const DEFAULT_SECURITY_LEVEL: usize = 2;

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


#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DilithiumConfig {
    pub k: usize,
    pub l: usize,
    pub eta: u32,
    pub tau: usize,
    pub beta: u32,
    pub gamma1: u32,
    pub gamma2: u32,
    pub d: usize,
}

impl Default for DilithiumConfig {
    fn default() -> Self {
        Self {
            k: 4,
            l: 4,
            eta: 2,
            tau: 39,
            beta: 78,
            gamma1: 131072,
            gamma2: 95232,
            d: 13,
        }
    }
}

impl DilithiumConfig {
    pub fn new(security_level: usize) -> Self {
        match security_level {
            2 => Self::default(),
            3 => Self {
                k: 6,
                l: 5,
                eta: 4,
                tau: 49,
                beta: 196,
                gamma1: 524288,
                gamma2: 261888,
                d: 13,
            },
            5 => Self {
                k: 8,
                l: 7,
                eta: 2,
                tau: 60,
                beta: 120,
                gamma1: 524288,
                gamma2: 261888,
                d: 13,
            },
            _ => panic!("Invalid security level"),
        }
    }
}

/// Validate threshold config
pub fn validate_threshold_config(
    threshold: usize,
    participants: usize,
) -> bool {
    // Threshold must be at least 2 (no single point of failure)
    // Threshold cannot exceed number of participant
    // Reasonable upper limit on participants. Use as a practical limit for the number of participants
    // Participants must be at least 2 (otherwise threshold >= 2 would be impossible)
    (2..=participants).contains(&threshold) && (2..N).contains(&participants)
}



#[cfg(test)]
mod tests {
    use super::*;

    // Tests for DilithiumConfig
    #[test]
    fn test_dilithium_config_default() {
        let config = DilithiumConfig::default();
        assert_eq!(config.k, 4);
        assert_eq!(config.l, 4);
        assert_eq!(config.eta, 2);
        assert_eq!(config.tau, 39);
        assert_eq!(config.beta, 78);
        assert_eq!(config.gamma1, 131072);
        assert_eq!(config.gamma2, 95232);
        assert_eq!(config.d, 13);
    }

    #[test]
    fn test_dilithium_config_security_level_2() {
        let config = DilithiumConfig::new(2);
        // Security level 2 should match default
        assert_eq!(config.k, 4);
        assert_eq!(config.l, 4);
        assert_eq!(config.eta, 2);
        assert_eq!(config.tau, 39);
        assert_eq!(config.beta, 78);
        assert_eq!(config.gamma1, 131072);
        assert_eq!(config.gamma2, 95232);
        assert_eq!(config.d, 13);

        // Verify it equals default
        assert_eq!(config, DilithiumConfig::default());
    }

    #[test]
    fn test_dilithium_config_security_level_3() {
        let config = DilithiumConfig::new(3);
        assert_eq!(config.k, 6);
        assert_eq!(config.l, 5);
        assert_eq!(config.eta, 4);
        assert_eq!(config.tau, 49);
        assert_eq!(config.beta, 196);
        assert_eq!(config.gamma1, 524288);
        assert_eq!(config.gamma2, 261888);
        assert_eq!(config.d, 13);
    }

    #[test]
    fn test_dilithium_config_security_level_5() {
        let config = DilithiumConfig::new(5);
        assert_eq!(config.k, 8);
        assert_eq!(config.l, 7);
        assert_eq!(config.eta, 2);
        assert_eq!(config.tau, 60);
        assert_eq!(config.beta, 120);
        assert_eq!(config.gamma1, 524288);
        assert_eq!(config.gamma2, 261888);
        assert_eq!(config.d, 13);
    }

    #[test]
    #[should_panic(expected = "Invalid security level")]
    fn test_dilithium_config_invalid_security_level_0() {
        DilithiumConfig::new(0);
    }

    #[test]
    #[should_panic(expected = "Invalid security level")]
    fn test_dilithium_config_invalid_security_level_1() {
        DilithiumConfig::new(1);
    }

    #[test]
    #[should_panic(expected = "Invalid security level")]
    fn test_dilithium_config_invalid_security_level_4() {
        DilithiumConfig::new(4);
    }

    #[test]
    #[should_panic(expected = "Invalid security level")]
    fn test_dilithium_config_invalid_security_level_6() {
        DilithiumConfig::new(6);
    }

    #[test]
    #[should_panic(expected = "Invalid security level")]
    fn test_dilithium_config_invalid_security_level_large() {
        DilithiumConfig::new(100);
    }

    #[test]
    fn test_dilithium_config_copy() {
        let config1 = DilithiumConfig::new(5);
        let config2 = config1; // Copy semantics
        assert_eq!(config1, config2);
    }

    #[test]
    fn test_dilithium_config_debug() {
        let config = DilithiumConfig::default();
        let debug_str = format!("{:?}", config);
        assert!(debug_str.contains("DilithiumConfig"));
        assert!(debug_str.contains("k: 4"));
        assert!(debug_str.contains("l: 4"));
    }

    #[test]
    fn test_dilithium_config_partial_eq() {
        let config1 = DilithiumConfig::new(2);
        let config2 = DilithiumConfig::new(2);
        let config3 = DilithiumConfig::new(3);

        assert_eq!(config1, config2);
        assert_ne!(config1, config3);
        assert_ne!(config2, config3);
    }

    #[test]
    fn test_default_security_level_constant() {
        assert_eq!(DEFAULT_SECURITY_LEVEL, 2);
        // Verify that default security level creates default config
        let config_from_const = DilithiumConfig::new(DEFAULT_SECURITY_LEVEL);
        let config_default = DilithiumConfig::default();
        assert_eq!(config_from_const, config_default);
    }

    // Tests for validate_threshold_config
    #[test]
    fn test_valid_threshold_configs() {
        // Valid configurations where 2 <= threshold <= participants and 2 <= participants < N
        assert!(validate_threshold_config(2, 2)); // Minimum valid case
        assert!(validate_threshold_config(2, 3));
        assert!(validate_threshold_config(3, 3)); // threshold == participants
        assert!(validate_threshold_config(2, 10));
        assert!(validate_threshold_config(5, 10));
        assert!(validate_threshold_config(10, 10));
        assert!(validate_threshold_config(2, 100));
        assert!(validate_threshold_config(50, 100));
        assert!(validate_threshold_config(100, 100));
    }

    #[test]
    fn test_threshold_too_small() {
        // Threshold must be at least 2
        assert!(!validate_threshold_config(0, 5));
        assert!(!validate_threshold_config(1, 5));
        assert!(!validate_threshold_config(0, 2));
        assert!(!validate_threshold_config(1, 2));
    }

    #[test]
    fn test_threshold_exceeds_participants() {
        // Threshold cannot exceed participants
        assert!(!validate_threshold_config(3, 2));
        assert!(!validate_threshold_config(11, 10));
        assert!(!validate_threshold_config(101, 100));
    }

    #[test]
    fn test_participants_too_small() {
        // Participants must be at least 2
        assert!(!validate_threshold_config(2, 0));
        assert!(!validate_threshold_config(2, 1));
        // Even if threshold is invalid, participants being < 2 should fail
        assert!(!validate_threshold_config(0, 0));
        assert!(!validate_threshold_config(1, 1));
    }

    #[test]
    fn test_edge_cases() {
        // Test boundary conditions
        assert!(validate_threshold_config(2, 2)); // Both at minimum

        // Test invalid combinations
        assert!(!validate_threshold_config(2, 1)); // participants < 2
        assert!(!validate_threshold_config(0, 0)); // both invalid

        // Note: Large values like 1000 may fail if N < 1000
        // The actual upper limit depends on the value of N from math::prelude::N
        // Only test with smaller values that are likely below N
        assert!(validate_threshold_config(50, 50));
        assert!(validate_threshold_config(49, 50));
    }

    #[test]
    fn test_typical_threshold_scenarios() {
        // Common threshold scenarios in practice

        // 2-of-3 multisig
        assert!(validate_threshold_config(2, 3));

        // 3-of-5 multisig
        assert!(validate_threshold_config(3, 5));

        // Majority threshold (n/2 + 1)
        assert!(validate_threshold_config(6, 11)); // 6 of 11
        assert!(validate_threshold_config(51, 100)); // 51 of 100

        // Super-majority (2n/3 + 1)
        assert!(validate_threshold_config(7, 10)); // 7 of 10
        assert!(validate_threshold_config(67, 100)); // 67 of 100
    }
}

