use rand::prelude::*;
use sha2::digest::core_api::XofReaderCoreWrapper;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256, Shake256ReaderCore,
};


pub fn get_randomness(randomness: Option<&[u8]>) -> Vec<u8> {
    match randomness {
        Some(r) => r.to_vec(),
        None => {
            let mut rng = rand::thread_rng();
            let mut bytes = vec![0u8; 32];
            rng.fill_bytes(&mut bytes);
            bytes
        }
    }
}

pub fn get_hash_reader(
    message: &[u8],
) -> XofReaderCoreWrapper<Shake256ReaderCore> {
    let mut hasher = Shake256::default();
    hasher.update(message);
    hasher.finalize_xof()
}

pub fn hash_message(message: &[u8]) -> Vec<u8> {
    let mut reader = get_hash_reader(message);
    let mut output = vec![0u8; 64];
    reader.read(&mut output);
    output
}


#[cfg(test)]
mod tests {
    use super::*;

    mod get_randomness_tests {
        use super::*;

        #[test]
        fn test_get_randomness_with_provided_bytes() {
            let input = vec![1, 2, 3, 4, 5];
            let result = get_randomness(Some(&input));
            assert_eq!(result, input);
        }

        #[test]
        fn test_get_randomness_with_empty_slice() {
            let input = vec![];
            let result = get_randomness(Some(&input));
            assert_eq!(result, vec![]);
        }

        #[test]
        fn test_get_randomness_without_input() {
            let result1 = get_randomness(None);
            let result2 = get_randomness(None);

            // Should return 32 bytes
            assert_eq!(result1.len(), 32);
            assert_eq!(result2.len(), 32);

            // Should be different (extremely unlikely to be the same)
            assert_ne!(result1, result2);
        }

        #[test]
        fn test_get_randomness_deterministic_with_seed() {
            let seed = vec![42; 32];
            let result1 = get_randomness(Some(&seed));
            let result2 = get_randomness(Some(&seed));

            // Same input should give same output
            assert_eq!(result1, result2);
            assert_eq!(result1, seed);
        }
    }

    // Integration tests
    mod integration_tests {
        use super::*;

        #[test]
        fn test_randomness_and_hashing_together() {
            // Test that randomness generation and hashing work well together
            let random_bytes = get_randomness(None);
            let hash1 = hash_message(&random_bytes);
            let hash2 = hash_message(&random_bytes);

            // Same input should give same hash
            assert_eq!(hash1, hash2);

            // Different random bytes should give different hash
            let other_random = get_randomness(None);
            let other_hash = hash_message(&other_random);
            assert_ne!(hash1, other_hash);
        }
    }
}
