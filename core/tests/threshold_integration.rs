use dilithium_core::basic::PublicKey;
use dilithium_core::dilithium::params::{K, L};
use dilithium_core::dilithium::threshold::{
    PartialSignature, ThresholdKeyShare, ThresholdSignature,
};
use dilithium_core::dilithium::utils::{
    derive_challenge_polynomial, hash_message,
};
use math::{field_element::FieldElement, poly::Polynomial};

const MESSAGE: &[u8] = b"threshold integration test message";
const RANDOMNESS: [u8; 32] = [0x42; 32];

fn generate_partials(
    ts: &ThresholdSignature,
    shares: &[ThresholdKeyShare<'static, FieldElement>],
    message: &[u8],
    randomness: &[u8],
) -> Vec<PartialSignature<'static, FieldElement>> {
    shares
        .iter()
        .map(|share| {
            ts.partial_sign::<FieldElement>(message, share, Some(randomness))
                .expect("partial signing should succeed")
        })
        .collect()
}

fn verify_combined(
    ts: &ThresholdSignature,
    shares: &[ThresholdKeyShare<'static, FieldElement>],
    message: &[u8],
) -> (
    PublicKey<'static, FieldElement>,
    Vec<PartialSignature<'static, FieldElement>>,
) {
    let public_key = shares[0].public_key.clone();
    let partials = generate_partials(
        ts,
        &shares[..ts.get_threshold_info().threshold],
        message,
        &RANDOMNESS,
    );
    (public_key, partials)
}

fn expected_challenge(message: &[u8]) -> Polynomial<'static, FieldElement> {
    let mu = hash_message(message);
    let mut seed = Vec::with_capacity(mu.len() + b"challenge".len());
    seed.extend_from_slice(&mu);
    seed.extend_from_slice(b"challenge");
    derive_challenge_polynomial::<FieldElement>(&seed)
}

fn combinations<T: Clone>(items: &[T], k: usize) -> Vec<Vec<T>> {
    fn helper<T: Clone>(
        items: &[T],
        k: usize,
        start: usize,
        current: &mut Vec<T>,
        result: &mut Vec<Vec<T>>,
    ) {
        if current.len() == k {
            result.push(current.clone());
            return;
        }
        for idx in start..items.len() {
            current.push(items[idx].clone());
            helper(items, k, idx + 1, current, result);
            current.pop();
        }
    }

    let mut result = Vec::new();
    let mut current = Vec::new();
    helper(items, k, 0, &mut current, &mut result);
    result
}

#[test]
fn complete_threshold_workflow() {
    let ts = ThresholdSignature::default();
    let shares = ts
        .distributed_keygen::<FieldElement>()
        .expect("distributed key generation succeeds");
    assert_eq!(shares.len(), ts.get_threshold_info().participants);

    let (public_key, partials) = verify_combined(&ts, &shares, MESSAGE);
    let signature = ts
        .combine_signatures::<FieldElement>(&partials, &public_key)
        .expect("combination succeeds");

    assert_eq!(signature.c, expected_challenge(MESSAGE));
    assert_eq!(signature.z.len(), L);
    assert_eq!(signature.h.len(), K);
    for share in &shares {
        assert_eq!(share.public_key, public_key);
    }
}

#[test]
fn different_threshold_configurations_succeed() {
    let configs = [(2, 3), (3, 5), (4, 6)];
    for (threshold, participants) in configs {
        let ts = ThresholdSignature::new(threshold, participants)
            .expect("configuration should be valid");
        let shares = ts
            .distributed_keygen::<FieldElement>()
            .expect("distributed key generation succeeds");
        let public_key = shares[0].public_key.clone();
        let partials =
            generate_partials(&ts, &shares[..threshold], MESSAGE, &RANDOMNESS);
        let signature = ts
            .combine_signatures::<FieldElement>(&partials, &public_key)
            .expect("combination succeeds");
        assert_eq!(signature.c, expected_challenge(MESSAGE));
    }
}

#[test]
fn insufficient_signatures_fail() {
    let ts = ThresholdSignature::default();
    let shares = ts
        .distributed_keygen::<FieldElement>()
        .expect("distributed key generation succeeds");
    let public_key = shares[0].public_key.clone();
    let threshold = ts.get_threshold_info().threshold;
    let partials =
        generate_partials(&ts, &shares[..threshold - 1], MESSAGE, &RANDOMNESS);
    let err = ts.combine_signatures::<FieldElement>(&partials, &public_key);
    assert!(err.is_err(), "expected error for insufficient shares");
}

#[test]
fn different_participant_combinations_work() {
    let ts = ThresholdSignature::default();
    let shares = ts
        .distributed_keygen::<FieldElement>()
        .expect("distributed key generation succeeds");
    let public_key = shares[0].public_key.clone();
    let threshold = ts.get_threshold_info().threshold;
    for combo in combinations(&shares, threshold) {
        let partials = generate_partials(&ts, &combo, MESSAGE, &RANDOMNESS);
        let signature = ts
            .combine_signatures::<FieldElement>(&partials, &public_key)
            .expect("combination succeeds");
        assert_eq!(signature.c, expected_challenge(MESSAGE));
    }
}

#[test]
fn varied_messages_are_supported() {
    let ts = ThresholdSignature::default();
    let shares = ts
        .distributed_keygen::<FieldElement>()
        .expect("distributed key generation succeeds");
    let public_key = shares[0].public_key.clone();
    let messages = [
        b"short msg".as_ref(),
        b"",
        b"\x00\x01\x02\x03\x04",
        b"A much longer message to exercise the pipeline".as_ref(),
        "unicode: αβγδε".as_bytes(),
    ];

    for message in messages {
        let partials = generate_partials(
            &ts,
            &shares[..ts.get_threshold_info().threshold],
            message,
            &RANDOMNESS,
        );
        let signature = ts
            .combine_signatures::<FieldElement>(&partials, &public_key)
            .expect("combination succeeds");
        assert_eq!(
            signature.c,
            expected_challenge(message),
            "challenge mismatch for message {:?}",
            message
        );
    }
}

#[test]
fn cross_message_verification_fails() {
    let ts = ThresholdSignature::default();
    let shares = ts
        .distributed_keygen::<FieldElement>()
        .expect("distributed key generation succeeds");
    let public_key = shares[0].public_key.clone();
    let threshold = ts.get_threshold_info().threshold;

    let message_a = b"first message";
    let message_b = b"second message";

    let partials_a =
        generate_partials(&ts, &shares[..threshold], message_a, &RANDOMNESS);
    let signature = ts
        .combine_signatures::<FieldElement>(&partials_a, &public_key)
        .expect("combination succeeds");
    assert_eq!(signature.c, expected_challenge(message_a));
    assert!(
        signature.c != expected_challenge(message_b),
        "signature unexpectedly matched challenge for different message"
    );
}

#[test]
fn threshold_info_matches_configuration() {
    let threshold = 4;
    let participants = 6;
    let ts = ThresholdSignature::new(threshold, participants)
        .expect("configuration should be valid");
    let info = ts.get_threshold_info();
    assert_eq!(info.threshold, threshold);
    assert_eq!(info.participants, participants);
    assert_eq!(info.min_signers, threshold);
    assert_eq!(info.max_participants, participants);
}

#[test]
fn single_share_cannot_form_signature() {
    let ts = ThresholdSignature::default();
    let shares = ts
        .distributed_keygen::<FieldElement>()
        .expect("distributed key generation succeeds");
    let public_key = shares[0].public_key.clone();
    let single_partial =
        generate_partials(&ts, &shares[..1], MESSAGE, &RANDOMNESS);
    let err =
        ts.combine_signatures::<FieldElement>(&single_partial, &public_key);
    assert!(err.is_err(), "single share should not succeed");
}
