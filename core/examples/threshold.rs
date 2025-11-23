use dilithium_core::dilithium::threshold::ThresholdSignature;
use math::field_element::FieldElement;

const MESSAGE: &[u8] = b"Original message";
const WRONG_MESSAGE: &[u8] = b"Wrong message";
const RANDOMNESS: [u8; 32] = [0x5A; 32];

/// Demonstrate partial signing and verification for a threshold instance.
/// The partial signature should verify the intended message and fail for a
/// mismatched message.
fn main() {
    let ts = ThresholdSignature::default();
    let shares = ts
        .distributed_keygen::<FieldElement>()
        .expect("distributed key generation should succeed");

    let partial_sig = ts
        .partial_sign::<FieldElement>(MESSAGE, &shares[0], Some(&RANDOMNESS))
        .expect("partial signature should succeed");

    let is_valid = ts.verify_partial_signature(MESSAGE, &partial_sig);
    let is_wrong_message_valid =
        ts.verify_partial_signature(WRONG_MESSAGE, &partial_sig);

    assert!(
        is_valid,
        "expected partial signature to verify for its message"
    );
    assert!(
        !is_wrong_message_valid,
        "partial signature must not verify a different message"
    );

    println!("Partial signature verified successfully.");
    println!(
        "Cross-message verification rejected as expected: {}",
        is_wrong_message_valid
    );
}
