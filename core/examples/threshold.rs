use dilithium_core::basic::PublicKey;
use dilithium_core::dilithium::threshold::{
    PartialSignature, ThresholdKeyShare, ThresholdSignature,
};
use dilithium_core::dilithium::utils::{
    derive_challenge_polynomial, hash_message,
};
use math::{field_element::FieldElement, poly::Polynomial};

const MESSAGE: &[u8] = b"example threshold Dilithium message";
const RANDOMNESS: [u8; 32] = [0x5Au8; 32];

fn collect_partials(
    ts: &ThresholdSignature,
    shares: &[ThresholdKeyShare<'static, FieldElement>],
    randomness: &[u8],
) -> Vec<PartialSignature<'static, FieldElement>> {
    shares
        .iter()
        .map(|share| {
            ts.partial_sign::<FieldElement>(MESSAGE, share, Some(randomness))
                .expect("partial signature should succeed")
        })
        .collect()
}

fn main() {
    let ts = ThresholdSignature::default();
    let shares = ts
        .distributed_keygen::<FieldElement>()
        .expect("distributed key generation should succeed");

    let ThresholdSignatureInfo {
        public_key,
        partial_signatures,
    } = build_threshold_signature(&ts, &shares);

    println!(
        "Threshold signature created with {} partials",
        partial_signatures.len()
    );
    println!(
        "Public key matrix has dimensions {}x{}",
        public_key.a.rows(),
        public_key.a.cols()
    );
}

struct ThresholdSignatureInfo {
    public_key: PublicKey<'static, FieldElement>,
    partial_signatures: Vec<PartialSignature<'static, FieldElement>>,
}

fn build_threshold_signature(
    ts: &ThresholdSignature,
    shares: &[ThresholdKeyShare<'static, FieldElement>],
) -> ThresholdSignatureInfo {
    let threshold = ts.get_threshold_info().threshold;
    let public_key = shares[0].public_key.clone();
    let partial_signatures =
        collect_partials(ts, &shares[..threshold], &RANDOMNESS);

    let signature = ts
        .combine_signatures::<FieldElement>(&partial_signatures, &public_key)
        .expect("threshold combination should succeed");

    let expected_challenge = derive_expected_challenge(MESSAGE);
    assert_eq!(signature.c, expected_challenge);

    println!(
        "Combined signature challenge matches derived value ({} coefficients)",
        signature.c.coefficients().len()
    );
    ThresholdSignatureInfo {
        public_key,
        partial_signatures,
    }
}

fn derive_expected_challenge(
    message: &[u8],
) -> Polynomial<'static, FieldElement> {
    let mu = hash_message(message);
    let mut seed = Vec::with_capacity(mu.len() + b"challenge".len());
    seed.extend_from_slice(&mu);
    seed.extend_from_slice(b"challenge");
    derive_challenge_polynomial::<FieldElement>(&seed)
}
