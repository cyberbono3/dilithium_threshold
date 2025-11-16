use dilithium_core::basic::PublicKey;
use dilithium_core::dilithium::threshold::{
    PartialSignature, ThresholdKeyShare, ThresholdSignature,
};
use math::field_element::FieldElement;

const MESSAGE: &[u8] = b"example threshold Dilithium message";

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
    let randomness = [0x5Au8; 32];
    let partial_signatures =
        collect_partials(ts, &shares[..threshold], &randomness);

    let signature = ts
        .combine_signatures::<FieldElement>(&partial_signatures, &public_key)
        .expect("threshold combination should succeed");

    assert!(
        public_key.verify(MESSAGE, &signature),
        "combined signature should verify with the public key"
    );

    ThresholdSignatureInfo {
        public_key,
        partial_signatures,
    }
}
