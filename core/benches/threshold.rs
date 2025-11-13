use criterion::{Criterion, black_box, criterion_group, criterion_main};
use dilithium_core::dilithium::threshold::ThresholdSignature;
use math::field_element::FieldElement;

const MESSAGE: &[u8] = b"benchmark threshold Dilithium workflow";

fn bench_threshold_workflow(c: &mut Criterion) {
    let ts = ThresholdSignature::default();
    let info = ts.get_threshold_info();

    c.bench_function("threshold_workflow_default", |b| {
        b.iter(|| {
            // 1. Distributed key generation
            let shares = ts
                .distributed_keygen::<FieldElement>()
                .expect("distributed key generation succeeds");
            let public_key = shares[0].public_key.clone();

            // 2. Partial signing + verification
            let mut partials = Vec::with_capacity(info.threshold);
            for share in shares.iter().take(info.threshold) {
                let partial = ts
                    .partial_sign::<FieldElement>(
                        black_box(MESSAGE),
                        share,
                        None,
                    )
                    .expect("partial signing succeeds");
                assert!(
                    ts.verify_partial_signature::<FieldElement>(
                        black_box(MESSAGE),
                        &partial,
                    ),
                    "partial signature must verify"
                );
                partials.push(partial);
            }

            // 3. Combine signatures
            ts.combine_signatures::<FieldElement>(&partials, &public_key)
                .expect("signature combination succeeds");
        });
    });
}

criterion_group!(benches, bench_threshold_workflow);
criterion_main!(benches);
