use core::dilithium::threshold::ThresholdSignature;
use criterion::{Criterion, black_box, criterion_group, criterion_main};
use math::field_element::FieldElement;

const MESSAGE: &[u8] = b"benchmark threshold Dilithium workflow";
const CONFIGURATIONS: &[(usize, usize)] = &[(2, 3), (3, 5), (4, 7)];

fn bench_threshold_workflow(c: &mut Criterion) {
    let mut group = c.benchmark_group("threshold_workflow");

    for &(threshold, participants) in CONFIGURATIONS {
        let ts = ThresholdSignature::new(threshold, participants)
            .expect("valid threshold configuration");
        group.bench_function(format!("{threshold}-of-{participants}"), |b| {
            b.iter(|| {
                // 1. Distributed key generation
                let shares = ts
                    .distributed_keygen::<FieldElement>()
                    .expect("distributed key generation succeeds");
                assert_eq!(shares.len(), participants);

                let public_key = shares[0].public_key.clone();
                for share in &shares {
                    assert_eq!(
                        share.public_key.a.shape(),
                        public_key.a.shape(),
                        "matrix shapes must match"
                    );
                    assert_eq!(
                        share.public_key.t, public_key.t,
                        "public vectors must match"
                    );
                }

                // 2. Partial signing + verification
                let mut partials = Vec::with_capacity(threshold);
                for share in shares.iter().take(threshold) {
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

    group.finish();
}

criterion_group!(benches, bench_threshold_workflow);
criterion_main!(benches);
