
# dilithium-threshold 

[![Tests](https://github.com/cyberbono3/dilithium_threshold/actions/workflows/tests.yml/badge.svg)](https://github.com/cyberbono3/dilithium_threshold/actions/workflows/tests.yml)
[![Benchmarks](https://github.com/cyberbono3/dilithium_threshold/actions/workflows/benchmarks.yml/badge.svg)](https://github.com/cyberbono3/dilithium_threshold/actions/workflows/benchmarks.yml)

Unofficial rust implementation of [Dilithium post-quantum digital signature scheme](https://repository.ubn.ru.nl/bitstream/handle/2066/191703/191703.pdf)

## Disclamer: Do not use it in production. This is the experimental implementation.

Dilithium was submitted to [NIST's Post-Quantum Cryptography standardization process](https://csrc.nist.gov/projects/post-quantum-cryptography/post-quantum-cryptography-standardization) and was selected as one of the algorithms to be standardized. It is one of the digital signature algorithms chosen by NIST, alongside FALCON and SPHINCS+. 

## Project structure
- `core/`: Dilithium implementation with basic and threshold signing APIs, examples, benches.
- `math/`: Field, polynomial, NTT, and matrix utilities used by the core crate.
- `scripts/`: Helper scripts (e.g., benchmark comparison).
- `.github/workflows/`: CI for tests, benchmarks, and coverage.
- `Makefile`: Convenience targets for examples, benches, and tests.
- `clippy.toml`: Lint configuration shared in CI and local runs.

## Library contents
- `basic` (in `core/src/basic`): keygen/sign/verify and signature serialization helpers.
- `threshold` (in `core/src/dilithium/threshold`): distributed keygen, partial signing, combination, and verification utilities.
- `dilithium/utils` (in `core/src/dilithium/utils.rs`): randomness, hashing, and polynomial helpers used across algorithms.
- `math` crate: finite field arithmetic, polynomials, NTT, matrix ops, and supporting traits.
- Examples: `core/examples/basic.rs` (single signer) and `core/examples/threshold.rs` (partial signing and verification).

## Security levels
- Default: Level 2 (ML-DSA-44 / Dilithium-2) via `DEFAULT_CONFIG` in `core/src/dilithium/params.rs`; exported constants like `K`, `L`, `BETA`, `GAMMA1`, `GAMMA2`, and `ETA` are derived from it and drive keygen, signing, and verification.
- Available sets: Levels 2, 3, and 5 are defined as `DILITHIUM_L{2,3,5}_CONFIG` and can be selected with `SecurityLevel` or `DilithiumConfig::new(level)` (returns an error if an unsupported level is requested).
- Tuning code: use `DilithiumConfig::for_level(SecurityLevel::Level3)` (or `.new(3)`) to retrieve parameters; consumers must thread the chosen config through any code that should differ from the defaults.
- Tests guard that `DEFAULT_CONFIG` equals the Level 2 constants; `DEFAULT_SECURITY_LEVEL` is kept for compatibility but not consumed outside tests.

## Examples

### 1. Basic scanario

```rust
let message = b"example Dilithium signing message";
let keypair = dilithium_core::basic::keygen::<FieldElement>().expect("key generation should succeed");
let signature = keypair.sign(message).expect("signing should succeed");
assert!(keypair.verify(message, &signature));
```

Run:
```bash
make examples-basic
```

### 2. Threshold scenario
```rust
 let threshold_sig = dilithium_core::threshold::sign::ThresholdSignature::default();
let shares = threshold_sig.distributed_keygen::<FieldElement>().unwrap();
let message = b"verified_test";
let partial_sig = threshold_sig.partial_sign(message, &shares[0], Some(&create_test_seed(1))).expect("Partial sign has to succeed");
assert!(threshold_sig.verify_partial_signature(message, &partial_sig));
```
Run:
```bash
make examples-threshold
```

## Benchmarks

### 1. Basic scenario
```bash
make bench-basic
```
### 2. Threshold scenario
```bash
make bench-threshold
```

## Testing
```bash
make test
```

## License 
Some code in `math` crate has been adopted from [twenty-first](https://github.com/Neptune-Crypto/twenty-first)  library under GPL-2.0 license

