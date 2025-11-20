Experimental implementatino of [Dilithium Threshold sginature scheme](https://repository.ubn.ru.nl/bitstream/handle/2066/191703/191703.pdf)
Work in progress

[![Tests](https://github.com/cyberbono3/dilithium_threshold/actions/workflows/tests.yml/badge.svg)](https://github.com/cyberbono3/dilithium_threshold/actions/workflows/tests.yml)
[![Coverage](https://github.com/cyberbono3/dilithium_threshold/actions/workflows/coverage.yml/badge.svg)](https://github.com/cyberbono3/dilithium_threshold/actions/workflows/coverage.yml)
[![Benchmarks](https://github.com/cyberbono3/dilithium_threshold/actions/workflows/benchmarks.yml/badge.svg)](https://github.com/cyberbono3/dilithium_threshold/actions/workflows/benchmarks.yml)




## Examples

Demonstrate the basic Dilithium API:

```bash
cargo run -p dilithium-core --example basic
# or via the Makefile target
make examples-basic
```

Demonstrate the threshold signing workflow:

```bash
cargo run -p dilithium-core --example threshold
# or via the Makefile target
make examples-threshold
```

Run both examples together:

```bash
make examples
```
