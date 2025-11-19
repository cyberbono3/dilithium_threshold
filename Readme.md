Experimental implementatino of [Dilithium Threshold sginature scheme](https://repository.ubn.ru.nl/bitstream/handle/2066/191703/191703.pdf)
Work in progress

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

## Selecting Dilithium security levels

The crate compiles one ML-DSA parameter set at a time via Cargo features. Use the
`xtask` helper CLI to build, test, or bench with a specific level:

```bash
# run tests with the default level (level2)
cargo run -p xtask -- level2 test

# run benches under level3
cargo run -p xtask -- level3 bench -- --bench

# build and run the threshold example using level5
cargo run -p xtask -- level5 run --example threshold
```

The CLI forwards everything after the subcommand directly to `cargo`, so the
last example is equivalent to manually executing
`cargo run -p dilithium-core --example threshold --no-default-features --features level5`.
