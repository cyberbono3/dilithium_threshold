SHELL := /bin/bash

.PHONY: test bench bench-basic bench-threshold examples-basic examples-threshold examples

test:
	cargo test --all-targets

bench:
	cargo bench -p core basic threshold

bench-basic:
	cargo bench -p core basic

bench-threshold:
	cargo bench -p core threshold

examples: examples-basic examples-threshold

examples-basic:
	cargo run -p dilithium-core --example basic

examples-threshold:
	cargo run -p dilithium-core --example threshold
