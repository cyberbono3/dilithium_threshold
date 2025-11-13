SHELL := /bin/bash

.PHONY: test bench bench-basic bench-threshold

test:
	cargo test --all-targets

bench:
	cargo bench -p core basic threshold

bench-basic:
	cargo bench -p core basic

bench-threshold:
	cargo bench -p core threshold
