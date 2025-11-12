SHELL := /bin/bash

.PHONY: test bench bench-basic

test:
	cargo test --all-targets

bench:
	cargo bench -p core basic

bench-basic:
	cargo bench -p core basic
