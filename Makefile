SHELL := /bin/bash

.PHONY: test bench

test:
	cargo test --all-targets

bench:
	cargo bench -p core sign

