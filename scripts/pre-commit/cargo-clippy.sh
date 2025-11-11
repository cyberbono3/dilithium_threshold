#!/usr/bin/env bash
set -euo pipefail

echo "Running cargo clippy --workspace --all-targets --all-features -- -D warnings"
cargo clippy --workspace --all-targets --all-features -- -D warnings
