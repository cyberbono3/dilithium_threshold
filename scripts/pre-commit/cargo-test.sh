#!/usr/bin/env bash
set -euo pipefail

echo "Running cargo test --workspace --all-features"
cargo test --workspace --all-features
