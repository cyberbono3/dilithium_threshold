#!/usr/bin/env bash
set -euo pipefail

echo "Running cargo fmt --all -- --check"
cargo fmt --all -- --check
