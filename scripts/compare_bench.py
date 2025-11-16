#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path

LINE_RE = re.compile(
    r"^test\s+(?P<name>\S+)\s+\.\.\.\s+bench:\s+(?P<value>[0-9.]+)\s+ns/iter(?:.*)$"
)

def parse_results(path: Path) -> dict[str, float]:
    results: dict[str, float] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        match = LINE_RE.match(line)
        if not match:
            continue
        results[match.group("name")] = float(match.group("value"))
    if not results:
        raise SystemExit(f"No benchmark results parsed from {path}")
    return results


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: compare_bench.py <baseline> <current>")

    baseline_path = Path(sys.argv[1])
    current_path = Path(sys.argv[2])

    baseline = parse_results(baseline_path)
    current = parse_results(current_path)

    regressions: list[tuple[str, float]] = []
    for name, base_val in baseline.items():
        if name not in current:
            continue
        cur_val = current[name]
        change = (cur_val - base_val) / base_val * 100.0
        print(
            f"{name}: baseline={base_val:.2f}ns current={cur_val:.2f}ns change={change:.2f}%"
        )
        if change > 5.0:
            regressions.append((name, change))

    if regressions:
        details = ", ".join(f"{name} ({chg:.2f}%)" for name, chg in regressions)
        raise SystemExit(f"Benchmark regressions detected: {details}")


if __name__ == "__main__":
    main()
