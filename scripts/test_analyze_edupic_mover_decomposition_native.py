#!/usr/bin/env python3
"""Focused tests for native mover-decomposition aggregation."""

from analyze_edupic_mover_decomposition_native import relative_range


def main() -> None:
    members = [{"metric": 0.95}, {"metric": 1.0}, {"metric": 1.05}]
    assert abs(relative_range(members, "metric") - 0.1) < 1e-15
    assert relative_range([{"metric": 0.0}] * 3, "metric") == 0.0


if __name__ == "__main__":
    main()
    print("native mover-decomposition analyzer tests passed")
