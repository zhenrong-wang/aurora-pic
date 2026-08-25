#!/usr/bin/env python3
"""Focused tests for threshold-crossing replication calculations."""

from analyze_edupic_threshold_crossing_replication import relative_range


def main() -> None:
    members = [{"rate": 96.0}, {"rate": 104.0}]
    assert abs(relative_range(members, "rate") - 0.08) < 1.0e-15
    assert relative_range([{"rate": 0.0}, {"rate": 0.0}], "rate") == 0.0


if __name__ == "__main__":
    main()
    print("threshold-crossing replication analyzer tests passed")
