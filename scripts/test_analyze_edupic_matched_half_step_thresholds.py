#!/usr/bin/env python3
"""Focused decision tests for matched-half-step threshold analysis."""

from analyze_edupic_matched_half_step_thresholds import relative_range


def main() -> None:
    metric = "field_push_promotions_per_million_pushes"
    passing = [{metric: 8.0}, {metric: 8.64}]
    failing = [{metric: 8.0}, {metric: 8.7}]
    assert relative_range(passing, metric) <= 0.08
    assert relative_range(failing, metric) > 0.08


if __name__ == "__main__":
    main()
    print("matched-half-step threshold analyzer tests passed")
