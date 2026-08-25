#!/usr/bin/env python3
"""Focused tests for the field-push threshold analyzer."""

import csv
from pathlib import Path
import tempfile

from analyze_edupic_field_push_thresholds import aggregate, relative_range


def main() -> None:
    assert abs(relative_range([{"x": 96.0}, {"x": 104.0}], "x") - 0.08) < 1e-15
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "field.csv"
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=(
                "phase_fraction", "region", "field_push_macro_observations",
                "field_push_promotions", "field_push_demotions"))
            writer.writeheader()
            writer.writerow({
                "phase_fraction": 0.25, "region": "critical",
                "field_push_macro_observations": 200,
                "field_push_promotions": 2, "field_push_demotions": 1})
            writer.writerow({
                "phase_fraction": 0.75, "region": "critical",
                "field_push_macro_observations": 999,
                "field_push_promotions": 999, "field_push_demotions": 999})
        value = aggregate(path, {"critical"}, 0.125, 0.5)
        assert value["field_push_macro_observations"] == 200
        assert value["field_push_promotions_per_million_pushes"] == 10000.0


if __name__ == "__main__":
    main()
    print("field-push threshold analyzer tests passed")
