#!/usr/bin/env python3
"""Focused unit tests for field-push promotion localization helpers."""

from __future__ import annotations

from localize_edupic_field_push_promotions import aggregate, relative_range


def main() -> int:
    rows = [
        {
            "phase_fraction": "0.2", "region_id": "2",
            "field_push_macro_observations": "100",
            "field_push_promotions": "2", "field_push_demotions": "1",
        },
        {
            "phase_fraction": "0.3", "region_id": "2",
            "field_push_macro_observations": "300",
            "field_push_promotions": "6", "field_push_demotions": "3",
        },
        {
            "phase_fraction": "0.2", "region_id": "1",
            "field_push_macro_observations": "999",
            "field_push_promotions": "999", "field_push_demotions": "999",
        },
    ]
    value = aggregate(rows, (2,), 0.125, 0.375)
    assert value["rows"] == 2
    assert value["field_push_macro_observations"] == 400
    assert value["field_push_promotions_per_million_pushes"] == 20000.0
    assert value["field_push_demotions_per_million_pushes"] == 10000.0
    assert abs(relative_range([8.0, 12.0]) - 0.4) < 1.0e-12
    print("field-push promotion localization tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
