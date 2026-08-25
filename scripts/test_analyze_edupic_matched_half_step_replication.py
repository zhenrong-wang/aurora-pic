#!/usr/bin/env python3
"""Focused raw-count pooling regression."""

from analyze_edupic_matched_half_step_replication import pool


def main() -> None:
    result = pool([
        {"field_push_macro_observations": 100, "field_push_promotions": 2,
         "field_push_demotions": 3},
        {"field_push_macro_observations": 300, "field_push_promotions": 6,
         "field_push_demotions": 5},
    ])
    assert result["field_push_macro_observations"] == 400
    assert result["field_push_promotions"] == 8
    assert result["field_push_demotions"] == 8
    assert result["field_push_promotions_per_million_pushes"] == 20000.0
    assert result["field_push_demotions_per_million_pushes"] == 20000.0


if __name__ == "__main__":
    main()
    print("matched-half-step replication analyzer tests passed")
