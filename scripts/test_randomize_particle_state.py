#!/usr/bin/env python3
"""Regression tests for constrained particle-state randomization."""

import math

from randomize_particle_state import (
    cell_velocity_fingerprint, nodal_number, randomize_state,
)


def records(offset: float) -> list[tuple[float, ...]]:
    return [
        (0.001 + offset, 0.0, 0.0, 1.0, 2.0, 3.0),
        (0.004 + offset, 0.0, 0.0, 4.0, 5.0, 6.0),
        (0.006 + offset, 0.0, 0.0, 7.0, 8.0, 9.0),
        (0.011 + offset, 0.0, 0.0, 10.0, 11.0, 12.0),
        (0.014 + offset, 0.0, 0.0, 13.0, 14.0, 15.0),
    ]


def main() -> None:
    source = {"electrons": records(0.0), "ions": records(0.0002)}
    first, report = randomize_state(source, 0.025, 6, 51949)
    repeated, _ = randomize_state(source, 0.025, 6, 51949)
    second, _ = randomize_state(source, 0.025, 6, 63059)
    assert first == repeated
    assert first != second
    assert first != source
    for name in source:
        assert len(first[name]) == len(source[name])
        assert cell_velocity_fingerprint(first[name], 0.025, 6) == \
            cell_velocity_fingerprint(source[name], 0.025, 6)
        before = nodal_number(source[name], 0.025, 6)
        after = nodal_number(first[name], 0.025, 6)
        assert all(math.isclose(a, b, abs_tol=1.0e-14)
                   for a, b in zip(before, after))
        assert report[name]["nodal_number_relative_l1_error"] <= 1.0e-14
        assert report[name]["cellwise_velocity_tuple_multisets_preserved"] is True
    try:
        randomize_state(source, 0.025, 6, -1)
    except Exception:
        pass
    else:
        raise AssertionError("negative seed must fail")


if __name__ == "__main__":
    main()
