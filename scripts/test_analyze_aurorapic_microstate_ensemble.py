#!/usr/bin/env python3
"""Regression tests for constrained-microstate ensemble helpers."""

import math

from analyze_aurorapic_microstate_ensemble import profile_scatter


def main() -> None:
    members = [
        {"value": [1.0] * (200 * 400)},
        {"value": [1.1] * (200 * 400)},
        {"value": [0.9] * (200 * 400)},
    ]
    ensemble, spatial, scatter = profile_scatter(members, "value")
    assert len(ensemble) == 200 * 400
    assert len(spatial) == 400
    assert all(math.isclose(value, 1.0) for value in ensemble)
    assert math.isclose(scatter[0], 0.0, abs_tol=1.0e-15)
    assert math.isclose(scatter[1], 0.1)
    assert math.isclose(scatter[2], 0.1)


if __name__ == "__main__":
    main()
