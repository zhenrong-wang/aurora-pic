#!/usr/bin/env python3
"""Regressions for coincident-node mesh comparison."""

from analyze_aurorapic_mesh_refinement import coincident, phase_space_average


def main() -> None:
    refined = list(range(16 * 5))
    assert coincident(refined, 3, 5)[:6] == [0, 2, 4, 5, 7, 9]
    x = [0.0, 0.5, 1.0]
    values = [2.0, 2.0, 2.0] * 16
    assert phase_space_average(x, values) == 2.0


if __name__ == "__main__":
    main()
