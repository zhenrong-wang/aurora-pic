#!/usr/bin/env python3
"""Regressions for same-grid particle-refinement comparison."""

from analyze_aurorapic_particle_refinement import require_same_grid


def main() -> None:
    require_same_grid([0.0, 0.5, 1.0], [0.0, 0.5, 1.0])
    try:
        require_same_grid([0.0, 1.0], [0.0, 0.9])
    except ValueError:
        pass
    else:
        raise AssertionError("different grids were accepted")


if __name__ == "__main__":
    main()
