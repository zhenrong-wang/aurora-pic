#!/usr/bin/env python3
"""Focused tests for discrete Poisson source attribution."""

from analyze_edupic_poisson_source_attribution import (
    add_components, recover_field, shapley_attribution, solve_dirichlet,
)


def main() -> None:
    rho = [0.0, 0.0, 0.0, 0.0, 0.0]
    potential = solve_dirichlet(rho, 4.0, 0.0, 0.25, 2.0)
    assert all(abs(value - expected) < 1e-14 for value, expected in zip(
        potential, [4.0, 3.0, 2.0, 1.0, 0.0], strict=True))
    field = recover_field(potential, rho, 0.25, 2.0)
    assert all(abs(value - 4.0) < 1e-14 for value in field)

    def matrix(value: float) -> list[list[float]]:
        return [[value]]

    aurora = {factor: matrix(value) for factor, value in zip(
        ("boundary_drive", "ion_space_charge", "electron_space_charge"),
        (1.0, 2.0, 3.0), strict=True)}
    native = {factor: matrix(value) for factor, value in zip(
        ("boundary_drive", "ion_space_charge", "electron_space_charge"),
        (2.0, 4.0, 6.0), strict=True)}
    allocation, start, end, residual = shapley_attribution(
        aurora, native, lambda parts: add_components(parts)[0][0] ** 2)
    assert start == 36.0 and end == 144.0
    assert abs(sum(allocation.values()) - 108.0) < 1e-12
    assert abs(residual) < 1e-12
    print("Poisson source-attribution analyzer tests passed")


if __name__ == "__main__":
    main()
