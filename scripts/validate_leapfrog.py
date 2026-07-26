#!/usr/bin/env python3
"""Standalone numerical checks for AuroraPIC's electrostatic leapfrog pusher.

This script mirrors the public pusher state contract without depending on the C++
test binary. It is intentionally dependency-free so it can run as part of the
smoke verification script on a minimal developer machine.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Iterable, Tuple


@dataclass
class Particle1D:
    x: float
    v: float
    v_half: float = 0.0


@dataclass
class Vec2:
    x: float
    y: float


@dataclass
class Particle2D:
    position: Vec2
    velocity: Vec2
    velocity_half: Vec2


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def require_near(actual: float, expected: float, tolerance: float, message: str) -> None:
    if abs(actual - expected) > tolerance:
        raise AssertionError(f"{message}: actual={actual:.17g} expected={expected:.17g} tolerance={tolerance:g}")


def initialize_half_step_1d(particle: Particle1D, electric: float, charge_to_mass: float, dt: float) -> None:
    particle.v_half = particle.v - 0.5 * charge_to_mass * electric * dt


def kick_1d(particle: Particle1D, electric: float, charge_to_mass: float, dt: float) -> None:
    particle.v_half += charge_to_mass * electric * dt


def drift_1d(particle: Particle1D, dt: float) -> None:
    particle.x += particle.v_half * dt


def synchronize_1d(particle: Particle1D, electric: float, charge_to_mass: float, dt: float) -> None:
    particle.v = particle.v_half + 0.5 * charge_to_mass * electric * dt


def initialize_half_step_2d(particle: Particle2D, electric: Vec2, charge_to_mass: float, dt: float) -> None:
    particle.velocity_half.x = particle.velocity.x - 0.5 * charge_to_mass * electric.x * dt
    particle.velocity_half.y = particle.velocity.y - 0.5 * charge_to_mass * electric.y * dt


def kick_2d(particle: Particle2D, electric: Vec2, charge_to_mass: float, dt: float) -> None:
    particle.velocity_half.x += charge_to_mass * electric.x * dt
    particle.velocity_half.y += charge_to_mass * electric.y * dt


def drift_2d(particle: Particle2D, dt: float) -> None:
    particle.position.x += particle.velocity_half.x * dt
    particle.position.y += particle.velocity_half.y * dt


def synchronize_2d(particle: Particle2D, electric: Vec2, charge_to_mass: float, dt: float) -> None:
    particle.velocity.x = particle.velocity_half.x + 0.5 * charge_to_mass * electric.x * dt
    particle.velocity.y = particle.velocity_half.y + 0.5 * charge_to_mass * electric.y * dt


def validate_1d_constant_acceleration() -> None:
    dt = 0.125
    acceleration = 2.5
    particle = Particle1D(x=0.2, v=-0.4)
    initialize_half_step_1d(particle, acceleration, 1.0, dt)
    for _ in range(8):
        kick_1d(particle, acceleration, 1.0, dt)
        drift_1d(particle, dt)
        synchronize_1d(particle, acceleration, 1.0, dt)
    t = 8.0 * dt
    require_near(particle.x, 0.2 + (-0.4) * t + 0.5 * acceleration * t * t, 1e-14,
                 "1D constant-acceleration position regression failed")
    require_near(particle.v, -0.4 + acceleration * t, 1e-14,
                 "1D constant-acceleration velocity regression failed")


def validate_2d_constant_acceleration() -> None:
    dt = 0.1
    acceleration = Vec2(1.5, -0.75)
    particle = Particle2D(position=Vec2(0.25, 0.75), velocity=Vec2(-0.2, 0.3), velocity_half=Vec2(0.0, 0.0))
    initialize_half_step_2d(particle, acceleration, 1.0, dt)
    for _ in range(6):
        kick_2d(particle, acceleration, 1.0, dt)
        drift_2d(particle, dt)
        synchronize_2d(particle, acceleration, 1.0, dt)
    t = 6.0 * dt
    require_near(particle.position.x, 0.25 + (-0.2) * t + 0.5 * acceleration.x * t * t, 1e-14,
                 "2D constant-acceleration x-position regression failed")
    require_near(particle.position.y, 0.75 + 0.3 * t + 0.5 * acceleration.y * t * t, 1e-14,
                 "2D constant-acceleration y-position regression failed")
    require_near(particle.velocity.x, -0.2 + acceleration.x * t, 1e-14,
                 "2D constant-acceleration x-velocity regression failed")
    require_near(particle.velocity.y, 0.3 + acceleration.y * t, 1e-14,
                 "2D constant-acceleration y-velocity regression failed")


def harmonic_oscillator_max_energy_error(dt: float, steps: int) -> float:
    particle = Particle1D(x=1.0, v=0.0)
    initialize_half_step_1d(particle, -particle.x, 1.0, dt)
    initial_energy = 0.5 * (particle.x * particle.x + particle.v * particle.v)
    max_energy_error = 0.0
    for _ in range(steps):
        kick_1d(particle, -particle.x, 1.0, dt)
        drift_1d(particle, dt)
        synchronize_1d(particle, -particle.x, 1.0, dt)
        energy = 0.5 * (particle.x * particle.x + particle.v * particle.v)
        max_energy_error = max(max_energy_error, abs(energy - initial_energy))
    require(isfinite(max_energy_error), "harmonic oscillator energy error is not finite")
    return max_energy_error


def validate_harmonic_oscillator() -> Tuple[Tuple[float, float], ...]:
    primary_error = harmonic_oscillator_max_energy_error(dt=0.02, steps=5000)
    require(primary_error < 1e-3,
            f"1D harmonic oscillator energy regression failed: max_error={primary_error:.17g}")

    sweep = tuple((dt, harmonic_oscillator_max_energy_error(dt=dt, steps=round(100.0 / dt)))
                  for dt in (0.04, 0.02, 0.01))
    # The time-centered diagnostic energy should improve approximately quadratically
    # as dt decreases for this smooth oscillator. Use loose monotonic checks to avoid
    # overfitting the exact floating-point trajectory.
    for (coarse_dt, coarse_error), (fine_dt, fine_error) in zip(sweep, sweep[1:]):
        require(fine_error < coarse_error,
                f"harmonic oscillator energy did not improve from dt={coarse_dt} to dt={fine_dt}")
    return sweep


def format_sweep(sweep: Iterable[Tuple[float, float]]) -> str:
    return ", ".join(f"dt={dt:g}:max_energy_error={error:.12g}" for dt, error in sweep)


def main() -> int:
    validate_1d_constant_acceleration()
    validate_2d_constant_acceleration()
    sweep = validate_harmonic_oscillator()
    print("leapfrog validation passed")
    print(format_sweep(sweep))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
