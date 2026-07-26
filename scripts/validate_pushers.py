#!/usr/bin/env python3
"""Standalone numerical checks for AuroraPIC's particle pushers.

This script mirrors the public pusher state contract without depending on the C++
test binary. It is intentionally dependency-free so it can run as part of the
smoke verification script on a minimal developer machine.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import atan, cos, isfinite, sin, sqrt
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
class Vec3:
    x: float
    y: float
    z: float


@dataclass
class Particle2D:
    position: Vec2
    velocity: Vec2
    velocity_half: Vec2


@dataclass
class Particle3D:
    position: Vec3
    velocity: Vec3
    velocity_half: Vec3


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def require_near(actual: float, expected: float, tolerance: float, message: str) -> None:
    if abs(actual - expected) > tolerance:
        raise AssertionError(f"{message}: actual={actual:.17g} expected={expected:.17g} tolerance={tolerance:g}")


def add3(a: Vec3, b: Vec3) -> Vec3:
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z)


def scale3(v: Vec3, factor: float) -> Vec3:
    return Vec3(factor * v.x, factor * v.y, factor * v.z)


def dot3(a: Vec3, b: Vec3) -> float:
    return a.x * b.x + a.y * b.y + a.z * b.z


def norm3(v: Vec3) -> float:
    return sqrt(dot3(v, v))


def cross3(a: Vec3, b: Vec3) -> Vec3:
    return Vec3(a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x)


def to_vec3(v: Vec2) -> Vec3:
    return Vec3(v.x, v.y, 0.0)


def to_vec2(v: Vec3) -> Vec2:
    return Vec2(v.x, v.y)


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


def boris_advance(velocity: Vec3, electric: Vec3, magnetic: Vec3, charge_to_mass: float, dt: float) -> Vec3:
    half_qm_dt = 0.5 * charge_to_mass * dt
    v_minus = add3(velocity, scale3(electric, half_qm_dt))
    t = scale3(magnetic, half_qm_dt)
    t2 = dot3(t, t)
    s = scale3(t, 2.0 / (1.0 + t2))
    v_prime = add3(v_minus, cross3(v_minus, t))
    v_plus = add3(v_minus, cross3(v_prime, s))
    return add3(v_plus, scale3(electric, half_qm_dt))


def initialize_boris_half_step_2d(particle: Particle2D, electric: Vec2, magnetic_z: float,
                                  charge_to_mass: float, dt: float) -> None:
    particle.velocity_half = to_vec2(boris_advance(to_vec3(particle.velocity), to_vec3(electric),
                                                  Vec3(0.0, 0.0, magnetic_z), charge_to_mass, -0.5 * dt))


def kick_boris_2d(particle: Particle2D, electric: Vec2, magnetic_z: float, charge_to_mass: float, dt: float) -> None:
    particle.velocity_half = to_vec2(boris_advance(to_vec3(particle.velocity_half), to_vec3(electric),
                                                  Vec3(0.0, 0.0, magnetic_z), charge_to_mass, dt))


def synchronize_boris_2d(particle: Particle2D, electric: Vec2, magnetic_z: float,
                         charge_to_mass: float, dt: float) -> None:
    particle.velocity = to_vec2(boris_advance(to_vec3(particle.velocity_half), to_vec3(electric),
                                             Vec3(0.0, 0.0, magnetic_z), charge_to_mass, 0.5 * dt))


def initialize_boris_half_step_3d(particle: Particle3D, electric: Vec3, magnetic: Vec3,
                                  charge_to_mass: float, dt: float) -> None:
    particle.velocity_half = boris_advance(particle.velocity, electric, magnetic, charge_to_mass, -0.5 * dt)


def kick_boris_3d(particle: Particle3D, electric: Vec3, magnetic: Vec3, charge_to_mass: float, dt: float) -> None:
    particle.velocity_half = boris_advance(particle.velocity_half, electric, magnetic, charge_to_mass, dt)


def synchronize_boris_3d(particle: Particle3D, electric: Vec3, magnetic: Vec3,
                         charge_to_mass: float, dt: float) -> None:
    particle.velocity = boris_advance(particle.velocity_half, electric, magnetic, charge_to_mass, 0.5 * dt)


def boris_rotation_angle(magnetic_magnitude: float, charge_to_mass: float, dt: float) -> float:
    return -2.0 * atan(0.5 * charge_to_mass * magnetic_magnitude * dt)


def rotate_about_axis(v: Vec3, axis: Vec3, angle: float) -> Vec3:
    c = cos(angle)
    s = sin(angle)
    axis_cross_v = cross3(axis, v)
    axis_dot_v = dot3(axis, v)
    return add3(add3(scale3(v, c), scale3(axis_cross_v, s)), scale3(axis, axis_dot_v * (1.0 - c)))


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


def validate_2d_boris_cyclotron() -> None:
    dt = 0.04
    steps = 37
    charge_to_mass = 1.25
    magnetic_z = 1.7
    initial_velocity = Vec2(0.8, -0.35)
    particle = Particle2D(position=Vec2(0.0, 0.0), velocity=Vec2(initial_velocity.x, initial_velocity.y),
                          velocity_half=Vec2(0.0, 0.0))
    electric = Vec2(0.0, 0.0)
    initialize_boris_half_step_2d(particle, electric, magnetic_z, charge_to_mass, dt)
    for _ in range(steps):
        kick_boris_2d(particle, electric, magnetic_z, charge_to_mass, dt)
    synchronize_boris_2d(particle, electric, magnetic_z, charge_to_mass, dt)

    angle = steps * boris_rotation_angle(abs(magnetic_z), charge_to_mass, dt)
    expected_x = initial_velocity.x * cos(angle) - initial_velocity.y * sin(angle)
    expected_y = initial_velocity.x * sin(angle) + initial_velocity.y * cos(angle)
    require_near(particle.velocity.x, expected_x, 1e-13, "2D Boris cyclotron x-velocity regression failed")
    require_near(particle.velocity.y, expected_y, 1e-13, "2D Boris cyclotron y-velocity regression failed")
    initial_speed = sqrt(initial_velocity.x * initial_velocity.x + initial_velocity.y * initial_velocity.y)
    final_speed = sqrt(particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y)
    require_near(final_speed, initial_speed, 1e-13, "2D Boris cyclotron speed conservation failed")


def validate_3d_boris_arbitrary_axis() -> None:
    dt = 0.03
    steps = 29
    charge_to_mass = -0.75
    magnetic = Vec3(0.4, -0.8, 1.1)
    initial_velocity = Vec3(0.6, -0.25, 0.9)
    particle = Particle3D(position=Vec3(0.0, 0.0, 0.0), velocity=Vec3(initial_velocity.x, initial_velocity.y, initial_velocity.z),
                          velocity_half=Vec3(0.0, 0.0, 0.0))
    electric = Vec3(0.0, 0.0, 0.0)
    initialize_boris_half_step_3d(particle, electric, magnetic, charge_to_mass, dt)
    for _ in range(steps):
        kick_boris_3d(particle, electric, magnetic, charge_to_mass, dt)
    synchronize_boris_3d(particle, electric, magnetic, charge_to_mass, dt)

    magnetic_magnitude = norm3(magnetic)
    axis = scale3(magnetic, 1.0 / magnetic_magnitude)
    angle = steps * boris_rotation_angle(magnetic_magnitude, charge_to_mass, dt)
    expected = rotate_about_axis(initial_velocity, axis, angle)
    require_near(particle.velocity.x, expected.x, 1e-13, "3D Boris arbitrary-axis x-velocity regression failed")
    require_near(particle.velocity.y, expected.y, 1e-13, "3D Boris arbitrary-axis y-velocity regression failed")
    require_near(particle.velocity.z, expected.z, 1e-13, "3D Boris arbitrary-axis z-velocity regression failed")
    require_near(norm3(particle.velocity), norm3(initial_velocity), 1e-13,
                 "3D Boris arbitrary-axis speed conservation failed")
    require_near(dot3(particle.velocity, axis), dot3(initial_velocity, axis), 1e-13,
                 "3D Boris arbitrary-axis parallel velocity conservation failed")


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
    validate_2d_boris_cyclotron()
    validate_3d_boris_arbitrary_axis()
    sweep = validate_harmonic_oscillator()
    print("pusher validation passed")
    print(format_sweep(sweep))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
