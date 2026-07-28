# Particle pusher validation

AuroraPIC uses time-centered particle pushers for the 1D, 2D, and 3D simulation paths. 1D remains electrostatic leapfrog-only. Multidimensional runs use the same electrostatic leapfrog path by default and switch to a Boris rotation/kick when a uniform configured magnetic field is nonzero. This note records the state contract and the regression targets used to guard the implementation.

## State contract

- `Particle::v`, the combined `Particle2D::{velocity,velocity_z}` state, and `Particle3D::velocity` remain the time-centered velocities used by diagnostics, kinetic-energy accounting, and particle CSV output.
- `Particle::v_half`, the combined `Particle2D::{velocity_half,velocity_half_z}` state, and `Particle3D::velocity_half` store the active half-step velocity used for kicks/rotation and position drift.
- Initialization converts a time-centered velocity to the previous half step:
  - `v_half = v - 0.5 * (q / m) * E(x) * dt`
- Each electrostatic step performs:
  1. interpolate the current electric field at the particle position;
  2. kick `v_half` by `(q / m) * E * dt`;
  3. drift position with the updated `v_half`;
  4. apply particle boundaries to position and half-step velocity;
  5. redeposit charge and solve fields;
  6. synchronize the diagnostic velocity with `v = v_half + 0.5 * (q / m) * E_new * dt`.

For multidimensional reflecting particle boundaries, the normal component of the half-step velocity is reversed during boundary handling. The synchronized diagnostic velocity is then derived from the reflected half-step state after fields are resolved.

## Regression targets

The automated CTest suite covers the core pusher behavior directly:

- 1D constant acceleration matches the analytic position and velocity to roundoff.
- 2D constant acceleration matches analytic x/y positions and velocities to roundoff.
- A 1D harmonic oscillator regression keeps the maximum time-centered energy error below `1e-3` for `dt = 0.02` over 5000 steps.
- 2D3V and 3D Boris cyclotron regressions verify constant-B rotation angle, speed conservation, and arbitrary-axis equivalence between planar 2D3V and full 3D velocity updates.

The standalone validation script `scripts/validate_pushers.py` repeats the leapfrog checks without depending on the C++ test binary, adds dependency-free 2D/3D Boris rotation checks, and reports harmonic-oscillator energy drift over a small timestep sweep. It is included in `scripts/verify.sh` so smoke verification fails if the baseline pusher contract regresses.

## Current scope

Boris support currently uses prescribed uniform `magnetic_field_x/y/z` in both planar 2D3V and 3D while retaining the electrostatic Poisson field solve. In 2D3V, particle position and electric field remain planar, but all three velocity components participate in magnetic rotation. Self-consistent electromagnetic field updates are still a later roadmap item.
