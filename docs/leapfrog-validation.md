# Leapfrog pusher validation

AuroraPIC uses an electrostatic, time-centered leapfrog particle pusher for the 1D and 2D simulation paths. This note records the state contract and the regression targets used to guard the implementation.

## State contract

- `Particle::v` and `Particle2D::velocity` remain the time-centered velocities used by diagnostics, kinetic-energy accounting, and particle CSV output.
- `Particle::v_half` and `Particle2D::velocity_half` store the active leapfrog half-step velocity used for kicks and position drift.
- Initialization converts a time-centered velocity to the previous half step:
  - `v_half = v - 0.5 * (q / m) * E(x) * dt`
- Each electrostatic step performs:
  1. interpolate the current electric field at the particle position;
  2. kick `v_half` by `(q / m) * E * dt`;
  3. drift position with the updated `v_half`;
  4. apply particle boundaries to position and half-step velocity;
  5. redeposit charge and solve fields;
  6. synchronize the diagnostic velocity with `v = v_half + 0.5 * (q / m) * E_new * dt`.

For 2D reflecting particle boundaries, the normal component of `velocity_half` is reversed during boundary handling. The synchronized diagnostic velocity is then derived from the reflected half-step state after fields are resolved.

## Regression targets

The automated CTest suite covers the core pusher behavior directly:

- 1D constant acceleration matches the analytic position and velocity to roundoff.
- 2D constant acceleration matches analytic x/y positions and velocities to roundoff.
- A 1D harmonic oscillator regression keeps the maximum time-centered energy error below `1e-3` for `dt = 0.02` over 5000 steps.

The standalone validation script `scripts/validate_leapfrog.py` repeats these checks without depending on the C++ test binary and also reports harmonic-oscillator energy drift over a small timestep sweep. It is included in `scripts/verify.sh` so smoke verification fails if the pusher contract regresses.

## Current scope

The pusher is electrostatic only. Magnetic-field rotation/Boris integration is intentionally not implemented yet; adding it should preserve the time-centered diagnostic contract or explicitly document any state-model change.
