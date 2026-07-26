# AuroraPIC methodology

AuroraPIC v0.1 implements a deliberately bounded, auditable electrostatic `1D1V` Particle-in-Cell model for first-version plasma dynamics research and engineering studies.

## Model

Particles carry position `x`, velocity `vx`, charge, mass, and macro-particle weight. Charge is deposited to a one-dimensional mesh with cloud-in-cell weighting. In 1D CLI configuration, `weight` is the direct macro-particle weight; alternatively, if `weight` is omitted, `density` is converted to `weight = density * initialization_width / particles` so the represented line density is explicit. The electrostatic field is obtained from Poisson's equation,

```text
d²phi/dx² = -rho / eps0,      E = -dphi/dx
```

in normalized units with `eps0 = 1`.

## Field solvers

- Periodic domains use a direct spectral Poisson solve. The zero mode is removed, enforcing global quasi-neutral compatibility.
- Dirichlet domains use a tridiagonal finite-difference Poisson solve with prescribed endpoint potentials.

## Particle advance and steady state

AuroraPIC advances particles with a time-centered electrostatic leapfrog update. Each step deposits charge at particle positions, solves Poisson's equation for `E^n`, kicks stored half-step velocities, drifts positions to the next time level, reapplies particle boundaries, redeposits charge, resolves fields, and synchronizes the public velocity fields used by diagnostics/output. In 2D, reflecting particle boundaries reverse the normal half-step velocity so subsequent diagnostics remain consistent after synchronization. Steady-state mode is an engineering stop condition based on the relative change between adjacent total-energy diagnostic windows, not a proof of physical equilibrium; use conservative tolerances/windows and inspect the emitted diagnostics for oscillatory systems.

## Collisions

Optional Monte-Carlo BGK-like velocity reset collisions model scattering against a prescribed neutral bath. The collision probability is `1 - exp(-nu dt)`.

## Stability guidance

For credible runs, choose:

- cell size small enough to resolve the Debye length;
- time step below electron plasma and transit-time scales;
- enough particles per cell to control PIC noise;
- output and steady-state windows long enough to distinguish relaxation from oscillation.

## Verification included

The automated test suite checks the periodic spectral Poisson solve against an analytic sinusoidal charge distribution and runs a short neutral two-species PIC simulation end-to-end.

## Extension path

AuroraPIC should extend from 1D to 2D first, then 3D, while keeping the current 1D model as a regression and validation target. The multidimensional plan is documented in `docs/multidimensional-roadmap.md`.
