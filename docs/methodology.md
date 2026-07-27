# AuroraPIC methodology

AuroraPIC v0.1 implements deliberately bounded, auditable electrostatic `1D1V`, structured `2D2V`, and structured `3D3V` Particle-in-Cell models for plasma dynamics research and engineering studies.

## Model

Particles carry position `x`, velocity `vx`, charge, mass, and macro-particle weight. Charge is deposited to a one-dimensional mesh with cloud-in-cell weighting. In 1D CLI configuration, `weight` is the direct macro-particle weight; alternatively, if `weight` is omitted, `density` is converted to `weight = density * initialization_width / particles` so the represented line density is explicit. The electrostatic field is obtained from Poisson's equation,

```text
d²phi/dx² = -rho / eps0,      E = -dphi/dx
```

in normalized units with `eps0 = 1`.

## Field solvers

- Periodic domains use a direct spectral Poisson solve. The zero mode is removed, enforcing global quasi-neutral compatibility.
- Dirichlet domains use a tridiagonal finite-difference Poisson solve with prescribed endpoint potentials.
- Imported finite-element domains support label-wise constant Dirichlet values and Neumann outward normal derivatives. For `-laplacian(phi) = rho / epsilon_0`, the Neumann term `dphi/dn` enters the weak-form right-hand side with a positive sign; therefore `E dot n = -dphi/dn`. At least one Dirichlet label is required to remove the constant-potential nullspace.
- Validated imported 2D domains use triangle/bilinear-quadrilateral finite-element stiffness assembly with lumped nodal charge, physical-label Dirichlet constraints, CSR storage, and Jacobi-preconditioned conjugate gradients. The solve reports its residual and iteration count, and projects element electric fields back to nodes. This API is not yet connected to the CLI simulation loop.

## Particle advance and steady state

AuroraPIC advances particles with a time-centered electrostatic leapfrog update. Each step deposits charge at particle positions, solves Poisson's equation for `E^n`, kicks stored half-step velocities, drifts positions to the next time level, reapplies particle boundaries, redeposits charge, resolves fields, and synchronizes the public velocity fields used by diagnostics/output. Reflecting multidimensional particle boundaries reverse the normal half-step velocity so subsequent diagnostics remain consistent after synchronization. All dimensions support transient and steady-state execution. Steady-state mode is an engineering stop condition based on the relative change between adjacent total-energy diagnostic windows, not a proof of physical equilibrium; use conservative tolerances/windows and inspect the emitted diagnostics for oscillatory systems.

The C++ imported-domain runtime follows the same cycle. It samples particles uniformly by decomposing validated cells into area-weighted triangles, tracks the earliest outward intersection of each drift segment with tagged geometry, and applies label-specific absorption or specular reflection. This segment-based check is required for concave domains where an endpoint can lie inside even though its path crossed outside. Absorbed counts are reported by physical label.

Imported runs are reachable from the CLI with `mesh = imported`. Their checkpoints include a deterministic signature over node coordinates, element connectivity, physical tags, and labels; a restart is rejected if the configured geometry differs.

## Collisions

Optional Monte-Carlo BGK-like velocity reset collisions model scattering against a prescribed neutral bath. The collision probability is `1 - exp(-nu dt)`.

## Stability guidance

For credible runs, choose:

- cell size small enough to resolve the Debye length;
- time step below electron plasma, cyclotron when using prescribed magnetic fields, and transit-time scales;
- enough particles per cell to control PIC noise;
- output and steady-state windows long enough to distinguish relaxation from oscillation.

The bounded smoke/performance envelope for the checked-in examples is documented in `docs/performance-envelope.md`; use it as an operability baseline, not as proof that arbitrary larger plasma cases are converged.

## Verification included

The automated test suite checks the periodic spectral Poisson solve against an analytic sinusoidal charge distribution, checks the imported finite-element solve against constant-potential, symmetric-source, and exact mixed-boundary linear solutions, verifies that repeated imported solves reuse one assembled operator without changing the numerical result, checks particle-location cache population, reuse, and cross-cell fallback, rejects singular or incomplete imported boundary specifications, exercises imported-domain multi-bounce reflection and label-attributed absorption end to end, and runs a short neutral two-species structured PIC simulation.

## Extension path

AuroraPIC should extend from 1D to 2D first, then 3D, while keeping the current 1D model as a regression and validation target. The multidimensional plan is documented in `docs/multidimensional-roadmap.md`.
