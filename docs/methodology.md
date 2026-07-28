# AuroraPIC methodology

AuroraPIC v0.1 implements deliberately bounded, auditable electrostatic `1D1V`, planar structured/imported `2D3V`, and structured `3D3V` Particle-in-Cell models for plasma dynamics research and engineering studies.

## Model

Particles carry position `x`, velocity `vx`, charge, mass, and macro-particle weight. Charge is deposited to a one-dimensional mesh with cloud-in-cell weighting. In 1D CLI configuration, `weight` is the direct macro-particle weight; alternatively, if `weight` is omitted, `density` is converted to `weight = density * initialization_width / particles` so the represented line density is explicit. The electrostatic field is obtained from Poisson's equation,

```text
d²phi/dx² = -rho / eps0,      E = -dphi/dx
```

with the configured homogeneous permittivity. Normalized mode uses a base
permittivity of `1`; SI mode uses `8.8541878128e-12 F/m`. See
`docs/units.md` for reduced-dimensional weight and energy conventions.

## Field solvers

- Periodic domains use a direct spectral Poisson solve. The zero mode is removed, enforcing global quasi-neutral compatibility.
- Dirichlet domains use a tridiagonal finite-difference Poisson solve with prescribed endpoint potentials.
- Imported finite-element domains support label-wise constant Dirichlet values and Neumann outward normal derivatives. For `-laplacian(phi) = rho / epsilon_0`, the Neumann term `dphi/dn` enters the weak-form right-hand side with a positive sign; therefore `E dot n = -dphi/dn`. At least one Dirichlet label is required to remove the constant-potential nullspace.
- Validated imported 2D domains use triangle/bilinear-quadrilateral finite-element stiffness assembly with lumped nodal charge, physical-label Dirichlet constraints, CSR storage, and Jacobi-preconditioned conjugate gradients. The solve reports its residual and iteration count, projects element electric fields back to nodes, and is integrated with the imported-mesh CLI simulation loop.

## Particle advance and steady state

AuroraPIC advances particles with a time-centered electrostatic leapfrog update. Each step deposits charge at particle positions, solves Poisson's equation for `E^n`, kicks stored half-step velocities, drifts positions to the next time level, reapplies particle boundaries, redeposits charge, resolves fields, and synchronizes the public velocity fields used by diagnostics/output. Reflecting multidimensional particle boundaries reverse the normal half-step velocity so subsequent diagnostics remain consistent after synchronization. All dimensions support transient and steady-state execution. Steady-state mode is an engineering stop condition based on the relative change between adjacent total-energy diagnostic windows, not a proof of physical equilibrium; use conservative tolerances/windows and inspect the emitted diagnostics for oscillatory systems.

The C++ imported-domain runtime follows the same cycle. It samples particles uniformly by decomposing validated cells into area-weighted triangles, tracks the earliest outward intersection of each drift segment with tagged geometry, and applies label-specific absorption or specular reflection. This segment-based check is required for concave domains where an endpoint can lie inside even though its path crossed outside. Absorbed counts are reported by physical label.

Imported runs are reachable from the CLI with `mesh = imported`. Their checkpoints include a deterministic signature over node coordinates, element connectivity, physical tags, and labels; a restart is rejected if the configured geometry differs.

Imported boundary sources inject a configured integer number of macro-particles at the beginning of each active timestep. Boundary faces sharing a physical label are selected by length, then sampled uniformly along the selected face and inset into the adjacent domain. The configured normal velocity is inward; thermal normal speed uses an inward half-range Gaussian magnitude, while tangential and out-of-plane thermal velocities remain signed Gaussian. New particles are initialized into the same leapfrog/Boris time staggering as initial particles before joining the push. Dead storage slots are reused deterministically. Imported checkpoint v6 records the source identity, schedule, three-component velocity parameters, cumulative counts, dynamic particle state, RNG state, and unit contract so continuation reproduces uninterrupted injection.

Absorbing impacts are recorded in parallel and then sorted by incident species and particle ID. This deterministic reduction accumulates species/tag-resolved macro-particle count, represented physical-particle count, charge, full three-velocity incident kinetic energy, last-step rate, and rate per tagged-boundary length. Configured secondary-emission rules are then evaluated serially. Their physical yield is converted through the incident/emitted macro weights, with stochastic rounding for fractional macro-particles and explicit per-impact/storage limits. Emitted velocities use the same inward half-range normal and signed tangential/out-of-plane distributions as sources, and emitted particles enter the pusher at the boundary-hit position inset into the domain. Imported checkpoint v6 preserves emission definitions, cumulative emitted counts, flux state, particle state, RNG state, unit metadata, and optional MCC state.

Imported-mesh quality reporting computes cell-area and edge-length extrema, the minimum cell-corner angle, and the maximum within-cell edge-length ratio. These inexpensive metrics complement the mandatory finite, nondegenerate, convex, consistently oriented, manifold, and exactly tagged boundary validation. The biased-probe integration mesh pins explicit angle, edge-ratio, area, topology, and physical-label envelopes.

## Collisions

The historical optional BGK velocity-reset model remains available for
compatibility. The tabulated 1D and imported 2D3V MCC paths evaluate
`nu_i(E) = neutral_density * sigma_i(E) * speed` and use exponential
null-collision candidate times with a strictly enforced user-supplied maximum
frequency. Elastic events conserve particle kinetic energy; excitation events
remove a configured threshold energy. The 1D path randomizes velocity sign;
the imported path samples an isotropic three-dimensional direction. Named
interval/cumulative counts are written to `collisions.csv`. Checkpoint v3 in
1D and v6 for imported geometry fingerprint effective tables and preserve
counters. Imported 2D3V ionization removes its threshold, equally partitions
the remaining primary/secondary electron energy, samples independent isotropic
directions, and creates a stationary opposite-charge ion with bounded,
preflighted storage. See `docs/collisions.md` for the gas metadata, scaling,
reactive-species constraints, and limitation contract.

## Stability guidance

For credible runs, choose:

- cell size small enough to resolve the Debye length;
- time step below electron plasma, cyclotron when using prescribed magnetic fields, and transit-time scales;
- enough particles per cell to control PIC noise;
- output and steady-state windows long enough to distinguish relaxation from oscillation.

The bounded smoke/performance envelope for the checked-in examples is documented in `docs/performance-envelope.md`; use it as an operability baseline, not as proof that arbitrary larger plasma cases are converged.

## Verification included

The automated test suite checks the periodic spectral Poisson solve against an analytic sinusoidal charge distribution, checks the imported finite-element solve against constant-potential, symmetric-source, and exact mixed-boundary linear solutions, verifies that repeated imported solves reuse one assembled operator without changing the numerical result, checks particle-location cache population, reuse, and cross-cell fallback, rejects singular or incomplete imported boundary specifications, exercises imported-domain multi-bounce reflection and label-attributed absorption, verifies bounded boundary injection, weight-aware secondary emission, species-resolved physical fluxes, deterministic serial/OpenMP behavior, source/emission restart, and runs a short neutral two-species structured PIC simulation.

## Extension path

AuroraPIC should extend from 1D to 2D first, then 3D, while keeping the current 1D model as a regression and validation target. The multidimensional plan is documented in `docs/multidimensional-roadmap.md`.
