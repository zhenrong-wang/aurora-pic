# AuroraPIC methodology

AuroraPIC v0.1 implements deliberately bounded, auditable electrostatic
`1D1V`/`1D3V`, planar structured/imported `2D3V`, and structured `3D3V`
Particle-in-Cell models for plasma dynamics research and engineering studies.

## Model

Particles carry position `x`, velocity `vx`, charge, mass, and macro-particle
weight in the backward-compatible 1D1V mode. With
`velocity_dimensions = 3`, they additionally carry `vy` and `vz`. Position,
charge deposition, and electrostatic acceleration remain one-dimensional;
the transverse components participate in initialization, kinetic-energy
diagnostics, BGK relaxation, isotropic MCC scattering, and restart. Charge is
deposited to a one-dimensional mesh with cloud-in-cell weighting. In 1D CLI
configuration, `weight` is the direct macro-particle weight; alternatively,
if `weight` is omitted, `density` is converted to
`weight = density * initialization_width / particles` so the represented line
density is explicit. The electrostatic field is obtained from Poisson's
equation,

```text
d²phi/dx² = -rho / eps0,      E = -dphi/dx
```

with the configured homogeneous permittivity. Normalized mode uses a base
permittivity of `1`; SI mode uses `8.8541878128e-12 F/m`. See
`docs/units.md` for reduced-dimensional weight and energy conventions.

The 1D runtime optionally assigns a positive integer
`timestep_multiplier = N` to each species. That species is pushed,
boundary-checked, and collision-processed on pre-step indices `0, N, 2N, ...`
with an effective interval `N*dt`; its charge distribution is held at the last
advanced position on intervening field solves. Newly created particles inherit
their product species interval. This is a multiple-timestep approximation, not
an automatic accuracy guarantee: transit distance, acceleration, and MCC
probability must be acceptable at the enlarged interval. The default is one.

## Initial-value problem

A new transient run constructs its electrostatic initial-value problem in four
ordered operations: load each particle population, deposit its charge, solve
Poisson's equation using the configured field boundaries, and initialize the
leapfrog or Boris half-step velocity. The initial electric field is therefore
consistent with the deposited charge and boundary values; AuroraPIC does not
accept an unrelated arbitrary electric field that could silently violate
Gauss's law. A checkpoint restart is different: it restores an already evolved
particle state, half-step state, time, and RNG engine rather than constructing a
new state at `t = 0`.

Species initialization has a separately versioned
`initialization_version = 1` contract. `loading = random` preserves the
historical independent uniform-position and Gaussian-velocity sampling.
`loading = quiet_start` places each structured-coordinate marginal at the
center of a deterministic equal-measure stratum and generates thermal
velocities in antithetic pairs. Every even-sized quiet-start population
therefore has its configured mean drift exactly, up to floating-point roundoff,
while retaining a sampled thermal spread. Imported 2D quiet-start loading
stratifies area-weighted cell selection and triangle coordinates, so all
particles remain in the validated geometry.

`density_profile = uniform` retains those historical position paths exactly.
`gaussian` defines a relative density proportional to
`exp(-0.5 * sum(((x_i-center_i)/scale_i)^2))` in physical/configured
coordinates and requires a center and positive scale on every active spatial
axis. `sinusoidal` defines a relative density
`1 + amplitude*cos(2*pi*sum(mode_i*xi_i) + phase)`, where each `xi_i` is
normalized over the active initialization envelope, modes are non-negative
integers with at least one nonzero component, phase is in radians, and
`|amplitude| <= 1` keeps density non-negative. The requested macro-particle
count and weight still fix the total represented population; profiles are
normalized spatial shapes rather than new weight conventions.

Nonuniform profiles use acceptance sampling from the existing uniform geometry
measure. Random loading uses pseudorandom candidates. Quiet-start uses a
deterministic low-discrepancy sequence for position and acceptance coordinates,
while retaining antithetic velocity pairs. Every species has an explicit
`max_profile_sampling_attempts` total work budget. A profile that cannot load
the requested population within that budget fails without returning a partial
or silently biased state.

For imported meshes, `initialization_region` selects a named Gmsh
dimension-two physical group. AuroraPIC builds a separate cumulative area
distribution for every cell label and samples only triangles belonging to the
requested label. This preserves area-uniform random and quiet-start loading
within disconnected or non-rectangular physical regions. A named region cannot
be combined with rectangular initialization bounds; unknown labels and
zero-area selections fail before the run starts.

Initialization acceptance gates are opt-in and run after generated-state or
checkpoint moments are collected but before scalar diagnostics or time
integration. `initialization_max_relative_charge_imbalance` bounds
`|sum(Q_i)| / sum(|Q_i|)`. The current gate bounds
`|sum(Q_i * mean(v_i))| / sum(|Q_i| * |mean(v_i)|)` over the active velocity
components. A zero numerator and denominator has zero residual. Explicit
`initialization_charge_pairs = first:second,...` require opposite represented
charge signs and bound the relative mismatch between their represented-charge
magnitudes with `initialization_max_relative_pair_imbalance`. All tolerances
are finite dimensionless values in `[0, 1]`; gates are disabled when their
keys are omitted. Every run writes `initialization_acceptance.csv`, including
disabled, passing, and failing states. A failed enabled gate leaves that audit
and `initialization.csv` in the output directory, then aborts before stepping.

`thermal_velocity` remains the backward-compatible isotropic Gaussian standard
deviation. Optional `thermal_velocity_x`, `thermal_velocity_y`, and
`thermal_velocity_z` values override it component by component; 1D1V accepts
only the x override, while 1D3V, 2D3V, and 3D3V accept all three. Values are
velocities in the configured unit system, not temperatures. Converting a
physical temperature requires the species mass and the documented unit
contract.

As a strict SI convenience, structured 2D species and its pair/current source
controls may specify temperatures in electron-volts. The loader resolves each
to the isotropic one-component Maxwellian standard deviation
`sqrt(e*T_eV/m)` using the configured species mass. It rejects component
overrides combined with a species temperature, velocity/temperature pairs on
a source, and all eV forms in normalized units.

Every run emits `initialization.csv` after generated initialization, external
particle-state loading, or checkpoint restoration. Its versioned rows report
the state source, species, loading and density-profile models or
`external`/`restart`, selected region where meaningful, live
macro-particle count and weight, represented physical-particle number and
charge, mean position and velocity, and realized component position and
velocity standard deviations.
This is an audit of the actually loaded numerical state rather than a copy of
requested configuration values.

The current imported quiet-start implementation intentionally rejects a
rectangular `init_x_*`/`init_y_*` clip instead of silently degrading to random
or biased rejection sampling. Random loading continues to support those
bounds. The versioned external `.aps` path validates time-centered particle
records and then rebuilds the field-consistent half step. Its deterministic
semantic signature can be pinned in configuration, and external runs record
the resolved source plus realized signature in `initial_state_metadata.txt`.
The public writer preserves that signature across a text round trip. Runtime
ingestion verifies integrity before delivering records directly into
simulation-owned arrays; reader-side auxiliary memory scales with species
count instead of particle count. The portable text backend uses repeated
species scans to preserve canonical signatures without buffering the
population. Physical-temperature inputs, general tabulated profiles, and a
chunked openPMD/HDF5 particle-state backend remain subsequent
initial-condition milestones.

## Field solvers

- Periodic domains use a direct spectral Poisson solve. The zero mode is removed, enforcing global quasi-neutral compatibility.
- Dirichlet domains use a tridiagonal finite-difference Poisson solve with
  prescribed endpoint potentials. The 1D path supports independent
  sinusoidal endpoint drives; the solve after particle drift evaluates them
  at `t^(n+1)`, matching the new charge-density time level.
- Structured 2D domains with exactly one periodic and one Dirichlet axis use a
  direct spectral-tridiagonal solve. The periodic transform uses radix-2,
  mixed-radix, or Bluestein algorithms according to its extent; each Fourier
  mode then produces an independent complex tridiagonal Dirichlet system.
  Electric fields use the same centered periodic and centered/one-sided
  Dirichlet finite differences as particle gathering.
- Imported finite-element domains support label-wise constant Dirichlet values and Neumann outward normal derivatives. For `-laplacian(phi) = rho / epsilon_0`, the Neumann term `dphi/dn` enters the weak-form right-hand side with a positive sign; therefore `E dot n = -dphi/dn`. At least one Dirichlet label is required to remove the constant-potential nullspace.
- Validated imported 2D domains use triangle/bilinear-quadrilateral finite-element stiffness assembly with lumped nodal charge, physical-label Dirichlet constraints, CSR storage, and Jacobi-preconditioned conjugate gradients. The solve reports its residual and iteration count, projects element electric fields back to nodes, and is integrated with the imported-mesh CLI simulation loop.

## Particle advance and steady state

AuroraPIC advances particles with a time-centered electrostatic leapfrog update. At the default species multiplier of one, each step deposits charge at particle positions, solves Poisson's equation for `E^n`, kicks stored half-step velocities, drifts positions to the next time level, reapplies particle boundaries, redeposits charge, resolves fields at `t^(n+1)`, and synchronizes the public velocity fields used by diagnostics/output. A subcycled 1D species performs those particle operations only on its scheduled base steps while participating in every charge deposition at its held position. A time-dependent 1D electrode is evaluated at the same field time level; its realized left/right values are recorded with every scalar sample, and restart recovers phase from the stored time. The 1D collision split is explicit: `collision_velocity_sampling = time_centered` retains the synchronized-velocity default, while `leapfrog_half_step` applies collisions to the drift velocity and synchronizes only afterward. The checkpoint collision identity protects this choice. When phase-EEDF particle history is enabled, the threshold ledger also records field-push promotions and demotions between the post-collision state entering the mover and the pre-collision state leaving it; the energy representation follows the configured collision-velocity sampling contract. These counts isolate mover-induced threshold changes, whereas consecutive synchronized interstep transitions span the preceding collision and following push. A configurable subthreshold band additionally accumulates promotion probability and signed, positive, and negative mover work. Checkpoint v22 preserves the restart-safe field-push counts and band-work accumulators and rejects unsupported continuation from older checkpoints unless spatial diagnostics are explicitly reset. Reflecting multidimensional particle boundaries reverse the normal half-step velocity so subsequent diagnostics remain consistent after synchronization. For 1D Dirichlet absorption, `boundary_losses.csv` records cumulative species/side macro-particle count, represented signed charge, and full-velocity kinetic energy evaluated from the crossing half-step longitudinal velocity and current transverse velocities. `power_transfer.csv` separately accumulates each species' represented kinetic-energy change caused only by the electric push, including the final push of an absorbed particle. Dividing a work difference by its time window gives the discrete, domain-integrated counterpart of the time-averaged `J_s · E` power per area; collision and wall-carried energy are not folded into it. Ordered per-worker reductions keep both diagnostics race-free and deterministic. Checkpoint v7 preserves the power and wall counters and their origins; legacy restarts explicitly begin unsupported counter coverage at their restart step. All dimensions support transient and steady-state execution. Steady-state mode is an engineering stop condition based on the relative change between adjacent total-energy diagnostic windows, not a proof of physical equilibrium; use conservative tolerances/windows and inspect the emitted diagnostics for oscillatory systems.

Checkpoint v23 extends the promotion band with origin total and longitudinal
energy and the exact identity
`delta K = m v_x delta v_x + 0.5 m delta v_x^2`. The signed linear term
measures velocity--increment alignment, while the nonnegative quadratic term
is proportional to the squared particle-sampled electric field under half-step
sampling. Restart validates both `signed = positive - negative` closures and
`signed total = linear + quadratic`. A v22 state remains usable when spatial
diagnostics are explicitly reset, but cannot continue the new accumulators.

The C++ imported-domain runtime follows the same cycle. It samples particles uniformly by decomposing validated cells into area-weighted triangles, tracks the earliest outward intersection of each drift segment with tagged geometry, and applies label-specific absorption or specular reflection. This segment-based check is required for concave domains where an endpoint can lie inside even though its path crossed outside. Absorbed counts are reported by physical label.

Imported runs are reachable from the CLI with `mesh = imported`. Their checkpoints include a deterministic signature over node coordinates, element connectivity, physical tags, and labels; a restart is rejected if the configured geometry differs.

Imported boundary sources inject a configured integer number of macro-particles at the beginning of each active timestep. Boundary faces sharing a physical label are selected by length, then sampled uniformly along the selected face and inset into the adjacent domain. The configured normal velocity is inward; thermal normal speed uses an inward half-range Gaussian magnitude, while tangential and out-of-plane thermal velocities remain signed Gaussian. New particles are initialized into the same leapfrog/Boris time staggering as initial particles before joining the push. Dead storage slots are reused deterministically. Imported checkpoint v6 records the source identity, schedule, three-component velocity parameters, cumulative counts, dynamic particle state, RNG state, and unit contract so continuation reproduces uninterrupted injection.

Structured 2D volumetric pair sources run at the same beginning-of-step point.
Every event samples one position inside its configured rectangle and creates
both species there. Uniform, Gaussian, and sinusoidal profile envelopes reuse
the bounded rejection sampler used by initialization; these profiles are
normalized probability shapes and therefore do not change the configured
integral source rate. The species have independent three-component drift and
isotropic thermal velocities, but must carry opposite equal charge and equal
macro weight. Thus deposited source charge cancels exactly at creation and
one macro-pair has one well-defined represented physical-pair count.

A source selects an integer macro-pair count, a total represented-pair rate,
or a peak volumetric rate. The peak form is multiplied by the analytic
normalized-profile area and global extrusion depth. For rate forms, each
active step adds
`rate * dt / macro_weight` to a deterministic fractional accumulator, creates
the integer part, and carries the remainder. Capacity is checked for both
species before storage changes; dead slots are reused deterministically.
Diagnostics report cumulative macro/represented pairs, the fractional
remainder, sampled full-3V kinetic energy, and configured rate. Structured
checkpoint v10 fingerprints extrusion depth and the complete source definition and preserves the
accumulator, diagnostics, particle state, and RNG. This is a prescribed
source, not a collision, neutral-depletion, recoil, or reaction-network model.

Structured 2D current regulation accumulates represented charge crossing one
absorbing boundary, resolved by species rather than by macro-particle count.
The new charge since the previous control update is divided by the configured
emitted species charge and macro weight. A positive integer count is emitted
uniformly along an inset plane. Cumulative mode carries signed fractional
surplus or debt for generic unequal-weight regulation. Timestep-local mode
clears reverse demand every step and carries only positive fractional
remainder; with equal species macro weights this is the LANDMARK
electron-minus-ion anode-loss rule.
`current_source.csv` separates cumulative negative and positive monitored
charge and records controller updates, reverse-demand steps and fraction,
cumulative reverse demand, and the largest single-step reverse demand in
macro-particle equivalents. It additionally bins reverse events into one-,
two-, and multi-macroparticle requests, reports demand mean/RMS, and accumulates
the negative and positive monitored boundary charge on those timesteps.
Checkpoint v10 preserves these counters. When a v1-v8 checkpoint is loaded,
the reverse-total start step is the restart step; when a v1-v9 checkpoint is
loaded, the distribution start step is the restart step. Partial coverage
therefore cannot be mistaken for whole-run statistics.

Every structured-2D scalar-output window also writes species-resolved loss
rows to `boundary_flux.csv`. Each row identifies the window start/end,
species, and boundary and contains incremental and cumulative macroparticle
loss, represented particles and charge, and their rates. On a fresh or
restarted run the first rows have zero duration and zero increment; subsequent
windows are differences of the checkpoint-preserved cumulative species loss
counters. Thus restart does not blend pre-restart losses into the first rate.

An optional x/y line-average potential reference has two modes. Gauge mode
applies a constant potential offset after each Poisson solve and leaves the
electric field invariant. Affine mode preserves the zero-coordinate
electrode, linearly reaches the requested internal-plane target, and adds the
exact constant correction to the corresponding electric-field component.
Checkpoint v10 preserves both mode choices, boundary counters, controller
accumulator, diagnostics, and RNG state.

Opt-in structured-2D resolved diagnostics reduce fields along one profile axis
using nodal control-width quadrature and deposit species number plus first and
second three-velocity moments to the same profile nodes with one-dimensional
CIC weights. Number density uses the node control volume, including the
configured extrusion depth. Mean velocity, component thermal speeds, scalar
SI temperature, and charge-current density are derived from those deposited
moments. Time averages trapezoidally integrate the field profiles and raw
density-weighted particle moments over the actual sample times before deriving
the averaged observables.

On the distinct periodic mode axis, field coefficients use nodal-area
quadrature. Species number and current are deposited once to the periodic
diagnostic grid and transformed there, avoiding a particle scan per requested
mode. The mode history retains complex coefficients as well as one-sided
amplitudes so offline analysis can recover frequency and propagation phase.
The reduced Hall smoke guards exact manufactured sine/cosine coefficients and
the complete output schema, but its short averaging window has no discharge
physics significance.

Absorbing impacts are recorded in parallel and then sorted by incident species and particle ID. This deterministic reduction accumulates species/tag-resolved macro-particle count, represented physical-particle count, charge, full three-velocity incident kinetic energy, last-step rate, and rate per tagged-boundary length. Configured secondary-emission rules are then evaluated serially. Their physical yield is converted through the incident/emitted macro weights, with stochastic rounding for fractional macro-particles and explicit per-impact/storage limits. Emitted velocities use the same inward half-range normal and signed tangential/out-of-plane distributions as sources, and emitted particles enter the pusher at the boundary-hit position inset into the domain. Imported checkpoint v6 preserves emission definitions, cumulative emitted counts, flux state, particle state, RNG state, unit metadata, and optional MCC state.

Imported-mesh quality reporting computes cell-area and edge-length extrema, the minimum cell-corner angle, and the maximum within-cell edge-length ratio. These inexpensive metrics complement the mandatory finite, nondegenerate, convex, consistently oriented, manifold, and exactly tagged boundary validation. The biased-probe integration mesh pins explicit angle, edge-ratio, area, topology, and physical-label envelopes.

## Collisions

The historical optional BGK velocity-reset model remains available for
compatibility; in 1D3V it redraws all three velocity components. The tabulated
1D and imported 2D3V MCC paths evaluate
`nu_i(E) = neutral_density * sigma_i(E) * relative_speed` and use exponential
null-collision candidate times with a strictly enforced user-supplied maximum
frequency. Positive-temperature SI runs sample a Maxwellian neutral velocity
bounded at eight component standard deviations and enforce a conservative
rate majorant over the reachable relative-speed interval. Elastic events with
gas mass metadata use
two-body kinematics that conserve projectile-plus-neutral momentum and energy;
excitation events remove a configured threshold energy. The 1D1V path
randomizes velocity sign at fixed speed; the 1D3V and imported paths sample
an isotropic three-dimensional relative direction. Named
interval/cumulative counts are written to `collisions.csv`. Named 1D collision
models may target distinct species simultaneously; their products are staged
until all targets finish. Equal-weight electron/ion pairs from 1D3V ionization
therefore conserve macro-charge and cannot collide in their birth timestep.
Capacity is preflighted before product insertion. Checkpoint v4 in
1D and v6 for imported geometry fingerprint effective tables and preserve
counters. 1D3V and imported 2D3V ionization remove their threshold, equally partition
the remaining primary/secondary electron energy, samples independent isotropic
directions, and creates an opposite-charge ion at the sampled target-neutral
velocity with bounded,
preflighted storage. Imported runs can load versioned external `.gas`
manifests; packaged physics and provenance are kept separate from
simulation-specific density, temperature, rate bounds, and reactive species
mappings. Effective dataset metadata and channel settings are emitted to
`collision_data.txt` and included in restart compatibility. See
`docs/collisions.md` for the gas metadata, scaling, reactive-species
constraints, and limitation contract.

Imported resonant charge exchange requires ion and neutral masses to match and
maps the tracked fast ion onto the sampled neutral velocity. This preserves
ion count and charge while transferring the fast neutral product to the
implicit reservoir.

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
