# AuroraPIC

AuroraPIC is a C++20 starting point for scientific plasma dynamics simulation. The current codebase implements electrostatic `1D1V` and `1D3V`, planar `2D3V` (structured and imported geometry), and structured `3D3V` Particle-in-Cell (PIC) paths with configurable species, periodic or Dirichlet boundaries, transient fixed-step and steady-state convergence modes, scalar diagnostics, and text checkpoint/restart files. The 1D path additionally provides restart-safe, species-resolved spatial density averaging with optional whole-RF-cycle alignment. The 1D baseline provides the historical BGK relaxation model plus tabulated null-collision MCC; 1D3V supports simultaneous named collision targets, isotropic elastic/excitation scattering, and bounded charge-conservative electron-impact ionization products. Imported 2D3V runs additionally support bounded-Maxwellian finite-mass neutrals, electron attachment, and resonant ion-neutral charge exchange. The multidimensional paths provide prescribed uniform or one-coordinate tabulated magnetic-field Boris pushes, VTK field output, side-specific particle boundaries, and optional particle inspection CSVs. A separate homogeneous electron-swarm runner scans reduced electric field with the same three-velocity collision kernel before a gas package is used in a device geometry.

## Why this methodology

PIC is a standard kinetic method for low-collisionality plasma simulation because it evolves macro-particles while solving fields on a mesh. AuroraPIC keeps the 1D path as the regression baseline while adding structured 2D building blocks before moving to geometry import, 3D, electromagnetic fields, and accelerated parallel backends.

For the recommended multidimensional expansion strategy, geometry/mesh format choices, and staged implementation plan, see `docs/multidimensional-roadmap.md`.

The first end-to-end nontrivial geometry case is documented in [`docs/real-case-validation.md`](docs/real-case-validation.md). It includes a Gmsh-authored chamber with an internal circular biased probe, a committed regenerable mesh, mesh-quality gates, tagged-boundary physics, and a deterministic simulation acceptance envelope. It is an integration-grade real geometry case, not yet an experimentally validated device model.

The quantitative system-level physics cases are documented in
[`docs/kinetic-validation.md`](docs/kinetic-validation.md). Deterministic,
resource-bounded 1D Landau-damping, 1D two-stream, and orthogonal 2D/3D Langmuir
runs measure damping rate, frequency, instability growth, nonlinear turnover,
directional symmetry, and total-energy drift against analytic or published
Vlasov-Poisson behavior.

The next device-physics benchmark is defined in
[`docs/ccp-validation.md`](docs/ccp-validation.md). It pins the published
Turner helium CCP parameters, records completed RF-electrode, 1D3V,
multispecies-MCC, and ionization support, and lists the remaining ion-scattering and
statistical-diagnostic gates that prevent a premature validation claim.

The Hall-effect-thruster verification and validation ladder is defined in
[`docs/hall-thruster-validation.md`](docs/hall-thruster-validation.md). It pins
the public LANDMARK axial-azimuthal and radial-azimuthal PIC cases, identifies
public WarpX reference output and NASA HERMeS measurements, provides a guarded
single-thread runtime qualification before costly campaigns, separates
code-to-code verification from experimental validation, and records the
physics, diagnostics, parallelism, and resource gates still required. The
installed `compare_hall.py` tool performs checksum-pinned, uncertainty-aware
profile and complex-mode comparisons against local reference data;
`preflight_hall.py` estimates production memory, storage, output rows, and
particle-update work without launching a simulation.
`qualify_hall_runtime.py` runs only a capped one-thread micro/workstation slice
and records binary/deck provenance plus initial-population cost projections.
`prepare_hall_campaign.py` requires an explicit production-cost
acknowledgement before it writes—but never launches—a full Case 2 candidate
deck.
`analyze_hall_pilot.py` checks bounded-tier diagnostic integrity, including
controllable cathode remainder and separately unserved reverse demand, and always reports
`physics_claim = none`.
`lock_hall_source.py` plans external acquisition from the committed source
registry and stream-hashes local artifacts without downloading them.
`verify_turner_source.py` checksum-locks and structurally validates a local
copy of the Turner helium CCP publisher supplement while keeping its
rights-restricted tables outside the repository.
`normalize_turner_source.py` then produces local, audit-hashed AuroraPIC gas
manifests and all four original/refined reference profiles without resampling
or committing the restricted values. Named 1D collision models can load these
manifests directly with `gas_data_file`.
`compare_turner.py` implements the paper's ion-density `X²` statistic against
the population standard deviation, locks the local reference through the
normalization audit, and applies the published 95% and 99% case ranges. It can
also compare the published electron-density distribution descriptively while
explicitly suppressing the ion-only acceptance ranges.
`analyze_turner_density_blocks.py` hash-checks consecutive reset 32-cycle
post-benchmark profiles and quantifies density drift, adjacent profile motion,
lag-one correlation, and an AR(1) effective block count without treating a
diagnostic continuation as a published benchmark pass.
`analyze_turner_amplitude_uncertainty.py` is a separate post-protocol design
tool. It fits stationary AR(1) errors to the linear-detrended amplitude from a
completed density-block audit and uses a fixed-seed parametric null ensemble to
ask whether the observed absolute drift is unusual for that correlated process.
It reports residual
correlation and explicit model/decision boundaries; it cannot replace a locked
gate, rescue a failed campaign, or establish stationarity.
`prepare_turner_ensemble.py` atomically creates three or more distinct-seed,
full-duration Turner decks and preflight reports behind a separate aggregate
cost acknowledgement. It never launches a run and records a sequential-only
workstation policy plus aggregate resource floors.
`attach_turner_ensemble_result.py` admits a completed seed only after proving
semantic deck equivalence, independently recomputing the locked Turner
comparison, and checksum-recording its profile, metadata, and final checkpoint.
`analyze_turner_ensemble.py` then requires every prepared seed, revalidates all
artifact hashes and comparisons, and emits a descriptive aggregate without
inventing a post-hoc ensemble acceptance threshold.
`audit_turner_credibility.py` checksum-locks the published-duration ensemble,
numerical-sensitivity, energy-closure, phase, and EEDF evidence into one
fail-closed scientific status. It distinguishes numerical integrity and
descriptive scalar agreement from the unresolved published ion-density
failure, and never promotes post-benchmark diagnostics into formal passes.
The prospective Turner subcycle-policy invariance control additionally proves
that the held-density correction is inactive for this unit-cadence helium
case: six physical diagnostic outputs remain byte-identical, so the existing
three-seed density discrepancy is preserved rather than silently reclassified.
The subsequent specification audit identified a distinct mismatch: Turner
permits one Bernoulli collision opportunity per particle step, while the
general AuroraPIC default uses a Poisson clock. An explicit fingerprinted
`single_bernoulli` mode now passes exposure-normalized statistical tests and
repeatedly moves helium ionization and density toward the Turner reference;
the full corrected equilibrium comparison remains outstanding.
`compare_edupic_phase_space.py` checksum-locks the public eduPIC Figure 11 raw
matrices and directly compares AuroraPIC over position and RF phase. It covers
potential, electric field, electron/ion density, current density, ohmic power,
and mean particle energy without fitted phase shifts or spatial reflection.
The current result is descriptive because the AuroraPIC measurement window is
still transient; it is not a cross-code acceptance or validation claim. The
first locked result is recorded in
`benchmarks/ccp/edupic-argon-phase-space-cycle80-20260812.json`.
Fresh 1D spatial-average windows also write `spatial_collision_rate.csv` and
`spatial_phase_collision_rate.csv`. These deposit represented collision events
with the same shape function as collision energy and expose SI volumetric event
rates, including an eduPIC-compatible phase-resolved ionization-rate field.
The first complete phase-space result, including ionization, is recorded in
`benchmarks/ccp/edupic-argon-phase-space-cycle80-v15-20260812.json`.
It also independently integrates the checksum-locked ionization cross section
over the measured EEDF and compares that prediction with the collision-event
ledger, keeping this post-diagnostic closure separate from acceptance claims.
The later checksum-locked
`benchmarks/ccp/edupic-ionization-exposure-synthesis-20260825.json` links an
independent native phase-region EEDF fold to native cycle histories: their
AuroraPIC/native ionization ratios are `0.8530` and `0.8693`, respectively.
This close descriptive agreement localizes the remaining deficit to
phase-space ionizing exposure while preserving the stated non-validation
claim boundary.
The matched-half-step field-push ledgers further localize reduced energetic
promotion to all three critical RF-phase octants, strongest at phase
`0.25--0.375` and in `x/L=0.2--0.4`; the reproducible post hoc result is
`benchmarks/ccp/edupic-field-push-promotion-localization-20260825.json`.
The phase-EEDF history diagnostic now also supports a configurable
subthreshold promotion band with restart-safe signed/positive/negative mover
work, enabling the next near-threshold causal comparison without changing the
particle trajectory.
The native counterpart passes a byte-exact one-cycle passivity smoke, and the
production comparison is preregistered in
`benchmarks/ccp/edupic-promotion-band-work-rule-20260825.json` before any
four-cycle work output is observed.
The completed prospective comparison is recorded in
`benchmarks/ccp/edupic-promotion-band-work-result-20260825.json`: all gates
pass, both microstates support reduced positive mover work and reduced
conditional promotion conversion, while the stricter band-supply-deficit rule
is not met. This localizes one matched CCP discrepancy; it is not an
experimental validation or a claim of general solver correctness.
Checkpoint v23 and a passive native transform now decompose that mover work
into exact linear alignment and quadratic sampled-field terms while recording
the band's origin energy partition. The next comparison is prospectively
locked in `benchmarks/ccp/edupic-mover-decomposition-rule-20260825.json`; no
production interpretation exists yet.
The three-member native baseline is complete and passes every locked integrity,
closure, population, resource, and repeatability gate in
`benchmarks/ccp/edupic-mover-decomposition-native-result-20260825.json`.
The completed prospective result is recorded in
`benchmarks/ccp/edupic-mover-decomposition-result-20260826.json`. Both locked
microstates support weaker particle-sampled field strength and less favorable
positive alignment; origin energy and longitudinal energy partition remain at
parity. This is a one-case discrepancy localization, not experimental or
general PIC validation.
The separately preregistered grid-field discriminator then finds an AuroraPIC
grid mean-square field ratio of `0.6710` and particle-sampling-factor ratio of
`1.1957`; all gates pass. Weaker grid fields are therefore supported, while
differential avoidance of strong-field locations is not. This remains a
bounded diagnosis in one evolved CCP state.
An exact discrete-Poisson reconstruction subsequently closes to relative RMS
errors below `2.0e-7` in both codes and places the boundary drive at `0.9973`
parity, localizing the mature field gap to net space charge. Independent
electron/ion attribution is strongly cancellation-conditioned and is not
treated as a physical causal separation.
A subsequent three-threshold net-charge/sheath comparison passes all integrity
and repeatability gates but supports none of its preregistered gross-mechanism
rules. The remaining mature difference is therefore associated with subtler
phase/spatial net-charge organization, not a large sheath-width or integrated
positive-charge mismatch.
A cycle-resolved follow-up then shows that the phase-0.3 regional field-energy
deficit is already present in the first measured cycle and persists in all four
cycles for both AuroraPIC microstates (`0.658--0.795` relative to the native
ensemble). Integrity and cycle-stability gates pass, but the formal result is
inconclusive because the preregistered pointwise repeatability gate fails at
several phase-0.2 cells and one charge-density cell. The narrower phase-0.3
field comparison is repeatable in both codes; it is a useful localization, not
a substitute for the failed joint gate or a validation claim.
A collision-free common-particle-state trace now closes the initial charge and
field profiles to relative RMS errors of `2.6e-14` and `2.7e-9`. The localized
field metric crosses its locked 2% divergence band after one step, while
populations remain identical through 100 steps. Source inspection identifies
ion-density refresh staggering as the leading mechanism candidate: eduPIC
holds the pre-push ion density between 20-step ion updates, whereas AuroraPIC
redeposits the moved ions. A prospective lag-control branch is still required
before treating that association as causal.
The prospectively declared later-window discriminator is summarized in
`benchmarks/ccp/edupic-argon-heating-trend-cycle80-to116-20260812.json`.
The subsequent strict continuation and cycle-148 horizon-sufficiency decision
are recorded in `benchmarks/ccp/edupic-argon-post-trend-horizon-cycle128-20260812.json`.
`prepare_turner_sensitivity.py` creates a non-launching, staged refinement
matrix with predeclared interpretation thresholds after a systematic Turner
density discrepancy is established.
`normalize_hall_reference.py` verifies local raw-table and case hashes,
performs explicit unit conversion and multi-code envelope reduction, and
produces comparator-ready reference artifacts without interpolation.
`prepare_hall_ensemble.py` writes three or more independent seeded decks
atomically; `aggregate_hall_ensemble.py` applies a conservative 95% Student-t
acceptance test to their checksum-bound comparison reports.
`prepare_hall_convergence.py` writes a cost-gated, five-stage workstation
population/duration ladder without launching it; `analyze_hall_convergence.py`
checks axial profiles and azimuthal spectra for decreasing three-level change.
Use `aurorapic_cli --validate-only <config.cfg>` to parse and validate any
deck without constructing a simulation or taking a timestep.
The CLI also blocks configurations exceeding 100 million estimated initial
particle updates unless the operator supplies the documented
`--allow-large-run` acknowledgement.

The dimensional contract is defined in [`docs/units.md`](docs/units.md). Configurations may select `units = normalized` or `units = si` plus a positive homogeneous `relative_permittivity`. Legacy omission remains normalized; maintained examples are explicit. SI 1D and imported 2D retain per-unit omitted measures; structured planar 2D uses an explicit extrusion depth and reports total particle/field energy over that volume.

The external initial-state contract is defined in [`docs/external-particle-state.md`](docs/external-particle-state.md). A strict versioned `.aps` interchange can initialize time-centered particle positions and velocities in structured 1D/2D/3D or imported 2D geometry while retaining configuration-owned species charge, mass, and constant macro weight. Runtime initialization uses a validated bounded-memory record consumer, so the text reader does not retain a second particle-sized record container. It remains a portable preprocessing/interoperability backend rather than a substitute for the planned high-volume openPMD/HDF5 backend.

The collision contract is defined in [`docs/collisions.md`](docs/collisions.md). The cross-section MCC path uses strict two-column tables, relative path resolution, explicit column scales, conservative maximum-frequency enforcement, bounded Maxwellian SI neutral motion, named-channel diagnostics, and restart fingerprints. Its committed data are synthetic validation inputs, not material data.

## Production milestone baseline

Production readiness is now tracked as explicit milestones instead of an open-ended roadmap narrative. The pinned milestone ladder and evidence expectations live in `docs/multidimensional-roadmap.md#production-readiness-milestone-ladder`; `scripts/validate_milestones.py` is part of the smoke suite and fails if those milestone IDs or README linkage drift. The current M6 baseline includes release-engineering and operability mitigation: configs may declare `config_version = 1`, public examples do so explicitly, unsupported future config versions fail with a clear diagnostic, CI workflow coverage exercises serial/OpenMP build variants across Linux/macOS runners, CPack can produce a `TGZ` package, and `docs/performance-envelope.md` documents the verified smoke/performance envelope. It also preserves the completed M5 higher-fidelity physics slice, M4 runtime-scaling interface slice, M3 VTK XML structured-grid output compatibility, and the M2 tagged 2D Gmsh v2 ASCII importer (`ImportedMesh2D`) for externally meshed planar domains.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 1
```

OpenMP support is enabled by default when CMake finds a C++ OpenMP toolchain. Disable it explicitly with `-DAURORA_ENABLE_OPENMP=OFF` to force serial-only builds.

To install the current build into a local prefix, create the `TGZ` package,
and smoke-test the installed simulation and swarm CLIs plus downstream CMake
package metadata, run:

```sh
python3 scripts/verify_install_package.py build
```

Installed downstream projects can consume the library target with:

```cmake
find_package(AuroraPIC CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE AuroraPIC::aurorapic)
```

## Verify

The full smoke suite builds the project, validates milestone and
release-engineering artifacts, tests local gas import and the homogeneous
swarm CLI, validates the synthetic Hall comparison/preflight workflow, runs
the CTest regression executable, runs the standalone pusher
validation script (leapfrog plus Boris checks), runs the quantitative Landau
damping, two-stream, and 2D/3D Langmuir kinetic benchmarks, runs isolated CLI smoke tests
for the included 1D/2D/3D examples, and runs the install/package smoke test
for the installed CLIs, CPack `TGZ`, and downstream
`find_package(AuroraPIC CONFIG REQUIRED)` consumer:

```sh
scripts/verify.sh
```

Local verification defaults to one compiler job, one CTest job, and one
implicit OpenMP thread so it remains responsive on shared workstations.
Dedicated build hosts can opt in to higher limits with
`AURORA_BUILD_JOBS`, `AURORA_TEST_JOBS`, and `AURORA_OPENMP_THREADS`.

The example smoke tests copy each config to a temporary `test_output_aurorapic_verify/` directory, rewrite only `output_dir`, run `aurorapic_cli`, and assert the expected scalar, field, VTK, and particle-inspection files are structurally valid. Temporary smoke outputs are removed on success; set `KEEP_VERIFY_OUTPUTS=1` or pass `--keep-output` to `scripts/verify_examples.py` to retain them for debugging:

```sh
python3 scripts/verify_examples.py build/aurorapic_cli --keep-output
```

## 2D status

2D runs are selected by `dimension = 2` in the config file and are reachable from both the CLI/config loader and the C++ API. The 2D path provides:

- `Mesh2D`: rectangular node-centered mesh with independent periodic or Dirichlet field topology on each coordinate axis, plus side boundary tags/potentials for electrode-style Dirichlet axes. The legacy single `boundary` value still applies to both axes.
- `FieldSolver` automatically selects a fully periodic spectral solve, a
  mixed periodic/Dirichlet spectral-tridiagonal solve, or the fully Dirichlet
  SOR path. Mixed topologies transform the periodic axis with radix-2,
  mixed-radix, or Bluestein FFT algorithms and directly solve the Dirichlet
  axis mode by mode.
- `ImportedMesh2D`: validated topology/label model for tagged planar Gmsh v2 ASCII imports; bounded parsing rejects non-finite coordinates, duplicate entities, degenerate or non-convex cells, inconsistent orientation, non-manifold edges, and incomplete boundary closure. The model exposes cell area/centroid and boundary-length metrics without exposing solver code to `.msh` details.
- `UnstructuredMesh2D`: computational state over validated imported topology, with lumped nodal control areas, spatially accelerated triangle/quad point location, checkpoint-independent per-particle cell-location caches, conservative element-shape charge deposition, and nodal electric-field interpolation.
- `UnstructuredPoissonSolver2D` / `solve_unstructured_poisson`: triangle/quad finite-element stiffness assembly, strict physical-label Dirichlet/Neumann mapping, CSR storage, Jacobi-preconditioned conjugate gradients, convergence reporting, and nodal electric-field recovery. The reusable solver binds to one topology and caches quadrature, boundary contributions, the constrained sparse operator, and its Jacobi diagonal; the free functions remain as one-shot convenience APIs.
- `UnstructuredSimulation2D`: imported-geometry runtime with area-uniform domain seeding, optional rectangular bounds or named physical cell-region initialization, deterministic tagged-boundary particle injection, weight-aware secondary emission, cached electrostatic field solves, electrostatic/Boris particle advance, earliest-crossing geometry-aware absorbing or reflecting boundary policies, species-resolved impact flux diagnostics, transient/steady execution, topology-checked restart, particle samples, and unstructured `.vtu` field output.
- `Species2D`: explicit `Particle2D` storage with planar position and three-component velocity initialization, CIC deposition, full-velocity kinetic-energy accounting, and live-particle accounting.
- `Simulation2D`: deposit -> solve -> particle push/drift -> redeposit/resolve loop using the existing 2D Poisson solvers, with per-side particle boundary policies (`auto`, `absorbing`, `reflecting`, `periodic`). The default zero magnetic field uses the electrostatic leapfrog pusher; any nonzero uniform component or tabulated magnetic profile switches particles to the 2D3V Boris rotation/kick.
- `Diagnostics2D`: scalar time histories in `scalars.csv`, cumulative absorbed-particle counts by side, optional sampled particle CSV files, and opt-in resolved profile/moment/mode diagnostics for structured meshes.
- `write_legacy_vtk` / `write_vtk_xml`: structured-grid VTK writers for `rho`, `phi`, and electric-field vectors on `Mesh2D`.

When `vtk_output = true`, 2D runs write field snapshots under `output_dir` for ParaView or VisIt. `vtk_format` selects `legacy` (`fields_<step>.vtk`, the default), `xml`/`vts` (`fields_<step>.vts`), or `both`.

Imported 2D runs use `mesh = imported` and `mesh_file = <planar-v2-ascii.msh>`. Each physical boundary must have a matching `[boundary.<physical-name>]` section and an independent `absorbing` or `reflecting` particle policy. Field conditions are either `field = dirichlet` with a constant `potential`, or `field = neumann` with a constant outward `normal_derivative = dphi/dn`; omitting `field` preserves the legacy Dirichlet interpretation. At least one Dirichlet boundary is required to fix the electrostatic gauge. Since `E = -grad(phi)`, the specified Neumann value corresponds to `E dot n = -normal_derivative`. Species use `[species.<name>]` sections.

An imported species may set `initialization_region = <physical-cell-name>` to
load uniformly by area from only the Gmsh 2D physical group with that exact
name. Unknown or empty regions fail during construction. Region selection
supports both random and quiet-start loading and is mutually exclusive with
the rectangular `init_x_*`/`init_y_*` clip.

All dimensions accept `density_profile = uniform | gaussian | sinusoidal`.
Uniform is the backward-compatible default and rejects profile parameters.
Gaussian profiles require `profile_center_<axis>` and positive
`profile_scale_<axis>` for every active spatial axis, in physical/configured
coordinates. Sinusoidal profiles require `profile_amplitude` with
`|amplitude| <= 1`, at least one nonzero integer `profile_mode_<axis>`, and an
optional radian `profile_phase`; modes are evaluated over the active
initialization envelope. `max_profile_sampling_attempts` is a total
per-species work budget and must be at least the requested particle count.
Exhaustion fails the run instead of continuing with a biased or incomplete
population.

Optional `[source.<name>]` sections inject a fixed number of macro-particles per step from a tagged boundary. They require `species`, `boundary`, and positive `particles_per_step`; optional controls are `start_step`, exclusive `end_step` (`0` means unlimited), non-negative inward `normal_velocity`, signed `tangential_velocity`, signed `out_of_plane_velocity`, and non-negative `thermal_velocity`. Boundary segments are sampled in proportion to length. The normal thermal component is half-range inward, while the tangent and out-of-plane thermal components are signed Gaussians; positive tangent is defined as a clockwise rotation of the inward normal. The represented source rate is `particles_per_step * species_weight / dt`. Dead slots are reused before particle storage grows, and `max_particles_per_species` bounds both growth and checkpoint loading. See [`examples/imported_plasma_2d.cfg`](examples/imported_plasma_2d.cfg) and its companion mesh for the complete strict syntax. The CLI automatically dispatches these configs to the imported runtime.

Optional `[emission.<name>]` sections attach secondary emission to an absorbing tagged boundary. They require `boundary`, `incident_species`, `emitted_species`, and a positive physical `yield`. The expected emitted macro-particle count per impact is `yield * incident_weight / emitted_weight`; deterministic integer production and stochastic rounding preserve that expectation. `max_particles_per_impact` bounds weight-ratio expansion, while `normal_velocity`, `tangential_velocity`, `out_of_plane_velocity`, and `thermal_velocity` use the same convention as boundary sources. Impacts are sorted by species and particle ID before emission, making RNG use reproducible across serial and OpenMP execution. Diagnostics report cumulative emitted counts and, for every species/absorbing-boundary pair, cumulative macro-particles, represented physical particles, charge, full 3V incident kinetic energy, last-step physical-particle rate, and rate per boundary length.

For a larger imported example, run `examples/biased_probe_2d.cfg`. Its source geometry is `examples/biased_probe_2d.geo`; `scripts/generate_real_case_mesh.py` regenerates the checked-in v2 ASCII mesh with Gmsh 4.x. Imported-run startup reports node/cell/boundary counts plus minimum corner angle and maximum cell edge ratio so mesh provenance and basic quality are visible in logs.

Every run writes `initialization.csv` under `output_dir`. One row per species
records the IC schema version, generated/restart source, dimension, loading
and density-profile models, imported region when applicable, live macro-particle count and weight,
represented physical particles and charge, and realized mean and standard
deviation of each position and velocity component. Restart reports label the loading model
as `restart` and do not claim an original imported region.

All 2D runs write `scalars.csv` with:
```text
step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles,absorbed_left,absorbed_right,absorbed_bottom,absorbed_top,live_particles_<species>...
```

The absorbed-particle columns are cumulative counts of particles removed by absorbing particle boundary policies on each side. Reflecting and periodic particle boundaries keep particles live and do not increment these counters.

If `particle_output = true`, the run also writes sampled particle inspection files named `particles_<step>.csv` with:

```text
species_id,species,x,y,vx,vy,vz,alive
```

Particle-output controls:

- `particle_output`: enable/disable sampled particle CSV output; default `false`.
- `particle_output_interval`: particle CSV interval; `0` inherits `output_interval`.
- `particle_output_stride`: write every Nth particle per species traversal; must be positive.
- `particle_sample_count`: maximum rows across all species for each file; `0` writes all stride-selected particles.

Structured 2D can additionally emit resolved kinetic diagnostics:

```ini
resolved_diagnostics = true
resolved_diagnostic_interval = 10
resolved_diagnostic_start_step = 1000
resolved_profile_axis = x
resolved_mode_axis = y
resolved_max_mode = 32
```

The profile and mode axes must be distinct; the mode axis must be periodic,
and `resolved_max_mode` cannot exceed its mesh Nyquist limit. An interval of
zero inherits `output_interval`. `resolved_field_profiles.csv` contains the
transverse mean potential, electric field, and charge density at every profile
node. `resolved_species_profiles.csv` contains CIC-consistent represented
number density, three mean velocities, three thermal speeds, scalar SI
temperature in eV, and all three current-density components for each species.
Normalized runs write `nan` for temperature because they have no implicit
electron-volt scale.

`resolved_modes.csv` records the real, imaginary, and one-sided amplitude of
periodic-axis Fourier coefficients for charge density, both electric-field
components, and each species' number density and three currents. Particles
are accumulated onto the periodic diagnostic grid once per sample before the
mode transform, so cost scales with particles plus mesh points times requested
modes, rather than particles times modes.

The final `resolved_field_time_average.csv` and
`resolved_species_time_average.csv` use trapezoidal integration over the
actual resolved-sample times. Species averages integrate density-weighted
first and second velocity moments before deriving mean velocity and
temperature. Their start/end time, duration, and sample count make the
averaging window explicit. A restarted run creates a new output segment and
therefore a segment-local average; production campaign tooling must combine
segments using their recorded durations rather than averaging averages
equally.

## 3D status

3D runs are selected by `dimension = 3` in the config file and are reachable from both the CLI/config loader and the C++ API. The 3D path provides:

- `Mesh3D`: node-centered Cartesian mesh with periodic or grounded Dirichlet field boundary mode.
- `Species3D`: explicit `Particle3D` storage with 3D position/velocity initialization, trilinear CIC deposition, kinetic-energy accounting, and live-particle accounting.
- `Simulation3D`: deposit -> solve -> particle push/drift -> redeposit/resolve loop using the 3D Poisson solvers, with per-side particle boundary policies (`auto`, `absorbing`, `reflecting`, `periodic`). The default zero magnetic field uses the electrostatic leapfrog pusher; a nonzero uniform component or tabulated magnetic profile switches particles to the Boris rotation/kick.
- `Diagnostics3D`: scalar time histories in `scalars.csv`, cumulative absorbed-particle counts by side, and optional sampled particle CSV files.
- `write_legacy_vtk` / `write_vtk_xml`: structured-grid VTK writers for `rho`, `phi`, and electric-field vectors on `Mesh3D`.

When `vtk_output = true`, 3D runs write field snapshots under `output_dir` for ParaView or VisIt. `vtk_format` selects `legacy` (`fields_<step>.vtk`, the default), `xml`/`vts` (`fields_<step>.vts`), or `both`.

All 3D runs write `scalars.csv` with:

```text
step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles,absorbed_left,absorbed_right,absorbed_bottom,absorbed_top,absorbed_back,absorbed_front,live_particles_<species>...
```

If `particle_output = true`, sampled particle files are named `particles_<step>.csv` with:

```text
species_id,species,x,y,z,vx,vy,vz,alive
```

The 3D API also exposes `Vec3`, `Particle3D`, `deposit_charge_cic(Mesh3D&, ...)`, `FieldSolver::solve(Mesh3D&)`, and `interpolate_electric(const Mesh3D&, Vec3)` for lower-level multidimensional development.

## Configuration format

AuroraPIC uses a strict INI-like format. The optional `config_version` key currently supports only `config_version = 1`; omitted versions are treated as version 1 for backward compatibility, while unsupported future versions are rejected with a clear diagnostic. 1D configs may omit `dimension` or set `dimension = 1`:

```ini
config_version = 1
velocity_dimensions = 1  # 1 (default) or 3
nx = 128
length = 1.0
dt = 0.002
steps = 300
boundary = periodic        # periodic or dirichlet
# For a 1D Dirichlet electrode, phi is the static offset and the
# optional drive is offset + amplitude*sin(2*pi*frequency*time + phase).
# phi_left = 0
# phi_right = 0
# phi_right_amplitude = 450
# phi_right_frequency = 13.56e6
# phi_right_phase = 0
mode = transient           # transient or steady_state
output_interval = 25
output_dir = output/run

[species]
name = electrons
charge = -1
mass = 1
weight = 0.01             # macro-particle weight; required to be positive
# density = 100.0         # optional alternative: if weight is omitted,
#                         # weight = density * initialization_width / particles
particles = 10000
drift_velocity = 0
thermal_velocity = 0.1
# Optional 1D species timestep: push, boundary handling, and collisions
# occur every N base steps with an interval of N*dt. Default: 1.
# timestep_multiplier = 20
# In 1D3V, x is the only spatial coordinate and electrostatic E_x is the
# only field component; vy/vz are retained for distributions and collisions.
# drift_velocity_y = 0
# drift_velocity_z = 0
# Optional versioned initial-condition controls:
# initialization_version = 1
# loading = quiet_start    # random (default) or quiet_start
# thermal_velocity_x = 0.1 # overrides thermal_velocity
# thermal_velocity_y = 0.1 # accepted only with velocity_dimensions = 3
# thermal_velocity_z = 0.1 # accepted only with velocity_dimensions = 3
# density_profile = gaussian
# profile_center_x = 0.5    # physical coordinate
# profile_scale_x = 0.1     # positive Gaussian standard deviation
# max_profile_sampling_attempts = 1000000

[collisions]
enabled = false
model = bgk
frequency = 0.0
neutral_temperature_velocity = 0.0
```

The default `velocity_dimensions = 1` preserves the historical 1D1V state and
random-number sequence. Selecting `3` evolves the same one-dimensional
electrostatic position and `E_x` push while carrying `vy` and `vz` through
initialization, kinetic-energy diagnostics, BGK relaxation, isotropic
elastic/excitation MCC, and restart. External `.aps` initialization is
currently restricted to 1D1V; 1D3V external-state ingestion is a subsequent
interoperability slice.

In 1D, `timestep_multiplier = N` provides deterministic species subcycling.
The species advances on pre-step indices `0, N, 2N, ...`, using `N*dt` for
its leapfrog push and MCC interval; its deposited charge remains fixed on
intermediate field solves. Boundary handling follows each active push. The
step-zero phase matches eduPIC's ion schedule. The default `N = 1` preserves
the historical update path. Users must assess particle transit, acceleration,
and collision-frequency accuracy against the enlarged species timestep.

The sinusoidal drive keys are currently restricted to transient 1D Dirichlet
domains. Driven `steady_state` mode is rejected until convergence is evaluated
over complete RF cycles.
Frequency is cycles per simulation-time unit (hertz in SI), and phase is in
radians. Nonzero amplitude requires positive frequency. The post-drift field
solve evaluates the waveform at the new field time level, and restart
reconstructs its phase from the checkpoint time. All 1D runs write:

```text
step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles,phi_left,phi_right,live_particles_<species>...
```

The potential columns contain the values actually applied by a Dirichlet
field solve and are zero for periodic domains, which have no electrodes.
One trailing `live_particles_<species>` column is written for every configured
species, in configuration order.
For driven systems, particle-plus-field energy is not conserved because the
external electrode supplies or removes energy.

Every 1D run also writes `boundary_losses.csv`. For each species and electrode
it records cumulative absorbed macro-particle count, represented signed
charge, and full-velocity kinetic energy at impact. Differences between
successive rows give output-window current density and kinetic power per area
in SI runs. Every 1D run additionally writes `power_transfer.csv`, containing
cumulative species-resolved kinetic-energy change caused by the electric
particle push. A row difference divided by elapsed time is the mean electrical
power per area delivered to that species over the window (`W m^-2` in SI).
`counter_origin_step` is zero for a complete fresh/v7 history; a legacy
restart reports its restart step so partial diagnostic coverage cannot be
mistaken for whole-run data.

For tabulated MCC, select `model = null_collision`, name the target `species`,
set `neutral_density` and a conservative `max_frequency`, then add one or more
`[collision.<name>]` elastic/excitation sections. This legacy singular schema
remains backward compatible. A 1D3V discharge that needs independent collision
targets uses `[collisions.<model>]` and
`[collisions.<model>.channel.<channel>]`. Named models must target distinct
species. Their diagnostic columns are qualified as `<model>.<channel>`.
Ionization channels name `secondary_species` and `ion_species`; the target,
secondary, and ion must satisfy the equal-weight charge-conservation contract.
Products are capacity-preflighted against `max_particles_per_species` and are
staged until every collision target finishes, so newborn particles cannot
collide in their birth timestep. See
[`examples/mcc_ionization_1d.cfg`](examples/mcc_ionization_1d.cfg).
Excitation and ionization retain the legacy `heavy_target` transform by
default; v2 gas manifests and inline channels can opt into
`inelastic_transform = finite_mass_center_of_mass` when a positive neutral
mass is available. This applies an explicit center-of-mass drift and mass
factor to post-threshold projectile velocities without dynamically tracking
the residual neutral/ion recoil product.
Imported 2D3V MCC also
requires `gas`, positive `neutral_mass`, non-negative `neutral_temperature`,
and `data_provenance`, and can use ionization sections naming secondary and
ion product species, attachment sections naming a negative-ion product, or
resonant `charge_exchange` channels. Alternatively,
`gas_data_file` loads a reusable,
versioned `.gas` manifest containing gas identity, mass, dataset/version,
provenance, citation, retrieval date, license, and channel tables while the
simulation retains operating conditions and product-species mappings. See
[`examples/mcc_relaxation.cfg`](examples/mcc_relaxation.cfg),
[`examples/mcc_ionization_1d.cfg`](examples/mcc_ionization_1d.cfg),
[`examples/imported_mcc_2d.cfg`](examples/imported_mcc_2d.cfg), and the
[`examples/imported_ionization_2d.cfg`](examples/imported_ionization_2d.cfg)
and
[`examples/imported_attachment_2d.cfg`](examples/imported_attachment_2d.cfg)
and
[`examples/imported_charge_exchange_2d.cfg`](examples/imported_charge_exchange_2d.cfg).
The collision contract documents the reactive kinematic and macro-weight
constraints. All bundled MCC datasets are synthetic validation
inputs, not gas-property data. The
[local real-gas workflow](docs/gas-data-workflow.md) converts a user-supplied
LXCat/BOLSIG+ export, records its SHA-256 and provenance, audits its collision
envelope, and produces a unit-safe versioned manifest without downloading or
vendoring the source dataset. The resulting package can be checked outside a
device geometry with the
[homogeneous electron-swarm runner](docs/swarm-validation.md), which scans
E/N using the production collision kernel and emits traceable transport and
channel-rate diagnostics. Positive-temperature SI runs derive the neutral
Maxwellian component speed from gas mass and temperature and enforce a bounded
thermal collision-rate majorant. Its local comparison tool maps user-supplied
measured or evaluated coefficients to those outputs and produces a hashed,
uncertainty-aware acceptance report without bundling the reference data.
The serialized campaign runner applies that contract to every resolution,
checks physics identity and uncertainty-aware numerical convergence, and
keeps all runs single-threaded and sequential by default.
The Turner helium CCP workflow additionally provides a hash-audited,
non-launching production Case 1 preparer. It locks the exact grid, initial
population, RF duration, restart-safe averaging window, and guarded collision
majorants while keeping the publisher tables under ignored `tmp/`; see
[`docs/ccp-validation.md`](docs/ccp-validation.md).
The Turner-only normalization audit pins the paper's required 2006 CODATA
eV-to-joule conversion rather than silently inheriting current SI constants.
An explicitly acknowledged bounded qualifier measures the exact initial
population for only a few serial steps under hard work and timeout limits
before any longer CCP execution is considered.
A separate pinned eduPIC argon workflow can generate a complete one-period
startup screen with `--steps 4000 --startup-diagnostics`; its exact population,
global energy, spatial collision-energy, and RF-phase accounting remain
explicitly below the cross-code validation claim boundary.
The checksum-gated `run_aurorapic_edupic_pilot.py` continuation advances that
state no farther than cycle 4, one low-priority process per cycle, and stops on
hard population-growth, particle-cap, field, closure, memory, or timeout gates.
`extend_aurorapic_edupic_horizon.py` continues only in hash-chained four-cycle
blocks. A safe pilot report starts the chain and each safe horizon report can
start the next block; unsafe or hash-mismatched reports are rejected. The tool
keeps its internal stationarity result distinct from safe execution and from
any external validation claim.
The original cycle-16 ceiling can be exceeded only with the repository's
checksum-approved comparison-readiness rule. That rule permits one block per
invocation through cycle 64 and requires two consecutive stationary blocks
before a separate measurement campaign becomes eligible.
The completed comparison-readiness horizon reached cycle 64 with every hard
gate passing, but population still grew at 1.159% per cycle and ionization
substantially exceeded wall loss. Measurement therefore remains locked; the
serial MCC hot path has since been accelerated with bitwise-equivalent cycle
64 evidence. Restart-safe, count-preserving species/electrode wall-impact
spectra have also passed a real-state cycle-64 pilot. The next prerequisite is
a prospectively declared production-scale equilibration campaign.
A separate checkpointed startup ladder advances one RF cycle at a time,
retaining species/collision balances and phase-matched field metrics while
hard-limiting each local horizon.
Hash-chained horizon reports feed a separate pre-benchmark stationarity
screen, so the published Turner density comparison remains unavailable until
population, ionization, field, and energy trends have settled.
The runner also offers an explicit bounded branching mode: ionization and
attachment respectively increase and decrease represented electron weight
while systematic resampling holds computational population fixed, enabling
temporal avalanche growth, rate-balance effective ionization, and
growth-over-flux-drift Townsend diagnostics without unbounded host load.
An optional finite-distance history experiment independently fits the slope
of steady signed electron flux across interior planes and writes the complete
profile, block uncertainty, fit R², and bounded-work evidence.
Elastic 3V channels may additionally provide a validated energy-dependent
mean-cosine table for Henyey-Greenstein anisotropic scattering; isotropic
behavior remains the explicit default. Elastic and resonant charge-exchange
tables can explicitly use projectile or center-of-mass lookup energy, and a
3V elastic channel can select exact backward center-of-mass scattering. These
generic controls cover the two-component He+-He collision law required by the
Turner CCP benchmark without changing the default collision contract.

2D configs must set `dimension = 2` and use `nx`/`ny`, `length_x`/`length_y`, 2D velocity keys, and 2D initialization bounds. `boundary = dirichlet` may also provide side electrode potentials (`phi_left`, `phi_right`, `phi_bottom`, `phi_top`) and side tags (`boundary_left_tag`, `boundary_right_tag`, `boundary_bottom_tag`, `boundary_top_tag`):

```ini
config_version = 1
dimension = 2
nx = 64
ny = 64
length_x = 1.0
length_y = 1.0
dt = 0.002
steps = 100
mode = transient
boundary = dirichlet
# Optional 2D overrides; omission preserves the global boundary on each axis.
boundary_x = dirichlet
boundary_y = periodic
phi_left = -5.0
phi_right = 5.0
phi_bottom = 0.0
phi_top = 0.0
boundary_left_tag = cathode
boundary_right_tag = anode
boundary_bottom_tag = grounded_wall
boundary_top_tag = grounded_wall
# Particle boundaries are independent from field/electrode potentials.
# particle_boundary sets the default for all sides; side-specific keys override it.
# auto resolves independently: periodic on a periodic coordinate axis and
# absorbing on a Dirichlet coordinate axis.
particle_boundary = absorbing
particle_boundary_right = reflecting
# Optional uniform B field. Any nonzero component uses the 2D3V Boris pusher.
magnetic_field_x = 0.0
magnetic_field_y = 0.0
magnetic_field_z = 0.0
output_interval = 10
output_dir = output/electrode_2d
checkpoint_output = true
checkpoint_interval = 25
# restart_path = output/electrode_2d/checkpoint_50.apc
vtk_output = true
# Optional field snapshot format: legacy (default), xml/vts, or both.
vtk_format = both
particle_output = true
particle_output_interval = 10
particle_output_stride = 5
particle_sample_count = 200

[species.electrons]
charge = -1
mass = 1
density = 100
particles = 10000
drift_velocity_x = 0.1
drift_velocity_y = 0.0
drift_velocity_z = 0.0
thermal_velocity = 0.02
# Optional anisotropic/quiet-start initialization. Per-axis values override
# thermal_velocity; omitted axes retain the scalar fallback.
initialization_version = 1
loading = quiet_start
thermal_velocity_x = 0.02
thermal_velocity_y = 0.01
thermal_velocity_z = 0.03
init_x_min = 0.0
init_x_max = 1.0
init_y_min = 0.0
init_y_max = 1.0
```

3D configs must set `dimension = 3` and use `nx`/`ny`/`nz`, `length_x`/`length_y`/`length_z`, 3D velocity keys, and 3D initialization bounds:

```ini
config_version = 1
dimension = 3
nx = 32
ny = 32
nz = 32
length_x = 1.0
length_y = 1.0
length_z = 1.0
dt = 0.001
steps = 100
mode = transient
boundary = periodic
particle_boundary = auto
# Optional uniform B field. Any nonzero component uses the Boris pusher.
magnetic_field_x = 0.0
magnetic_field_y = 0.0
magnetic_field_z = 0.0
output_interval = 10
output_dir = output/plasma_3d
checkpoint_output = true
checkpoint_interval = 25
# restart_path = output/plasma_3d/checkpoint_50.apc
vtk_output = true
# Optional field snapshot format: legacy (default), xml/vts, or both.
vtk_format = both
particle_output = true
particle_output_interval = 10
particle_output_stride = 10
particle_sample_count = 200

[species.electrons]
charge = -1
mass = 1
density = 100
particles = 10000
drift_velocity_x = 0.1
drift_velocity_y = 0.0
drift_velocity_z = 0.0
thermal_velocity = 0.02
initialization_version = 1
loading = quiet_start
thermal_velocity_x = 0.02
thermal_velocity_y = 0.01
thermal_velocity_z = 0.03
init_x_min = 0.0
init_x_max = 1.0
init_y_min = 0.0
init_y_max = 1.0
init_z_min = 0.0
init_z_max = 1.0
```

## Checkpoint/restart controls

All 1D, 2D, and 3D runs support text `.apc` checkpoints for deterministic restart and regression debugging. Structured 2D checkpoint v10 records extrusion depth, volumetric and current-regulated source state, controller/correction modes, reverse-demand statistics and distribution moments, species-resolved boundary losses, potential-reference configuration, unit metadata, RNG state, and full 3V particle state. Structured 2D v1-v9 remain readable when their controller contract is compatible; reverse-demand totals begin at the restart step for pre-v9 inputs and distribution statistics begin there for pre-v10 inputs. Current 1D v16 records whether spatial diagnostics sample before or after collisions, extending the v15 collision-event-aware format. A changed sampling order is rejected as part of the averaging contract; older checkpoints imply the historical `post_collision` order. A changed spectrum contract is rejected; enabling spectra from a pre-v14 checkpoint starts an origin-labeled window at the restart step. A changed species schedule is rejected; pre-v13 checkpoints require every multiplier to remain 1. By default a changed averaging contract is rejected. Explicit `spatial_average_reset_on_restart = true` discards stored sums and starts a new window, which must begin after the checkpoint step so no sample is missed. Older checkpoints retain the state they support and begin newer origin-labeled counters at the restart step; pre-v11 collision-localization and pre-v12 EEDF windows must be reset because those checkpoint versions lack the corresponding state. Structured 3D v2 and imported 2D v6 retain their documented compatibility rules. Imported checkpoints additionally fingerprint mesh topology, coordinates, and tags.

- `checkpoint_output`: enable/disable checkpoint writes during `run()`; default `false`.
- `checkpoint_interval`: checkpoint interval in steps; `0` inherits `output_interval` when `checkpoint_output = true`.
- `checkpoint_path`: optional fixed checkpoint file path. If omitted, checkpoints are written as `output_dir/checkpoint_<step>.apc`. If provided, each checkpoint write updates that same path.
- `restart_path`: optional checkpoint file to load before the run loop starts. The run resumes from the checkpoint step/time and continues until the configured `steps`/`max_steps` limit.
- `initial_state_path`: optional versioned `.aps` particle initial state, resolved relative to the config file. It starts at step/time zero, validates dimension, units, species, counts, finite records, and geometry, and rebuilds the field and half-step state. It cannot be combined with `restart_path`.
- `initial_state_signature`: optional decimal or `0x`-prefixed canonical 64-bit semantic signature for `initial_state_path`. A mismatch rejects changed input before loading. External runs write the resolved source and realized/expected signatures to `initial_state_metadata.txt`.

## Run modes and termination

All 1D, 2D, and 3D simulations support the same run controls:

- `mode = transient`: run until `steps` is reached.
- `mode = steady_state`: run until convergence is detected or `max_steps` is reached.
- `steady_window`: number of emitted diagnostic samples in each adjacent comparison window.
- `steady_tolerance`: maximum relative change between the mean total energies of those windows.
- `max_steps`: hard safety cap for a steady-state run.

Steady-state convergence is an engineering termination criterion, not by itself proof of physical equilibrium. It is evaluated only when scalar diagnostics are sampled, so `output_interval`, `steady_window`, and `dt` jointly determine the physical duration represented by a convergence window. Run summaries and CLI output distinguish a converged steady run from one that exhausted `max_steps`. When checkpoint output is enabled, convergence forces a final checkpoint even when the regular checkpoint interval has not been reached.

Periodically driven 1D cases use the stronger opt-in RF-cycle controller:

```ini
periodic_convergence = true
periodic_convergence_reset_on_restart = false
periodic_convergence_rf_frequency = 13.56e6
periodic_convergence_cycles_per_block = 32
periodic_convergence_minimum_blocks = 16
periodic_convergence_minimum_effective_blocks = 8
periodic_convergence_maximum_absolute_projected_fractional_drift = 0.01
periodic_convergence_maximum_absolute_split_half_fractional_change = 0.01
periodic_convergence_maximum_relative_standard_error = 0.01
```

The RF period must contain an integer number of timesteps and its frequency
must match every active sinusoidal electrode. AuroraPIC samples represented
species populations and total energy at a consistent RF phase, groups them
into complete-cycle blocks, and requires every observable to pass nominal
block-count, autocorrelation-adjusted effective-sample, drift, split-half, and
standard-error gates. This controller permits `mode = steady_state` for driven
1D cases. Its complete history and thresholds are stored in checkpoint v25;
the CSV decision record is written to the configured output directory.
Strict restart restoration is the default. Set
`periodic_convergence_reset_on_restart = true` only to begin a declared fresh
statistical epoch from an older checkpoint or to discard an existing epoch.
The checkpoint step must lie exactly on the configured RF phase; AuroraPIC
rejects between-phase resets and never treats pre-checkpoint samples as part of
the new convergence decision.

## Runtime controls

All 1D, 2D, and 3D configs accept runtime controls for the M4 scaling interface:

- `runtime_backend`: `serial` (default) or `openmp`; `mpi` and `gpu` are reserved placeholders and are rejected until those backends are implemented.
- `runtime_threads`: positive thread count. Serial requires `1`; OpenMP may use values greater than `1` only when AuroraPIC was built with OpenMP support.

Structured particle initialization/synchronization loops and the 1D particle advance use `RuntimePolicy` through deterministic static scheduling. The imported 2D runtime additionally parallelizes particle push, geometry-aware boundary handling, synchronization, and charge deposition. Imported deposition uses per-worker nodal buffers followed by an ordered reduction, avoiding shared-node races while preserving charge conservation. Each imported particle retains its last valid element and shape coordinates in a runtime sidecar: interpolation and deposition first validate that cell, then use the spatial index only after a crossing or cache miss. The regression suite compares serial and OpenMP structured and imported 2D runs when OpenMP is available, while keeping serial behavior as the portability baseline.

2D/3D particle-boundary controls:

- `particle_boundary`: default particle policy for all sides; one of `auto`, `absorbing`, `reflecting`, or `periodic`.
- `particle_boundary_left`, `particle_boundary_right`, `particle_boundary_bottom`, `particle_boundary_top`: per-side overrides for the default policy; 3D also supports `particle_boundary_back` and `particle_boundary_front`.
- Structured 2D may set `boundary_x` and `boundary_y` independently to `periodic` or `dirichlet`. Omitted axis keys inherit the legacy global `boundary`. Endpoint-free spacing, node areas, CIC deposit/gather, Poisson neighbors, and electric gradients follow the selected axis topology.
- `auto`: in structured 2D resolves per axis to `periodic` or `absorbing`; legacy 1D/3D behavior continues to follow the global field boundary.
- `absorbing`: removes particles that leave the domain and increments the corresponding `absorbed_*` scalar diagnostic.
- `reflecting`: mirrors escaped particles back into the domain and reverses the normal velocity component.
- `periodic`: wraps escaped particles across that coordinate direction.

2D/3D magnetic-field controls:

- Planar structured/imported 2D3V and structured 3D support uniform `magnetic_field_x`, `magnetic_field_y`, and `magnetic_field_z`. They default to `0.0`; any nonzero component activates the full three-velocity Boris pusher.
- Alternatively, `magnetic_field_profile_file` and `magnetic_field_profile_axis = x | y | z` select a strict whitespace table with `coordinate Bx By Bz` columns. Coordinates must be finite and strictly increasing, values are linearly interpolated at each particle position, and the table must cover the full domain on its selected axis. Relative paths resolve from the configuration file.
- Uniform magnetic components and a profile are mutually exclusive. Extrapolation is rejected instead of silently clamped. [`examples/hall_field_profile_smoke.cfg`](examples/hall_field_profile_smoke.cfg) exercises a published Hall-benchmark profile as one consumer of this generic interface; it is explicitly not a discharge result.
- The current field solve remains electrostatic Poisson. These controls add prescribed magnetic rotation to particle pushes, not a self-consistent electromagnetic field update or an arbitrary 2D/3D field-map importer.

Structured 2D configs may also declare bounded profiled volumetric pair
sources:

```ini
[source.channel_pair_seed]
first_species = electrons
second_species = ions
peak_volumetric_pair_rate = 5.23e23
start_step = 0
end_step = 0
x_min = 0.005
x_max = 0.020
y_min = 0
y_max = 0.0128
first_drift_velocity_y = 10000
first_thermal_velocity = 0
second_thermal_velocity = 0
density_profile = sinusoidal
profile_amplitude = -1
profile_mode_x = 1
```

Structured 2D uses a positive global `out_of_plane_depth` (default `1`) as its extrusion measure. Density-derived macro weights include this depth, deposition divides total macro weight by it to recover volume charge density, and field/particle energies are totals over the extruded volume.

Exactly one rate is required. `pairs_per_step` requests a fixed integer macro-pair count. `represented_pair_rate` specifies the total represented rate over the source region. `peak_volumetric_pair_rate` derives that total as peak rate times the analytic normalized-profile area times `out_of_plane_depth`. Fractional macro-pairs carry deterministically across steps.

Every created pair shares one sampled position and has independent isotropic Gaussian velocity distributions around the two configured drift vectors. The optional `density_profile = uniform | gaussian | sinusoidal` and associated `profile_*` keys use the same normalized rejection-sampling contract as species initialization. The profile determines relative position probability inside the source rectangle. It does not alter an explicitly configured total `represented_pair_rate`, but its analytic integral is part of the conversion from `peak_volumetric_pair_rate` to total rate. For example, sinusoidal amplitude `-1`, x mode `1`, and zero phase produce a centered `sin²(pi*x_normalized)` envelope. `end_step` is exclusive and zero means unlimited.

The species must have opposite equal charge and equal macro weight, making every event exactly charge balanced. `sources.csv` reports cumulative macro and represented pairs, fractional remainder, injected energy, effective profile area, extrusion depth, peak rate, and derived total rate. Structured 2D checkpoint v10 preserves this contract. The machine-validated [`examples/hall_landmark_axial_azimuthal.case`](examples/hall_landmark_axial_azimuthal.case) pins the public benchmark derivation and reduced-run limits; `scripts/validate_hall_case.py` rejects drift. These remain prescribed sources, without neutral depletion, reaction closure, recoil, or arbitrary tabulated profiles.

Structured 2D can additionally regulate a single-species source from represented charge lost at an absorbing boundary:

```ini
current_source_species = electrons
current_source_monitor_boundary = left
current_source_emission_boundary = right
current_source_emission_inset = 0.001
current_source_temperature_ev = 10
potential_reference_axis = x
potential_reference_coordinate = 0.024
potential_reference_target = 0
```

The controller sums `absorbed_count * charge * macro_weight` over every species at the monitored boundary. It converts the newly observed charge to the emitted species’ macro weight and creates the non-negative integer part at a uniformly sampled emission plane. `current_source_control_mode = cumulative` retains signed fractional surplus or debt for backward-compatible generic regulation; `timestep_local` clears reverse demand each step and retains only positive fractional remainder, matching the LANDMARK electron-minus-ion loss rule for equal macro weights. `current_source.csv` distinguishes raw cumulative mismatch, unserved reverse charge, and controllable remainder. It also reports one-, two-, and multi-macroparticle reverse-event bins, demand mean/RMS, and negative/positive monitored charge accumulated specifically on reverse timesteps. Structured 2D runs additionally emit `boundary_flux.csv`: for every species, boundary, and scalar-output window it records macro and represented losses, represented charge, particle rate, and charge rate. The initial zero-duration rows establish a restart-safe baseline and never masquerade as a physical rate.

The optional potential reference sets the interpolated transverse mean potential at the configured x or y coordinate to its target. `potential_reference_correction = gauge` applies the backward-compatible spatially constant shift and leaves the electric field unchanged. `affine` preserves the zero-coordinate electrode, applies a linear potential correction, and adjusts the corresponding electric-field component analytically. The Hall campaign explicitly selects timestep-local current control and affine potential correction.

For SI structured-2D configurations, species may use `temperature_ev` instead
of scalar or component `thermal_velocity` keys. Pair sources similarly accept
`first_temperature_ev` and `second_temperature_ev`, while the regulated source
accepts `current_source_temperature_ev`. AuroraPIC converts each value to the
one-component Maxwellian standard deviation `sqrt(e * temperature_ev / mass)`.
Temperature and velocity forms are mutually exclusive, and eV inputs are
rejected in normalized mode because normalized particle mass has no implicit
kilogram scale.

The parser is intentionally strict: unsupported `config_version` or species `initialization_version` values, unknown sections/keys, invalid initial loading models, density profiles, sampling budgets, component thermal velocities, or external particle-state metadata/records, invalid unit systems or relative permittivities, invalid enum values, invalid particle-boundary values, invalid booleans, non-finite numbers, non-positive `dt`/`output_interval`, invalid checkpoint intervals when checkpoint output is enabled, invalid electrode drive amplitude/frequency/phase combinations, non-positive particle limits/output strides, malformed collision channels/tables or unsafe collision-rate bounds, empty 2D boundary tags, nonzero electrode potentials on a periodic coordinate axis, non-finite magnetic-field values, invalid source schedules/velocities/references, invalid emission yields/limits/references, and invalid species initialization intervals are rejected instead of silently falling back to defaults. Emission rules must target an absorbing boundary, and unsafe macro-particle expansion is rejected during construction. For structured species definitions, provide either an explicit positive `weight` or omit `weight` and provide a positive `density`; the loader converts density to macro-particle weight over the configured initialization interval or area. With a nonuniform profile, this density fixes the total represented population (equivalently the volume-average density); the profile fixes its normalized relative spatial shape.

Optional global initialization gates reject inconsistent generated or restarted
particle states before time integration:

```ini
initialization_max_relative_charge_imbalance = 1e-10
initialization_max_relative_current_imbalance = 1e-8
initialization_max_relative_pair_imbalance = 1e-10
initialization_charge_pairs = electrons:ions
```

Charge imbalance is normalized by total absolute represented charge. Current
imbalance uses the norm of charge-weighted species mean velocities, normalized
by the corresponding absolute contributions. Named pairs must have opposite
represented-charge signs and matching charge magnitudes within the pair
tolerance. Tolerances are dimensionless values in `[0, 1]`; omitted gates are
disabled. Every run writes `initialization_acceptance.csv`, and a failed gate
leaves both initialization audit files before aborting.

## Performance and validation envelope

The verified smoke/performance envelope is documented in `docs/performance-envelope.md`, and the quantitative Landau-damping, two-stream, and 2D/3D Langmuir cases in `docs/kinetic-validation.md`. The staged Hall-thruster targets and strict production-only resource policy are documented in `docs/hall-thruster-validation.md`. Imported scalar diagnostics expose cumulative particle, deposition, and field-solve timings plus location-cache hits and spatial searches, and `scripts/benchmark_unstructured.py` reports repeat medians for a chosen imported config. The checked-in examples prove that the documented 1D/2D/3D CLI paths, diagnostics, VTK output, particle samples, prescribed uniform/profiled magnetic-field Boris activation, and checkpoint-style text outputs remain structurally valid at small CI-friendly sizes. The kinetic cases additionally verify published damped and unstable collisionless responses plus analytic 2D/3D directional plasma oscillations; neither suite proves convergence for arbitrary plasma regimes or validates a real device. Before using larger runs, document resolution, timestep, particles-per-cell/noise, output cadence, boundary model, and convergence checks against mesh/time/particle refinements.

This is a serious first version, not a final plasma platform. Key known gaps are: no MPI/GPU backend yet, OpenMP remains a shared-memory particle-path implementation rather than a domain-decomposed whole-solver model, external initial-state loading does not yet have a chunked openPMD/HDF5 high-volume backend, MCC thermal neutrals have fixed temperature and zero bulk flow, excitation and ionization do not dynamically track the residual target recoil product (although an explicit finite-mass projectile transform is available), charge exchange is limited to the resonant mass-matched case, and there is no neutral depletion, gas heating, or general reaction network; prescribed magnetic fields are uniform or one-coordinate profiles only (no arbitrary multidimensional map or self-consistent electromagnetic field solve yet), imported electrostatic field conditions are limited to label-wise constant Dirichlet/Neumann data, and the imported runtime has not been performance-qualified on production-scale meshes. No authoritative He/Ar/Kr/Xe cross-section set is bundled yet. High-volume particle dumps are intentionally deferred to an openPMD/HDF5-style format in a later phase; current `.aps` initial states, text checkpoints, and particle CSV output are for preprocessing, restart, inspection, and regression/debug workflows.
