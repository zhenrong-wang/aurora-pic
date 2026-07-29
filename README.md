# AuroraPIC

AuroraPIC is a C++20 starting point for scientific plasma dynamics simulation. The current codebase implements electrostatic `1D1V` and `1D3V`, planar `2D3V` (structured and imported geometry), and structured `3D3V` Particle-in-Cell (PIC) paths with configurable species, periodic or Dirichlet boundaries, transient fixed-step and steady-state convergence modes, scalar diagnostics, and text checkpoint/restart files. The 1D baseline provides the historical BGK relaxation model plus tabulated null-collision MCC; 1D3V supports simultaneous named collision targets, isotropic elastic/excitation scattering, and bounded charge-conservative electron-impact ionization products. Imported 2D3V runs additionally support bounded-Maxwellian finite-mass neutrals, electron attachment, and resonant ion-neutral charge exchange. The multidimensional paths provide prescribed uniform magnetic-field Boris pushes, VTK field output, side-specific particle boundaries, and optional particle inspection CSVs. A separate homogeneous electron-swarm runner scans reduced electric field with the same three-velocity collision kernel before a gas package is used in a device geometry.

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
public WarpX reference output and NASA HERMeS measurements, separates
code-to-code verification from experimental validation, and records the
physics, diagnostics, parallelism, and resource gates still required.

The dimensional contract is defined in [`docs/units.md`](docs/units.md). Configurations may select `units = normalized` or `units = si` plus a positive homogeneous `relative_permittivity`. Legacy omission remains normalized; maintained examples are explicit. SI reduced-dimensional runs report per-unit omitted measure (`J/m²` in 1D and `J/m` in planar 2D).

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
swarm CLI, runs the CTest regression executable, runs the standalone pusher
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

- `Mesh2D`: rectangular node-centered mesh with periodic or Dirichlet field boundary mode plus side boundary tags/potentials for electrode-style Dirichlet domains.
- `ImportedMesh2D`: validated topology/label model for tagged planar Gmsh v2 ASCII imports; bounded parsing rejects non-finite coordinates, duplicate entities, degenerate or non-convex cells, inconsistent orientation, non-manifold edges, and incomplete boundary closure. The model exposes cell area/centroid and boundary-length metrics without exposing solver code to `.msh` details.
- `UnstructuredMesh2D`: computational state over validated imported topology, with lumped nodal control areas, spatially accelerated triangle/quad point location, checkpoint-independent per-particle cell-location caches, conservative element-shape charge deposition, and nodal electric-field interpolation.
- `UnstructuredPoissonSolver2D` / `solve_unstructured_poisson`: triangle/quad finite-element stiffness assembly, strict physical-label Dirichlet/Neumann mapping, CSR storage, Jacobi-preconditioned conjugate gradients, convergence reporting, and nodal electric-field recovery. The reusable solver binds to one topology and caches quadrature, boundary contributions, the constrained sparse operator, and its Jacobi diagonal; the free functions remain as one-shot convenience APIs.
- `UnstructuredSimulation2D`: imported-geometry runtime with area-uniform domain seeding, optional rectangular bounds or named physical cell-region initialization, deterministic tagged-boundary particle injection, weight-aware secondary emission, cached electrostatic field solves, electrostatic/Boris particle advance, earliest-crossing geometry-aware absorbing or reflecting boundary policies, species-resolved impact flux diagnostics, transient/steady execution, topology-checked restart, particle samples, and unstructured `.vtu` field output.
- `Species2D`: explicit `Particle2D` storage with planar position and three-component velocity initialization, CIC deposition, full-velocity kinetic-energy accounting, and live-particle accounting.
- `Simulation2D`: deposit -> solve -> particle push/drift -> redeposit/resolve loop using the existing 2D Poisson solvers, with per-side particle boundary policies (`auto`, `absorbing`, `reflecting`, `periodic`). The default zero magnetic field uses the electrostatic leapfrog pusher; any nonzero `magnetic_field_x/y/z` component switches particles to the 2D3V Boris rotation/kick.
- `Diagnostics2D`: scalar time histories in `scalars.csv`, cumulative absorbed-particle counts by side, and optional sampled particle CSV files.
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

## 3D status

3D runs are selected by `dimension = 3` in the config file and are reachable from both the CLI/config loader and the C++ API. The 3D path provides:

- `Mesh3D`: node-centered Cartesian mesh with periodic or grounded Dirichlet field boundary mode.
- `Species3D`: explicit `Particle3D` storage with 3D position/velocity initialization, trilinear CIC deposition, kinetic-energy accounting, and live-particle accounting.
- `Simulation3D`: deposit -> solve -> particle push/drift -> redeposit/resolve loop using the 3D Poisson solvers, with per-side particle boundary policies (`auto`, `absorbing`, `reflecting`, `periodic`). The default zero magnetic field uses the electrostatic leapfrog pusher; nonzero `magnetic_field_x/y/z` switches particles to the Boris rotation/kick.
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

The sinusoidal drive keys are currently restricted to transient 1D Dirichlet
domains. Driven `steady_state` mode is rejected until convergence is evaluated
over complete RF cycles.
Frequency is cycles per simulation-time unit (hertz in SI), and phase is in
radians. Nonzero amplitude requires positive frequency. The post-drift field
solve evaluates the waveform at the new field time level, and restart
reconstructs its phase from the checkpoint time. All 1D runs write:

```text
step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles,phi_left,phi_right
```

The potential columns contain the values actually applied by a Dirichlet
field solve and are zero for periodic domains, which have no electrodes.
For driven systems, particle-plus-field energy is not conserved because the
external electrode supplies or removes energy.

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
behavior remains the explicit default.

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
# auto maps to periodic when boundary = periodic, and absorbing when boundary = dirichlet.
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

All 1D, 2D, and 3D runs support a text `.apc` checkpoint format intended for deterministic restart and regression debugging. Checkpoints include the simulation dimension, unit contract, step/time, RNG engine state, per-species particle positions/velocities/leapfrog half-step state, live flags, and 2D/3D absorbed-particle counters. Imported 2D checkpoints additionally store a deterministic topology/coordinate/tag fingerprint and refuse restart against a different mesh. The current 1D v4 format records the configured velocity dimensionality, all active velocity components, collision fingerprints, and collision counters; legacy 1D v1-v3 remain readable only by 1D1V runs. Structured 2D v3, structured 3D v2, and imported v6 record unit metadata. The current 2D formats preserve all three velocity components; legacy 2D checkpoints load with zero out-of-plane velocity. Imported v6 records source/emission and boundary-flux state plus optional MCC identity, gas metadata, effective-table fingerprint, RNG state, and collision counters. Imported v1–v5 remain readable for collision-free runs, while imported MCC restart requires v6. Imported v1–v3 and structured v1 checkpoints remain readable only with the historical normalized unit contract; 1D v1/v2 cannot restart null-collision MCC. Physical flux history begins at zero when loading imported v1 or v2.

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

## Runtime controls

All 1D, 2D, and 3D configs accept runtime controls for the M4 scaling interface:

- `runtime_backend`: `serial` (default) or `openmp`; `mpi` and `gpu` are reserved placeholders and are rejected until those backends are implemented.
- `runtime_threads`: positive thread count. Serial requires `1`; OpenMP may use values greater than `1` only when AuroraPIC was built with OpenMP support.

Structured particle initialization/synchronization loops and the 1D particle advance use `RuntimePolicy` through deterministic static scheduling. The imported 2D runtime additionally parallelizes particle push, geometry-aware boundary handling, synchronization, and charge deposition. Imported deposition uses per-worker nodal buffers followed by an ordered reduction, avoiding shared-node races while preserving charge conservation. Each imported particle retains its last valid element and shape coordinates in a runtime sidecar: interpolation and deposition first validate that cell, then use the spatial index only after a crossing or cache miss. The regression suite compares serial and OpenMP structured and imported 2D runs when OpenMP is available, while keeping serial behavior as the portability baseline.

2D/3D particle-boundary controls:

- `particle_boundary`: default particle policy for all sides; one of `auto`, `absorbing`, `reflecting`, or `periodic`.
- `particle_boundary_left`, `particle_boundary_right`, `particle_boundary_bottom`, `particle_boundary_top`: per-side overrides for the default policy; 3D also supports `particle_boundary_back` and `particle_boundary_front`.
- `auto`: resolves to `periodic` for `boundary = periodic` field solves and `absorbing` for `boundary = dirichlet` field solves.
- `absorbing`: removes particles that leave the domain and increments the corresponding `absorbed_*` scalar diagnostic.
- `reflecting`: mirrors escaped particles back into the domain and reverses the normal velocity component.
- `periodic`: wraps escaped particles across that coordinate direction.

2D/3D magnetic-field controls:

- Planar 2D3V supports uniform `magnetic_field_x`, `magnetic_field_y`, and `magnetic_field_z`. They default to `0.0`; any nonzero component activates the full three-velocity Boris pusher while position and electrostatic fields remain planar.
- 3D supports uniform `magnetic_field_x`, `magnetic_field_y`, and `magnetic_field_z`. They default to `0.0`; any nonzero component activates the Boris pusher for 3D particles.
- Magnetic-field values must be finite. The current field solve remains electrostatic Poisson; these controls add prescribed uniform magnetic rotation to particle pushes, not a self-consistent electromagnetic field update.

The parser is intentionally strict: unsupported `config_version` or species `initialization_version` values, unknown sections/keys, invalid initial loading models, density profiles, sampling budgets, component thermal velocities, or external particle-state metadata/records, invalid unit systems or relative permittivities, invalid enum values, invalid particle-boundary values, invalid booleans, non-finite numbers, non-positive `dt`/`output_interval`, invalid checkpoint intervals when checkpoint output is enabled, invalid electrode drive amplitude/frequency/phase combinations, non-positive particle limits/output strides, malformed collision channels/tables or unsafe collision-rate bounds, empty 2D boundary tags, non-finite magnetic-field values, invalid source schedules/velocities/references, invalid emission yields/limits/references, and invalid species initialization intervals are rejected instead of silently falling back to defaults. Emission rules must target an absorbing boundary, and unsafe macro-particle expansion is rejected during construction. For structured species definitions, provide either an explicit positive `weight` or omit `weight` and provide a positive `density`; the loader converts density to macro-particle weight over the configured initialization interval or area. With a nonuniform profile, this density fixes the total represented population (equivalently the volume-average density); the profile fixes its normalized relative spatial shape.

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

The verified smoke/performance envelope is documented in `docs/performance-envelope.md`, and the quantitative Landau-damping, two-stream, and 2D/3D Langmuir cases in `docs/kinetic-validation.md`. The staged Hall-thruster targets and strict production-only resource policy are documented in `docs/hall-thruster-validation.md`. Imported scalar diagnostics expose cumulative particle, deposition, and field-solve timings plus location-cache hits and spatial searches, and `scripts/benchmark_unstructured.py` reports repeat medians for a chosen imported config. The checked-in examples prove that the documented 1D/2D/3D CLI paths, diagnostics, VTK output, particle samples, prescribed uniform-B Boris activation, and checkpoint-style text outputs remain structurally valid at small CI-friendly sizes. The kinetic cases additionally verify published damped and unstable collisionless responses plus analytic 2D/3D directional plasma oscillations; neither suite proves convergence for arbitrary plasma regimes or validates a real device. Before using larger runs, document resolution, timestep, particles-per-cell/noise, output cadence, boundary model, and convergence checks against mesh/time/particle refinements.

This is a serious first version, not a final plasma platform. Key known gaps are: no MPI/GPU backend yet, OpenMP remains a shared-memory particle-path implementation rather than a domain-decomposed whole-solver model, external initial-state loading does not yet have a chunked openPMD/HDF5 high-volume backend, MCC thermal neutrals have fixed temperature and zero bulk flow, excitation and ionization omit neutral recoil, charge exchange is limited to the resonant mass-matched case, and there is no neutral depletion, gas heating, or general reaction network; prescribed uniform magnetic fields only (no self-consistent electromagnetic field solve yet), imported field conditions are limited to label-wise constant Dirichlet/Neumann data, and the imported runtime has not been performance-qualified on production-scale meshes. No authoritative He/Ar/Kr/Xe cross-section set is bundled yet. High-volume particle dumps are intentionally deferred to an openPMD/HDF5-style format in a later phase; current `.aps` initial states, text checkpoints, and particle CSV output are for preprocessing, restart, inspection, and regression/debug workflows.
