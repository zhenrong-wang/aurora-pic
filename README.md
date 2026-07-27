# AuroraPIC

AuroraPIC is a C++20 starting point for scientific plasma dynamics simulation. The current codebase implements electrostatic `1D1V`, structured `2D2V`, and structured `3D3V` Particle-in-Cell (PIC) paths with configurable species, periodic or Dirichlet boundaries, transient fixed-step and steady-state convergence modes, scalar diagnostics, and text checkpoint/restart files. The 1D baseline also provides optional Monte-Carlo collisions. The multidimensional paths provide prescribed uniform magnetic-field Boris pushes, VTK field output, side-specific particle boundaries, and optional particle inspection CSVs.

## Why this methodology

PIC is a standard kinetic method for low-collisionality plasma simulation because it evolves macro-particles while solving fields on a mesh. AuroraPIC keeps the 1D path as the regression baseline while adding structured 2D building blocks before moving to geometry import, 3D, electromagnetic fields, and accelerated parallel backends.

For the recommended multidimensional expansion strategy, geometry/mesh format choices, and staged implementation plan, see `docs/multidimensional-roadmap.md`.

## Production milestone baseline

Production readiness is now tracked as explicit milestones instead of an open-ended roadmap narrative. The pinned milestone ladder and evidence expectations live in `docs/multidimensional-roadmap.md#production-readiness-milestone-ladder`; `scripts/validate_milestones.py` is part of the smoke suite and fails if those milestone IDs or README linkage drift. The current M6 baseline includes release-engineering and operability mitigation: configs may declare `config_version = 1`, public examples do so explicitly, unsupported future config versions fail with a clear diagnostic, CI workflow coverage exercises serial/OpenMP build variants across Linux/macOS runners, CPack can produce a `TGZ` package, and `docs/performance-envelope.md` documents the verified smoke/performance envelope. It also preserves the completed M5 higher-fidelity physics slice, M4 runtime-scaling interface slice, M3 VTK XML structured-grid output compatibility, and the M2 tagged 2D Gmsh v2 ASCII importer (`ImportedMesh2D`) for externally meshed planar domains.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

OpenMP support is enabled by default when CMake finds a C++ OpenMP toolchain. Disable it explicitly with `-DAURORA_ENABLE_OPENMP=OFF` to force serial-only builds.

To install the current build into a local prefix, create the `TGZ` package, and smoke-test the installed CLI plus downstream CMake package metadata, run:

```sh
python3 scripts/verify_install_package.py build
```

Installed downstream projects can consume the library target with:

```cmake
find_package(AuroraPIC CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE AuroraPIC::aurorapic)
```

## Verify

The full smoke suite builds the project, validates milestone and release-engineering artifacts, runs the CTest regression executable, runs the standalone pusher validation script (leapfrog plus Boris checks), runs isolated CLI smoke tests for the included 1D/2D/3D examples, and runs the install/package smoke test for the installed CLI, CPack `TGZ`, and downstream `find_package(AuroraPIC CONFIG REQUIRED)` consumer:

```sh
scripts/verify.sh
```

The example smoke tests copy each config to a temporary `test_output_aurorapic_verify/` directory, rewrite only `output_dir`, run `aurorapic_cli`, and assert the expected scalar, field, VTK, and particle-inspection files are structurally valid. Temporary smoke outputs are removed on success; set `KEEP_VERIFY_OUTPUTS=1` or pass `--keep-output` to `scripts/verify_examples.py` to retain them for debugging:

```sh
python3 scripts/verify_examples.py build/aurorapic_cli --keep-output
```

## 2D status

2D runs are selected by `dimension = 2` in the config file and are reachable from both the CLI/config loader and the C++ API. The 2D path provides:

- `Mesh2D`: rectangular node-centered mesh with periodic or Dirichlet field boundary mode plus side boundary tags/potentials for electrode-style Dirichlet domains.
- `ImportedMesh2D`: validated topology/label model for tagged planar Gmsh v2 ASCII imports; bounded parsing rejects non-finite coordinates, duplicate entities, degenerate or non-convex cells, inconsistent orientation, non-manifold edges, and incomplete boundary closure. The model exposes cell area/centroid and boundary-length metrics without exposing solver code to `.msh` details.
- `UnstructuredMesh2D`: computational state over validated imported topology, with lumped nodal control areas, spatially accelerated triangle/quad point location, conservative element-shape charge deposition, and nodal electric-field interpolation.
- `solve_unstructured_poisson`: triangle/quad finite-element stiffness assembly, strict physical-label Dirichlet mapping, CSR storage, Jacobi-preconditioned conjugate gradients, convergence reporting, and nodal electric-field recovery. Runtime/config integration remains pending.
- `Species2D`: explicit `Particle2D` storage with 2D position/velocity initialization, CIC deposition, kinetic-energy accounting, and live-particle accounting.
- `Simulation2D`: deposit -> solve -> particle push/drift -> redeposit/resolve loop using the existing 2D Poisson solvers, with per-side particle boundary policies (`auto`, `absorbing`, `reflecting`, `periodic`). The default zero magnetic field uses the electrostatic leapfrog pusher; nonzero `magnetic_field_z` switches particles to the Boris rotation/kick.
- `Diagnostics2D`: scalar time histories in `scalars.csv`, cumulative absorbed-particle counts by side, and optional sampled particle CSV files.
- `write_legacy_vtk` / `write_vtk_xml`: structured-grid VTK writers for `rho`, `phi`, and electric-field vectors on `Mesh2D`.

When `vtk_output = true`, 2D runs write field snapshots under `output_dir` for ParaView or VisIt. `vtk_format` selects `legacy` (`fields_<step>.vtk`, the default), `xml`/`vts` (`fields_<step>.vts`), or `both`.

All 2D runs write `scalars.csv` with:
```text
step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles,absorbed_left,absorbed_right,absorbed_bottom,absorbed_top,live_particles_<species>...
```

The absorbed-particle columns are cumulative counts of particles removed by absorbing particle boundary policies on each side. Reflecting and periodic particle boundaries keep particles live and do not increment these counters.

If `particle_output = true`, the run also writes sampled particle inspection files named `particles_<step>.csv` with:

```text
species_id,species,x,y,vx,vy,alive
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
nx = 128
length = 1.0
dt = 0.002
steps = 300
boundary = periodic        # periodic or dirichlet
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

[collisions]
enabled = false
frequency = 0.0
neutral_temperature_velocity = 0.0
```

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
# Optional uniform out-of-plane B field. Nonzero values use the Boris pusher.
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
thermal_velocity = 0.02
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
init_x_min = 0.0
init_x_max = 1.0
init_y_min = 0.0
init_y_max = 1.0
init_z_min = 0.0
init_z_max = 1.0
```

## Checkpoint/restart controls

All 1D, 2D, and 3D runs support a text `.apc` checkpoint format intended for deterministic restart and regression debugging. Checkpoints include the simulation dimension, step/time, RNG engine state, per-species particle positions/velocities/leapfrog half-step state, live flags, and 2D/3D absorbed-particle counters. Loading a checkpoint validates the format dimension plus species count/name metadata against the active config before repopulating particles and recomputing mesh fields.

- `checkpoint_output`: enable/disable checkpoint writes during `run()`; default `false`.
- `checkpoint_interval`: checkpoint interval in steps; `0` inherits `output_interval` when `checkpoint_output = true`.
- `checkpoint_path`: optional fixed checkpoint file path. If omitted, checkpoints are written as `output_dir/checkpoint_<step>.apc`. If provided, each checkpoint write updates that same path.
- `restart_path`: optional checkpoint file to load before the run loop starts. The run resumes from the checkpoint step/time and continues until the configured `steps`/`max_steps` limit.

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

The particle initialization/synchronization loops and the 1D particle advance use `RuntimePolicy` through deterministic static scheduling. The regression suite compares serial and OpenMP 2D smoke runs when OpenMP is available, while keeping serial behavior as the portability baseline.

2D/3D particle-boundary controls:

- `particle_boundary`: default particle policy for all sides; one of `auto`, `absorbing`, `reflecting`, or `periodic`.
- `particle_boundary_left`, `particle_boundary_right`, `particle_boundary_bottom`, `particle_boundary_top`: per-side overrides for the default policy; 3D also supports `particle_boundary_back` and `particle_boundary_front`.
- `auto`: resolves to `periodic` for `boundary = periodic` field solves and `absorbing` for `boundary = dirichlet` field solves.
- `absorbing`: removes particles that leave the domain and increments the corresponding `absorbed_*` scalar diagnostic.
- `reflecting`: mirrors escaped particles back into the domain and reverses the normal velocity component.
- `periodic`: wraps escaped particles across that coordinate direction.

2D/3D magnetic-field controls:

- 2D supports `magnetic_field_z`, a uniform out-of-plane magnetic field. It defaults to `0.0`; nonzero values activate the Boris pusher for 2D particles.
- 3D supports uniform `magnetic_field_x`, `magnetic_field_y`, and `magnetic_field_z`. They default to `0.0`; any nonzero component activates the Boris pusher for 3D particles.
- Magnetic-field values must be finite. The current field solve remains electrostatic Poisson; these controls add prescribed uniform magnetic rotation to particle pushes, not a self-consistent electromagnetic field update.

The parser is intentionally strict: unsupported `config_version` values, unknown sections/keys, invalid enum values, invalid particle-boundary values, invalid booleans, non-finite numbers, non-positive `dt`/`output_interval`, invalid checkpoint intervals when checkpoint output is enabled, non-positive `particle_output_stride`, empty 2D boundary tags, non-finite magnetic-field values, and invalid species initialization intervals are rejected instead of silently falling back to defaults. For species definitions, provide either an explicit positive `weight` or omit `weight` and provide a positive `density`; the loader converts density to macro-particle weight over the configured initialization interval or area.

## Performance and validation envelope

The verified smoke/performance envelope is documented in `docs/performance-envelope.md`. In short, the checked-in examples prove that the documented 1D/2D/3D CLI paths, diagnostics, VTK output, particle samples, prescribed uniform-B Boris activation, and checkpoint-style text outputs remain structurally valid at small CI-friendly sizes. They do not prove convergence for arbitrary plasma regimes. Before using larger runs, document resolution, timestep, particles-per-cell/noise, output cadence, boundary model, and convergence checks against mesh/time/particle refinements.

This is a serious first version, not a final plasma platform. Key known gaps are: no MPI/GPU backend yet, OpenMP is limited to safe particle-loop slices rather than a whole-solver scaling model, simplified collision model, prescribed uniform magnetic fields only (no self-consistent electromagnetic field solve yet), and no complete unstructured simulation/config path around the current imported particle-grid and electrostatic field-solve components. High-volume particle dumps are intentionally deferred to an openPMD/HDF5-style format in a later phase; current text checkpoint and particle CSV output are for restart, inspection, and regression/debug workflows. These extension points are documented in `docs/methodology.md`, `docs/multidimensional-roadmap.md`, and `docs/performance-envelope.md`.
