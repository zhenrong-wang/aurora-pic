# AuroraPIC

AuroraPIC is a C++20 starting point for scientific plasma dynamics simulation. The current codebase implements an electrostatic `1D1V` Particle-in-Cell (PIC) baseline with configurable species, periodic or Dirichlet boundaries, optional Monte-Carlo collisions, transient fixed-step simulation, and steady-state convergence mode. It also includes a structured `2D2V` electrostatic path (`Mesh2D`, `Species2D`, and `Simulation2D`) for periodic/Dirichlet rectangular domains, VTK field output, scalar histories, and optional particle inspection CSVs.

## Why this methodology

PIC is a standard kinetic method for low-collisionality plasma simulation because it evolves macro-particles while solving fields on a mesh. AuroraPIC keeps the 1D path as the regression baseline while adding structured 2D building blocks before moving to geometry import, 3D, electromagnetic fields, and accelerated parallel backends.

For the recommended multidimensional expansion strategy, geometry/mesh format choices, and staged implementation plan, see `docs/multidimensional-roadmap.md`.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Verify

```sh
ctest --test-dir build --output-on-failure
# or run the full smoke suite:
scripts/verify.sh
```

The full smoke suite builds the project, runs the CTest regression executable, and runs the included 1D/2D examples:

```sh
./build/aurorapic_cli examples/two_stream.cfg
./build/aurorapic_cli examples/sheath_steady.cfg
./build/aurorapic_cli examples/plasma_2d.cfg
./build/aurorapic_cli examples/electrode_2d.cfg
```

## 2D status

2D runs are selected by `dimension = 2` in the config file and are reachable from both the CLI/config loader and the C++ API. The 2D path provides:

- `Mesh2D`: rectangular node-centered mesh with periodic or Dirichlet field boundary mode plus side boundary tags/potentials for electrode-style Dirichlet domains.
- `Species2D`: explicit `Particle2D` storage with 2D position/velocity initialization, CIC deposition, kinetic-energy accounting, and live-particle accounting.
- `Simulation2D`: deposit -> solve -> leapfrog kick/drift -> redeposit/resolve loop using the existing 2D Poisson solvers, with per-side particle boundary policies (`auto`, `absorbing`, `reflecting`, `periodic`).
- `Diagnostics2D`: scalar time histories in `scalars.csv`, cumulative absorbed-particle counts by side, and optional sampled particle CSV files.
- `write_legacy_vtk`: structured-grid VTK writer for `rho`, `phi`, and electric-field vectors on `Mesh2D`.

When `vtk_output = true`, 2D runs write legacy VTK structured-grid snapshots (`fields_0.vtk`, interval snapshots, and the final `fields_<step>.vtk`) under `output_dir` for ParaView or VisIt.

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

## Configuration format

AuroraPIC uses a strict INI-like format. 1D configs may omit `dimension` or set `dimension = 1`:

```ini
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
dimension = 2
nx = 64
ny = 64
length_x = 1.0
length_y = 1.0
dt = 0.002
steps = 100
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
output_interval = 10
output_dir = output/electrode_2d
vtk_output = true
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

2D particle-boundary controls:

- `particle_boundary`: default particle policy for all sides; one of `auto`, `absorbing`, `reflecting`, or `periodic`.
- `particle_boundary_left`, `particle_boundary_right`, `particle_boundary_bottom`, `particle_boundary_top`: per-side overrides for the default policy.
- `auto`: resolves to `periodic` for `boundary = periodic` field solves and `absorbing` for `boundary = dirichlet` field solves.
- `absorbing`: removes particles that leave the domain and increments the corresponding `absorbed_*` scalar diagnostic.
- `reflecting`: mirrors escaped particles back into the domain and reverses the normal velocity component.
- `periodic`: wraps escaped particles across that coordinate direction.

The parser is intentionally strict: unknown sections/keys, invalid enum values, invalid particle-boundary values, invalid booleans, non-finite numbers, non-positive `dt`/`output_interval`, non-positive `particle_output_stride`, empty 2D boundary tags, and invalid species initialization intervals are rejected instead of silently falling back to defaults. For species definitions, provide either an explicit positive `weight` or omit `weight` and provide a positive `density`; the loader converts density to macro-particle weight over the configured initialization interval or area.

## Current limitations

This is a serious first version, not a final plasma platform. Key known gaps are: no MPI/OpenMP backend yet, no 3D solver/CLI yet, no checkpoint/restart format yet, simplified collision model, and no magnetic-field/Boris rotation yet. High-volume particle dumps are intentionally deferred to an openPMD/HDF5-style format in a later phase; current particle CSV output is for inspection and regression/debug workflows. These extension points are documented in `docs/methodology.md` and `docs/multidimensional-roadmap.md`.
