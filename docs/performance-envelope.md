# AuroraPIC performance and validation envelope

AuroraPIC is currently a bounded electrostatic PIC prototype, not a general-purpose plasma production platform. This document defines the envelope that the automated smoke suite and public examples are intended to cover, plus the scaling and physics assumptions users should check before trusting larger studies.

## Verified smoke envelope

The checked-in examples are intentionally small so they run in CI and on developer laptops. They validate integration paths and file structure, not high-resolution physics fidelity.

| Example | Dimension | Mesh | Particles | Steps | Main coverage |
| --- | ---: | ---: | ---: | ---: | --- |
| `examples/two_stream.cfg` | 1D1V | 128 cells | 12,000 | 300 | Periodic electrostatic transient with multiple species. |
| `examples/sheath_steady.cfg` | 1D1V | 96 cells | 6,000 | up to 2,000 | Dirichlet boundaries, absorbing-wall loss, collisions, steady-state stop condition. |
| `examples/plasma_2d.cfg` | 2D2V | 32 x 32 nodes | 200 | 20 | Periodic 2D field solve, VTK output, particle samples, prescribed uniform-B Boris activation. |
| `examples/electrode_2d.cfg` | 2D2V | 32 x 24 nodes | 160 | 10 | Dirichlet electrode fields and mixed particle-boundary policies. |
| `examples/plasma_3d.cfg` | 3D3V | 8 x 8 x 8 nodes | 128 | 3 | Structured 3D CLI path, VTK legacy/XML output, particle samples. |

`scripts/verify.sh` builds the code and runs these examples through `scripts/verify_examples.py`, which checks that scalar histories, field snapshots, VTK files, and sampled particle files are structurally valid. Passing this suite means the documented smoke envelope works; it does not establish convergence for arbitrary plasma regimes.

## Practical scaling expectations

The current implementation is aimed at correctness and regression coverage before whole-solver performance. Expect memory and runtime to scale approximately with:

- particles advanced per step: `O(total_particles)`;
- structured field storage: `O(nx * ny * nz)` for the active dimension;
- direct structured Poisson solve cost depending on boundary mode and dimension;
- output volume proportional to written field nodes, particle sample count, and output frequency.

Use the serial backend as the portability baseline. Optional OpenMP currently covers safe particle-loop slices and uses deterministic static scheduling; it is not yet a full MPI/GPU or domain-decomposed scaling model. Treat `runtime_backend = mpi` and `runtime_backend = gpu` as reserved future options that intentionally fail fast.

Transient and steady-state execution are available in all structured dimensions. Steady-state termination compares adjacent windows of emitted total-energy diagnostics and always remains bounded by `max_steps`. This is an operational convergence signal only; production studies must also establish problem-specific field, charge, flux, and distribution-function convergence.

## Before using larger runs

For credible physical studies, document these checks with the run configuration and diagnostics:

1. **Resolution:** cell size resolves the relevant Debye length and geometry scale.
2. **Time step:** `dt` is below electron plasma, cyclotron when using prescribed `B`, and particle transit-time limits.
3. **Noise:** particles per cell are high enough for the desired charge/current noise tolerance.
4. **Convergence:** repeat with smaller `dt`, finer mesh, and/or more particles; compare scalar histories and field snapshots.
5. **Output cadence:** output intervals are short enough to detect transients but not so frequent that I/O dominates.
6. **Boundary model:** particle and field boundaries match the intended physical problem; current imported Gmsh support is topology/tag import only, not an unstructured field solve.
7. **Physics scope:** current fields are electrostatic Poisson fields plus optional prescribed uniform magnetic rotation. There is no self-consistent electromagnetic field update yet.

## Release-engineering envelope

The M6 baseline now includes:

- strict `config_version = 1` compatibility checks;
- a CI workflow that exercises serial/OpenMP configuration variants across Linux and macOS runners;
- CPack `TGZ` packaging rules for the CLI, library, headers, examples, documentation, and installable CMake package metadata;
- an install/package smoke test that installs the built tree, runs the installed CLI, builds a downstream `find_package(AuroraPIC CONFIG REQUIRED)` consumer, and inspects the generated `TGZ` package;
- this documented performance envelope.

The remaining hardening work is to add measured benchmark history, publish signed release artifacts, expand compiler/platform coverage, and implement deeper physics validation for electromagnetic fields and improved collision models when those features exist.
