# AuroraPIC performance and validation envelope

AuroraPIC is currently a bounded electrostatic PIC prototype, not a general-purpose plasma production platform. This document defines the envelope that the automated smoke suite and public examples are intended to cover, plus the scaling and physics assumptions users should check before trusting larger studies.

## Verified smoke envelope

The checked-in examples are intentionally small so they run in CI and on developer laptops. They validate integration paths and file structure, not high-resolution physics fidelity.

| Example | Dimension | Mesh | Particles | Steps | Main coverage |
| --- | ---: | ---: | ---: | ---: | --- |
| `examples/two_stream.cfg` | 1D1V | 128 cells | 12,000 | 300 | Periodic electrostatic transient with multiple species. |
| `examples/sheath_steady.cfg` | 1D1V | 96 cells | 6,000 | up to 2,000 | Dirichlet boundaries, absorbing-wall loss, collisions, steady-state stop condition. |
| `examples/mcc_relaxation.cfg` | 1D1V | 32 cells | 1,000 | 50 | Synthetic tabulated elastic/excitation null-collision MCC with channel diagnostics. |
| `examples/plasma_2d.cfg` | 2D3V | 32 x 32 nodes | 200 | 20 | Periodic 2D field solve, VTK output, particle samples, prescribed uniform-B Boris activation. |
| `examples/electrode_2d.cfg` | 2D3V | 32 x 24 nodes | 160 | 10 | Dirichlet electrode fields and mixed particle-boundary policies. |
| `examples/imported_plasma_2d.cfg` | 2D3V | 6 nodes / 3 mixed cells | 64 initial / 70 final | 3 | Imported Gmsh CLI path, mixed-boundary FEM solve, tagged reflection/injection, VTU, particle samples, checkpoint. |
| `examples/imported_mcc_2d.cfg` | 2D3V | 6 nodes / 3 mixed cells | 64 | 6 | Synthetic stationary-neutral isotropic elastic MCC, provenance metadata, diagnostics, and checkpoint v6. |
| `examples/imported_ionization_2d.cfg` | 2D3V | 6 nodes / 3 mixed cells | 32 initial / 38 final | 4 | External synthetic `.gas` manifest, electron-impact ionization, paired charged products, threshold-energy accounting, bounded storage, and checkpoint v6. |
| `examples/biased_probe_2d.cfg` | 2D3V | 725 nodes / 1,342 triangles / internal probe hole | 800 initial / 404 final in the pinned run | 20 | Gmsh-authored real geometry, local refinement, mixed probe fields, injection, collection, secondary emission, flux diagnostics, VTU, checkpoint. |
| `examples/plasma_3d.cfg` | 3D3V | 8 x 8 x 8 nodes | 128 | 3 | Structured 3D CLI path, VTK legacy/XML output, particle samples. |

`scripts/verify.sh` builds the code and runs these examples through `scripts/verify_examples.py`, which checks that scalar histories, field snapshots, VTK files, and sampled particle files are structurally valid. Passing this suite means the documented smoke envelope works; it does not establish convergence for arbitrary plasma regimes.

The local verification entry point is deliberately resource-conservative:
one compiler job, one CTest job, and one implicit OpenMP thread by default.
Use `AURORA_BUILD_JOBS`, `AURORA_TEST_JOBS`, or
`AURORA_OPENMP_THREADS` only to opt in to greater concurrency on a dedicated
machine. Individual runtime configs can still request a tested explicit
thread count.

## Practical scaling expectations

The current implementation is aimed at correctness and regression coverage before whole-solver performance. Expect memory and runtime to scale approximately with:

- particles advanced per step: `O(total_particles)`;
- structured field storage: `O(nx * ny * nz)` for the active dimension;
- direct structured Poisson solve cost depending on boundary mode and dimension;
- imported-mesh particle lookup depending on spatial-bin occupancy; the FEM quadrature, mixed-boundary contributions, constrained CSR operator, and Jacobi diagonal are assembled once per simulation, while charge right-hand-side formation, preconditioned-CG iterations, and field recovery remain per-step costs;
- output volume proportional to written field nodes, particle sample count, and output frequency.

Use the serial backend as the portability baseline. Optional OpenMP currently covers safe particle-loop slices and uses deterministic static scheduling; it is not yet a full MPI/GPU or domain-decomposed scaling model. Treat `runtime_backend = mpi` and `runtime_backend = gpu` as reserved future options that intentionally fail fast.

Imported 2D diagnostics include cumulative `particle_seconds`, `deposition_seconds`, and `field_solve_seconds` columns, plus cumulative `location_cache_hits` and `location_searches`. Timings use a monotonic wall clock and are operational measurements, not simulation state; timings and location caches are intentionally excluded from checkpoints and numerical convergence decisions. Each live particle first validates its last element and recomputes shape coordinates locally. Only a cache miss or element crossing invokes spatial point location; restart reconstructs the sidecars during its required redeposition. Imported charge deposition uses one dense nodal buffer per active worker and reduces those buffers in worker order. This avoids atomics and data races, gives repeatable results for a fixed worker count, and uses additional memory proportional to `active_threads * mesh_nodes`.

Boundary injection can grow particle storage, so production runs must set `max_particles_per_species` from a conservative memory budget. Absorbed/dead slots are reused before growth. The configured macro-particle injection rate and expected residence time should be used to estimate the live population; the hard limit remains a safety bound rather than a population-control model.

Run a repeat-median benchmark with:

```sh
python3 scripts/benchmark_unstructured.py build/aurorapic_cli \
  examples/imported_plasma_2d.cfg --repeats 5
```

The benchmark disables field, particle, and checkpoint output, preserves the configured mesh, particle population, timestep count, and runtime backend, and reports median end-to-end and cumulative phase timings and location-cache counters. Record the compiler, build type, CPU, thread affinity, mesh size, live-particle count, step count, backend, and thread count alongside results. The checked-in smoke invocation validates the measurement path; its tiny example is not a scaling claim.

Transient and steady-state execution are available in all structured dimensions. Steady-state termination compares adjacent windows of emitted total-energy diagnostics and always remains bounded by `max_steps`. This is an operational convergence signal only; production studies must also establish problem-specific field, charge, flux, and distribution-function convergence.

## Before using larger runs

For credible physical studies, document these checks with the run configuration and diagnostics:

1. **Resolution:** cell size resolves the relevant Debye length and geometry scale.
2. **Time step:** `dt` is below electron plasma, cyclotron when using prescribed `B`, and particle transit-time limits.
3. **Noise:** particles per cell are high enough for the desired charge/current noise tolerance.
4. **Convergence:** repeat with smaller `dt`, finer mesh, and/or more particles; compare scalar histories and field snapshots.
5. **Output cadence:** output intervals are short enough to detect transients but not so frequent that I/O dominates.
6. **Boundary model:** particle and field boundaries match the intended physical problem; imported Gmsh domains are checked for manifold topology and exact tagged-boundary closure. Their particle-grid coupling uses cached element-local shapes with spatial fallback, their electrostatic solve supports strict mixed Dirichlet/Neumann labels, and the runtime applies segment-based absorbing/reflecting and length-weighted injection policies by physical label.
7. **Mesh convergence:** use the reported minimum angle and maximum edge ratio as quality guards, then repeat the physical study on at least three systematically refined meshes. The biased-probe smoke case verifies one mesh artifact; it does not establish mesh-independent physics.
8. **Units:** require an explicit `units` declaration for maintained studies, archive `units.txt`, and apply the 1D/2D per-unit omitted-measure convention when comparing energy, current, or flux with physical data.
9. **Physics scope:** current fields are electrostatic Poisson fields plus optional prescribed uniform magnetic rotation. There is no self-consistent electromagnetic field update yet.

## Release-engineering envelope

The M6 baseline now includes:

- strict `config_version = 1` compatibility checks;
- a CI workflow that exercises serial/OpenMP configuration variants across Linux and macOS runners;
- CPack `TGZ` packaging rules for the CLI, library, headers, examples, documentation, and installable CMake package metadata;
- an install/package smoke test that installs the built tree, runs the installed CLI, builds a downstream `find_package(AuroraPIC CONFIG REQUIRED)` consumer, and inspects the generated `TGZ` package;
- this documented performance envelope.

The remaining hardening work is to add measured benchmark history, publish signed release artifacts, expand compiler/platform coverage, and implement deeper physics validation for electromagnetic fields plus reactive collision models and authoritative gas datasets.
