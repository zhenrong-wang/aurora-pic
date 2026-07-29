# Multidimensional roadmap

AuroraPIC should grow from the current `1D1V` electrostatic code into 2D first, then 3D. Do not invent a bespoke 3D geometry language as the primary interface. Use open mesh, geometry, and diagnostic formats at the boundaries, and keep a small internal abstraction that the solver can own and test.

## Recommendation

1. **Keep 1D as the regression baseline.** Every multidimensional refactor should preserve the existing 1D examples and tests.
2. **Implement 2D electrostatic PIC before full 3D.** 2D exposes the important algorithmic changes--vector fields, multidimensional deposition/interpolation, boundary tagging, and larger memory footprints--without the full cost of 3D.
3. **Use open formats for geometry and data exchange.** Start with structured Cartesian meshes generated from configuration, then add importers for external meshes rather than a custom geometry interface.
4. **Design internal APIs around mesh topology and boundary labels, not around file formats.** File readers should populate the same `Mesh`, `BoundaryCondition`, and field containers used by generated grids.

## Format strategy

AuroraPIC needs three kinds of external data. They should be treated separately:

| Need | Recommended format | Role |
| --- | --- | --- |
| Mesh topology and physical groups | Gmsh `.msh` | Good first unstructured mesh import target; supports 2D/3D elements and named physical boundary groups. |
| CAD exchange | STEP, via a meshing pre-process | Keep CAD out of the PIC core; convert to a mesh with tagged boundaries before simulation. |
| Surface-only geometry | STL/OBJ, via meshing pre-process | Useful for imported surfaces, but insufficient alone for volume PIC solve domains. |
| Visualization output | VTK legacy or VTK XML (`.vti`, `.vtu`) | Easy inspection in ParaView/VisIt; better than CSV for vector/tensor fields. |
| Particle/field checkpoints | openPMD/HDF5, later phase | Portable high-volume simulation data and restart/checkpoint support. |
| Multiphysics mesh/data exchange | CGNS, later phase if needed | Stronger fit for CFD-style interoperability than for the first PIC extension. |

The first production-facing geometry workflow is:

```text
CAD/surface description -> external mesher -> tagged Gmsh v2 ASCII mesh -> AuroraPIC mesh importer
```

For the M2 importer slice, AuroraPIC reads tagged 2D Gmsh v2 ASCII `.msh` files only. CAD (`STEP`) and surface (`STL`/`OBJ`) inputs stay outside the core: prepare them in an external CAD/meshing tool, assign physical names to 1D boundaries and 2D regions, export a planar v2 ASCII mesh, then let `ImportedMesh2D` translate nodes, triangle/quadrilateral cells, and line boundary faces into AuroraPIC-owned labels. Import counts and per-element tag counts are bounded. Validation requires finite coordinates, unique entities, consistently oriented nondegenerate cells, convex quadrilaterals, manifold cell edges, and tagged faces that exactly close the domain boundary. Binary Gmsh files, Gmsh v4 files, non-planar nodes, and 3D elements are intentionally rejected or ignored until later milestones.

This avoids making AuroraPIC responsible for CAD kernels, boolean operations, or meshing quality.

## Internal architecture targets

### Mesh and geometry

Introduce a dimension-aware mesh layer independent of the current 1D `Grid` implementation:

- `Mesh1D`, `Mesh2D`, and eventually `Mesh3D`, or a templated `Mesh<Dim>` after the 2D design stabilizes.
- Cell-centered and node-centered coordinate access.
- Boundary tags such as `left`, `right`, `wall`, `inlet`, `outlet`, `absorbing`, or user-defined physical names imported from mesh files.
- Structured Cartesian meshes remain the runtime baseline; validated unstructured element connectivity and nodal computational state are now available for the imported 2D path.

Do not expose solver code directly to Gmsh, VTK, or CAD concepts. Importers should translate external files into AuroraPIC's internal mesh and boundary-tag model.

### Particles

Replace scalar position and velocity storage with vector-valued state:

- position `x[Dim]`;
- velocity `v[VelDim]` (the planar runtimes now use `2D3V`);
- species charge, mass, macro-particle weight, and live/dead state;
- boundary interaction policy for reflection, absorption, emission, or reinjection.

Prefer a structure-of-arrays layout once particle count becomes large, but keep the first 2D implementation simple and testable.

### Deposition and interpolation

Generalize cloud-in-cell operations:

- 1D: two-node weighting;
- 2D structured: four-node bilinear weighting;
- 3D structured: eight-node trilinear weighting;
- unstructured: barycentric triangle and isoparametric bilinear-quadrilateral shape functions, with a spatial candidate index for particle point location.

Charge conservation and boundary behavior should be covered by tests before adding more physics.

### Field solve

Stage the field solver work:

1. 2D structured Cartesian electrostatic Poisson solver with periodic and Dirichlet cases.
2. 2D boundary-tag support for walls and electrodes.
3. 3D structured electrostatic solver.
4. Unstructured finite-element Poisson path with physical-label Dirichlet constraints, CSR storage, preconditioned iteration, convergence reporting, and nodal field recovery.
5. Electromagnetic fields and Boris/leapfrog pusher as a separate physics milestone.

### Diagnostics and restart

CSV is acceptable for 1D scalar histories, but multidimensional fields need richer formats:

- keep `scalars.csv` for time histories;
- add VTK output for fields, charge density, and mesh inspection;
- add particle sampling output for debugging;
- consider openPMD/HDF5 once checkpointing and large particle dumps are needed.

## Production-readiness milestone ladder

The historical multidimensional phases below remain useful context, but AuroraPIC now has enough 1D/2D/3D coverage that production work should be pinned to explicit, verifiable milestones. Each milestone must preserve the existing `scripts/verify.sh` baseline unless the milestone itself intentionally updates that baseline and documents the change.

| ID | Status | Milestone | Exit evidence |
| --- | --- | --- | --- |
| M0 | Complete | Regression-preserving multidimensional PIC core | 1D/2D/3D CLI examples smoke successfully; CTest covers Poisson, pusher, checkpoint/restart, strict config validation, multidimensional diagnostics, and Boris activation; standalone pusher validation passes. |
| M1 | Complete | Validation and benchmark suite | Added named analytic/symmetry benchmarks for deposition, interpolation, Poisson solves, particle boundaries, restart determinism, and representative 2D/3D examples with documented tolerances. |
| M2 | Complete | Geometry and mesh import workflow | Added bounded tagged 2D Gmsh v2 ASCII import, geometric/topological validation, cell and boundary metrics, and internal boundary labels without exposing solver code to file-format details; documented CAD/surface-to-mesh preprocessing. |
| M3 | Complete | Scalable data and restart formats | Added VTK XML structured-grid field output with legacy/XML/both format selection and compatibility tests as the first scalable data-format step; later openPMD/HDF5 work remains planned. |
| M4 | Complete | Runtime scaling backend | Added OpenMP/MPI/GPU-ready runtime policy interfaces with optional OpenMP static-schedule execution, serial default behavior, deterministic single-rank comparisons, and scaling smoke tests. |
| M5 | Complete | Higher-fidelity physics | Added prescribed uniform-B Boris integration for 2D/3D electrostatic runs and tabulated collision models for 1D plus imported 2D3V isotropic elastic/excitation scattering and bounded ionization product creation, with CTest and CLI examples; self-consistent electromagnetic fields and broader reaction models remain later M5 extensions. |
| M6 | Current baseline | Release engineering and operability | Added config_version=1 configuration compatibility checks, clearer failure diagnostics for unsupported config versions, CI matrix coverage, CPack TGZ packaging rules, installable CMake package metadata, an install/package smoke test, and documented performance envelopes guarded by CTest, release-artifact validation, and explicit example configs. |

### Immediate coding target

The current baseline now includes bounded M6 release-engineering and operability mitigations: every loader accepts an optional `config_version = 1`, the public examples declare that version explicitly, unsupported future versions fail with a clear diagnostic before any silent fallback, `.github/workflows/ci.yml` defines a CI matrix for Linux/macOS plus OpenMP-on/off build variants, CMake installs the CLI/library/headers/examples/docs and can produce a CPack `TGZ` package with installable `find_package(AuroraPIC CONFIG REQUIRED)` metadata, and `docs/performance-envelope.md` records the documented performance envelopes for the small verified smoke examples. The install/package smoke test installs the built tree, runs the installed CLI, builds a downstream CMake consumer, and inspects the generated `TGZ` package so release artifacts are exercised as runnable software instead of only statically checked. The completed M5 higher-fidelity physics slice remains part of the baseline: prescribed uniform `magnetic_field_x/y/z` activates full three-velocity Boris rotation in planar 2D3V and 3D while retaining the electrostatic Poisson solve. The completed M4 runtime-scaling slice also remains part of the baseline: configs accept `runtime_backend`/`runtime_threads`, serial remains the default deterministic path, optional OpenMP builds use static-schedule particle loops for safe single-rank scaling smoke tests, and MPI/GPU are explicit future backends that fail fast instead of silently degrading. The M3 VTK XML slice remains part of the baseline: 2D/3D structured-grid field output can be written as legacy `.vtk`, XML `.vts`, or both while keeping legacy VTK as the default for existing configs. The M2 geometry-import slice also remains visible: tagged 2D Gmsh v2 ASCII meshes are translated into `ImportedMesh2D` nodes, cells, dimension-aware physical names, and boundary labels while keeping solver code independent of file-format details. The machine-checkable milestone contract (`scripts/validate_milestones.py`) and release-artifact guard (`scripts/validate_release_artifacts.py`) keep the smoke suite pinned to the production ladder. Remaining hardening work is to publish signed release artifacts, collect measured benchmark history across platforms, and implement deeper physics validation for electromagnetic fields and improved collision models when those features exist.

Post-M6 production hardening has unified transient and steady-state run semantics across 1D, 2D, and 3D. An explicit normalized/SI contract now propagates homogeneous permittivity through every field solver, energy diagnostic, run metadata file, and versioned checkpoint; reduced-dimensional SI weight and energy conventions are documented. Imported 2D topology is validated as a finite, nondegenerate, consistently oriented, manifold domain with exact tagged-boundary closure, and startup exposes basic mesh-quality metrics. `UnstructuredMesh2D` provides indexed nodal state, spatially accelerated point location, conservative particle-grid coupling, finite-element Poisson solution, and nodal field recovery. `UnstructuredSimulation2D` closes the transient/steady loop with area-uniform initialization, deterministic length-weighted tagged-boundary injection, segment-based absorbing/reflecting interactions, weight-aware secondary emission, species/tag-resolved physical impact fluxes, `.vtu` output, particle samples, and topology-fingerprinted deterministic restart. Its reusable Poisson workspace supports constant, label-wise mixed Dirichlet/Neumann field data and reuses assembled operators. Imported particle loops use the runtime backend with race-free deposition, while checkpoint-independent location sidecars avoid repeated spatial searches. Cumulative phase timings, cache/search counters, and a repeat-median benchmark driver make scaling measurable. Strict `[boundary.<label>]`, `[source.<name>]`, `[emission.<name>]`, and `[species.<name>]` config sections plus `mesh = imported` make this path CLI-accessible. The Gmsh-authored biased-probe case now exercises an internal curved electrode hole, local refinement, mixed probe fields, injection, collection, emission, fluxes, VTU, and checkpoint end to end. The deterministic quantitative cases in `docs/kinetic-validation.md` now verify the integrated collisionless PIC path against published Landau damping/frequency and two-stream growth, nonlinear turnover, and analytic orthogonal 2D/3D Langmuir frequency and symmetry. The Turner helium CCP ladder in `docs/ccp-validation.md` now has restart-safe sinusoidal 1D electrode forcing, a backward-compatible 1D3V path, simultaneous named electron/ion MCC, and bounded charge-conservative ionization products. The next prerequisite is benchmark-specific two-component He+-He scattering, followed by pinned collision data and phase-averaged diagnostics. Source physics still needs energy/angle-dependent material emission and distribution-based inflows; performance work continues with production-scale benchmarks, sparse/thread-memory deposition, fixed-storage locations, and PCG workspace reuse; geometry work continues with spatially varying data and higher-order/3D import.

The Hall-effect-thruster path is now pinned in
`docs/hall-thruster-validation.md`. Public LANDMARK axial-azimuthal and
radial-azimuthal PIC cases provide the code-verification targets, public WarpX
output provides an external reference corpus, and NASA HERMeS performance,
ion-velocity, and plume measurements provide later experimental-validation
targets. The first generic Hall-driven prerequisite is now implemented:
strict one-coordinate tabulated magnetic profiles share interpolation and
domain-coverage validation across structured 2D/3D and imported 2D while
retaining the uniform-field interface. Independent periodic/Dirichlet
structured 2D axes are now implemented through mesh spacing, node measures,
CIC deposit/gather, Poisson neighbors, field gradients, initialization, and
default particle policies, with analytic tests in both orientations. The
correctness-first mixed SOR field solve still needs a production
FFT-periodic/direct-axial implementation. Structured 2D now also has bounded
volumetric pair sources with shared positions, uniform/Gaussian/sinusoidal
normalized profiles, mutually exclusive fixed-macro and total represented
rates, deterministic fractional accumulation, injected-energy diagnostics,
storage caps, and restart. The next bounded HET target is a versioned
LANDMARK manifest that converts its published peak volumetric source into the
integral represented-rate contract with an explicit out-of-plane measure. Full
LANDMARK runs require exact source/current-control contracts,
statistical spectral diagnostics, MPI and scalable output; real HET validation
additionally requires measured geometry/field inputs, neutral and wall
physics, cathode/facility conditions, and uncertainty-aware comparison.
The first initial-value-problem hardening slice now adds a shared,
species-versioned loading contract to structured 1D/2D/3D and imported 2D.
It preserves historical random loading, adds deterministic stratified
quiet-start positions and antithetic velocity pairs, and supports per-component
thermal-velocity overrides for anisotropic 2D3V/3D3V distributions. A focused
low-cost regression covers every runtime path, strict parsing, geometry
containment, exact mean drift, and unsupported-combination rejection. The
second slice adds area-uniform selection by named imported cell physical group
and writes a versioned `initialization.csv` audit for generated and restarted
states, including represented charge and realized velocity moments. Analytic
density profiles now form the third slice: bounded Gaussian and sinusoidal
relative-density loading is shared by structured 1D/2D/3D and imported 2D,
with random or deterministic low-discrepancy candidates, strict
dimension-specific parameters, explicit work exhaustion, realized-profile
audit metadata, and distribution-moment regressions. Explicit charge/current
acceptance gates now form the fourth slice: generated and restarted states can
enforce normalized global charge and charge-weighted drift-current residuals,
plus explicit opposite-charge species-pair balance, with a machine-readable
audit written before any time integration. Portable external particle state
now forms the fifth slice: a strict `.aps` version 1 backend loads
time-centered records across structured 1D/2D/3D and imported 2D, validates
units/species/counts/geometry, rebuilds the self-consistent half step, and
feeds the same audits and acceptance gates. Deterministic writing, canonical
semantic signatures, optional signature pinning, and explicit source
provenance now harden that interchange contract. A validated record-consumer
path now populates simulation-owned arrays without retaining a second
particle-sized state, establishing the bounded-memory adapter boundary for a
future chunked openPMD/HDF5 backend. Native chunked storage remains the next
I/O scalability extension to this contract.
The collision baseline now includes strict tabulated cross-section ingestion,
elastic/excitation null-collision sampling, conservative rate-bound
enforcement, channel diagnostics, and checkpoint fingerprints while retaining
1D BGK compatibility. Imported 2D3V adds isotropic finite-mass neutral
scattering, required gas/provenance metadata, bounded electron-impact
ionization with paired charged products, and charge-conservative electron
attachment to a kinetic negative-ion product. Stationary-target finite-mass
elastic recoil and resonant mass-matched ion-neutral charge exchange are also
covered. Positive-temperature SI operation samples a bounded Maxwellian
neutral at each null-collision candidate, applies relative-rate and
center-of-mass kinematics, transfers the target velocity to reactive heavy
products, and enforces a conservative thermal rate majorant while preserving
the exact zero-temperature RNG path.
A strict external `.gas` manifest
interface now separates reusable table physics and licensing metadata from
simulation operating conditions and fingerprints both for restart. A local
LXCat/BOLSIG+ converter now generates SI manifests and audit reports from
user-supplied complete sets while rejecting unit, threshold, mass-ratio, and
coverage ambiguity. A homogeneous electron-swarm runner now
uses the production 3V MCC kernel for reproducible E/N scans with enforced
energy coverage, block uncertainty estimates, diffusion diagnostics, and
per-channel rates. Its optional branching/resampled contract tracks
ionization and attachment through conserved statistical weight at a fixed
computational-particle count and reports temporal growth, rate-balance
effective ionization, and clearly labeled Townsend approximations. A separate
bounded finite-distance history experiment now measures signed steady electron
flux across interior planes and fits a spatial effective Townsend coefficient
with history-block uncertainty, R², termination, population, and work gates.
A strict local
reference-comparison contract now maps
user-supplied measured or evaluated swarm coefficients to simulation columns,
combines declared numerical/reference uncertainty with explicit tolerances,
and emits a hashed per-value acceptance report without redistributing the
reference data. A strict serialized campaign manifest now runs bounded
resolution studies one at a time, applies that independent reference contract
to every run, enforces a common collision-model identity and E/N grid, checks
uncertainty-aware convergence against a designated finest run, and emits a
hashed aggregate audit. Energy-dependent anisotropic elastic scattering is now
available through a strict mean-cosine table and Henyey-Greenstein phase
function, with total-cross-section semantics, 3V-only enforcement,
energy-coverage checks, and restart fingerprinting. Executing and publishing
campaigns with authoritative gas-specific reference data, plus full
differential cross-section models, remain subsequent work. Charge exchange
beyond the resonant case, neutral flow/heating/depletion, richer ionization
kinematics, and curated measured material datasets remain future work.

## Historical implementation phases

### Phase 0: preserve and isolate 1D

- Keep existing examples and tests passing.
- Split algorithm concepts in names and docs: mesh, deposition, field solve, particle pusher, boundary interaction, diagnostics.
- Add tests around charge conservation and boundary validation.

### Phase 1: structured 2D electrostatic core

- Add `Vec<Dim>` or equivalent small vector type.
- Add `Mesh2D` with rectangular domain, spacing, and boundary labels.
- Implement 2D cloud-in-cell deposition and field interpolation.
- Implement a 2D Poisson solver for structured grids.
- Add a small 2D example that has an analytic or symmetry-based verification target.

### Phase 2: 2D boundaries and geometry import

- Add boundary-condition policies: periodic, Dirichlet electrode, absorbing wall, reflecting wall.
- Add a Gmsh `.msh` importer for boundary tags and simple triangular/quadrilateral domains.
- Keep external meshing outside AuroraPIC; document the mesh preparation workflow.

### Phase 3: 3D structured PIC

- Extend the vector, deposition, interpolation, and field containers to 3D.
- Add 3D structured Poisson support and small regression examples.
- Revisit memory layout and parallelization because 3D will stress cache, memory bandwidth, and particle count.

### Phase 4: high-volume data and acceleration

- Add VTK XML or openPMD/HDF5 output for larger runs.
- Add OpenMP/MPI/GPU backends only after the numerical interfaces are stable.
- Add restart/checkpoint support.

## Non-goals for the first production hardening step

- A custom CAD or geometry modeling language.
- Direct CAD boolean operations inside AuroraPIC.
- Full electromagnetic PIC before electrostatic multidimensional deposition and boundary handling are verified.
- Unstructured 3D before tagged 2D import and validation are stable.
