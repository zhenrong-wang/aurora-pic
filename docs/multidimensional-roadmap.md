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

The first production-facing geometry workflow should be:

```text
CAD/surface description -> external mesher -> tagged mesh file -> AuroraPIC mesh importer
```

This avoids making AuroraPIC responsible for CAD kernels, boolean operations, or meshing quality.

## Internal architecture targets

### Mesh and geometry

Introduce a dimension-aware mesh layer independent of the current 1D `Grid` implementation:

- `Mesh1D`, `Mesh2D`, and eventually `Mesh3D`, or a templated `Mesh<Dim>` after the 2D design stabilizes.
- Cell-centered and node-centered coordinate access.
- Boundary tags such as `left`, `right`, `wall`, `inlet`, `outlet`, `absorbing`, or user-defined physical names imported from mesh files.
- Structured Cartesian mesh support first; unstructured element connectivity later.

Do not expose solver code directly to Gmsh, VTK, or CAD concepts. Importers should translate external files into AuroraPIC's internal mesh and boundary-tag model.

### Particles

Replace scalar position and velocity storage with vector-valued state:

- position `x[Dim]`;
- velocity `v[VelDim]` (`2D3V` should be allowed after `2D2V` basics work);
- species charge, mass, macro-particle weight, and live/dead state;
- boundary interaction policy for reflection, absorption, emission, or reinjection.

Prefer a structure-of-arrays layout once particle count becomes large, but keep the first 2D implementation simple and testable.

### Deposition and interpolation

Generalize cloud-in-cell operations:

- 1D: two-node weighting;
- 2D structured: four-node bilinear weighting;
- 3D structured: eight-node trilinear weighting;
- unstructured: element-local shape functions in a later phase.

Charge conservation and boundary behavior should be covered by tests before adding more physics.

### Field solve

Stage the field solver work:

1. 2D structured Cartesian electrostatic Poisson solver with periodic and Dirichlet cases.
2. 2D boundary-tag support for walls and electrodes.
3. 3D structured electrostatic solver.
4. Unstructured finite-element or finite-volume Poisson path only after mesh import and boundary tags are stable.
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
| M0 | Current baseline | Regression-preserving multidimensional PIC core | 1D/2D/3D CLI examples smoke successfully; CTest covers Poisson, pusher, checkpoint/restart, strict config validation, multidimensional diagnostics, and Boris activation; standalone pusher validation passes. |
| M1 | Next | Validation and benchmark suite | Add named analytic/symmetry benchmarks for deposition, interpolation, Poisson solves, particle boundaries, restart determinism, and representative 2D/3D examples with documented tolerances. |
| M2 | Planned | Geometry and mesh import workflow | Import tagged 2D Gmsh meshes into internal boundary labels without exposing solver code to file-format details; document CAD/surface-to-mesh preprocessing. |
| M3 | Planned | Scalable data and restart formats | Add production-oriented field/particle output and checkpoint formats (VTK XML and later openPMD/HDF5) with compatibility tests and migration notes from CSV/text checkpoints. |
| M4 | Planned | Runtime scaling backend | Introduce OpenMP/MPI/GPU-ready interfaces only after validation remains stable; include deterministic single-rank comparisons and scaling smoke tests. |
| M5 | Planned | Higher-fidelity physics | Extend beyond prescribed electrostatic/uniform-B operation with self-consistent electromagnetic fields and improved collision models, each guarded by conservation and benchmark tests. |
| M6 | Planned | Release engineering and operability | Add packaged builds, CI matrix coverage, versioned configuration compatibility, documented performance envelopes, and clearer failure diagnostics for invalid production inputs. |

### Immediate coding target

The current M1 batch adds named physics-facing benchmark cases for analytic periodic Poisson solves across 1D/2D/3D, exact CIC shape-function deposition, affine electric-field interpolation, analytic one-step 2D particle-boundary policies, and checkpoint/restart determinism across 1D/2D/3D. The machine-checkable milestone contract (`scripts/validate_milestones.py`) keeps the smoke suite pinned to the production ladder. Future M1 work should continue adding physics-facing benchmark cases rather than only documentation checks.

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
