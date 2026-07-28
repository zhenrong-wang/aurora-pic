# Real-case geometry and validation path

## Current integration case

`examples/biased_probe_2d.geo` is the first reproducible geometry-to-simulation
case that is more than a rectangular solver fixture. It defines a planar
0.12 by 0.08 chamber with a circular internal biased probe, local probe
refinement, and named `inlet`, `outlet`, `wall`, `probe`, and `plasma`
physical groups. Regenerate its committed portable mesh with:

```bash
python3 scripts/generate_real_case_mesh.py
```

The command requires Gmsh 4.x and explicitly requests v2 ASCII output, which
is the imported format currently supported by AuroraPIC. The committed mesh
was generated with Gmsh 4.12.1; regeneration with that version is
byte-identical.

The imported artifact contains 725 nodes, 1,342 triangular cells, and 108
boundary faces. Its acceptance gate requires the circular probe center to be
outside the plasma domain, all four boundary labels to survive import, total
area between 0.0090 and 0.0093, minimum corner angle above 30 degrees, and
maximum per-cell edge ratio below 2.0. AuroraPIC reports these quality metrics
when starting any imported-mesh run.

`examples/biased_probe_2d.cfg` then exercises:

- mixed Dirichlet/Neumann electrostatic data on external and internal faces;
- electron and ion boundary injection;
- absorbing outlet, wall, and probe interactions;
- ion-induced, weight-aware secondary-electron emission from the probe;
- species/tag-resolved particle, charge, energy, rate, and flux diagnostics;
- prescribed magnetic field, VTU field output, particle sampling, and
  checkpoint v3.

The smoke acceptance envelope requires 20 completed timesteps, 350–450 live
macro-particles, at least 250 outlet impacts, at least 80 probe impacts,
nonzero wall impacts and secondary emission, exactly 40 injected
macro-particles per source, and a final Poisson residual below `1e-6`. These
are deterministic integration bounds, not agreement with an experiment.

## What “real” means at each level

The project now satisfies a **real meshing and software-integration case** for
planar electrostatic PIC: geometry is authored in Gmsh, meshed externally,
imported through physical tags, simulated through the public CLI, and checked
end to end.

It does not yet satisfy a **validated engineering plasma case**. That requires
the following work:

1. Choose a published experiment or device with accessible geometry,
   operating conditions, material properties, and diagnostic data.
2. Pin an explicit SI/normalized-unit contract and verify dimensional
   consistency of every input and diagnostic.
3. Add the dominant device physics, likely neutral background gas,
   electron-neutral and ion-neutral collisions, energy/angle-dependent wall
   yields, and possibly an external circuit or RF waveform.
4. Run particle-count, timestep, mesh, domain, and random-seed convergence
   studies with uncertainty bars.
5. Compare fields, currents, densities, energy distributions, and wall fluxes
   against analytic limits and experimental measurements.
6. For devices that are not faithfully planar, add axisymmetric support or
   imported 3D topology and scalable parallel field/particle infrastructure.

The current biased-probe case is therefore the correct gateway test before
selecting and claiming a validated device simulation.
