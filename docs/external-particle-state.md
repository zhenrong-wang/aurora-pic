# External particle initial state

AuroraPIC particle-state (`.aps`) version 1 is a strict, portable interchange
format for supplying a time-zero particle population from a mesh generator,
kinetic preprocessor, analytic script, or another simulation code. It is an
initial-condition format, not a deterministic restart: simulation time,
half-step velocities, RNG state, fields, collisions, sources, and accumulated
diagnostics are not stored.

Configure it with a global path resolved relative to the config file:

```ini
initial_state_path = external_state_1d.aps
```

`initial_state_path` and `restart_path` are mutually exclusive. Every
configured species must occur in the state file, no unknown species are
allowed, and each record count must exactly match that species' configured
`particles` value. Charge, mass, and constant macro-particle weight remain
species configuration properties. Per-particle weights are not supported in
version 1.

## Version 1 grammar

The file is whitespace-delimited and has a fixed order:

```text
AuroraPIC-particle-state-v1
dimension 2
units normalized
weighting species_constant
velocity_staggering time_centered
particle_count 2
records
particle electrons 0.2 0.3 0 1.0 2.0 3.0
particle electrons 0.8 0.7 0 -1.0 -2.0 -3.0
end
```

Each particle row is:

```text
particle species x y z vx vy vz
```

The declared dimension must match the simulation. A 1D1V file requires
`y = z = vy = vz = 0`; a planar 2D3V file requires `z = 0` while retaining all
three velocity components. All numbers must be finite. `units` is exactly
`normalized` or `si` and must match the run. Only `species_constant` weighting
and time-centered velocities are accepted.

All structured positions are validated against their species' configured
initialization interval/rectangle/box and the domain; periodic domain upper
endpoints are rejected so the input has one canonical coordinate
representation. Imported 2D positions must locate inside the validated mesh
and, when configured, inside the species' rectangular initialization bounds
or named physical region.

After loading, AuroraPIC deposits the imported charge, solves the configured
Poisson problem, and derives leapfrog or Boris half-step velocities from the
time-centered input velocities and self-consistent electric field. The
realized state is then subject to the same initialization acceptance gates as
generated and restarted states. `initialization.csv` identifies the source,
loading, and density profile as `external`.

## Scope and high-volume path

The text backend is deliberately bounded by the exact configured population
and is intended for auditable small/medium preprocessing and interoperability
tests. The in-memory loader separates record validation from the storage
backend so a future optional openPMD/HDF5 implementation can map particle
position and velocity records onto the same runtime contract. Version 1 is not
an openPMD file and should not be presented as one. Large production
populations still require the planned chunked HDF5/openPMD backend.
