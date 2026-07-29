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
# Optional decimal or 0x-prefixed canonical semantic signature:
# initial_state_signature = 123456789
```

`initial_state_path` and `restart_path` are mutually exclusive. Every
configured species must occur in the state file, no unknown species are
allowed, and each record count must exactly match that species' configured
`particles` value. Charge, mass, and constant macro-particle weight remain
species configuration properties. Per-particle weights are not supported in
version 1.

`initial_state_signature` pins the canonical 64-bit FNV-1a semantic signature
of the declared metadata, sorted species populations, and every binary64
position/velocity value. Decimal and `0x`-prefixed values are accepted. A
changed record, species name, population, dimension, or unit contract fails
before particle loading when a signature is pinned. The signature is a
deterministic integrity mechanism, not a cryptographic authenticity
guarantee.

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
loading, and density profile as `external`. External runs also write
`initial_state_metadata.txt` with the resolved source path, realized signature,
optional expected signature, dimension, units, and particle count.

The public `write_external_particle_state` API emits deterministic version 1
files in sorted species order with round-trip-stable binary64 text precision.
It refuses to replace an existing file unless its caller explicitly enables
overwrite. `external_particle_state_signature` lets preprocessors compute the
same semantic signature before writing a run configuration.

## Scope and high-volume path

The public API supports both an in-memory state object for preprocessors and
`load_validated_external_particle_state_bounded` for runtime ingestion. The
bounded path verifies metadata, species counts, the canonical signature, and
file stability before delivering any records. Its consumer receives the
configured species index, the record index within that species, and one
time-centered record. Structured and imported simulation paths use this API
to populate their own particle arrays directly, avoiding a second
particle-sized record container.

Canonical signatures group records by sorted species name, while version 1
permits interleaved records. The bounded text reader therefore scans the file
once per configured species to reproduce the canonical signature, then once
to deliver validated records. Auxiliary reader memory is proportional to the
species count, at the cost of repeated text parsing. This tradeoff is
appropriate for portable preprocessing and moderate data sets, not
high-volume production I/O.

The consumer boundary is deliberate: a future optional openPMD/HDF5 reader
can deliver chunks into the same simulation-owned arrays without changing
geometry validation or particle conversion. Version 1 is not an openPMD file
and should not be presented as one. Large production populations still
require the planned chunked HDF5/openPMD backend.
