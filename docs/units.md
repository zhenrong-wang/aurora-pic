# Units contract

AuroraPIC supports two explicit global unit systems:

```ini
units = normalized
relative_permittivity = 1.0
```

or:

```ini
units = si
relative_permittivity = 1.0
```

Omitting `units` preserves the historical `normalized` behavior for legacy
configuration files. All maintained examples declare it explicitly.
`relative_permittivity` must be positive and finite and describes one
homogeneous domain-wide material. Spatially varying dielectric coefficients
are not yet supported.

## Normalized mode

Normalized mode uses base permittivity `1`; the Poisson permittivity is
therefore `relative_permittivity`. Every supplied quantity must belong to one
self-consistent normalization chosen by the user. AuroraPIC does not infer or
convert reference length, time, mass, charge, voltage, magnetic field, or
density scales.

The common examples use `relative_permittivity = 1` implicitly. This exactly
preserves all solver and diagnostic results from before the units contract.

## SI mode

SI mode uses vacuum permittivity
`8.8541878128e-12 F/m`, multiplied by `relative_permittivity`.

Inputs have these meanings:

| Quantity | SI unit |
| --- | --- |
| position and mesh coordinates | m |
| timestep and time | s |
| particle velocity and thermal velocity standard deviation | m/s |
| potential | V |
| electric field | V/m |
| prescribed magnetic field | T |
| species charge | C |
| species mass | kg |
| collision frequency | 1/s |
| MCC neutral density | particles/m³ |
| MCC neutral temperature | K |
| MCC cross section | m² |
| MCC table and threshold energy | J |
| 3D number density | particles/m³ |

`thermal_velocity` is a velocity-distribution standard deviation, not a
temperature in kelvin or electron-volts. Convert a physical temperature to
the desired thermal-speed convention before supplying it. In contrast,
positive SI MCC `neutral_temperature` is in kelvin and the collision kernel
derives its one-component neutral standard deviation as
`sqrt(k_B * neutral_temperature / neutral_mass)`.

Structured 2D provides a strict SI convenience form: species
`temperature_ev`, pair-source `first_temperature_ev` /
`second_temperature_ev`, and `current_source_temperature_ev` are converted to
the one-component Maxwellian standard deviation
`sqrt(ELEMENTARY_CHARGE_SI * temperature_ev / mass_kg)`. These keys are
mutually exclusive with their velocity forms and are rejected under
`units = normalized`.

Reduced-dimensional electrostatic PIC represents the omitted dimensions per
unit measure:

| Runtime | Macro-particle weight | Charge density | Energy diagnostic |
| --- | --- | --- | --- |
| 1D | particles/m² | C/m³ | J/m² |
| Structured 2D planar | particles over configured depth | C/m³ | J |
| Imported 2D planar | particles/m | C/m³ | J/m |
| 3D | particles | C/m³ | J |

Accordingly, density-derived weights remain dimensionally consistent:
`density * initialization_length / particles` in 1D,
`density * initialization_area * out_of_plane_depth / particles` in
structured 2D, `density * initialization_area / particles` in imported 2D, and
`density * initialization_volume / particles` in 3D.

## Solver, diagnostics, and restart behavior

The selected permittivity is propagated through periodic and Dirichlet
structured Poisson solvers, the mixed-boundary unstructured finite-element
solver, and field-energy diagnostics. CLI startup prints the selected unit
system and effective permittivity.

Each run writes `units.txt` beside `scalars.csv`. It records the unit system,
relative and effective permittivity, spatial dimension, and the dimensional
basis of macro-particle weight and energy.

For SI structured-2D resolved diagnostics, species profile number density is
in `m^-3`, current density is in `A/m^2`, velocity and thermal speed are in
`m/s`, and scalar temperature is in eV. In normalized mode the same columns
retain normalized density/current/velocity values and `temperature_ev` is
written as `nan`; AuroraPIC does not invent an SI temperature scale.

Current 1D checkpoint v14, structured 2D checkpoint v10, structured 3D
checkpoint v2, and imported checkpoint v6 record and validate the unit
contract. The 1D v4 format also records velocity dimensionality; 1D v1-v3
can initialize only 1D1V. Historical structured v1
and imported v1–v3 checkpoints remain readable only with normalized units and
`relative_permittivity = 1`, because those formats contain no unit metadata.

## Current limitations

- No automatic normalized-to-SI conversion is performed.
- Permittivity is homogeneous and scalar.
- Planar structured 2D uses an explicit positive extrusion depth, defaulting
  to one. It is not an axisymmetric metric.
- Axisymmetric metric factors are not implemented.
- Configuration cannot prove that user-provided charge, mass, voltage,
  density, and velocity values are mutually realistic; validation studies
  must document their provenance.
