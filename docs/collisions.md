# Collision models

AuroraPIC provides a bounded 1D BGK compatibility model and tabulated
null-collision MCC for 1D1V and imported planar 2D3V runs. Collision processing
occurs after electrostatic/Boris velocity synchronization at each timestep.

## BGK compatibility model

```ini
[collisions]
enabled = true
model = bgk
frequency = 0.5
neutral_temperature_velocity = 0.03
```

Each live particle has probability `1 - exp(-frequency * dt)` of having its
velocity redrawn from a zero-mean Gaussian with the configured standard
deviation. This is a relaxation model, not a cross-section-based binary
collision model.

## Tabulated null-collision MCC

The MCC slice supports stationary-heavy-neutral elastic and excitation
channels for one named kinetic species:

```ini
[collisions]
enabled = true
model = null_collision
species = electrons
neutral_density = 5.0
max_frequency = 4.0
max_candidates_per_particle = 64

[collision.elastic]
type = elastic
cross_section_file = elastic.dat

[collision.excitation]
type = excitation
cross_section_file = excitation.dat
threshold_energy = 0.5
```

Cross-section paths are resolved relative to the configuration file. Each
table is whitespace-separated text with exactly two numeric columns:

```text
# energy  cross_section
0.0       0.0
1.0       2.0e-20
10.0      8.0e-21
```

Blank lines and `#`/`;` comments are accepted. Energies must be strictly
increasing; energy and cross section must be finite and non-negative. Linear
interpolation is used between rows, while values outside the table range use
the nearest endpoint. Tables should span the complete simulated energy range
and end with a physically justified high-energy value.

The table columns are interpreted in the active simulation's unit system.
Optional positive scales convert source-table columns:

```ini
energy_scale = 1.602176634e-19
cross_section_scale = 1.0
```

For example, this `energy_scale` converts electron-volts to joules in an SI
run. In SI mode, `neutral_density` is in `m^-3`, cross section is in `m^2`,
particle mass is in kg, and velocity is in `m/s`, producing a rate in `s^-1`.
Normalized inputs must form the corresponding self-consistent normalized rate.

For particle energy `E = m v^2 / 2`, channel frequency is
`nu_i = neutral_density * sigma_i(E) * abs(v)`. Candidate times are sampled
from an exponential distribution with `max_frequency`; a candidate is
accepted into a channel in proportion to `nu_i / max_frequency`, otherwise it
is a null collision. Multiple candidates per timestep are supported.

`max_frequency` is a user-supplied conservative bound on the sum of channel
frequencies. The run fails when an evaluated total exceeds it.
`max_candidates_per_particle` protects against an unreasonable
`max_frequency * dt`; exceeding it also fails instead of truncating collision
history.

Elastic events preserve kinetic energy and randomize the sign of the 1D
velocity. Excitation events remove exactly `threshold_energy`, then randomize
the sign. Imported 2D3V events instead sample an isotropic direction over the
unit sphere while preserving the post-collision speed. A channel below its
threshold has zero rate regardless of its table.

## Imported 2D3V gas metadata

Imported runs use the same channel sections and additionally require:

```ini
[collisions]
enabled = true
model = null_collision
species = electrons
gas = argon
neutral_density = 2.4e20
neutral_mass = 6.6335209e-26
neutral_temperature = 300.0
data_provenance = dataset name, version, and citation
max_frequency = 1.0e8
```

`gas` is an explicit identity label. In SI, `neutral_mass` is kg,
`neutral_temperature` is K, and `neutral_density` is `m^-3`; normalized runs
must use a self-consistent normalized contract. Gas identity, mass,
temperature, provenance, rate controls, target species, and effective tables
are checkpoint-fingerprinted.

The present kinematics use the stationary-heavy-neutral approximation.
Therefore neutral mass and temperature are required provenance and
forward-compatibility metadata but do not yet add thermal neutral velocity,
finite-mass recoil, gas heating, or depletion.

`examples/imported_mcc_2d.cfg` exercises the complete imported parser,
isotropic scatter, diagnostics, and v6 restart path using deliberately
synthetic constant cross sections. It is not an Argon or other material model.

## Diagnostics and restart

An enabled collision model writes `collisions.csv`, containing interval and
cumulative candidate, null-collision, and named-channel counts at the scalar
output cadence.

1D checkpoint v3 and imported checkpoint v6 record collision model identity, a
fingerprint of effective cross-section tables and MCC parameters, cumulative
collision counts, and RNG state. They reject restart with changed collision
data or model parameters. Historical 1D v1/v2 and imported v1-v5 checkpoints
cannot restart null-collision MCC because they contain no compatible MCC
provenance. Historical 1D v3 MCC signatures remain compatible when the new
optional gas metadata is absent.

## Current limitations

- Collision sampling is currently serial to preserve deterministic RNG order.
- Neutrals are stationary and infinitely heavy in the kinematics; thermal
  motion, recoil, depletion, and gas heating are absent.
- Structured 2D and structured 3D do not yet expose MCC configuration.
- Ionization, attachment, charge exchange, Coulomb collisions, and secondary
  particle creation are not implemented.
- Cross-section licensing and provenance are the user's responsibility. The
  checked-in MCC tables are synthetic software-validation data, not He, Ar,
  Kr, Xe, or other material data.
