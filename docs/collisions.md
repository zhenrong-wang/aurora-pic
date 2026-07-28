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

The MCC slice supports stationary-neutral elastic and excitation channels for
one named kinetic species. Imported 2D3V additionally supports bounded
electron-impact ionization product creation and resonant ion-neutral charge
exchange:

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

[collision.ionization]
type = ionization
cross_section_file = ionization.dat
threshold_energy = 1.0
secondary_species = electrons
ion_species = ions

[collision.charge_exchange]
type = charge_exchange
cross_section_file = charge_exchange.dat
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

Without gas mass metadata, legacy elastic events preserve projectile kinetic
energy and randomize the sign of the 1D velocity. With positive
`neutral_mass`, elastic events use stationary-target two-body center-of-mass
kinematics. Imported 2D3V samples the post-collision relative direction
isotropically; projectile plus implicit neutral recoil conserve momentum and
total kinetic energy, while the tracked projectile can gain or lose energy.
Excitation removes exactly `threshold_energy` and retains the heavy-neutral
approximation. A channel below its threshold has zero rate regardless of its
table.

Ionization is available only through the imported 2D3V interface. Each
accepted event removes `threshold_energy`, divides the remaining incident
electron energy equally between the scattered primary and one secondary
electron, samples both directions independently and isotropically, and creates
one stationary ion at the event position. The target and secondary species
must have identical nonzero charge, mass, and macro weight; the ion must have
the opposite charge and the same macro weight. These constraints make each
macro-event charge conservative and make the implemented electron energy
partition well-defined. Product storage is preflighted against
`max_particles_per_species`, and new particles do not collide again during
their birth timestep.

This is a deliberately bounded ionization model: ion recoil, differential
scattering data, neutral depletion, metastables, and multi-ionization are not
yet represented. The bundled `examples/imported_ionization_2d.cfg` and
`mcc_2d3v_ionization.dat` are deterministic software-validation inputs, not
material data.

Resonant charge exchange represents the identity swap
`A+_fast + A_slow -> A_fast + A+_slow`. It requires projectile mass equal to
`neutral_mass`, a charged target species, zero threshold, and imported 2D3V
configuration. Because the present neutral background is stationary, an
accepted event resets the kinetic ion velocity to zero without changing
particle count or charge. The outgoing fast neutral is absorbed into the
untracked neutral reservoir.

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

For reusable external data, replace the inline `gas`, `neutral_mass`,
`data_provenance`, and channel physics with a versioned manifest:

```ini
# argon.gas
gas_data_version = 2
units = si
gas = argon
neutral_mass = 6.6335209e-26
dataset_id = provider.argon.electron-neutral
dataset_version = 2026-01
data_provenance = provider dataset name and source URL
citation = citation requested by the data contributor
retrieved = 2026-01-15
license = applicable dataset terms

[collision.elastic]
type = elastic
cross_section_file = argon_elastic.dat

[collision.ionization]
type = ionization
cross_section_file = argon_ionization.dat
threshold_energy = 2.524e-18
```

The simulation supplies operating conditions and reactive species mappings:

```ini
[collisions]
enabled = true
model = null_collision
species = electrons
gas_data_file = argon.gas
neutral_density = 2.4e20
neutral_temperature = 300.0
max_frequency = 1.0e8

[collision.ionization]
secondary_species = electrons
ion_species = argon_ions
```

Manifest table paths are resolved relative to the manifest, so one package can
be reused by simulations in different directories. `gas_data_version = 2`,
an explicit `si` or `normalized` unit contract, gas identity, positive neutral
mass, stable dataset ID and version,
provenance, citation, a valid `YYYY-MM-DD` retrieval date, license text, and at
least one valid channel are mandatory. Unknown or duplicate keys, malformed
tables, invalid thresholds, and missing files fail during loading. Simulation
files cannot override packaged channel type, table, scales, or threshold; they
only map ionization products to configured kinetic species. The public
`pic::load_gas_dataset` API exposes the same validated manifest contract to
embedding applications.

Version 1 manifests remain readable as legacy normalized data. Version 2 is
required for newly converted physical datasets and is rejected when its unit
system differs from the simulation. The complete local conversion and audit
workflow is documented in
[`gas-data-workflow.md`](gas-data-workflow.md).

The neutral background remains stationary. Neutral mass is active in elastic
recoil and resonant charge-exchange validation; neutral temperature remains
provenance and forward-compatibility metadata and does not yet add thermal
neutral velocity. Gas heating and depletion are not modeled.

`examples/imported_mcc_2d.cfg` exercises the complete imported parser,
isotropic scatter, diagnostics, and v6 restart path using deliberately
synthetic constant cross sections. It is not an Argon or other material model.
`examples/imported_ionization_2d.cfg` similarly exercises paired reactive
product creation, capacity enforcement, energy accounting, and deterministic
restart.

## Diagnostics and restart

An enabled collision model writes `collisions.csv`, containing interval and
cumulative candidate, null-collision, and named-channel counts at the scalar
output cadence. Imported runs also write `collision_data.txt`, recording the
resolved gas metadata, operating state, model signature, effective channel
settings, table paths, and product mappings used by that run.

1D checkpoint v3 and imported checkpoint v6 record collision model identity, a
fingerprint of effective cross-section tables and MCC parameters, cumulative
collision counts, and RNG state. They reject restart with changed collision
data, external dataset metadata, or model parameters. Historical 1D v1/v2 and
imported v1-v5 checkpoints
cannot restart null-collision MCC because they contain no compatible MCC
provenance. Historical 1D v3 MCC signatures remain compatible when the new
optional gas metadata is absent. Checkpoints made with the earlier
infinite-neutral-mass imported MCC signature are intentionally rejected after
enabling finite-mass recoil.

## Current limitations

- Collision sampling is currently serial to preserve deterministic RNG order.
- Neutrals are stationary. Elastic recoil is finite-mass, but thermal motion,
  excitation/ionization recoil, depletion, and gas heating are absent.
- Structured 2D and structured 3D do not yet expose MCC configuration.
- Ionization is limited to the equal-sharing imported 2D3V model above.
  Charge exchange is limited to the resonant mass-matched model above.
  Attachment, Coulomb collisions, and general reaction networks are not
  implemented.
- Cross-section licensing and provenance are the user's responsibility. The
  checked-in MCC tables are synthetic software-validation data, not He, Ar,
  Kr, Xe, or other material data.

For a real gas, use an authoritative threshold and a separately licensed,
cited cross-section dataset. As one reference point, the
[NIST Atomic Spectra Database](https://physics.nist.gov/cgi-bin/ASD/ie.pl?at_num_out=1&biblio=1&e_out=0&el_name_out=1&ion_charge_out=1&spectra=argon&unc_out=1&units=1)
reports the first Argon ionization energy; this does not provide the required
energy-dependent cross section. Dataset redistribution and use rights must be
checked at the source—for example, the
[LXCat terms](https://nl.lxcat.net/data/preselect.php?a=181&b=8&d=26&t=swrm)
place citation, retrieval-date, redistribution, and commercial-use conditions
on contributed data. AuroraPIC therefore does not vendor those tables.
