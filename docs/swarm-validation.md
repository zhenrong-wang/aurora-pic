# Electron-swarm validation

`aurorapic_swarm` exercises the same three-velocity null-collision kernel used
by imported-geometry simulations in a homogeneous prescribed electric field.
It isolates gas transport from mesh, wall, and space-charge effects and scans
one or more reduced electric fields:

```bash
aurorapic_swarm examples/synthetic_swarm.swarm
```

The checked-in example and its cross sections are synthetic interface data,
not argon or another material. For physical work, set `gas_data_file` to a
locally imported version-2 SI gas package and use the collision-frequency
recommendation from its `audit.json`.

## Configuration

```ini
swarm_config_version = 1
gas_data_file = local-gases/argon/Ar.gas
neutral_density = 2.4e20
reduced_fields_td = 1, 2, 5, 10, 20, 50, 100
max_frequency = 2.0e8
timestep = 2.5e-10
steps = 20000
burn_in_steps = 10000
particles = 10000
population_model = fixed_population_no_avalanche
uncertainty_blocks = 10
work_item_limit = 1500000000
initial_mean_energy_ev = 0.1
max_energy_ev = 100
seed = 271828
output_file = argon_swarm.csv
```

The runner requires `max_frequency * timestep <= 0.1`. `max_energy_ev` must
not exceed the upper energy of any collision table, and the run fails as soon
as an electron exceeds it. This turns energy coverage into an enforced
precondition rather than constant extrapolation beyond measured data. The
sampling-step count (`steps - burn_in_steps`) must divide evenly into
`uncertainty_blocks`. A scan is also bounded to 10 million allocated
particles and, by default, 100 million particle-step-field work items. Raising
`work_item_limit` is an explicit acknowledgement that a larger run belongs on
an appropriately monitored compute host.

`population_model` selects one of two explicit contracts:

- `fixed_population_no_avalanche` is the default transport mode. Ionization
  events are counted and change the primary electron energy, but their
  secondaries do not join the ensemble.
- `branching_resampled` adds each ionization secondary with its parent's
  statistical weight, then systematically resamples back to `particles`
  computational electrons after every step while preserving total represented
  electron weight. This keeps computational cost fixed while the physical
  population grows.

Branching runs may set `population_limit` to a value greater than `particles`.
It caps the temporary pre-resampling population; if omitted, a conservative
bounded value is selected. The temporary population plus the resampling
target may not exceed 10 million allocations. This is a memory fail-safe, not
a convergence control: production studies still need particle-count and
timestep refinement.

The CSV embeds the dataset identity, version, citation, provenance, retrieval
date, license text, gas-manifest path, numerical controls, per-field seed, and
collision-model signature. The signature fingerprints cross-section and
angular tables as loaded by the kernel.
For each E/N it reports:

- signed mean electron velocity along the electric field and the conventional
  electron drift velocity opposite the field;
- reduced mobility `N * drift_velocity / electric_field`;
- mean electron energy;
- longitudinal and transverse endpoint diffusion estimates in fixed mode;
- maximum observed energy;
- null-collision candidate counts and per-channel rates;
- block standard errors for drift and mean energy, plus Poisson counting
  errors for channel rates.

Branching mode additionally reports initial/final represented electron
weight, the fixed final computational-particle count, and a temporal growth
rate fitted to `log(total electron weight)`. When the conventional electron
drift is positive it also reports `growth rate / flux drift velocity` as a
clearly named Townsend approximation. `diffusion_available` and
`townsend_available` distinguish unavailable quantities from numerical zero;
diffusion columns are empty in branching mode because resampling breaks the
independent endpoint-lineage estimator.

The reduced-mobility unit is `1 / (V m s)`. One Townsend is
`1e-21 V m^2`.

## Comparison with reference coefficients

Reference data remain local, just like imported cross sections. Place the
measured or evaluated coefficients in a strict CSV with
`reduced_field_td` and one or more observable columns, then describe the
source and acceptance contract in a `.swarm-reference` manifest:

```ini
[reference]
swarm_reference_version = 1
data_file = argon-reference.csv
reference_id = laboratory.argon.swarm
reference_version = 2026-01
gas = Ar
population_model = fixed_population_no_avalanche
coefficient_convention = flux_fixed_population
provenance = laboratory or evaluated database and source URL
citation = citation requested by the contributor
retrieved = 2026-01-15
license = applicable dataset terms
field_absolute_tolerance_td = 1e-12
field_relative_tolerance = 1e-12

[observable.drift]
simulation_column = electron_drift_velocity_m_s
reference_column = drift_velocity_m_s
simulation_uncertainty_column = mean_velocity_x_standard_error_m_s
reference_uncertainty_column = drift_velocity_standard_uncertainty_m_s
relative_tolerance = 0.05
absolute_tolerance = 0
uncertainty_multiplier = 2

[observable.mean_energy]
simulation_column = mean_energy_ev
reference_column = mean_energy_ev
simulation_uncertainty_column = mean_energy_standard_error_ev
reference_uncertainty_column = mean_energy_standard_uncertainty_ev
relative_tolerance = 0.05
absolute_tolerance = 0
uncertainty_multiplier = 2
```

Run the local comparison with:

```bash
python3 scripts/compare_swarm.py argon_swarm.csv \
  argon.swarm-reference --output argon-comparison.json
```

For each observable and E/N point, acceptance is:

```text
abs(simulation - reference)
  <= absolute_tolerance
   + relative_tolerance * abs(reference)
   + uncertainty_multiplier
     * hypot(simulation_uncertainty, reference_uncertainty)
```

Uncertainty columns are optional; omitting one contributes zero to the
combined uncertainty. Every reference E/N must match exactly one simulated
point within the declared field tolerance. The comparator does not silently
interpolate reference data. Extra simulated points are retained in the report
but do not affect acceptance. All simulation rows must also share one
cross-section dataset ID/version and match the manifest's gas and
`population_model`. `coefficient_convention` is recorded explicitly because
bulk and flux swarm coefficients are not interchangeable; the current
fixed-population mean-velocity estimator is described as
`flux_fixed_population`.

The JSON records the reference metadata, acceptance rule, residuals,
uncertainties, per-value decisions, unmatched simulated fields, and SHA-256
of the simulation CSV, reference CSV, and manifest. Exit status is zero when
all values pass, one when a valid comparison misses its criteria, and two for
an invalid or ambiguous input. Existing reports require `--overwrite`.

## Model boundary

The fixed mode intentionally excludes multiplication. The
`branching_resampled` mode provides bounded electron-impact avalanche
multiplication, but its growth-over-flux-drift result is not a spatial
steady-state bulk Townsend coefficient. It does not yet include attachment,
photoionization, space charge, or a dedicated steady Townsend experiment.
Ionization uses the engine's current equal-sharing energy model. Both modes
assume stationary zero-temperature neutrals. Elastic scattering is isotropic
unless the gas package explicitly supplies a validated energy-dependent
Henyey-Greenstein mean-cosine table; configured neutral mass remains active
in either elastic recoil path.

Consequently, this benchmark can validate the current MCC implementation's
drift, mean energy, fixed-mode diffusion trend, collision rates, and bounded
transient avalanche growth. It cannot yet claim high-accuracy transport for
datasets requiring a full differential angular cross section beyond the
mean-cosine phase-function approximation, thermal neutral motion, non-equal
ionization energy sharing, or spatial steady-state Townsend coefficients.

## Production study checklist

1. Use one complete and internally consistent collision set.
2. Compare against independent measured or evaluated swarm coefficients
   using a traceable `.swarm-reference` contract.
3. Repeat with smaller timesteps and higher particle counts.
4. Extend burn-in until block means no longer show a transient trend.
5. Confirm the maximum observed energy remains comfortably below
   `max_energy_ev`.
6. Record the imported package and `audit.json` alongside the CSV.
7. Do not tune cross sections against device results before the homogeneous
   swarm comparison is understood.
8. For branching runs, repeat with larger computational populations and verify
   both the temporal growth rate and its block uncertainty stabilize.
