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
neutral_temperature = 300
reduced_fields_td = 1, 2, 5, 10, 20, 50, 100
max_frequency = 2.0e8
timestep = 2.5e-10
steps = 20000
burn_in_steps = 10000
particles = 10000
population_model = branching_resampled
uncertainty_blocks = 10
work_item_limit = 1500000000
initial_mean_energy_ev = 0.1
max_energy_ev = 100
seed = 271828
output_file = argon_swarm.csv

# Optional finite-distance steady-flux experiment:
population_limit = 20000
spatial_histories = 1000
spatial_length_m = 0.02
spatial_bins = 20
spatial_fit_begin_bin = 4
spatial_fit_end_bin = 16
spatial_max_steps = 200000
spatial_work_item_limit = 100000000
spatial_min_r_squared = 0.95
spatial_profile_file = argon_spatial_profile.csv
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
  secondaries do not join the ensemble. This mode rejects attachment because
  retaining a consumed electron would violate its fixed-population contract.
- `branching_resampled` adds each ionization secondary and removes each
  attached primary with its statistical weight, then systematically resamples
  back to `particles` computational electrons after every step while
  preserving total represented electron weight. This keeps computational
  cost fixed while the physical population grows or decays.

Branching runs may set `population_limit` to a value greater than `particles`.
It caps the temporary pre-resampling population; if omitted, a conservative
bounded value is selected. The temporary population plus the resampling
target may not exceed 10 million allocations. This is a memory fail-safe, not
a convergence control: production studies still need particle-count and
timestep refinement.

## Spatial steady-flux Townsend experiment

Setting positive `spatial_histories` enables a separate finite-distance
experiment for every configured E/N point. It is not derived from the
temporal-growth result. Each independent history injects one flux-weighted
half-Maxwellian electron at distance zero, with its initial velocity directed
into the drift domain. Electrons then move through the same prescribed field
and 3V collision kernel until they:

- cross the upstream source boundary;
- reach the downstream collector at `spatial_length_m`; or
- are consumed by attachment.

`initial_mean_energy_ev` defines the underlying isotropic Maxwellian
temperature; flux weighting changes the mean energy of particles actually
injected through the source plane.

Ionization secondaries remain in the same history. If a history's active
population exceeds `particles`, systematic resampling returns it to that
target while preserving represented weight. `population_limit` bounds the
temporary pre-resampling population. A history that does not terminate within
`spatial_max_steps`, or a complete scan that exceeds
`spatial_work_item_limit` particle updates, fails instead of returning a
truncated coefficient. The retained history-by-plane sample matrix is capped
at one million values.

The domain contains `spatial_bins` equally spaced interior planes. Each
forward crossing contributes positive statistical weight and each backward
crossing contributes negative weight, producing net electron flux per
injected electron. Bin indices are zero-based;
`spatial_fit_begin_bin` is inclusive and `spatial_fit_end_bin` is exclusive.
An omitted/zero end selects all bins through the final plane. At least three
fit bins are required.

The spatial effective Townsend coefficient is the slope of
`log(net electron flux)` versus distance. Independent history blocks provide
its standard error, while the aggregate fit supplies R². Every fit plane and
every uncertainty block must have positive net flux.
`spatial_min_r_squared` optionally makes fit quality an enforced acceptance
gate. The companion `spatial_profile_file` records every plane's distance,
net crossings, standard error, and fit-selection flag.

This experiment represents the linear, prescribed-field steady response to a
continuous source by superposition of independent histories. It does not
include space-charge distortion, electrode sheaths, or secondary wall
processes. Those belong in an imported-geometry device simulation.

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
rate fitted to `log(total electron weight)`. Aggregate ionization,
attachment, and net creation rates are reported with weight-aware counting
uncertainties. When the conventional electron drift is positive, the output
includes both `growth rate / flux drift velocity` and the rate-balance
effective coefficient `(ionization rate - attachment rate) / flux drift
velocity`. These are clearly named Townsend approximations.
`diffusion_available` and
`townsend_available` distinguish unavailable quantities from numerical zero;
diffusion columns are empty in branching mode because resampling breaks the
independent endpoint-lineage estimator.

When the spatial experiment is enabled, the main CSV also reports the
spatial flux coefficient, its history-block standard error, fit R², completed
history count, maximum active history population, and actual particle-update
count.

The reduced-mobility unit is `1 / (V m s)`. One Townsend is
`1e-21 V m^2`.

## Comparison with reference coefficients

Reference data remain local, just like imported cross sections. Place the
measured or evaluated coefficients in a strict CSV with
`reduced_field_td` and one or more observable columns, then describe the
source and acceptance contract in a `.swarm-reference` manifest:

```ini
[reference]
swarm_reference_version = 2
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
neutral_temperature_k = 300
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
`population_model`. Version 2 additionally requires
`neutral_temperature_k` and rejects simulation rows at a different
temperature. Version 1 remains readable for legacy stationary-neutral
contracts but carries no temperature compatibility check.
`coefficient_convention` is recorded explicitly because
bulk and flux swarm coefficients are not interchangeable; the current
fixed-population mean-velocity estimator is described as
`flux_fixed_population`.

The JSON records the reference metadata, acceptance rule, residuals,
uncertainties, per-value decisions, unmatched simulated fields, and SHA-256
of the simulation CSV, reference CSV, and manifest. Exit status is zero when
all values pass, one when a valid comparison misses its criteria, and two for
an invalid or ambiguous input. Existing reports require `--overwrite`.

## Serialized validation campaigns

A physical cross-section package is not validated by one stochastic run.
Use `run_swarm_campaign.py` to execute two or more resolution studies
serially and test numerical convergence against a designated finest run. If
`reference_manifest` is present, every result is also compared with that same
independent reference contract:

```ini
[campaign]
swarm_campaign_version = 1
campaign_id = laboratory.argon.transport
campaign_version = 2026-01
provenance = local validation campaign and source-control revision
retrieved = 2026-01-15
reference_manifest = argon.swarm-reference
run_order = baseline, timestep_refined, particle_refined
reference_run = particle_refined
field_absolute_tolerance_td = 1e-12
field_relative_tolerance = 1e-12

[run.baseline]
config_file = argon-baseline.swarm
result_file = results/argon-baseline.csv

[run.timestep_refined]
config_file = argon-timestep-refined.swarm
result_file = results/argon-timestep-refined.csv

[run.particle_refined]
config_file = argon-particle-refined.swarm
result_file = results/argon-particle-refined.csv

[observable.drift]
simulation_column = electron_drift_velocity_m_s
uncertainty_column = mean_velocity_x_standard_error_m_s
relative_tolerance = 0.03
absolute_tolerance = 0
uncertainty_multiplier = 2

[observable.mean_energy]
simulation_column = mean_energy_ev
uncertainty_column = mean_energy_standard_error_ev
relative_tolerance = 0.03
absolute_tolerance = 0
uncertainty_multiplier = 2
```

Each `result_file` must be the output selected by its corresponding swarm
configuration. Run the campaign with:

```bash
python3 scripts/run_swarm_campaign.py argon.swarm-campaign \
  --swarm-executable build/aurorapic_swarm \
  --output argon-campaign.json
```

Runs are launched one at a time with `OMP_NUM_THREADS=1`,
`OMP_DYNAMIC=FALSE`, and `OMP_MAX_ACTIVE_LEVELS=1`; a manifest is limited to
16 runs. Existing result, comparison, or aggregate report files are rejected
unless `--overwrite` is explicit. Every run is first evaluated by
`compare_swarm.py` when a reference manifest is supplied. The campaign always
requires identical gas, dataset
ID/version, population model, collision-model signature, and E/N points
across resolutions. It applies the same uncertainty-aware residual form to
each non-reference run relative to `reference_run`.

The aggregate JSON records SHA-256 values for the campaign manifest, every
configuration, and every result. When applicable it also hashes every per-run
reference report and the reference manifest, plus the comparator. Captured
command output is retained per run. Per-run reference reports are stored in
`<report-stem>.artifacts/`. Exit status is zero for a full pass, one for a
well-formed campaign that misses reference or convergence criteria, and two
for invalid inputs or failed simulations.

For an acquisition-stage numerical study, `reference_manifest` may be
omitted. The runner then skips external comparisons, records
`external_reference_available = false` and
`reference_validation_passed = null`, and evaluates only convergence against
`reference_run`. The report carries an explicit claim boundary that such a
pass does not validate the gas data or physical transport model. Adding an
independent reference later restores the combined external-reference and
convergence gate without changing the run schema.

The runner deliberately performs no download and assigns no license. Keep
the original gas package, importer audit, reference data, and their terms
beside the campaign. A convergence pass demonstrates stability under the
declared refinements; it does not prove the collision data or model are
physically correct.

## Model boundary

The fixed mode intentionally excludes multiplication. The
`branching_resampled` mode provides bounded electron-impact avalanche
multiplication, but its growth-over-flux-drift result is not a spatial
steady-state bulk Townsend coefficient. The optional spatial experiment
provides an independent steady-flux coefficient over a declared finite fit
range, but it is still a prescribed-field, linear swarm benchmark rather than
a self-consistent discharge. Attachment removes electron weight
according to its tabulated channel, but the homogeneous runner does not track
the negative-ion product. It does not yet include detachment, photoionization,
space charge, or a self-consistent discharge-level Townsend experiment.
Ionization uses the engine's current equal-sharing energy model. A positive
`neutral_temperature` with an SI gas package activates the same bounded
Maxwellian neutral sampling used by imported geometry; zero is the default.
Elastic scattering is isotropic
unless the gas package explicitly supplies a validated energy-dependent
Henyey-Greenstein mean-cosine table; configured neutral mass remains active
in either elastic recoil path.

Consequently, this benchmark can validate the current MCC implementation's
drift, mean energy, fixed-mode diffusion trend, collision rates, and bounded
transient avalanche growth. It cannot yet claim high-accuracy transport for
datasets requiring a full differential angular cross section beyond the
mean-cosine phase-function approximation, non-equal ionization energy sharing,
neutral bulk flow, or gas-specific spatial Townsend accuracy without
independent measured or evaluated validation.

## Production study checklist

1. Use one complete and internally consistent collision set.
2. Compare against independent measured or evaluated swarm coefficients
   using a traceable `.swarm-reference` contract.
3. Use a serialized `.swarm-campaign` to repeat with smaller timesteps and
   higher particle counts.
4. Extend burn-in until block means no longer show a transient trend.
5. Confirm the maximum observed energy remains comfortably below
   `max_energy_ev`.
6. Record the imported package and `audit.json` alongside the CSV.
7. Do not tune cross sections against device results before the homogeneous
   swarm comparison is understood.
8. For branching runs, repeat with larger computational populations and verify
   both the temporal growth rate and its block uncertainty stabilize.
9. For spatial runs, vary the source-to-fit distance, fit range, plane count,
   history count, and downstream length; require coefficient stability rather
   than accepting one high R² value in isolation.
