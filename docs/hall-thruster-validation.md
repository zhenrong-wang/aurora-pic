# Hall-effect-thruster verification and validation

AuroraPIC will use a staged Hall-effect-thruster (HET) program rather than
claiming that one visually plausible discharge proves the solver. The first
device-level target is the public LANDMARK family of electrostatic PIC
benchmarks. Experimental HET measurements become useful only after the
corresponding geometry, neutral, material, magnetic-field, cathode, and
facility assumptions can be represented and disclosed.

This document distinguishes three activities:

- **verification** checks algorithms against analytic motion, conservation,
  manufactured fields, and independently implemented code benchmarks;
- **validation** compares a fully specified model with experimental
  measurements and their uncertainty;
- **calibration** estimates uncertain model parameters from a declared subset
  of measurements, which must then be separated from validation data.

Agreement with another PIC code is verification, not experimental validation.
The word "validated" must not be used for a HET result until the experimental
configuration and an uncertainty-aware acceptance test are both recorded.

## Public reference cases

### LANDMARK axial-azimuthal benchmark

Charoy et al. published a collisionless 2D3V axial-azimuthal benchmark for
partially magnetized E x B plasmas in *Plasma Sources Science and Technology*
28, 105010 (2019), DOI
[10.1088/1361-6595/ab46c5](https://doi.org/10.1088/1361-6595/ab46c5).
An [open accepted manuscript](https://www.osti.gov/servlets/purl/1572876)
contains the complete model. Seven independently developed PIC codes compared
time- and azimuth-averaged profiles and the dominant instability. This is a
verification benchmark representative of Hall physics, not an exact model of
a manufactured thruster.

The pinned published inputs are:

| Quantity | Value |
| --- | ---: |
| Axial by azimuthal domain | `2.5 cm` by `1.28 cm` |
| Structured cells | `500` by `256` |
| Cell size | `50 um` square |
| Time step and final time | `5e-12 s`, `20 us` (`4,000,000` steps) |
| Discharge voltage | `200 V` |
| Initial plasma density | `5e16 m^-3` |
| Initial electron / Xe+ temperature | `10 eV` / `0.5 eV` |
| Initial particles/cell/species | Case 1: `150`; Case 2: `75`; Case 3: `300` |
| Radial magnetic field | prescribed axial Gaussian-like profile; `6 mT` at the anode, `10 mT` at `x=0.75 cm`, `1 mT` downstream |
| Pair source | prescribed cosine profile from `x=0.25 cm` to `1.0 cm`, peak `5.23e23 m^-3 s^-1` |
| Boundaries | periodic azimuth; absorbing axial particle boundaries; `200 V` anode and cathode-side emission model |

Electron-Xe collisions and neutral transport are deliberately omitted. The
source creates colocated electron-Xe+ pairs at `10 eV` and `0.5 eV`. Cathode
electrons are emitted from `x=2.4 cm` according to current continuity, and the
published potential correction holds the mean emission-plane voltage drop at
200 V. Those source, current-control, and potential-correction rules are part
of the benchmark; replacing them with generic inflow is a different case.

The primary comparisons are time- and azimuth-averaged axial electric field,
ion density, and electron temperature, followed by azimuthal mode wavelength,
frequency, propagation, and electron cross-field current. The paper reports
remaining inter-code differences of about 5% and substantial sensitivity to
random seeds. Five percent is therefore context for designing an ensemble
acceptance interval, not a tolerance to copy without reproducing the published
averaging and convergence study.

[University of Michigan Deep Blue Data](https://deepblue.lib.umich.edu/data/concern/data_sets/vq27zp39n?locale=en)
also publishes WarpX output for the 75-particle/cell Case 2a run: 20 us,
four million iterations, output every 5,000 iterations, and approximately
32 GB in AMReX plotfile form. This is a valuable independent reference corpus,
but it must be downloaded explicitly outside the source tree and must never be
an ordinary CI dependency.

### LANDMARK radial-azimuthal benchmark

Villafana et al. published a complementary 2D3V radial-azimuthal benchmark in
*Plasma Sources Science and Technology* 30, 075002 (2021), DOI
[10.1088/1361-6595/ac0a4a](https://doi.org/10.1088/1361-6595/ac0a4a).
The [open manuscript and record](https://www.osti.gov/pages/biblio/1826078)
describe seven independent PIC implementations, plasma-wall interaction,
electron cyclotron drift instability (ECDI), modified two-stream instability
(MTSI), convergence, and supplementary verification outputs.

Its principal inputs are a `1.28 cm` by `1.28 cm` periodic-azimuthal,
grounded-radial-wall domain with `256` by `256` cells, `50 um` spacing,
`1.5e-11 s` timestep, `30 us` final time, xenon, `5e16 m^-3`, `10 eV`
electrons, `0.5 eV` ions, and 100 initial particles/cell/species. It imposes a
uniform `20 mT` radial magnetic field, a `10 kV/m` out-of-plane axial electric
field, and a virtual `1 cm` axial domain. It is collisionless, absorbs
particles at the grounded radial walls, injects a prescribed radial pair
source, and uses a virtual-axial replacement model.

This benchmark is the stronger wall/instability test, but not the first HET
implementation target: AuroraPIC must first support periodic topology combined
with radial walls, three-component prescribed electric fields, and the exact
virtual-axial replacement/source contract.

### Experimental validation targets

NASA publishes useful HET measurements through the Technical Reports Server.
The [extended HERMeS performance characterization](https://ntrs.nasa.gov/citations/20210024578)
covers a 12.5 kW magnetically shielded thruster over a broad throttle range.
The [HERMeS ion-velocity study](https://ntrs.nasa.gov/citations/20180006879)
uses laser-induced fluorescence in the discharge channel and near field
specifically to support model validation, and the
[facility-effect campaign](https://ntrs.nasa.gov/citations/20170000953)
reports thrust, probe, spectrometer, camera, and plume comparisons over
background pressure.

These data can eventually test thrust, discharge current, efficiency,
oscillation spectra, ion velocity, plume current density, and plume
divergence. They are not presently a reproducible AuroraPIC case: public
papers do not constitute one versioned package containing every internal
dimension, magnetic-field map, cathode boundary condition, wall property,
neutral inlet, and facility condition needed for an unambiguous full PIC
input deck. Facility pressure must be treated as a modeled condition, not
silently ignored.

## PIC, fluid, and hybrid roles

Full kinetic PIC is appropriate when electron distribution functions,
sheaths, ECDI/MTSI, anomalous cross-field transport, and ion velocity
distributions are the quantities of interest. It is also expensive because an
explicit electron PIC model must resolve electron plasma and Debye scales.

Whole-device engineering studies often use multifluid or hybrid models:
kinetic ions and neutrals with fluid electrons, or related reduced models.
For a weakly collisional, partially magnetized HET, this is generally more
appropriate language than ideal single-fluid MHD. Such models are far cheaper
for breathing modes, plume coupling, magnetic optimization, lifetime sweeps,
and operating maps, but require closures for electron transport and wall
physics.

AuroraPIC should keep one geometry/mesh, field-map, gas-data, boundary,
diagnostic, and provenance layer, then allow PIC and future hybrid/fluid
backends to share cases. A reduced model does not replace kinetic verification;
overlapping cases should compare moments, currents, and conservation between
the models.

## AuroraPIC readiness and gaps

| Capability | Current state | Required action |
| --- | --- | --- |
| Electrostatic 2D3V PIC and Boris push | Available on structured and imported 2D meshes | Preserve analytic Larmor and E x B drift gates |
| Prescribed magnetic field | Uniform vectors and strict one-coordinate tabulated profiles are available across structured 2D/3D and imported 2D | Add profile provenance fingerprints and arbitrary sampled-map import |
| Mixed topology | Structured 2D automatically uses a direct spectral-tridiagonal Poisson solve for either periodic/Dirichlet orientation; mixed-radix FFT and Bluestein paths cover composite and prime periodic sizes | Measure the serial production grid, then distribute the transform and axial mode solves with MPI before full LANDMARK campaigns |
| Volumetric pair source | Structured 2D has normalized profiles, explicit extrusion depth, analytic peak-volumetric-to-total conversion, SI eV thermal loading, fractional accumulation, diagnostics, restart, and a versioned reduced LANDMARK manifest | Qualify source statistics at production population and duration |
| Cathode/current control | Selectable generic cumulative regulation and published timestep-local electron-minus-ion loss injection are available with explicit flux/reverse-demand diagnostics and checkpoint v9; the exact-control workstation pilot passed | Qualify long-duration response and cathode-temperature sensitivity |
| Cathode potential correction | Selectable gauge and affine internal-plane corrections are available; affine mode preserves the anode and changes the axial field analytically; the exact-control workstation pilot passed | Compare converged profiles against checksum-pinned published data |
| Radial benchmark virtual axis | Not available | Implement bounded virtual-axis replacement and audit its energy/particle flux |
| HET diagnostics | Structured 2D emits transverse field/species profiles, density-weighted three-velocity moments, all current components, trapezoidal time averages, complex periodic-axis Fourier histories, checksum-pinned reference comparisons, seeded ensemble statistics, a guarded population/duration convergence workflow, and checkpoint-chained horizon stages | Execute the guarded horizon ladder, add long-run segment aggregation, and qualify with native reference data |
| Xenon material data | No authoritative bundled package | Keep LANDMARK collisionless; separately provenance and validate Xe collision/wall data for real devices |
| Scale-out runtime | Serial/OpenMP only; a bounded, provenance-pinned runtime qualifier measures a selected host without permitting a production launch | Add MPI domain decomposition before production-size LANDMARK runs |
| High-volume output | VTK/XML, CSV, and text restart | Add openPMD/HDF5 or equivalent parallel, chunked output before large campaigns |
| Real HET geometry | Tagged planar Gmsh import exists | Add spatial field maps, axisymmetric/3D geometry, dielectric/material walls, neutral flow, and cathode/facility models |

The axial-azimuthal benchmark is rectangular and therefore does not itself
verify external geometry import. Geometry and HET physics are independent
gates: both the LANDMARK benchmark and an imported-mesh manufactured-field
case must pass before attempting a real thruster.

## Execution ladder and acceptance

HET work proceeds in this order:

1. **H0, analytic and manufactured verification:** bounded single-particle
   gyro-orbit and E x B drift tests, spatial magnetic-profile interpolation,
   periodic wrapping, absorbing axial flux, and colocated neutral pair-source
   charge conservation.
2. **H1, reduced LANDMARK integration:** the exact physical contract on a
   deliberately coarse grid and short time, used only for invariants,
   restart determinism, diagnostics, and resource behavior. It cannot make a
   physics-agreement claim.
3. **H2, axial-azimuthal Case 2:** convergence tiers in grid, timestep,
   particles/cell, averaging window, and random-seed ensemble; compare axial
   profiles and instability spectra with the paper and public WarpX data.
4. **H3, radial-azimuthal case:** verify grounded-wall losses, ECDI/MTSI
   spectra, source balance, virtual-axis behavior, and statistical
   convergence.
5. **H4, imported HET precursor:** use a public or publishable geometry and
   measured magnetic map; verify mesh refinement, wall fluxes, neutral
   balance, and integrated current before comparison to experiment.
6. **H5, experimental validation:** freeze calibration and validation
   subsets, propagate measurement/input uncertainty, and compare performance,
   ion velocity, oscillation, and plume observables at multiple operating
   points.

Every H2-H5 report must pin source revision, compiler/backend, case manifest,
input-data hashes and licenses, random seeds, mesh and timestep, particles per
cell, averaging interval, convergence results, and diagnostic scripts. A fixed
seed is required for restart regression; multiple seeds are required for a
turbulent physics claim.

## External reference comparison

AuroraPIC does not redistribute the approximately 29 GB WarpX corpus or
derived paper data. The public repository may also be inaccessible from some
hosts. Acquisition is therefore separate from normalization. The committed
`examples/hall_landmark_case2.sources` registry distinguishes the original
500×256 LANDMARK supplement from the CC0 WarpX 512×256 AMReX corpus and pins
their landing records, DOI, license, artifact identity, and acquisition
method.

The non-downloading source tool first writes an acquisition plan:

```sh
python3 scripts/lock_hall_source.py \
  examples/hall_landmark_case2.sources \
  --source warpx_deepblue \
  --output local/warpx-acquisition-plan.json
```

Deep Blue directs this 29 GB artifact to external Globus transfer. After that
transfer, the same tool stream-hashes the local archive. Files above 64 MiB
require the explicit
`I_UNDERSTAND_THIS_MAY_HASH_A_VERY_LARGE_FILE` acknowledgement. Supplying a
checksum independently obtained from repository metadata records
`repository_checksum_verified = true`; omitting it records that only local
post-acquisition byte identity has been established. The tool never downloads,
extracts, or loads the full archive into memory.

An operator then extracts the required profile/mode tables with a recorded
procedure and creates a local `.hall-source` table lock that pins those exact
raw bytes and their column/unit map:

```ini
[source]
hall_source_version = 1
case_id = landmark-axial-azimuthal-2019
case_variant = case-2-multicode-supplement
case_manifest_sha256 = <lowercase SHA-256>
profile_file = published-profiles.csv
profile_sha256 = <lowercase SHA-256>
mode_file = published-modes.csv
mode_sha256 = <lowercase SHA-256>
source_url = <direct artifact or repository record URL>
source_artifact_id = <repository identifier or supplement filename>
provenance = <how the two raw tables were obtained>
citation = Charoy et al. 2019, doi:10.1088/1361-6595/ab46c5
retrieved = YYYY-MM-DD
license = <terms applying to the source tables>

[profile]
coordinate_column = x_cm
coordinate_scale_to_m = 0.01
electric_field_columns = code_a_ex,code_b_ex
electric_field_scale_to_v_m = 1
ion_density_columns = code_a_ni_cm3,code_b_ni_cm3
ion_density_scale_to_m3 = 1e6
electron_temperature_columns = code_a_te,code_b_te
electron_temperature_scale_to_ev = 1

[mode]
mode_column = azimuthal_mode
frequency_columns = code_a_frequency_khz,code_b_frequency_khz
frequency_scale_to_hz = 1000
comparison_mode = 16

[acceptance]
coordinate_absolute_tolerance_m = 1e-12
relative_tolerance = 0.05
uncertainty_multiplier = 1
```

Column lists contain independent published codes or replicates. Normalization
uses their midpoint as the reference and half their range as uncertainty.
One column is valid and produces zero range uncertainty. Values are converted
only by the declared multiplicative scales. Source row order and coordinates
are preserved; the normalizer refuses duplicate/non-increasing coordinates
and performs no interpolation or extrapolation:

```sh
python3 scripts/normalize_hall_reference.py \
  local/case2.hall-source \
  --case-manifest examples/hall_landmark_axial_azimuthal.case \
  --output-dir local/case2-reference
```

The command verifies the raw-table and case-manifest hashes before writing an
atomic output directory containing canonical profile/mode CSVs,
`reference.hall-reference`, and `normalization.json`. Existing output is never
overwritten. The audit pins every input and output hash, source artifact
identity, the `midpoint_and_half_range` envelope method, and the
`native_no_interpolation` coordinate policy.

The generated `.hall-reference` manifest pins the exact normalized profile and
mode CSV bytes:

```ini
[reference]
hall_reference_version = 1
case_id = landmark-axial-azimuthal-2019
case_variant = case-2a-warpx
case_manifest_sha256 = <lowercase SHA-256>
profile_data_file = landmark-case2a-profiles.csv
profile_data_sha256 = <lowercase SHA-256>
mode_data_file = landmark-case2a-modes.csv
mode_data_sha256 = <lowercase SHA-256>
profile_axis = x
mode_axis = y
coordinate_column = coordinate_m
coordinate_absolute_tolerance = 1e-12
provenance = conversion procedure and source artifact identity
citation = Charoy et al. 2019 and public WarpX dataset
retrieved = YYYY-MM-DD
license = terms applying to the local reference files

[profile.axial_field]
simulation_source = field
simulation_column = electric_x
reference_column = electric_x_v_m
reference_uncertainty_column = electric_x_uncertainty_v_m
relative_tolerance = 0.05
absolute_tolerance = 0
uncertainty_multiplier = 2

[profile.ion_density]
simulation_source = species
simulation_species = ions
simulation_column = number_density
reference_column = ion_density_m3
relative_tolerance = 0.05

[profile.electron_temperature]
simulation_source = species
simulation_species = electrons
simulation_column = temperature_ev
reference_column = electron_temperature_ev
relative_tolerance = 0.05

[mode.dominant_frequency]
simulation_quantity = electric_y
# Omit simulation_species for a mesh-field quantity.
mode = 16
metric = frequency_hz
reference_column = frequency_hz
reference_uncertainty_column = frequency_uncertainty_hz
relative_tolerance = 0.05
uncertainty_multiplier = 2
```

The reference profile CSV has exactly one row per reference coordinate.
Coordinates must match exactly one AuroraPIC profile point within the declared
absolute tolerance; the comparator never silently interpolates or
extrapolates. The mode reference CSV has exactly one row per integer `mode`.
Each `[profile.*]` or `[mode.*]` section declares its own acceptance rule:

```text
abs(simulation - reference)
<= absolute_tolerance
 + relative_tolerance * abs(reference)
 + uncertainty_multiplier * reference_uncertainty
```

`scripts/compare_hall.py` verifies both reference hashes and the simulation
case ID, verifies that field and species files share one nonzero averaging
window, and writes an atomic JSON report containing every input hash and
residual. For mode histories it unwraps the complex coefficient phase and
supports `frequency_hz`, `signed_frequency_hz`, `phase_velocity_m_s`,
`growth_rate_s`, and `rms_amplitude`. At least three unique samples with
nonzero amplitude are required. Sampling must be frequent enough that the
phase advance between samples remains below pi; the comparator cannot recover
a frequency already aliased by the output cadence.

```sh
python3 scripts/compare_hall.py production-output \
  local/case2a.hall-reference \
  --case-manifest examples/hall_landmark_axial_azimuthal.case \
  --output production-output/hall-comparison.json
```

Exit status is zero for agreement, one for a valid comparison outside its
criteria, and two for malformed, ambiguous, unpinned, or mismatched inputs.
Existing reports are never overwritten without `--overwrite`. The bounded
regression uses synthetic data only; it tests the comparison machinery, not a
LANDMARK result.

### Published Figure 6 profile screening

The public accepted manuscript contains the seven-code axial-profile envelope
for Case 2 in Figure 6. It is averaged from 16 to 20 microseconds and covers
axial electric field, ion density, and electron temperature. This is useful
before the native supplement or the 29 GB WarpX archive is available, but a
digitized plot is not native validation data.

AuroraPIC therefore treats this source as a version 2
`digitized_profile_screening` contract. Its report always records
`physics_claim = none`, both averaging windows, and whether those windows
match. It contains no fabricated mode reference. Version 1 profile-plus-mode
validation contracts remain strict and unchanged.

Keep the manuscript and derived tables outside the repository. After obtaining
the public accepted manuscript, reproduce the checksum-pinned screening data:

```sh
python3 scripts/digitize_charoy_figure6.py \
  /external/reference/charoy-2019-accepted-manuscript.pdf \
  --case-manifest examples/hall_landmark_axial_azimuthal.case \
  --output-dir /external/reference/charoy-figure6-screening
```

Then compare an existing resolved run:

```sh
python3 scripts/compare_hall.py \
  /external/run/output \
  /external/reference/charoy-figure6-screening/reference.hall-reference \
  --case-manifest examples/hall_landmark_axial_azimuthal.case \
  --output /external/reference/charoy-figure6-screening/comparison.json
```

The digitizer verifies the source PDF hash, page-vector structure, curve count,
and color ordering before emitting profiles. It samples the seven curves on
the diagnostic 0.2 mm grid and records their midpoint and half-range. That
half-range measures inter-code spread; it is not experimental uncertainty.
Use the original authors' numerical supplement, or the checksum-pinned WarpX
archive, for a validation claim. Screening-only reports are not accepted by
the seeded ensemble aggregator.

### Guarded time-horizon extension

The workstation duration-convergence ladder ends at 10,000 steps, or 50 ns.
That is useful for early-transient sensitivity but cannot approach the
published 16–20 microsecond averaging window. Time-horizon work therefore uses
one checkpoint-chained stage at a time. A stage is prepared only from a
completed serial, single-thread prior run whose final scalar step, physical
time, and checkpoint header agree.

For the completed 5,000-step workstation pilot, prepare the 20,000-step
(100 ns) stage without launching it:

```sh
python3 scripts/prepare_hall_horizon_stage.py \
  examples/hall_landmark_axial_azimuthal.case \
  --prior-deck /external/pilot/case.cfg \
  --prior-output /external/pilot/output \
  --output-dir /external/campaign/horizon-100ns \
  --acknowledge-cost \
  I_UNDERSTAND_THIS_EXTENDS_A_HALL_RUN_FROM_A_PINNED_CHECKPOINT
```

The generated `horizon.json` pins the prior deck, scalar history, and final
checkpoint hashes. It records measured species growth, a four-times growth
projection, 25% capacity headroom, added particle-update cost, and an
11-sample average over the final 20% of the target horizon. The generated deck
remains serial with one thread and is still protected by the CLI large-run
acknowledgement. The preparer never launches it.

After execution, integrity checks and published-profile screening must be
repeated. The stage analyzer binds both comparison reports, verifies restart
and final checkpoints, checks population capacity and potential correction,
and requires every profile L2/Linf error to improve:

```sh
python3 scripts/analyze_hall_horizon_stage.py \
  /external/campaign/horizon-100ns/horizon.json \
  --current-comparison /external/campaign/horizon-100ns/comparison.json \
  --prior-comparison /external/pilot/comparison.json \
  --report /external/campaign/horizon-100ns/credibility.json
```

It separately blocks the next extension when cumulative unserved reverse
charge grows in the one-way cathode controller. Prepare another stage only if
capacity remains safe, the profile trend is credible, and that controller gate
is resolved. Every stage retains `physics_claim = none`;
reaching 20 microseconds, population/time-step/grid convergence, multiple
seeds, and native-reference agreement remain separate gates.

The completed 100 ns continuation is recorded in
`benchmarks/hall/landmark-horizon-100ns-20260730-seed24680.json`. It ran from
the pinned 5,000-step checkpoint to step 20,000 in 850.56 seconds on one
low-priority CPU core, used 32.6 MiB peak resident memory, and passed numerical
integrity checks. Relative L2 error against the digitized Figure 6 envelope
fell by 48.7% for axial field, 5.0% for ion density, and 21.8% for electron
temperature. This is encouraging time evolution, not reference agreement.
The cumulative one-way cathode demand that could not be served increased from
158 to 849 macro-electron charges, or from 1.25% to 4.53% of cumulative
emission. The analyzer therefore holds a 400 ns extension for review of the
reduced early-time cathode/anode flux balance.

Checkpoint v9 adds restart-safe, per-update observability to distinguish that
cumulative diagnostic from a controller instability. The guarded
1,000-update probe recorded in
`benchmarks/hall/landmark-controller-probe-v9-20260730-seed24680.json`
covered steps 20,001 through 21,000 (5 ns): reverse demand occurred in 43
updates (4.3%), totaled 26 macro-electron charges, and never exceeded one
macro-electron charge in an update. The run completed on one pinned,
low-priority CPU core in 55.08 seconds with 35.7 MiB peak resident memory.
This rules out runaway reverse demand over the observed interval and instead
identifies a sparse, population-quantized edge condition. It does not prove
population convergence or justify bypassing the controller gate. The next
physics-facing test is a paired population-sensitivity probe with the same
physical inputs, timestep, seed policy, and observation window.

## Seeded ensemble comparison

A turbulent Hall result cannot be accepted from one favorable seed.
`prepare_hall_ensemble.py` requires 3–64 unique unsigned 32-bit seeds and
atomically writes one deck per seed plus `ensemble.json`. It never launches
the decks, and workstation/production generation retains the same explicit
cost acknowledgement as a single deck:

```sh
python3 scripts/prepare_hall_ensemble.py \
  examples/hall_landmark_axial_azimuthal.case \
  --tier production \
  --seeds 104729,130363,155921 \
  --output-dir campaign/case2-ensemble \
  --acknowledge-cost I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_RUN
```

Each completed run is compared separately. Passing `--runtime-config` binds
the seed and exact deck hash into its comparison report:

```sh
python3 scripts/compare_hall.py \
  campaign/case2-ensemble/results/seed_104729 \
  local/case2-reference/reference.hall-reference \
  --case-manifest examples/hall_landmark_axial_azimuthal.case \
  --runtime-config campaign/case2-ensemble/case_seed_104729.cfg \
  --output campaign/case2-ensemble/comparisons/seed_104729.json
```

After all declared reports exist:

```sh
python3 scripts/aggregate_hall_ensemble.py \
  campaign/case2-ensemble/ensemble.json \
  --output campaign/case2-ensemble/ensemble-comparison.json
```

Aggregation rejects changed decks, seed/path mismatches, inconsistent
reference hashes or averaging windows, different comparison shapes, and
identical field/species/mode artifact hashes reused across nominally
independent seeds. For every profile and mode observable it computes the
ensemble mean, sample standard deviation, standard error, and two-sided 95%
Student-t interval. The conservative error bound
`abs(ensemble_mean - reference) + confidence_half_width` must fit inside the
original comparison threshold, and at least two thirds of individual reports
must pass by default. A passing reduced tier records
`physics_claim_eligible = false`; production eligibility still does not by
itself constitute experimental validation.

## Workstation convergence campaign

The same-seed workstation convergence contract separates population and
duration sensitivity while holding geometry, timestep, fields, sources, and
control models fixed. Its five unique stages are:

| Axis | Coarse | Baseline | Fine |
|---|---:|---:|---:|
| particles/cell/species | 8 | 16 | 32 |
| steps | 2,500 | 5,000 | 10,000 |

The shared baseline is run once, so the complete plan contains five runs and
7.68 billion lower-bound initial-particle updates. Each stage retains 11
samples over its final 20% and scales particle capacity with population.
Planning requires an explicit acknowledgement and never launches a solver:

```sh
python3 scripts/prepare_hall_convergence.py \
  examples/hall_landmark_axial_azimuthal.case \
  --output-dir campaign/case2-convergence \
  --acknowledge-cost \
  I_UNDERSTAND_THIS_IS_AN_OPT_IN_HALL_CONVERGENCE_PLAN
```

Every generated deck remains independently protected by the CLI large-run
acknowledgement. The population stages also preserve the checkpoint-v9
controller observability contract from step zero. After separately executing
all five decks, analyze their time-averaged axial field/species profiles,
time-averaged azimuthal spectra, and cathode reverse-demand scaling:

```sh
python3 scripts/analyze_hall_convergence.py \
  campaign/case2-convergence/convergence.json \
  --report campaign/case2-convergence/convergence-report.json
```

For each observable, the report records absolute and baseline-normalized L2
and Linf changes from coarse to baseline and fine to baseline. Acceptance
requires the fine change to satisfy the pinned tolerances and not exceed the
coarse change. The controller population gate compares represented reverse
charge per update and the maximum represented reverse impulse, rather than
raw macro-particle counts. The 32-particle/cell stage may not exceed 1.25
times the baseline reverse charge per update, its maximum represented impulse
must fall to at most 0.75 times baseline, and its largest reverse request must
remain at or below 1.5 macroparticles. Reverse-event frequency is reported but
is deliberately not required to decrease: smaller macro weights can resolve
more sub-baseline events even while their represented physical charge falls.
All five runs require complete v9 diagnostics; a legacy restart with partial
coverage is rejected. This is a fixed-grid, same-seed sensitivity gate with
`physics_claim = none`; grid, timestep, seed-ensemble, and external-reference
agreement remain separate requirements.

## Resource policy

No published-scale HET case runs in `scripts/verify.sh`, CTest, or ordinary CI.
H0 and H1 tests must default to one process and one thread and have explicit
particle, cell, step, memory, and output limits. H2-H5 are opt-in production
campaigns with a preflight resource estimate, explicit output directory,
checkpoint cadence, storage quota, and wall-time limit. Reference datasets
remain external and checksum-pinned; they are never vendored into the Git
repository.

Run the non-launching resource gate before creating a production deck:

```sh
python3 scripts/preflight_hall.py \
  examples/hall_landmark_axial_azimuthal.case \
  --tier production \
  --report hall-case2a-preflight.json
```

The manifest defines three monotonic tiers:

| Tier | Cells | Initial particles/species | Steps | Lower-bound updates | Generation |
|---|---:|---:|---:|---:|---|
| `micro` | 32×16 | 2,048 | 200 | 819,200 | unlocked |
| `workstation` | 125×64 | 128,000 | 5,000 | 1.28 billion | cost acknowledgement |
| `production` | 500×256 | 9.6 million | 4 million | 76.8 trillion | cost acknowledgement |

Every tier retains the same physical domain, timestep, prescribed magnetic
field, source, cathode controller, and potential-reference model. The reduced
tiers make no physics claim; they test integration, population growth,
diagnostics, and scaling only.

The production Case 2 contract contains 128,000 cells,
19.2 million initial macro-particles, and at least 76.8 trillion particle
updates. Its actual configured hard capacity of 80 million particles per
species, a 96-byte particle assumption, two retained checkpoints, 161 resolved
samples over 16–20 microseconds, and modes 0–128 estimate about 14.32 GiB of
resident memory and 28.76 GiB of diagnostics plus checkpoints. These are
capacity-based planning estimates, not expected live-population measurements;
sources can increase population, while VTK, particle dumps, external data,
replication, and temporary files are excluded. At an explicitly supplied
measured rate of 100 million particle updates/s, the lower-bound push time
alone is 768,000 seconds (8.9 days).

### Runtime qualification before execution

Case 2 is still the smallest suitable first full-PIC HET target with a public
multi-code axial-azimuthal comparison. Cases 1 and 3 use two and four times the
Case 2 initial particle population. The radial-azimuthal LANDMARK case reduces
the initial-population update count to about 26.2 trillion, but it is not an
immediate substitute: AuroraPIC does not yet implement its virtual-axial
replacement model. The one-dimensional LANDMARK Hall case is a fluid/hybrid
benchmark and therefore cannot verify this kinetic PIC path.

Before any longer local pilot, run a bounded workstation-resolution timing
slice:

```sh
taskset -c 0 nice -n 19 python3 scripts/qualify_hall_runtime.py \
  build/aurorapic_cli \
  examples/hall_landmark_axial_azimuthal.case \
  --tier workstation \
  --steps 40 \
  --max-initial-updates 11000000 \
  --timeout-seconds 60 \
  --acknowledge-cost I_UNDERSTAND_THIS_IS_A_BOUNDED_HALL_PROBE \
  --report local/hall-runtime-qualification.json
```

The tool supports only the `micro` and `workstation` tiers, forces the serial
one-thread backend, refuses more than 25 million initial-particle updates,
applies a hard timeout, disables resolved/particle/VTK/checkpoint output, and
never passes the CLI large-run acknowledgement. It records the exact case,
deck, and executable hashes, CPU affinity, measured wall time, and
initial-population-only workstation/production projections. Existing reports
are rejected before a simulation is launched.

The projection is a screening estimate, not a scheduling promise. It includes
short-run startup overhead but excludes growth from the pair and cathode
sources. Published Case 2 data show growth from 75 to roughly 290
particles/cell, so a projection based only on the initial 76.8 trillion
updates is an optimistic lower bound for a complete discharge.

The first complete 5,000-step workstation pilot is recorded in
`benchmarks/hall/landmark-workstation-20260729-seed24680.json`. It completed
on one low-priority CPU core in 267.56 seconds with 32.2 MiB peak resident
memory and 121.5 MiB of artifacts. All integration checks passed, but the run
exposed timestep-local cathode and affine-potential differences that motivated
the selectable exact modes implemented afterward. The committed record pins the executable, deck,
manifest, analyzer, scalar, resolved-diagnostic, and checkpoint hashes while
retaining `physics_claim = none`.

The exact-control rerun is recorded in
`benchmarks/hall/landmark-workstation-exact-20260729-seed24680.json`. With the
same seed, grid, initial population, timestep, and duration, it completed in
270.37 seconds on one low-priority CPU core with 32.1 MiB peak resident
memory. The timestep-local controller retained zero negative debt and reduced
the controllable charge residual by 98.5%; its one-way unserved reverse demand
is reported separately. The affine correction changed the averaged axial-field
profile by 13.6 kV/m RMS relative to the earlier gauge run. This is still an
integration comparison with `physics_claim = none`, not published-case
agreement.

The preflight distinguishes the paper's 500 by 256 cells from AuroraPIC's
501 by 256 structured nodes, writes all assumptions and arithmetic to JSON,
returns one when
a declared memory or storage budget is exceeded, never launches AuroraPIC, and
always records `launch_authorized = false`. The current case remains
`reduced_integration_only`. A production-candidate deck can only be generated
after spelling out an explicit cost acknowledgement:

```sh
python3 scripts/prepare_hall_campaign.py \
  examples/hall_landmark_axial_azimuthal.case \
  --tier production \
  --output campaign/case2.cfg \
  --acknowledge-cost I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_RUN
```

This copies the checksum-verified magnetic profile beside the deck and writes
the original 500 by 256 cell, 75-particle-per-cell-per-species, 4-million-step
contract. It starts resolved averaging at step 3.2 million (16 microseconds),
uses one serial thread, disables VTK and particle dumps, and sets an explicit
particle capacity. The tool never runs the generated deck. Measured solver
throughput, adequate capacity from a pilot campaign, MPI decomposition,
convergence tiers, scheduler quotas, and explicit execution authorization
remain separate gates.

The generated deck can be parsed and semantically validated without allocating
its particles or starting a timestep:

```sh
aurorapic_cli --validate-only campaign/case2.cfg
```

The CLI independently blocks configurations above 100 million initial
particle updates. After scheduler and resource review, a deliberate large run
requires `--allow-large-run I_UNDERSTAND_THIS_IS_A_LARGE_RUN`. This second
guard prevents an already-generated workstation or production deck from being
launched accidentally.

The `micro` tier is the only Hall discharge tier exercised by ordinary
verification. After it completes, `scripts/analyze_hall_pilot.py` checks
sampling cadence, finite diagnostics, particle capacity, prescribed pair rate,
cathode emitted-charge/residual/remainder identities, one-way actuator debt,
potential-reference accuracy,
profile-average shape, and complete azimuthal mode coverage. Its report always
records `physics_claim = none`; passing it is an integration result, not
agreement with the published discharge.

The first H0 field slice is complete: a shared strict `coordinate Bx By Bz`
profile supports linear interpolation and full-domain coverage checks across
structured 2D/3D and imported 2D, while preserving uniform-field
compatibility. The committed LANDMARK magnetic profile is exercised only as a
four-step integration smoke and cannot make a discharge claim. Independent
structured 2D axis topology is also complete:
Dirichlet axial plus periodic azimuthal and its transposed orientation have
analytic vacuum-field, charge-conservation, spacing, and default particle
policy regressions. Mixed cases now transform the periodic axis and solve one
complex tridiagonal Dirichlet system per Fourier mode. Composite sizes use a
mixed-radix FFT and prime sizes use Bluestein convolution, avoiding the former
quadratic transform fallback; manufactured discrete-Poisson regressions cover
both axis orientations and non-power-of-two sizes. A field-only regression
also solves the original published 500 axial by 256 azimuthal cell Case 2
grid (501 by 256 AuroraPIC nodes) with vacuum charge and
nonzero axial electrodes; it contains no particles or timesteps and is not a
throughput claim. The implementation remains single-rank and requires measured
production-grid qualification. A generic
structured-2D volumetric source layer is now also
complete: named scheduled sources create equal-weight, opposite-charge pairs
at shared positions drawn from normalized uniform, Gaussian, or sinusoidal
profiles; fixed macro rates, total represented physical rates, and peak
volumetric rates are mutually exclusive; fractional production is
deterministic across steps and restart; and diagnostics include rate,
remainder, effective area, extrusion depth, and injected energy. The Hall
smoke exercises the centered cosine-family profile and pinned peak
volumetric-rate path without claiming to be a resolved discharge.
The versioned `examples/hall_landmark_axial_azimuthal.case` now derives
`2.5104e19 s^-1` from the published peak, normalized profile integral, and an
explicit `1 m` reduced extrusion depth. Its checksum, budgets, runtime linkage,
and no-claim limitations are enforced by `scripts/validate_hall_case.py`.
The same manifest pins timestep-local charge-regulated electron emission at
`x = 2.4 cm`, the published `10 eV` velocity scale, and the affine zero-mean
potential/field correction at that plane. Focused regressions preserve the
generic unequal-weight cumulative and gauge modes while independently
verifying the LANDMARK timestep-local and affine modes plus deterministic
checkpoint continuation.
The manifest also pins the published `10 eV` electron and `0.5 eV` Xe+
temperatures for both initial loading and pair creation; strict SI conversion
and sampled three-component moment regressions guard the contract. The reduced
runtime now emits five bounded axial profile/moment/current samples,
trapezoidal field and density-weighted species averages, and complex azimuthal
mode histories through mode three. Manufactured regressions pin the
normalization and complex Fourier convention, while the manifest and example
smoke pin the artifact schema. These diagnostics enable comparison but do not
constitute one. The external comparator now pins reference/profile/mode/case
hashes, performs uncertainty-aware residual checks, and derives modal
frequency or growth from complex histories; the non-launching preflight pins
resource arithmetic and refuses exceeded budgets. Both are guarded with
synthetic data only. The next H2-enabling slice is execution and trend
analysis of the guarded 100 ns checkpoint stage, followed by long-run
diagnostic segment aggregation and MPI decomposition required for the
published run.
