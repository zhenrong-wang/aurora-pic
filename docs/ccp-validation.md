# Turner helium CCP validation target

AuroraPIC targets the four one-dimensional helium capacitively coupled plasma
(CCP) cases published by Turner et al., *Physics of Plasmas* 20, 013507
(2013), DOI
[10.1063/1.4775084](https://doi.org/10.1063/1.4775084). These are
code-to-code physics benchmarks with quantified statistical agreement, not
experimental validation of a particular reactor.

## Published contract

The model is a 6.7 cm planar gap containing uniform 300 K helium. One
electrode is grounded and the other is driven by

```text
phi(t) = V*sin(2*pi*13.56e6*t),
```

so the applied voltage is zero at time zero. Electrons and singly charged
helium monomer ions are absorbed at both electrodes, with no secondary
emission. The initial electron temperature is 30,000 K and the initial ion
temperature is 300 K.

| Case | Neutral density (m^-3) | Voltage amplitude (V) | Initial plasma density (m^-3) | Cells | Particles/cell/species | Steps/RF period | RF periods |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `9.64e20` | `450` | `2.56e14` | 128 | 512 | 400 | 1280 |
| 2 | `3.21e21` | `200` | `5.12e14` | 256 | 256 | 800 | 5120 |
| 3 | `9.64e21` | `150` | `5.12e14` | 512 | 128 | 1600 | 5120 |
| 4 | `3.21e22` | `120` | `3.84e14` | 512 | 64 | 3200 | 15360 |

The table uses the paper's cell count. AuroraPIC's 1D Dirichlet `nx` is the
number of endpoint-inclusive nodes, so the corresponding configurations use
`nx = cells + 1` (for example, `nx = 129` for Case 1).

The final 32 RF periods are averaged. The primary comparison is the
time-averaged ion-density profile, evaluated using the published mean,
standard-deviation, and chi-squared data. Useful secondary quantities include
mid-plane ion density and electron temperature, electron/ion power per area,
and electrode ion current.

For Case 1, the published chi-squared acceptance ranges are `55–303` at 95%
and `48–405` at 99%. Passing one stochastic realization is evidence of
statistical consistency, not a guarantee for arbitrary discharges; repeated
runs should leave those ranges at approximately their declared frequencies.

The prescribed collision set is Biagi 7.1 electron-He elastic momentum
transfer, two excitation channels, and ionization, all with isotropic
center-of-mass scattering. Ionization divides residual energy equally between
the primary and secondary electrons. He+-He scattering uses separate
isotropic and backward components. Intermediate cross sections are linearly
interpolated and values above the table range use the final tabulated value.

## AuroraPIC readiness

The first five bounded prerequisites are complete:

- 1D Dirichlet electrodes accept independent static offsets plus sinusoidal
  amplitude, frequency, and phase;
- the field solve applies the voltage at the new field time level;
- restart reconstructs the waveform phase from the stored simulation time;
- `scalars.csv` records the actual Dirichlet `phi_left` and `phi_right` plus
  live macro-particle counts for each species;
- `boundary_losses.csv` records restart-safe cumulative macro-particle count,
  represented charge, and impact kinetic energy by species and electrode;
- `examples/rf_electrode_1d.cfg` checks the zero, quarter-cycle, and
  half-cycle values with a bounded normalized 1D3V run;
- `velocity_dimensions = 3` retains transverse velocity through initialization,
  energy diagnostics, BGK, isotropic elastic/excitation MCC, and deterministic
  velocity-aware checkpoint/restart while preserving 1D1V as the default;
- `examples/mcc_relaxation.cfg` exercises the 1D3V MCC command-line path.
- named 1D collision models can target electrons and ions simultaneously with
  qualified diagnostics and a combined restart fingerprint;
- 1D3V ionization creates equal-weight electron/ion pairs, removes the
  configured threshold energy, preflights bounded storage, and defers newborn
  collisions until the next timestep;
- `examples/mcc_ionization_1d.cfg` exercises simultaneous electron/ion MCC,
  reactive products, diagnostics, capacity, and checkpoint integration.
- elastic tables can explicitly select projectile or center-of-mass lookup
  energy, and 1D3V supports separate finite-mass isotropic and exactly backward
  channels. For equal-mass He+-He, the backward channel exchanges projectile
  and neutral velocities while conserving momentum and energy. Both choices
  are restart-fingerprinted.

The configuration is:

```ini
boundary = dirichlet
phi_left = 0
phi_right = 0
phi_right_amplitude = 450
phi_right_frequency = 13.56e6
phi_right_phase = 0
```

In SI mode, voltage is in volts, frequency is in hertz, time is in seconds,
and phase is in radians. In normalized mode, frequency is cycles per
normalized time unit. A nonzero amplitude requires a positive frequency and
Dirichlet boundaries. Driven cases require `mode = transient`; the existing
instantaneous energy-window steady-state test is intentionally rejected until
cycle-averaged convergence is implemented.

## Blocking capabilities

AuroraPIC must not claim a Turner result until all of these are complete:

1. retain and audit the locally normalized, checksum-pinned collision tables
   without changing their interpolation contract, while keeping the
   publisher-derived files local unless redistribution permission is
   established;
2. finish volume power-transfer diagnostics. Species-resolved electrode
   charge/current and wall kinetic power are derivable from restart-safe
   cumulative boundary losses, ionization source is recorded by the collision
   diagnostics, and spatial-density time averaging is complete;
3. connect the implemented chi-squared comparator to a statistically bounded,
   run-contract-validated campaign;
4. implement whole-RF-cycle convergence and phase/time averaging;
5. run the full case only through an explicit production profile. Case 1
   alone requires 512,000 steps at the published resolution, so it must never
   become an ordinary laptop/CI regression.

Open-access implementations such as the
[WarpX capacitive-discharge example](https://warpx.readthedocs.io/en/latest/usage/examples/capacitive_discharge/README.html)
are useful independent integration references. LXCat is an open-access
platform, but access does not erase dataset-specific attribution, version,
or redistribution conditions. AuroraPIC's existing local import workflow
therefore remains the route for user-supplied data until an exact
redistributable benchmark package is identified and reviewed.

## Publisher supplement lock

The exact publisher supplement was manually acquired on 30 July 2026. The
committed `examples/turner_ccp.sources` registry pins the 133,969-byte archive
as SHA-256
`a0a5fe93900d7d7b213157f1eab664e06aab6e718f2189910e65f23bd699d661`
and separately pins all four members. It contains the Biagi 7.1 electron
table, the 101-row centre-of-mass He+-He table, original results for all four
published grids, and refined results. Case 1 has the required 129 nodes over
the 6.7 cm gap.

The electron table carries an embedded `All rights reserved` notice.
Consequently the archive and extracted tables remain under ignored `tmp/`;
AuroraPIC commits only identity and semantic metadata. Verify a local copy
without modifying or redistributing it with:

```sh
python3 scripts/verify_turner_source.py \
  examples/turner_ccp.sources \
  tmp/013507_1_supplements.zip \
  --output tmp/turner-source-verification.json
```

Normalize the verified bytes into local AuroraPIC manifests, channel tables,
reference CSV files, and a transformation audit with:

```sh
python3 scripts/normalize_turner_source.py \
  examples/turner_ccp.sources \
  tmp/013507_1_supplements.zip \
  --output-dir tmp/turner-normalized-v3
```

Electron thresholds are converted from eV to joules using the 2006 CODATA
elementary charge required by the paper. Normalization audit version 2 records
that constant edition explicitly. Ion energy remains tabulated in centre-of-mass eV and ion
cross sections remain in their source `1e-20 m2` scale; the gas manifest
performs both SI conversions at load time. No table is resampled. Named 1D
collision models load either generated manifest through `gas_data_file`, with
only ionization product-species mappings supplied by the simulation deck.

## Exact Case 1 campaign preflight

[`examples/turner_helium_ccp_case1.case`](../examples/turner_helium_ccp_case1.case)
locks the published Case 1 physics and numerical contract independently of
the restricted local tables. After normalization, generate the production
deck and its machine-readable preflight report with:

```sh
python3 scripts/prepare_turner_case.py \
  examples/turner_helium_ccp_case1.case \
  tmp/turner-normalized-v3 \
  --output tmp/turner-case1-campaign/turner_case1.cfg \
  --acknowledge-cost I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_RUN
```

The preparer verifies the normalization audit and every normalized file,
derives the exact 65,536 macro-particles per species and weight, locks the
512,000-step/1,280-cycle duration and final 12,800-sample window, and audits
thermal-neutral collision majorants through the declared 10 keV projectile
envelopes. The scan follows the production kernel's eight-standard-deviation
thermal bound, applies a safety factor, and records both sampled peaks and
configured majorants. During execution, the collision kernel independently
checks every encountered particle and fails rather than biasing the result if
a majorant is exceeded.

Generation never launches the calculation and the report explicitly carries
`physics_claim = none`. The exact campaign has a floor of 67,108,864,000
initial-particle updates before collision work and particle creation, so it
must first receive a bounded, single-core runtime qualification. Use
`aurorapic_cli --validate-only` to parse the generated deck without allocating
or advancing its particles.

Run that qualification explicitly with:

```sh
python3 scripts/qualify_turner_runtime.py \
  build/aurorapic_cli \
  examples/turner_helium_ccp_case1.case \
  tmp/turner-normalized-v3 \
  --steps 4 \
  --max-initial-updates 1000000 \
  --timeout-seconds 60 \
  --acknowledge-cost I_UNDERSTAND_THIS_IS_A_BOUNDED_TURNER_PROBE \
  --report tmp/turner-case1-runtime-qualification.json
```

The qualifier regenerates and audits the exact production deck, changes only
duration, output, checkpoint, averaging, and serial-runtime controls for the
probe, then requests lower process priority and one CPU. It will not exceed
its built-in two-million-initial-update ceiling, refuses to overwrite a
report, and never forwards the CLI's large-run authorization. Its projection
is an initial-population-only planning estimate, not a promise: ionization,
electrode losses, checkpoints, storage, and long-run hardware behavior can
all change the production cost.

The next bounded rung executes exactly one RF cycle in two half-cycle stages:

```sh
python3 scripts/run_turner_startup.py \
  build/aurorapic_cli \
  examples/turner_helium_ccp_case1.case \
  tmp/turner-normalized-v3 \
  --work-dir tmp/turner-case1-startup \
  --report tmp/turner-case1-startup/report.json \
  --max-initial-updates 60000000 \
  --timeout-seconds 120 \
  --acknowledge-cost I_UNDERSTAND_THIS_IS_A_ONE_CYCLE_TURNER_STARTUP
```

This stage advances 52,428,800 initial-particle updates, checkpoints at the
half-cycle, restarts to the full cycle, and retains both checkpoints and both
diagnostic directories under ignored `tmp/`. Its report requires exact
integer population and collision-counter continuity, floating diagnostic
continuity to roundoff, the prescribed electrode waveform, finite energy and
charge histories, charge-paired ionization/population balance, and a
129-node final field. It reports an early boundary-to-bulk field indicator
but explicitly does not interpret it as a stationary sheath.

Extend that retained one-cycle checkpoint through cycles two to four with:

```sh
python3 scripts/extend_turner_horizon.py \
  build/aurorapic_cli \
  examples/turner_helium_ccp_case1.case \
  tmp/turner-normalized-v3 \
  --prior-work-dir tmp/turner-case1-startup \
  --prior-report tmp/turner-case1-startup/report.json \
  --work-dir tmp/turner-case1-horizon \
  --report tmp/turner-case1-horizon/report.json \
  --additional-cycles 3 \
  --max-initial-updates 160000000 \
  --timeout-seconds 120 \
  --acknowledge-cost I_UNDERSTAND_THIS_IS_A_BOUNDED_TURNER_HORIZON
```

The horizon runner verifies the prior report and checkpoint hashes, advances
one whole cycle per low-priority serial process, and writes a distinct
checkpoint and phase-matched diagnostic set after every cycle. A completed
horizon report/work directory can be supplied as the next `--prior-*` pair,
so bounded blocks chain without repeating earlier cycles. Per-cycle
reports contain collision deltas, charge-paired population balances, inferred
electrode losses, energy and charge ranges, waveform error, restart
continuity, and boundary-to-bulk field structure. Four added cycles are the
built-in maximum; this is a startup trend screen, not a stationarity test.

Screen a contiguous hash-chained report history before authorizing any
published-profile comparison:

```sh
python3 scripts/analyze_turner_stationarity.py \
  tmp/turner-case1-startup/report.json \
  tmp/turner-case1-horizon/report.json \
  tmp/turner-case1-horizon-next/report.json \
  --window-cycles 4 \
  --max-population-change 0.005 \
  --max-observable-span 0.05 \
  --output tmp/turner-case1-stationarity.json
```

The default engineering screen requires both species to change by no more
than 0.5% per cycle throughout the final four cycles and requires the
ionization count, phase-zero boundary field, and phase-zero total energy each
to span no more than 5% across that window. These are conservative AuroraPIC
startup gates, not thresholds published by Turner et al. Passing only permits
a longer confirmation window; failing explicitly keeps the published
chi-squared comparison inapplicable.

The case manifest separately pins the published Case 1 summary from Table III:
mid-plane ion density `1.40e14 m^-3`, electron temperature `9.36 eV`,
electron and ion power `34.3` and `90.6 W m^-2`, ion current
`0.219 A m^-2`, and approximately 31,900 total macro-particles. The
stationarity report includes the current-to-reported macro-particle ratio as
context, but does not turn that reported stochastic population into an
acceptance gate.

## Restart-safe density averaging

The primary Turner observable is now produced by a generic 1D spatial
averager. By default it samples after the collision operator. The optional
`spatial_average_sampling_order = pre_collision` mode samples the same
completed-step state immediately before collisions, without changing particle
evolution; this exists for matched external-code diagnostic protocols. It
deposits each species' represented particle number with
the same linear particle-grid shape and node-volume convention used for
charge, then accumulates nodal density independently of the ordinary output
interval. For Case 1, the exact final-32-cycle contract is:

```ini
spatial_average = true
spatial_average_interval = 1
spatial_average_start_step = 499201
spatial_average_end_step = 512000
spatial_average_rf_frequency = 13.56e6
spatial_average_rf_cycles = 32
spatial_average_phase_bins = 16
```

This selects 12,800 post-collision samples. Configuration validation requires an
integer number of timesteps per RF cycle, a whole-cycle window ending at the
time-zero drive phase, an interval that divides the cycle, and agreement with
every active electrode-drive frequency. `spatial_average.csv` uses long-form
species/node density rows. `spatial_kinetic_energy.csv` reports the
density-weighted mean total kinetic energy and the explicitly labeled
effective kinetic temperature `2 <K> / velocity_dimensions`; it includes
directed energy and is not silently presented as a thermodynamic temperature.
`spatial_field_average.csv` reports nodal mean potential, mean electric field,
and RMS electric field. These are sheath-localization observables, not an
automatic sheath-edge definition. `spatial_average_metadata.json` records the
window, timestep, density and moment sample counts, definitions, species, and
completeness gates.

Optional `spatial_average_phase_bins` must divide the number of sampled steps
per RF cycle. It writes `spatial_phase_moments.csv` with species-resolved
density, mean three-velocity, mean kinetic energy, and drift-separated
temperature in each phase bin, plus `spatial_phase_fields.csv` with matching
potential and electric-field statistics. The drift-separated definition is
`2/d * (<K> - m|<v>|^2/2)`, where `d` is the configured velocity dimension.
Metadata explicitly records `sampling_order`. Checkpoint v17 preserves that
contract in addition to every phase-bin accumulator and sample count; bounded
regressions require byte-identical continuous and checkpoint-split products.

1D checkpoint v5 and later stores the density averaging contract, sample
count, and every nodal density sum; v6 adds species/side wall count and impact
energy, v7 adds species electric work, v8 adds kinetic-energy, potential,
mean-field, and squared-field sums, and v9 adds phase-binned velocity and field
moments. A changed averaging window is rejected by default.
When a pre-v8 checkpoint continues an active density window, density remains
restart-correct but the new moment products remain header-only and
`moments_complete` is false; use an explicit post-checkpoint reset window to
obtain complete moment profiles.
Setting `spatial_average_reset_on_restart = true` explicitly discards stored
sums and permits a new averaging contract; its first sample must be after the
checkpoint step. Metadata records the reset, and the Turner comparator exposes
a separate post-benchmark mode that never applies the published-duration
acceptance gate. A legacy v1-v4 checkpoint can otherwise restart only when
spatial averaging is disabled. Older restarts retain the state they support
and start each newer, explicitly origin-labeled counter at the restart step.
Bounded regressions prove byte-identical continuous/checkpoint-split and
restart-reset profiles, represented-number conservation, exact wall
accounting, power-work closure, and counter restart continuity.

The baseline statistical comparison uses the original-grid reference, not the
numerically refined profile. For every mesh node it computes

```text
X^2 = sum((candidate_ion_mean - reference_ion_mean)^2
          / reference_ion_population_stddev^2).
```

The population standard deviation is the denominator; the standard deviation
of the reference mean is not. Once a candidate final-32-period profile exists,
run:

```sh
python3 scripts/compare_turner.py \
  --case 1 \
  --reference tmp/turner-normalized-v3/turner_case1_benchmark.csv \
  --candidate tmp/turner-case1-run/spatial_average.csv \
  --candidate-metadata \
    tmp/turner-case1-run/spatial_average_metadata.json \
  --species ions \
  --normalization-audit tmp/turner-normalized-v3/audit.json \
  --output tmp/turner-case1-comparison.json
```

The tool requires ordered correspondence to the prescribed uniform mesh and
does not interpolate. Candidate coordinates must match that mesh to numerical
roundoff. The publisher supplement prints coordinates with limited precision,
so reference coordinates receive a narrowly bounded `2.5e-4`-cell rounding
allowance; the report records both maximum coordinate errors. The tool rejects
partial data unless metadata proves the exact case timestep and final 32 RF
cycles. A report still makes no overall physics claim because density metadata
alone cannot prove the prescribed particle count, collision model, or prior
stationary state.

## First exact Case 1 campaign result

The first exact Case 1 campaign completed all 512,000 steps on 30 July 2026
using the 2006-CODATA-normalized local supplement. The final averaging metadata
records all 12,800 required samples over steps 499,201--512,000 and the final
checkpoint is present. The retained local artifact identities are:

| Artifact | SHA-256 |
| --- | --- |
| Generated deck | `b39c9e70b92757f561c88f7d1a2a03bd0c617861a393d6fce12e53c445d13927` |
| Normalization audit | `29ed7dae41e329d8187e51a13b850f999d13e88ea291c56a2d633d253deca740` |
| Final checkpoint | `ba8383387194671ab77001bd8892027a170d67ed7e46a5c23ed7b6ab786c2092` |
| Averaged profile | `13f7dbfdd12f81986c28251676a895b8b67a2fcd664d7e265fa8e44cb750a702` |
| Averaging metadata | `3a8ad14babf4fc23a510e788c7a08c5cf6f6e7276241505087f1c47fda064278` |
| Comparison report | `5f9f5fb7a50c7e438552e74276f5161438698cde785f2bb1702fb22531f9bc4f` |

The unmodified published comparison gives `X² = 574.399`, outside both the
Case 1 95% range (`55--303`) and 99% range (`48--405`). The relative profile
L2 error is 2.874%, the maximum pointwise relative error is 4.811%, and the
mid-plane ion density is 2.211% above the reference mean. The final
macro-particle population is 32,280 and the mean of the final 32
cycle-boundary samples is 32,157.6, close to the approximately 31,900 reported
by Turner et al.

A post-hoc weighted amplitude fit multiplies the candidate density by
`0.980363` and reduces `X²` to `121.834`; an unweighted amplitude removal
leaves 1.036% relative shape L2 error. These fitted values are diagnostic and
are not part of the published test, so they do not change the failed campaign
classification. The candidate symmetry error is 1.429% in relative L2, and
the live population rises 1.465% across the final 32 cycle-boundary samples.
The result proves the complete benchmark execution and comparison path, but
does **not** yet establish Turner code-to-code verification. The next
credibility work must isolate late-time stochastic drift and the roughly 2%
density-amplitude bias before another production campaign is justified.

## Post-benchmark source/wall diagnostic

Checkpoint v6 wall accounting was introduced after the exact 1,280-cycle run,
so a diagnostic continuation started from its final v5 checkpoint. One cycle
first established the counter transition; a subsequent 32-cycle window covered
steps 512,400--525,200. These steps are deliberately outside the published
duration and cannot change the failed density comparison. The retained local
artifact identities are:

| Artifact | SHA-256 |
| --- | --- |
| Solver executable (commit `d11875a`) | `c05bf308f18a21a32ca23e520b1410ac7267d329a7f23025395f9f64863f015c` |
| Original step-512,000 v5 checkpoint | `ba8383387194671ab77001bd8892027a170d67ed7e46a5c23ed7b6ab786c2092` |
| One-cycle diagnostic deck | `230a17ed38f02b269780692d09051eb51229b7bdae21180b1608e75841aeab03` |
| 32-cycle diagnostic deck | `5c91148fc5a5a2fb6b59f5f832050eea2da407a18370a31d3ba19cc7902334c5` |
| Final step-525,200 v6 checkpoint | `23c82ef6ce77a8acf9b4465d9c38aab0da38621da0dfd89a1c9e3fc54061cd95` |
| 32-cycle balance report | `e224564a6d9aa06a517a86090703920dce1e7e449b09408e68c3ff6ce7a813e3` |

The 32-cycle window recorded 24,657 ionizations. Electron absorption was
12,721 left plus 12,038 right, giving the exact observed population change of
-102. Ion absorption was 12,325 left plus 12,478 right, giving the exact
observed change of -146. Both integer balance residuals are zero.

The left and right ion-current magnitudes are `0.218999` and
`0.221717 A m^-2`; their mean is `0.220358 A m^-2`, 0.620% above the Turner
Table III value `0.219 A m^-2`. This has no published single-window acceptance
range, but it is an encouraging independent observable. The corresponding
mean electron absorbed-current magnitude is `0.219967 A m^-2`, within 0.18%
of the ion mean. The reported wall kinetic powers are particle energy carried
to each electrode, not the volume electron/ion power transfer in Table III.

The total live population fell from 32,298 to 32,050 during this continuation,
moving toward the approximately 31,900 reported value. Together with exact
source/loss closure and close electrode current, this argues against a gross
wall-loss implementation error. It does not resolve whether the original
density excess is a seed fluctuation, residual finite-duration drift, or a
smaller collision/heating discrepancy.

### Species electrical-power continuation

For 1D collision-enabled runs, `collisions.csv` also reports a signed tracked
particle kinetic-energy change for every collision channel, both per output
interval and cumulatively. Positive values add kinetic energy to simulated
charged particles; negative values remove it. SI columns use `J_m-2` because
a 1D run represents energy per unit transverse area; normalized columns are
explicitly labeled. Ionization includes the post-collision primary and the
created electron and ion, so its entry is the net tracked charged-particle
kinetic-energy change. Neutral kinetic energy, internal excitation energy, and
ionization potential are outside that tracked-particle ledger and appear as a
signed channel loss. Checkpoint v11 preserves the cumulative channel totals
and the spatial/phase collision-energy state.

When spatial averaging is enabled, `spatial_collision_power.csv` deposits each
signed channel transfer onto the same nodal control volumes and with the same
linear shape function used for particle moments. Its mean power density uses
every physical timestep in the declared window, independently of
`spatial_average_interval`. Consequently, integrating a channel over nodal
volume exactly recovers its global collision ledger over the same window.
`spatial_phase_collision_power.csv` performs the same accounting in the
configured RF bins and records the physical timestep count and duration of
each bin. Summing phase-bin energy densities recovers the unbinned spatial
field. SI energy and power densities are `J m^-3` and `W m^-3`; normalized
runs label both quantities explicitly. The diagnostic is empty when no
collision channels are configured; when spatial averaging is disabled it
allocates no deposition state and the event hook returns immediately.

A one-cycle Case 1 smoke window at steps 757,201--757,600 exercised all six
helium collision channels with 16 RF bins of exactly 25 timesteps each. The
largest difference between a mesh-integrated spatial channel and its global
ledger was `3.98e-20 J m^-2`; the largest difference between summed phase bins
and the unbinned spatial channel was `2.54e-21 J m^-2`. This verifies the real
Case 1 deposition, phase partition, pre-v11 restart-reset path, and v11 output
checkpoint. It remains an implementation smoke test rather than a statistically
converged physics comparison. Evidence and artifact hashes are recorded in
[`turner-case1-spatial-collision-smoke-20260804.json`](../benchmarks/ccp/turner-case1-spatial-collision-smoke-20260804.json).

The first bounded application was a four-RF-cycle pilot continuing Case 1
from step 755,600 through 757,200. The independent tracked-particle balance

```text
dK/dt = electric work + collision kinetic transfer - wall kinetic loss
```

closed with a residual of `1.02e-11 W m^-2` (`8.06e-14` relative). Electric
work was `127.061 W m^-2`, collisions removed `81.689 W m^-2`, walls removed
`46.290 W m^-2`, and the measured kinetic-energy rate was
`-0.917814 W m^-2`. Every configured electron and ion channel was exercised
and had the physically expected negative tracked-particle kinetic transfer.
This is a strong implementation-level conservation result, but the four-cycle
channel powers are not a new Turner acceptance comparison; production means
still require a predeclared longer window. The reusable calculation is
`scripts/analyze_1d_energy_budget.py`, and the checksum-pinned evidence is
[`turner-case1-energy-ledger-pilot-20260804.json`](../benchmarks/ccp/turner-case1-energy-ledger-pilot-20260804.json).

The predeclared 32-cycle localization window then covered steps
757,601--770,400 with 16 bins of 800 physical timesteps. Its global kinetic
ledger closed to `1.18e-11 W m^-2` (`9.41e-14` relative); maximum
spatial-to-global and phase-to-spatial channel residuals were respectively
`1.25e-19` and `2.98e-19 J m^-2`. Integer electron and ion source/loss
balances were exact.

The post-benchmark comparison is especially encouraging but remains outside
the formal published-duration gate. Ion-density `X^2 = 218.4` lies inside
Turner's published 95% range, electron and ion electrical powers differ from
Table III by `-0.17%` and `+0.52%`, and mean ion current differs by `+0.44%`.
Midplane electron temperature remains `6.50%` low while midplane ion density
is `3.43%` high. Thus global power deficit is not a plausible explanation for
the temperature discrepancy.

Electron inelastic losses were `4.713 W m^-2` (triplet excitation),
`8.404 W m^-2` (singlet excitation), and `10.670 W m^-2` (ionization).
Approximately 98% of each electron-channel loss occurred outside the two
predeclared outer-10% geometric regions. Ion-neutral losses were more edge
weighted: the outer regions contained 50.9% of backward-channel and 41.8% of
isotropic-channel transfer. Ionization and singlet loss were strongly RF
modulated (60.9% and 49.3% conditional peak-to-trough over mean) while their
half-cycle asymmetry remained only 2.55% and 3.35%. These observations point
next to phase-selected bulk/sheath EEDFs and threshold-resolved moments, not
another global power diagnostic. The reusable analyzer is
`scripts/analyze_1d_spatial_collision.py`; the full claim boundaries and
checksums are in
[`turner-case1-spatial-collision-32cycle-20260804.json`](../benchmarks/ccp/turner-case1-spatial-collision-32cycle-20260804.json).

Checkpoint v12 adds configurable phase-resolved regional EEDFs. The first
32-cycle campaign used 0.5 eV bins through 500 eV in the central 20% and both
outer 20% regions. Overflow and histogram-normalization residual were zero,
and the independent kinetic ledger remained closed to `2.10e-11 W m^-2`.
The bulk mean energy was `13.562 eV`; its drift-separated temperature ranged
from `8.275` to `9.458 eV`, and its half-cycle histogram total variation was
only `0.79%`. Fractions above triplet, singlet, and ionization thresholds were
22.38%, 21.50%, and 15.57%. The populated outer-region phases were hotter,
with approximately 22.4--22.5% above ionization threshold, while complementary
phase bins were completely electron depleted near opposite electrodes. These
are well-closed AuroraPIC diagnostics, not a published EEDF validation. The
next credible discriminator is an independently generated EEDF with identical
regions and phase bins. Evidence is recorded in
[`turner-case1-phase-eedf-32cycle-20260804.json`](../benchmarks/ccp/turner-case1-phase-eedf-32cycle-20260804.json).

### Independent cross-code discriminator

The independent comparison now has a code-neutral interchange rather than an
AuroraPIC-specific CSV contract. `export_phase_eedf.py` converts native
phase/region histograms to probability mass on explicit energy-bin edges;
`compare_phase_eedf.py` integrates total variation exactly for the resulting
piecewise-uniform distributions even when the two solvers use different
energy grids. Exact mean energy, drift-separated temperature, overflow, and
optional energetic-tail fractions are reported separately. Formal pass/fail
is disabled unless all acceptance limits were declared on the command line.
The format and claim boundary are documented in
[`phase-eedf-interchange.md`](phase-eedf-interchange.md).

The first external target is eduPIC 1.0, pinned for assessment at upstream
commit `32050728c961a317d6d6acd6bc86d026da403326`. It is an independently
implemented, GPL-3.0, transparent 1D3V electrostatic PIC/MCC code whose
published default argon CCP is intended to reproduce the paper's reference
figures. Its default contract is 25 mm, 10 Pa argon at 350 K, 250 V amplitude,
13.56 MHz, 400 nodes, and 4,000 electron steps per RF cycle. It records the
central 10% EEPF, cycle-averaged density, 200-bin phase-space diagnostics, wall
flux/ion energy, and species `j E` power.

This is intentionally a second benchmark, not a replacement for Turner helium
Case 1. The unmodified eduPIC reference is argon, uses analytic Phelps-based
cross sections and its own collision/scattering conventions, and normally
requires about 1,500 equilibration plus at least 1,000 measurement cycles.
Running it immediately would be both physically unmatched to Turner and too
expensive for the workstation safety policy. The next campaign gate is
therefore: reproduce eduPIC's unmodified reference in conservative,
checkpointed blocks; lock its actual stationary window and output hashes;
then run AuroraPIC under that same argon/collision/boundary contract and feed
both outputs through the neutral comparator. No cross-code agreement claim is
made yet.

A one-cycle initialization-only pilot of that exact upstream commit completed
on one affinity-pinned, low-priority CPU in 0.56 s with 31.6 MiB peak resident
memory and no swap. The population changed from 1,000 electrons plus 1,000
ions to 835 electrons plus 3,318 ions. This proves only that the upstream code
builds and that its checkpoint path works locally. Because the published
reference approaches roughly 110,000 particles, the startup timing must not
be extrapolated linearly to the stationary campaign. The next run remains
gated behind measured checkpoint blocks. The checksum-bearing runtime record
is
[`edupic-reference-init-pilot-20260804.json`](../benchmarks/ccp/edupic-reference-init-pilot-20260804.json).

`run_edupic_stage.py` now enforces the external equilibration ladder. Each
stage copies rather than overwrites its input state, validates the binary
checkpoint layout against its exact byte count, cross-checks cycle and species
counts with `conv.dat`, limits requested cycles and initial particle-timesteps,
pins one CPU, applies low priority and a hard timeout, captures output, and
writes a hash-bearing report only after exact cycle coverage is proved. The
synthetic regression also verifies that the input checkpoint is unchanged.

The first guarded continuation advanced only cycles 2--3. It completed in
0.85 s under a 10 s timeout and a 40 million initial-particle-timestep limit;
the total population increased from 4,153 to 6,788. This remains rapidly
evolving initialization, not stationarity. Its purpose is to prove the safe
stage mechanism before increasing any block size. Evidence is
[`edupic-reference-equilibration-stage-0001-0003-20260804.json`](../benchmarks/ccp/edupic-reference-equilibration-stage-0001-0003-20260804.json).

The runner was then strengthened to require predeclared executable and input
checkpoint SHA-256 values. A one-cycle exercise of those gates advanced cycle
3 to 4 in 0.70 s; population rose again, from 6,788 to 8,109. The preserved
cycle-4 checkpoint is the only authorized input to the next stage. The result
also demonstrates why block size must grow from measured state rather than
startup timing: throughput already fell from 2.36 to 1.42 cycles/s as the
population increased. Evidence is
[`edupic-reference-equilibration-stage-0003-0004-20260804.json`](../benchmarks/ccp/edupic-reference-equilibration-stage-0003-0004-20260804.json).

Two further two-cycle stages reached cycle 8 under respective 70 and 90
million initial-particle-timestep ceilings. Population rose from 8,109 to
10,685 and then 13,129; throughput decreased to 1.13 cycles/s. The complete
eight-cycle trace still has a total-population relative linear slope of
`+0.148` per cycle and an endpoint change of `+1.032` relative to its mean.

`analyze_edupic_convergence.py` now makes the equilibration gate objective and
predeclared. It cannot become eligible before cycle 1,500 and requires the
last 100 cycles split into four 25-cycle blocks. Electron, ion, and total
population must each have absolute relative linear slope no greater than
`1e-4` per cycle and block-mean range no greater than 2% of the window mean.
Passing this gate authorizes a measurement window; it is not itself a physics
validation. Current status is correctly ineligible and nonstationary. Evidence
is
[`edupic-reference-equilibration-through-cycle8-20260804.json`](../benchmarks/ccp/edupic-reference-equilibration-through-cycle8-20260804.json).

`advance_edupic_equilibration.py` now chains these immutable stages without
manual checkpoint handling. It verifies the initial binary and checkpoint
hashes, sizes every stage from the current population and a particle-timestep
budget, limits individual and aggregate wall time, delegates to the one-core
runner, and atomically updates a recovery manifest after every completed
stage. A failed or interrupted stage cannot overwrite its input or erase the
last completed checkpoint; a new bounded campaign can restart from that stage.

Longer campaigns may additionally declare `--max-host-load-per-cpu`,
`--min-available-memory-mib`, and `--max-swap-io-pages-per-stage`. The
coordinator samples only the enabled Linux host metrics before and after every
immutable stage, records each sample and any violation in the campaign
manifest, and stops before launching more work when a threshold is crossed.
These guards complement rather than replace CPU affinity, reduced priority,
particle-work limits, stage timeouts, and the aggregate wall-time ceiling.

The first adaptive campaign used four two-cycle stages to advance cycle 8 to
16 in 10.75 s under a 20 s aggregate ceiling. Population rose from 13,129 to
22,618 while per-stage wall time rose from 2.12 to 3.13 s. The relative total
population slope over all 16 cycles decreased from the cycle-8 diagnostic but
remains strongly positive at `+0.0902` per cycle. It is therefore still early
equilibration and ineligible for measurement. Evidence is
[`edupic-reference-adaptive-through-cycle16-20260804.json`](../benchmarks/ccp/edupic-reference-adaptive-through-cycle16-20260804.json).

The next adaptive advance exposed and safely contained a scheduling edge case.
Seven stages reached cycle 31, after which the residual aggregate wall budget
reduced the final stage timeout to 3 s, slightly below its required runtime.
The runner timed out without changing the cycle-31 input. The coordinator now
predicts the next stage from the preceding measured per-cycle time with a 1.5
safety factor and stops cleanly before launch unless both aggregate wall time
and the stage timeout cover that prediction. A regression locks this behavior.

A new hash-locked recovery campaign then advanced the preserved cycle-31 state
to cycle 32 in 3.08 s. Population is 40,431 and the full-history total
population slope has fallen to `+0.0509` per cycle, but this is still rapid
growth and remains ineligible for measurement. The failed attempt, scheduler
correction, recovery chain, and final hashes are retained in
[`edupic-reference-adaptive-through-cycle32-20260804.json`](../benchmarks/ccp/edupic-reference-adaptive-through-cycle32-20260804.json).

The following campaign advanced cycle 32 to 44 in seven stages totaling 43.62
s of solver wall time. Population increased from 40,431 to 52,332; the
full-history total-population slope decreased to `+0.0380` per cycle but is
still far from stationary. Stage sizing automatically changed from two cycles
to one as the particle-work ceiling became active.

Observation of the final stage occurred after its immutable report and
checkpoint existed but before the campaign manifest was known to contain it.
The coordinator now supports `--resume-existing`: it replays and hashes every
recorded stage from the original locked checkpoint, then reconciles exactly
one contiguous completed-but-unrecorded stage directory. Changed, ambiguous,
or discontinuous chains are rejected. A synthetic interruption regression
proves incorporation without rerunning the completed stage. The manifest also
separates accumulated solver-stage wall time from individual coordinator
invocation overhead. Evidence is
[`edupic-reference-adaptive-through-cycle44-20260804.json`](../benchmarks/ccp/edupic-reference-adaptive-through-cycle44-20260804.json).

Five further two-cycle stages advanced cycle 44 to 54 in 41.48 s of
coordinator wall time. Population increased from 52,332 to 61,136, while the
total relative slope decreased to `+0.0311` per cycle. Fifty-four samples now
contain two complete 25-cycle blocks, whose total-population means are 18,882
and 46,322. Their range is 79.3% of the full-history mean, far above the
predeclared 2% limit. Electron and ion relative slopes are independently
positive at `+0.0341` and `+0.0288` per cycle. This quantitatively confirms
that the default eduPIC discharge is still forming, not merely fluctuating
around a stationary population. Evidence is
[`edupic-reference-adaptive-through-cycle54-20260804.json`](../benchmarks/ccp/edupic-reference-adaptive-through-cycle54-20260804.json).

Another five two-cycle stages reached cycle 64 in 48.56 s. Population rose to
69,145. The full-history total slope is `+0.0262` per cycle and the two
complete 25-cycle block means still have a 69.5% relative range.

The convergence analyzer now also reports a 25-cycle provisional recent
window, explicitly labeled descriptive and excluded from the stationarity
gate. Over cycles 40--64, electron, ion, and total-population relative slopes
are respectively `+0.0155`, `+0.0136`, and `+0.0145` per cycle; total
population changes by 34.9% relative to that window's mean. This is a better
measure of current equilibration speed than the entire startup history, while
leaving the cycle-1,500 and final-100-cycle contract unchanged. Both views say
the discharge is still substantially transient. Evidence is
[`edupic-reference-adaptive-through-cycle64-20260804.json`](../benchmarks/ccp/edupic-reference-adaptive-through-cycle64-20260804.json).

Five more two-cycle stages reached cycle 74 in 53.46 s. Population is 76,476.
The full-history total relative slope is `+0.0225` per cycle and the first two
25-cycle block means retain a 62.3% range relative to the full-history mean.
The descriptive recent window, now cycles 50--74, has mean population 67,413,
relative slope `+0.0116` per cycle, and an endpoint increase equal to 27.8% of
its mean. The recent rate is declining, but it remains more than two orders of
magnitude above the formal `1e-4`-per-cycle limit. Evidence is
[`edupic-reference-adaptive-through-cycle74-20260804.json`](../benchmarks/ccp/edupic-reference-adaptive-through-cycle74-20260804.json).

The next invocation advanced cycle 74 to 82 in four two-cycle stages, then
stopped cleanly because the safety-adjusted next-stage prediction exceeded the
remaining stage timeout. A separate hash-locked two-cycle invocation completed
cycle 84. Population is 83,006. The first three complete 25-cycle block means
are 18,882, 46,322, and 68,191, showing the long startup trajectory directly.
The descriptive cycles 60--84 window has mean 74,849, relative slope
`+0.00947` per cycle, and endpoint change 22.7% of its mean. Current growth is
now below 1% per cycle but remains approximately 95 times the formal slope
limit. Evidence is
[`edupic-reference-adaptive-through-cycle84-20260804.json`](../benchmarks/ccp/edupic-reference-adaptive-through-cycle84-20260804.json).

Five further two-cycle stages reached cycle 94 in 66.77 s. Population is
89,372. The descriptive cycles 70--94 window has mean population 81,671;
electron, ion, and total relative slopes are `+0.00839`, `+0.00759`, and
`+0.00797` per cycle. Its endpoint increase is 19.2% of the recent mean.
Equilibration is consistently decelerating, but the total recent slope remains
about 80 times the formal limit and no measurement is authorized. Evidence is
[`edupic-reference-adaptive-through-cycle94-20260804.json`](../benchmarks/ccp/edupic-reference-adaptive-through-cycle94-20260804.json).

The next five two-cycle stages reached cycle 104 in 67.81 s. Population is
95,355. This is the first checkpoint with four complete 25-cycle population
blocks: their total-population means are 18,882, 46,322, 68,191, and 85,535.
The block range is 118% of the 104-cycle mean, and the full-history relative
slope is `+0.0155` per cycle. The descriptive cycles 80--104 window has mean
population 88,040; electron, ion, and total relative slopes are `+0.00743`,
`+0.00677`, and `+0.00708` per cycle, with an endpoint increase equal to 16.9%
of its mean. The decreasing recent slope is encouraging for equilibration,
but it remains about 71 times the formal limit. Cycle 104 is also far below
the predeclared cycle-1,500 eligibility boundary, so the four-block result is
diagnostic only and neither measurement nor cross-code comparison is
authorized. Evidence is
[`edupic-reference-adaptive-through-cycle104-20260804.json`](../benchmarks/ccp/edupic-reference-adaptive-through-cycle104-20260804.json).

Five more two-cycle stages reached cycle 114 in 71.75 s. Population is
100,799. The full-history total relative slope decreased to `+0.0140` per
cycle. Over the descriptive cycles 90--114 window, mean population is 94,015;
electron, ion, and total relative slopes are `+0.00644`, `+0.00588`, and
`+0.00614` per cycle. Its endpoint increase is 14.9% of the recent mean. The
recent rate continues to decline, but remains about 61 times the formal limit.
The run is therefore still equilibration rather than an authorized measurement
window. Evidence is
[`edupic-reference-adaptive-through-cycle114-20260805.json`](../benchmarks/ccp/edupic-reference-adaptive-through-cycle114-20260805.json).

The next bounded advance reached cycle 124 in five two-cycle stages and 87.00
s. Population is 105,966. The full-history total relative slope is `+0.0128`
per cycle. The descriptive cycles 100--124 window has mean population 99,619;
electron, ion, and total relative slopes are `+0.00563`, `+0.00518`, and
`+0.00539` per cycle, and its endpoint increase is 12.9% of the recent mean.
This is the third consecutive decrease in recent-window growth, but the total
rate remains about 54 times the formal limit. Equilibration and the explicit
no-measurement boundary therefore remain in force. Evidence is
[`edupic-reference-adaptive-through-cycle124-20260805.json`](../benchmarks/ccp/edupic-reference-adaptive-through-cycle124-20260805.json).

The first host-guarded campaign then reached cycle 134 in five two-cycle
stages and 82.85 s. All ten pre/post-stage checks passed: normalized one-minute
load peaked at 0.289 against a 0.5 limit, available memory bottomed at 4,879
MiB against a 4,096 MiB floor, and swap I/O remained zero. Population is
110,577. The new fifth complete 25-cycle block has mean population 100,155,
still well above the preceding block's 85,535. The descriptive cycles 110--134
window has mean 104,794 and total relative slope `+0.00478` per cycle, about
48 times the formal limit. The campaign establishes the operational guard but
does not change the no-measurement physics boundary. Evidence is
[`edupic-reference-guarded-through-cycle134-20260805.json`](../benchmarks/ccp/edupic-reference-guarded-through-cycle134-20260805.json).

A larger guarded window advanced cycle 134 to 154 in ten two-cycle stages and
171.57 s. All 20 host checks passed: normalized load stayed below 0.273,
available memory stayed above 6,074 MiB, and swap I/O stayed zero. Population
is 119,430. The sixth complete 25-cycle block has mean population 112,395,
compared with 100,155 in the fifth. The descriptive cycles 130--154 window has
mean population 114,195 and total relative slope `+0.00390` per cycle, about
39 times the formal limit. Safe campaign scaling is demonstrated, but the
continued block-to-block growth retains the no-measurement boundary. Evidence
is
[`edupic-reference-guarded-through-cycle154-20260805.json`](../benchmarks/ccp/edupic-reference-guarded-through-cycle154-20260805.json).

The next guarded campaign reached cycle 174 in 186.88 s. The particle-work
guard automatically reduced stage size from two cycles to one after cycle 168
as population crossed 125,000; all 13 stages remained below one billion
initial particle-steps. All 26 host checks passed, with normalized load below
0.170, available memory above 5,395 MiB, and zero swap I/O. Population is
127,760. The descriptive cycles 150--174 window has mean population 122,895
and total relative slope `+0.00347` per cycle, about 35 times the formal limit.
Adaptive stage reduction worked as intended, while the scientific state
remains equilibration rather than measurement. Evidence is
[`edupic-reference-guarded-through-cycle174-20260805.json`](../benchmarks/ccp/edupic-reference-guarded-through-cycle174-20260805.json).

The long-horizon guarded campaign subsequently reached cycle 250 through 76
immutable one-cycle stages and three coordinator invocations. Two invocations
stopped after valid stages when available memory crossed the 4,096 MiB floor;
neither launched more work, and the checkpoint chain contains no failed or
partial stage. Population is 151,437. The ninth and tenth complete 25-cycle
block means are 140,671 and 148,071. The descriptive cycles 226--250 window
has total relative slope `+0.00191` per cycle, about 19 times the formal
limit. This materially advances equilibration and exercises recovery under
real host pressure, but remains ineligible for measurement. Evidence is
[`edupic-reference-guarded-through-cycle250-20260805.json`](../benchmarks/ccp/edupic-reference-guarded-through-cycle250-20260805.json).

Continued guarded execution reached cycle 528 through 354 immutable one-cycle
stages and 14 coordinator invocations. Eight memory and two load violations
were detected only at stage boundaries and prevented subsequent launches; no
failed or partial stage entered the chain. Population is 195,803. The last
four complete 25-cycle block means are 186,453, 189,555, 192,252, and 194,633.
The descriptive cycles 504--528 window has total relative slope
`+0.000403` per cycle and endpoint change 0.952% of its mean. This is a major
reduction from cycle 250 and only about four times the formal slope limit, but
cycle 528 remains below the cycle-1,500 eligibility boundary and is not an
authorized measurement window. Evidence is
[`edupic-reference-guarded-through-cycle528-20260805.json`](../benchmarks/ccp/edupic-reference-guarded-through-cycle528-20260805.json).

The guarded campaign then reached cycle 763 through 589 immutable one-cycle
stages and 20 coordinator invocations. Population is 210,644. The last four
complete 25-cycle block means are 205,754, 207,234, 208,480, and 209,622. The
descriptive cycles 739--763 window has total relative slope `+0.000197` per
cycle and endpoint change 0.460% of its mean. This is about half the cycle-528
rate and roughly twice the formal slope limit. The trend is approaching
stationarity, but it is still above threshold and cycle 763 remains below the
cycle-1,500 eligibility boundary. Evidence is
[`edupic-reference-guarded-through-cycle763-20260805.json`](../benchmarks/ccp/edupic-reference-guarded-through-cycle763-20260805.json).

The campaign next reached cycle 1,022 through 848 immutable one-cycle stages
and 27 coordinator invocations. Population is 216,805. The descriptive cycles
998--1,022 window has electron, ion, and total relative slopes of
`+8.32e-5`, `+7.45e-5`, and `+7.88e-5` per cycle. This is the first checkpoint
where all three descriptive slopes lie below the formal `1e-4` limit, and its
total endpoint change is 0.261% of the recent mean. It is encouraging evidence
that equilibration is approaching completion, but the diagnostic is only 25
cycles, cycle 1,022 remains below the cycle-1,500 eligibility boundary, and
the required final 100 cycles have not yet been split into and tested as four
25-cycle blocks. No measurement is authorized. Evidence is
[`edupic-reference-guarded-through-cycle1022-20260805.json`](../benchmarks/ccp/edupic-reference-guarded-through-cycle1022-20260805.json).

The guarded reference campaign has now completed the predeclared cycle-1,500
stationarity milestone. From the locked cycle-174 input, 1,326 immutable
one-cycle stages completed across 41 bounded coordinator invocations. The
accepted chain contains no failed or recovered-unrecorded stage. Its 2,653
host-health checks recorded zero swap I/O. Eight low-memory and two high-load
checks stopped further work only at stage boundaries; no partial result
entered the chain.

The formal cycles 1,401--1,500 window is eligible and passes all six locked
criteria. Electron, ion, and total-population relative slopes are
`1.12e-5`, `9.52e-6`, and `1.03e-5` per cycle, respectively, against the
`1e-4` absolute limit. Their four 25-cycle block-mean ranges are `0.0851%`,
`0.0714%`, and `0.0781%`, respectively, against the `2%` limit. This is the
first formal authorization to begin the upstream eduPIC measurement phase.
It establishes population stationarity only: it does not yet validate
eduPIC observables, AuroraPIC physics, or cross-code agreement. Checksum-bearing
evidence is
[`edupic-reference-stationarity-cycle1500-20260805.json`](../benchmarks/ccp/edupic-reference-stationarity-cycle1500-20260805.json).

The first native measurement-mode pilot then advanced the immutable
cycle-1,500 checkpoint by exactly one RF cycle. The external process was
pinned to one low-priority CPU, bounded by 886,876,000 initial
particle-timesteps and a 60 s timeout, and completed in 23.38 s. It produced
all expected density, EEPF, IFED, phase-resolved field, density, current,
power, energy, and ionization outputs. Every table has its upstream shape and
finite values. EEPF normalization is `0.99999997`; powered and grounded IFED
normalizations are `1.000003` and `0.999990`.

All four upstream stability checks pass: plasma-frequency times electron
timestep is `0.090`, grid spacing over central Debye length is `0.723`, and
maximum electron and ion collision-frequency products are `0.013` and
`0.020`. The reported one-cycle density, flux, energy, and power values are
descriptive runtime/output qualification only; one RF cycle cannot supply the
recommended measurement statistics or support cross-code agreement. The
upstream accumulators are process-local and absent from `picdata.bin`, so the
next safety gate is an immutable measurement-block coordinator plus an
explicit block aggregation/uncertainty contract, rather than an uncheckpointed
multi-hour process. Evidence is
[`edupic-reference-measurement-pilot-cycle1501-20260806.json`](../benchmarks/ccp/edupic-reference-measurement-pilot-cycle1501-20260806.json).

The new measurement coordinator has qualified real immutable resume through
cycle 1,508. Two independent four-cycle measurement processes form one
contiguous checkpoint chain from the original stationary cycle-1,500 state.
They completed in 93.36 and 92.82 s, compared with a conservative 140.31 s
prediction, and each invocation stopped after its declared single-stage cap.
Every stage is replayed against its checkpoint, convergence history, native
diagnostic hashes, binary hash, and input-state hash before resume.

Four pre/post-stage host checks all passed. Available memory remained above
4,735 MiB, normalized one-minute load remained below 0.362, and swap I/O was
zero. The two blocks retain valid EEPF and powered/grounded IFED normalization.
Their central electron densities differ by 0.264% relative to their mean and
total power densities by 1.19%. These two-block ranges are encouraging
operational consistency, not measurement convergence or uncertainty. The
upstream process-normalized EEPF, IFED, conditional moments, and derived power
arrays also require an explicit aggregation definition; they must not be
silently treated as raw pooled accumulators. Evidence is
[`edupic-reference-measurement-blocks-through-cycle1508-20260806.json`](../benchmarks/ccp/edupic-reference-measurement-blocks-through-cycle1508-20260806.json).

The predeclared short-horizon qualification subsequently reached cycle 1,516:
four equal four-cycle measurement blocks, four invocations, and 369.98 s of
solver time. Eight host checks passed with available memory above 4,735 MiB,
normalized load below 0.362, and zero swap. A new analyzer replays the complete
hash chain and makes the upstream aggregation semantics explicit. Density,
potential, field, species density/current, and ionization-rate arrays are
linear duration averages. EEPF and IFED are equal-time mixtures of already
normalized block distributions rather than native pooled histograms.
Conditional mean energies and products of separately averaged current and
field remain blockwise because the necessary raw weights are absent.

Across the four blocks, electron density profiles differ from their exact
duration aggregate by 0.31--0.41% relative L2 and ion profiles by 1.04--1.14%.
EEPF total variation from the normalized mixture is 0.42--0.57%. Total power
density has a 1.21% range around its block mean. In contrast, powered IFED
variation is 9.79--12.19% and grounded IFED variation is 10.21--10.96%, so
wall-energy statistics clearly require a longer horizon. Naive block standard
errors are retained but explicitly not corrected for autocorrelation. This is
an eligible short-horizon consistency analysis, not cross-code or physical
validation. Evidence is
[`edupic-reference-measurement-block-analysis-cycle1516-20260806.json`](../benchmarks/ccp/edupic-reference-measurement-block-analysis-cycle1516-20260806.json).

Before collecting any later native measurement block, the next sampling rule
was frozen in
[`edupic-native-measurement-stationarity-rule-20260806.json`](../benchmarks/ccp/edupic-native-measurement-stationarity-rule-20260806.json).
It requires 16 contiguous four-cycle blocks (64 measured cycles), at least
eight AR(1)-effective blocks for each species' line-integrated density, no
more than 1% projected or split-half density drift, and no adjacent density
profile movement above 2.5% relative L2. These are internal sampling-readiness
gates aligned with the existing Turner density-series audit, not published
eduPIC thresholds. EEPF and IFED remain descriptive because the upstream
files discard the raw observation counts needed for a defensible pooled
uncertainty estimate. The analyzer can now join separately immutable campaign
directories only when their source-binary and checkpoint hashes form a
continuous chain.

The first predeclared continuation block then reached cycle 1,520 in 91.29 s.
Its single-core host checks retained at least 5,740 MiB available memory,
normalized load at or below 0.269, and zero swap I/O. Joined with the original
four blocks, the five-block density series remains well inside the numerical
drift and profile-movement thresholds, but it is explicitly ineligible because
the 16-block horizon is incomplete. The current negative lag-one estimates cap
the nominal AR(1) count at five and are not reliable at this sample count.
Latest-block EEPF total variation is 0.67%, while powered and grounded IFED
variation remains 11.22% and 12.81%. Progress evidence is
[`edupic-native-measurement-continuation-cycle1520-20260806.json`](../benchmarks/ccp/edupic-native-measurement-continuation-cycle1520-20260806.json).

Three further single-stage invocations reached the eight-block midpoint at
cycle 1,532. Density drift remains small: projected series changes are 0.026%
for electrons and 0.033% for ions, while maximum adjacent-profile movement is
0.59% and 1.93%, respectively. The longer series exposes positive serial
correlation, however. Lag-one values of 0.303 and 0.344 reduce the nominal
AR(1) counts to only 4.28 electron and 3.91 ion blocks. Thus treating the
eight block means as independent would materially overstate confidence.
Latest-block EEPF variation is 0.54%; IFED remains near 11.5%.

The cycle-1,532 stage itself completed and passed all output and checkpoint
checks. Its after-stage host check then measured normalized load 0.627, above
the frozen 0.5 limit, while unrelated host workloads were active. The
coordinator stopped with `host_load_above_maximum`; no process was terminated
and no ninth block was attempted. Midpoint evidence is
[`edupic-native-measurement-continuation-cycle1532-20260807.json`](../benchmarks/ccp/edupic-native-measurement-continuation-cycle1532-20260807.json).

After host load returned below policy, three guarded invocations advanced the
series through cycle 1,544 and 11 total blocks. All six surrounding host
checks passed, with at least 4,764 MiB available memory, normalized load no
higher than 0.225, and zero swap. Electron and ion lag-one correlations are
now 0.282 and 0.256, giving 6.16 and 6.52 AR(1)-effective blocks. Projected
density changes remain only -0.015% and -0.009%, and maximum adjacent-profile
movement remains below 0.69% and 1.94%. The result is encouraging but still
ineligible: five blocks remain, and neither species has reached the frozen
eight-effective-block gate. IFED variation remains 12.2--13.9%. Evidence is
[`edupic-native-measurement-continuation-cycle1544-20260807.json`](../benchmarks/ccp/edupic-native-measurement-continuation-cycle1544-20260807.json).

Five final guarded stages completed the frozen horizon at cycle 1,564: 16
contiguous four-cycle blocks and 64 measured cycles. A dedicated evaluator now
applies the committed rule without manual threshold interpretation. The
horizon, drift, split-half, and adjacent-profile gates all pass. Electron and
ion effective counts are only 6.24 and 7.53, however, below the required eight
for both species. The formal classification is therefore
`internal_density_stationarity_screen_failed`.

This is a useful negative result rather than a solver failure. Projected
density changes are only -0.015% and -0.009%; split-half changes are -0.021%
and -0.017%; maximum adjacent-profile movement is 0.69% and 1.93%. The failed
effective-count gates instead show that four-cycle block means remain
correlated. EEPF block variation is 0.49--0.75%, while powered and grounded
IFED ranges remain 10.4--14.1% and 11.6--14.5%; no IFED acceptance threshold
is applied. The correct next step is a predeclared longer horizon, not a
physical-validation claim. Complete evidence is
[`edupic-native-measurement-stationarity-cycle1564-20260807.json`](../benchmarks/ccp/edupic-native-measurement-stationarity-cycle1564-20260807.json).

Before collecting cycle 1,565 or later, a follow-up rule was committed in
[`edupic-native-measurement-extension-rule-20260807.json`](../benchmarks/ccp/edupic-native-measurement-extension-rule-20260807.json).
It retains the four-cycle blocks and every density threshold, but extends the
checkpoint to 24 total blocks and 96 measured cycles through cycle 1,596. At
the observed electron lag-one correlation of 0.439, the AR(1) approximation
requires about 21 total blocks to reach eight effective blocks; 24 provides
some margin without being a guarantee of passage. The original 64-cycle
failure remains unchanged regardless of the extension outcome.

The first three extension blocks reached cycle 1,576 and 19 total blocks.
Every host check passed, with at least 4,519 MiB available memory, normalized
load no higher than 0.331, and zero swap. Contrary to the extension's planning
estimate, lag-one correlation increased to 0.671 for electron density and
0.631 for ion density. Nominal effective counts consequently fell to 3.74 and
4.30 even though three blocks were added. Projected density changes remain
only 0.038% and 0.050%, and adjacent-profile movement remains below 0.69% and
2.06%. This is evidence of longer correlation, not observed density drift.
The 24-block rule remains incomplete and is not assumed to pass. Evidence is
[`edupic-native-measurement-extension-cycle1576-20260807.json`](../benchmarks/ccp/edupic-native-measurement-extension-cycle1576-20260807.json).

The high AR(1) correlation was then triaged as a measurement-analysis issue,
not assumed to be a solver defect. Across 19 blocks the line-integrated
density varies by less than +/-0.08%. AR(1) correction raises the relative
standard error from 0.0080% to 0.0179% for electrons and from 0.0086% to
0.0180% for ions. Non-overlapping eight-cycle batches reduce lag-one
correlation to 0.160 and 0.062, with relative standard errors of 0.0101% and
0.0105%. Thus correlation matters, but absolute density-mean uncertainty is
still very small. The analyzer now reports AR(1)-corrected uncertainty,
variance inflation, and power-of-two reblocking diagnostics. This additional
post-hoc view does not change either frozen gate. Evidence is
[`edupic-native-measurement-reblocking-diagnostic-cycle1576-20260807.json`](../benchmarks/ccp/edupic-native-measurement-reblocking-diagnostic-cycle1576-20260807.json).

No faulty field equation, particle mover, collision model, or boundary model
has been identified by this result. Changing one would be scientifically
unjustified. The remaining external IFED limitation is different: eduPIC's
normalized output omits raw observation counts, so defensible pooled
uncertainty cannot be reconstructed afterward. AuroraPIC's native phase-EEDF
diagnostics already preserve macro and represented observation counts; the
external IFED comparison remains descriptive until count-preserving source
output is available.

The extension subsequently completed its frozen endpoint at cycle 1,596: 24
four-cycle blocks and 96 measured cycles. It fails the effective-count gates
again, now with raw AR(1) counts of 1.00 electron and 1.57 ion blocks. Every
other frozen gate passes. Projected density changes are 0.181% and 0.175%,
split-half changes are 0.084% and 0.083%, and maximum adjacent-profile
movement is 0.75% and 2.06%. AR(1)-corrected density-mean standard errors
remain below 0.077% and 0.058%.

Reblocking confirms that the long component is not removed by merely pairing
four-cycle windows: even six non-overlapping 16-cycle means retain lag-one
correlations of 0.453 and 0.503. Therefore another arbitrary short extension
would be methodologically weak. The next measurement design must use
substantially longer batches and a production-scale horizon informed by
eduPIC's approximately 1,000-cycle recommendation. The completed negative
result is recorded in
[`edupic-native-measurement-extension-cycle1596-20260807.json`](../benchmarks/ccp/edupic-native-measurement-extension-cycle1596-20260807.json).

The production-scale replacement was frozen before cycle 1,597 in
[`edupic-native-production-measurement-rule-20260807.json`](../benchmarks/ccp/edupic-native-production-measurement-rule-20260807.json).
It starts a statistically separate 1,024-cycle segment at the cycle-1,596
checkpoint, satisfying the upstream at-least-1,000-cycle measurement guidance.
Sixty-four native 16-cycle stages support exact density aggregation and
power-of-two reblocking through 128-cycle means. Prior short windows retain
checkpoint provenance but are excluded from the formal production statistics.
At the qualified rate the campaign costs approximately 6.65 serial solver
hours and 1.3 GiB of retained evidence.

The measurement coordinator now also samples free filesystem capacity before
and after every stage. The production contract requires at least 32 GiB free,
in addition to the existing one-core, 4 GiB available-memory, normalized-load
0.5, zero-swap, 600 s timeout, and one-stage-per-invocation guards. This closes
a production-safety gap on the currently 88%-used workspace filesystem.

The first production invocation was refused before execution when available
memory fell to 3,655 MiB; no simulation state changed. After the unrelated
host workload released memory, the identical immutable contract completed
cycles 1,597--1,612 in 362.11 s. The admitted stage retained at least 5,234
MiB available memory and 47,889 MiB free disk, with normalized load no higher
than 0.453 and zero swap. The observed 22.63 s/cycle projects to about 6.44
solver hours for the full campaign.

After 25 valid stages (400 production cycles), the operational available-memory
floor was reduced from 4,096 to 3,072 MiB. This was an audited safety-policy
amendment, not a change to the checkpoint chain or scientific analysis contract.
The solver remained pinned to one CPU with approximately 45 MiB resident memory,
zero swap I/O, and a 600 s stage timeout. Across those stages, available memory
bottomed at 3,701 MiB and the largest observed within-stage decline was about
1,180 MiB. The revised 3 GiB floor therefore avoids repeated false-positive
stops while retaining margin for unrelated host-memory variability.

The single-block preview has valid EEPF and powered/grounded IFED
normalizations, but is correctly ineligible for statistics. An analyzer edge
case discovered by this preview was fixed: single-block variance,
autocorrelation, and standard-error fields are now `null`, and no batch series
is invented. Evidence is
[`edupic-native-production-measurement-cycle1612-20260807.json`](../benchmarks/ccp/edupic-native-production-measurement-cycle1612-20260807.json).

The production campaign subsequently reached its exact cycle-2,620 target:
64 contiguous 16-cycle blocks and 1,024 standalone measurement cycles completed
in 22,851.00 s of single-core solver time. Every checkpoint and native output
hash revalidated, the final checkpoint SHA-256 is
`3a914d71c9740ae12c67daebb45d6199b32f42e45941ed27af1433a7b28f9dc4`,
and swap I/O remained zero. Four host-memory checks stopped coordination safely
under the original 4 GiB floor; no memory stop occurred after the audited 3 GiB
amendment.

The complete horizon passes the predeclared drift, split-half, and adjacent
profile-movement gates for both densities. Electron and ion projected drift
magnitudes are 0.0800% and 0.0775%, and their autocorrelation-corrected relative
standard errors are 0.0716% and 0.0680%. Nevertheless, adjacent 16-cycle density
blocks have lag-one correlations of 0.9003 and 0.8998, leaving only 3.36 and
3.37 AR(1)-effective blocks against the required eight. The correct frozen
classification is therefore `internal_density_stationarity_screen_failed`.
Small corrected uncertainty does not authorize changing a predeclared gate
after observing the result.

This was an eduPIC execution, not an AuroraPIC execution. It does not validate
AuroraPIC or establish cross-code agreement. The checksum-pinned result is
[`edupic-native-production-measurement-cycle2620-20260810.json`](../benchmarks/ccp/edupic-native-production-measurement-cycle2620-20260810.json).
The next cross-code prerequisite is an AuroraPIC argon contract matching the
reference cross sections, ionization energy partition, scattering, boundary
drive, timestep/subcycling, and observable definitions.

That prerequisite has now advanced to an executable contract preflight.
`import_edupic_cross_sections.py` consumes the one-million-row table generated
by the pinned upstream `test_cross_sections()` routine without embedding its
GPL formulas in AuroraPIC. It validates the exact `0.001 eV` grid and writes
separate, checksum-audited local electron and ion gas manifests. The electron
manifest selects isotropic elastic/excitation and the Opal-style `10 eV`
ionization distribution. The ion manifest selects center-of-mass isotropic and
backward scattering. Generated tables remain ignored and are not redistributed.

The locked case manifest is
[`edupic_argon_ccp_reference.case`](../examples/edupic_argon_ccp_reference.case).
Its preparer validates every local package hash and can generate at most one RF
cycle; it explicitly cannot authorize production. A two-step, one-core SI
preflight loaded all five million table rows and advanced the exact 400-node,
4,000-step-per-cycle geometry/drive contract in 7.06 s with a 177,880 KiB peak
resident set and zero swap. This is integration evidence only. The
checksum-bearing record is
[`edupic-argon-aurorapic-contract-preflight-20260810.json`](../benchmarks/ccp/edupic-argon-aurorapic-contract-preflight-20260810.json).

The cross-section evaluator now provides a fingerprinted `lower_bin` mode,
selected by all five generated eduPIC channels. This matches the reference's
integer selection on its `0.001 eV` source grid while preserving linear
interpolation as the default for every existing gas package.
Thermal-neutral frequency validation uses an exact segment-tree range maximum,
reducing each table query from a million-row scan to logarithmic work. The
index raises this five-table preflight's peak memory by about 80 MiB, but the
optimized interpolation-only run was bit-for-bit identical to its
pre-optimization run and reduced its measured wall time from 11.38 s to
7.09 s.

The generated AuroraPIC deck now gives ions a timestep multiplier of 20. Ion
pushes, boundary checks, and MCC calls occur at pre-step indices divisible by
20 and use `20*dt`, including the reference-compatible update at index zero;
ion charge is held between those updates. Checkpoint v14 retains the v13
species-schedule validation. The new two-step preflight therefore exercises
one long ion update and completes in 7.06 s with zero swap.

The previously explicit inelastic-electron transform difference is now closed:
the generated v2 gas manifest selects the finite-mass center-of-mass transform
for excitation and ionization, and the importer audit and case preparer enforce
that choice. Unit regressions recover the post-threshold relative velocity for
excitation and both Opal electron velocities after the transform. This closes
an input/kinematic contract gap; it is not yet a discharge-scale validation.
Independent random streams and microscopic initial states are expected for a
black-box statistical comparison and are not defects. The failed
external effective-block gate also remains a limitation on formal uncertainty,
not something this preflight repairs.

The preparer now has a strict `--startup-diagnostics` mode. It accepts only the
complete 4,000-step RF period, emits 40 scalar samples, accumulates spatial
collision energy every step, and assigns exactly 250 steps to each of 16 phase
bins. The first one-cycle serial screen completed in 10.08 s with a 171,728 KiB
peak resident set and zero swap. It sampled 75,134 collision candidates,
including 1,703 excitations and 2,381 ionizations. Particle accounting closed
exactly: the 2,381 created pairs and wall losses led from 1,000 particles per
species to 863 electrons and 3,369 ions. The global kinetic-energy ledger
closed to `3.30e-13 W m^-2` (`1.18e-15` relative), while global-to-spatial and
spatial-to-phase collision-energy residuals stayed below
`6.36e-21 J m^-2`.

This is positive internal numerical evidence, not equilibrium or external
validation. The ion inventory and field energy were still increasing at the
cycle boundary, so the discharge was plainly forming. The checksum-bearing
record is
[`edupic-argon-aurorapic-startup-cycle1-20260810.json`](../benchmarks/ccp/edupic-argon-aurorapic-startup-cycle1-20260810.json).
The next authorized development slice is a checkpointed multi-cycle pilot with
per-cycle population, field, collision, energy-closure, and resource gates.

That pilot is now implemented by `run_aurorapic_edupic_pilot.py`. It accepts
only the checksum-locked cycle-1 deck, executable, and checkpoint; runs one
low-priority serial process per cycle; resets the complete spatial/phase window
at each restart; and has non-relaxable limits of cycle 4, 60 seconds per cycle,
512 MiB available-memory launch floor, 4x per-cycle particle growth, 25% of the
configured particle cap, `10^7 V m^-1`, `10^-10` relative energy residual, and
`10^-15 J m^-2` spatial/phase residual. Every stage preserves its input/output
checkpoint hashes and population/event balance.

The cycle-2 through cycle-4 continuation passed every gate. Total-particle
growth moderated from 1.277x to 1.237x to 1.195x; electron counts at the cycle
boundaries were 1,329, 1,886, and 2,483, while ion counts were 4,077, 4,801,
and 5,506. Ionization settled near 766--792 events per cycle after the 2,381
event first-cycle avalanche. Peak sampled field still rose from 31.2 to 33.6
to 35.2 kV/m, and field energy rose each cycle. Relative energy closure stayed
below `2.8e-14`, spatial/phase residuals below `9.8e-21 J m^-2`, peak resident
memory near 172 MiB, and each continuation completed in under 12 seconds.

This is a bounded, credible equilibration trend—not stationarity. The retained
record is
[`edupic-argon-aurorapic-multicycle-pilot-cycle4-20260810.json`](../benchmarks/ccp/edupic-argon-aurorapic-multicycle-pilot-cycle4-20260810.json).
The next step must be a separately authorized blockwise horizon extension with
cycle-boundary stationarity gates; equilibrium observable comparison remains
premature.

`extend_aurorapic_edupic_horizon.py` now provides that extension in immutable
four-cycle blocks, with an absolute cycle-16 ceiling. In addition to the pilot
safety gates, it freezes an internal stationarity screen before execution:
absolute normalized slopes below 1% per cycle for total population, field
energy, and peak field; below 2% for ionization; and ionization coefficient of
variation below 5%. The tool verifies the prior report and input checkpoint
hashes, resets each full-cycle diagnostic window, and emits a new checkpoint
chain. A failed stationarity screen is retained as evidence and does not
relabel safe execution as equilibrium.

The cycles 5--8 block passed every hard execution/accounting gate but correctly
failed stationarity. Ionization passed both trend gates, with 1.27% normalized
slope and 1.87% coefficient of variation. Total population still had an 11.84%
normalized slope per cycle, field energy 2.97%, and peak field 2.32%. At cycle
8 the simulation contained 4,841 electrons and 8,108 ions; peak field was
38.66 kV/m. Each cycle stayed below 16 seconds and about 172 MiB, relative
energy residuals below `3.5e-14`, and spatial/phase residuals below
`1.8e-20 J m^-2`. The checksum-bearing negative stationarity result is
[`edupic-argon-aurorapic-horizon-cycle8-20260810.json`](../benchmarks/ccp/edupic-argon-aurorapic-horizon-cycle8-20260810.json).
The horizon runner now accepts either the completed safe pilot or a completed
safe prior horizon report, while rejecting an unsafe report. This permits an
immutable multi-block checkpoint chain without weakening any execution gate.
For cumulative step counts that trigger the CLI's conservative initial-update
estimate, the runner supplies the explicit large-run acknowledgement. Actual
incremental work remains limited to one 4,000-step, low-priority, single-core
cycle at a time under the runner's non-relaxable timeout and resource gates.

The chained cycles 9--12 block also passed every hard execution and accounting
gate and correctly failed stationarity. The normalized total-population slope
moderated from 11.84% to 7.97% per cycle, field-energy slope from 2.97% to
1.72%, and peak-field slope from 2.32% to 1.86%. All three remain above their
frozen 1% thresholds. Ionization passed both of its gates: its normalized
slope was -1.93% per cycle and its coefficient of variation was 3.32%. At
cycle 12 the simulation contained 7,183 electrons and 10,667 ions, with a
41.57 kV/m peak sampled field. Each cycle completed in under 19.4 seconds at
about 172 MiB peak resident memory; relative energy residuals remained below
`3.8e-14` and spatial/phase residuals below `3.2e-20 J m^-2`.

These improving trends support continued bounded equilibration, but they do
not establish equilibrium or agreement with eduPIC. The checksum-bearing
negative result is
[`edupic-argon-aurorapic-horizon-cycle12-20260810.json`](../benchmarks/ccp/edupic-argon-aurorapic-horizon-cycle12-20260810.json).
The next and final permitted block is cycles 13--16 under the unchanged frozen
thresholds.

Cycles 13--16 completed that predeclared horizon with every hard gate passing
and with another correct stationarity failure. The normalized population slope
moderated to 6.00% per cycle, field-energy slope to 1.36%, and peak-field slope
to 1.30%; each remains above the frozen 1% limit. Ionization passed with a
0.96% normalized slope and 1.39% coefficient of variation. At cycle 16 the
simulation contained 9,502 electrons and 13,205 ions, and the peak sampled
field was 43.85 kV/m. Per-cycle runtime remained below 23.4 seconds and peak
resident memory near 172 MiB. Relative energy closure stayed below `7.1e-14`
and spatial/phase closure below `1.2e-19 J m^-2`.

The failed stationarity result means that no equilibrium observable or eduPIC
comparison is yet eligible. It also exhausts the original absolute cycle-16
ceiling. The next step is to freeze a new, longer blockwise equilibration
horizon before examining more endpoints, retaining the same thresholds and
resource gates. The evidence is
[`edupic-argon-aurorapic-horizon-cycle16-20260810.json`](../benchmarks/ccp/edupic-argon-aurorapic-horizon-cycle16-20260810.json).

Before observing cycle 17, the follow-up equilibration and comparison-readiness
contract was frozen in
[`edupic-argon-aurorapic-equilibration-extension-rule-20260810.json`](../benchmarks/ccp/edupic-argon-aurorapic-equilibration-extension-rule-20260810.json).
It extends the absolute ceiling to cycle 64 only when its exact checksum is
provided. Execution remains limited to one low-priority four-cycle block per
invocation with all existing safety and stationarity thresholds unchanged.
Two consecutive blocks must pass the internal stationarity screen before a
separate measurement campaign becomes eligible. That future campaign excludes
all equilibration samples and targets 1,024 cycles in 16-cycle blocks to match
the duration of the retained native eduPIC reference measurement.

The comparison contract requires electron and ion density profiles, EEPF, and
powered/grounded ion-impact energy distributions. Density and electron-energy
diagnostics already have AuroraPIC output paths. A count-preserving 1D ion
wall-impact spectrum remains an implementation prerequisite. The external
eduPIC reference's failed effective-block gate also remains visible, so the
first eligible cross-code result will be descriptive rather than a formal
acceptance claim.

The first prospectively governed block, cycles 17--20, passed every hard gate
but did not start the required stationarity streak. Population slope moderated
to 4.64% per cycle and remains the dominant failure. Field-energy slope was
1.044%, narrowly above its 1% gate, while peak-field slope passed for the first
time at 0.984%. Ionization slope and coefficient of variation passed at 0.171%
and 0.632%. Cycle 20 ended with 11,734 electrons, 15,603 ions, and a 45.67 kV/m
peak sampled field. The checksum-bearing result is
[`edupic-argon-aurorapic-horizon-cycle20-20260810.json`](../benchmarks/ccp/edupic-argon-aurorapic-horizon-cycle20-20260810.json).
Measurement remains ineligible; the next block is cycles 21--24.

Cycles 21--24 also passed every hard gate but did not start the stationarity
streak. Population slope moderated further to 3.93% per cycle. Field-energy
slope passed at 0.838%, while peak-field slope narrowly failed at 1.022%.
Ionization coefficient of variation passed at 4.16%, but its normalized slope
failed at 3.68% because the four cycle counts rose from 738 to 821. Cycle 24
ended with 13,971 electrons, 18,003 ions, and a 47.56 kV/m peak sampled field.
The retained result is
[`edupic-argon-aurorapic-horizon-cycle24-20260810.json`](../benchmarks/ccp/edupic-argon-aurorapic-horizon-cycle24-20260810.json).
Measurement remains ineligible; the next block is cycles 25--28.

The post-processing wrapper originally allowed 15 seconds per analyzer. During
the later comparison-readiness horizon, unrelated host CPU and I/O contention
caused a valid energy analysis to take 12.34 seconds on retry and one attempt
to cross that narrow wrapper limit. The analyzer allowance is now 45 seconds.
This is an operational robustness change only: the solver remains capped at
60 seconds per cycle and no execution, physics, accounting, or stationarity
threshold changed.

The prospectively declared extension subsequently completed all ten remaining
four-cycle blocks through cycle 64. Every solver, resource, population-cap,
field, energy-closure, spatial-closure, and checkpoint-chain gate passed. From
cycle 28 to cycle 64, normalized population slope declined from 3.30% to
1.159% per cycle. At the final endpoint, field-energy slope was 0.135%, peak-
field slope 0.286%, ionization slope 0.042%, and ionization coefficient of
variation 1.07%; all four pass their frozen gates. Population alone remains
above its 1% threshold, so no block passed the complete stationarity screen
and the required two-block streak is zero.

This is physically more than a marginal statistical miss. Cycle 64 recorded
715 ionizations but only 336 electron and 308 ion wall losses, producing net
gains of 379 electrons and 407 ions. The discharge therefore remains in a
filling transient. It ended with 32,549 electrons, 37,229 ions, 45.41
`microJ m^-2` field energy, and a 54.78 kV/m peak sampled field. A side-by-side
measurement must not begin from this state.

The late campaign also exposed its performance boundary: single-core cycle
time reached 56.39 seconds at 69,778 particles, close to the unchanged
60-second guard, while RSS remained only about 171 MiB. One contended-core
cycle-60 attempt timed out; a rerun on a quieter single core reproduced the
cycle 57--59 checkpoint hashes and completed without changing the limit. The
predeclared horizon is exhausted. Complete evidence is
[`edupic-argon-aurorapic-equilibration-cycle64-20260811.json`](../benchmarks/ccp/edupic-argon-aurorapic-equilibration-cycle64-20260811.json).

The next credible step is not another short extension. AuroraPIC first needs a
profile-guided 1D3V performance pass and count-preserving ion wall-impact
spectra; a new production-scale equilibration horizon should then be declared
before further cycles, informed by eduPIC's approximately 1,500-cycle startup.

The first performance prerequisite is now complete. Reusable per-model MCC
workspaces remove per-particle vector allocation, and ionization products use
a monotonic dead-slot cursor rather than repeated full storage scans. An exact
single-core replay from the cycle-63 checkpoint completed cycle 64 in 45.85
seconds, down from the certified 56.39 seconds. Its final checkpoint SHA-256
remained `f99b58a0b39a04c190e5ee9b4d5b98d2a65f0cdb9bf42f6165cf1d745541d47c`,
and ten emitted diagnostic products were byte-identical. This advanced the
runtime prerequisite without changing the failed stationarity conclusion; at
that milestone, ion wall-impact spectra remained the next prerequisite.

That diagnostic prerequisite is now implemented. The opt-in keys are:

```ini
wall_impact_spectrum = true
wall_impact_reset_on_restart = true
wall_impact_energy_bins = 200
wall_impact_energy_max = 500
```

Energy bounds are eV in SI runs and normalized energy otherwise. The runtime
accumulates separate left/right histograms for every species, recording macro
counts and represented counts independently. Overflow is never folded into
the last bin. `wall_impact_spectrum_summary.csv` closes every species/side
macro count and represented kinetic-energy total against the pre-existing
boundary-loss ledger; a mismatch is fatal. Checkpoint v14 retains histogram,
overflow, closure baseline, and origin state. Enabling the diagnostic while
restarting a pre-v14 checkpoint begins a new origin-labelled window without
discarding the older cumulative boundary ledger.

`wall_impact_reset_on_restart = true` explicitly begins a fresh histogram
window at a v14 checkpoint while retaining the cumulative boundary ledger as
the new closure baseline. It is intended for prospectively declared
measurement windows; leave it false on later restarts that continue the same
window.

An exact cycle-64 replay from the cycle-63 checkpoint exercised the diagnostic
with 200 bins over 0--500 eV. It captured all 336 electron impacts and all 308
ion impacts with no overflow. Mean ion impact energies were 29.56 eV at the
left electrode and 26.20 eV at the right; these are nonstationary pilot values,
not comparison measurements. All four count and energy closures passed, the
ten pre-existing diagnostic products stayed byte-identical, and the particle
plus RNG checkpoint tail stayed byte-identical. The retained evidence is
[`edupic-argon-aurorapic-cycle64-wall-impact-pilot-20260811.json`](../benchmarks/ccp/edupic-argon-aurorapic-cycle64-wall-impact-pilot-20260811.json).

With performance and wall-impact diagnostics complete, the next scientific
step is to predeclare a longer equilibration campaign. Its measurement gate
remains unchanged: two consecutive stationary blocks must pass before any
AuroraPIC-to-eduPIC distribution comparison begins.

Commit `df8765d` added restart-safe, species-resolved electric-work
accounting. The diagnostic records represented kinetic-energy change caused
only by the electric particle push; collision-energy transfer and kinetic
energy carried into a wall remain separate. A collisionless no-loss regression
requires this work to equal the particle kinetic-energy change. Ordered
per-worker reductions preserve deterministic execution, checkpoint v7
preserves cumulative values and their origin, and a v6 restart begins the new
counter at its restart step.

A second 32-cycle continuation covered steps 525,200--538,000. Its quantities
are directly comparable in definition to the line-integrated,
time-averaged species `J · E` powers reported in
[Turner et al., Table III](https://arxiv.org/abs/1211.5246), but this window is
still beyond the published 1,280-cycle duration and has no published
single-window uncertainty band. The retained local artifact identities are:

| Artifact | SHA-256 |
| --- | --- |
| Solver executable (commit `df8765d`) | `4ba77160ef71c4f2630139f42c2aa0e77662a2d8f3c11c5aa11d03419849a668` |
| Input step-525,200 v6 checkpoint | `23c82ef6ce77a8acf9b4465d9c38aab0da38621da0dfd89a1c9e3fc54061cd95` |
| 32-cycle power deck | `b05efe0002434e3d5e6f9a1c137b57178c7f0783ef1ac8b937f07ce197040659` |
| Final step-538,000 v7 checkpoint | `b7f67ed11e41611025edb51ab537b6ea6008f721a82fcf9062c25af1b8d3f5ff` |
| 32-cycle power/balance report | `d010cbcc22a04bd5591a7d705a914b06b80634548aad3f281a52f2f8a8b68b35` |

The measured electron power is `34.2044 W m^-2`, 0.279% below Turner's
`34.3 W m^-2`. The measured ion power is `91.5628 W m^-2`, 1.063% above
Turner's `90.6 W m^-2`. The mean two-electrode ion-current magnitude is
`0.220216 A m^-2`, 0.555% above `0.219 A m^-2`. The window also has exact
integer source/loss closure for both species: 24,822 ionizations, electron
population change -14, ion population change +35, and zero balance residuals.

Agreement of three independent global observables at about 1% or better is
strong evidence that the discharge's field heating, ion acceleration,
ionization source, and absorbing-wall flux are mutually consistent. It is
materially stronger than the density comparison alone, but it is not relabeled
as a formal published-duration pass: the original density-profile `X²` gate
remains failed and this continuation is diagnostic evidence.

### Post-benchmark density blocks

Commit `826097a` added an explicit
`spatial_average_reset_on_restart = true` contract. It discards checkpointed
averaging sums only when requested, rejects a new window whose first sample is
not after the checkpoint step, records the reset in metadata, and reproduces
an uninterrupted late window byte-for-byte in regression. The Turner
comparator's `--post-benchmark-window` mode requires a reset, complete,
32-cycle window after the published duration and sets
`published_acceptance_applicable = false`.

Two consecutive 32-cycle blocks were then measured: cycles 1,314--1,345
(steps 525,201--538,000) and cycles 1,346--1,377
(steps 538,001--550,800). They follow one trajectory and are non-overlapping,
but are not statistically independent. Key retained artifact identities are:

| Artifact | SHA-256 |
| --- | --- |
| Solver executable (commit `826097a`) | `ed6c8c31f73754cbc5851678ade9ef85d2445dd456cefe271cb9e8c0390ef072` |
| Block 1 deck | `d7fd52a4900f4bc37e4729fd5ef9ef7453b44637672474323ec591ce329ada19` |
| Block 1 density profile | `f02335735ab7e6f2501daf1aa352d2c0da35d89a2859e00ff9496455b1dfc139` |
| Block 1 comparison report | `d3aea6cb4e1626444dbf2fb5d923511687f97a40319e409510adc714340f9194` |
| Step-538,000 v7 checkpoint | `51f26168508bc090274af85b2bd8c3d7e38e09b63ebf1774c9bacdd064ffbaed` |
| Block 2 deck | `9fa4b7cfb1f361d471e883f462c85ac5c52565cd498f9a9fafffb85d00588eed` |
| Block 2 density profile | `b648f2496e18adc905ad86310b3fee5df172801b15ea29a5d5ed3856ddfaec2f` |
| Block 2 comparison report | `58fadc3444f1efdad922c93948f8459cf2524e1355caa10f19964e9bd4543cd0` |
| Final step-550,800 v7 checkpoint | `d81837888681e7822a22d02d912888e92321a1a1d25bea01102a35c8e20d83b9` |

| 32-cycle window | X² | Relative L2 | Maximum pointwise error | Published gate applicable |
| --- | ---: | ---: | ---: | --- |
| Published-duration final block | 574.399 | 2.874% | 4.811% | Yes: failed 99% |
| Post-benchmark block 1 | 365.343 | 2.362% | 4.638% | No |
| Post-benchmark block 2 | 849.284 | 3.478% | 5.717% | No |

Block 1 happens to lie inside the published 99% numerical range and block 2
lies outside it, but neither is an acceptance result. Relative profile
movement is 1.573% L2 from the published-duration block to block 1 and 1.827%
from block 1 to block 2. Their line-integrated ion densities change by -0.689%
and +1.261%, respectively. Thus the original density excess is not a stable,
monotonic offset over these late windows: stochastic sampling and/or a
slowly varying discharge state is comparable in size to the reported
discrepancy. This weakens the case for a gross systematic physics error, but
does not distinguish correlated noise from insufficient stationarity. A
formal resolution requires either a multi-seed published-duration ensemble or
enough consecutive blocks to estimate autocorrelation and effective sample
size.

The restart-safe sequential-block audit is now implemented by
`scripts/analyze_turner_density_blocks.py`. It rejects profile or metadata
hash drift, mixed cases/species, incomplete or non-reset windows, changed
spatial grids, and any gap or overlap between blocks. A minimum of eight
consecutive 32-cycle blocks is the declared floor before the series is even
eligible for stationarity interpretation; this is an internal diagnostic
floor, not a threshold published by Turner et al. The tool always leaves
`published_acceptance_applicable = false`.

Applied initially to the first two consecutive blocks, the analyzer confirmed
a 1.827% relative profile movement and a 1.261% rise in line-integrated ion
density. That two-block audit report has SHA-256
`77071ae543d658b967d7adb6e273989f4a84243f3f3857fe01b28bcbfac189a9`.

A third low-priority serial continuation then covered cycles 1,378--1,409
(steps 550,801--563,600). Its comparison gives `X² = 484.785`, 3.335%
relative profile L2 error, and 6.532% maximum pointwise error. The
line-integrated density falls 0.526% from block 2 and the adjacent profile
movement is 1.268%. Across all three blocks the fitted end-to-end integrated
density drift is +0.723%; the preliminary lag-one correlation is -0.446.
With only three samples neither value is a reliable stationarity estimate,
and the analyzer correctly retains `insufficient_consecutive_blocks`.

The third window also has exact integer source/loss closure: 24,744
ionizations, electron and ion population changes of -74 and -85, and zero
balance residuals. Its two-electrode mean ion-current magnitude is
`0.220589 A m^-2`, 0.726% above Turner; electron and ion electrical powers
are `34.1306 W m^-2` and `91.6173 W m^-2`, respectively 0.494% below and
1.123% above Turner. These remain diagnostic global-observable comparisons,
not new acceptance gates.

| Third-block artifact | SHA-256 |
| --- | --- |
| Input deck | `c57fc271e8135469aef68a4990494578547e2443c081cdd1f91e23eb471f0d2c` |
| Density profile | `903cc6b8b608e8cb8f0bbd72d60776c4e171b767839ec248cebe71b9c0cf0145` |
| Averaging metadata | `5e42146a09e6845b263cd19d5a78a6b8ca1bf90512125f7cbdd4ceaa4dbf0bf0` |
| Final step-563,600 checkpoint | `bc018fc6ee68db4d8d96908e04454142622a21f6cbfd0b52ad6b6993952d1331` |
| Density comparison | `df7a898799b5760f95b7ddd8a318065f39063413306a63bf03eb359e305c0810` |
| Three-block audit | `4e60cd3ce3ce50109d528b1f2ed19d02c714bda7c8dc1757c4b7b3dd1522289a` |
| Source/wall/power balance | `92ce097d68dcae2121a367aa307fbf54ef509d7d59ac1ca84ff7fc40f41ba67e` |

```sh
python3 scripts/analyze_turner_density_blocks.py \
  diagnostic-block-1-comparison.json \
  diagnostic-block-2-comparison.json \
  --output density-block-analysis.json
```

Five further low-priority serial continuations completed the declared
eight-block floor through cycle 1,569 (step 627,600). One attempt at the final
block was externally interrupted after 400 steps; its partial output was
excluded and the block was rerun cleanly from the intact step-614,800
checkpoint. Every new complete block has exact integer source/loss closure.

| 32-cycle block | X² | Relative L2 | Integrated-density change | Gate |
| --- | ---: | ---: | ---: | --- |
| 1,314--1,345 | 365.343 | 2.362% | baseline | diagnostic only |
| 1,346--1,377 | 849.284 | 3.478% | +1.261% | diagnostic only |
| 1,378--1,409 | 484.785 | 3.335% | -0.526% | diagnostic only |
| 1,410--1,441 | 423.871 | 2.796% | -0.225% | diagnostic only |
| 1,442--1,473 | 482.358 | 3.170% | +0.123% | diagnostic only |
| 1,474--1,505 | 257.454 | 1.934% | -1.037% | diagnostic only |
| 1,506--1,537 | 283.680 | 1.926% | +0.308% | diagnostic only |
| 1,538--1,569 | 198.202 | 1.683% | -0.366% | diagnostic only |

The complete series is now eligible for interpretation, but does not
demonstrate stationarity. The line-integrated density has lag-one correlation
`0.261`, leaving an AR(1) effective count of only `4.69` blocks, and its fitted
change across the series is -1.124%. Adjacent profile movement reaches 1.908%.
Individual windows span `X² = 198--849`; three happen to lie inside the
published 95% interval, demonstrating why selecting one favorable
post-benchmark window would be misleading.

The final block again supports the global-physics result: its ion current,
electron power, and ion power differ from Turner by +0.247%, +0.355%, and
+0.727%, respectively. The committed nonrestricted summary evidence is
[`benchmarks/ccp/turner-case1-density-blocks-8-20260731.json`](../benchmarks/ccp/turner-case1-density-blocks-8-20260731.json)
(SHA-256
`23001a2086cae8df032bcd970963e0c34fbb2528646eba4eec79cdf0b33a1d7c`).
The original published-duration `X² = 574.399` failure remains unchanged.

Before observing block 9, the longer stationarity horizon and stopping rule
were predeclared in
[`benchmarks/ccp/turner-case1-stationarity-rule-20260731.json`](../benchmarks/ccp/turner-case1-stationarity-rule-20260731.json)
(SHA-256
`75591319e7e90cc4cbcb35409b348e6b4ef6e413e293a284d5e64df2d7757fdd`).
The horizon is 16 consecutive blocks. All of these internal gates must pass:

- at least 16 total blocks and 8 AR(1)-effective blocks;
- no more than 1% absolute fitted integrated-density drift;
- no more than 1% absolute first-half/second-half integrated-density shift;
- no adjacent density-profile movement above 2.5% relative L2.

The 1% practical-equivalence scale follows the 1.03% median pointwise
ion-density population scatter in the locked Turner Case 1 reference. These
are AuroraPIC diagnostic gates, not thresholds published by Turner et al.
Passing them would establish a sufficiently stationary continuation for
interpretation; it would not reclassify the original published-duration
failure. If the 16-block screen passes, the stronger next test is an
independent-seed published-duration ensemble, which costs roughly three full
1,280-cycle trajectories.

Block 9 (cycles 1,570--1,601, steps 627,601--640,400) was collected only
after that rule was committed. Its `X² = 513.199`, relative profile L2 error
is 3.190%, and line-integrated density rises 1.205% from block 8. Exact
source/loss balance closes; ion current, electron power, and ion power differ
from Turner by +0.649%, +0.283%, and +0.993%. Across nine blocks, absolute
fitted drift is 0.655%, the split-half shift is 0.692%, maximum adjacent
profile movement is 2.030%, and the AR(1) effective count is 7.11. Thus the
three numerical-equivalence gates currently pass, while the predeclared
16-block horizon and 8-effective-block gates remain incomplete.

| Block 9 artifact | SHA-256 |
| --- | --- |
| Input deck | `68c828c639b72ccee51ef253008e70ce08b6ef4b1bbda831e18a6afd369bf0ad` |
| Density profile | `d4f5d067485e0a98ae2cf4309b3d4fdc6a464d696b436eca0007a0b91c0e3226` |
| Averaging metadata | `2ebdd28be609bbbcb17df04482bb16f8a08383c230ef16aace368385cdeb318a` |
| Final step-640,400 checkpoint | `042b2ab4d88ec7a675f803bfc0ec77a4322a77029d681f34eb175f554439e688` |
| Density comparison | `a24dc3517d018ee58b525d57bd222023b3fe7da4ea1209cde306d2f1131509fc` |
| Source/wall/power balance | `e119c8f8a581b6c238446ee4a2f735bd68e25dfd1a8745efce1832af03323077` |
| Nine-block analysis | `ac991bbb48605c746e07a3f66a409a1220ac301b7570af87db199577620eadd6` |

Block 10 (cycles 1,602--1,633, steps 640,401--653,200) gives
`X² = 329.914`, 2.517% relative profile L2 error, and a -0.572%
line-integrated-density change. Its exact source/loss balance closes and its
global-observable differences remain within 1.09%. Across ten blocks,
absolute fitted drift is 0.618%, the split-half shift is 0.649%, maximum
adjacent profile movement remains 2.030%, and the AR(1) effective count rises
to 8.29. All numerical stationarity gates now pass; only the predeclared
16-block horizon is incomplete.

| Block 10 artifact | SHA-256 |
| --- | --- |
| Input deck | `6dfe5911b37542646cd6e1bc6a3e33cf423284dc46d9493c54096b941f772089` |
| Density profile | `b70a3ab23f591b35ea36b4c02a6d9e2ee039b77e17664915bc67dc5e9a9d506f` |
| Averaging metadata | `1d825788588fd127415535a4d78d82b5d34d25a240c074acd28762879459d72c` |
| Final step-653,200 checkpoint | `87501112f4ab336ba2e32c32d86543055aaf728c75b18d0a294b99ec36f84dc2` |
| Density comparison | `49a6cca9fbadcd9c9b655bc975342d548bb130ae8f9d92ed77b04da007260e18` |
| Source/wall/power balance | `18a4ba5fd2d01aaf2f7a198de7d1d4d2114142f43b47e070262090c04f01d4e8` |
| Ten-block analysis | `e3a4c9a26256ccce3fb9c0a76e535a952272cea6d596ab70be91a7f54f3bd2af` |

Block 11 (cycles 1,634--1,665, steps 653,201--666,000) gives
`X² = 411.450`, 3.149% relative profile L2 error, and a +0.465%
line-integrated-density change. Exact source/loss balance closes and all
three global-observable differences remain below 0.99%. Across eleven blocks,
absolute fitted drift is 0.373%, the split-half shift is 0.443%, maximum
adjacent profile movement remains 2.030%, and the AR(1) effective count is
9.53. All numerical gates remain satisfied; the 16-block horizon is the only
incomplete gate.

| Block 11 artifact | SHA-256 |
| --- | --- |
| Input deck | `4997f0d76f4d96204fed94b6be0defffdcbede3725488d2d26c2db65ef7deab1` |
| Density profile | `2e1e77079db46c5115e5152cf1c8b143d567e99165ec43c63d6249ca873a4f5f` |
| Averaging metadata | `0015cad03402f5738ad21f8d648e32346a62e181bee20b6e3a29a3edfe81112c` |
| Final step-666,000 checkpoint | `517c1779c6d4d9f2a5db23e09c97af554a79a870969383b48dfec73190f00882` |
| Density comparison | `efcefb05bf42816c55b87817967e008e90ceb206050e7e8635ddf6d4d2987965` |
| Source/wall/power balance | `84bcf8574f3ced778b6ea55d935bf12e640dc86a4f5b712c23819d085abb3adc` |
| Eleven-block analysis | `8e6785fb3801b4f711f04829533921eb9bcea492e3ca23c2e79023c4987ff669` |

Block 12 (cycles 1,666--1,697, steps 666,001--678,800) gives
`X² = 241.053`, 2.296% relative profile L2 error, and a -0.706%
line-integrated-density change. Exact source/loss balance closes and all
global-observable differences remain below 0.76%. Across twelve blocks,
absolute fitted drift is 0.495%, the split-half shift is 0.316%, maximum
adjacent profile movement remains 2.030%, and the AR(1) effective count is
11.26. All numerical gates remain satisfied; four blocks remain in the
predeclared horizon.

| Block 12 artifact | SHA-256 |
| --- | --- |
| Input deck | `f86d70d0b0d0a3fdd08f20f639acb753286e94078983087ef336765084a96c1a` |
| Density profile | `cea539da10f3a5ee024fd86049c14cd0bd7f54a3f6ab909bb3543659bf778588` |
| Averaging metadata | `1fe14ed70bc0091649d7d5e523bf65e8c13abfe252e0ec41252816dc59a954ea` |
| Final step-678,800 checkpoint | `e507d5f161b15ac62c01f3423db0fca1fe18275c6425f3e00009be414826d1cd` |
| Density comparison | `ce504d5fbe51b770b19117f10689ee38f1bf8187e89c7d269b03b3bfa7331778` |
| Source/wall/power balance | `c02c3e73eaf59fd6f67dfec6c6948bce79eff7bd5bc21eb94e219dfc86ca2cb6` |
| Twelve-block analysis | `c8086be567055b5f4aefc17b707b4e8440947253a5ae0f0afbc6795c648d387c` |

Block 13 (cycles 1,698--1,729, steps 678,801--691,600) gives
`X² = 249.082`, 1.866% relative profile L2 error, and a -0.243%
line-integrated-density change. Exact source/loss balance closes and all
global-observable differences remain below 1.15%. Across thirteen blocks,
absolute fitted drift is 0.674%, the split-half shift is 0.355%, maximum
adjacent profile movement remains 2.030%, and the AR(1) effective count is
10.58. All numerical gates remain satisfied; three blocks remain in the
predeclared horizon.

| Block 13 artifact | SHA-256 |
| --- | --- |
| Input deck | `5cc50612ba7e2bfa8a9f96e9ec165473bcadd2c2a1431a72f3edc9c828131f71` |
| Density profile | `b8f8033cae4a03263886cac043ad1649cc5d8a1db6299e769fd1254dec38ab42` |
| Averaging metadata | `7ad30653507bc6f6c3d000edadf45d04b4611dae09608ad8454b235b5ff5e41a` |
| Final step-691,600 checkpoint | `0d4c7e9255616b20cc865ef106ce7dcff5b27e72a57fc59f5ee80aa9a5274cd4` |
| Density comparison | `6229d38d94adb784d3561645be86b47d8cafbf4553cc17fccb0aba1df94d1c4e` |
| Source/wall/power balance | `3ad725933d4265278781fb6fcfd41b64a6b444efbeabe92860addcc7eb7b9e5e` |
| Thirteen-block analysis | `19278776bc1e951dea28e82f11c6025bae9f3759521588e10acaf953255f2718` |

Block 14 (cycles 1,730--1,761, steps 691,601--704,400) gives
`X² = 159.130`, 1.621% relative profile L2 error, and a -0.341%
line-integrated-density change. Exact source/loss balance closes and all
global-observable differences remain below 0.54%. Across fourteen blocks,
absolute fitted drift is 0.926%, the split-half shift is 0.385%, maximum
adjacent profile movement remains 2.030%, and the AR(1) effective count is
9.08. All numerical gates still pass, but the drift is close to its 1% limit
and two blocks remain in the predeclared horizon.

| Block 14 artifact | SHA-256 |
| --- | --- |
| Input deck | `c1dfdbdbf6d009eb37584c8858547a5c443848bd6bbbea266423d8e85349e449` |
| Density profile | `8c5138efd21bf12e58fc563760c246927debd7f6b03b381f64f66e00f010e364` |
| Averaging metadata | `9fc6b11655a67d8d8d655b640ee25a23d565028573f00f788b1cc150b0e37906` |
| Final step-704,400 checkpoint | `f802a08fa2963732f2d65f60e3fa137f6b361965167feb26f2866cbe04163d17` |
| Density comparison | `6aaab56126232e8306a0df22498bb1fed239d066ab5a8d4bb6b6d3fda58b2e9c` |
| Source/wall/power balance | `1bc86c3253e33877c55a9a459bb4f12a50eccf1d58fdbae2a7a6d646f3693a92` |
| Fourteen-block analysis | `bf5233ef79905d03016721b10ef1f04bb9d1d9e9d4a11d1f787d3d241bc1ddb8` |

Block 15 (cycles 1,762--1,793, steps 704,401--717,200) gives
`X² = 259.494`, 1.847% relative profile L2 error, and a +0.440%
line-integrated-density change. Exact source/loss balance closes and all
global-observable differences remain below 0.87%. Across fifteen blocks,
absolute fitted drift is 0.957%, the split-half shift is 0.352%, maximum
adjacent profile movement remains 2.030%, and the AR(1) effective count is
8.30. The drift and effective-count gates remain narrowly satisfied; one
block remains in the fixed horizon.

| Block 15 artifact | SHA-256 |
| --- | --- |
| Input deck | `41309d525348c3b66fa81276835e33bce0132359215fc1b8192d863ffe0aaefa` |
| Density profile | `63e989cfd9c9f70f173abb36e12acbf502c108d6b93759c74746e1cf1c2611cf` |
| Averaging metadata | `31d5cb28789d7a480a51cbc7db9dd5d655d5ec00cd4fbdcd1fa1fcab52b4a7cb` |
| Final step-717,200 checkpoint | `9ec744e5feca4745497b369250cd2503f4c74aef5e8f8468de5c6f451fe70f29` |
| Density comparison | `3fa6c31b546d2f68a0e1ffa9a043686564cfb8b6d581db749432bdfb7c121938` |
| Source/wall/power balance | `ecd038908a0efd79d68f7b1ea358ae5b09db2b79f16abd8d7dccef36ee2efc1d` |
| Fifteen-block analysis | `d5c08866f0a6345bf0b1053955faa9cd71d33ac8c13893f7441f4d037ee480da` |

### Predeclared 16-block stationarity result

Block 16 (cycles 1,794--1,825, steps 717,201--730,000) gives
`X² = 799.512`, 3.388% relative profile L2 error, and a +1.486%
line-integrated-density change. Exact source/loss balance closes. Its global
ion-current, electron-power, and ion-power differences are +1.537%, +1.479%,
and +1.833%; these remain close but are the largest of the recent windows and
have no single-window published acceptance gates.

Despite the noisy final window, the predeclared full-series stationarity
screen **passes all five gates**:

| Gate | Threshold | 16-block value | Result |
| --- | ---: | ---: | --- |
| Total blocks | at least 16 | 16 | pass |
| AR(1) effective blocks | at least 8 | 11.84 | pass |
| Absolute fitted density drift | at most 1% | 0.484% | pass |
| Absolute split-half density shift | at most 1% | 0.094% | pass |
| Maximum adjacent profile L2 | at most 2.5% | 2.030% | pass |

This establishes that the late continuation is stationary at the resolution
of the predeclared AuroraPIC diagnostic. It also quantifies substantial
32-cycle sampling variability: `X²` spans 159--849, only 7 of 16 windows lie
inside the published 95% numerical interval, and 9 of 16 lie inside the 99%
interval. Therefore the result does not retroactively make the original
published-duration `X² = 574.399` comparison pass. It instead supports the
interpretation that the approximately 2--3% density discrepancy is comparable
to correlated finite-window noise on a single trajectory.

The committed, nonrestricted result is
[`benchmarks/ccp/turner-case1-density-blocks-16-20260731.json`](../benchmarks/ccp/turner-case1-density-blocks-16-20260731.json)
(SHA-256
`f9b92a2b0232cd192910254a998f8030422266f4c3fe9f0f3c1f7e8ed7a0846b`).
The next credibility test is an independent-seed, published-duration ensemble;
that is the appropriate way to distinguish seed variability from a residual
systematic density bias.

## Independent-seed production ensemble

`scripts/prepare_turner_case.py` now accepts an audited unsigned-32-bit
`--seed` override and records it in both the generated deck and checksum-bound
preflight contract. `scripts/prepare_turner_ensemble.py` builds on that
contract to atomically prepare 3--16 unique full-duration seed decks without
launching them. It refuses overwrite, requires a separate aggregate cost
acknowledgement, and fixes workstation concurrency at one run.

```sh
python3 scripts/prepare_turner_ensemble.py \
  examples/turner_helium_ccp_case1.case \
  tmp/turner-normalized-v3 \
  --output-dir /external/campaign/turner-case1-ensemble \
  --seeds 13507,24680,97531 \
  --acknowledge-cost \
    I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_ENSEMBLE
```

The ensemble manifest records every deck and preflight SHA-256, result and
comparison locations, the aggregate initial particle-update floor, and an
explicit `launched = false` claim boundary. Preparation does not authorize
parallel launch or establish any physics result. The existing completed
seed-13507 trajectory may be attached later only through a checksum-verified
result-ingestion step; the preparer does not silently treat it as an ensemble
member.

The first real campaign preparation uses seeds 13,507, 24,680, and 97,531.
All three full 512,000-step decks pass `aurorapic_cli --validate-only`; none
was launched. The aggregate preflight floor is 201,326,592,000 initial
particle updates and 1,572,864 aggregate capacity slots, with only one run
authorized concurrently. Based on the completed local seed-13,507 trajectory,
planning should allow roughly 2.2 hours and about 0.5 GiB per run on this host;
these are observations, not portable performance guarantees.

| Prepared ensemble artifact | SHA-256 |
| --- | --- |
| Ensemble manifest | `5ebf633d827403c8fab97c7fabbb04505249bc90fde1dc412f57c929ec338358` |
| Seed-13,507 deck | `1cb7abc71c6e38aec2c0cd51dfe77d22158889923463ed06d45a587a46ba4a0f` |
| Seed-13,507 preflight | `3e289b7efc2cd67fbfd56acf0d4b95cd5a01a92a009192b3df62c67098e393b4` |
| Seed-24,680 deck | `552cb91fe121fc408f72b6c328722153c91e0de96bfaf5f31e549b5d1c012838` |
| Seed-24,680 preflight | `1b4d4f4d9223ae452e5edf91f793ec6ee1d56fa2fca695e59c5cdca84f60d47e` |
| Seed-97,531 deck | `0153288faf1ebb0baaee0be3a4d6e23395175cc22987cf39dbb7d9f518f83f53` |
| Seed-97,531 preflight | `9af247d3e97de6bcc071f505908175a0f558e283406b60f6bd23c9f8c06a3100` |

That attachment gate is now implemented by
`scripts/attach_turner_ensemble_result.py`. It leaves the immutable preparation
manifest unchanged, permits the executed deck to differ only in `output_dir`,
recomputes the complete published-duration comparison from the locked local
reference, candidate profile, and averaging metadata, and checksum-records the
final checkpoint. A changed waveform or any other physics/numerical deck field
is rejected in regression.

All three prepared runs are now complete and attached. The independently
recomputed individual results are:

| Seed | `X²` | 95% result | 99% result | Integrated ion-density bias |
| ---: | ---: | :---: | :---: | ---: |
| 13,507 | 574.399 | fail | fail | +2.483% |
| 24,680 | 320.601 | fail | pass | +1.730% |
| 97,531 | 604.587 | fail | fail | +2.497% |

`scripts/analyze_turner_ensemble.py` requires exactly one attachment for every
prepared seed, revalidates the manifest, executed configs, checkpoints,
profiles, averaging metadata, normalization audit, and stored comparisons,
then recomputes every comparison before aggregating. It deliberately declares
`formal_ensemble_acceptance_rule = none_predeclared`: the published intervals
apply to individual realizations, and no aggregate threshold was selected
before observing these runs.

The three-seed mean integrated density bias is +2.237%, with a 0.439%
between-seed sample standard deviation; all three biases have the same sign.
Only one member passes the published 99% interval and none passes the 95%
interval. Under an explicitly diagnostic assumption of independent 1% failure
probability, observing at least two 99% failures in three trials has probability
`0.000298`. This is not a new formal pass/fail test, but it makes a purely
chance explanation implausible and elevates the positive density-amplitude
offset to a systematic discrepancy requiring isolation.

The nonrestricted aggregate is
[`benchmarks/ccp/turner-case1-seed-ensemble-3-20260801.json`](../benchmarks/ccp/turner-case1-seed-ensemble-3-20260801.json).
The complete local analysis report SHA-256 is
`423b420efdd10ed4bb899ab82ad173a4f94a290e23f2c129d6a85367f49b9254`.
The next credibility action is therefore a predeclared physics sensitivity
matrix, not additional unstructured production seeds: isolate collision-table
interpretation, scattering/ionization choices, and numerical resolution while
tracking density amplitude, current, and species power together.

## Predeclared discrepancy-isolation matrix

Numerical convergence is tested before changing any prescribed collision
physics. `scripts/prepare_turner_sensitivity.py` consumes one checksum-verified
ensemble deck and creates four same-seed, full-1,280-cycle diagnostic variants
without launching them:

| Stage | Variant | Change | Question |
| ---: | --- | --- | --- |
| 1 | `particles_2x` | 1,024 particles/cell/species | Is the density excess driven by macro-particle noise or discrete heating? |
| 1 | `timestep_2x` | 800 steps/RF cycle | Is the excess sensitive to pushing, RF integration, or MCC time splitting? |
| 2 | `grid_2x_fixed_particles` | 256 cells, fixed total particles | Isolate field-grid resolution from total particle count. |
| 2 | `grid_2x_same_ppc` | 256 cells and 512 particles/cell/species | Check joint grid/particle convergence. |

All variants retain the physical duration and final-32-cycle average, but they
change the published numerical contract. Their published `X²` values may be
reported descriptively and must not be classified as Turner benchmark passes.
`compare_turner.py --numerical-sensitivity` enforces this distinction. It
derives an integer RF timestep-refinement ratio from the averaging metadata,
accepts only integer-refined grids (mapping coincident nodes without
interpolation), and always emits
`published_acceptance_applicable = false`.
Only one low-priority serial run is permitted at a time; stage 2 is deferred
until stage 1 is interpreted.

The primary sensitivity observable is the change in integrated ion-density
bias relative to the paired seed-13,507 baseline (`+2.483%`). Before running a
variant, an absolute shift of at least 0.75 percentage points is declared
material, at most 0.50 percentage points is practical equivalence, and the
interval between them is ambiguous. The 0.75-point threshold is approximately
three standard errors of the completed three-seed baseline mean, rounded
conservatively. Current and species power must be tracked as corroborating
observables; no variant will be selected merely because it produces a favorable
`X²`.

The stage-1 `particles_2x` run completed all 512,000 steps with 1,024
particles/cell/species. Its integrated ion-density bias is `+2.409%`, versus
`+2.483%` for the paired seed-13,507 baseline: a shift of only `-0.075`
percentage points. This is inside the predeclared 0.50-point practical
equivalence boundary. The relative profile L2 error is essentially unchanged
(`2.865%` versus `2.874%`). Its descriptive `X² = 428.312` is outside the
published 99% interval, but cannot be classified as a published Turner result
because the particle count changed.

The result argues against ordinary macro-particle count noise or discrete
particle heating as the primary source of the approximately 2.2% ensemble
density excess. It does not prove infinite-particle convergence from one
refinement ratio. The retained nonrestricted evidence is
[`benchmarks/ccp/turner-case1-particles-2x-sensitivity-20260802.json`](../benchmarks/ccp/turner-case1-particles-2x-sensitivity-20260802.json).
The next stage-1 probe is `timestep_2x`; grid variants remain deferred.

The `timestep_2x` run then completed 1,024,000 steps at 800 steps/RF cycle.
Its integrated ion-density bias is `+3.117%`, a `+0.634` percentage-point
shift from the paired baseline. This lies in the predeclared ambiguous interval
and moves away from, rather than toward, the reference. Its relative profile
L2 error is `3.725%`; the descriptive `X²` is `784.852`. Thus straightforward
timestep error does not explain the positive density offset, but one refinement
ratio does not establish an asymptotic temporal trend.

The retained nonrestricted evidence is
[`benchmarks/ccp/turner-case1-timestep-2x-sensitivity-20260803.json`](../benchmarks/ccp/turner-case1-timestep-2x-sensitivity-20260803.json).
Because stage 1 did not isolate the discrepancy, the predeclared next run is
`grid_2x_fixed_particles`.

The `grid_2x_fixed_particles` run completed 512,000 steps on 256 cells while
retaining 65,536 particles/species. Its integrated ion-density bias is
`+2.018%`, a `-0.465` percentage-point shift from the paired baseline and just
inside the predeclared practical-equivalence boundary. The coincident-node
relative profile L2 error improves to `2.480%`; its descriptive `X²` is
`426.639`. The result suggests a modest spatial-resolution contribution, but
does not classify the density excess as materially grid-sensitive.

Holding total particles fixed halves particles/cell from 512 to 256, coupling
grid refinement to particle statistics. The retained nonrestricted evidence is
[`benchmarks/ccp/turner-case1-grid-2x-fixed-particles-sensitivity-20260803.json`](../benchmarks/ccp/turner-case1-grid-2x-fixed-particles-sensitivity-20260803.json).
The final predeclared variant, `grid_2x_same_ppc`, separates that effect by
restoring 512 particles/cell on the refined grid.

That final variant completed with integrated ion-density bias `+2.334%`, only
`-0.149` percentage points from the paired baseline. Its coincident-node
relative L2 error is `2.687%`, and its descriptive `X²` is `429.713`. It is
therefore practically equivalent under the predeclared rule. The retained
nonrestricted evidence is
[`benchmarks/ccp/turner-case1-grid-2x-same-ppc-sensitivity-20260803.json`](../benchmarks/ccp/turner-case1-grid-2x-same-ppc-sensitivity-20260803.json).

The numerical sensitivity matrix is now complete. None of the particle,
timestep, fixed-particle grid, or fixed-particles-per-cell grid refinements
materially reduces the density bias. All four final windows have exact source
and wall-loss balance. Their represented ionization rates span only -0.270% to
+0.235% about the four-variant mean; ion current differs from Turner by
+0.803% to +1.363%, electron power by +0.455% to +1.025%, and ion power by
+1.128% to +1.708%. This combination argues against ordinary resolution or a
gross ionization/heating error as the primary explanation. It redirects the
next diagnostic toward spatial transport, sheath structure, collision-model
conventions, and independent cross-code comparison.

The complete nonrestricted matrix summary is
[`benchmarks/ccp/turner-case1-numerical-sensitivity-matrix-20260803.json`](../benchmarks/ccp/turner-case1-numerical-sensitivity-matrix-20260803.json).

### Published electron-density diagnostic

The publisher supplement also contains electron-density mean and population
scatter profiles, although Turner et al. give formal `X²` acceptance ranges for
the ion-density statistic. `compare_turner.py --species electrons` therefore
uses the same locked coordinate and uncertainty calculation but always reports
`published_acceptance_applicable = false`.

Across the three exact seeds, integrated electron-density bias is +3.103% on
average, versus +2.237% for ions. The electron bias exceeds the ion bias in
every exact and sensitivity run, by 0.623--1.359 percentage points. The same
seed ordering appears in both species, while the four sensitivity runs retain
closely consistent ionization, current, and power. This descriptive evidence
localizes the remaining discrepancy more toward electron/sheath transport and
density-profile formation than toward a gross total ionization-source error.
It is not a formal electron-density benchmark pass/fail conclusion.

The retained nonrestricted evidence is
[`benchmarks/ccp/turner-case1-electron-density-diagnostic-20260803.json`](../benchmarks/ccp/turner-case1-electron-density-diagnostic-20260803.json).

### Electron-energy and sheath-structure localization

Checkpoint v8 now carries the final-window species kinetic-energy density,
mean potential, mean electric field, and squared electric field. A 32-cycle
reset window at steps 730001--742800 continued the stationary seed-13507
trajectory in serial mode. Source/loss balance remained exact. Ion current was
`+1.302%` above Turner, while electron and ion electrical power were `+1.180%`
and `+1.492%` above the published Table III values.

The new profiles are strongly left-right consistent: electron energy differs
from its mirror by `0.512%` relative L2, RMS field by `0.169%`, and the mean
field satisfies antisymmetry to `0.217%`. The outer 20% regions contain
`85.35%` of integrated RMS-field squared. Their density-weighted effective
electron temperatures are `11.290` and `11.325 eV`, versus `9.528 eV` over the
central half. At the mid-plane the effective temperature is `8.835 eV`,
`-5.609%` from Turner's `9.36 eV`; mid-plane ion density is `+4.439%` high.

This narrows, but does not close, the credibility gap. Symmetry, exact balance,
edge-localized field, and near-reference global power argue against a gross
solver asymmetry or global heating failure. The post-benchmark ion profile
still has `X² = 603.861`, `3.542%` relative L2 error, and `+2.887%` integrated
bias. Published acceptance is inapplicable outside the prescribed duration,
and this result would lie above its 99% range in any case. The next localization
target is phase-resolved, drift-separated electron moments and spatially
resolved collision-channel energy loss, followed by matched cross-code output.

Run the deterministic analysis with:

```sh
python3 scripts/analyze_turner_spatial_structure.py \
  --density diagnostic-output/spatial_average.csv \
  --kinetic-energy diagnostic-output/spatial_kinetic_energy.csv \
  --field diagnostic-output/spatial_field_average.csv \
  --metadata diagnostic-output/spatial_average_metadata.json \
  --output spatial-structure.json
```

The retained checksum-bearing evidence is
[`benchmarks/ccp/turner-case1-spatial-structure-20260803.json`](../benchmarks/ccp/turner-case1-spatial-structure-20260803.json).

A following 16-bin phase-resolved window completed at step 755600 with 800
samples in every bin. Mid-plane drift-separated electron temperature spans
`8.016--9.296 eV`; coherent mean motion contributes only `1.96%` on average
and at most `3.34%` of mean kinetic energy. Thus RF bulk drift does not explain
the approximately `6%` low cycle-averaged mid-plane temperature. Half-cycle
temperature symmetry is `0.480%` relative L2, peak-field symmetry is `0.136%`,
and the applied voltage satisfies half-cycle antisymmetry to `1.6e-10`
relative L2. The field maximum correctly transfers between electrodes as the
voltage polarity reverses.

The density discrepancy persists: the later block has ion `X² = 781.004` and
`+3.281%` integrated ion bias. Global power/current agreement and exact balance
remain intact. This moves the next diagnostic toward spatial collision-channel
energy transfer and phase-selected velocity distributions rather than coherent
RF drift. The retained evidence is
[`benchmarks/ccp/turner-case1-phase-structure-20260804.json`](../benchmarks/ccp/turner-case1-phase-structure-20260804.json).

```sh
python3 scripts/prepare_turner_sensitivity.py \
  tmp/turner-case1-ensemble-v1/ensemble.json \
  --baseline-seed 13507 \
  --baseline-density-bias-percent 2.4834268915580937 \
  --executable build/aurorapic_cli \
  --output-dir /external/campaign/turner-case1-sensitivity-v1 \
  --acknowledge-cost \
    I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_SENSITIVITY
```

Generate a checksum-bearing balance report for any fully covered SI diagnostic
window with:

```sh
python3 scripts/analyze_turner_balance.py \
  --scalars diagnostic-output/scalars.csv \
  --collisions diagnostic-output/collisions.csv \
  --boundary-losses diagnostic-output/boundary_losses.csv \
  --power-transfer diagnostic-output/power_transfer.csv \
  --expected-steps 12800 \
  --output diagnostic-balance.json
```

The analyzer rejects mismatched step windows, partial wall-counter coverage,
partial power-counter coverage, decreasing source/wall counters, non-finite
values, and nonzero integer source/loss residuals. Its reports retain
`physics_claim = none` because a post-benchmark continuation is diagnostic
evidence rather than the published test.

For a final window embedded in a full numerical-sensitivity output, add exact
`--window-start-step` and `--window-end-step` endpoints plus
`--scope published_duration_numerical_sensitivity_window`. The analyzer
selects matching rows from every diagnostic, requires all endpoints and counter
coverage, and records that the published numerical contract changed.

## eduPIC normalization audit

The twelve-cycle matched-heating comparison left a stable `+11.898%`
AuroraPIC electron-density difference and a `-13.560%` power-per-electron
difference. Before launching more kinetic trajectories, the complete 1D
normalization path was audited against the pinned eduPIC C implementation and
the published 400-by-200 phase matrices.

The two codes are normalization-equivalent. eduPIC's superparticle weight of
`7.0e4` over its fictive `1.0e-4 m2` electrode area is exactly AuroraPIC's
`7.0e8 m-2` line weight. Both use linear CIC deposition. eduPIC's factor of two
at either boundary is the same half-cell control volume that AuroraPIC uses
directly. Both phase diagnostics are arithmetic means, both use pre-collision
moments with leapfrog-centered longitudinal velocity, and both form electron
current as `-e n u`.

The archived result provides a direct conservation check. Trapezoidal spatial
integration of every published density phase implies a mean of `106463.3074`
eduPIC electron macroparticles; the twelve-cycle AuroraPIC phases imply
`119130.4915`. Their ratio, `1.118981688364513`, equals the independently
reported physical-density ratio to `2.22e-16` relative error. Replacing the
physical half-cell quadrature with a simple 400-node average changes that
cross-code ratio by only `1.54e-7`. The published eduPIC heating matrix also
reconstructs from its current and electric-field matrices to `2.30e-7`
relative L2, consistent with text-output rounding.

This rules out macro weight, electrode area, endpoint volume, phase averaging,
current units, and spatial quadrature as explanations for the density and
heating discrepancies. It does not show that the kinetic states are equal.
The next discriminator is an independent-initial-realization ensemble from
the same macroscopic state.

The checksum-bearing evidence is
[`benchmarks/ccp/edupic-argon-normalization-audit-20260819.json`](../benchmarks/ccp/edupic-argon-normalization-audit-20260819.json).
The audit can be reproduced locally with:

```sh
PYTHONPATH=scripts python3 scripts/audit_edupic_normalization.py \
  candidate-run-root reference-raw-data eduPIC.cc \
  --output normalization-audit.json
```

### Constrained initial-microstate ensemble

The next test is predeclared in
[`benchmarks/ccp/edupic-argon-microstate-ensemble-rule-20260820.json`](../benchmarks/ccp/edupic-argon-microstate-ensemble-rule-20260820.json).
Two deterministic conditional randomizations preserve the source state much
more tightly than an unconstrained bootstrap: cell occupancy and cellwise
velocity-tuple multisets are exact, while each cell's summed CIC fraction—and
therefore its nodal number deposition—is preserved to roundoff. Subcell
positions are redrawn under pair-sum constraints and velocity tuples are
independently re-paired with positions inside each cell.

For seeds 51,949 and 63,059, more than 99.8% of particle positions change,
with approximately `18 micrometres` RMS displacement. Species nodal-number
relative L1 errors are between `9.3e-16` and `1.1e-15`. These are controlled
independent microstates conditional on the empirical cell populations and
velocity samples; they are not independent draws from an unknown continuous
distribution. The two branches will use the same collision RNG seed, discard
two RF cycles, and measure four fresh cycles with the established 200-bin
pre-collision protocol. Runs remain serial with a 256 MiB RSS ceiling.

Both branches subsequently completed and passed every solver, diagnostic, and
resource gate. Across the locked source member and the two constrained
microstates, the volume-density range is only `0.0859%`, the power range is
`2.75%`, and maximum member-to-ensemble cycle-average spatial-current scatter
is `4.01%`; all pass their prospective `3%`, `8%`, and `8%` limits.

The ensemble does not close the cross-code differences. Mean density remains
`3.3396e15 m-3`, corresponding to an `11.87%` spatial-profile L2 mismatch.
Mean power per electron is `0.8643` of eduPIC, a `13.57%` deficit. The
cycle-average spatial-current mismatch is `36.95%` L2 with `0.9955`
correlation, versus at most `4.01%` member scatter. Density-profile scatter is
at most `0.232%`. The electric-field phase-space mismatch is only `3.05%`.
Thus the persistent density and current-amplitude/localization differences are
not explained by the audited constrained microscopic realization.

The retained result is
[`benchmarks/ccp/edupic-argon-microstate-ensemble-20260820.json`](../benchmarks/ccp/edupic-argon-microstate-ensemble-20260820.json).
The next prospective discriminator is electron current per represented
particle together with ionization and wall-loss balance, separating transport
and heating from population regulation.

### Phase-resolution and particle-balance audit

That discriminator is now recorded in
[`benchmarks/ccp/edupic-argon-transport-balance-audit-20260820.json`](../benchmarks/ccp/edupic-argon-transport-balance-audit-20260820.json).
It is a post-hoc localization audit, not a new acceptance test. The analyzer
validates the three matched branch contracts, verifies the locked continuous
reports, all 64 native eduPIC stage reports, and their relevant diagnostic
hashes, and records direct hashes for the AuroraPIC ledgers. AuroraPIC's
ionization, wall-loss, and live-particle counters close exactly for electrons
and ions.

Conservatively averaging both codes from 200 to 16 RF phase bins changes the
matched-ensemble electron-current phase-space relative L2 error from `17.90%`
to `5.33%`, while correlation rises to `0.99936`. The field error is `3.03%`
at 16 bins. Density is effectively unchanged at `11.87%`. This reproduces the
earlier 16-bin current result and shows that the large 200-bin local-current
error is dominated by the short candidate window's per-bin sampling noise;
it is not evidence by itself of an incorrect current definition or particle
push. Power-density error similarly falls from `20.73%` to `12.80%`.

The particle balance isolates the remaining concern. Across the three
four-cycle matched members, mean electron source-minus-wall-loss is `-9.97%`
of the ionization source. The independently run twelve-cycle window remains
negative at `-7.59%`. In contrast, the checksum-bound 1,024-cycle native
eduPIC campaign has mean ionization source `4.6222e18 m-2 s-1`, electron wall
flux `4.6088e18 m-2 s-1`, and mean imbalance `+0.280%` across its 64 contiguous
blocks. Those blocks are descriptive and are not treated as independent
samples.

The strongest justified conclusion is therefore narrower and more useful:
AuroraPIC reproduces the published electric-field and electron-current RF
waveforms closely at a common noise-reduced phase resolution, while its present
candidate state has an unresolved density/source-loss stationarity error. The
candidate event counts are finite and the constrained members share a
collision seed, so this audit does not assign formal confidence intervals. The
next prospective run should qualify source/loss balance over longer independent
blocks before density agreement is claimed; changing the pusher or current
normalization is not supported by this evidence.

The audit can be reproduced without advancing either simulation:

```sh
PYTHONPATH=scripts python3 scripts/analyze_edupic_transport_balance.py \
  microstate-rule baseline branch-51949 branch-63059 \
  continuous-cycle24 continuous-cycle28 continuous-cycle32 long-window \
  reference-raw-data native-campaign-root --output transport-audit.json
```

The follow-up source/loss campaign is prospectively locked in
[`benchmarks/ccp/edupic-argon-balance-stationarity-rule-20260820.json`](../benchmarks/ccp/edupic-argon-balance-stationarity-rule-20260820.json).
It continues the full cycle-14 checkpoint in one low-priority serial 16-cycle
block per invocation. Before any solver step completed, strict restart checks
showed that the full checkpoint requires compatible spatial and wall-impact
contracts. The locked operational amendment therefore resets spatial sampling
at the ordinary 400-step output interval with ten phase bins and retains the
stored 200-bin wall-impact shape; physical inputs, thresholds, and resource
limits are unchanged. Scalar, collision, boundary-loss, field, resource, and
final-checkpoint evidence remains checksum-bound. A block must close both
species ledgers
exactly, keep each source/loss imbalance within `3%`, keep each normalized
species-population slope within `5e-4` per cycle, and pass field-energy and
ionization-variation gates. Two consecutive passing blocks are required before
a fresh density measurement; eight blocks are targeted before blockwise
uncertainty is assessed. Contiguous blocks will not be assumed independent.

Block 1 subsequently completed 16 cycles safely. Its checksum-bearing report
is
[`benchmarks/ccp/edupic-argon-balance-stationarity-block01-20260820.json`](../benchmarks/ccp/edupic-argon-balance-stationarity-block01-20260820.json).
Both species particle ledgers close exactly, peak RSS is `183956 KiB`, and all
hard-safety gates pass. Four of six stationarity gates pass: electron and ion
population slopes are `-2.55e-4` and `-2.42e-4` per cycle, field-energy slope
is `+2.47e-4` per cycle, and ionization-count CV is `4.86%`. Electron and ion
source-minus-wall-loss remain `-6.80%` and `-6.39%`, failing the prospective
`3%` limits. The imbalance is therefore persistent across a longer window,
even though macroscopic slopes are slow; block 2 remains warranted to test
whether balance is approaching zero or has settled at a biased plateau.

Block 2 and the chained campaign analysis are retained in
[`benchmarks/ccp/edupic-argon-balance-stationarity-block02-20260820.json`](../benchmarks/ccp/edupic-argon-balance-stationarity-block02-20260820.json)
and
[`benchmarks/ccp/edupic-argon-balance-stationarity-campaign-20260820.json`](../benchmarks/ccp/edupic-argon-balance-stationarity-campaign-20260820.json).
Block 2 again passes every hard, population-slope, field-energy, and
ionization-variation gate, but electron and ion balance worsen to `-10.34%`
and `-10.05%`. Across both blocks, `15230` ionizations compete with `16531`
electron and `16477` ion wall losses, producing combined imbalances of
`-8.54%` and `-8.19%`.

The two-block trend localizes the deficit further. From block 1 to block 2,
ionization events per cycle fall `3.10%`, and ionization per live electron per
cycle falls `2.53%`; electron wall losses per cycle change only `+0.109%`.
Thus the quiet candidate state is not merely obscured by short-window noise:
its ionization source is weakening relative to a stable wall-loss sink. More
measurement averaging cannot repair that state. Additional expensive balance
blocks are paused while effective ionization probability per electron, EEDF,
and collision sampling are audited against native eduPIC. The two contiguous
blocks do not identify a unique causal defect or support independent-sample
confidence intervals.

### Ionization-path closure

The requested follow-up is retained in
[`benchmarks/ccp/edupic-argon-ionization-path-audit-20260821.json`](../benchmarks/ccp/edupic-argon-ionization-path-audit-20260821.json).
It checksum-locks the native eduPIC source, AuroraPIC gas manifest and
ionization table, long-window EEDF and counters, the earlier cycle-32
phase-space comparison, and the prospectively gated region-matched collision
audit. Both implementations use the same 0.001 eV energy grid, 15.8 eV
threshold, laboratory collision energy, cold argon target, finite-mass
center-of-mass transform, and 10 eV Opal--Beaty--Peterson ejected-energy law.
eduPIC selects the nearest table bin while AuroraPIC selects the lower bin.

Folding the 12-cycle full-gap AuroraPIC EEDF through the locked
`n_g sigma_ion(E) v(E) dt` kernel predicts `5541.48` ionization events at the
0.25 eV histogram-bin centers. Conservative minimum and maximum kernels within
each histogram bin give `5283.74` to `5810.04`; the exact counter records
`5652`, inside those bounds and only `1.99%` above the center estimate. The
independent cycle-32 window closed at `2.23%`, and the earlier prospective
region audit closed its full-gap ionization rate within `2.70%`.

The direct-collision algorithm difference is also bounded. At the configured
maximum null frequency, `nu_max dt = 0.01844`; the Poisson mean used by
AuroraPIC exceeds eduPIC's one-collision-per-step Bernoulli probability by at
most `0.916%`, and in the direction of slightly more AuroraPIC opportunities.
It cannot explain an `8--10%` source deficit.

The long-window measured effective ionization frequency is `53.61 kHz`, versus
the locked eduPIC value `62.30 kHz`, or a ratio of `0.8605`. This independently
reproduces the cycle-32 ratio `0.8566`; the region-matched EEDF fold gave
`0.8737`. The collision sampler and particle creation ledger are therefore
credible for this case. The deficit is carried by the energetic electron
distribution. The next discriminator is a spatially and phase-resolved audit
of the ionizing EEDF tail and electron power deposition; this result does not
yet identify the unique kinetic cause.

### Spatial and RF-phase ionization localization

The next discriminator is retained in
[`benchmarks/ccp/edupic-argon-spatial-ionization-audit-20260821.json`](../benchmarks/ccp/edupic-argon-spatial-ionization-audit-20260821.json).
It reuses the checksum-bound 12-cycle window and the published 200-phase by
400-node eduPIC matrices; no additional solver run was required. Candidate and
reference ionization rate, electron density, and signed electron `J E` power
are compared without a fitted phase shift or spatial reflection. Whole-gap
integrals use trapezoidal nodal weights, and the non-overlapping spatial and
phase contributions add exactly to the reported net source difference.

The total AuroraPIC ionization rate is `96.287%` of eduPIC while its electron
density is `111.898%`, reproducing the effective-ionization-frequency ratio
`0.86049` from the independent ionization-path audit. The total electron power
ratio is `0.96725`, but power per electron is only `0.86440`. Thus the large
per-electron discrepancy is not an event-counter normalization error: a denser
candidate population nearly hides a substantially less ionizing population.

The source difference is strongly localized. The `0.2--0.4` and `0.4--0.6`
gap bands account for `42.82%` and `51.89%` of the net deficit, respectively;
together they account for `94.71%`. The `0.8--0.9` band instead has `9.62%`
more ionization than the reference and offsets `6.30%` of the net gap. The
electrode-adjacent `0--0.1` and `0.9--1.0` bands together contribute only
`3.14%` of the net gap, even though the final band has a large relative error
where the absolute ionization contribution is tiny.

The RF-phase result is similarly selective: `92.49%` of the net source gap is
accumulated in phase fractions `0--0.5`. The largest octant contributions are
`0.125--0.25` (`34.80%` of the net gap) and `0.25--0.375` (`40.23%`). In the
opposite half-cycle the total source is much closer to the reference. This
combination favors a phase-specific creation and transport difference in the
energetic electron tail that reaches the discharge interior, rather than a
collision sampler defect localized at either electrode.

The signed local power ratios do not map monotonically to local ionization:
for example, the `0.4--0.6` band has more cycle-averaged signed power per
electron but lower ionization per electron. Signed `J E` includes alternating
heating and cooling and does not resolve energy-selective tail production.
The next targeted diagnostic should therefore resolve the above-threshold EEDF
or ionization kernel jointly in space and RF phase, especially over
`x/L=0.2--0.6` and phase `0.125--0.5`, rather than extending the same scalar
balance campaign. These results remain descriptive cross-simulation evidence;
they have no prospective spatial acceptance threshold and are not experimental
validation or a unique diagnosis of the kinetic cause.

### Regional ionizing-tail closure

The follow-up diagnostic was prospectively locked in
[`benchmarks/ccp/edupic-argon-ionizing-tail-rule-20260821.json`](../benchmarks/ccp/edupic-argon-ionizing-tail-rule-20260821.json)
before execution. It continued the checksum-locked cycle-14 state for exactly
four RF cycles on one low-priority CPU. Seven contiguous gap bands, 200 RF
phase bins, and 320 energy bins from 0 to 80 eV were sampled every second
timestep. The run completed in `363.31 s` with `196124 KiB` peak RSS, no energy
overflow, a maximum sampled field of `66.70 kV/m`, and `243436` maximum live
macro-particles.

The execution record is retained in
[`benchmarks/ccp/edupic-argon-ionizing-tail-block-20260822.json`](../benchmarks/ccp/edupic-argon-ionizing-tail-block-20260822.json).
The global prospective observation gate failed and is retained as a failed
gate: the electrode-adjacent 10% bands contain no electrons during 102 of 200
phase bins as the sheaths expand. This is physical depletion rather than a
missing accumulator. All safety, finite-output, shape, sampling, and resource
gates passed. The prospectively named interior `x/L=0.2--0.4` and `0.4--0.6`
bands have at least `1181067` macro-particle observations in every phase bin,
so the declared critical scope is statistically populated despite the honest
global failure.

The checksum-bound result is
[`benchmarks/ccp/edupic-argon-ionizing-tail-audit-20260822.json`](../benchmarks/ccp/edupic-argon-ionizing-tail-audit-20260822.json).
Over the predeclared `x/L=0.2--0.6`, phase `0.125--0.5` window, the measured
AuroraPIC effective ionization frequency is `66.60 kHz`, versus `81.54 kHz`
from the published eduPIC density and ionization matrices: a ratio of
`0.81669`. Independently folding the sampled AuroraPIC EEDF through the locked
argon `n_g sigma_ion(E) v(E)` kernel predicts `67.10 kHz`, or `0.82286` of
eduPIC. Candidate measurement divided by its EEDF prediction is `0.99250`.
Thus the localized EEDF accounts for the cross-code ionization deficit to
within `0.75%` in the combined critical window.

Across the full cycle in the same interior region, the measured candidate to
eduPIC frequency ratio is `0.86710`, the EEDF-folded ratio is `0.85206`, and
the candidate measurement closes against its fold at `1.01765`. The strongest
measured octant deficits are phase `0.25--0.375` (`0.80496` of eduPIC) and
`0.375--0.5` (`0.79584`), independently confirming the phase localization
found from collision-rate matrices.

This is direct evidence that the remaining ionization discrepancy is encoded
in AuroraPIC's energetic electron distribution rather than introduced by its
ionization-event sampler. It still does not identify which kinetic mechanism
creates that distribution difference. The next discriminator must separate
phase-dependent sheath energization, energetic-electron transport into the
bulk, and inelastic cooling. The comparison remains against published
simulation output, not experiment, and no retrospective cross-code acceptance
threshold is assigned.

### Interior electron-energy pathway

The same checksum-bound continuation contains phase-resolved electric power,
tracked collision-energy exchange, density, and mean electron energy. The
post-hoc decomposition is retained in
[`benchmarks/ccp/edupic-argon-energy-pathway-audit-20260822.json`](../benchmarks/ccp/edupic-argon-energy-pathway-audit-20260822.json).
It uses the same predeclared `x/L=0.2--0.6`, phase `0.125--0.5` scope and
introduces no new simulation or fitted alignment.

AuroraPIC's integrated electric power in this window is `56.59 W/m2`, versus
`58.91 W/m2` for eduPIC, a seemingly close ratio of `0.96063`. The candidate
electron column density is `1.12073` times the reference, however, so electric
power per electron is only `0.85714` of eduPIC. The independently measured
effective-ionization ratio is `0.81669`; ionization divided by the power-per-
electron ratio is `0.95281`. Lower energy input per electron therefore explains
most, but not all, of the localized ionization deficit.

The candidate tracked collisional kinetic-power change is `-31.20 W/m2`:
`-20.86 W/m2` from excitation, `-10.23 W/m2` from ionization, and
`-0.106 W/m2` from elastic collisions. The electron kinetic-energy storage rate
is `+2.46 W/m2`. Applying the phase-averaged kinetic-energy equation leaves an
inferred outward energy-flux divergence of `+22.92 W/m2`. This final term is a
residual, not a directly measured flux; the ionization collision channel also
contains the very small newborn-ion kinetic contribution. It cannot yet be
used as an exact species-separated transport measurement.

Seven of eight phase octants have an ionization ratio within about `3%` of
their electric-power-per-electron ratio. The clear exception is phase
`0.375--0.5`: electric power per electron remains `0.93081` of eduPIC while
ionization falls to `0.79584`, giving an ionization-to-power ratio of `0.85499`.
This isolates an energy-selective redistribution or transport effect during
that octant rather than a uniform heating-amplitude error.

### Direct internal electron-energy transport

The required generic, checkpoint-safe surface diagnostic is documented in
[`phase-surface-flux.md`](phase-surface-flux.md). Its first CCP use was locked
before execution in
[`benchmarks/ccp/edupic-argon-surface-flux-rule-20260822.json`](../benchmarks/ccp/edupic-argon-surface-flux-rule-20260822.json).
The four-cycle continuation placed surfaces at `x/L=0.2` and `0.6`, retained
200 RF phase bins and 0.25 eV energy bins, and ran serially at nice level 10.
It completed in `379.35 s` with `198284 KiB` peak RSS and at most `243214`
live macro-particles. All surface-flux shape, finite-value, histogram-closure,
crossing-sufficiency, field, particle, memory, and sampling gates passed. The
surface totals contain `118462` and `1004429` macro-crossings. The retained
global EEDF observation gate again failed only because an electrode-adjacent
region is empty in some sheath phases; this does not affect either interior
surface and is recorded in
[`benchmarks/ccp/edupic-argon-surface-flux-block-20260822.json`](../benchmarks/ccp/edupic-argon-surface-flux-block-20260822.json).

The checksum-bound same-block closure is
[`benchmarks/ccp/edupic-argon-surface-flux-audit-20260822.json`](../benchmarks/ccp/edupic-argon-surface-flux-audit-20260822.json).
Over phase `0.125--0.5`, the directly measured outward electron kinetic-energy
flux divergence is `20.12 W/m2`, versus `21.29 W/m2` independently inferred
from electric power, collision exchange, and kinetic-energy storage: a `5.51%`
closure error. In the exceptional `0.375--0.5` octant, direct transport is
`110.47 W/m2` versus `112.22 W/m2` inferred, closing to `1.55%`. Across all
eight octants the largest absolute difference is only `2.52 W/m2`.

The exceptional octant also carries a positive approximate `16.56 W/m2`
outward divergence in electrons above the `15.8 eV` ionization threshold; the
tail energy uses 0.25 eV histogram-bin centers and had zero overflow. This
shows that the former balance residual is a real particle-transport signal,
not an event-ledger artifact, and that energetic-electron transport materially
contributes during the phase with the extra ionization deficit. It does not
yet prove that AuroraPIC transport is excessive relative to eduPIC: the
published reference has no matching internal crossing spectrum. The result is
a strong same-code conservation and mechanism-localization milestone, not a
cross-code flux validation or experimental validation.

### Direct-transport timestep sensitivity

The direct transport result was next subjected to a prospectively locked 2:1
timestep refinement from the same portable cycle-32 particle state. The rule
is retained in
[`benchmarks/ccp/edupic-argon-surface-flux-timestep-rule-20260822.json`](../benchmarks/ccp/edupic-argon-surface-flux-timestep-rule-20260822.json).
Each branch discarded two equilibration cycles and measured two fresh cycles
with identical 200-bin RF, internal-surface, and energy-spectrum contracts.
The baseline used 4000 steps per RF cycle and the refined branch used 8000.
They ran sequentially on one nice-level-10 CPU; both branches passed energy,
field, EEDF, surface-shape, finite-value, histogram-closure, crossing-count,
and resource gates. Their immutable reports are
[`baseline`](../benchmarks/ccp/edupic-argon-surface-flux-dt-baseline-20260822.json)
and
[`half-dt`](../benchmarks/ccp/edupic-argon-surface-flux-dt-half-20260822.json).

The paired result is retained in
[`benchmarks/ccp/edupic-argon-surface-flux-timestep-result-20260822.json`](../benchmarks/ccp/edupic-argon-surface-flux-timestep-result-20260822.json).
All five prospective gates pass. Halving the timestep changes the directly
measured `0.125--0.5` transport divergence by `4.69%` (limit `15%`) and the
exceptional `0.375--0.5` divergence by only `0.067%` (limit `10%`). The
exceptional above-15.8-eV contribution changes by `2.67%` (limit `20%`).
Baseline/refined direct-versus-inferred closure errors are `3.67%`/`3.99%`
over `0.125--0.5` and `0.095%`/`0.686%` in the exceptional octant, within the
declared `10%` and `6%` ceilings. Neither spectrum has overflow above 80 eV.

This rules out material ordinary timestep error as the source of the direct
transport finding at the declared resolution and sampling horizon. It does not
establish formal temporal order or asymptotic convergence. Absolute values in
this paired common-state experiment need not equal the adjacent cycle-18--22
continuation because the state, seed, equilibration history, and measurement
window differ; only the within-pair timestep change carries the refinement
claim.

### Direct-transport mesh sensitivity

The next prospective rule continued the already paired 400-node and 799-node
fixed-particle checkpoints for two fresh RF cycles, avoiding a redundant
equilibration replay while resetting every measurement accumulator. The rule,
branch reports, and paired result are retained in
[`mesh rule`](../benchmarks/ccp/edupic-argon-surface-flux-mesh-rule-20260822.json),
[`400-node branch`](../benchmarks/ccp/edupic-argon-surface-flux-mesh-baseline-20260822.json),
[`799-node branch`](../benchmarks/ccp/edupic-argon-surface-flux-mesh-refined-20260822.json),
and
[`paired mesh result`](../benchmarks/ccp/edupic-argon-surface-flux-mesh-result-20260822.json).
Both branches completed in about `222 s` on one low-priority CPU with peak RSS
below `209 MiB`; every surface, finite-value, histogram, sampling, field,
particle, and resource gate passed.

All five prospective paired gates pass. The directly measured transport over
phase `0.125--0.5` changes from `20.85 W/m2` to `24.40 W/m2`, a `14.55%`
relative change against the declared `20%` limit. The exceptional
`0.375--0.5` result is much more stable: `112.12 W/m2` versus `113.55 W/m2`,
or `1.27%` against a `12%` limit. Its approximate above-15.8-eV component
changes from `17.68 W/m2` to `20.27 W/m2`, or `12.76%` against a `25%` limit.
Direct-versus-inferred closure remains within `4.29%` in the broad critical
window and `0.65%` in the exceptional octant; both spectra have zero overflow.

The result rules out material ordinary field-grid error for the exceptional
transport signal at these declared tolerances. It also reveals that the
broader phase-window difference is appreciably more grid-sensitive and should
not be described as fully mesh independent. This comparison holds total
particle count fixed, so the next direct-transport discriminator is the
same-grid particle-count refinement; neither result establishes asymptotic
spatial order by itself.

### Direct-transport particle-count sensitivity

The fixed-799-node particle refinement was prospectively locked in the
[`particle rule`](../benchmarks/ccp/edupic-argon-surface-flux-particle-rule-20260823.json).
It reuses the completed standard-particle surface branch and continues the
previously equilibrated doubled-particle state for two fresh measurement
cycles. Halving macro weight from `7e8` to `3.5e8` doubles represented phase-
space sampling. The doubled branch completed serially in `436.22 s`, peaked at
`220704 KiB` RSS and `488227` live macroparticles, and passed every declared
surface, histogram, crossing-count, field, particle, and resource gate. Its
immutable report is the
[`doubled-particle branch`](../benchmarks/ccp/edupic-argon-surface-flux-particle-double-20260823.json).

The checksum-bound
[`paired particle result`](../benchmarks/ccp/edupic-argon-surface-flux-particle-result-20260823.json)
passes all five prospective physics gates. Doubling particle count changes the
direct `0.125--0.5` transport from `24.40` to `25.57 W/m2`, or `4.57%`
against a `20%` limit. Exceptional `0.375--0.5` transport changes from
`113.55` to `119.28 W/m2`, or `4.80%` against a `15%` limit, while its
approximate above-15.8-eV component changes from `20.27` to `20.93 W/m2`, or
`3.14%` against a `25%` limit. Direct-versus-inferred closure errors are
`1.80%`/`1.65%` in the broad window and `0.64%`/`0.20%` in the exceptional
octant; neither spectrum overflows its 80 eV range.

Ordinary macroparticle-count sensitivity is therefore not a material
explanation of the direct transport finding at these declared tolerances.
The doubled particles originated from a split state and were collisionally
decorrelated during the preceding paired equilibration, so this is stronger
than an immediate duplicated-sample comparison but weaker than an independent-
seed ensemble. It is not proof of asymptotic particle convergence, cross-code
flux agreement, or experimental validation.

### Direct-transport continuation-seed uncertainty

The next prospective test reused three completed 400-node states whose
collision-RNG continuations used seeds `13507`, `24601`, and `35713`. Each
state already had six RF cycles of independently seeded evolution before two
fresh surface-flux cycles were measured. The immutable
[`seed rule`](../benchmarks/ccp/edupic-argon-surface-flux-seed-rule-20260824.json)
locked all checkpoint and prior-report hashes, execution limits, diagnostics,
and acceptance thresholds before any branch was launched. All three serial
branches passed their surface, histogram, crossing-count, field, particle,
and resource gates; their reports are retained for
[`13507`](../benchmarks/ccp/edupic-argon-surface-flux-seed-13507-20260824.json),
[`24601`](../benchmarks/ccp/edupic-argon-surface-flux-seed-24601-20260824.json),
and
[`35713`](../benchmarks/ccp/edupic-argon-surface-flux-seed-35713-20260824.json).

The checksum-bound
[`ensemble result`](../benchmarks/ccp/edupic-argon-surface-flux-seed-result-20260824.json)
is deliberately a mixed outcome: four of five prospective gates pass. The
exceptional `0.375--0.5` direct transport spans `112.12--122.07 W/m2`, an
`8.49%` relative range against a `15%` limit. Its approximate above-15.8-eV
component spans `17.68--19.73 W/m2`, an `11.06%` range against a `30%` limit.
The ensemble means are `117.16 W/m2` and `18.54 W/m2`, respectively. The
maximum direct-versus-inferred closure error is `0.93%` in this exceptional
octant, well inside the `6%` ceiling. This strengthens the conclusion that the
exceptional energetic-electron transport mechanism is not an artifact of one
continuation RNG history.

The broader `0.125--0.5` direct transport spans `20.85--26.80 W/m2`, a
`25.38%` range that fails the prospectively locked `20%` limit, even though
all member closure errors remain below `6.93%` and pass the `10%` conservation
ceiling. This signal must therefore be described as statistically under-
resolved at a two-cycle horizon; it needs longer blockwise sampling before a
seed-robust magnitude is claimed. The arms share one original particle
realization, so the experiment quantifies continuation-seed scatter rather
than unrestricted initial-condition uncertainty, cross-code agreement, or
experimental validation.

### Four-cycle direct-transport seed uncertainty

The failed two-cycle critical-window gate motivated a longer-window test, not
a relaxed tolerance. The
[`four-cycle rule`](../benchmarks/ccp/edupic-argon-surface-flux-seed-long-rule-20260824.json)
retained the same `20%`, `15%`, `30%`, `10%`, and `6%` limits and continued
all three members from their step-32000 checkpoints. Measurement accumulators
were reset and four entirely fresh RF cycles were collected. All runtime and
diagnostic gates passed for
[`13507`](../benchmarks/ccp/edupic-argon-surface-flux-seed-long-13507-20260824.json),
[`24601`](../benchmarks/ccp/edupic-argon-surface-flux-seed-long-24601-20260824.json),
and
[`35713`](../benchmarks/ccp/edupic-argon-surface-flux-seed-long-35713-20260824.json).
Each branch used one low-priority CPU, completed in `446--456 s`, and stayed
below `191 MiB` peak RSS.

The checksum-bound
[`four-cycle ensemble result`](../benchmarks/ccp/edupic-argon-surface-flux-seed-long-result-20260824.json)
passes all five prospective gates. The broad `0.125--0.5` direct transport now
spans `22.56--25.45 W/m2`, an `11.86%` relative range against the retained
`20%` limit; its ensemble mean is `24.40 W/m2`. The exceptional `0.375--0.5`
transport spans `114.19--117.68 W/m2`, only `3.01%` against `15%`, with a
mean of `115.97 W/m2`. Its approximate above-15.8-eV component spans
`17.43--20.16 W/m2`, or `14.41%` against `30%`, with a mean of `18.92 W/m2`.
Maximum direct-versus-inferred closure errors are `6.01%` in the broad window
and `0.85%` in the exceptional octant.

The critical-window failure at two cycles therefore behaves as finite-window
sampling noise and is resolved at the prospectively tested four-cycle horizon.
The exceptional transport and energetic-tail findings remain substantially
more seed stable. This establishes repeatability across three decorrelated
continuation-RNG histories at the declared horizon; it still does not provide
unrestricted initial-condition uncertainty, asymptotic statistical
convergence, matching eduPIC internal flux data, or experimental validation.

### Direct-transport constrained-microstate uncertainty

The remaining shared-original-state limitation was tested with the locked
source particle state and two independently randomized conditional
microstates. The existing randomization preserves cell occupancy, CIC nodal
density to roundoff, and each cell's empirical velocity-tuple multiset while
redrawing constrained subcell positions and position--velocity pairing. All
three states had six RF cycles of evolution before the prospectively locked
[`microstate rule`](../benchmarks/ccp/edupic-argon-surface-flux-microstate-rule-20260824.json)
measured four fresh cycles. Runtime and diagnostic gates passed for the
[`source state`](../benchmarks/ccp/edupic-argon-surface-flux-microstate-locked-20260824.json),
[`microstate 51949`](../benchmarks/ccp/edupic-argon-surface-flux-microstate-51949-20260824.json),
and
[`microstate 63059`](../benchmarks/ccp/edupic-argon-surface-flux-microstate-63059-20260824.json).

The checksum-bound
[`microstate result`](../benchmarks/ccp/edupic-argon-surface-flux-microstate-result-20260824.json)
passes all five prospective gates. Broad `0.125--0.5` direct transport spans
`23.44--24.23 W/m2`, only `3.32%` against a `20%` limit, with an ensemble
mean of `23.94 W/m2`. Exceptional `0.375--0.5` transport spans
`115.03--116.38 W/m2`, or `1.17%` against `15%`, with a mean of
`115.86 W/m2`. Its approximate above-15.8-eV component spans
`15.64--18.87 W/m2`, or `18.98%` against `30%`, with a mean of
`16.97 W/m2`. Maximum direct-versus-inferred closure errors are `5.86%` in
the broad window and `0.89%` in the exceptional octant.

At the declared four-cycle horizon, the direct transport magnitude is much
less sensitive to constrained particle microstate than the earlier two-cycle
estimate was to sampling duration. This independently strengthens both the
broad and exceptional transport findings. The ensemble remains conditional on
the empirical cell populations, CIC first moments, and velocity tuples; it is
not unrestricted distribution-function uncertainty, cross-code internal-flux
validation, or experimental validation.

### Collision-velocity staggering audit

A direct source audit subsequently found one unrecorded numerical-contract
difference in the eduPIC argon comparison. eduPIC moves particles with its
stored leapfrog velocity, applies boundaries, and passes that same drift
velocity to MCC. AuroraPIC historically synchronized a time-centered velocity
at the new position, applied MCC to it, and then rebuilt the half step. These
schemes converge toward the same continuous-time split as the timestep is
reduced, but they are not identical at finite timestep and the case manifest
should not have classified them as matched.

AuroraPIC now exposes the checkpoint-protected global
`collision_velocity_sampling` choice. The default `time_centered` retains all
legacy deck and checkpoint behavior; `leapfrog_half_step` matches eduPIC's
move-boundary-collision ordering and is selected by the eduPIC case generator.
Configuration, pusher-state, product-staggering, and incompatible-restart
regressions cover the new path.

An exploratory same-state, same-seed one-cycle A/B smoke then ran both choices
serially at low priority. Both completed in about 81 seconds below 191 MiB RSS.
Relative to the legacy branch, `leapfrog_half_step` changed total energy by
`+0.246%`, electron electrical work by `-1.124%`, ionization events by
`-2.292%`, and the final electric-field profile by `1.329%` L2. The ionization
counts were only `480` and `469`, so the event change is comparable to
short-window sampling noise; the smoke does not establish whether agreement
with eduPIC improves. The checksum-bearing result is
[`benchmarks/ccp/edupic-argon-collision-velocity-staggering-smoke-20260824.json`](../benchmarks/ccp/edupic-argon-collision-velocity-staggering-smoke-20260824.json).

The prospectively locked two-cycle discard plus four-cycle, 200-phase
discriminator has now completed from the same portable particle state. It ran
serially at low priority in 576 seconds with a 205 MiB peak resident set. The
maximum instantaneous field was `66.4 kV/m`; the whole-run energy-ledger
residual was `0.228%`, below the declared `0.5%` gate. All execution and
diagnostic gates passed.

The external result rejects the proposed remedy. The candidate/reference
density ratio moved from `1.12069` to `1.12021`, only a `0.00048` reduction in
distance to unity against the predeclared `0.02` materiality threshold. The
phase-binned electron-power-per-particle ratio moved from `0.87515` to
`0.87031`, a small worsening below its `0.03` threshold. The ionization-source
ratio moved from `0.98673` to `0.94258`; its `0.04415` increase in distance to
unity is a prospectively material worsening. Thus collision-velocity
staggering neither explains nor repairs the remaining density/heating gap.
The new option is still required to state and reproduce the finite-timestep
algorithm contract honestly.

The checksum-bound
[`discriminator result`](../benchmarks/ccp/edupic-argon-collision-velocity-staggering-result-20260824.json)
preserves this negative finding. It is an elimination result, not an eduPIC
acceptance pass or experimental validation. The next credibility step should
target the remaining ionization/heating discrepancy with an independently
checkable collision-kinematics or swarm observable, rather than tuning the CCP
case until its global outputs happen to agree.

### Native eduPIC internal-transport cross-code comparison

The missing external discriminator is now measured directly. A deterministic,
passive transform instruments the exact pinned `C/eduPIC.cc` implementation
that produced the 2620-cycle reference checkpoint. It records electron
crossings at `x/L = 0.2` and `0.6` with the same post-kick drift velocity,
200 phase bins, 0.25 eV energy bins, and represented-particle normalization as
AuroraPIC. Three four-cycle continuations used prospectively fixed collision
seeds. Every run was serial, low priority, hard-limited to 256 MiB virtual
memory, and completed in `96.29--98.92 s` at `47,424--47,824 KiB` peak RSS.

All predeclared measurement gates pass. Native eduPIC's exceptional
`0.375--0.5` phase-window outward kinetic-energy divergence is
`116.39 W/m2`; AuroraPIC's constrained-microstate mean is `115.86 W/m2`.
The absolute ratio is `0.9954`, an unexpectedly close `0.46%` difference.
For the approximate above-15.8-eV component, eduPIC gives `20.16 W/m2` and
AuroraPIC `16.97 W/m2`, an absolute ratio of `0.8420`. Native replicate
relative ranges are `2.90%` for exceptional total transport and `6.37%` for
its ionizing tail, both well inside the `15%` and `30%` gates.

The prospectively declared excessive-transport hypothesis is rejected. After
dividing AuroraPIC flux by its independently measured `1.12069`
AuroraPIC/eduPIC electron-density ratio, the exceptional total and tail ratios
are `0.8882` and `0.7513`, respectively. AuroraPIC therefore does not lose too
much electron energy through these interior surfaces per reference-density
equivalent; its ionizing-tail transport is lower. The remaining retained-tail,
heating-per-electron, and ionization discrepancy must be sought in energy-space
formation or redistribution, collision kinematics/rates, sheath interaction,
or another mechanism—not in excessive outward interior transport.

The first attempted diagnostic targeted the separate `Cpp/eduPIC.cpp`
implementation and encountered its incompatible checkpoint ABI before any
measurement began. The kernel OOM-killed that process after it misread the
double-scalar C checkpoint count as roughly 1.08 billion electrons. The
amended prospective rule records this failure, pins the correct implementation,
and adds the hard address-space cap. No scientific output was observed before
the amendment, and the input checkpoint remained hash-identical.

The checksum-bound
[`rule`](../benchmarks/ccp/edupic-native-surface-flux-crosscode-rule-20260824.json),
[`result`](../benchmarks/ccp/edupic-native-surface-flux-crosscode-result-20260824.json),
and [`execution record`](../benchmarks/ccp/edupic-native-surface-flux-crosscode-execution-20260824.json)
preserve the full finding. This is a direct local cross-code diagnostic under
one nominal CCP case, not a published crossing-spectrum comparison,
experimental validation, or proof of general PIC correctness.

### Direct native regional phase-EEDF comparison

The next prospectively locked discriminator directly samples the native
eduPIC energetic distribution rather than inferring it from the published
density and ionization matrices. The exact `C/eduPIC.cc` implementation was
instrumented with the same seven spatial regions, 200 phase bins, every-second-
step cadence, 0.25 eV bins, 80 eV ceiling, and pre-collision synchronized-
velocity convention used by the AuroraPIC regional diagnostic. Three
four-cycle deterministic continuations completed serially in
`101.25--101.83 s` at `49,036--49,344 KiB` peak RSS. Every critical phase-
region bin contains at least `1,034,024` native macro observations, overflow is
zero, and all repeatability and execution gates pass.

The direct comparison confirms the energetic-tail discrepancy. Across
`x/L=0.2--0.6`, phase `0.125--0.5`, AuroraPIC's histogram-folded ionization
frequency is `67.10 kHz` versus native eduPIC's `78.67 kHz`, a ratio of
`0.85296`. The above-15.8-eV population fraction ratio is `0.87182`; the
above-30-eV ratio falls further to `0.78090`. Native three-seed relative ranges
are only `5.55%` for the folded ionization frequency and `3.55%` for the
ionizing-tail fraction. This independently reproduces the earlier EEDF-folded
deficit without relying on eduPIC's published ionization matrix.

The discrepancy is a shape effect rather than a scalar temperature error.
AuroraPIC's mean energy in the same critical scope is `3.62%` higher while its
ionizing population and ionization kernel are lower; the probability-
distribution total-variation distance is `4.82%`. The three critical phase
octants independently give folded-ionization ratios of `0.8595`, `0.8521`,
and `0.8492`, so the result is not caused by one isolated phase bin.

The prospective spatial-development clause is rejected. The deficit is
already present in `x/L=0.1--0.2`, where the folded-ionization ratio is
`0.84197`, rather than appearing only deeper in the `0.2--0.6` interior. This
weakens bulk transport and downstream inelastic cooling as primary causes and
points more strongly toward sheath energization timing, velocity-space
anisotropy, or energetic-particle residence. The earlier observation that
density-normalized energetic crossing flux near `x/L=0.2` is close while the
local EEDF tail fraction is low is not contradictory: crossing flux weights
speed, whereas an EEDF population fraction weights residence.

The EEDF and surface-flux transforms add no random draws. For all three shared
seeds their independently instrumented binaries produce byte-identical final
particle checkpoints and stdout, providing a direct passivity cross-check.
The checksum-bound
[`rule`](../benchmarks/ccp/edupic-native-phase-eedf-crosscode-rule-20260824.json),
[`result`](../benchmarks/ccp/edupic-native-phase-eedf-crosscode-result-20260824.json),
and [`execution record`](../benchmarks/ccp/edupic-native-phase-eedf-crosscode-execution-20260824.json)
retain the evidence. This narrows the mechanism under one matched CCP case; it
does not establish experimental or general PIC validation.

### Direct native velocity-anisotropy comparison

The next prospective test decomposes that tail by component temperature,
longitudinal kinetic-energy share, signed longitudinal population, and mean
longitudinal speed. The AuroraPIC block and three native four-cycle
continuations use the same critical `x/L=0.2--0.6`, phase `0.125--0.5`
window and a locked 15.8 eV threshold. All scientific measurement gates pass:
the candidate has `795,489` critical tail macro observations, each native
member has at least `811,852`, and all three native diagnostic checkpoints are
byte-identical to prior independently instrumented runs.

This comparison produces a useful elimination result. AuroraPIC retains only
`87.13%` of native eduPIC's above-threshold population, yet its tail
longitudinal-energy fraction is `0.38023` versus `0.38246` (ratio `0.99417`),
and its signed directional imbalance is `0.30080` versus `0.30257` (ratio
`0.99418`). The longitudinal, first-transverse, and second-transverse
temperature shares differ by only `-0.00038`, `-0.00022`, and `+0.00060`.
Across the three predeclared phase octants, the largest absolute directional-
imbalance difference is `0.00906`, far below the prospective `0.05` threshold.
Native relative ranges are `0.65%` for longitudinal tail share and `0.56%`
for directional imbalance.

Gross velocity anisotropy and directional sheath timing are therefore not the
leading explanation for the energetic-tail population deficit. The next
discriminator should measure energetic-particle residence/collision history,
for example age since the last wall encounter and cumulative elastic versus
inelastic events for electrons crossing 15.8 eV. This conclusion is local to
one argon CCP case and does not validate the code against experiment.

The checksum-bound
[`rule`](../benchmarks/ccp/edupic-native-phase-anisotropy-rule-20260825.json),
[`result`](../benchmarks/ccp/edupic-native-phase-anisotropy-result-20260825.json),
and [`execution record`](../benchmarks/ccp/edupic-native-phase-anisotropy-execution-20260825.json)
retain the evidence. The AuroraPIC solver completed normally, but a legacy
postprocessor key mismatch occurred before its in-memory peak-RSS record was
written. The record leaves that value null rather than imputing it; the
identical prior physics run measured `196,124 KiB`. The schema adapter is now
covered by the runner test.

### Direct native energetic-particle history comparison

The next prospective discriminator follows the identity of every sampled
electron from the beginning of the same four-cycle continuation. It records
tracked age, time above 15.8 eV, current energetic streak, threshold entries,
elastic/excitation/ionization collisions, and whether the particle was born
during the measurement window. Histories update every electron timestep while
the phase EEDF remains sampled every second step. Native array compaction copies
the history with the moved particle, ionization products receive fresh history,
and no random draws or particle-state mutations are added.

Diagnostic passivity is exact for all three native seeds: final checkpoints,
stdout trajectories, phase histograms, and anisotropy moments are byte-identical
to the earlier independent runs. Native continuations completed in
`104.18--104.81 s` at `80,612--80,836 KiB` peak RSS. AuroraPIC completed in
`376.15 s` at `206,052 KiB`; all execution, sampling, population, overflow,
finite-value, and three-seed repeatability gates pass.

This is another useful elimination result. In the critical `x/L=0.2--0.6`,
phase `0.125--0.5` window, AuroraPIC/native energetic-duty fractions are
`0.15372/0.15380`, a ratio of `0.99950`. Their current-streak fractions differ
by only `0.61%`, tail-entry rates per 1000 tracked steps by `-0.23%`, elastic
collision exposure by `+1.65%`, and excitation exposure by `-1.06%`.
Born-during-window fractions differ by only `0.00095` absolute. None crosses
its prospectively declared mechanism threshold.

AuroraPIC ionization-collision exposure among observed tail particles is
`9.49%` lower. That is directionally consistent with the energetic-population
deficit but remains below the locked `15%` materiality threshold; it is a
subthreshold clue, not a selected mechanism. Overall, electrons that reach the
tail have remarkably similar finite-window persistence, turnover, and collision
histories in both codes even though AuroraPIC has about `13%` fewer energetic
particles as a fraction of its sampled electron population. The next controlled
discriminator should therefore target
promotion into the tail—particularly collision scheduling and ionization-
product handling—rather than energetic-particle retention.

These histories are observation-conditioned and left-censored for particles
already alive at the window origin, so they are not unbiased lifetime
distributions. The checksum-bound
[`rule`](../benchmarks/ccp/edupic-native-phase-history-rule-20260825.json),
[`result`](../benchmarks/ccp/edupic-native-phase-history-result-20260825.json),
and [`execution record`](../benchmarks/ccp/edupic-native-phase-history-execution-20260825.json)
retain the evidence. This remains a one-case code-to-code mechanism diagnostic,
not experimental validation or a proof of general PIC correctness.

### Unconditional energetic-threshold traffic

The next diagnostic removes the tail-conditioned survivor bias by counting
every live electron timestep in the critical phase/space window. It records
energetic occupancy, net transitions between consecutive synchronized
pre-collision states, accepted-collision transitions by process, and energetic
versus subthreshold ionization births. The native transform remains exactly
passive for all three seeds. Native runs completed in `116.18--116.76 s` at
`80,644--80,784 KiB`; AuroraPIC completed in `393.26 s` at `206,220 KiB`, and
all solver/output gates passed.

The raw observation is strong but the prospectively locked campaign is
formally inconclusive. AuroraPIC/native ratios are `0.8713` for energetic
occupancy, `0.8266` for interstep promotion rate, `0.8603` for interstep
demotion rate, `0.8966` for excitation-collision demotion, `0.8396` for
ionization-collision demotion, and `0.8456` for subthreshold birth rate. Thus
the data point toward reduced energetic-population traffic, with promotion
falling more than demotion; they do not support enhanced loss.

However, the global native repeatability gate also included elastic threshold
demotions. Only `4`, `4`, and `0` such events occurred, producing a `150%`
relative range against the predeclared `30%` limit. The high-population primary,
excitation, ionization, birth, and occupancy repeatability gates all pass, but
the locked analyzer correctly sets `interpretation_allowed` to false. The
promotion-limited reading is therefore a strong exploratory observation, not a
confirmed prospective mechanism result.

The checksum-bound
[`rule`](../benchmarks/ccp/edupic-native-threshold-crossing-rule-20260825.json),
[`result`](../benchmarks/ccp/edupic-native-threshold-crossing-result-20260825.json),
and [`execution record`](../benchmarks/ccp/edupic-native-threshold-crossing-execution-20260825.json)
preserve both the signal and the failed gate. A clean confirmation must use
independent AuroraPIC microstates and a population-qualified repeatability rule
declared before those outputs are observed; the original gate will not be
silently weakened.

### Confirmatory threshold-traffic microstate replication

That clean confirmation is now complete. Before observing either new ledger,
the replication rule locked two pre-existing constrained AuroraPIC microstates,
the same critical `x/L=0.2--0.6`, phase `0.125--0.5` window, an 8% candidate
repeatability limit, population floors, and the already observed three-seed
native ensemble by hash. Both four-cycle AuroraPIC continuations passed every
execution, resource, sampling, shape, finite-value, closure, and population
gate. They completed serially in `403.44 s` and `402.70 s`, with peak RSS of
`205,888 KiB` and `205,744 KiB`.

The prospective mechanism result passes. The two AuroraPIC/native promotion-
rate ratios are `0.81265` and `0.80066`, both below the locked `0.90` boundary.
The corresponding demotion-rate ratios are `0.84982` and `0.84261`, both below
the `1.10` enhanced-loss boundary. Candidate relative ranges are `1.69%` for
energetic occupancy, `1.49%` for promotion, `0.85%` for demotion, `4.06%` for
excitation demotion, `1.12%` for ionization demotion, and `0.94%` for
subthreshold births; all are below 8%. Excitation and ionization comparisons
also exceed the locked 500-event floor in every candidate and native member.
Elastic threshold demotions remain descriptive because their tiny native count
does not qualify for relative comparison.

Within this matched argon CCP, the energetic-tail deficit is therefore robust
to the two constrained particle realizations and is associated with reduced
promotion traffic, not enhanced demotion. This is a useful numerical-physics
finding: the next diagnostic should resolve the promotion interval into field-
push work and preceding collision/ionization-product effects, then test the
responsible implementation choice directly. Interstep transitions are not yet
a field-only attribution because each interval spans the previous accepted
collision and the following push.

The checksum-bound
[`rule`](../benchmarks/ccp/edupic-threshold-crossing-replication-rule-20260825.json),
[`result`](../benchmarks/ccp/edupic-threshold-crossing-replication-result-20260825.json),
and [`execution record`](../benchmarks/ccp/edupic-threshold-crossing-replication-execution-20260825.json)
retain the evidence. The native ensemble was reused and already observed, so
this confirms candidate-microstate robustness rather than an independent
native replication. It remains code-to-code evidence for one case—not
experimental validation or proof of general PIC correctness.

A subsequent diagnostic audit found and corrected an AuroraPIC first-sample
classification bug: history age was incremented before checking whether a
previous synchronized state existed. This added `348` and `311` false
promotions across all regions of phase bin 0 in the two runs. It does not alter
the locked result because the critical analyzer selects phase fractions
`0.125--0.5`; phase bin 0 contributes no selected rows. The
[`first-observation audit`](../benchmarks/ccp/edupic-threshold-crossing-first-observation-audit-20260825.json)
records the affected and corrected solver hashes, counts, scope proof, and
regression-test status.

### Direct field-push threshold attribution

The next prospective campaign isolates the mover-stage interval rather than
inferring it from consecutive pre-collision states. For each surviving
electron push it compares the post-collision velocity entering the mover with
the pre-collision velocity leaving it, attributes the event at the post-drift
position, and counts threshold promotions and demotions. AuroraPIC preserves
these accumulators in checkpoint v21. The independent native transform adds no
random draws or particle-state mutations; all three final native checkpoints
remain byte-identical to prior passive runs.

All locked gates pass. The two AuroraPIC branches completed in `408.19 s` and
`395.95 s` at about `206 MiB` peak RSS. Native seeds completed in
`105.58--112.15 s` at `80,620--80,876 KiB`. Critical-window populations are
approximately `525 million` candidate pushes and `475 million` native pushes
per member, with more than 2,500 demotions and 4,100 promotions in every
candidate and more than 2,700 demotions and 4,600 promotions in every native
member. Candidate relative ranges are `1.56%` for promotions and `2.61%` for
demotions; native ranges are `3.04%` and `4.21%`, within the prospective `8%`
and `15%` limits.

The field-push promotion-deficit outcome is supported. AuroraPIC/native
promotion-rate ratios are `0.80585` and `0.79340`, remarkably close to the
earlier whole-interstep ratios of `0.81265` and `0.80066`. Field-push demotion
ratios are also lower, `0.83234` and `0.81092`; enhanced mover-stage loss is
not supported. Thus, under the currently compared algorithm contracts, the
missing energetic-tail traffic is already present across the mover-stage
interval rather than being introduced primarily by the accepted-collision
stage.

This result deliberately stops short of blaming field interpolation or the
leapfrog mover. These constrained AuroraPIC checkpoints use time-centered MCC
velocity, whereas native eduPIC applies MCC to its half-step velocity. Each
ledger isolates its own configured mover interval, but their staggered
before/after states are not identical. The next decisive experiment is an
AuroraPIC matched-half-step continuation with the same field-push diagnostic;
only persistence there would localize the discrepancy beyond collision-
velocity staggering.

The checksum-bound
[`rule`](../benchmarks/ccp/edupic-field-push-threshold-rule-20260825.json),
[`result`](../benchmarks/ccp/edupic-field-push-threshold-result-20260825.json),
and [`execution record`](../benchmarks/ccp/edupic-field-push-threshold-execution-20260825.json)
retain the evidence and the two non-scientific execution corrections. This is
one-case code-to-code mechanism evidence, not experimental validation or proof
of general PIC correctness.

### Matched-half-step continuation

The predeclared follow-up removed the velocity-staggering ambiguity without
overriding checkpoint identity. Live positions and time-centered velocities
from each locked step-24000 checkpoint were exported to checksum-bound APS v2
states. AuroraPIC rebuilt the adjacent half-step velocity from the
self-consistent field, used `leapfrog_half_step` collision sampling, relaxed
for one complete RF cycle, and then measured the same four-cycle field-push
window. Both serial branches passed every runner, population, safety, and
provenance gate with zero EEDF overflow and about `206 MiB` peak RSS.

The directional evidence does not support velocity staggering as the missing
promotion mechanism: matched-half-step AuroraPIC/native promotion-rate ratios
are `0.76039` and `0.82876`, both below the prospective persistence boundary
of `0.90` and farther from parity than the time-centered ratios. Demotion-rate
ratios are `0.76548` and `0.84844`. This is not yet a formal persistence result,
however. Promotion and demotion relative ranges across the two microstates are
`8.60%` and `10.28%`, narrowly exceeding the locked `8%` repeatability limit.
The analyzer therefore correctly forbids interpretation and records only a
directional persistence signal. A prospectively pooled second four-cycle block
is the smallest follow-up that can reduce finite-window uncertainty without
changing the diagnostic, critical window, or decision boundaries.

The checksum-bound
[`rule`](../benchmarks/ccp/edupic-matched-half-step-threshold-rule-20260825.json),
[`result`](../benchmarks/ccp/edupic-matched-half-step-threshold-result-20260825.json),
and [`execution record`](../benchmarks/ccp/edupic-matched-half-step-threshold-execution-20260825.json)
preserve the gated outcome. This remains one-case code-to-code mechanism
evidence, not published-observable or experimental validation.

The prospectively declared replication then continued both half-step
checkpoints through a second non-overlapping four-cycle block and pooled raw
counts within each microstate. Every second-block runner and population gate
passed. The pooled windows contain approximately `1.047 billion` critical
pushes per microstate, with `8,289` and `8,755` promotions. Promotion and
demotion relative ranges fell to `5.46%` and `6.50%`, passing the unchanged
`8%` repeatability gate. Pooled AuroraPIC/native promotion-rate ratios are
`0.79137` and `0.83575`; both pass the prospective persistence criterion of at
most `0.90` and neither approaches the `0.95--1.05` parity interval. The formal
outcome is therefore that the field-push energetic-promotion deficit persists
when AuroraPIC uses matched leapfrog-half-step collision sampling. Velocity
staggering is excluded as the explanation under this case contract; the next
mechanism audit can focus on field interpolation, push timing, and the
self-consistent field/particle-state coupling.

The pooled replication
[`rule`](../benchmarks/ccp/edupic-matched-half-step-replication-rule-20260825.json),
[`result`](../benchmarks/ccp/edupic-matched-half-step-replication-result-20260825.json),
and [`execution record`](../benchmarks/ccp/edupic-matched-half-step-replication-execution-20260825.json)
are checksum-bound. Pooling extends the observation window but does not add
independent microstates, so this is a mechanism-localization milestone rather
than general or experimental validation.

### Mover-contract audit

A checksum-locked source audit now compares the pinned eduPIC 1.0 mover with
AuroraPIC's 1D electrostatic path. The bulk contracts agree: both use centered
interior nodal fields, linear cloud-in-cell interpolation, the same electron
kick algebra, and kick-then-drift leapfrog ordering. No bulk mover-formula
mismatch was found.

Two narrower differences remain. eduPIC advances time and solves the field
before each kick, while AuroraPIC kicks with the retained current-time field
and solves at the end of the step. At 4,000 steps per RF cycle this is a phase
offset of `0.0015707963 rad`, or at most `0.3927 V` on the `250 V` waveform.
That scale is too small to assume it explains a `16--21%` promotion deficit,
but it is cleanly testable by advancing AuroraPIC's configured phase by one
timestep. Second, eduPIC's two electrode-node fields include the charged
half-cell Gauss correction `rho dx/(2 epsilon)`; AuroraPIC currently uses only
the one-sided potential gradient. This directly affects wall-adjacent
interpolation cells and deserves a separate prospective solver branch if the
phase-aligned test cannot restore parity.

The reproducible
[`audit`](../benchmarks/ccp/edupic-mover-contract-audit-20260825.json)
is a static source-contract and scale result, not a dynamic equivalence or
validation claim.

The predeclared phase-aligned control advanced only `phi_left_phase` by the
audited one-timestep offset, reloaded each locked step-36000 half-step
checkpoint, relaxed one complete RF cycle, and measured four complete cycles.
Both serial branches passed every provenance, population, sampling, field,
and resource gate. They collected approximately `520.5 million` critical
pushes each with zero EEDF overflow and about `210 MiB` peak RSS. Promotion
rates are `8.01771` and `7.96845` per million pushes, only `1.28%` above and
`4.69%` below their respective pooled baselines. Their native eduPIC ratios
are `0.80148` and `0.79656`, both below the locked `0.90` persistence boundary
and far outside the `0.95--1.05` parity interval. Promotion and demotion
relative ranges are `0.62%` and `3.39%`, comfortably passing the `8%`
repeatability gate.

The formal result is therefore that the one-electron-timestep RF field-time
convention does not explain the energetic-promotion deficit. Combined with
the source audit, this excludes the bulk mover algebra, MCC velocity
staggering, and the RF phase offset under this case contract. The next isolated
solver mechanism is the electrode-node charged half-cell Gauss correction;
that requires a prospectively locked code branch because it changes the field
seen in the wall-adjacent cells.

The checksum-bound phase-alignment
[`rule`](../benchmarks/ccp/edupic-phase-aligned-mover-rule-20260825.json),
[`result`](../benchmarks/ccp/edupic-phase-aligned-mover-result-20260825.json),
and [`execution record`](../benchmarks/ccp/edupic-phase-aligned-mover-execution-20260825.json)
preserve this mechanism result and the one pre-timestep runner correction. It
remains one-case code-to-code evidence, not experimental validation or proof
of general PIC correctness.

The electrode-node discrepancy has now also been corrected in the production
1D Dirichlet solver. Its endpoint fields integrate Gauss's law across the
half-volume represented by each electrode node, while the interior potential
and field formulas remain unchanged. A focused solver test verifies both
correction signs and magnitudes and guards against endpoint charge leaking
into the interior Dirichlet potential solution.

The prospectively locked dynamic control loaded the same step-36000 particle
and RNG states with the corrected solver, relaxed for one RF cycle, and
measured four cycles. Both branches passed every gate with approximately
`520.7 million` critical pushes each, zero EEDF overflow, and about `210 MiB`
peak RSS. Corrected promotion rates are `8.27668` and `8.06922` per million;
their native eduPIC ratios are `0.82737` and `0.80663`. Both remain below the
locked `0.90` persistence boundary. The correction changed the two finite
windows by `+4.55%` and `-3.48%` relative to their pooled baselines, while
promotion and demotion relative ranges of `2.54%` and `5.21%` pass the `8%`
repeatability gate.

Thus the physically required endpoint Gauss correction is retained as a
production numerical fix, but it does not explain the energetic-promotion
deficit under this CCP contract. The checksum-bound endpoint-control
[`rule`](../benchmarks/ccp/edupic-endpoint-gauss-control-rule-20260825.json),
[`result`](../benchmarks/ccp/edupic-endpoint-gauss-control-result-20260825.json),
and [`execution record`](../benchmarks/ccp/edupic-endpoint-gauss-control-execution-20260825.json)
preserve the formal result.

### Frozen-state operator equivalence

The mover investigation now separates operator correctness from differences
in the evolved plasma distribution. A compiled deterministic test evaluates a
401-node frozen Poisson state using an independent transcription of the
checksum-pinned eduPIC Thomas solve, endpoint fields, CIC interpolation,
electron kick, and drift. The same field and 20,000 identical particle states
are evaluated through AuroraPIC's production solver and pusher. One hundred
samples lie in the two wall-adjacent interpolation cells, and the velocity set
deliberately exercises both directions across the `15.8 eV` threshold.

All equivalence gates pass. Maximum discrepancies are `2.84e-14 V` in
potential, `2.27e-10 V/m` in nodal field, `2.29e-10 V/m` in interpolated
field, `9.31e-10 m/s` after the kick, and `3.47e-18 m` after drift. Both
implementations classify exactly `10,000` promotions and `10,000` demotions.
The test runs in less than a second and is now a permanent CTest regression.

This proves formula-level frozen-state equivalence of the compared 1D
electrostatic Poisson/interpolation/kick/drift operator to roundoff. It also
changes the interpretation of the dynamic promotion deficit: the deficit is
not evidence of an incorrect instantaneous AuroraPIC mover. It must arise
from state evolution before the measured push, such as collision-channel
sampling, collision kinematics, particle creation/loss history, or the
self-consistent distribution those mechanisms produce. A frozen collision-
operator comparison is therefore the next localization step.

The checksum-bound frozen-operator
[`result`](../benchmarks/ccp/edupic-frozen-operator-result-20260825.json) and
[`execution record`](../benchmarks/ccp/edupic-frozen-operator-execution-20260825.json)
retain the numerical errors, exact threshold counts, source and binary hashes,
and the pre-result tolerance correction. This is an operator verification,
not a dynamic discharge or experimental validation result.

### Frozen collision-rate localization

The collision investigation next compares the checksum-pinned native source,
AuroraPIC collision implementation, and all five local generated argon tables
without evolving either plasma. The source contracts agree on electron and
ion energy frames, threshold losses, isotropic/ backward angular models,
finite-mass center-of-mass transforms, and the Opal-Beaty-Peterson ionization
partition with a `10 eV` scale. Three known algorithm differences remain:
eduPIC uses nearest `0.001 eV` bins, AuroraPIC uses lower bins; eduPIC permits
one Bernoulli event per species step, while AuroraPIC's null-collision clock
permits repeated Poisson opportunities; and the random engines differ.

An upper-half-bin worst-case scan over `0--80 eV` finds a maximum electron
event-probability difference of `3.07e-6` per electron step. The ion maximum is
larger, `1.96e-4` per ion step at `0.00175 eV`, where the argon-ion cross
section is steep. The isolated worst point is not representative of the
particle population, so both kernels were then folded over every live
half-step velocity in the two locked step-36000 checkpoints. Ion-neutral
thermal motion at `350 K` is integrated deterministically with a 27-node
three-dimensional normal quadrature.

The actual-state result is tightly replicated. AuroraPIC lower-bin/native
nearest-bin total first-event ratios are `0.999869` and `0.999871` for
electrons, and `1.006676` and `1.006678` for ions. Ionization-channel ratios
are `0.999755` and `0.999775`. Including AuroraPIC's repeated-opportunity
Poisson mean gives total/native ratios of `1.00273--1.00278` for electrons
and `1.00719` for ions. Thus lookup and opportunity-clock differences alter
direct frozen traffic by less than one percent in both species and cannot
directly supply a `16--20%` collision-acceptance deficit.

The checksum-bound frozen collision
[`kernel audit`](../benchmarks/ccp/edupic-frozen-collision-kernel-audit-20260825.json),
[`actual-state fold`](../benchmarks/ccp/edupic-frozen-collision-state-fold-20260825.json),
and [`execution record`](../benchmarks/ccp/edupic-frozen-collision-execution-20260825.json)
retain the result without redistributing the locally generated gas tables.
Collision-product velocity distributions remain the next statistical
operator comparison; these frozen rate results do not prove that the evolved
EEDFs must agree.

### Frozen collision-product moments

The next bounded test drives AuroraPIC's compiled collision implementation
and compares isolated product distributions with the closed-form moments of
the checksum-pinned eduPIC laws. Each stochastic channel uses `12,000`
accepted single-collision samples. The test covers finite-mass isotropic
electron elastic and excitation events, finite-mass Opal-Beaty-Peterson
ionization, and equal-mass argon-ion isotropic and backward scattering.

The ideal isotropic moments are `E[cos(chi)] = 0` and
`E[cos(chi)^2] = 1/3`. AuroraPIC gives `(0.00217, 0.33629)` for elastic,
`(-0.00938, 0.33587)` for excitation, and `(0.00464, 0.33281)` for ion
isotropic scattering. Equal-mass backward scattering performs the expected
velocity exchange. The sampled Opal mean ejected energy is `1.03333` against
the analytic `1.03504`, a `0.165%` difference; electron-pair energy and
paired-angle momentum invariants close below `6.3e-15` and `1.4e-15`.

The checksum-bound
[`collision-product result`](../benchmarks/ccp/edupic-collision-product-moments-20260825.json)
therefore excludes a large direct product-kinematics error in the tested
laws. Combined with the actual-state rate fold, the remaining dynamic deficit
is localized to accumulated stochastic evolution and its self-consistent
particle creation, loss, and field feedback, rather than the frozen mover,
collision traffic, or isolated product distributions. This is not yet a
whole-discharge or experimental validation result.

### Coupled population-balance audit

A checksum-verified post hoc audit reuses the two completed, consecutive
matched-half-step blocks to test whether ordinary AuroraPIC microstate
variation can explain the persistent cross-code deficit. Only files already
bound by each runner report are admitted: scalar states, collision counts,
and boundary losses. The unbound historical power-transfer files are
deliberately excluded.

Species accounting closes exactly. In microstate 51949, `1,853`
ionizations minus `2,058` electron wall losses gives the observed `-205`
electron change. In microstate 63059, `1,826 - 2,067 = -241`. Ion balances
also close exactly, including the unequal electron and ion wall losses in the
second branch.

Gross coupled evolution is substantially more repeatable than the cross-code
energetic-tail deficit. Cross-branch relative ranges are `1.47%` for
ionization, `0.44%` for electron wall loss, and `2.64%` for ion wall loss.
At five common RF-phase endpoints, mean kinetic and field energies differ by
only `0.25%` and `0.28%`; their final differences are `0.98%` and `0.50%`.
The prospectively pooled promotion-rate range remains `5.46%`, below its
locked `8%` gate and well below the `16--20%` code-to-code deficit.

The [`coupled-evolution result`](../benchmarks/ccp/edupic-coupled-evolution-audit-20260825.json)
therefore does not support within-Aurora microstate variability as the
explanation. The next comparison needs time-resolved native state variables
or a deliberately matched state-transition experiment; extending the same
two AuroraPIC branches alone has diminishing diagnostic value.

### Native cycle-history contract

A new passive native transform now records one checksum-ready state row after
boundary handling and collisions at the end of every RF cycle. It composes
the existing threshold ledger and adds species populations, energetic
population above `15.8 eV`, electron collision channels, ion collision total,
wall losses, kinetic and field energies, charge L1, and maximum field. The
transform adds no random draws and does not mutate particle state.

The one-cycle seed-13507 smoke test is trajectory-exact: deterministic control
and instrumented checkpoints have the identical SHA-256
`cb489136...21cc2ce`. Its first state row also closes population accounting:
`1000 + 2262 - 2458 = 804` electrons and
`1000 + 2262 - 14 = 3248` ions. The
[`smoke record`](../benchmarks/ccp/edupic-cycle-history-smoke-20260825.json)
retains the source, instrumenter, binary, checkpoint, and output hashes.

Before observing any production-state cycle-history output, a
[`prospective rule`](../benchmarks/ccp/edupic-cycle-history-rule-20260825.json)
locks the common cycle-2620 native checkpoint, three seeds, their previously
observed passive final-checkpoint hashes, four-cycle horizon, serial resource
limits, accounting gates, repeatability thresholds, and interpretation
boundaries. This comparison can identify differing time-resolved observables;
because the native and AuroraPIC states evolved independently, it cannot by
itself assign causal onset to one algorithm.

All three preregistered native members subsequently completed serially in
`101.97--102.74 s` with at most `80,728 KiB` RSS. Every final checkpoint is
byte-identical to its previously locked passive run, and all twelve per-cycle
species balances and electron-channel sums close exactly. Native relative
ranges are `2.90%` for ionization, `1.54%` for electron wall loss, `2.01%`
for ion wall loss, and `2.84%` for final energetic fraction, passing the
prospective `8%` gate.

The [`cycle-history result`](../benchmarks/ccp/edupic-cycle-history-result-20260825.json)
reveals a robust distribution-level difference. Normalized by starting
electrons, AuroraPIC/native four-cycle ionization is `0.8693`, whereas the
electron-wall-loss ratio is `0.9685`. AuroraPIC therefore replenishes about
`13.1%` fewer electrons through ionization; its declining population is not
caused primarily by excessive wall loss. At the same time its total kinetic
energy is `12.2%` higher than native, its field energy is within `1.6%`, and
its global final energetic fraction is `6.4%` lower. This combination points
to an EEDF-shape or energy-partition problem: AuroraPIC has more total kinetic
energy but less ionizing-tail supply. Lower ionization is both a consequence
and a feedback amplifier of that tail difference, so this result still does
not identify the first causal divergence.

A deterministic decomposition of the five locked final checkpoints then
qualified the total-energy statement. AuroraPIC has `9.29%` more electrons
and `8.95%` more ions than the independently evolved native ensemble, which
explains most of its larger population-weighted kinetic energy. Its mean
electron energy is only `4.17%` higher and the increase is nearly isotropic:
the per-particle x, y, and z energy ratios are `1.045`, `1.040`, and `1.040`.
Ion mean energy is `1.6%` lower, excluding ion overheating and a directional
heating defect at this endpoint.

The [`checkpoint energy partition`](../benchmarks/ccp/edupic-checkpoint-energy-partition-20260825.json)
also finds that total ionizing-band energy per area agrees within `1.6%`, even
though the fraction of electrons above `15.8 eV` is `6.4%` lower in
AuroraPIC. Therefore the `13.1%` normalized ionization deficit is not explained
by a missing global endpoint ionizing-energy reservoir. The next localization
must resolve where that reservoir resides in phase, space, and time, and how
often its electrons sample the ionization cross section.

### Ionization-exposure synthesis

That phase-space localization had already been completed independently under
the earlier native phase-EEDF rule. Its prospectively locked regional EEDF
fold, using the common checksum-bound ionization cross section, predicts an
AuroraPIC/native ionization-frequency ratio of `0.8530` over `x/L=0.2--0.6`
and RF phase `0.125--0.5`. The later cycle histories independently measure
`0.8693` realized ionizations per starting electron. These differ by only
`0.0164` in absolute ratio, or `1.92%` relative.

The checksum-locked
[`ionization-exposure synthesis`](../benchmarks/ccp/edupic-ionization-exposure-synthesis-20260825.json)
connects those results with the endpoint partition. It shows why a global
energy scalar was misleading: ionizing-band energy per area is near parity at
`1.015`, yet the endpoint ionizing population fraction is `0.936`, and the
phase-region cross-section-weighted exposure is `0.853`. Meanwhile electron
wall loss is much closer to parity at `0.969`. Thus the observed event deficit
is quantitatively tracked by where and when electrons sample the ionization
kernel, not by excessive wall loss or a missing global endpoint energy
reservoir.

This is meaningful cross-diagnostic evidence, not a new formal acceptance
gate. The EEDF and cycle histories come from independently evolved microstates
and the synthesis is post hoc, so their close agreement does not prove causal
onset, experimental validity, or general PIC correctness. It does sharply
narrow the remaining problem: locate the earliest divergence in energetic
electron production and transport before the common collision kernel is
sampled.

### Matched-half-step promotion localization

The accepted matched-half-step campaign already isolated the reduced
above-`15.8 eV` traffic to the field-push interval, but its formal result
pooled the entire critical phase/space window. A fail-closed postprocessor now
reuses the four checksum-bound AuroraPIC ledgers and three native ledgers to
localize that existing result. It first exactly reproduces the prospectively
accepted pooled promotion ratios, `0.79137` and `0.83575`, before examining any
subdivision.

The deficit spans all three critical RF-phase octants. It is strongest over
phase `0.25--0.375`, where the two AuroraPIC/native promotion-rate ratios are
`0.7549` and `0.7916`; the adjacent octants give `0.7783/0.8354` and
`0.8519/0.8977`. Over phase `0.125--0.5`, the upstream `x/L=0.1--0.2` ratios
are `0.8902` and `0.9038`, then fall to `0.7729` and `0.8011` in
`x/L=0.2--0.4`. Candidate relative ranges in that latter region are only
`3.58%` for promotion and `3.69%` for demotion. The farther `x/L=0.4--0.6`
subdivision has `15.9%` candidate promotion dispersion and is not used for a
strong localized inference.

The checksum-locked
[`promotion-localization result`](../benchmarks/ccp/edupic-field-push-promotion-localization-20260825.json)
therefore shows that the mover-stage discrepancy is distributed in phase and
amplifies immediately beyond the upstream region; it is not generated by one
isolated RF bin. These subdivisions were selected after the aggregate result
was known, so they remain descriptive. The next causal test should compare the
field work delivered to near-threshold electrons across the `x/L=0.2`
interface during phase `0.25--0.375`, rather than launch another broad
whole-discharge continuation.

The required causal diagnostic is now implemented in both codes. AuroraPIC
checkpoint v22 records a configurable subthreshold band's observation count,
promotion count, and signed/positive/negative mover work without changing the
particle trajectory. A composed passive transform adds the same quantities to
the pinned native C implementation for the argon excitation-to-ionization band
`11.5--15.8 eV`.

The native one-cycle
[`promotion-band smoke`](../benchmarks/ccp/edupic-promotion-band-work-smoke-20260825.json)
is trajectory-exact against the previously verified field-push binary: final
checkpoint, stdout, phase EEDF, threshold ledger, and existing field-push
ledger hashes are identical. The new ledger contains `4,749,306` band
observations and `3,902` promotions, and its signed work equals positive minus
negative work to `6.4e-12 eV`. Peak RSS remains about `80.4 MiB` under the
`256 MiB` address-space cap.

Before observing any production-state work output, the
[`prospective rule`](../benchmarks/ccp/edupic-promotion-band-work-rule-20260825.json)
locks two existing matched-half-step AuroraPIC continuations, three native
seeds, all solver/instrumenter/checkpoint hashes, the `x/L=0.2--0.4`, phase
`0.25--0.375` scope, population and repeatability gates, and three possible
mechanisms: reduced supply into the band, reduced positive work per band
electron, or reduced conversion across `15.8 eV`. No physics outcome exists
until every locked member and gate is complete.

The three native production members have now completed serially. Every final
checkpoint matches its preregistered passive hash, each critical member has
more than `1.06 million` band observations and `1,850` promotions, and peak
RSS is at most `80,628 KiB`. Native relative ranges are `0.45%` for band
supply, `1.84%` for promotion probability, and `3.55%` for mean positive
work, passing the locked native limits. The
[`native result`](../benchmarks/ccp/edupic-promotion-band-work-native-result-20260825.json)
records the exact hashes and quantities. This completes only the native half;
cross-code interpretation remains forbidden until both AuroraPIC members run.

Both locked AuroraPIC members have now completed, serially and without a guard
failure. Each contributes more than `1.08 million` critical band observations
and `1,619` promotions; peak RSS is below `216 MiB`. All integrity,
population, work-closure, and repeatability gates pass. The
[`prospective cross-code result`](../benchmarks/ccp/edupic-promotion-band-work-result-20260825.json)
therefore permits the declared interpretation. Relative to the three-seed
native mean, the two AuroraPIC mean-positive-work ratios are `0.8795` and
`0.8870`, while the conditional band-promotion-probability ratios are `0.8294`
and `0.8615`. Both mechanisms satisfy the locked `<=0.90` rule. The band-supply
ratios are `0.9217` and `0.8971`; because both do not satisfy the rule, a band
supply deficit is not supported. The ensemble-mean ratios are `0.8832`,
`0.8454`, and `0.9094`, respectively.

This is a meaningful localization result: in the strongest previously
identified region/phase window, AuroraPIC supplies nearly as many electrons to
the `11.5--15.8 eV` band but delivers less positive mover work per band
observation and converts fewer of those observations across `15.8 eV`. The
conditional promotion probability also depends on the within-band energy
distribution, so these two supported deficits are not claimed as independent
causes. This one-case solver-to-solver comparison neither compares experiment
nor establishes a published benchmark pass or general PIC correctness.

The next discriminator is now implemented but has no production outcome yet.
For every band observation, checkpoint v23 accumulates origin total and
longitudinal energy and decomposes mover work exactly as
`delta K = m v_x delta v_x + 0.5 m delta v_x^2`. The linear term probes
velocity--field alignment; the quadratic term probes particle-sampled field
strength and is proportional to `E^2` under the locked half-step convention.
The composed native transform's
[`one-cycle smoke`](../benchmarks/ccp/edupic-mover-decomposition-smoke-20260825.json)
is byte-exact for the trajectory and all prior diagnostics. Its two algebraic
closures are below `8.2e-12 eV`, with `81,088 KiB` peak RSS.

Before any production decomposition output is observed, the
[`prospective mover-decomposition rule`](../benchmarks/ccp/edupic-mover-decomposition-rule-20260825.json)
locks the two AuroraPIC states, three native seeds, binaries, instrumenter,
scope, population/repeatability gates, and four possible associations: sampled
field strength, favorable alignment, longitudinal energy partition, and origin
energy distribution. Multiple associations may be supported because the
quantities are decomposition terms or covariates, not independent causal
variables. No physics interpretation is allowed until every locked member and
gate completes.

The three native production members have now completed serially and are sealed
in the
[`native mover-decomposition result`](../benchmarks/ccp/edupic-mover-decomposition-native-result-20260825.json).
Every final checkpoint and prior promotion-work ledger is byte-exact, all
population and closure gates pass, and peak RSS is at most `81,056 KiB`.
Relative ranges are `0.15%` for mean origin energy, `2.40%` for longitudinal
energy fraction, `3.55%` for positive linear work, and `3.87%` for quadratic
work. The native half is therefore qualified and repeatable. Cross-code
interpretation remains forbidden until both locked AuroraPIC members complete.

Both AuroraPIC members have now completed and the
[`prospective mover-decomposition result`](../benchmarks/ccp/edupic-mover-decomposition-result-20260826.json)
passes every joint integrity, population, closure, and repeatability gate. The
two AuroraPIC/native quadratic-work ratios are `0.7796` and `0.8251`, so the
locked sampled-field-strength deficit is supported. The positive-linear-work
ratios are `0.8799` and `0.8873`, so the favorable-alignment deficit is also
supported. In contrast, origin-energy ratios `0.9967/0.9971` and longitudinal
energy-fraction ratios `0.9944/0.9936` are in the declared parity interval;
their deficit rules are not supported. AuroraPIC relative ranges are at most
`5.66%` for the four declared metrics.

This sharpens the previous localization: the near-threshold population enters
the mover with nearly the same mean energy and longitudinal partition, but
experiences a smaller squared velocity increment and less favorable positive
linear work. Under the locked half-step convention these correspond to weaker
particle-sampled field strength and weaker energizing velocity--field
alignment. They are simultaneous associations in an already evolved
self-consistent discharge, not independent causal effects and not proof of
where the two solvers first diverge. The result remains one solver-to-solver
CCP comparison, not experimental or general PIC validation.

### Grid-field versus particle-sampling discriminator

The quadratic mover term establishes a particle-sampled field-strength
difference, but by itself cannot distinguish a weaker grid field from a
different distribution of band electrons within that field. Before reducing
the existing field outputs, the
[`grid-field sampling rule`](../benchmarks/ccp/edupic-grid-field-sampling-rule-20260826.json)
locked both AuroraPIC members, all three native members, exact hashes, the
`x/L=0.2--0.4` and phase `0.25--0.375` window, spatial boundary interpolation,
phase reduction, repeatability gates, and `0.90` decision thresholds.

The
[`prospective result`](../benchmarks/ccp/edupic-grid-field-sampling-result-20260826.json)
passes every hash, shape, positivity, and repeatability gate. The two
AuroraPIC/native grid-mean-square-field ratios are `0.6653` and `0.6766`, and
the ensemble ratio is `0.6710`; the corresponding RMS fields are approximately
`2.00/2.02 kV/m` versus `2.44--2.47 kV/m` natively. Both members therefore
satisfy the locked grid-field-deficit rule. In contrast, their conditional
particle-sampling-factor ratios are `1.1719` and `1.2195`, with an ensemble
ratio of `1.1957`, so the locked differential-sampling-deficit rule is not
supported. AuroraPIC relative ranges are `1.69%` and `3.98%`, and native
ranges are `1.74%` and `2.13%`, respectively.

This is a meaningful diagnostic localization: the approximately `0.8024`
ensemble quadratic-work ratio is associated with a substantially weaker
self-consistent phase-mean grid field, not with AuroraPIC band electrons
preferentially avoiding strong-field positions inside the already localized
window. Indeed, conditional sampling is modestly higher in AuroraPIC and
partly offsets the grid-field difference. This does not yet identify when or
why the grid fields first diverge. The next prospective discriminator should
compare the charge-density/Poisson inputs and sheath-edge evolution at earlier
matched checkpoints and phases, before the mature field deficit is present.
As before, this one-case solver comparison is neither experimental validation
nor proof of general PIC correctness.

### Discrete Poisson source attribution

Both implementations use the same nodal tridiagonal Dirichlet Poisson solve,
centered interior electric field, and half-control-volume Gauss-law endpoint
field. Before evaluating the density outputs, the
[`Poisson source-attribution rule`](../benchmarks/ccp/edupic-poisson-source-attribution-rule-20260826.json)
locked the phase-space files and hashes, exact operator, critical window,
reconstruction tolerance, boundary-drive parity interval, three-factor
substitution, and Shapley decision rule.

The
[`prospective result`](../benchmarks/ccp/edupic-poisson-source-attribution-result-20260826.json)
passes every declared gate. Reconstructing the full phase-resolved electric
field from electron density, ion density, and endpoint potentials gives
relative RMS errors of about `2.91e-10` for both AuroraPIC members and
`1.95--1.97e-7` for the three native members. The native error is consistent
with its six-digit text output. This provides unusually direct evidence that
the density deposition, potential, and electric-field diagnostics obey the
same discrete field equation in both production calculations.

The critical-window boundary-drive mean-square-field ratio is `0.99727` for
both AuroraPIC members, excluding the prescribed voltage as the mature field
gap. Under the preregistered Shapley rule, electron-space-charge fractions are
`69.17` and `72.23`, ion-space-charge fractions are `-68.18` and `-71.24`, and
boundary fractions are approximately `0.008`; the allocations close to the
field-energy gap within about `1.1e-7 V^2/m^2`. The formal declared outcome is
therefore `electron_space_charge_dominant`.

That label requires an important qualification. The absolute attribution sums
are `137.36` and `143.47` times the small net gap, a conditioning diagnostic
added transparently after the locked outcome. Electron and ion space-charge
fields almost cancel, as expected in a quasineutral plasma; replacing either
species independently creates a large, non-evolved charge imbalance. Thus the
formal result must not be presented as an independent electron-causation
claim. The robust findings are narrower: the operator closes, the external
drive is at parity, and the mature field difference resides in the net
space-charge distribution. The next discriminator should use net-charge and
sheath-edge modes directly and then seek their earliest temporal divergence.

### Net-charge and sheath-structure discriminator

The next
[`prospective rule`](../benchmarks/ccp/edupic-net-charge-sheath-rule-20260826.json)
therefore avoids separate electron/ion substitutions. It locks the bounded
charge separation `(ni-ne)/(ni+ne)`, the direct net-space-charge field,
drive--space-charge cancellation, and left-sheath density-ratio crossings at
`ne/ni = 0.8`, `0.9`, and `0.95`. All three thresholds are reported because a
single density-ratio crossing is not a unique physical sheath edge.

The
[`result`](../benchmarks/ccp/edupic-net-charge-sheath-result-20260826.json)
passes every hash, shape, crossing, positivity, and repeatability gate, but
none of the four preregistered association rules passes. The two AuroraPIC
critical charge-separation-RMS ratios are `0.7945/0.7981`, opposite the
declared larger-separation hypothesis. Drive--space-charge cancellation ratios
are `1.0883/1.0839`, below the locked `1.10` support threshold. Sheath-width
ratios span `0.9607--0.9719`, rather than indicating a wider sheath, while
positive-sheath-charge ratios span `1.0206--1.0277`, rather than indicating
missing positive charge. AuroraPIC relative ranges are below `0.46%` and
native ranges below `2.01%` across all metrics.

This disciplined null result excludes a large width, integrated positive
charge, or bounded charge-separation-amplitude explanation under the declared
criteria. The modestly stronger cancellation is real and repeatable but does
not cross the prospective support threshold. Combined with the exact Poisson
closure, the remaining field gap points to the detailed phase and spatial
organization of net charge rather than a gross sheath scalar. Establishing
onset now requires matched earlier phase-resolved charge/field histories; the
mature averages cannot reveal which evolution step diverged first.

### Four-cycle phase-snapshot onset test

The subsequent
[`prospective onset rule`](../benchmarks/ccp/edupic-phase-snapshot-onset-rule-20260826.json)
locks ten instantaneous grid snapshots per RF cycle, three native seeds, two
AuroraPIC microstates, and a four-cycle window. Passive native instrumentation
is confirmed by byte-identical final checkpoints and all eight pre-existing
diagnostics. Every source, executable, snapshot-set, shape, coordinate,
resource, and passivity check in the
[`result`](../benchmarks/ccp/edupic-phase-snapshot-onset-result-20260826.json)
passes.

At the primary phase `0.3`, each AuroraPIC member remains below the native
ensemble in every cycle. The eight field-energy ratios span `0.6581--0.7946`;
the deficit is therefore already present in cycle one and does not emerge
during the locked window. Four-cycle relative ranges at this phase are below
`0.161` for every individual trajectory, and the cross-member field
repeatability limits pass at phase `0.3` in both codes.

The preregistered *joint* repeatability gate nevertheless fails. AuroraPIC has
five failing cycle/phase/metric cells (maximum relative range `0.2582`), and
the native ensemble has three, all in the phase-`0.2` field metric (maximum
`0.5049`). The four-cycle phase-neighborhood rule also fails because one
AuroraPIC member's phase-`0.2` ratio is `1.0150`, despite phase-`0.3` and
phase-`0.4` deficits. The formal outcome is consequently
`inconclusive_failed_joint_gate`. The defensible observation is narrower: the
previously localized phase-`0.3` field deficit is stable, repeatable at that
phase, and predates this measurement window. This test cannot identify its
first divergence from a common state or establish published-benchmark or
experimental validity.

### Collision-free common-state divergence

APS v3 adds an explicitly signed `leapfrog_half_step` interchange contract for
structured electrostatic 1D states. A real 222,591-particle eduPIC checkpoint
was converted and independently loaded with matching semantic signature while
using less than 27 MiB RSS. The
[`prospective common-state rule`](../benchmarks/ccp/edupic-common-state-divergence-rule-20260826.json)
then disables collisions in both codes and compares logarithmically spaced
charge and field profiles from the same mature particle coordinates and
half-step velocities. The checksum-bound
[`result`](../benchmarks/ccp/edupic-common-state-divergence-result-20260826.json)
and [`execution record`](../benchmarks/ccp/edupic-common-state-divergence-execution-20260826.json)
preserve the reductions, passive-instrumentation check, and resource envelope.

The initial profiles agree extremely closely: charge relative RMS is
`2.63e-14`, field relative RMS is `2.70e-9`, and the critical regional
field-energy ratio is `1.00000015`. This is direct dynamic evidence that the
state conversion, macro-weight normalization, charge deposition, boundary
drive, Poisson solve, endpoint Gauss correction, and field reduction agree at
the common initial state.

After one step, full-domain field relative RMS remains only `5.37e-4`, but the
critical regional field-energy ratio becomes `0.96687` and stays outside the
locked `[0.98,1.02]` band for the next samples. The preregistered outcome is
therefore `one_step_mover_or_boundary_mismatch`. Electron and ion populations
are exactly equal through horizon 100; the first four-particle electron-loss
difference occurs only at horizon 200. The field discrepancy consequently
precedes the boundary-loss bifurcation.

Source ordering provides a specific post-hoc mechanism candidate. eduPIC
deposits ion density immediately before each 20-electron-step ion push, moves
the ions, and then reuses the pre-push density until the next ion update.
AuroraPIC moves due ions and redeposits all live species at their new positions
for the following field solve. That difference acts immediately after the
first common-state ion move and is consistent with the observed localized
departure. It is not yet a causal result: the next prospectively locked control
must make AuroraPIC retain the pre-push ion density and test whether the trace
collapses toward eduPIC.

That causal test is now prospectively fixed in the
[`held-density control rule`](../benchmarks/ccp/edupic-ion-density-refresh-control-rule-20260826.json).
It locks the cache time level, boundary-removal behavior, restart contract,
early-horizon improvement thresholds, resource envelope, and claim boundary
before the control implementation or any control output exists.

The completed [`control result`](../benchmarks/ccp/edupic-ion-density-refresh-control-result-20260826.json)
passes every integrity and initial-parity gate and gives
`strong_ion_density_refresh_mechanism_support`. At horizons 1, 2, and 5,
the critical-field log-error is only `9.4e-6`, `1.4e-5`, and `2.9e-5` of
the corresponding baseline error, while full-field relative RMS falls from
approximately `5e-4` to `5e-9--2e-8`. Populations match native at every
locked early horizon. The stricter through-horizon-20 explanation gate does
not pass: the regional field-energy ratio becomes `1.03438` precisely at the
next ion-cache refresh boundary. Thus the pre-push density time level is a
causally supported explanation for the immediate discrepancy, while a
one-base-step refresh-cadence alignment remains to be prospectively tested.
The checksum and resource envelope are preserved in the
[`execution record`](../benchmarks/ccp/edupic-ion-density-refresh-control-execution-20260826.json).
The follow-up cadence hypothesis and its non-retroactive thresholds are fixed
in the [`upcoming-due refresh rule`](../benchmarks/ccp/edupic-ion-density-refresh-cadence-rule-20260826.json)
before implementation or execution.

The prospective [`cadence result`](../benchmarks/ccp/edupic-ion-density-refresh-cadence-result-20260826.json)
passes every gate and yields
`full_collision_free_common_state_trace_closure`. Electron and ion
populations match eduPIC exactly at all 15 horizons. No material field flag
occurs through horizon 3999. At the final horizon, full-field relative RMS is
`0.00551` and the critical regional field-energy ratio is `0.999477`; at
horizon 20 those values are `5.78e-8` and `1.00000356`. This establishes that
the prior deterministic collision-free difference came from the held-density
time level and its refresh cadence, within this common-state experiment. The
[`execution record`](../benchmarks/ccp/edupic-ion-density-refresh-cadence-execution-20260826.json)
preserves the binary, runner, report, reductions, and resource hashes. It does
not by itself validate collisions or establish agreement with Turner's
published steady-state ensemble.

### Collision-enabled common-state ensemble

The next prospective layer starts both implementations from that same mature
state, enables the pinned eduPIC argon electron and ion collision laws, and
runs five independent RNG members per code to the aligned one-period endpoint.
The generators differ (`mt19937` versus `mt19937_64`), so the
[`rule`](../benchmarks/ccp/edupic-collision-enabled-common-state-ensemble-rule-20260826.json)
compares ensemble statistics rather than falsely pairing seed trajectories.
The native instrumentation adds no random draws and only counts collision
branches already selected by eduPIC. The checksum-bound
[`execution record`](../benchmarks/ccp/edupic-collision-enabled-common-state-ensemble-execution-20260826.json)
contains all ten endpoint observations, binary hashes, and resource results.

The prospective
[`result`](../benchmarks/ccp/edupic-collision-enabled-common-state-ensemble-result-20260826.json)
passes every integrity and physics gate, yielding
`one_period_collision_enabled_stochastic_consistency_supported`. Symmetric
relative differences between mean accepted counts are `0.000288` for electron
elastic, `0.00330` for excitation, `0.00207` for ionization, `0.0193` for ion
isotropic, and `0.0108` for ion backward collisions. All four electrode-loss
means pass their 10% gates. Mean final populations differ by 0.6 electron and
7.4 ions out of initial populations above 108,000, with matching net-change
signs. The ensemble-mean endpoint electric-field profiles differ by `0.00621`
relative RMS, and the AuroraPIC-to-eduPIC mean field-energy proxy ratio is
`0.999083`.

This is an end-to-end collision-enabled dynamic cross-code result: it jointly
exercises collision selection and products, ion creation, particle absorption,
subcycled charge deposition, Poisson solution, and the mover for one RF
period. It materially strengthens numerical and implementation credibility.
It remains a bounded pilot for one mature 1D argon state—not a steady-state
Turner benchmark pass, experimental validation, convergence proof, or evidence
for other gases and dimensions.

The prospectively locked four-period extension then repeated three seeds per
code through 15,999 pushes. Its
[`execution record`](../benchmarks/ccp/edupic-collision-enabled-common-state-four-period-execution-20260826.json)
completed all six members serially in 954.15 member-seconds. Peak RSS was
38.5 MiB for instrumented eduPIC and 188.7 MiB for AuroraPIC, below the frozen
256 MiB reporting gate.

The
[`four-period result`](../benchmarks/ccp/edupic-collision-enabled-common-state-four-period-result-20260826.json)
passes the collision, wall-loss, ion-population, and field gates. Relative
ensemble-mean differences are `0.00230` for elastic collisions, `0.00955` for
excitation, `0.00119` for ionization, `0.00985` for ion-isotropic collisions,
and `0.00223` for ion-backward collisions. Wall-loss differences range from
`0.00102` to `0.0151`. The final mean-field relative RMS is `0.00860`, and
the field-energy ratio is `0.999834`.

The frozen electron-population sign clause alone fails. eduPIC's mean change
is -3.67 electrons and AuroraPIC's is +9.67 electrons from an initial 108,586;
the mean endpoints differ by only 13.33 particles (`0.0123%`). Those changes
are smaller than the member scatter (endpoint standard deviations 37.1 and
12.7 particles), while the separately observed ionization and both electron
wall-loss components all pass. The formal outcome must therefore remain
`localized_collision_enabled_common_state_discrepancy`; tolerances are not
retuned after observation. Physically, this localizes the unresolved question
to a near-zero residual of stochastic electron creation and loss, not to a
material collision-rate, field, or density disagreement. A higher-power
prospective balance test is required to distinguish a real small bias from
sampling noise.

That higher-power follow-up is now complete. The corrected
[`five-member rule`](../benchmarks/ccp/edupic-collision-enabled-common-state-five-member-balance-rule-20260826.json)
reuses the three parent members only after exact endpoint and field-hash
verification and adds seeds 51949 and 63059 in each code. An initial execution
attempt exposed a JSON-generation defect that rounded the unsigned 64-bit APS
signature. AuroraPIC rejected the state before initialization; the attempt,
cause, quarantine, and unchanged-threshold disposition are preserved in the
[`preflight failure record`](../benchmarks/ccp/edupic-collision-enabled-common-state-five-member-balance-preflight-failure-20260826.json).
The exact signature was restored and the execution was relocked before all four
new members were rerun.

The checksum-bound
[`five-member result`](../benchmarks/ccp/edupic-collision-enabled-common-state-five-member-balance-result-20260826.json)
yields `four_period_electron_balance_equivalence_supported`. The AuroraPIC
minus eduPIC final-electron mean difference is 36.6 particles. Its conservative
two-sided 90% interval, using the preregistered t(4) critical value, is
`[-7.24, 80.44]` particles and lies wholly inside the practical-equivalence
band `[-108.586, 108.586]`, or +/-0.1% of the initial population. The combined
five-member collision, wall-loss, ion-population, field-profile, and
field-energy gates also pass. The final field-profile relative RMS is
`0.00722`, and the field-energy ratio is `1.000148`. The
[`execution record`](../benchmarks/ccp/edupic-collision-enabled-common-state-five-member-balance-execution-20260826.json)
also proves that all six reused parent members matched their locked endpoints
and fields.

This prospective equivalence result resolves the earlier sign-only failure at
the declared practical scale without altering or erasing that original formal
outcome. It supports four-period integrated stochastic consistency for this
state; steady-state published-profile validation remains the next distinct
scientific claim.

### Corrected-cadence phase-resolved EEDF closure

An older four-cycle argon comparison found AuroraPIC/native-eduPIC ratios of
`0.85296` for the EEDF-folded ionization frequency and `0.87182` for the
electron fraction above 15.8 eV in the critical interior phase window. That
run predated the held-ion-density cadence correction and began from an
independently evolved AuroraPIC trajectory. It therefore could not distinguish
a collision/heating defect from the already identified field-staggering error.

The new
[`prospective rule`](../benchmarks/ccp/edupic-corrected-cadence-phase-eedf-rule-20260826.json)
exports the exact native cycle-2620 checkpoint as a half-step APS state and
runs three corrected AuroraPIC continuations for the same four-cycle window.
Each member samples 200 RF phase bins, seven spatial regions, and 320 energy
bins every two electron steps. The
[`execution record`](../benchmarks/ccp/edupic-corrected-cadence-phase-eedf-execution-20260826.json)
shows more than 1.03 million macro-observations in every critical
region/phase bin, zero histogram overflow, approximately 206 MiB peak RSS,
and all three members completed within the frozen resource envelope.

The checksum-bound
[`result`](../benchmarks/ccp/edupic-corrected-cadence-phase-eedf-result-20260826.json)
yields `strong_corrected_cadence_phase_eedf_closure`. The corrected critical
folded-ionization ratio is `1.01714`, the above-15.8-eV tail ratio is
`0.995780`, and the three critical phase-slice folded-ionization ratios are
`1.04277`, `1.00525`, and `1.01737`. Critical histogram total-variation
distance is only `0.001764`, compared with approximately `0.0482` previously.
Candidate folded-ionization and tail relative ranges are `0.0297` and
`0.0178`, so the closure is repeatable across the three independent random
streams.

This resolves the old regional energetic-tail deficit as overwhelmingly a
state/charge-refresh staggering artifact rather than a missing AuroraPIC
collision-heating mechanism. Together with collision traffic, population,
wall-loss, and field closure, it substantially strengthens the numerical and
physical implementation case for this 1D argon CCP. It does not turn the
separate helium Turner density discrepancy into a pass; that published-profile
claim remains a separate unresolved result.

### Turner subcycle-policy invariance control

Inspection of the exact Turner Case 1 deck shows why the argon cadence fix
cannot be transferred to the helium result: neither Turner species declares a
`timestep_multiplier`, so both advance on every electron step. The
[`prospective invariance rule`](../benchmarks/ccp/turner-case1-subcycle-policy-invariance-rule-20260826.json)
therefore continued one identical late Turner checkpoint for one RF cycle
under `current_position` and `pre_push_held`, changing only that policy and the
output directory.

The checksum-bound
[`result`](../benchmarks/ccp/turner-case1-subcycle-policy-invariance-result-20260826.json)
classifies the control as `turner_subcycle_policy_invariance_established`.
Final fields, scalar history, collision counters, boundary losses, power
transfer, and the one-cycle spatial average are byte-identical. Both serial
branches completed in approximately 3.8 seconds with about 12 MiB peak RSS.
This formally closes the cadence ambiguity without wasting three new
512,000-step runs: the existing three-seed Turner mean density bias of
`+2.237%` is unchanged and remains unresolved.

The control also exposed an overly strict restart guard. Checkpoints predating
the held-charge-cache record were rejected under `pre_push_held` even when all
species had unit cadence and no cache was required. The loader now permits
that physically equivalent case while continuing to reject legacy restarts
with any genuinely subcycled species; both paths have regression coverage.
The next Turner discriminator must target helium collision/sheath-transport
conventions or an independent matched implementation, not the subcycle charge
policy.

## Bounded execution ladder

Case 1 remains the smallest whole-discharge target, but shortening it changes
the statistical benchmark. Work therefore proceeds without relabeling a
reduced run as a pass:

1. **C0, collision-law verification:** exact center-of-mass energy lookup,
   finite-mass isotropic recoil, backward velocity exchange, electron
   equal-sharing ionization, and restart fingerprints. This stage is complete.
2. **C1, source and diagnostic lock:** the electronic supplement is acquired
   outside the repository and its archive/member hashes, structure, and usage
   constraint are pinned. Exact tables are locally normalized without
   resampling and load through provenance-bearing gas manifests. Restart-safe
   whole-cycle spatial-density averaging and the statistical comparator are
   complete. Restart-safe source, wall-current, and species-power observables
   now cover the published global Case 1 quantities.
3. **C2, bounded startup screening:** run the exact Case 1 grid, timestep,
   particle weight, waveform, and collision model for a small declared number
   of RF periods. Check invariants, collision balance, sheath formation,
   restart, and resource behavior. This cannot use the steady-state
   chi-squared acceptance range. The exact-population, few-step runtime
   qualification and its hard resource gates are complete. A checkpoint-split
   whole-RF-period startup now proves restart continuity, species/collision
   balance, waveform timing, finite diagnostics, and early boundary-field
   formation. The bounded four-cycle horizon adds phase-matched population,
   collision, and field trends; longer blockwise stationarity remains.
4. **C3, checkpointed steady-state comparison:** continue only through
   explicit low-priority blocks, establish whole-cycle stationarity, average
   the final 32 periods, and apply the published chi-squared test.

This ladder gives useful bounded evidence early while preserving one honest
definition of a Turner benchmark pass.
