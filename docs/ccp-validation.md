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

The primary Turner observable is now produced by a generic 1D post-step
spatial averager. It deposits each species' represented particle number with
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

This selects 12,800 post-step samples. Configuration validation requires an
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
Checkpoint v9 preserves every phase-bin accumulator and sample count; bounded
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
