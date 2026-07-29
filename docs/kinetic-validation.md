# Kinetic verification

AuroraPIC includes deterministic, quantitative linear Landau-damping,
two-stream-instability, and orthogonal 2D and 3D Langmuir benchmarks. These
system-level physics cases check that particle loading, charge deposition, the
periodic Poisson solve, field interpolation, leapfrog advancement, and
diagnostics reproduce known collisionless kinetic responses together. This is
stronger than a structural smoke test, but it is not experimental validation
of a device or gas model.

## Linear Landau-damping case

The benchmark uses the normalized 1D1V Vlasov-Poisson problem on
`0 <= x < 4*pi`:

```text
f_e(x, v, 0) = (1 + 0.01*cos(0.5*x))*exp(-v^2/2)/sqrt(2*pi)
```

Electrons have unit negative charge, unit mass, and unit thermal velocity.
An equal-weight, spatially uniform ion population with mass `1e6` provides an
effectively stationary neutralizing background. The generated state contains
32,768 particles of each species. Electron positions are obtained by
inverting the perturbed spatial cumulative distribution, and normal velocity
quantiles are assigned with a bit-reversal permutation. This quiet,
deterministic loading avoids a committed multi-megabyte particle artifact.

The run uses 64 periodic cells, `dt = 0.05`, and 240 steps. For `k = 0.5`,
the linear Vlasov-Poisson dispersion relation gives electric-field damping
rate `gamma = -0.1533` and angular frequency `omega = 1.4156`. These reference
values are tabulated for the standard linear case by Crouseilles,
Mehrenberger, and Vecil,
[ESAIM: Proceedings 32 (2011), 211-230](https://doi.org/10.1051/proc/2011022).
The independent
[PICLas Landau-damping tutorial](https://piclas.readthedocs.io/en/stable/userguide/tutorials/pic-poisson-Landau-damping/pic-poisson-Landau-damping.html)
uses the same periodic domain and wavenumber and documents the `-0.153`
electric-field envelope.

The validator extracts the magnitude of the `k = 0.5` electric-field Fourier
mode from every field snapshot. It fits the logarithm of four or more
successive envelope peaks between normalized times 1 and 10.5, and obtains
the angular frequency from their half-period spacing.

The committed acceptance envelope is deliberately wider than one platform's
last digits:

| Quantity | Reference | Acceptance |
| --- | ---: | ---: |
| Initial electric-field mode amplitude | `alpha/k = 0.020` | `0.019` to `0.021` |
| Electric-field damping rate | `-0.1533` | `-0.18` to `-0.13` |
| Angular frequency | `1.4156` | `1.34` to `1.48` |
| Maximum relative total-energy drift | `0` | at most `1e-3` |

Run it after building:

```sh
OMP_NUM_THREADS=1 python3 scripts/validate_kinetic_benchmarks.py \
  build/aurorapic_cli --benchmark landau \
  --report build/landau-report.json
```

The validator uses only the Python standard library, emits a machine-readable
JSON report, and removes its temporary state and 241 small field snapshots on
success. Pass `--keep-output` to retain them after a failure or for inspection.
The Landau case's bounded envelope is 65,536 total particles, 64 cells, 240
steps, one AuroraPIC process, and the caller's OpenMP thread limit. The
repository's `scripts/verify.sh` fixes that limit to one by default.

## Linear two-stream instability and nonlinear turnover

The second benchmark uses the normalized symmetric warm-beam distribution

```text
f_e(x, v, 0) =
  (1 + 0.001*cos(0.2*x))
  * (exp(-(v-2.4)^2/2) + exp(-(v+2.4)^2/2))
  / (2*sqrt(2*pi))
```

on `0 <= x < 10*pi`. A spatially uniform ion population of mass `1e6`
provides the effectively stationary neutralizing background. Each electron
beam contains 16,384 particles and the ion background contains 32,768, for
65,536 particles total. Beam particles at each quiet-start position have
exactly opposite velocities, and bit-reversed normal quantiles represent the
unit thermal spread without random sampling noise.

For this standard case, the linear Vlasov-Poisson electric-field growth rate
is `lambda = 0.2258`. It is tabulated by Roberts et al.,
[Computers & Mathematics with Applications 154 (2024), 103-119](https://doi.org/10.1016/j.camwa.2023.11.014).
The setup also matches the independently documented symmetric warm-beam
benchmark in the
[Sandia DPG report](https://www.osti.gov/servlets/purl/1891588).

AuroraPIC runs 128 periodic cells with `dt = 0.05` through normalized time
50. The validator fits the logarithm of the `k = 0.2` electric-field mode
over `14 <= t <= 28`, after the startup transient and before nonlinear
saturation. It then locates the global mode-amplitude peak and requires a
subsequent reduction, demonstrating turnover into the nonlinear saturated
regime. The peak amplitude is a numerical regression observable, not a
universal theoretical constant.

| Quantity | Reference | Acceptance |
| --- | ---: | ---: |
| Initial electric-field mode amplitude | `alpha/k = 0.005` | `0.0048` to `0.0052` |
| Electric-field growth rate | `0.2258` | `0.19` to `0.26` |
| Linear-fit coefficient of determination | `1` | at least `0.97` |
| Nonlinear peak time | numerical regression | `33` to `43` |
| Peak/initial mode amplification | numerical regression | `80` to `180` |
| Minimum post-peak/peak amplitude | nonlinear turnover | at most `0.95` |
| Maximum relative total-energy drift | `0` | at most `1e-3` |

Run only this case with:

```sh
OMP_NUM_THREADS=1 python3 scripts/validate_kinetic_benchmarks.py \
  build/aurorapic_cli --benchmark two-stream \
  --report build/two-stream-report.json
```

## Orthogonal 2D Langmuir oscillations

The first multidimensional quantitative case initializes a cold electron
density perturbation

```text
n_e(s, 0) = 1 + 0.01*cos(s),  s = x or y
```

in a periodic `2*pi` by `2*pi` domain. Uniform ions with mass `1e6` provide an
effectively stationary background. In normalized units, `q_e = -1`,
`m_e = 1`, `n_e = 1`, and `epsilon = 1`, so the cold electron plasma
frequency is

```text
omega_p = sqrt(n_e*q_e^2/(epsilon*m_e)) = 1.
```

This is the standard plasma-oscillation relation documented by the
[PICLas plasma-wave tutorial](https://piclas.readthedocs.io/en/latest/userguide/tutorials/pic-poisson-plasma-wave/pic-poisson-plasma-wave.html).
WarpX likewise publishes analytic
[1D, 2D, and 3D Langmuir-wave examples](https://warpx.readthedocs.io/en/latest/usage/examples/langmuir/README.html).

AuroraPIC runs separate x- and y-directed modes on the same 32 by 32 mesh.
Each run uses a deterministic 64 by 64 particle lattice for 4,096 electrons
and 4,096 ions, `dt = 0.05`, and 320 steps. Only scalar diagnostics are
written. The validator obtains the mode frequency from five half-period
field-energy peaks, tests amplitude retention, and compares the orthogonal
results. For `alpha = 0.01` and `k = 1`, the analytic initial field
amplitude is `alpha/k = 0.01`; integrating its squared sinusoid over the
domain gives initial field energy `pi^2*1e-4`.

| Quantity | Reference | Acceptance |
| --- | ---: | ---: |
| x-directed angular frequency | `1` | `0.97` to `1.03` |
| y-directed angular frequency | `1` | `0.97` to `1.03` |
| Relative x/y frequency difference | `0` | at most `0.005` |
| Relative x/y initial-field difference | `0` | at most `0.001` |
| Initial field energy, each direction | `pi^2*1e-4` | `0.00095` to `0.00102` |
| Last/first field-amplitude peak, each direction | `1` | `0.98` to `1.02` |
| Maximum relative total-energy drift, each direction | `0` | at most `0.003` |

Run only the multidimensional case with:

```sh
OMP_NUM_THREADS=1 python3 scripts/validate_kinetic_benchmarks.py \
  build/aurorapic_cli --benchmark langmuir-2d \
  --report build/langmuir-2d-report.json
```

## Orthogonal 3D Langmuir oscillations

The structured 3D3V case extends the same cold, normalized Langmuir problem
to separate x-, y-, and z-directed perturbations in a periodic
`(2*pi)^3` domain. It therefore exercises trilinear charge deposition, the
3D periodic Poisson solve, vector-field interpolation, 3D particle
advancement, and scalar diagnostics along every coordinate axis.

Each direction uses a 16 by 16 by 16 mesh and deterministic particle lattice
with 4,096 electrons plus 4,096 effectively stationary ions. The timestep is
`0.05`, and each direction runs for 320 steps with field and particle files
disabled. The deliberately bounded mesh keeps the full default regression
safe on developer machines. For a sinusoidal field of amplitude
`alpha/k = 0.01`, the analytic initial field energy is
`alpha^2*(2*pi)^3/4`.

The acceptance envelope records the accuracy expected from this coarse
regression rather than claiming a converged 3D solution. In particular, its
2% energy-drift gate is looser than the 2D case and must not be reused as a
production-study convergence criterion.

| Quantity | Reference | Acceptance |
| --- | ---: | ---: |
| Angular frequency, each direction | `1` | `0.97` to `1.03` |
| Maximum relative directional frequency spread | `0` | at most `1e-6` |
| Maximum relative directional initial-field spread | `0` | at most `1e-6` |
| Initial field energy, each direction | `alpha^2*(2*pi)^3/4` | `0.0059` to `0.0063` |
| Last/first field-amplitude peak, each direction | `1` | `0.98` to `1.02` |
| Maximum relative total-energy drift, each direction | `0` | at most `0.02` |

Run only this case with:

```sh
OMP_NUM_THREADS=1 python3 scripts/validate_kinetic_benchmarks.py \
  build/aurorapic_cli --benchmark langmuir-3d \
  --report build/langmuir-3d-report.json
```

The default invocation runs all four kinetic benchmarks sequentially. Use
`--benchmark landau`, `--benchmark two-stream`, `--benchmark langmuir-2d`,
or `--benchmark langmuir-3d` to select one.

## What remains

Passing these cases verifies damped and unstable collisionless electrostatic
kinetics in the current 1D path, including nonlinear two-stream turnover, and
cold directional plasma oscillations in structured 2D3V and 3D3V. It does
not validate collision cross sections, material boundaries, imported
geometry, warm multidimensional dispersion or damping, steady-state
convergence, electromagnetic fields, or a real thruster/discharge.

The next verification and validation ladder is:

1. reproduce the published Turner helium capacitively coupled plasma
   benchmark under the staged readiness contract in
   [ccp-validation.md](ccp-validation.md);
2. add an imported-geometry probe current-voltage comparison;
3. follow the pinned LANDMARK and experimental Hall-thruster ladder in
   [hall-thruster-validation.md](hall-thruster-validation.md), beginning with
   analytic magnetized-particle and spatial prescribed-field verification.
