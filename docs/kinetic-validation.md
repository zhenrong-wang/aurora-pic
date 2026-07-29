# Kinetic verification

AuroraPIC includes a deterministic, quantitative linear Landau-damping
benchmark. It is the first system-level physics verification case: it checks
that particle loading, charge deposition, the periodic Poisson solve,
field interpolation, leapfrog advancement, and diagnostics reproduce a known
collisionless kinetic response together. This is stronger than a structural
smoke test, but it is not experimental validation of a device or gas model.

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
  build/aurorapic_cli --report build/landau-report.json
```

The validator uses only the Python standard library, emits a machine-readable
JSON report, and removes its temporary state and 241 small field snapshots on
success. Pass `--keep-output` to retain them after a failure or for inspection.
Its bounded envelope is 65,536 total particles, 64 cells, 240 steps, one
AuroraPIC process, and the caller's OpenMP thread limit. The repository's
`scripts/verify.sh` fixes that limit to one by default.

## What remains

Passing this case verifies linear collisionless electrostatic kinetics in the
current 1D path. It does not validate collision cross sections, material
boundaries, imported geometry, multidimensional mode propagation, steady-state
convergence, electromagnetic fields, or a real thruster/discharge.

The next verification and validation ladder is:

1. add a quantitative two-stream growth-rate and saturation benchmark;
2. add 2D and 3D Langmuir-mode dispersion benchmarks;
3. reproduce the published Turner helium capacitively coupled plasma
   benchmark with an authoritative open collision dataset;
4. add an imported-geometry probe current-voltage comparison;
5. target the LANDMARK Hall-thruster benchmark after the required
   magnetic-field, source, collision, and parallel-runtime capabilities exist.

