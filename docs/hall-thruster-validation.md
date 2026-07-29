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
| Mixed topology | Structured periodic or Dirichlet modes exist, but the complete LANDMARK side contract is not implemented | Support periodic azimuth plus independent axial/radial field and particle boundaries |
| Volumetric pair source | Collision ionization and tagged-boundary sources exist | Add deterministic-rate spatial pair profiles with equal-position creation |
| Cathode/current control | Not available | Implement the exact emission-plane current balance and potential correction |
| Radial benchmark virtual axis | Not available | Implement bounded virtual-axis replacement and audit its energy/particle flux |
| HET diagnostics | Scalar and field snapshots exist | Add transverse/time averages, species moments, current components, Fourier spectra, and ensemble statistics |
| Xenon material data | No authoritative bundled package | Keep LANDMARK collisionless; separately provenance and validate Xe collision/wall data for real devices |
| Scale-out runtime | Serial/OpenMP only | Add MPI domain decomposition before production-size LANDMARK runs |
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

## Resource policy

No published-scale HET case runs in `scripts/verify.sh`, CTest, or ordinary CI.
H0 and H1 tests must default to one process and one thread and have explicit
particle, cell, step, memory, and output limits. H2-H5 are opt-in production
campaigns with a preflight resource estimate, explicit output directory,
checkpoint cadence, storage quota, and wall-time limit. Reference datasets
remain external and checksum-pinned; they are never vendored into the Git
repository.

The first H0 field slice is complete: a shared strict `coordinate Bx By Bz`
profile supports linear interpolation and full-domain coverage checks across
structured 2D/3D and imported 2D, while preserving uniform-field
compatibility. The committed LANDMARK magnetic profile is exercised only as a
four-step integration smoke and cannot make a discharge claim. The next
bounded slice is mixed periodic/nonperiodic structured topology plus analytic
periodic wrapping and absorbing axial-flux verification. That capability is
generic to crossed-field devices and open plasma domains; it is not embedded
Hall-specific solver logic.
