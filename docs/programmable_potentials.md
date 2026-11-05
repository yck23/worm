# Programmable site potentials and numerical speckle

The Bose-Hubbard solver accepts an arbitrary **static** real potential at every
site. Compile with `-DPROGRAMMABLE_POTENTIAL=ON` (which selects `UNIFORM=OFF`)
and set `onsite_potential_file` in the INI or on the command line. Positive
potential repels particles:

```text
H = -t_hop sum_<ij> (b_i^dagger b_j + h.c.)
    + U_on/2 sum_i n_i(n_i-1) + V_nn sum_<ij> n_i*n_j
    + sum_i V_i*n_i - mu sum_i n_i
```

The existing scalar couplings remain configurable. `V_nn` is the nearest-neighbor
interaction, distinct from the site potential `V_i`. `Total_Energy` includes
`-mu*N`, as before; add `mu*<N>` to obtain the energy without that term.
`Potential_Energy` is the full diagonal energy, not just the imposed disorder.
All energies, including the file values, use the same units as `t_hop`; `beta`
uses inverse energy units. There is no automatic conversion from Hz or metres.

## Choose particle number and spatial potential independently

**Yes: arbitrary site potentials and speckle work with fixed total particle
number.** Two independent choices control a fresh run:

| Choice | INI setting | Windows runner selection |
| --- | --- | --- |
| Fixed total N (canonical) | `canonical = N`, e.g. `1024` | Add `-Canonical` |
| Fluctuating N (grand canonical) | `canonical = -1` and choose `mu` | Omit `-Canonical` |
| Use a potential map | `onsite_potential_enabled = 1` and `onsite_potential_file = "path/map.dat"` | Add `-ProgrammablePotential`, or pass `-PotentialFile path/map.dat` |
| Turn the added site potential off | `onsite_potential_enabled = 0` | Or add `-DisablePotential` to override the INI |

`-Canonical` selects the executable; it does **not** set the numerical value of
N. The `.ini` file sets N. The `.dat` file contains only the potential values;
it does not select the ensemble, lattice size, or executable. With the runner,
prefer `-PotentialFile`: its path is resolved from your current directory.
An INI-only map path must be absolute or relative to `-OutputDirectory`
(`results` by default).

"Uniform" describes the values across space. "Static" means unchanged during
the run. All supported maps are static, including spatially nonuniform speckle.
The potentials do not refresh on each sweep, and MPI replicas use the same map.

| Desired added site potential | What to put in the `.dat` file |
| --- | --- |
| Uniform offset `V0` everywhere | Exactly `Lx*Ly` copies of `V0` |
| Spatially nonuniform potential | One value per site, varying across sites |
| Zero added site potential | All zeros, or disable the potential without a file |

One programmable executable handles all these choices; changing the map values
requires no rebuild. `UNIFORM=OFF` is the implementation that permits different
site values, and it also accepts a constant or zero map. The separately compiled
`UNIFORM=ON` executable does not load an enabled map.

## Build and run on Windows

Run from the checkout containing this guide. For the current workspace:

```powershell
Set-Location '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
```

The executables are already built here. Skip the first command below unless
you need to rebuild; generate a map and then start the example:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/build_windows.ps1 -ProgrammablePotential -Jobs 32
python scripts/generate_speckle.py results/speckle.dat --nx 32 --ny 32 --dx 1 --correlation-length 2 --mean-strength 0.5 --seed 42 --metadata-json results/speckle.json
powershell -ExecutionPolicy Bypass -File scripts/run_bose_square.ps1 -ProgrammablePotential -PotentialFile results/speckle.dat -ParameterFile parameter_files/BoseHubbard_speckle.ini -Lx 32 -Ly 32 -Sweeps 10000 -Thermalization 1000 -RuntimeLimit 60 -JobName speckle
```

The separate package is `dist/windows-square-programmable`. It contains ordinary
and canonical serial/MPI executables, just like the original package. To reuse
an already-built compatible dependency, pass `-AlpsInstall <prefix>` to the build
script; see `WINDOWS_BUILD.md`. The runner resolves `-PotentialFile` before
changing to its output directory. For a direct executable invocation, relative
potential paths are relative to the process's working directory, not the INI.
Use an absolute path when in doubt. The generic runner overrides dimensions and
run lengths with its switches, so pass dimensions matching the generated map.
This supplied run is an integration example, not a converged physics result.

For other platforms using an installed ALPSCore dependency:

```sh
cmake -S src -B build/programmable -DLATTICE=square -DPROGRAMMABLE_POTENTIAL=ON -DCMAKE_PREFIX_PATH=/path/to/alps-install
cmake --build build/programmable --parallel
./build/programmable/qmc_worm parameter_files/BoseHubbard_speckle.ini
```

Add `-DCWINDOW=ON` in a separate build and set `canonical=N` for fixed particle
number. A uniform executable rejects an enabled potential file instead of ignoring it.
Site-potential input currently applies to Bose-Hubbard, not the XXZ field.

## Fixed-N example: speckle, uniform, or off

The supplied [canonical INI](../parameter_files/BoseHubbard_speckle_canonical.ini)
uses a 32 by 32 lattice and exactly 1024 total particles (mean filling 1):

```ini
canonical = 1024
canonical_window = 0.5
onsite_potential_enabled = 1
onsite_potential_file = "results/speckle.dat"
```

Edit `canonical` to choose your total N; it must be an integer from 0 to
`Lx*Ly*nmax`. `initial_particle_number` and `initial_occupancy` do not choose N
in canonical mode. Keep `0 < canonical_window < 1`: this permits the auxiliary
open-worm excursion while every measured closed configuration has exactly N.
The site occupations can redistribute in response to the potential.

### Run with nonuniform speckle

From this repository root, generate and run:

```powershell
python scripts/generate_speckle.py results/speckle.dat --nx 32 --ny 32 --dx 1 --correlation-length 2 --mean-strength 0.5 --seed 42
powershell -ExecutionPolicy Bypass -File scripts/run_bose_square.ps1 -Canonical -PotentialFile results/speckle.dat -ParameterFile parameter_files/BoseHubbard_speckle_canonical.ini -Lx 32 -Ly 32 -Sweeps 10000 -Thermalization 1000 -RuntimeLimit 60 -JobName canonical_speckle
```

The result is `results/canonical_speckle.out.h5`. To use MPI, add `-Processes 2`
(or your desired process count) to the run command and use a new `-JobName`.
Every rank has N=1024; N is not divided among ranks. The runner selects
`qmc_worm_canonical.exe` for one process or `qmc_worm_canonical_mpi.exe` for MPI
from `dist/windows-square-programmable`.

### Run with a spatially uniform potential from a `.dat` file

Create 1024 identical values, then select that file with the same canonical INI:

```powershell
python -c "import numpy as np; from scripts.generate_speckle import write_potential; write_potential('results/uniform.dat', np.full((32,32), 0.5))"
powershell -ExecutionPolicy Bypass -File scripts/run_bose_square.ps1 -Canonical -PotentialFile results/uniform.dat -ParameterFile parameter_files/BoseHubbard_speckle_canonical.ini -Lx 32 -Ly 32 -Sweeps 10000 -Thermalization 1000 -RuntimeLimit 60 -JobName canonical_uniform
```

`-PotentialFile` overrides the INI's map path. Replace `0.5` with any finite
constant; use `0.0` for a zero map. At fixed N, adding a constant `V0` shifts
reported energy by `V0*N` and leaves the equilibrium relative weights unchanged.
Here the shift is `0.5*1024=512`. A nonuniform map can change the density profile
while its sum stays N. Constant offsets and `mu` can still affect worm sampling
efficiency; `mu` does not tune particle number in a canonical run.

### Turn the site potential off without deleting the `.dat` setting

For a fresh run with zero added site potential:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/run_bose_square.ps1 -Canonical -PotentialFile results/speckle.dat -DisablePotential -ParameterFile parameter_files/BoseHubbard_speckle_canonical.ini -Lx 32 -Ly 32 -Sweeps 10000 -Thermalization 1000 -RuntimeLimit 60 -JobName canonical_no_potential
```

Or edit `onsite_potential_enabled = 0` in the INI and omit `-DisablePotential`.
The configured map file is then ignored and need not exist. The added `V_i` is
zero; hopping, interactions, `mu`, and fixed N remain as configured. This does
not turn off the whole Hamiltonian or the measured `Potential_Energy`.
To enable the map again, set the INI flag to `1` and rerun the speckle or uniform
command above with `-PotentialFile` and a new job name, omitting
`-DisablePotential`. Existing INIs without this flag default
to enabled; an empty or omitted `onsite_potential_file` still means a zero map.
Omitting only `-PotentialFile` does not disable a map already named in the INI.

Use `scripts/run_bose_square.ps1` for these commands. The older
`run_rb87_bkt64.ps1` workflow selects the uniform package and has no map/off
switches. For a grand-canonical map or off run, use
[BoseHubbard_speckle.ini](../parameter_files/BoseHubbard_speckle.ini) and omit
`-Canonical`; all potential choices work the same way.

Changing N, the map, or the enabled flag changes the run setup. Use a **fresh
run and a new job name**. Checkpoint restarts retain the saved map and enabled
state; an explicit map/enabled override on restore is rejected. Editing the
`.dat` or INI does not change a checkpoint's Hamiltonian.

## Controlling numerical disorder

The generator preserves the supplied random-phase, Gaussian Fourier-filter,
FFT-intensity construction. The experimental waist, empirical normalization
factor, inclusive `linspace` grid, and hardcoded random seed are removed.
It needs NumPy; the simulation does not need Python once the file is generated.

- `--nx`, `--ny`: lattice dimensions. `--ny` defaults to `--nx`.
- `--dx`, `--dy`: numerical grid spacing. Defaults are `dx=1`, `dy=dx`.
- `--correlation-length ell`: intensity correlation length, in the same length
  units as the spacing. With `dx=dy=1`, `ell=2` means two lattice spacings.
- `--mean-strength V`: mean of the nonnegative potential, in simulation energy
  units. Each finite map is rescaled to this requested mean exactly up to
  floating-point rounding. It is not an exactly prescribed RMS fluctuation.
- `--seed`: reproducible disorder realization. It is separate from the INI's
  Monte Carlo seed. Reusing it with a different strength rescales the same map.
- `--subtract-mean`: subtract the mean after normalization for a zero-mean
  potential. Strength still refers to the positive intensity before subtraction.
- `--metadata-json`: save generation parameters and measured map statistics.

The correlation convention is the connected **intensity** correlation

```text
C(r) / C(0) approximately exp[-r^2 / (2*ell^2)].
Fourier amplitude filter A(f) = exp[-2*pi^2*ell^2*(fx^2 + fy^2)].
```

Thus `sqrt(2)*ell` is the intensity 1/e radius, and the full width at half maximum
along an axis is `2*sqrt(2*ln(2))*ell`. This Gaussian relation is the continuum,
many-mode limit; fixed-amplitude random-phase modes, finite box size, and sampling
alter the measured correlation. A single finite realization does not have an
exact prescribed fitted width or an exactly exponential intensity distribution.
The FFT generates a periodic map, with box lengths `nx*dx`, `ny*dy`. Resolve the
correlation length with several grid points and use a box much larger than it.
The generator reports sampling warnings for poorly resolved or small boxes.

For a disorder-strength scan, keep seed, spacing, size, and correlation length
fixed and vary only strength. For disorder averaging, generate different seeds
and run each as a separate job. MPI ranks within a job sample the same potential;
they are Monte Carlo replicas, not different disorder realizations. Centering
the potential shifts the effective global chemical potential: grand-canonical
results change unless `mu` is adjusted consistently. At fixed N it changes the
energy offset without changing equilibrium relative weights.

## Different disorder configurations on different workers

Use `run_disorder_ensemble.py` to run independent realizations concurrently.
Each worker starts a separate serial simulation with its own map, Monte Carlo
seed, thermalization, measurements, checkpoint, and output directory. It works
with fixed N and grand-canonical sampling. The operating system schedules these
processes across cores; the script limits concurrency rather than pinning CPUs.
No additional C++ build is needed beyond the programmable package above.
The launcher requires NumPy; the summary also requires h5py. Install missing
packages with `python -m pip install numpy h5py`. On another platform or with
another build directory, pass `--executable /path/to/qmc_worm` (or the canonical
serial executable) to the launcher.

For 24 different maps running concurrently at fixed N=1024:

```powershell
python scripts/run_disorder_ensemble.py --parameter-file parameter_files/BoseHubbard_speckle_canonical.ini --canonical --realizations 24 --workers 24 --correlation-length 2 --mean-strength 0.5 --base-disorder-seed 4200 --base-mc-seed 9100 --output-dir results/canonical_disorder_24
python scripts/summarize_disorder_ensemble.py results/canonical_disorder_24/manifest.json
```

For a smaller trial use `--realizations 4 --workers 2` and a new output directory.
For 96 realizations on 24 workers, use `--realizations 96 --workers 24`: each
worker starts another independent run when it finishes the previous one.
A new map is therefore used for every new *simulation*, with fresh equilibration.
Every canonical realization keeps the INI's total N; N is not split among
workers. For grand-canonical calculations, use `BoseHubbard_speckle.ini` and
omit `--canonical`.

The launcher reads `Lx`, `Ly`, `beta`, `canonical`, hopping, interactions,
`nmax`, `thermalization`, `sweeps`, and `runtimelimit` from the supplied INI.
Set those physical and sampling parameters there before launching. It overrides
the map, enables it, assigns distinct seeds, and chooses separate filenames.
An existing map/off setting in that input does not select this ensemble's maps;
use the original single-map runner when you want the potential disabled.
Use `--dx`, `--dy`, and `--subtract-mean` for the same numerical speckle options
as the generator. Seeds are assigned by realization index, so changing worker
count does not change the generated ensemble. `--prepare-only` writes all maps
and input files without starting simulations, for inspection or cluster use.

Choose a new output directory for each ensemble. The launcher refuses to replace
existing results. Its `manifest.json` records seeds, potential hashes, paths,
and success/failure for every requested realization. Subdirectories `r0000`,
`r0001`, ... each contain `parameters.ini`, `potential.dat`, `potential.json`,
`run.log`, `result.out.h5`, and `state.clone.h5` after running. Individual
checkpoints retain their own map. A failed or unfinished realization is not
silently discarded when computing the summary.

### Why not redraw the map after every worm?

In this code, inserting and closing a worm is one Monte Carlo update of the
existing many-particle worldline configuration, not a fresh equilibrium sample.
The configuration equilibrated under map A is generally not equilibrated under
map B. Unconditionally replacing A with B and making one update does not, in
general, sample the equilibrium distribution associated with B. Repeating that
procedure does not generally produce the desired quenched-disorder average:

```text
quenched mean = average over maps V of the equilibrated thermal mean <O>_V
```

Generating another map is valid when the new simulation is equilibrated before
its measurements. Sampling disorder itself as a thermal variable instead needs
joint acceptance rules and usually targets an annealed ensemble, a different
physical problem. The current solver implements equilibrium sampling with
static disorder, not physical-time fluctuating disorder. Independent fixed-map
QMC runs followed by disorder averaging are also used in the original
[Lin-Sorensen-Ceperley disordered Bose-Hubbard study, Section II.1](https://arxiv.org/pdf/1107.3107).

The existing `qmc_worm_mpi.exe` workflow combines Monte Carlo samples from ranks
sharing the same Hamiltonian. It still rejects different maps across those
ranks. Independent disorder workers use separate serial processes and outputs,
so disorder variation is retained instead of being hidden inside a pooled
Monte Carlo error bar. `--workers` here is distinct from the single-map runner's
`-Processes` option.

### Reading a disorder average

`summary.json` gives every realization equal weight, regardless of its number
of measurements or estimated error. For realization means m_r and R maps:

```text
mean = sum_r(m_r) / R
standard_error = sample_standard_deviation(m_r) / sqrt(R)
thermal_standard_error = sqrt(sum_r(thermal_error_r**2)) / R
```

The first error estimate includes thermal noise in the individual means as well
as variation between disorder maps. The thermal estimate is a diagnostic
component: do not add it to the first estimate in quadrature. With only one
map, the disorder-ensemble standard error is undefined and reported as `null`.
Nonfinite or unavailable per-map thermal errors are also represented by `null`.
A successful run does not prove equilibration; check thermalization, independent
starts, and convergence with the number of disorder realizations.

For each observable the JSON retains `realization_means` and measurement
`counts`, so you can inspect the distribution and difficult samples. It also
reports `Number_variance`, computing `<N^2>_r - <N>_r^2` separately within each
realization before averaging. Computing this from disorder-averaged moments
would incorrectly include differences of mean N between maps as thermal number
fluctuations. Its thermal error is not propagated without the covariance of
those moments. The canonical number and this number variance should remain
exact N and zero, respectively.

## Arbitrary maps and notebook use

The script is importable from a notebook launched in the repository root:

```python
import numpy as np
from scripts.generate_speckle import generate_speckle, write_potential

V = generate_speckle(32, 32, dx=1.0, correlation_length=2.0,
                     mean_strength=0.5, seed=42)
y, x = np.indices(V.shape)
V += 0.002 * ((x - 15.5)**2 + (y - 15.5)**2)  # optional numerical trap
write_potential("results/speckle_plus_trap.dat", V)
```

A file contains exactly one finite real value per site, separated by whitespace;
blank lines and `#` comments are accepted. No coordinate or index columns are
present. For a square grid, `V[y,x]` is stored at `site=x+Lx*y`, so use
`V.ravel(order="C")`. For example, a 3 by 2 map is ordered
`(0,0), (1,0), (2,0), (0,1), (1,1), (2,1)`. General lattice indexing is
`basis + N_BASIS*(x + Lx*(y + Ly*z))`. A missing parameter means a zero map;
when the potential is enabled, a specified but missing, malformed, nonfinite,
or incorrectly sized file is an error. A disabled potential does not read the file. Keep the generator dimensions consistent with the INI dimensions;
comments are descriptive metadata, not simulation parameters.

## Output and restarts

Results store the resolved map under `/hamiltonian/onsite_potential`, with model,
lattice, and enabled-state metadata. A disabled potential stores an all-zero map. Checkpoints carry the corresponding Hamiltonian under
`/simulation/realizations/0/clones/0/checkpoint/hamiltonian`. Restarts use this embedded map, so the source file can
be moved or changed without altering the restored Hamiltonian. Start a fresh run
for a new potential; `--reset-statistics` continues the saved Hamiltonian and is
not a way to switch it. Old nonuniform checkpoints without stored Hamiltonian
metadata are rejected because their parameters were previously hardcoded.

`Density_Distribution` reports site occupations in file order, including for
open boundaries. Fully periodic lattices support `Winding_number_squared`.
Translationally binned one-body and density correlations remain disabled in
nonuniform builds. The scalar `Equal_time_susceptibility` is available for fully
periodic single-basis Bose-Hubbard lattices, including static disorder: it sums
the off-site equal-time worm events before binning and adds physical density.
`Winding_number_squared_total` likewise records the component sum before binning.
For the canonical beta=0.35 scan, matching the clean-system comparison notebook,
see the [line-by-line build, execution and plotting guide](disorder_response_beta035.md). The pre-existing on-site `Density_Matrix` limitation of the
uniform solver is not addressed by this feature.

```python
import h5py
with h5py.File("results/speckle.out.h5") as f:
    potential = f["hamiltonian/onsite_potential"][()].reshape(32, 32)
    density = f["simulation/results/Density_Distribution/mean/value"][()].reshape(32, 32)
```

## Verification

Run generator and independent-disorder workflow checks with:

```powershell
python -m unittest discover -s tests -p test_generate_speckle.py
python -m unittest discover -s tests -p test_disorder_ensemble.py
```

The independent solver checks require NumPy and h5py:

```powershell
python scripts/test_programmable_potential.py --executable dist/windows-square-programmable/qmc_worm.exe --canonical-executable dist/windows-square-programmable/qmc_worm_canonical.exe --mpi-executable dist/windows-square-programmable/qmc_worm_mpi.exe --canonical-mpi-executable dist/windows-square-programmable/qmc_worm_canonical_mpi.exe --uniform-executable dist/windows-square/qmc_worm.exe
```

This compares a small nonuniform open lattice against exact diagonalization,
including zero hopping, finite occupancy interactions, and fixed N. It also
checks constant-map/shifted-mu cases, malformed input, rectangular periodic
geometry, density-sum consistency, and restarts without the source map.
Logs, HDF5 outputs, and numerical comparisons remain under
`build/potential-validation/`. These are implementation regressions;
production calculations still need sampling, size, and occupation-cutoff checks.

### Verified build and run on 19 October 2026

Both `dist/windows-square` and `dist/windows-square-programmable` were built,
including ordinary/canonical serial/MPI executables, using the compatible
ALPSCore install from the neighboring original checkout. All 32 hashes in each
package manifest matched. Compiler diagnostics were the existing Boost template
and unused-variable warnings; there were no build errors.

The 14 generator unit tests passed. The complete integration command above
passed all twelve exact-diagonalization comparisons, missing/changed-map restarts,
explicit map-override and checkpoint-metadata rejection, malformed-map rejection,
rectangular periodic mapping/winding, fixed-N restart, and ordinary/canonical
MPI restart and mixed-rank-map rejection. The constant-map and shifted-mu cases
also produced identical energy estimates for the tested seed.

The enabled/disabled tests additionally checked existing and missing map files,
zero-map physics in ordinary/canonical/uniform builds, persisted disabled state,
restore-override rejection, and disabled-map runs on two MPI ranks. A real
checkpoint created before the enabled flag was introduced also restored with
its original map and the default enabled state. Its log is retained locally in
`build/potential-validation/legacy-toggle-ojtr_hhd/`.

The numerical report and per-case logs are retained locally in
`build/potential-validation/check-jmuvahfd/`. This directory is ignored by Git;
rerunning the test script creates a new independent directory. Build logs are
`build-repro/build-default-final.log` and
`build-repro/build-programmable-final.log`.

An example map is available locally as `results/speckle.dat` (32 by 32, dx=1,
ell=2, mean strength=0.5, seed=42), with `results/speckle.json`. A 200-measurement
runner smoke test produced `results/speckle_smoke.out.h5` and passed configuration
checks. The stored map matches the generated text file exactly. This very short
smoke run is not enough to estimate reliable statistical errors or establish
physical equilibration; use the longer example command and assess convergence
for actual calculations.

The canonical guide commands were also run as short 200-measurement checks on
32 by 32 sites with N=1024. Speckle, constant V=0.5, command-line off, and
INI-controlled off each reported exactly N=1024, N squared=1048576, and a site
density sum of 1024. Both off cases used a deliberately missing map path and
stored all-zero potentials with enabled=0. These logs and HDF5 files are under
`results/guide-checks/`. As with the other short smoke runs, these verify command
wiring and invariants rather than statistical convergence.

### Independent disorder workflow verified on 20 October 2026

All 16 ensemble unit tests passed. They cover reproducible map/seed assignment,
separate artifacts, overwrite protection, fixed-N configuration, worker failures,
equal weighting with unequal measurement counts, uncertainty calculations,
within-map number variance, and refusal of incomplete or altered results.

Real grand-canonical and canonical runs each completed four distinct speckle
realizations on two concurrent workers, in batches. Every saved potential matched
its generated file, and every run produced a separate checkpoint. The canonical
runs on 8 by 8 sites retained exactly N=64 and zero number variance. The summary
CLI and the observed two-process concurrency also passed. These short runs check
implementation and invariants, not convergence of the physical predictions.

Inputs, manifests, logs, HDF5 outputs, and summaries are retained locally under
`build/ensemble-validation/check-0hmfkzno/`, in the `grand_canonical` and
`canonical` subdirectories. Like other validation output, this is ignored by Git.
