# Disorder scan at beta = 0.35: build, launch, and reproduce

The **current combined dataset** contains 36 points and 288 runs, including
V/t=0.6, 0.7, 0.8 and 0.9 at all three sizes. Open
`results/disorder_beta035_20261026/disorder_response_simple.ipynb`.
Both top-level notebooks now include every collected strength. The dataset's
`README.md` describes the consolidated layout; the earlier commands below record
how the two collections were originally built and launched.

This is a test on `assess/site-programmable-potentials`. Nothing is merged into
`main`. The reference is `origin/main` commit `1089533`,
`results/cleanRbRDM/density_matrix_comparison.ipynb`. The scan stores the full
reference commit and a copy of that notebook in its own provenance directory.

## What this campaign measures

The canonical square lattice has periodic boundaries, beta=0.35, hopping t=1,
U=0.152, nearest-neighbor interaction zero, mu=-3.47, and nmax=5. The populations
match the clean notebook: N=1875 for L=32, N=7500 for L=64, N=30000 for L=128.
All three have density 1.8310546875. In canonical sampling mu affects sampling
and the energy convention, while the total measured particle number stays fixed.

The mean nonnegative speckle strengths are V/t = 0, 0.5, 1, 2, 3, 4, 6, 8.
The numerical intensity correlation length is two lattice spacings, dx=dy=1.
Eight independently seeded realizations are measured at each of the 24 points:
192 simulations in total. Within each size, the same unit-strength maps are
rescaled at different strengths; Monte Carlo seeds differ across strengths.
Thus different strengths share disorder samples, with correlated errors across
strengths. The zero-strength point uses eight independent MC chains of the clean
system. Nonnegative maps are used directly; adding a constant to a potential has
no effect on canonical equilibrium relative weights.

The new `Equal_time_susceptibility` result is exactly the notebook observable:

```text
chi_s = (1/L^2) sum_ij <b_i^dagger b_j>
      = sum_displacements G - G(0,0) + N/L^2.
```

There is no extra beta factor or site-count division. The off-site worm events
are summed before binning, and physical density supplies the on-site term.
This works for a fixed nonuniform map without assuming translation invariance
of that map. It records a scalar with covariance across all off-site events;
adding the marginal errors of a matrix would discard that covariance. The
programmable build still does not record the full displacement matrix.
`Winding_number_squared` records both components as before;
`Winding_number_squared_total` records their sum before estimating its error.
The new partial susceptibility histogram is saved/restored with checkpoints and
cleared by thermalization completion and `--reset-statistics`.

Each realization has 50,000 thermalization worm excursions followed by production.
`Nmeasure=1`; `Nmeasure2=100` bins 100 excursions into each scalar susceptibility
sample. The sweep ceiling is 1,000,000,000, so wall time normally limits the run.
The wall budgets, including thermalization, are 1,800 s for L32, 3,600 s for L64,
and 10,800 s for L128. Initialization, configuration tests, and HDF5 I/O add
elapsed time. Up to three points, each with eight serial processes, run together:
24 processes total. The configured total budget is about 288 CPU hours, or roughly
12 hours at full 24-process utilization, with queue imbalance and overhead adding
time. It is a bounded exploratory campaign, not a guarantee of convergence.

A point must have production susceptibility, winding and density measurements
in all eight successful outputs before it appears in the plots. If a wall limit
is reached during thermalization, that point is marked failed for analysis;
its logs and checkpoints remain, and the queue continues to other points.
The controller does not silently discard a failed disorder sample.

This retains the clean notebook's nmax=5 for direct comparison. Strong disorder
can increase local occupations; occupation-cutoff, thermalization, production
length and disorder-sample convergence must be checked before physical claims.

## Commands from this desktop, line by line

Open PowerShell. Run each numbered block in order. These are the same build and
execution commands used for this test; use a new output-directory name for a new
campaign, because existing scans and results are protected against overwriting.

1. Select this checkout and confirm the branch.

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
git switch assess/site-programmable-potentials
git status --short --branch
```

2. Select Python and check the analysis dependencies. If imports succeed, no
installation command is needed.

```powershell
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
& $python -c "import numpy, h5py, matplotlib; print('Python analysis dependencies OK')"
```

If a package is missing:

```powershell
& $python -m pip install numpy h5py matplotlib
```

3. Build the programmable ordinary/canonical serial/MPI package, reusing the
compatible ALPSCore installation in the original neighboring checkout.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/build_windows.ps1 -ProgrammablePotential -Jobs 8 -AlpsInstall '\\alfs1.physics.ox.ac.uk\al\kuo\projects\worm\build-repro\alps-install' *> build-repro/build-disorder-response-20261026.log
if ($LASTEXITCODE -ne 0) { throw 'Programmable build failed; inspect the build log' }
```

On a machine without that prebuilt ALPS install, omit `-AlpsInstall` and its path:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/build_windows.ps1 -ProgrammablePotential -Jobs 8
```

The script checks its required compiler/CMake inputs. See
[WINDOWS_BUILD.md](../WINDOWS_BUILD.md) for installation prerequisites.

4. For the scalar estimator regression, also build the uniform reference package.
This is needed for the density-matrix identity check, not for production scans.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/build_windows.ps1 -Jobs 8 -AlpsInstall '\\alfs1.physics.ox.ac.uk\al\kuo\projects\worm\build-repro\alps-install' *> build-repro/build-uniform-response-20261026.log
if ($LASTEXITCODE -ne 0) { throw 'Uniform reference build failed; inspect the build log' }
& $python scripts/test_disorder_response.py
if ($LASTEXITCODE -ne 0) { throw 'Scalar response regression failed' }
```

The test compares susceptibility on a periodic 3x3 canonical hard-core lattice
against exact diagonalization, including nonuniform and zero-hopping cases. It
also checks the clean scalar against the complete density-matrix sum, winding
component/total consistency, fixed N, and checkpoint reset. Evidence is written
to a new `build/response-validation/check-*` directory.

5. Read or edit the complete campaign configuration. All changes belong in this
file before preparing a new output directory.

```powershell
Get-Content parameter_files/disorder_response_beta035.json
```

`strengths` changes the x-axis points, `realizations` changes the number of maps,
`workers` bounds concurrency, and `runtime_seconds` controls the wall budget for
each size. `workers` must be a multiple of `realizations`; e.g. 24 and 8 gives
three concurrent ensembles. Change `correlation_length` to change numerical
speckle length. `parameters` contains the common Hamiltonian and measurement
settings; `particles` sets the exact total N independently for each L.

6. Prepare every map and input file without starting a solver. The default
configuration is the file displayed above. This command records the binary hash,
reference notebook, effective inputs, seeds, map hashes, source snapshots and diff.

```powershell
& $python scripts/disorder_response_scan.py prepare --output-dir results/disorder_beta035_20261026
if ($LASTEXITCODE -ne 0) { throw 'Scan preparation failed' }
```

To use an edited configuration saved elsewhere, add `--config path/to/config.json`.
To select another build, add `--executable path/to/qmc_worm_canonical.exe`.
Do not rebuild or replace the executable between preparing and launching: the
controller verifies its hash before starting.

7. Launch once as a detached, hidden process. It keeps running when this
PowerShell command finishes. Keep the computer awake for the campaign.

```powershell
$scanDir = Join-Path $repo 'results\disorder_beta035_20261026'
$controller = Join-Path $repo 'scripts\disorder_response_scan.py'
$scanFile = Join-Path $scanDir 'scan.json'
$scanArgs = '"{0}" run --scan "{1}"' -f $controller, $scanFile
$job = Start-Process -FilePath $python -ArgumentList $scanArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $scanDir 'controller.stdout.log') -RedirectStandardError (Join-Path $scanDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $scanDir 'controller.pid')
$job.Id
```

Launch exactly once. The controller refuses a scan already marked running or
completed. It automatically launches queued points and regenerates the tables,
plots and executed notebook after each completed point. No interactive supervision
or repeated assistant checks are required. A controller crash requires inspecting
its manifests/checkpoints before planning a continuation; automatic crash recovery
is deliberately not assumed.

For a foreground run instead, replace the detached launch with:

```powershell
& $python scripts/disorder_response_scan.py run --scan results/disorder_beta035_20261026/scan.json
```

8. Return later and inspect status. This command only reads the saved manifests.

```powershell
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta035_20261026/scan.json
Get-Content results/disorder_beta035_20261026/controller.stdout.log -Tail 20
Get-Content results/disorder_beta035_20261026/controller.stderr.log -Tail 20
```

`running` is the controller's persisted state. To check the process exists too:

```powershell
$controllerId = [int](Get-Content results/disorder_beta035_20261026/controller.pid)
Get-Process -Id $controllerId -ErrorAction SilentlyContinue
```

9. Plots are generated automatically. To regenerate them later from completed
outputs, without running another simulation:

```powershell
& $python scripts/disorder_response_scan.py analyze --scan results/disorder_beta035_20261026/scan.json
```

Open `results/disorder_beta035_20261026/disorder_response_comparison.ipynb` in
your notebook editor. It contains ordinary NumPy/Matplotlib plotting cells and
embedded figures. Its default root discovery works when opened in the scan
folder; otherwise set `root` in the first code cell to that absolute scan path.
The notebook reads saved CSV tables; `analyze` regenerates those tables from raw
HDF5 results. Do not run two analyzers concurrently in the same output folder.

## Files and interpretation

- `scan.json`: complete configuration, queue state, controller PID, binary hash.
- `inputs/*.ini`: effective source INIs for each size/strength point.
- `L32_V0`, `L64_V0p5`, etc.: one ensemble directory per point.
- Each point's `manifest.json`: per-realization seeds, hashes, paths and outcome.
- Each `r0000`, `r0001`, etc.: static `potential.dat`, `potential.json`, effective
  `parameters.ini`, `run.log`, `result.out.h5`, `state.clone.h5`, and private temp.
- `analysis/summary.csv`: one equal-disorder-weight mean and one SE for each
  complete point; susceptibility, scaled susceptibility, both winding components
  and total winding.
- `analysis/realizations.csv`: every individual whole-map mean, its stored MC
  error, measurement counts, seeds and raw output path.
- `analysis/analysis_status.json`: completed/requested point counts and all omitted
  incomplete points. A partial figure is explicitly a partial scan.
- `analysis/*.png`: susceptibility, scaled susceptibility, Wx squared, Wy squared,
  total winding versus disorder; susceptibility versus size and relative errors.
- `provenance/`: source snapshots, working diff and original reference notebook.

Dots show each realization. The mean weights maps equally, regardless of histogram
count. The plotted one-SE error is `SD(realization means, ddof=1)/sqrt(R)` and
already includes finite-MC noise. Stored per-map MC errors are diagnostics and
must not be added again. Scalar whole-map means preserve spatial and winding-
component covariance. The new error bars measure disorder-ensemble uncertainty;
the reference clean notebook used a rank jackknife for one clean potential.

Continue testing on this branch. Merging into `main` is a separate later step
once the measurements and numerical convergence have been assessed.

## Verification before the 26 October 2026 launch

Both packages rebuilt successfully. `scripts/test_disorder_response.py` passed:
zero-hopping chi=0.222222222222 (exact 2/9), nonuniform chi=0.594763 +/- 0.007913
(exact 0.591205), and clean chi=0.592150 +/- 0.007369 (exact 0.598232).
The clean scalar matched the corrected full density-matrix sum to roundoff;
component/total winding, fixed N, and checkpoint reset checks passed. Full evidence:
`build/response-validation/check-mszetzn6/`.

The complete unattended controller was exercised on two strengths with two
realizations each. Preparation, queueing, HDF5 validation, all seven PNG figures,
CSV export and executed notebook generation passed. Evidence:
`build/response-validation/scan-smoke-dena7v5i/`.
These are implementation checks; they do not establish production convergence.

The existing 16 independent-disorder launcher and summary unit tests also passed
before launch. The production input set is `results/disorder_beta035_20261026/`.

## Completed results and simple analysis - 27 October 2026

The campaign completed at 01:25:38 GMT on 27 October: all 24 points and 192
realizations succeeded, with no failed or skipped runs. Every run collected
production susceptibility. Per-run susceptibility-bin ranges were 2,078-3,696
for L32, 1,062-4,885 for L64, and 881-17,801 for L128; all stored susceptibility
MC errors were finite. This verifies measurement collection, not convergence.

Open the executed notebook:
`results/disorder_beta035_20261026/disorder_response_simple.ipynb`.
It now contains nine figures, one per cell, with raw-HDF5 loading and averaging
visible in ordinary Python, including all refinement points. No cached summary
input or custom analysis module is needed. Earlier collection-specific notebooks
are retained under `campaigns/`; both top-level notebooks show the combined data.

From this checkout's root, reproduce the same calculations and PNG files with:

```powershell
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
& $python scripts/analyze_disorder_response.py
```

The standalone script saves files without opening interactive windows. Its
`# %%` sections match the simple notebook's cells. Outputs are
`results/disorder_beta035_20261026/analysis/summary.csv` and nine PNGs
in that same directory, together with the 288-row `realizations.csv`. Each raw map hash, beta, shape, fixed N, positive
measurement count and winding sum is checked. The simple analysis was also
compared against both archived collections for all 36 means and standard errors.

The scaled susceptibility reverses size ordering between V/t=0.5 and 1, while
winding decreases strongly across the scan. This is an apparent finite-size
crossing in the measured curves; more disorder samples, a finer strength grid,
and convergence checks are needed before estimating a transition precisely.

## Additional strengths 0.6, 0.7, 0.8, 0.9 - 27 October 2026

The refinement adds all four strengths at L32, L64 and L128, with eight
realizations per point: 12 additional points and 96 runs. The physical parameters,
fixed particle numbers, correlation length, thermalization, size-specific wall
budgets and maximum 24-process concurrency match the original campaign.
Unit-strength disorder seeds are reused, so the new points rescale the same
maps. Fresh MC seed ranges start after all streams allocated to the original
strengths. The original results remain available at the dataset root; the added
scan is stored under `refinement_V06_V09_20261027/` within that dataset.

The same solver binary hash was checked against the completed original scan.
No rebuild is required. From this checkout's root, prepare with:

```powershell
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
& $python scripts/disorder_response_scan.py prepare --config parameter_files/disorder_response_beta035_refinement.json --output-dir results/disorder_beta035_20261026/refinement_V06_V09_20261027
if ($LASTEXITCODE -ne 0) { throw 'Refinement preparation failed' }
```

Launch once, using the same detached procedure:

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
$scanDir = Join-Path $repo 'results\disorder_beta035_20261026\refinement_V06_V09_20261027'
$controller = Join-Path $repo 'scripts\disorder_response_scan.py'
$scanFile = Join-Path $scanDir 'scan.json'
$scanArgs = '"{0}" run --scan "{1}"' -f $controller, $scanFile
$job = Start-Process -FilePath $python -ArgumentList $scanArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $scanDir 'controller.stdout.log') -RedirectStandardError (Join-Path $scanDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $scanDir 'controller.pid')
$job.Id
```

The controller manages all new points unattended and produces their own figures
and executed notebook. To inspect this refinement later:

```powershell
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta035_20261026/refinement_V06_V09_20261027/scan.json
```

The parent dataset's `additional_scans.json` records the extension location.
The additional collection is complete and is now incorporated in the combined
notebooks and tables. Its collection records are archived under
`campaigns/refinement_V06_V09_20261027/`, and its raw point directories are beside
the original points under `data/`.

## Consolidated layout and current reproduction command

All 36 size/strength points and 288 runs are indexed by the dataset root's
`scan.json`, with strengths sorted as
`0, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 6, 8`.
All raw point directories are now together under `data/`. The initial and
refinement scan records, effective source inputs, logs, provenance and previous
plots are retained under `campaigns/original/` and
`campaigns/refinement_V06_V09_20261027/`. Collection-time scan records are saved
as `scan.at_collection.json`; the archived `scan.json` files have updated paths.
The dataset `relocation.json` records the moves and verifies every ensemble
manifest plus 1,152 required artifact sizes. Raw inputs and HDF5 contents were
preserved.

Use the same single command to refresh **both** top-level notebooks, all tables,
and all nine figures from the consolidated raw results:

```powershell
python scripts/analyze_disorder_response.py
```

Current figures and tables are directly in
`results/disorder_beta035_20261026/analysis/`. The nine figures include the prior
seven comparisons and two detailed views over V/t=0.5 to 1. The uncertainty-only
panel spaces sampled strengths evenly for readability; all physical-observable
plots use the actual numerical disorder strength.

To check the complete combined dataset using the existing status tool:

```powershell
python scripts/disorder_response_scan.py status --scan results/disorder_beta035_20261026/scan.json
```

Analysis discovers any future completed collections registered in the root's
`additional_scans.json`. Duplicate size/strength points and incompatible physical
parameters are rejected, preventing additions from being silently omitted or
counted twice. There is no need to edit a strength list inside the notebook.

Validation after consolidation passed: all 36 means and SEs agree with their
original collection tables to rounding tolerance; all 288 distinct raw-output
paths resolve; all 288 MC seeds are distinct; both notebooks have nine embedded
figures and valid code cells. No simulation results were recalculated or changed
by the filesystem moves.

For the separate beta=0.60 sweep at V/t=0.6, 1, 2, 3, see
[the beta=0.60 build and execution commands](disorder_response_beta060.md).

## Paired beta=0.35 / beta=0.60 notebook

The main notebooks now place beta=0.60 immediately below beta=0.35 for every
observable. The beta=0.35 dataset remains complete at 36 points and 288 runs.
Partial beta=0.60 points are explicitly marked and excluded as whole ensembles
until every realization completes. See the beta=0.60 guide for its L128 checkpoint
continuation. The same `python scripts/analyze_disorder_response.py` command
regenerates both datasets' paired notebooks and per-temperature figures.
