# Beta=0.60 disorder sweep

This scan runs on `assess/site-programmable-potentials` as a separate dataset from
the completed beta=0.35 measurements. There are 12 points / 96 independent runs:
V/t = 0.6, 1.0, 2.0, 3.0, each at L=32,64,128 with eight realizations per point.

The settings match the earlier sweep except beta=0.60 and new MC seeds:
N=1875,7500,30000 respectively; t=1, U=0.152, V_nn=0, mu=-3.47, nmax=5,
canonical_window=0.5, speckle correlation length=2 sites. The same disorder seeds
are used, allowing comparison between temperatures on the same maps. MC seeds
use a fresh block beyond those allocated to the initial and refinement scans.

Thermalization is 50,000 excursions; Nmeasure=1, Nmeasure2=100. Per-run wall
budgets including thermalization are 1,800 s (L32), 3,600 s (L64), and 10,800 s
(L128). Up to 24 serial processes run simultaneously, grouped into three active
points with eight workers each. The nominal total budget is 144 CPU hours,
approximately six hours at full utilization plus queue imbalance and I/O.
Runtime limits do not establish physical or statistical convergence.

## Reproduce from build to execution

Run from the checkout containing these changes:

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
git switch assess/site-programmable-potentials
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
& $python -c "import numpy, h5py, matplotlib; print('Dependencies OK')"
```

The current programmable package is already built and its binary hash was checked
against the validated beta=0.35 solver. To rebuild on this desktop if needed:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/build_windows.ps1 -ProgrammablePotential -Jobs 8 -AlpsInstall '\\alfs1.physics.ox.ac.uk\al\kuo\projects\worm\build-repro\alps-install'
if ($LASTEXITCODE -ne 0) { throw 'Build failed' }
```

On another machine, follow [WINDOWS_BUILD.md](../WINDOWS_BUILD.md) for dependencies;
omit `-AlpsInstall` and its path to build the pinned dependency rather than reuse it.

Inspect the complete settings, then prepare once into a fresh output directory:

```powershell
Get-Content parameter_files/disorder_response_beta060.json
& $python scripts/disorder_response_scan.py prepare --config parameter_files/disorder_response_beta060.json --output-dir results/disorder_beta060_20261027
if ($LASTEXITCODE -ne 0) { throw 'Preparation failed' }
```

All raw points are prepared directly under `data/`. Existing results are protected
against overwriting. If repeating, select a new output directory and use it in
all commands below. Do not replace the binary between preparation and launch;
the controller verifies the recorded executable hash.

Launch once, as a detached hidden process:

```powershell
$scanDir = Join-Path $repo 'results\disorder_beta060_20261027'
$controller = Join-Path $repo 'scripts\disorder_response_scan.py'
$scanFile = Join-Path $scanDir 'scan.json'
$scanArgs = '"{0}" run --scan "{1}"' -f $controller, $scanFile
$job = Start-Process -FilePath $python -ArgumentList $scanArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $scanDir 'controller.stdout.log') -RedirectStandardError (Join-Path $scanDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $scanDir 'controller.pid')
$job.Id
```

The controller handles the queue and refreshes tables, seven plots and an executed
notebook after each completed point. Plot labels read beta from the dataset
configuration. No interactive supervision is required.

Inspect later, using:

```powershell
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta060_20261027/scan.json
Get-Content results/disorder_beta060_20261027/controller.stdout.log -Tail 20
Get-Content results/disorder_beta060_20261027/controller.stderr.log -Tail 20
```

The automatic notebook is
`results/disorder_beta060_20261027/disorder_response_comparison.ipynb`.
To rebuild its plots from completed raw results without starting a simulation:

```powershell
& $python scripts/disorder_response_scan.py analyze --scan results/disorder_beta060_20261027/scan.json
```

## Dataset files

- `scan.json`, `config.json`: queue state, effective settings and point paths.
- `data/L32_V0p6`, etc.: all 12 point directories together.
- Each `r0000` to `r0007`: potential, effective INI, log, HDF5 result/checkpoint.
- `analysis/summary.csv`, `analysis/realizations.csv`: averages and individual means.
- `analysis/*.png`: susceptibility, scaled susceptibility, both winding components,
  total winding, size dependence and relative errors.
- `provenance/`: source snapshots, reference notebook and parameter records.

Maps have equal weight. Plotted errors are SD of whole-map means divided by sqrt(8),
including MC noise. Failed or incomplete maps are not silently discarded from
an ensemble. The beta=0.35 results remain a separately indexed dataset.

## Progress check and checkpoint continuation - 28 October 2026

The initial collection completed all eight L32/L64 points. L128 at V/t=0.6 and 1
reached the original three-hour limit with some chains still below the unchanged
50,000-excursion thermalization threshold, and several others had very few bins.
The V/t=2 and 3 chains wrote results/checkpoints, but their completion records were
stale; an `Invalid argument` file-write error was recorded. No inferred exit code
was assigned to those chains.

All four L128 ensembles are continued from their saved states, RNGs and statistics,
with a further 7,200 seconds per chain and 16 simultaneous serial processes.
There are 32 continuations in two batches, nominally four hours plus overhead.
The maps, beta, N, Hamiltonian and thermalization threshold are unchanged.
This adds measurements to the same realizations rather than creating extra maps.
The original inputs, results, checkpoints and logs remain at their initial paths.
New work is under each case's `continuation01/`; the point manifest selects its
current output. Original manifests are `manifest.before_continuation.json`, with
`scan.before_continuation.json` at the dataset root.

Preparation copies checkpoints, modifies both the typed ALPS `runtimelimit` and
its serialized INI value, and preserves an immutable `start.clone.h5`. The solver
loads/saves its separate working `state.clone.h5`. Preparation was tested on a
copied L128 chain: the state advanced, its map and thermalization stayed unchanged,
and the source checkpoint hash was preserved. Manifest-write retries were tested
for transient errors and permanent failure; all 16 ensemble unit tests passed.

The preparation command used was:

```powershell
python scripts/resume_disorder_response.py prepare --scan results/disorder_beta060_20261027/scan.json --seconds 7200 --workers 16
```

Do not repeat preparation on these existing continuation directories. To launch
this prepared continuation once, use the same hidden-launch pattern with the
continuation controller:

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
$scanDir = Join-Path $repo 'results\disorder_beta060_20261027'
$controller = Join-Path $repo 'scripts\resume_disorder_response.py'
$scanFile = Join-Path $scanDir 'scan.json'
$scanArgs = '"{0}" run --scan "{1}"' -f $controller, $scanFile
$job = Start-Process -FilePath $python -ArgumentList $scanArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $scanDir 'continuation.stdout.log') -RedirectStandardError (Join-Path $scanDir 'continuation.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $scanDir 'controller.pid')
$job.Id
```

The controller validates real process exits and all eight results before marking
a point complete, then refreshes both temperatures' notebooks automatically.
It does not omit an unfinished map from a disorder average. Manifest writes now
use unique temporary files and bounded retries for transient network-file errors.

## Paired notebook layout

Both datasets' `disorder_response_simple.ipynb` and
`disorder_response_comparison.ipynb` show each beta=0.60 plot directly below its
beta=0.35 counterpart. The same simple code cell loops over `[0.35, 0.60]`.
The beta=0.60 figures are explicitly marked partial until all 12 points finish;
missing L128 points are listed in the notebook output and
`analysis/paired_analysis_status.json`. Raw HDF5 files are read separately for each
beta, with no mixing of measurements or uncertainty estimates between temperatures.

Refresh both datasets' paired notebooks manually with:

```powershell
python scripts/analyze_disorder_response.py
```

This produces nine plot pairs (18 embedded figures), retaining full-range,
refinement-range and uncertainty views. Existing single-temperature PNG paths in
each dataset's `analysis/` still contain that dataset's own beta only.

## Complete results and main notebook entry point

The continuation completed at 12:05:14 GMT on 28 October 2026. All 12 beta=0.60
points and 96 realizations are complete, including all four L128 points. Every
run has production susceptibility: 586-824 bins at L32, 75-185 at L64, and
435-777 at L128. This records collection success, not convergence.

Open `notebooks/disorder_beta035_beta060.ipynb` from the repository root.
It is the clearly named entry point for the complete comparison, with each
beta=0.60 plot beneath beta=0.35 and 18 embedded figures. The older notebook
paths under both datasets are updated too. If an editor still displays earlier
partial results, close that tab and open this file directly.

`python scripts/analyze_disorder_response.py` also regenerates this main entry
point, so subsequent plot refreshes cannot leave it behind.

## Additional V/t=4,6,8 measurements - 28 October 2026

The next collection adds nine points / 72 runs: V/t=4,6,8 at each of L32,L64,L128,
with eight realizations per point. Beta, fixed N, couplings, nmax, measurement
intervals, 50,000 thermalization excursions and correlation length remain the same.
Unit-strength disorder fields are reused and a fresh MC seed block begins after
the initial beta=0.60 strengths. To accommodate the previous L128 warmup time,
the L128 wall budget is 18,000 s (five hours) per chain; L32 stays at 1,800 s and
L64 at 3,600 s. Maximum concurrency is 16 processes, grouped into two ensembles
of eight. This is about 156 CPU hours / 9.75 hours at full utilization, with
queue imbalance and I/O overhead adding elapsed time.

The new collection lives inside the beta=0.60 dataset at
`extensions/V4_V6_V8_20261028/`, with its own `data/` and manifests. It is registered
in the parent `additional_scans.json`; no physical results are moved or replaced.
The paired notebook shows the additional nine points as pending until complete.
After every point finishes, the controller automatically regenerates the paired
notebooks, including `notebooks/disorder_beta035_beta060.ipynb`. Complete maps
are always required for a point's disorder mean.

From the repository root, the exact preparation command is:

```powershell
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
& $python scripts/disorder_response_scan.py prepare --config parameter_files/disorder_response_beta060_V4_V6_V8.json --output-dir results/disorder_beta060_20261027/extensions/V4_V6_V8_20261028
if ($LASTEXITCODE -ne 0) { throw 'Additional scan preparation failed' }
```

The configuration file records all settings, including `refresh_paired_notebook`
and the longer L128 budget. Use a fresh directory for any repeat. The executable
hash matches the previously validated solver; no rebuild was needed.

Launch once with:

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
$scanDir = Join-Path $repo 'results\disorder_beta060_20261027\extensions\V4_V6_V8_20261028'
$controller = Join-Path $repo 'scripts\disorder_response_scan.py'
$scanFile = Join-Path $scanDir 'scan.json'
$scanArgs = '"{0}" run --scan "{1}"' -f $controller, $scanFile
$job = Start-Process -FilePath $python -ArgumentList $scanArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $scanDir 'controller.stdout.log') -RedirectStandardError (Join-Path $scanDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $scanDir 'controller.pid')
$job.Id
```

Inspect the additional queue later with:

```powershell
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta060_20261027/extensions/V4_V6_V8_20261028/scan.json
```

The parent collection's status command reports its initial 12 points; the paired
analysis reads both the parent and registered extension. At completion beta=0.60
will contain 21 points / 168 runs, with strengths 0.6,1,2,3,4,6,8. Regenerate all
comparison plots at any later time using:

```powershell
& $python scripts/analyze_disorder_response.py
```

## Queued beta=0.45 and beta=0.52 sweeps - 28 October 2026

The follow-on queue waits for the beta=0.60 V/t=4,6,8 extension to complete
successfully, then runs beta=0.45 followed by beta=0.52. Each new temperature
has V/t=1,2,3,4 at L=32,64,128: 12 points / 96 runs per beta, or 192 new runs.
The canonical populations are N=1875,7500,30000. All other physical and
measurement settings match the extension above, including eight realizations,
correlation length two sites, and 50,000 thermalization excursions. Disorder
maps are paired with the existing sweeps, and MC seeds use distinct blocks.

Each scan uses at most 16 serial workers (two eight-realization ensembles).
The per-chain wall limits, including thermalization, are 1,800 / 3,600 / 18,000 s
for L32 / L64 / L128. The new sweeps budget about 208 CPU hours per temperature,
416 CPU hours total: at least 26 hours at full utilization, plus queue imbalance,
analysis/I/O time and the remaining beta=0.60 work. The morning check can show
partial progress; this sequence is not expected to finish overnight.

The existing programmable canonical solver was reused; its recorded SHA256 is
`1a3fefe2985d3623f83e8ead92a6112f7343b1205489408eef4facc0dda54ce9`.
Use the build commands near the start of this guide if rebuilding is needed.
The commands below were used after that build, from the repository root.
Preparation creates all inputs and maps without launching QMC. Run preparation
only once; repeats require fresh output directories and matching queue paths.

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
Get-Content parameter_files/disorder_response_beta045.json
Get-Content parameter_files/disorder_response_beta052.json
& $python scripts/disorder_response_scan.py prepare --config parameter_files/disorder_response_beta045.json --output-dir results/disorder_beta045_20261028
if ($LASTEXITCODE -ne 0) { throw 'Beta 0.45 preparation failed' }
& $python scripts/disorder_response_scan.py prepare --config parameter_files/disorder_response_beta052.json --output-dir results/disorder_beta052_20261028
if ($LASTEXITCODE -ne 0) { throw 'Beta 0.52 preparation failed' }
& $python -m unittest discover -s tests -p test_disorder_queue.py
if ($LASTEXITCODE -ne 0) { throw 'Queue checks failed' }
& $python scripts/run_disorder_queue.py prepare --queue results/disorder_queue_beta045_beta052_20261028/queue.json --after results/disorder_beta060_20261027/extensions/V4_V6_V8_20261028/scan.json --scan results/disorder_beta045_20261028/scan.json --scan results/disorder_beta052_20261028/scan.json
if ($LASTEXITCODE -ne 0) { throw 'Queue preparation failed' }
```

The dataset catalog `parameter_files/disorder_comparison_datasets.json` lists
beta=0.35,0.60,0.45,0.52 in that panel order. Both new configurations set
`refresh_paired_notebook=true`. The simple analysis script reads raw results and
refreshes the shared comparison after every completed point. Open
`notebooks/disorder_temperature_comparison.ipynb`; the older
`notebooks/disorder_beta035_beta060.ipynb` remains an updated alias. Beta=0.60 is
still directly below beta=0.35, with beta=0.45 and beta=0.52 following for each
observable. Queued/incomplete points are listed explicitly and never included
in an eight-map average. All four datasets retain their own analysis tables.

Launch the prepared sequence once as a hidden, detached process:

```powershell
$queueDir = Join-Path $repo 'results\disorder_queue_beta045_beta052_20261028'
$queueFile = Join-Path $queueDir 'queue.json'
$controller = Join-Path $repo 'scripts\run_disorder_queue.py'
$queueArgs = '"{0}" run --queue "{1}"' -f $controller, $queueFile
$job = Start-Process -FilePath $python -ArgumentList $queueArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $queueDir 'controller.stdout.log') -RedirectStandardError (Join-Path $queueDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $queueDir 'controller.pid')
$job.Id
```

The queue checks the dependency once per minute, leaving its running processes
and files alone. When it completes, each new scan gets its own controller log
inside its dataset directory. `queue.json` records the waiting/running/completed
state and stage times. An exclusive lock and prepared-state checks prevent
relaunching this queue. If the dependency fails or its controller disappears
for three checks, the queue stops with a recorded error. A normally finished
scan with failed points permits the next independent temperature to proceed;
an interrupted scan controller stops the sequence for inspection.

Check progress later with these read-only commands:

```powershell
& $python scripts/run_disorder_queue.py status --queue results/disorder_queue_beta045_beta052_20261028/queue.json
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta045_20261028/scan.json
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta052_20261028/scan.json
Get-Content results/disorder_queue_beta045_beta052_20261028/controller.stdout.log -Tail 20
Get-Content results/disorder_queue_beta045_beta052_20261028/controller.stderr.log -Tail 20
```

If a refresh is needed after the controllers finish:

```powershell
& $python scripts/analyze_disorder_response.py
```

No further disorder points are queued beyond these requested sweeps. All work
remains on `assess/site-programmable-potentials` for assessment.

Launch record: the detached queue started at 23:16:32 GMT on 28 October 2026,
PID 16148, initially waiting for the beta=0.60 extension (7/9 points complete).
Verification passed six queue tests, all 192 input/map hashes and seed/population
checks, and the notebook export with 36 figures and explicit queued panels.

## Recovery and resumed queue - 29 October 2026

The overnight extension stopped at 07:27 GMT: seven of nine points passed
validation, while L64,V/t=8 and L128,V/t=8 each had one chain with no production
measurements. All processes exited successfully, but those two chains had only
26,916 and 21,054 warmup excursions, below the required 50,000. The L64 chain
also substantially overran its nominal wall limit. Wall limits therefore do
not guarantee a hard process deadline. Seventy of the extension's 72 runs have
production data. The original waiting queue stopped on this failed dependency,
so neither beta=0.45 nor beta=0.52 had started.

On the morning check, only `L64_V8/r0001` and `L128_V8/r0004` were continued.
The selective continuation option retains all successful run entries and files.
The two source checkpoints were copied into their own `continuation01/` folders;
all HDF5 datasets were verified equal except the increased runtime limit.
State, RNG, disorder, physical parameters and any statistics are retained.
The new budget is 43,200 s (12 hours) per chain, with two workers. Existing
results and checkpoints remain in place. `recovery_audit_20261029.json` records
all 144 original artifact hashes, and `scan.before_continuation.json` plus each
`manifest.before_continuation.json` retain the old metadata. These continuations
are the same realizations, not additional independent disorder maps.

These are the exact recovery and execution commands used. Preparation is already
done; do not repeat it in the same directory or relaunch a running controller.

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
& $python -m unittest discover -s tests -p test_disorder_continuation.py
if ($LASTEXITCODE -ne 0) { throw 'Continuation checks failed' }
& $python scripts/resume_disorder_response.py prepare --scan results/disorder_beta060_20261027/extensions/V4_V6_V8_20261028/scan.json --seconds 43200 --workers 2 --only-missing-measurements
if ($LASTEXITCODE -ne 0) { throw 'Continuation preparation failed' }
$scanDir = Join-Path $repo 'results\disorder_beta060_20261027\extensions\V4_V6_V8_20261028'
$controller = Join-Path $repo 'scripts\resume_disorder_response.py'
$scanFile = Join-Path $scanDir 'scan.json'
$runArgs = '"{0}" run --scan "{1}"' -f $controller, $scanFile
$job = Start-Process -FilePath $python -ArgumentList $runArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $scanDir 'recovery_20261029.stdout.log') -RedirectStandardError (Join-Path $scanDir 'recovery_20261029.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $scanDir 'recovery_20261029.pid')
```

Recovery controller PID: 27168. A fresh dependency queue preserves the previous
failed queue as history and runs the unchanged prepared beta=0.45 and beta=0.52
scans once the recovered beta=0.60 scan passes validation:

```powershell
& $python scripts/run_disorder_queue.py prepare --queue results/disorder_queue_beta045_beta052_20261029_resume/queue.json --after results/disorder_beta060_20261027/extensions/V4_V6_V8_20261028/scan.json --scan results/disorder_beta045_20261028/scan.json --scan results/disorder_beta052_20261028/scan.json
if ($LASTEXITCODE -ne 0) { throw 'Recovery queue preparation failed' }
$queueDir = Join-Path $repo 'results\disorder_queue_beta045_beta052_20261029_resume'
$queueFile = Join-Path $queueDir 'queue.json'
$controller = Join-Path $repo 'scripts\run_disorder_queue.py'
$queueArgs = '"{0}" run --queue "{1}"' -f $controller, $queueFile
$job = Start-Process -FilePath $python -ArgumentList $queueArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $queueDir 'controller.stdout.log') -RedirectStandardError (Join-Path $queueDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $queueDir 'controller.pid')
```

The active queue is now `disorder_queue_beta045_beta052_20261029_resume`,
controller PID 29640, launched at 08:06:56 GMT. Use its path for future checks;
the queue dated 20261028 is historical and remains marked failed.

```powershell
& $python scripts/run_disorder_queue.py status --queue results/disorder_queue_beta045_beta052_20261029_resume/queue.json
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta060_20261027/extensions/V4_V6_V8_20261028/scan.json
Get-Content results/disorder_beta060_20261027/extensions/V4_V6_V8_20261028/recovery_20261029.stdout.log -Tail 20
Get-Content results/disorder_beta060_20261027/extensions/V4_V6_V8_20261028/recovery_20261029.stderr.log -Tail 20
```

The recovery controller refreshes the shared notebooks when each point finishes,
including `notebooks/disorder_temperature_comparison.ipynb`. Three continuation
regression tests passed, including preservation of completed entries, checkpoint
state/RNG, and refusal to alter an active scan or repeat preparation.


## Live progress clarification - 29 October 2026

Read-only debugger snapshots at 2026-10-29T14:16:14.591917+00:00 confirmed that both
recovered chains have passed the 50,000-excursion warmup and are advancing in
production. L64_V8_r0001: 972,165 -> 974,624 excursions. L128_V8_r0004: 310,408 -> 310,770 excursions. 
They remain active because their continuation budget is 12 hours from 08:06 GMT
(approximately 20:06, subject to a long worm delaying the stop callback), rather
than stopping immediately after warmup. The requested maximum is one billion
sweeps. This is now productive running, not evidence of stalled warmup.

Disk progress alone was misleading: Nsave=100,000,000 is counted after warmup,
so no periodic checkpoint is expected at these counters. Final results are saved
when the solver exits. Debugger snapshots and counter derivation are recorded in
`diagnostic_20261029_*.gdb.log` and `diagnostic_20261029_progress.json` in the
extension directory. Processes were immediately detached and left running.

## Queued beta=0.60 refinement V/t=6.5,7,7.5 - 29 October 2026

The additional beta=0.60 collection has nine points / 72 runs: V/t=6.5,7,7.5
at L=32,64,128, with eight realizations each. It runs after the existing queue
finishes beta=0.45 and beta=0.52. The new detached queue waits on that entire
queue, so failures are detected rather than leaving a prepared scan waiting
indefinitely. No running queue or simulation was restarted.

The new inputs retain canonical N=1875,7500,30000, the same physical parameters,
correlation length two sites and paired disorder seeds. A fresh MC seed block
starts at 102610644 before the lattice/realization offsets. Wall budgets are
1,800 / 3,600 / 18,000 seconds for L32 / L64 / L128, including the existing
50,000-excursion warmup; maximum concurrency is 16. As with other wall-limited
collections, every chain must actually reach production before a point is used.
The new nominal budget is 156 CPU hours, at least 9.75 hours at full utilization,
plus queue imbalance and overhead, after the preceding temperatures finish.

The collection is stored under the beta=0.60 dataset at
`extensions/V6p5_V7_V7p5_20261029/`. It is registered in `additional_scans.json`,
so the same analysis script and notebooks include its points when complete.
Beta=0.60 now requests 30 points / 240 runs, with strengths
0.6,1,2,3,4,6,6.5,7,7.5,8. Pending points remain explicit until all eight maps
have valid results. Each completed point automatically refreshes the comparison.

These preparation commands were used with the existing validated executable.
As before, run each preparation only once; repeats require fresh directories.

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
Get-Content parameter_files/disorder_response_beta060_V65_V7_V75.json
& $python scripts/disorder_response_scan.py prepare --config parameter_files/disorder_response_beta060_V65_V7_V75.json --output-dir results/disorder_beta060_20261027/extensions/V6p5_V7_V7p5_20261029
if ($LASTEXITCODE -ne 0) { throw 'Refinement preparation failed' }
& $python -m unittest discover -s tests -p test_disorder_queue.py
if ($LASTEXITCODE -ne 0) { throw 'Queue checks failed' }
& $python scripts/run_disorder_queue.py prepare --queue results/disorder_queue_beta060_V65_V7_V75_20261029/queue.json --after results/disorder_queue_beta045_beta052_20261029_resume/queue.json --scan results/disorder_beta060_20261027/extensions/V6p5_V7_V7p5_20261029/scan.json
if ($LASTEXITCODE -ne 0) { throw 'Refinement queue preparation failed' }
```

The extension registry entry records scan
`extensions/V6p5_V7_V7p5_20261029/scan.json`, strengths `[6.5,7.0,7.5]`,
sizes `[32,64,128]`, and `realizations_per_point=8`.
For a fresh reproduction, register the prepared extension before refreshing
plots with these commands (the current dataset is already registered):

```powershell
@'
from pathlib import Path
import json,sys
sys.path.insert(0,'scripts')
from run_disorder_ensemble import _write_json
p=Path('results/disorder_beta060_20261027/additional_scans.json')
r=json.loads(p.read_text())
entry=dict(scan='extensions/V6p5_V7_V7p5_20261029/scan.json',
           strengths=[6.5,7.0,7.5],sizes=[32,64,128],realizations_per_point=8)
if not any(e['scan']==entry['scan'] for e in r['scans']):
    r['scans'].append(entry)
    _write_json(p,r)
'@ | & $python
```
 It was added only after
all nine point manifests were fully prepared and validated. All 72 input and
potential hashes, particle numbers, runtimes, paired disorder seeds and distinct
MC seeds passed checks; the record is `prelaunch_validation.json` in the extension.

Launch the prepared queue once:

```powershell
$queueDir = Join-Path $repo 'results\disorder_queue_beta060_V65_V7_V75_20261029'
$controller = Join-Path $repo 'scripts\run_disorder_queue.py'
$queueFile = Join-Path $queueDir 'queue.json'
$queueArgs = '"{0}" run --queue "{1}"' -f $controller,$queueFile
$job = Start-Process -FilePath $python -ArgumentList $queueArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $queueDir 'controller.stdout.log') -RedirectStandardError (Join-Path $queueDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $queueDir 'controller.pid')
```

The successor controller started at 20:53:47 GMT, PID 60896, and is waiting on
the active beta=0.45/beta=0.52 queue. Nine queue regression tests passed, including
waiting for every stage of a parent queue and refusing a failed dependency.

```powershell
& $python scripts/run_disorder_queue.py status --queue results/disorder_queue_beta060_V65_V7_V75_20261029/queue.json
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta060_20261027/extensions/V6p5_V7_V7p5_20261029/scan.json
& $python scripts/analyze_disorder_response.py
```

Open `notebooks/disorder_temperature_comparison.ipynb` for the plots. The older
`notebooks/disorder_beta035_beta060.ipynb` alias is also regenerated.


## Queued beta=0.45 extra strengths - 30 October 2026

The additional beta=0.45 collection adds V/t=2.5,3.5,4.5,5,6 at L=32,64,128,
with eight realizations per point: 15 points / 120 runs. It follows the existing
beta=0.52 sweep and the queued beta=0.60 V/t=6.5,7,7.5 refinement. The original
beta=0.45 sweep finished successfully on 30 October, with all 12 points / 96 runs.

The new extension retains canonical N=1875,7500,30000, U/t=0.152, nmax=5,
correlation length two sites, 50,000 warmup excursions and the existing
measurement intervals. It reuses the paired disorder seeds and assigns a fresh
MC block starting at 102610708 before lattice/strength/realization offsets.
Per-chain budgets remain 1,800 / 3,600 / 18,000 seconds for L32 / L64 / L128,
including warmup. Maximum concurrency is 16 processes. The additional budget is
260 CPU hours: at least 16.25 hours at full utilization plus queue imbalance and
I/O, starting only after the preceding queues finish.

Inputs and eventual results are under the beta=0.45 dataset at
`extensions/V2p5_V3p5_V4p5_V5_V6_20261030/`. The registered combined beta=0.45
collection now requests 27 points / 216 runs, at V/t=1,2,2.5,3,3.5,4,4.5,5,6.
The analysis script reads this extension through the parent `additional_scans.json`.
All eight realizations must have production data before a point is plotted.

The existing validated executable is reused. These commands reproduce preparation;
the current preparation is already complete, and repeats require fresh directories.

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
Get-Content parameter_files/disorder_response_beta045_V25_V35_V45_V5_V6.json
& $python scripts/disorder_response_scan.py prepare --config parameter_files/disorder_response_beta045_V25_V35_V45_V5_V6.json --output-dir results/disorder_beta045_20261028/extensions/V2p5_V3p5_V4p5_V5_V6_20261030
if ($LASTEXITCODE -ne 0) { throw 'Beta 0.45 extension preparation failed' }
```

For a fresh reproduction, register the fully prepared extension before refreshing
plots. The following registration is already present in the current dataset:

```powershell
@'
from pathlib import Path
import json,sys
sys.path.insert(0,'scripts')
from run_disorder_ensemble import _write_json
p=Path('results/disorder_beta045_20261028/additional_scans.json')
r=json.loads(p.read_text()) if p.exists() else dict(schema_version=1,parent_scan='scan.json',scans=[])
entry=dict(scan='extensions/V2p5_V3p5_V4p5_V5_V6_20261030/scan.json',
           strengths=[2.5,3.5,4.5,5.0,6.0],sizes=[32,64,128],realizations_per_point=8)
if not any(e['scan']==entry['scan'] for e in r['scans']):
    r['scans'].append(entry)
    _write_json(p,r)
'@ | & $python
```

Prepare the successor queue once, then launch it detached:

```powershell
& $python scripts/run_disorder_queue.py prepare --queue results/disorder_queue_beta045_extra_20261030/queue.json --after results/disorder_queue_beta060_V65_V7_V75_20261029/queue.json --scan results/disorder_beta045_20261028/extensions/V2p5_V3p5_V4p5_V5_V6_20261030/scan.json
if ($LASTEXITCODE -ne 0) { throw 'Beta 0.45 successor queue preparation failed' }
$queueDir = Join-Path $repo 'results\disorder_queue_beta045_extra_20261030'
$controller = Join-Path $repo 'scripts\run_disorder_queue.py'
$queueFile = Join-Path $queueDir 'queue.json'
$queueArgs = '"{0}" run --queue "{1}"' -f $controller,$queueFile
$job = Start-Process -FilePath $python -ArgumentList $queueArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $queueDir 'controller.stdout.log') -RedirectStandardError (Join-Path $queueDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $queueDir 'controller.pid')
```

This queue waits on the entire preceding beta=0.60 refinement queue. A failed
parent queue is reported as an error. Existing controllers and simulations
continue independently; no active queue file is rewritten to add this extension.

Check progress or regenerate the shared notebook with:

```powershell
& $python scripts/run_disorder_queue.py status --queue results/disorder_queue_beta045_extra_20261030/queue.json
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta045_20261028/extensions/V2p5_V3p5_V4p5_V5_V6_20261030/scan.json
& $python scripts/analyze_disorder_response.py
```

Open `notebooks/disorder_temperature_comparison.ipynb` for the updated plots.
The controller automatically refreshes it after every completed point. The new
strengths stay explicitly pending until valid eight-realization results exist.

Launch record: successor queue PID 72396, started 2026-10-30T11:18:51.914008+00:00,
waiting behind the beta=0.60 refinement queue. Validation checked all 120 input
and potential hashes, particle numbers, runtimes and paired disorder seeds, and
confirmed that the new MC seeds are distinct from every registered earlier run.
The validation record is `prelaunch_validation.json` in the extension directory.
