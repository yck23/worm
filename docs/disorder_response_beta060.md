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


## Queued beta=0.60 at V/t=9 - 31 October 2026

Disorder 9 adds three points / 24 runs at beta=0.60: L=32,64,128, with eight
realizations each. It follows the queued beta=0.45 extra strengths. The preceding
beta=0.60 refinement and beta=0.45 extras keep their order and running controllers.

All physical settings and wall budgets match the beta=0.60 refinement: canonical
N=1875,7500,30000, U/t=0.152, nmax=5, correlation length two sites, 50,000 warmup
excursions, and 1,800 / 3,600 / 18,000 seconds per chain for L32 / L64 / L128.
Maximum concurrency is 16. Paired disorder seeds are retained, with a new MC
seed block starting at 102610772 before lattice/realization offsets. The nominal
budget is 52 CPU hours; the two-ensemble schedule takes about 5.5 hours plus
overhead after earlier work, with long worm excursions able to delay time limits.

The extension is `results/disorder_beta060_20261027/extensions/V9_20261031/`,
registered in its parent's `additional_scans.json`. Beta=0.60 now requests 33
points / 264 runs. New points appear in the shared notebook only after all eight
realizations have production data. Preparation checked all 24 input and potential
hashes, fixed populations, wall budgets, paired disorder seeds, distinct MC seeds,
and the executable hash; see the extension's `prelaunch_validation.json`.

The following commands reproduce preparation using the existing validated binary.
The current preparation and launch are already done; repeats need fresh directories.

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
Get-Content parameter_files/disorder_response_beta060_V9.json
& $python scripts/disorder_response_scan.py prepare --config parameter_files/disorder_response_beta060_V9.json --output-dir results/disorder_beta060_20261027/extensions/V9_20261031
if ($LASTEXITCODE -ne 0) { throw 'Disorder 9 preparation failed' }
```

Register the fully prepared extension before plotting. This entry is already
present in the current dataset:

```powershell
@'
from pathlib import Path
import json,sys
sys.path.insert(0,'scripts')
from run_disorder_ensemble import _write_json
p=Path('results/disorder_beta060_20261027/additional_scans.json')
r=json.loads(p.read_text())
entry=dict(scan='extensions/V9_20261031/scan.json',strengths=[9.0],
           sizes=[32,64,128],realizations_per_point=8)
if not any(e['scan']==entry['scan'] for e in r['scans']):
    r['scans'].append(entry)
    _write_json(p,r)
'@ | & $python
```

Prepare and launch the successor queue once:

```powershell
& $python scripts/run_disorder_queue.py prepare --queue results/disorder_queue_beta060_V9_20261031/queue.json --after results/disorder_queue_beta045_extra_20261030/queue.json --scan results/disorder_beta060_20261027/extensions/V9_20261031/scan.json
if ($LASTEXITCODE -ne 0) { throw 'Disorder 9 queue preparation failed' }
$queueDir = Join-Path $repo 'results\disorder_queue_beta060_V9_20261031'
$controller = Join-Path $repo 'scripts\run_disorder_queue.py'
$queueFile = Join-Path $queueDir 'queue.json'
$queueArgs = '"{0}" run --queue "{1}"' -f $controller,$queueFile
$job = Start-Process -FilePath $python -ArgumentList $queueArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $queueDir 'controller.stdout.log') -RedirectStandardError (Join-Path $queueDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $queueDir 'controller.pid')
```

Inspect or refresh later:

```powershell
& $python scripts/run_disorder_queue.py status --queue results/disorder_queue_beta060_V9_20261031/queue.json
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta060_20261027/extensions/V9_20261031/scan.json
& $python scripts/analyze_disorder_response.py
```

Open `notebooks/disorder_temperature_comparison.ipynb` for the current plots.
The controller refreshes the notebook after every completed point.

Launch record: disorder 9 queue PID 38792, started 2026-10-31T11:11:50.885067+00:00.
Progress at enqueue: beta=0.35 has 36/36 points complete; beta=0.45 has
12/27 complete (15 extra points queued); beta=0.52 has 12/12 complete.
Beta=0.60 has 29/33 complete: L128,V/t=7.5 is running with eight chains,
and the three disorder 9 points are queued. No failed points were recorded.
The timestamped machine-readable snapshot is `progress_at_enqueue.json` in the
disorder 9 queue directory.


## Queued beta=0.52 at V/t=4.5,5,6 - 31 October 2026

The beta=0.52 additions are exactly V/t=4.5,5,6 at L=32,64,128, with eight
realizations per point: nine points / 72 runs. The successor queue waits for
the previously queued beta=0.60 disorder 9 collection. The combined beta=0.52
dataset requests 21 points / 168 runs at strengths 1,2,3,4,4.5,5,6.

Physical settings, fixed populations and run lengths match the original beta=0.52
sweep: N=1875,7500,30000, U/t=0.152, nmax=5, correlation length two sites,
50,000 warmup excursions, and 1,800 / 3,600 / 18,000 seconds per chain for
L32 / L64 / L128. Maximum concurrency is 16 processes. Disorder fields remain
paired with earlier sweeps; MC seeds use a fresh block starting at 102610836
before lattice/strength/realization offsets. The nominal additional budget is
156 CPU hours, at least 9.75 hours at full utilization plus queue imbalance and
overhead, after all previously requested work.

The new collection lives under the beta=0.52 dataset at
`extensions/V4p5_V5_V6_20261031/`, registered in `additional_scans.json`.
Every point needs all eight valid production results before it is plotted.
The shared comparison notebook and dataset tables are refreshed automatically
when each point completes.

Using the existing validated executable, preparation was:

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
Get-Content parameter_files/disorder_response_beta052_V45_V5_V6.json
& $python scripts/disorder_response_scan.py prepare --config parameter_files/disorder_response_beta052_V45_V5_V6.json --output-dir results/disorder_beta052_20261028/extensions/V4p5_V5_V6_20261031
if ($LASTEXITCODE -ne 0) { throw 'Beta 0.52 extension preparation failed' }
```

For a fresh reproduction, register the fully prepared extension with:

```powershell
@'
from pathlib import Path
import json,sys
sys.path.insert(0,'scripts')
from run_disorder_ensemble import _write_json
p=Path('results/disorder_beta052_20261028/additional_scans.json')
r=json.loads(p.read_text()) if p.exists() else dict(schema_version=1,parent_scan='scan.json',scans=[])
entry=dict(scan='extensions/V4p5_V5_V6_20261031/scan.json',
           strengths=[4.5,5.0,6.0],sizes=[32,64,128],realizations_per_point=8)
if not any(e['scan']==entry['scan'] for e in r['scans']):
    r['scans'].append(entry)
    _write_json(p,r)
'@ | & $python
```

These commands prepare and launch the successor queue once. Current inputs and
registration are already prepared; repeats require fresh directories.

```powershell
& $python scripts/run_disorder_queue.py prepare --queue results/disorder_queue_beta052_extra_20261031/queue.json --after results/disorder_queue_beta060_V9_20261031/queue.json --scan results/disorder_beta052_20261028/extensions/V4p5_V5_V6_20261031/scan.json
if ($LASTEXITCODE -ne 0) { throw 'Beta 0.52 queue preparation failed' }
$queueDir = Join-Path $repo 'results\disorder_queue_beta052_extra_20261031'
$controller = Join-Path $repo 'scripts\run_disorder_queue.py'
$queueFile = Join-Path $queueDir 'queue.json'
$queueArgs = '"{0}" run --queue "{1}"' -f $controller,$queueFile
$job = Start-Process -FilePath $python -ArgumentList $queueArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $queueDir 'controller.stdout.log') -RedirectStandardError (Join-Path $queueDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $queueDir 'controller.pid')
```

Inspect or refresh later:

```powershell
& $python scripts/run_disorder_queue.py status --queue results/disorder_queue_beta052_extra_20261031/queue.json
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta052_20261028/extensions/V4p5_V5_V6_20261031/scan.json
& $python scripts/analyze_disorder_response.py
```

Open `notebooks/disorder_temperature_comparison.ipynb` for current plots. The new
points stay explicitly pending until their eight-realization averages are valid.

Launch record: beta=0.52 extra queue PID 63960, started 2026-10-31T11:24:29.477310+00:00,
waiting after the beta=0.60 disorder 9 queue. All 72 prepared runs passed input
and potential hash, canonical population, runtime, paired-disorder-seed and
distinct-MC-seed checks. The exact added strengths are 4.5,5,6.


## Queued beta=0.39 at V/t=0.6,1,2 - 31 October 2026

The beta=0.39 sweep adds nine points / 72 runs: V/t=0.6,1,2 at L=32,64,128,
with eight realizations each. It follows the existing beta=0.52 extra-point
queue, preserving all previously requested work. On the progress check, all nine
beta=0.60 refinement points were complete and the beta=0.45 extras had started
16 serial chains. No failed points or queue errors were recorded.

The new dataset is `results/disorder_beta039_20261031/`, with raw runs under
`data/`. Fixed populations are N=1875,7500,30000. Other settings match the recent
sweeps: U/t=0.152, nmax=5, correlation length two sites, 50,000 warmup excursions,
Nmeasure=1 and Nmeasure2=100. Per-chain wall budgets are 1,800 / 3,600 / 18,000 s
for L32 / L64 / L128, with maximum concurrency 16. Paired disorder fields are
reused and a distinct MC seed block starts at 102610900 before the lattice,
strength and realization offsets. The nominal budget is 156 CPU hours, at least
9.75 hours at full utilization plus queue imbalance and overhead, after the
preceding queues finish.

The existing executable was reused. Preparation commands are below; preparation
and launch have already been performed for the current dataset. Repeats need
fresh output directories and corresponding queue/catalog paths.

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
Get-Content parameter_files/disorder_response_beta039.json
& $python scripts/disorder_response_scan.py prepare --config parameter_files/disorder_response_beta039.json --output-dir results/disorder_beta039_20261031
if ($LASTEXITCODE -ne 0) { throw 'Beta 0.39 preparation failed' }
```

Register the fully prepared dataset in the comparison catalog. This is already
present in the current checkout:

```powershell
@'
from pathlib import Path
import json,sys
sys.path.insert(0,'scripts')
from run_disorder_ensemble import _write_json
p=Path('parameter_files/disorder_comparison_datasets.json')
r=json.loads(p.read_text())
entry=dict(beta=.39,directory='results/disorder_beta039_20261031')
if not any(e['beta']==entry['beta'] for e in r['datasets']):
    r['datasets'].append(entry)
    _write_json(p,r)
'@ | & $python
```

Prepare the successor queue once and launch it detached:

```powershell
& $python scripts/run_disorder_queue.py prepare --queue results/disorder_queue_beta039_20261031/queue.json --after results/disorder_queue_beta052_extra_20261031/queue.json --scan results/disorder_beta039_20261031/scan.json
if ($LASTEXITCODE -ne 0) { throw 'Beta 0.39 queue preparation failed' }
$queueDir = Join-Path $repo 'results\disorder_queue_beta039_20261031'
$controller = Join-Path $repo 'scripts\run_disorder_queue.py'
$queueFile = Join-Path $queueDir 'queue.json'
$queueArgs = '"{0}" run --queue "{1}"' -f $controller,$queueFile
$job = Start-Process -FilePath $python -ArgumentList $queueArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $queueDir 'controller.stdout.log') -RedirectStandardError (Join-Path $queueDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $queueDir 'controller.pid')
```

Check later or regenerate the notebooks with:

```powershell
& $python scripts/run_disorder_queue.py status --queue results/disorder_queue_beta039_20261031/queue.json
& $python scripts/disorder_response_scan.py status --scan results/disorder_beta039_20261031/scan.json
& $python scripts/analyze_disorder_response.py
```

The shared `notebooks/disorder_temperature_comparison.ipynb` now uses five
ordered beta panels per observable: 0.35,0.60,0.45,0.52,0.39 (45 embedded figures).
The beta=0.60 panel remains directly below beta=0.35. A simple analysis script is
also copied into each dataset's `analysis/simple_analysis.py`. The new panels
remain explicitly pending until complete eight-map points have production data.
Each completed point triggers an automatic refresh of all registered datasets.


## Recovery handoff and current active queue - 1 November 2026

All beta=0.45 points are complete, including the extra strengths (27 points /
216 runs). The beta=0.60 refinement at 6.5,7,7.5 is complete. Disorder 9 has
22 of 24 process results, but one L64 map and one L128 map are still executing
very long warmup worms. Read-only live snapshots showed 26,897 and 20,272
completed excursions, below the required 50,000; repeated snapshots were
unchanged while both processes consumed CPU. The solver checks its wall limit
after a worm closes, so these processes can substantially overrun their budgets.
They were left running; their state and the 22 completed results were preserved.

The beta=0.39 inputs were prepared on 31 October, but spend-cap approval failure
prevented that day's registration and launch. On 1 November all 72 inputs/maps,
fixed populations, runtimes, paired disorder seeds and fresh MC seeds were
validated. The dataset is now registered in the five-temperature catalog.

A replacement queue now owns the remaining work. It waits for the active
disorder 9 scan to finish and its controller to exit. If the scan fails validation
because chains have no production data, it copies their checkpoints, retains
state/RNG/maps/statistics, and continues only those missing chains with a
43,200-second budget and two workers. Valid runs remain untouched. It permits
one recovery attempt and records an error if that attempt fails. Once the scan
is complete it runs the beta=0.52 additions, then beta=0.39.

The two former waiting-only Python controllers (63960 and 18272) were stopped
after verifying their waiting states and absence of work children; their console
hosts were the only children. Their queue files were backed up as
`queue.before_recovery_handoff_20261101.json` and marked `superseded`, with a
`replacement_queue` pointer. Their manifests and prepared inputs were preserved.
The active QMC processes 54868 and 67308 were not stopped. Do not relaunch the
superseded queues alongside this replacement.

The recovery feature passed 15 queue regression tests and three checkpoint
continuation tests. The tests cover waiting for the original controller's exit,
continuing only missing chains, preserving completed data, a single recovery
attempt, failure handling, and successors waiting through the recovery state.

These commands prepared the replacement queue (already launched now; use fresh
queue output directories for any repeat):

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
& $python -m unittest discover -s tests -p test_disorder_queue.py
if ($LASTEXITCODE -ne 0) { throw 'Queue recovery checks failed' }
& $python -m unittest discover -s tests -p test_disorder_continuation.py
if ($LASTEXITCODE -ne 0) { throw 'Checkpoint continuation checks failed' }
& $python scripts/run_disorder_queue.py prepare --queue results/disorder_queue_recovery_V9_beta052_beta039_20261101/queue.json --after results/disorder_beta060_20261027/extensions/V9_20261031/scan.json --scan results/disorder_beta052_20261028/extensions/V4p5_V5_V6_20261031/scan.json --scan results/disorder_beta039_20261031/scan.json --recover-missing-dependency --recovery-seconds 43200 --recovery-workers 2
if ($LASTEXITCODE -ne 0) { throw 'Recovery queue preparation failed' }
$queueDir = Join-Path $repo 'results\disorder_queue_recovery_V9_beta052_beta039_20261101'
$controller = Join-Path $repo 'scripts\run_disorder_queue.py'
$queueFile = Join-Path $queueDir 'queue.json'
$queueArgs = '"{0}" run --queue "{1}"' -f $controller,$queueFile
$job = Start-Process -FilePath $python -ArgumentList $queueArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $queueDir 'controller.stdout.log') -RedirectStandardError (Join-Path $queueDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $queueDir 'controller.pid')
```

Use this active queue path for future checks:

```powershell
& $python scripts/run_disorder_queue.py status --queue results/disorder_queue_recovery_V9_beta052_beta039_20261101/queue.json
Get-Content results/disorder_queue_recovery_V9_beta052_beta039_20261101/controller.stdout.log -Tail 30
Get-Content results/disorder_queue_recovery_V9_beta052_beta039_20261101/controller.stderr.log -Tail 30
& $python scripts/analyze_disorder_response.py
```

The notebook remains `notebooks/disorder_temperature_comparison.ipynb`, with
45 embedded figures across five temperatures. Automatic updates continue after
each recovered or newly measured point. Long open worms have no intermediate
closed-configuration checkpoint; the recovery controller must wait for them to
exit before it can safely use their final checkpoints.

Active replacement controller PID: 11104; started 2026-11-01T15:12:03.033052+00:00.


## Live recovery progress - late 1 November 2026

The original long disorder 9 worms closed and their processes exited. Recovery
started at 21:38 GMT from retained checkpoints. Both active chains are now in
production, beyond the 50,000-excursion warmup. Repeated live counters advanced:
- L64_V9_r0001: 236,759 -> 242,581 excursions.
- L128_V9_r0000: 92,760 -> 94,366 excursions.

The third continuation, L128_V9/r0004, is queued with 48,728 saved excursions.
The active runs continue because of their 12-hour continuation budgets, whose
nominal limits are around 09:38 GMT on 2 November, rather than because warmup
remains incomplete. Sparse checkpoint settings hide this production from disk
status until save/exit. Subsequent beta=0.52 and beta=0.39 sweeps remain queued.


## Parallel progress and current controller - 2 November 2026

The disorder 9 recovery reached 23 of 24 completed runs. The last chain,
L128_V9/r0004, passed warmup (115,772 live excursions on the check) and remains
in production. Waiting for its full 12-hour budget had left the later sweeps
idle. The schedule now runs the beta=0.52 additions alongside that final chain,
then beta=0.39. Scientific parameters, measurements, disorder maps, seeds and
per-chain run budgets are unchanged.

The new controller retained a Windows handle to the original process, checked
its creation time, and signalled readiness before the former Python controller
was stopped. Only that Python controller (11104) was stopped. QMC PID 43900 was
never stopped or restarted. The commit transfers metadata ownership; Windows
still displays its original parent PID, which is expected. Its actual exit code
will be read from the retained handle before validating the saved HDF5 result.
The old queue is marked superseded and backed up as
`queue.before_parallel_20261102.json`.

Each future ensemble still contains eight realizations, with seven workers per
ensemble and two ensembles running at once. Thus up to 14 new chains run beside
the one recovery chain: at most 15 QMC processes, below the previous limit of16.
Manifest worker metadata was backed up as
`manifest.before_parallel_20261102.json`. No realization is dropped and no
unfinished map is included in an average. The prepared INIs and potential files
were not changed.

Shared figure and notebook exports now acquire a cross-process file lock,
preventing concurrent point completions from overwriting common notebook files.
The lock is released when the exporting process exits, including abnormal exits.
Five handoff/worker-preservation tests and one cross-process locking test passed.
An actual shared-file lock was checked before the handoff.

The current plan is `results/disorder_parallel_tail_20261102/plan.json`. These
commands describe its preparation and launch; they have already been executed.
Do not repeat them or relaunch the superseded controller. For a different handoff,
obtain the actual active QMC PID and use a fresh plan directory.

```powershell
$repo = '\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
$python = 'C:\Users\kuo\AppData\Local\miniforge3\python.exe'
& $python -m unittest discover -s tests -p test_parallel_disorder_tail.py
if ($LASTEXITCODE -ne 0) { throw 'Parallel handoff tests failed' }
& $python -m unittest discover -s tests -p test_disorder_analysis_lock.py
if ($LASTEXITCODE -ne 0) { throw 'Shared analysis lock test failed' }
& $python scripts/run_parallel_disorder_tail.py prepare --plan results/disorder_parallel_tail_20261102/plan.json --old-queue results/disorder_queue_recovery_V9_beta052_beta039_20261101/queue.json --pid 43900 --scan results/disorder_beta052_20261028/extensions/V4p5_V5_V6_20261031/scan.json --scan results/disorder_beta039_20261031/scan.json
if ($LASTEXITCODE -ne 0) { throw 'Parallel tail preparation failed' }
$planDir = Join-Path $repo 'results\disorder_parallel_tail_20261102'
$controller = Join-Path $repo 'scripts\run_parallel_disorder_tail.py'
$planFile = Join-Path $planDir 'plan.json'
$planArgs = '"{0}" run --plan "{1}"' -f $controller,$planFile
$job = Start-Process -FilePath $python -ArgumentList $planArgs -WorkingDirectory $repo -WindowStyle Hidden -RedirectStandardOutput (Join-Path $planDir 'controller.stdout.log') -RedirectStandardError (Join-Path $planDir 'controller.stderr.log') -PassThru
$job.Id | Set-Content -LiteralPath (Join-Path $planDir 'controller.pid')
```

The controller initially waits in `awaiting_handoff`. The handoff was committed
only after confirming `adoption_ready_utc`, the new controller's liveness, the
old controller's `recovering_dependency` state, and that its only work child was
QMC PID43900. These are the remaining commands used after those checks. They
are historical, already completed; the old queue is now superseded so the guard
rejects a repeat.

```powershell
$oldPath = 'results/disorder_queue_recovery_V9_beta052_beta039_20261101/queue.json'
$old = Get-Content -LiteralPath $oldPath -Raw | ConvertFrom-Json
$ready = Get-Content -LiteralPath $planFile -Raw | ConvertFrom-Json
if ($old.status -ne 'recovering_dependency' -or $ready.status -ne 'awaiting_handoff') { throw 'Handoff states changed' }
Get-Process -Id $ready.controller_pid -ErrorAction Stop | Out-Null
$children = @(Get-CimInstance Win32_Process -Filter ('ParentProcessId={0}' -f $old.controller_pid) | Where-Object { $_.Name -ne 'conhost.exe' })
if ($children.Count -ne 1 -or $children[0].ProcessId -ne $ready.adopted.pid -or $children[0].Name -ne 'qmc_worm_canonical.exe') { throw 'Unexpected work children' }
Stop-Process -Id $old.controller_pid -ErrorAction Stop
@'
from pathlib import Path
import sys
sys.path.insert(0,'scripts')
from run_disorder_ensemble import _write_json
from run_disorder_queue import read,process_alive,now
old_path=Path('results/disorder_queue_recovery_V9_beta052_beta039_20261101/queue.json')
old=read(old_path);plan=read('results/disorder_parallel_tail_20261102/plan.json')
assert not process_alive(old['controller_pid']) and process_alive(plan['adopted']['pid'])
_write_json(old_path.with_name('queue.before_parallel_20261102.json'),old)
old.update(status='superseded',replacement_queue='results/disorder_parallel_tail_20261102/plan.json',finished_utc=now())
_write_json(old_path,old)
Path('results/disorder_parallel_tail_20261102/handoff.commit').write_text(now()+'\n')
'@ | & $python
```

Use these read-only checks for the current controller:

```powershell
& $python scripts/run_parallel_disorder_tail.py status --plan results/disorder_parallel_tail_20261102/plan.json
Get-Content results/disorder_parallel_tail_20261102/controller.stdout.log -Tail 30
Get-Content results/disorder_parallel_tail_20261102/controller.stderr.log -Tail 30
```

Notebook refresh remains `python scripts/analyze_disorder_response.py`. Open
`notebooks/disorder_temperature_comparison.ipynb` for completed data. The first
beta=0.52 extra chains started at 10:43:55 GMT on 2 November. Beta=0.39 follows
automatically after the beta=0.52 collection; it no longer waits directly for
the last disorder9 recovery chain.

Active parallel controller PID: 54224; handoff committed 2026-11-02T10:43:53.302790+00:00.

Local scan table/notebook writes also take the shared export lock, then release
it before launching the paired exporter. Two additional integration tests passed
for release-before-child-export and release after an analysis error. This protects
the standalone beta=0.39 analysis directory when controllers finish concurrently.


## L=128, beta=0.60, disorder 9: equilibration audit (2 November 2026)

All eight original processes finished, including the two checkpoint recoveries.
This is collection completion, **not validated equilibration**. The original
50,000-excursion warmup has not been demonstrated sufficient. Susceptibility
relative Monte Carlo errors range from about 1.3% to 42% across the eight maps;
these errors and aggregate logarithmic-bin summaries do not establish stationarity.
Do not infer thermal drift by comparing different disorder maps.

`parameter_files/disorder_quality_flags.json` marks this point provisional.
The existing comparison notebook and beta=0.60 figure titles display that warning.
Original results remain visible for review and are preserved unchanged. Audit
results are kept separately and never silently replace the original point.

The audit starts from immutable copies of the latest eight checkpoints, retaining
each map, canonical N=30000, beta, configuration and RNG state. For each map it:

1. Discards 200,000 additional closed excursions (plus ten unused production excursions).
2. Resets all statistics once and collects window A: 200,000 production excursions.
3. Resets all statistics once again and collects window B: 200,000 production excursions.

Runs checkpoint in one-hour chunks; statistics are retained between chunks of
the same window. The solver checks time limits between closed worms, so a chunk
can exceed one hour. A 48-chunk limit per stage stops an unproductive audit for
inspection; it does not label it equilibrated. Only run controls are changed in
checkpoint copies. Source checkpoints and original production accumulators are
never overwritten. Two audit chains run concurrently while the existing sweep
pipeline is pending; the audit permits eight only after that entire pipeline is
successfully completed. Main sweeps retain their existing fourteen-worker ceiling.

These are the exact preparation and execution steps, using the already built
canonical executable described earlier in this guide. Preparation refuses to
reuse an existing audit directory; the commands below document the audit already
launched, and should not be repeated against that directory.

```powershell
$repo='\\alfs1.physics.ox.ac.uk\al\kuo\projects\duplicatesGit\worm'
Set-Location -LiteralPath $repo
$python='C:\Users\kuo\AppData\Local\miniforge3\python.exe'
& $python -m unittest discover -s tests -p test_disorder_equilibration_audit.py
& $python scripts/check_disorder_equilibration.py prepare --audit results/disorder_beta060_20261027/equilibration/L128_V9_20261102/audit.json --manifest results/disorder_beta060_20261027/extensions/V9_20261031/data/L128_V9/manifest.json --main-plan results/disorder_parallel_tail_20261102/plan.json
$auditRoot=$repo+'\results\disorder_beta060_20261027\equilibration\L128_V9_20261102'
$auditScript=$repo+'\scripts\check_disorder_equilibration.py'
$auditFile=$auditRoot+'\audit.json'
Start-Process -FilePath $python -ArgumentList @('-u',('"'+$auditScript+'"'),'run','--audit',('"'+$auditFile+'"')) -WorkingDirectory $repo -RedirectStandardOutput ($auditRoot+'\controller.log') -RedirectStandardError ($auditRoot+'\controller.err') -WindowStyle Hidden
```

Read-only progress checks and notebook regeneration:

```powershell
& $python scripts/check_disorder_equilibration.py status --audit results/disorder_beta060_20261027/equilibration/L128_V9_20261102/audit.json
Get-Content ($auditRoot+'\controller.err')
Get-Content ($auditRoot+'\r0000\discard\chunk0000\run.log') -Tail 20
& $python scripts/run_parallel_disorder_tail.py status --plan results/disorder_parallel_tail_20261102/plan.json
& $python scripts/analyze_disorder_response.py
```

Audit controller PID at launch: 40888, 22:51 GMT on 2 November. `audit.json`
records per-map stages, child PIDs, source hashes, chunk counters, resets and fresh
window statistics. `original_statistics_review.json` records the original scalar
susceptibility summaries. Code and test copies are in `provenance/`.
After all eight maps finish, `window_comparison.json` reports A and B means,
errors and per-map differences for susceptibility, winding, total energy and N.
The two windows are adjacent and correlated; they are not extra disorder samples.
Inspect paired changes and whether longer binning reaches stable errors.
Absence of detected drift does not prove equilibration. Independent initial-state
checks and longer windows may still be required before removing the provisional
flag. The audit never removes that flag automatically. A failed/interrupted audit
must be inspected before continuing; rerunning `run` refuses a partially used audit
rather than accidentally resetting a partial measurement window.


### Integration after extended thermalisation (3 November)

The additional discarded evolution continues for all eight maps. A background
integration process waits until every map completes both fresh windows, then
copies each final window B into the normal point directory:
`extensions/V9_20261031/data/L128_V9/r000i/extended_thermalisation_20261102/`.
It validates all eight maps and physical parameters before updating the existing
`manifest.json`. The old manifest is retained as
`manifest.before_extended_thermalisation_20261102.json`; original result files
and all audit windows remain unchanged. The ordinary scan analysis and comparison
notebook then read the eight fresh windows from their usual dataset structure.
The provisional equilibration flag remains until convergence is assessed.

The integration process was started once with the same hidden Start-Process
pattern above, substituting `scripts/integrate_extended_thermalisation.py` and
redirecting output to `integration.log` and errors to `integration.err` in the
audit directory. Its state is `integration.json`. It must not be launched twice.
Read-only status command:

```powershell
Get-Content results/disorder_beta060_20261027/equilibration/L128_V9_20261102/integration.json
```

The comparison notebook also displays the partial per-map audit at its end;
`scripts/watch_disorder_equilibration_notebook.py` refreshes that section as
stages finish. This section does not substitute a partial-map disorder average.


### Progress visibility repair (4 November 2026)

The simulations continued, but the notebook watcher stopped on a transient
network PermissionError when reading audit.json. It now retries these reads,
and refreshes on saved chunk progress as well as stage completion. The notebook
introduction shows per-map extra discarded, window A and window B counters.
Counters describe saved closed excursions; a currently running chunk can be
further ahead. The waiting integration controller also retries network reads.
Only these lightweight controllers were restarted; no QMC process was stopped.
Current watcher logs are notebook_watcher_20261104.log/.err; integration logs
are integration_20261104.log/.err, in the audit directory. Two progress-reporting
tests passed, including handling the counter reset before fresh production.
All requested beta=0.52 and beta=0.39 points are now complete. At 11:03 UTC,
four maps had completed both fresh windows, one was in window A, and three were
in extra warmup. Map r0005 had no new closed-state checkpoint since 05:02,
although its worker continued consuming a core: an open worm can overrun the
one-hour chunk limit. Its running state was preserved rather than restarted.
