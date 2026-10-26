# Winding and scaled susceptibility measurements

Open [response_measurements.ipynb](response_measurements.ipynb). It contains plain analysis cells, saved tables, and five figures: the full 27-point scan, transition-region view, baseline/pilot comparison, rank traces, and temporal blocking curves. It needs NumPy, Matplotlib and IPython; it does not need the large restart files or pandas.

The [main guide, Section 11.6](../../COMPLETE_BKT_GUIDE.md#116-equal-time-susceptibility-and-reliable-error-bars) explains the physics, estimator, cadence, resampling and convergence requirements.

## Dataset and observables

The baseline covers L = 32, 64, 128 (L² sites), beta = .28, .30, .32, .35, .39, .45, .52, .60, .70. Temperature is T/t = 1/beta with t = 1. All points have U/t = .152, V = 0, nmax = 5, canonical window .5 and density N/L² = 1.8310546875; N = 1875, 7500, 30000. This differs from the older N = 9948, nmax = 16 reference scan. Its cutoff convergence does not transfer to these data.

- `W2 = <Wx² + Wy²>` is the sum of squared winding components, without division by two. For the isotropic square system beta times stiffness is W2/2. The BKT infinite-size reference for W2 is 4/pi.
- `chi = N/L² + sum(G[r] for r != 0)` is the equal-time zero-momentum occupation. Exact density replaces the known biased on-site estimator. There is no extra beta, volume factor or radial multiplicity.
- `chi_scaled = chi/L^1.75` is the leading BKT scaling used in [Carleo et al.](https://arxiv.org/html/1305.3032v3). Finite-size and logarithmic corrections still matter.

[baseline_summary.csv](data/baseline_summary.csv) contains 27 means and one-standard-error estimates. [baseline_rank_observables.csv](data/baseline_rank_observables.csv) contains all 520 input rank records. Each family reproduces its selected final `prod03.h5` vector means and counts: 24 ranks at 19 original points, eight ranks at eight recollected points. The notebook recomputes every mean and error from these compact records.

Susceptibility is weighted by completed density-histogram counts; winding is weighted by its own diagonal measurement counts. Each delete-one-rank jackknife recomputes the weighted mean. Summing inside a rank preserves the covariance of spatial bins and winding directions. Whole ranks account for within-chain dependence without pretending individual intervals are independent. The validity still requires independently sampled, equilibrated rank streams. Eight ranks give a coarse uncertainty estimate; shared bias remains invisible.

This baseline is **final production only**. The earlier comparison notebook pools three production windows. These errors must not be attached to its different central values. All 27 baseline chi errors are below 5%; 26 winding errors are below 5%. The remaining near-zero W2 point meets the absolute SE target .02. These are precision statements, not proof of equilibrium; all baseline rows retain `conditional_unverified` and no established temporal plateau.

## Additional simulations actually run

Four eight-rank warm continuations initially ran for 600 seconds each on 25 October 2026. On 26 October, the L128 beta=.30 and .45 families each received another 1200 seconds concurrently, without resetting their statistics. These families now have 1800 seconds each; the other two retain 600 seconds. All four cumulative histories meet the defined numerical precision and blocking screens.

The cumulative 205767 completed scalar intervals cover 2057670 eligible updates, plus 148 partial-interval updates retained only in checkpoints. Original source checkpoints were read only, and every previously completed history record was preserved exactly. These are continued source RNG lineages, so they remain separate from the baseline.

| L | beta | Requested seconds | Completed intervals | Min intervals/rank | W2 +/- rank SE | chi/L^1.75 +/- rank SE | W2 blocking relative SE |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 128 | 0.28 | 600 | 109779 | 12559 | 0.01006 +/- 0.00210 | 0.25201 +/- 0.00497 | 20.0% (absolute target) |
| 128 | 0.30 | 1800 | 73287 | 8867 | 0.38384 +/- 0.01141 | 0.99939 +/- 0.01126 | 3.4% |
| 128 | 0.45 | 1800 | 13402 | 1477 | 2.56715 +/- 0.06022 | 3.46732 +/- 0.08436 | 2.7% |
| 64 | 0.70 | 600 | 9299 | 1083 | 4.66855 +/- 0.10204 | 3.98947 +/- 0.06981 | 2.5% |

All four pass the numerical screens, including both rank and blocking errors against the stated precision target. These remain **conditional sampling results**: independent-start/equilibration and nmax convergence review are required before a precision physical claim. No transition temperature is fitted. The L128 beta=.45 late-minus-early scaled susceptibility change is about 2.18 paired SE; the three-SE screening threshold does not prove absence of drift.

Before extension, beta=.30 missed the 5% winding target (blocked SE .02339, 6.0%) and beta=.45 lacked a demonstrated plateau. After extension, their winding blocking SEs are .01292 (3.4%) and .07036 (2.7%), respectively. The earlier run01 HDF5 snapshots remain archived; run02 outputs and history.npz at those two points are cumulative and must not be pooled with run01 as independent data.

The analyzer reblocks within each rank, then combines rank variances with squared count weights. It checks the three coarsest sufficiently populated blocking levels, requiring errors within 25%, at least 50 total blocks and at least four per rank, and the shortest of these levels to exceed four estimated autocorrelation times. These are preliminary numerical gates, not a certification. It separately requires at least eight ranks, 64 intervals per rank and reasonably balanced counts. The reported `_se` columns remain whole-rank jackknife errors; `_blocking_se` columns are the maximum error across the selected triple and are unavailable if that screen fails. Both must meet the precision target for the combined precision flag.

The displayed histories and baseline windows are correlated continuations. Visual differences cannot be tested by adding their errors in quadrature. Paired half-window diagnostics preserve rank pairing. No statistically independent chronology is constructed by concatenating different ranks.

## Run settings and sufficient errors

The initial pilot settings were:

```ini
Nmeasure = 1
Nmeasure2 = 100
Nresponse_measure = 10
Nresponse_max_blocks = 1000000
C_worm = 2.0
sweeps = 1000000000
thermalization = 50000
runtimelimit = 600
```

The two extensions used `runtimelimit=1200`; all other measurement settings were retained. The saved parameters.ini files describe the initial setup; per-invocation runtime and output selection are recorded in each manifest and HDF5 output.

`thermalization=50000` is inherited provenance, not a recommendation for a fresh start. These runs resumed warm states; the one-time statistics reset sets the counter at the threshold. Select fresh-run burn-in from multiple starts and stability of winding, coherence, density and energy. A short stable energy trace alone is insufficient.

`Nresponse_measure` controls the scalar recorder independently of `Nmeasure` and `Nmeasure2`. Ten means ten eligible closed update calls, not ten sweeps over every lattice site. Every interval includes zero-hit updates. Completed records pair the accumulated off-site susceptibility contribution with density and winding averaged over that same interval. Smaller intervals resolve correlations but do not increase the amount of independent information.

Keep `Nmeasure=1` for diagnostics. `Nmeasure2=100` or 1000 is a starting choice for spatial-array cost; it does not dictate the scalar time resolution. Preserve rank identity, seeds, chronology, actual update counts and partial intervals. `--reset-statistics` starts a new measurement window on the same chain, not a new independent simulation. Repeated cumulative outputs must not be counted repeatedly.

Use a relative SE target of 5% for chi and W2, or absolute W2 SE .02 near zero. After finding stable blocking and stationarity, estimate additional cost with `(current_SE/target_SE)^2`. Before a plateau, doubling the production length and checking again is more defensible than promising a fixed runtime. A practical stronger target is 50–100 long, approximately independent blocks **per rank**.

For B long blocks of k recorded intervals, each rank needs at least `B*k*Nresponse_measure` eligible updates after burn-in. Using the selected pilot block sizes gives conservative planning examples:

| Point | Selected coarse block size k | Updates/rank for 50 blocks | Current minimum completed updates/rank |
|---|---:|---:|---:|
| L128 beta 0.28 | 1024 | 512000 | 125590 |
| L128 beta 0.30 | 1024 | 512000 | 88670 |
| L128 beta 0.45 | 128 | 64000 | 14770 |
| L64 beta 0.70 | 128 | 64000 | 10830 |

All four now have a justified preliminary blocking scale under the numerical screen, but none meets this stronger 50-long-blocks-per-rank planning target at its selected coarse scale. These calculations are counts, not guaranteed wall times or proofs of decorrelation. Slow modes and common start bias can exceed estimates from a short history. Occupation-cutoff, continuum/lattice-spacing and finite-size uncertainties must be assessed separately.

## Reproduce the notebook and analysis

From the repository root, with NumPy, Matplotlib, IPython and h5py available:

```powershell
python scripts/prepare_response_dataset.py
python scripts/analyze_response_measurements.py --self-test
python scripts/analyze_response_measurements.py
python scripts/execute_response_notebook.py
```

Alternatively open the notebook in Jupyter and run all cells. The first script rebuilds the baseline from compact audit CSVs; the analyzer reads the four NPZ histories; the executor saves notebook outputs. The original source checkpoints are unnecessary for these commands. The notebook independently recomputes baseline/pilot means, errors and blocking curves without importing helper functions.

The implementation tests are:

```powershell
python scripts/test_response_measurements.py
python scripts/test_response_collection.py
```

The simulator tests require the built executable. Nine integration groups passed, including aligned scalar/vector identities, bitwise restart equivalence, partial intervals, resets, burn-in, legacy migration, default-off behavior and invalid-setting rejection. Six collection tests passed, including empty histories and rejecting stale or inconsistent checkpoints. All four production families also passed MPI merged-output versus checkpoint count/mean checks. See [validation](validation/README.md).

## Build and continue simulations locally

The response executable was built separately from the packaged release using the installed ALPSCore and MSYS2 dependencies:

```powershell
cmake -S V:/kuo/projects/worm/src -B V:/kuo/projects/worm/build-repro/worm-square-response -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DLATTICE=square -DCWINDOW=ON -DMATSUBARA=OFF -DALPSCore_DIR=V:/kuo/projects/worm/build-repro/alps-install/share/ALPSCore -DBOOST_ROOT=C:/msys64/ucrt64 -DHDF5_ROOT=C:/msys64/ucrt64 "-DCMAKE_CXX_FLAGS=-DALPS_UTILITY_NO_STACKTRACE -DALPS_NO_SIGNALS"
cmake --build V:/kuo/projects/worm/build-repro/worm-square-response --parallel 8
```

Adjust local dependency paths on another machine. Physics/measurement INI changes do not themselves require rebuilding. Source changes do. The collector currently targets Windows, Microsoft MPI and the repository's `dist/windows-square` DLLs.

The following extension was executed on 26 October 2026; repeating it would start another cumulative extension:

```powershell
python scripts/collect_response_measurements.py --points l128t030 l128t045 --run --seconds 1200 --jobs 2
python scripts/collect_response_measurements.py --points l128t028 l128t030 l128t045 l64t070 --verify
python scripts/analyze_response_measurements.py
python scripts/execute_response_notebook.py
```

Local restart checkpoints must be present for continuation. The first invocation in a new family uses copied source states and resets statistics; later invocations preserve and extend histories. Keep eight ranks and interval 10 in these families. Editing a nearby INI does not alter restored checkpoint parameters; the collector updates its own copied parameter records explicitly. It checks logs for checkpoint-save failures and compares exported rank counts/means to the latest merged HDF5 result before accepting a run.

`Nresponse_max_blocks` bounds retained history; reaching it fails explicitly. Archive a completed window and deliberately reset, preserving lineage and chronology for any later combined analysis, or increase the cap. Do not silently change the interval of saved histories. ALPS null-dataspace histories are normalized to empty numeric arrays by the exporter.

## Files and provenance

- `data/`: clean 27-point table and 520 scalar rank inputs, with schema notes.
- `audit/`: compact checkpoint-to-production audit and original drift screen.
- `provenance/`: source matching, hashes and baseline selection rules.
- `sampling/<point>/`: committed NPZ histories, parameters, manifests, history checks and merged production HDF5 outputs. The local multi-gigabyte restart files, temporary files and logs are excluded.
- `analysis/`: pilot summary, per-rank table and reblocking curves.
- `figures/`: standalone PNGs, also embedded in the notebook.
- `validation/`: compact passing reports; raw test runs remain local.
- `FILES_SHA256.csv`: checksums and sizes for this published result package, excluding itself.

Original absolute paths inside HDF5/INI files describe the collection machine. Analysis uses portable paths and the compact inputs in this folder. Old source archive paths in provenance need not exist to rerun the notebook.
