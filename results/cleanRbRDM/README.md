# cleanRbRDM — additional density-matrix collection

Started 20 October 2026. This folder collects new equal-time `Density_Matrix` G(r)
statistics for the points missing usable correlation errors in `clean_Rb87`.

## Completed collection and comparison notebook

Completed **20 October 2026, 17:00 BST**. All eight points finished three production
blocks: **24 files and 2,501 complete G(r) histograms**, excluding the pilot. Every
stored density-matrix mean and error is finite; the production arrays are distinct
and all production stderr logs are empty. Original source checkpoint sizes and
timestamps match the recorded provenance.

Open [density_matrix_comparison.ipynb](density_matrix_comparison.ipynb) and read down
from the top. The simplified version has ordinary NumPy, SciPy and Matplotlib lines,
no custom function definitions, no analysis-module imports, and one graph per cell:

1. Read the three productions and average them using their measurement counts.
2. Show the L32, L64 and L128 density matrices separately at an editable beta.
3. Show G(r) at every temperature, then compare the sizes at one temperature.
4. Fit one selected point to algebraic and exponential decay, with each fit visible.
5. Calculate chi-squared, degrees of freedom, reduced chi-squared and nominal p-values explicitly.
6. Repeat the fits across temperatures and plot the parameters and statistics.
7. Plot Wx squared, Wy squared and their sum for all sizes and temperatures, showing every selected production as a faint dot.
8. Calculate the equal-time susceptibility from the complete displacement matrix, replacing its biased on-site entry by the exact density.
9. Plot susceptibility, susceptibility divided by L^(7/4), and its system-size dependence; save the 27 means and 81 production values.

Both notebooks have been run and contain their plots. Edit `beta` for the matrix
plots, or `L` and `beta` at the individual-fit example. Every temperature axis is
T=1/beta. Fits use the same default interval as the previous analysis. Propagated
errors and chi-squared neglect covariance; nominal p-values are approximate.

The short [exponent_vs_temperature.ipynb](exponent_vs_temperature.ipynb) independently
reads the raw files and shows exponent, decay length, chi-squared and reduced
chi-squared against temperature with the same direct lines. It has no dependency
on exports from the main notebook.

The main notebook saves all 27 points, both fits and winding means together in
`analysis/simple_density_fits_and_winding.csv`. The previous notebooks are preserved
as `*.before_simplification_20261023.ipynb`; older figures and tables remain available.

These checks establish successful collection and file consistency. The highest-beta
L128 productions have only 4-5 histograms per rank and the output has one error-bin
level, so finite errors alone do not establish convergence or fit precision.

## Susceptibility and statistical interpretation - 23 October 2026

The final notebook section computes `chi_s = G.sum() - G[0, 0] + N / L**2`
for each production, then averages using the density-histogram counts. This is
the lattice version of the equal-time susceptibility in
[Carleo et al.](https://arxiv.org/html/1305.3032v3), or zero-momentum occupation.
No extra beta, site-count division, or normalization by G(0) is applied.

The main guide's [Section 11.6](../../COMPLETE_BKT_GUIDE.md#116-equal-time-susceptibility-and-reliable-error-bars)
explains the observables, measurement parameters, scalar-estimator implementation
still needed, and a production protocol for reliable uncertainties. The current
notebook provides central estimates and production dots, not validated susceptibility
error bars or a fitted transition temperature. The density matrices use nmax = 5;
this clean scan has not inherited cutoff convergence from the separate nmax = 16 audit.

Current outputs are:

- [susceptibility_summary.csv](analysis/susceptibility_summary.csv): 27 count-weighted means.
- [susceptibility_productions.csv](analysis/susceptibility_productions.csv): 81 productions,
  relative source paths, counts, raw on-site values and corrected susceptibility.
- [Susceptibility versus temperature](analysis/susceptibility_vs_temperature.png).
- [Scaled susceptibility versus temperature](analysis/susceptibility_scaled_vs_temperature.png).
- [Susceptibility versus size](analysis/susceptibility_vs_size.png).

The 24 `prod*.rank_means.npz` files at the eight recollected points retain eight
rank means and their counts per production. Sum each rank's whole matrix before
estimating errors; combine its three consecutive windows before a rank jackknife.
Those eight chains can provide an exploratory error estimate if equilibrated, but
there is no time series with which to demonstrate a blocking plateau. The other
19 points do not have these saved rank summaries.

## Reproducing the committed analysis

The committed selection is described in [results/README.md](../README.md).
It includes every raw production file used by both current notebooks, the new
rank summaries, parameter/provenance records, and the current output tables and
susceptibility figures. `analysis_inputs.csv` lists all 81 selected HDF5 paths,
physical parameters, counts and SHA-256 hashes. `../FILES_SHA256.csv` inventories
the complete published result selection. Notebook paths are resolved relative to
the checkout; absolute paths embedded in older provenance record the original run.

Open `density_matrix_comparison.ipynb` and run all cells from top to bottom.
The required analysis packages are NumPy, SciPy, Matplotlib and h5py, with a Jupyter
Python kernel. The companion `exponent_vs_temperature.ipynb` is also self-contained.
The first cells locate the repository from the current directory or its parents.

The older notebooks, scratch previews, older analysis subdirectories described
below, and large restart checkpoints remain local archives rather than part of
this committed selection. The collection commands below require those local
checkpoints and the built Windows executable; a clone alone reproduces the
analysis, not the original Monte Carlo continuation.

## Queue

In order: L128 beta=0.35, 0.39, 0.45, 0.52, 0.60, 0.70; L64 beta=0.60, 0.70.
The `t` directory token represents beta, not temperature.

- Eight MPI ranks, one point at a time, using the packaged canonical executable.
- `Nmeasure=1`, `Nmeasure2=1000`: one G histogram per 1,000 measurement sweeps per rank.
- First point: a 180-second setup check, saved separately as `rdm_pilot.h5`.
- Each point: three 900-second production blocks, `rdmprod01.h5` through `rdmprod03.h5`.
- Bounded queue: approximately six hours of production plus setup, pilot, checkpoint
  copying/loading/saving, and scheduler overhead. Completion time depends on the machine.

## Starting states and data integrity

The controller repacks copies of live configuration and RNG data from ranks 0–7 of
existing `clean_Rb87` checkpoints. All these source ranks have passed their configured
thermalization threshold. This is a warm continuation with fresh measurements, not
an independent-seed equilibration study. The source archives are opened read-only;
physical parameters and original files are preserved. Source paths, sizes, timestamps,
and the executable hash are recorded in `manifest.json`.

Every copied checkpoint starts with an empty measurement group. Before each block,
the controller removes measurement groups only in its own copied checkpoints and
uses `--reset-statistics`. Each output filename must be absent before launch. These
steps prevent ALPS append-mode saving from retaining old correlation groups. The
measurement interval and filenames are changed in both the typed HDF5 parameter
datasets and serialized raw INI attributes. ALPS reparses the raw values on restore;
keeping both synchronized ensures the constructor uses the new interval. INI changes alone would not suffice.

Eight ranks limit memory use and checkpoint I/O on the desktop. Original rank offsets
and RNG streams are retained. This subset can be continued with eight ranks in this
folder; do not mix in the other original rank files.

## Inspecting progress

From the project root, with the Conda Python environment used by the notebook:

```powershell
& 'C:\Users\kuo\AppData\Local\miniforge3\python.exe' scripts/collect_rdm.py --status
```

`active.json` records the active point, block, controller/MPI PIDs and start time.
`status.jsonl` includes 30-second heartbeats, completed-block histogram counts, and
error availability. Each point has stdout/stderr logs, effective parameter INIs,
new checkpoints, and completed HDF5 results. `*.metrics.json` verifies each rank's
histogram count equals `floor(diagonal_count/1000)` and the merged count equals
the rank sum. `*.rank_means.npz` retains per-rank G means and counts for each block.
These are block summaries, not the individual 1,000-sweep histogram time series.

After each point's three blocks finish, `analysis/gr` is refreshed with the audit,
per-production axis and shell curves, count-weighted means, raw arrays, and plots.
Queued or missing results remain visible. Production SD is descriptive; finite site
errors alone do not establish convergence or a blocking plateau. The pilot is a
separate diagnostic and is excluded from the three-production collation.

To finish the current block and stop before starting another:

```powershell
& 'C:\Users\kuo\AppData\Local\miniforge3\python.exe' scripts/collect_rdm.py --stop
```

The collection is run by `scripts/collect_rdm.py`. Do not launch a duplicate controller.
A stopped or failed queue should be inspected before restarting; checkpoint states
and existing output files are deliberately protected against silent overwrites.

## Earlier unweighted decay exponent analysis (archived)

The earlier notebook fitted an effective power-law exponent for L32, L64,
and L128 and plots all three against T=1/beta. New measurements are preferred; valid old productions
fill the remaining temperatures, giving nine temperatures at each size. L32 uses
all 27 existing valid production files from clean_Rb87 (three per temperature),
containing 8,718 complete histogram records. Old/new sources use different markers
and are not pooled at the same point.

That version used an unweighted log-log fit on max(3,L/16) <= r <= L/4:
r=3..8 for L32, r=4..16 for L64, and r=8..32 for L128. The absolute short-distance
cutoff therefore overrides L/16 for L32. Controls include axis/shell choice and
fit-window bounds. A second plot compares upper cutoffs L/8, L/4, and L/3.
The L32 L/8 test contains only two distances and is flagged insufficient rather
than fitted; valid comparisons require at least five distances. Another shows the fitted curve
and interval for a selectable beta. Tables are exported to `analysis/notebook`:
`decay_exponents_vs_temperature.csv`, `decay_exponents_by_production.csv`, and
`decay_exponents_by_fit_window.csv`.

The previous G(r) shading was the pooled curve plus/minus the unweighted sample SD
of three production curves (ddof=1); normalized curves used the SD of each production
normalized by its own G(1). This was descriptive spread, not a propagated HDF5 error
or an error on the pooled mean, and was not divided by sqrt(3). The simplified notebook uses no shading. The earlier optional exponent bars
showed production-fit SD, not a
confidence interval or formal regression uncertainty. Effective slopes at curved
or exponentially decaying points should not be interpreted as established algebraic
exponents or used alone to infer a transition temperature.

The earlier unweighted exponent tables remain as reference outputs. The focused
notebook now uses the new weighted two-model comparison described below.

## Algebraic versus exponential comparison - 23 October 2026

Both notebooks now include the requested two-model fit comparison at every one of
27 lattice/temperature points (L32, L64, L128; nine temperatures each):

- Algebraic: G(r) = A * r^(-eta).
- Exponential: G(r) = A * exp(-r/xi).
- Two fitted parameters per model, no constant offset, identical points and errors.
- Fits minimize squared residuals in raw G space with fixed absolute sigma values.
- Default r ranges: 3..8, 4..16, 8..32 for L32, L64, L128.

For each production, the spatial mean error is approximated by
sqrt(sum(stored_marginal_error_j^2)) / number_of_unique_directions. Production means
are pooled with complete-histogram-count weights w; the pooled sigma is approximated
by sqrt(sum(w_p^2 * sigma_p^2)). The resulting chi2_diag is sum(((G-model)/sigma)^2),
with dof=n_points-2 and reduced_chi2_diag=chi2_diag/dof. Errors are not rescaled to
make the reduced statistic one. This differs from the previous descriptive SD of
production means, which is not used as the fit-error model.

The full spatial, inter-distance, and inter-production covariances are unavailable.
The above propagation and fits neglect them. Consequently p-values are explicitly
nominal under independent Gaussian errors, and formal_p_value_valid is False for
every fit. Smaller diagonal chi-squared measures relative weighted misfit; it does
not identify a phase or make a poor absolute fit adequate. The models are nonnested,
so the difference is not assigned a likelihood-ratio chi-squared p-value. Stored
histogram errors also lack a demonstrated blocking plateau.

Both simplified notebooks recompute the fits directly from raw HDF5 arrays.
There are no MODEL_* controls or cached-table dependencies. The main notebook
shows the two formulas, each fit, the chi-squared sums and the plotting lines.
The earlier module-based implementation and its exported diagnostics are retained
for reference, with the same default-window numerical results.

The earlier detailed diagnostics remain under analysis/model_comparison:

- model_fits.csv: 54 fitted models, parameters, nominal parameter errors, chi2,
  reduced chi2, degrees of freedom, nominal p-values, and validity flags.
- model_comparison.csv: 27 paired comparisons, including exp-minus-algebraic chi2.
- model_fit_residuals.csv: every fitted distance, observed G, assumed sigma,
  model prediction, standardized residual and chi2 contribution.
- model_fit_windows.csv: both models under three upper cutoffs. The L32 L/8
  windows have only two distances, so all 18 such model/window entries are flagged
  insufficient_points and excluded from fitting and plots.
- model_fit_profiles.npz and model_fit_settings.json: pooled data/errors and the
  exact settings, model definitions, covariance assumptions, and references.
- PNGs: eta, xi and both amplitudes versus temperature; chi2 and reduced chi2
  versus temperature for all sizes; nominal p-values; delta chi2; window sensitivity;
  and nine-temperature fit galleries for each lattice size.

Default-window result: algebraic has lower diagonal weighted misfit at 26 points;
exponential is lower at L128 beta=0.28 (T=3.5714). At that latter point both reduced
statistics are large (about 230 and 61), so the relative preference is not a claim
that the exponential model passes an absolute goodness-of-fit test.

Method references:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chi2.html


## Susceptibility error bars in the main notebook

The susceptibility section of `density_matrix_comparison.ipynb` now computes whole-rank jackknife errors directly and uses `errorbar` for raw susceptibility, scaled susceptibility and the size dependence. These plots use matched final-prod03 means and one-standard-error bars. A percentage-error plot shows all 27 susceptibility uncertainties even when the physical-scale bars are small.

The calculation reads the compact rank summaries in `../response_measurements/data/baseline_rank_observables.csv`, sums off-site contributions inside each rank, adds exact density, and count-weights the rank means. No custom functions or analysis-module imports are used. The direct result is checked against the actual final-production density matrix at every point.

`analysis/susceptibility_final_with_errors.csv` records the 27 plotted means and errors, including their L^(7/4) scaling. The older `susceptibility_summary.csv` remains the historical pooled-three-production table; it is not the central value of the new error-bar plots. All errors remain conditional on independent, equilibrated rank streams. The separate `../response_measurements/response_measurements.ipynb` includes winding, rank/temporal errors and the new continuation histories.

## Squared winding error bars

The winding section in `density_matrix_comparison.ipynb` plots final-prod03 means of Wx^2, Wy^2 and their sum, each with +/- one whole-rank jackknife standard error. Winding measurement counts weight each rank; summing components before resampling preserves their covariance. These means replace the pooled-three-production means in the winding plots to match the error window. `analysis/winding_final_with_errors.csv` contains all 27 matching means and errors. The historical pooled CSV remains unchanged. Valid errors require independent, equilibrated rank streams.
