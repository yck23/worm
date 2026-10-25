# Final-production baseline data

`baseline_summary.csv` contains all 27 (L, beta) points. `baseline_rank_observables.csv` retains all 520 rank summaries and source paths; `baseline_scalar_ranks.npz` contains the same scalar measurements and their separate counts. No checkpoint configurations are needed to inspect or reanalyse these files.

Each point uses only its matched final `prod03` production. `chi` replaces G(0) with exact canonical density N/L^2. `chi_scaled` divides chi by L^(7/4). `W2` is Wx^2+Wy^2; it is not divided by two.

The `_se` columns are one-standard-error estimates from deleting entire ranks and recomputing count-weighted means. Chi uses `g_count`; W2 uses `w_count`. The rank SD columns describe rank-to-rank spread and are not the plotted SEs. The jackknife degrees of freedom are recorded for context, not used to claim calibrated confidence intervals.

All uncertainties are conditional on independent, equilibrated rank streams. Neither a scalar blocking plateau nor equilibrium has been established. `coarse_histories` marks points with fewer than 20 G histograms in at least one rank; this is a descriptive warning threshold.

Precision flags compare SE/abs(mean) with 5%. The optional winding criterion uses SE <= max(0.05*abs(W2), 0.02), so an absolute SE target is available when W2 is small. These flags do not establish convergence.

Rebuild with `python scripts/prepare_response_dataset.py`. The compact input audits under `../audit/` are required; hashes and formulas are recorded in `../provenance/baseline_manifest.json`. Source checkpoint paths are historical provenance.
