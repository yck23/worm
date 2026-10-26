# Committed clean-system results

For the clean **winding-sum and scaled-susceptibility dataset with rank error bars**, start with the new [response notebook](response_measurements/response_measurements.ipynb) and [collection notes](response_measurements/README.md). It includes 27 final-production points, 520 rank records and four supplemental scalar-history pilots. Its separate [checksum manifest](response_measurements/FILES_SHA256.csv) covers that package.

Start with [the main comparison notebook](cleanRbRDM/density_matrix_comparison.ipynb)
and [its collection and analysis notes](cleanRbRDM/README.md). The notebook contains
ordinary analysis and plotting cells, with saved outputs. Its final section evaluates
the equal-time susceptibility and its leading BKT size scaling.

The measurement definitions and the additional simulation work required for reliable
error bars are documented in [the main guide, Section 11.6](../COMPLETE_BKT_GUIDE.md#116-equal-time-susceptibility-and-reliable-error-bars).
The susceptibility plots in both notebooks now show final-production means with
matching whole-rank jackknife standard errors. Historical pooled-production CSVs
remain separate. These errors do not certify equilibrium or a precision BKT temperature.

## Included inputs and results

- 81 selected HDF5 production outputs at L = 32, 64 and 128, with nine beta values
  and three productions per point. These files contain the stored observables,
  marginal errors, counts and simulation parameters used by both notebooks.
- 57 outputs from `clean_Rb87` and 24 replacement outputs from `cleanRbRDM`.
  The notebook uses replacements whenever the corresponding recollection folder
  exists. Stale, missing or superseded density matrices are not selected.
- 24 new per-rank NPZ summaries and their measurement-count verification records;
  effective INI files and the original recollection manifest.
- The two current notebooks, current fit/winding CSV, susceptibility CSVs, and
  three rendered susceptibility plots. Earlier notebooks and scratch exports
  remain local archives.
- [analysis_inputs.csv](cleanRbRDM/analysis_inputs.csv), with portable relative paths,
  parameters and checksums for the selected HDF5 files, and
  [FILES_SHA256.csv](FILES_SHA256.csv), with sizes and hashes of the original published
  result snapshot except the checksum manifest itself. The response package has its own manifest.

The matched scan has t = 1, U = 0.152, V = 0, nmax = 5 and density
N/L^2 = 1.8310546875. N is 1875, 7500 and 30000 for the three sizes. This is a
separate dataset from the audited N = 9948, L = 64 reference in the main guide.
Folder tokens such as `t035` denote beta = 0.35, not temperature; plot temperatures
are T = 1/beta in hopping units. Cutoff, equilibration and blocking convergence
must be assessed before a precision physical claim.

## Reproduce the analysis

From the repository root, install the analysis packages in your preferred Python
environment, then open the main notebook:

```sh
python -m pip install numpy scipy matplotlib h5py jupyterlab
python -m jupyterlab results/cleanRbRDM/density_matrix_comparison.ipynb
```

Select that environment's kernel and run all cells in order. The notebooks read
HDF5 inputs without modifying them and regenerate their CSV outputs. Susceptibility
figures are embedded in the main notebook; PNG copies are also committed for
viewing without Jupyter. The recorded environment versions are:

- numpy 2.4.2
- scipy 1.17.1
- matplotlib 3.10.8
- h5py 3.15.1
- ipykernel 6.29.5

The original collector, `scripts/collect_rdm.py`, is included for provenance and
future local collection. Running it requires the original rank checkpoints,
Microsoft MPI, and a built executable; those are not required for notebook analysis.
Absolute paths in the old HDF5/JSON/INI provenance describe the original machine.
Use the manifest's relative paths for analysis in another checkout.

## Version-control boundary

Large `*.clone.h5*` restart files, builds, binaries, runtime logs, temporary files,
notebook backups and unrelated scans remain ignored. The selected result files
are explicitly tracked despite the general `/results/` ignore rule. Ordinary
`git add` records later changes to those tracked files; deliberately adding new
result files under this directory requires `git add -f` with explicit filenames.
