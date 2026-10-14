# High occupation cutoff: initialization memory fix

The uniform build previously allocated a dense four-dimensional hopping-weight
lookup table with `(nmax - nmin + 1)^4` doubles on every MPI rank. For the
Bose-Hubbard input `nmax = 1024`, this required 8.0313 TiB per rank, or about
192.75 TiB for 24 ranks. The resulting `std::bad_alloc` occurred during worm
construction, before simulation. Reducing the process count could not make
this particular cutoff practical.

Hopping weights now call the same model formula used to populate that table,
on demand. The two remaining occupation matrices require about 16.03 MiB per
rank at this cutoff, in addition to the simulation's other memory. The change
preserves the occupation cutoff, particle number, weight formulas, and
checkpoint format. It trades cached hopping lookups for formula evaluation;
production throughput has not been benchmarked.

## Validation on 13 October 2026

- Rebuilt and packaged square-lattice serial/MPI binaries for both ensembles.
- Four fixed-seed, 3,000-sweep serial comparisons against the previous package:
  canonical and grand-canonical Bose-Hubbard, XXZ spin 1/2, and XXZ spin 1.
  All 45 saved result datasets in each case matched exactly.
- A short 24-rank MPI run retained the physical parameters of
  `parameter_files/Rb87_L16_T028.ini`: L = 16, beta = 0.28, N = 30,000,
  and nmax = 1024. Only duration and diagnostic/measurement intervals were
  shortened in a separate test input. All ranks completed an update sweep,
  passed the configuration check, and saved checkpoints with exactly 30,000
  particles. The time limit was reached during thermalization; this was an
  initialization/update regression, not a production or convergence test.
- Test inputs, comparison scripts, logs, checkpoints, and previous executables
  are retained locally under `results/allocation_fix_validation/`.

The original input and `results/clean_Rb87/l16t028` outputs were not modified.
The rebuilt runtime is in `dist/windows-square`, with an updated checksum
manifest. The original Fresh command can be rerun.
