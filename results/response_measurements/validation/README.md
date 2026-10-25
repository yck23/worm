# Scalar response measurement validation

The compact [report.json](report.json) records the successful end-to-end regression run. Its executable SHA-256 identifies the tested build. Paths are relative to the repository root.

Run the tests from the repository root with a Python installation containing NumPy and h5py:

```powershell
python scripts/test_response_measurements.py
```

The script uses `build-repro/worm-square-response/qmc_worm.exe` by default; `--executable` selects another compatible canonical serial build. It prepends the packaged DLL directory, supplies a writable ALPS temporary directory, and launches child processes with hidden windows on Windows. `--output` selects the validation artifact directory. Each invocation creates a new timestamped directory.

The tests use a seeded 3 x 3 canonical Bose-Hubbard system with N = 13, nmax = 5, and beta = 0.3. They check:

- Corrected equal-time susceptibility against the complete displacement-matrix sum, and scalar winding against the sum of both components, over identical measurement windows.
- Default-off preservation of the RNG trajectory and worldline configuration.
- Legacy checkpoints without the new parameters and explicit activation of response measurements on copied legacy checkpoints.
- Bitwise equality of ordered histories, RNG state, and worldline configuration after a restart that interrupts an incomplete response block.
- Correct accounting for partial-block winding and density sums.
- Statistics reset clearing both completed histories and partial sums.
- Saving and restoring empty histories before thermalisation ends.
- Rejection of incompatible interval changes, silent disabling of existing histories, history-cap overflow, and unsupported open boundaries.

The final run passed all nine groups in approximately 11 seconds. Scalar/vector differences were 3.02e-14 for susceptibility and 1.11e-16 for the winding sum. These are bookkeeping and reproducibility checks; they do not establish equilibrium, a converged autocorrelation error, or a BKT transition. The tested executable is canonical and serial; grand-canonical and MPI execution need their own runtime checks.

ALPS writes empty vector histories as HDF5 null dataspaces. Readers must map `dataset.shape is None` to an empty numeric array. The regression reader includes this handling.

Only this README and the compact report are intended for version control. Raw timestamped runs contain diagnostic INIs, logs, results, and checkpoints and remain local. The earlier failed report was a test-reader issue with null dataspaces; it did not identify a simulator failure. No binaries or checkpoint files are needed in the commit.

## Completed-family export verification

[family_verification.json](family_verification.json) records validation of all four completed MPI families without launching simulations. Pooled rank scalar counts and means agree with the latest merged outputs. Completed block sweeps plus retained partial sweeps also equal both checkpoint and merged particle-number measurement counts.

Re-export and validate completed families with:

```powershell
python scripts/collect_response_measurements.py --points l128t028 l128t030 l128t045 l64t070 --verify
```

The collector rejects checkpoint-save warnings even when the MPI process returns success. It checks merged counts and means before publishing a history export. [collector_checks.json](collector_checks.json) records six synthetic HDF5 regression groups covering empty histories, omitted zero-count accumulators, stale counts, inconsistent means, missing nonempty accumulators, and save warnings. Reproduce those checks with `python scripts/test_response_collection.py`; this test launches no simulations.
