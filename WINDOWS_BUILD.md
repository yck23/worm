# Reproducible Windows build and run

This repository uses a pinned ALPSCore 2.3.3 commit and a small MinGW compatibility patch. The build creates a bundled runtime package in `dist/windows-square`; Windows system libraries and the Microsoft MPI runtime remain host prerequisites.

## Prerequisites

- CMake in `C:\Program Files\CMake`
- MSYS2 UCRT64 GCC, Boost, Eigen and HDF5 under `C:\msys64\ucrt64`
- Microsoft MPI SDK and runtime
- Git

## Build from source

From PowerShell in the repository root:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\build_windows.ps1 -Jobs 32
```

The script performs these steps:

1. Clones ALPSCore if it is not already present.
2. Checks out commit `3606edfbabc44e0fa7d05efa2a50c6a8a481340f`.
3. Applies `patches/alpscore-mingw.patch`.
4. Configures ALPSCore from a fresh CMake cache in Release mode with `Testing=OFF`.
5. Builds and installs ALPSCore locally.
6. Configures WORM from fresh CMake caches with `LATTICE=square`, both normally and with the canonical-window option.
7. Builds serial and MPI executables.
8. Collects the required MinGW and ALPSCore DLLs into one package.
9. Writes SHA-256 hashes to `dist/windows-square/SHA256SUMS.txt`.

## Core selection on this desktop

The audited processor is a 13th Gen Intel Core i9-13900 with 24 physical cores and 32 logical processors.

- Use `-Jobs 32` while compiling.
- Start production simulations with `-Processes 24` (one MPI rank per physical core).
- Benchmark `-Processes 32` for long runs if maximum throughput matters; Hyper-Threading may or may not improve Monte Carlo throughput.
- Use `-Processes 1` while checking a new parameter set.

Each MPI process is a simulation replica whose statistics are combined by the ALPSCore MPI adapter.

## Verified state

The independent recheck was completed on 2 October 2026 for WORM commit `bd0a6d9`. The pinned ALPSCore source and compatibility patch were checked again, and all 32 package hashes matched. No rebuild or environment activation is required for parameter-file edits. The copy-and-edit execution procedure is in [GRAND_CANONICAL_RUN_GUIDE.md, Section 8](GRAND_CANONICAL_RUN_GUIDE.md#8-run-one-editable-parameter-file).

On 2026-09-28 the build was repeated after the final parameter, checkpoint, invariant, and no-data-output guards. ALPSCore was confirmed at commit `3606edfbabc44e0fa7d05efa2a50c6a8a481340f`; its only two modified files are exactly those reproduced by `patches/alpscore-mingw.patch`, as verified by a clean reverse-patch check. The package contains ordinary and canonical serial/MPI executables. All 32 manifest entries matched their SHA-256 hashes, there were no unlisted payload files, and every non-system import was bundled; MPI executables additionally require the host `msmpi.dll`.

The final build completed with no compiler errors. Its warnings are repeated GCC diagnostics in Boost.Variant/standard-library templates, two existing Release-build `linkdir_found` variables used only by assertions, and benign CMake discovery notices. Final two-rank smoke tests then preserved exact canonical N = 38 and produced a grand-canonical mean N = 194.473 with Var(N) = 1148.45 after starting from N = 155. Both sampled nonzero winding. A deliberately measurement-free run wrote parameters and a checkpoint without a NaN result group; a premature statistics reset was rejected. Inputs with fewer than ten sweeps were also rejected before the former zero progress-interval operation.

Earlier exact-diagonalisation regressions found grand-canonical and canonical 3 x 3 energies within 1.79 estimated standard errors or better. The canonical MPI executable completed the independent 24-rank 64 x 64 checkpoint/production audit in `RB87_BKT_AUDIT.md`. The grand-canonical response identity was checked independently, and a 156-block L = 16, 24, and 32 fixed-mu campaign resolved high- and low-temperature winding regimes. Residual number-sector drift near the crossover limits that campaign to a diagnostic, not a precision Tc or an automatic validation of L = 64 grand-canonical sampling.

ALPSCore's upstream test targets are disabled because some are POSIX-specific on this Windows/MinGW toolchain. Verification here is therefore an end-to-end WORM integration test, not the complete upstream ALPSCore unit-test suite or a scientific convergence benchmark.

## Parameter file

The generic `run_bose_square.ps1` helper intentionally overrides lattice dimensions, periodic boundaries, run duration, output filename, and checkpoint filename from its own command-line switches. Other model values such as `t_hop`, `U_on`, `mu`, `beta`, and `nmax` come from its chosen parameter file.

The dedicated `run_rb87_bkt64.ps1` runner is different: it does not override physical values in the selected Rb-87 INI. It selects the executable, run mode, process count, parameter-file path, and output directory. Fresh stores the INI parameters in checkpoints; Resume and Production restore those stored values.

To use a different physical setup, copy the INI file, edit it, and pass it with `-ParameterFile`.

For the Rb-87 workflows, `parameter_files/Rb87_BKT_64.ini` is exact-N and `parameter_files/Rb87_BKT_64_grand_canonical.ini` is mu-controlled. The latter must be run with `scripts/run_rb87_bkt64.ps1 -Ensemble GrandCanonical`. The independent-start sweep helper is `scripts/run_gc_bkt_audit.ps1`, and `scripts/analyze_gc_bkt_audit.py` combines only numbered production blocks. See `GRAND_CANONICAL_RUN_GUIDE.md` for the required burn-in, convergence tests, and BKT interpretation.

## Where ALPSCore is

The build helper keeps the dependency below this repository:

| Path | Meaning |
|---|---|
| `build-repro/alps-src` | Pinned ALPSCore Git checkout |
| `build-repro/alps-build` | Temporary ALPSCore build tree |
| `build-repro/alps-install` | Local ALPSCore installation used to link WORM |
| `patches/alpscore-mingw.patch` | Reproducible Windows compatibility changes tracked by this repository |
| `dist/windows-square` | WORM runtime package, except Windows and Microsoft MPI system libraries |

You do not need to version-control the nested ALPSCore checkout separately. The tracked build script records the exact upstream commit and reapplies the tracked patch. The build and install trees are reproducible local products and remain ignored by Git.

## Run

Default 8x8 serial run:

```powershell
.\scripts\run_bose_square.ps1
```

Short smoke test:

```powershell
.\scripts\run_bose_square.ps1 -Lx 4 -Ly 4 -Sweeps 200 -Thermalization 20 -RuntimeLimit 5
```

Four MPI replicas:

```powershell
.\scripts\run_bose_square.ps1 -Processes 4
```

Custom parameter file:

```powershell
.\scripts\run_bose_square.ps1 -ParameterFile .\my_bose_hubbard.ini -JobName my_scan_point
```

Results are written to `results/` by default.
