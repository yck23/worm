# Rb-87 2D BKT build and convergence audit

Audit date: 28 September 2026

## Bottom line

The completed independent audit passes the Windows build and the core energy/winding physics used for an exact-N, homogeneous 64 x 64 calculation. A wholly new base-seed family was rebuilt and run through fresh initialization, extended burn-in, and three separately reset production blocks. Its two late blocks agree and reproduce the earlier strict-window seed family within conservative rank-level uncertainty. The 40 nK, L = 64 point is therefore independently replicated as a diagnostic equilibrium point for energy and winding. It is still not a determination of the BKT transition temperature: that requires a temperature grid, several sizes, repeated seeds at the fit points, and logarithmic finite-size scaling.

The audit is therefore a qualified pass for energy and winding calculations, not an unconditional pass for every observable. In particular, the on-site element of the current one-body density-matrix estimator fails an exact operator identity in two independent seeds, and one optional Matsubara histogram is incomplete across checkpoint restore. Neither issue enters the packaged BKT winding result because the BKT analysis uses energy and winding and the audited binaries have Matsubara measurements disabled.

The chemical-potential-controlled path also passes its implementation, checkpoint, fluctuation-response, start-erasure, and winding-isotropy checks. A 156-block L = 16, 24, and 32 fixed-mu scan resolves normal and superfluid regimes and gives size-ordered raw crossover diagnostics. It does not pass a precision transition-temperature standard because number and winding retain material block motion near the crossover. The code and workflow are suitable to commit; the present scan is suitable as validation evidence, not as a publication-quality Tc.

The production configuration is:

```text
L = 64, N = 9948, U/t = 0.152, beta*t = 0.279077157506
nmax = 16, canonical_window = 0.5, 24 MPI ranks
```

Use [parameter_files/Rb87_BKT_64.ini](parameter_files/Rb87_BKT_64.ini) and [scripts/run_rb87_bkt64.ps1](scripts/run_rb87_bkt64.ps1) for the audited exact-N result. The simple instructions are in [BASIC_RUN_GUIDE.md](BASIC_RUN_GUIDE.md); the complete physics-to-analysis treatment is in [COMPLETE_BKT_GUIDE.md](COMPLETE_BKT_GUIDE.md). The smaller-size fixed-mu path and its convergence limits are isolated in [GRAND_CANONICAL_RUN_GUIDE.md](GRAND_CANONICAL_RUN_GUIDE.md); it is not yet validated at L = 64.

## Physical sources and mapping

The Oxford source is Abel Beregi's 2024 DPhil thesis, [Probing universality of 2D quantum systems with bilayer Bose gases](https://ora.ox.ac.uk/objects/uuid:b2f4f0a1-8576-4528-bbd3-557d273cfbdd). Appendix D.8 uses Rb-87, dimensionless interaction strength 0.076 for 1 kHz axial confinement, 40 nK, and a 0.5 micrometre numerical lattice spacing. The thesis also describes an approximately 32 micrometre uniform region and quotes a uniform-system critical phase-space density near 8.5.

The atomic mass is 1.44316089500 x 10^-25 kg from Daniel Steck's [Rubidium 87 D Line Data](https://steck.us/alkalidata/rubidium87numbers.pdf), revision 2.3.4.

For a continuum field discretised on spacing a, the mapping used here is:

```text
t = hbar^2 / (2 m a^2)
U/t = 2 g_tilde
```

With a = 0.5 micrometres and g_tilde = 0.076:

| Quantity | Value |
|---|---:|
| Physical width, 64a | 32 micrometres |
| t/h | 232.6009778 Hz |
| t/kB | 11.1630863 nK |
| U/t | 0.152 |
| T/t at 40 nK | 3.583238445 |
| beta*t at 40 nK | 0.279077157506 |

The filling was inferred from the weak-gas critical phase-space density D_c = ln(380/g_tilde) = 8.51719. At 40 nK this gives 2.42864 particles per grid site, hence N = 9948 on 64 x 64 after integer rounding.

This is an approximate homogeneous lattice discretisation. The program does not implement the thesis's radial trap, finite box wall, bilayer coupling, or disorder. At 40 nK the thermal wavelength is only about 1.87 grid sites, so continuum discretisation error is a real systematic uncertainty. A smaller spacing and continuum extrapolation would be needed for quantitative comparison with experiment.

## Build audit

The reproducible build script pins ALPSCore commit 3606edfbabc44e0fa7d05efa2a50c6a8a481340f and applies the recorded MinGW compatibility patch. It builds and packages four square-lattice programs:

```text
qmc_worm.exe
qmc_worm_mpi.exe
qmc_worm_canonical.exe
qmc_worm_canonical_mpi.exe
```

The canonical-window executable is required for the independently replicated L = 64 production results in this document. Early short grand-canonical L = 64 runs became trapped in starting particle-number sectors and produced zero winding; those runs are rejected as unequilibrated. A later, separate grand-canonical audit passed at L = 16 and tested finite-temperature behavior at L = 32, without claiming that the L = 64 grand-canonical case is thereby converged.

The final fresh-cache rebuild completed with GCC 14.2.0, HDF5 1.14.3, Microsoft MPI, and the pinned ALPSCore source. All 32 payload entries matched the generated SHA-256 manifest. The only unbundled DLL imports are Windows system/runtime libraries and msmpi.dll, which is supplied by the required Microsoft MPI runtime.

A fresh two-rank canonical smoke run then recorded requested N = 155, measured N = 155, window 0.5, and 15448 winding measurements. Checkpoint restore, reset-statistics production mode, result archiving, and restoration of the caller's PowerShell directory were also exercised. A deliberate canonical_window = 1.0 run exited with the intended validation error.

The code audit also found and corrected the following issues:

- Fresh exact-N states placed excess particles on consecutive row-major sites. Their locations are now shuffled reproducibly from the rank seed, avoiding an artificial density stripe.
- The canonical open-worm window was a hard-coded 0.1. It is now a parameter and must satisfy 0 < canonical_window < 1, which preserves exact N in closed configurations.
- Restoring a checkpoint now recomputes the potential-energy measurement estimator.
- Serial and MPI thermalisation timeouts with zero retained measurements now report that state cleanly and omit a misleading NaN result group.
- Both binaries now reject a parameter file for the wrong ensemble instead of silently ignoring an incompatible canonical value.
- Production reset now refuses to bypass an unfinished thermalisation threshold.
- The end-of-run invariant test independently reconstructs boundary-slice particle number and measured potential energy; the canonical build also rechecks exact N on every closed exit.
- Invalid beta, energy offset, worm weight, measurement interval, update probability, time tolerance scale, and Bose-Hubbard occupation cutoff values now fail early.
- The package now includes both ordinary and canonical serial/MPI executables and a SHA-256 manifest.

ALPSCore's complete upstream unit-test suite is not run on this MinGW build because it includes POSIX-specific targets. Validation is therefore based on the pinned source, package dependency/hash checks, parser guards, and end-to-end serial/MPI WORM jobs.

The 23-24 September independent rerun additionally established all of the following:

- clean Release rebuilds of the ordinary and canonical serial/MPI targets completed;
- ALPSCore was exactly commit `3606edfbabc44e0fa7d05efa2a50c6a8a481340f`, with only the recorded MinGW patch applied;
- all 32 packaged payload hashes matched `dist/windows-square/SHA256SUMS.txt`;
- fresh 1.2-million-sample grand-canonical and canonical 3 x 3 energy tests agreed with exact diagonalisation within 1.79 estimated standard errors or better;
- a fresh two-rank canonical MPI job, checkpoint continuation, and statistics reset all completed with exact particle number; and
- an independent parser reconstructed particle number, hopping-vertex count, event pairing, chronology, periodic continuity, and raw winding from every worldline event in all 24 new 64 x 64 rank checkpoints after every long stage. All checks passed, every rank had N = 9948, and every reconstructed value matched its stored counter.

The grand-canonical regression's energies and density-density correlations passed, but its on-site density-matrix element was 7.0 reported standard errors from exact diagonalisation. In the canonical N = 9 regression, all tested energy, density-density, and density-matrix values were within 1.79 standard errors. A separate canonical N = 13 identity test, described below, independently confirmed that the exceptional on-site estimator branch is biased.

### Optional chemical-potential-controlled extension

On 24 September the ordinary and canonical variants were cleanly rebuilt after adding `initial_particle_number` for fresh grand-canonical starts. This parameter creates a randomized initial state near a guessed N but applies no later constraint; `mu` remains the thermodynamic number control. The ensemble-aware runner requires `canonical = -1` with `-Ensemble GrandCanonical` and rejects mismatched parameter files. All 32 packaged hashes passed after the rebuild.

A two-rank L = 8 smoke calculation at `mu=-3.47` started from N = 155, changed number sectors, restored both rank checkpoints, and completed a reset-statistics production block. That block contained 17,238 measurements with mean N = 241.6857, Var(N) = 543.2207, and winding sum 3.26036. This establishes that the new path is genuinely grand canonical and checkpoint-complete for the tested observables.

It is not a convergence result: the first block mean was about 191.94, so the later value was still drifting strongly. The earlier 64 x 64 starts at occupations 1, 2, and 3 gave incompatible means of about 3569, 5936, and 7643 with zero winding. They show why changing N is not by itself a convergence test.

The follow-up used independent low/central/high starts and separately reset blocks at L = 16. A final count-weighted reanalysis pools first and second moments inside each family and gives the three families equal weight. At beta = 0.2790771575 and mu = -3.58, 15 family-block results gave mean N = 517.943 plus or minus 0.885, Var(N) = 7114.00 plus or minus 124.12, and winding sum 0.77072 plus or minus 0.00581. The low/mid/high family means for N were 518.35, 519.23, and 516.25.

At mu = -3.53, nine family-block results gave mean N = 612.968 plus or minus 2.829, Var(N) = 6949.16 plus or minus 560.11, and winding sum 1.30067 plus or minus 0.01586. The low/mid/high family means were 609.50, 610.83, and 618.57. The finite difference gave delta mean N / delta mu = 1900.50 plus or minus 59.28, while beta times the endpoint-average number variance was 1962.35 plus or minus 80.05, a ratio of 0.9685. This 3.15 percent difference, over a finite 0.05 mu interval, independently verifies quantitative chemical-potential response within the descriptive start-family uncertainty. The calculation is reproducible with `scripts/analyze_gc_response.py`.

These checks validate the smaller-size grand-canonical workflow and do not alter the canonical L = 64 audit pass above. The detailed commands, energy conventions, and BKT interpretation are in [GRAND_CANONICAL_RUN_GUIDE.md](GRAND_CANONICAL_RUN_GUIDE.md).

The final fixed-mu campaign used beta = 0.24, 0.27, 0.29, and 0.32 at L = 16, 24, and 32, always with three independent starts. L = 16 and 24 retained five blocks per family; L = 32 retained three, for 156 retained HDF5 blocks in total. Its legacy manifests predate the hash column, but the analyzer checked every HDF5 block's stored geometry, model parameters, seed, and initial N. The current runner requires matching INI hashes and refuses legacy or changed manifests before launch. At beta = 0.24 the winding sums were 0.5593, 0.2576, and 0.1425; at beta = 0.32 they were 2.1454, 2.0244, and 2.0131. Every point had less than 2 percent x-y winding mismatch. Raw intersections with 4/pi occurred at beta approximately 0.27815, 0.28516, and 0.28961, or about 40.13, 39.15, and 38.55 nK under the stated Rb-87 grid mapping.

Those crossings are diagnostics only. At L = 24 and beta = 0.29, winding changed from 1.1988 in block 2 to about 1.514 in blocks 4 and 5 while density rose; the all-block mean of 1.3848 and the blocks-3-to-5 mean of about 1.4859 are materially different. This is direct evidence that the number-sector mode remains slow near the crossover. Longer blocks, more independent starts, a finer beta grid, selected nmax checks, and a logarithmic finite-size fit are still required for a publishable grand-canonical Tc.

## Worm-parameter tuning

The exact-N window was tested at L = 16 using two independent replicas and long thermalisation/measurement blocks:

| Window | Kinetic energy, replicas | Winding sum, replicas | Decision |
|---:|---:|---:|---|
| 0.1 | -2179, -2127 | 2.467, 2.107 | trapped; reject |
| 0.5 | -2019, -2025 | 1.365, 1.384 | exact N; use |
| 1.0 | -2022, -2018 | 1.340, 1.332 | mixes, but not strictly exact N |
| 2.0 | about -2021 | 1.385, 1.380 | mixes, but not strictly exact N |

The 0.5 window gives the large mixing improvement while keeping closed configurations at exactly N. The executable now rejects 1.0 and larger.

An occupation-cutoff comparison rejected nmax = 8. Results for 10, 12, and 16 were mutually compatible within residual chain scatter; nmax = 16 was retained conservatively.

The selected update probabilities are:

```text
p_moveworm = 0.4
p_insertkink = 0.3
p_deletekink = 0.2
p_glueworm = 0.1
```

The thesis's 500000 updates per site refers to a different classical-field Metropolis calculation and cannot be translated directly into this continuous-time worm code. Here equilibration was assessed from checkpoint sectors, energy, winding autocorrelation, successive reset measurement blocks, and replica comparisons.

## 64 x 64 convergence assessment at 40 nK

### Earlier strict-window family: seed 96040

The earlier family used one 15-minute, 24-rank fresh block followed by four 15-minute reset blocks. The rank checkpoints retained exactly N = 9948 and sampled both x and y winding sectors.

A later checkpoint/source audit found that the fresh block ended at roughly 30000-37000 sweeps per rank, below the configured 50000. The old reset routine then advanced the bookkeeping counter to the thermalisation threshold. This explains why the first reset blocks still drifted. The current source now rejects a reset at or below the threshold and tells the user to Resume first; the historical observation still shows why merely crossing that technical threshold does not prove equilibration.

| Block | Samples | Wx^2 | Wy^2 | Sum | Winding tau range | Kinetic energy |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1,224,312 | 0.633763 | 0.649312 | 1.283075 | 44.1-44.4 | -32216.6 |
| 2 | 1,244,273 | 0.611959 | 0.617730 | 1.229689 | 45.2-45.4 | -32203.1 |
| 3 | 1,259,354 | 0.587457 | 0.589747 | 1.177204 | 46.4-48.3 | -32202.3 |
| 4 | 1,231,369 | 0.597359 | 0.585063 | 1.182422 | 45.3-46.8 | -32209.5 |

Blocks 1 and 2 drift monotonically and must not be treated as equilibrated production. The combined means of blocks 3 and 4 agree. Their simple late-block mean is:

```text
<Wx^2 + Wy^2> = 1.1798
```

For block 4, ALPS reports a quadrature component error of about 0.0108. Treating the 24 rank means as independent units gives a more conservative standard error of about 0.014. The rank means span 1.057 to 1.301 in winding sum and -32271.7 to -32132.4 in kinetic energy. This is why pooled within-chain errors are not used as the only convergence diagnostic.

### Independent corrected-workflow family: seed 230926

The independent family used a different base seed and a new output directory. It ran `Fresh -> Resume -> Resume` as burn-in and then three reset production blocks, each for 900 seconds on 24 physical cores. Fresh stopped during thermalisation, as expected. After the second Resume, all rank sweep counters were between 119372 and 129292, all 24 checkpoint configurations passed the external reconstruction, and 18 ranks occupied nonzero winding sectors.

| Production block | Samples | Wx^2 | Wy^2 | Sum | Winding tau range | Kinetic energy | Potential energy |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1,206,814 | 0.601002 | 0.600731 | 1.201734 | 44.2-44.4 | -32218.0 | 37059.0 |
| 2 | 1,238,032 | 0.584875 | 0.571814 | 1.156689 | 45.4-45.7 | -32193.2 | 37067.7 |
| 3 | 1,238,948 | 0.575524 | 0.588271 | 1.163795 | 45.4-46.3 | -32212.1 | 37074.6 |

The first-to-second winding change is -0.04461 +/- 0.01862 when each rank is treated as one paired unit, a 2.4-standard-error shift. Block 1 is therefore still equilibration-contaminated and is discarded. The second-to-third change is +0.00729 +/- 0.01868; kinetic, potential, and total energies also agree within two paired rank-level standard errors. The count-weighted late result is:

```text
seed 230926, accepted blocks 2+3:
<Wx^2 + Wy^2> = 1.160243
rank-level standard error from 24 late rank means = 0.01175
```

This rerun shows that two 15-minute plain Resume calls were not enough for this particular L = 64 start, even though all ranks had exceeded the nominal 50000-sweep threshold. The safer operational default is now at least three plain Resume blocks after Fresh before the first reset Production block, followed by explicit stationarity checks. No fixed count replaces those checks.

The late count-weighted family summaries are:

| Base seed | Accepted late blocks | Winding sum | Kinetic energy | Potential energy | Total energy |
|---:|---:|---:|---:|---:|---:|
| 96040 | historical 3+4 | 1.179783 | -32205.8 | 37070.6 | 4864.8 |
| 230926 | corrected 2+3 | 1.160243 | -32202.7 | 37071.1 | 4868.4 |

The winding-family difference is 0.01954, about 1.1 combined conservative standard errors when the earlier final-block rank error and the new late-rank error are used. The late energies agree closely. An equal-family central diagnostic is about 1.170, but two families do not support a precision uncertainty estimate. The standard square-lattice universal-jump reference is 4/pi = 1.27324, so the replicated L = 64 point lies below it. This is evidence that 40 nK is on the high-temperature side of the finite-size critical region for L = 64, not a standalone measurement of T_BKT.

The external parser passed every new checkpoint stage with zero failures. It independently verified occupation bounds, event chronology and continuity, paired hopping events, exact N, stored vertex count, and stored winding. The final state had 20 of 24 ranks in nonzero winding sectors. Integrating every imaginary-time occupation segment gave exactly N/4096 = 2.4287109375, with only 2.7e-6 to 8.0e-6 probability at the n = 16 boundary across the audited stages. The cutoff is not visibly saturated.

Historical HDF5 files are retained in `results/rb87_audit`. The independent rerun, stage outputs, and rank tables are retained in `results/full_independent_audit_seed230926/bkt64_seed230926`.

## Temperature bracket

Strict exact-N L = 32 runs used N = 2487, nmax = 16, window 0.5, eight MPI ranks, and the same interaction and density mapping.

| Temperature | beta*t | Samples | Wx^2 + Wy^2 | Propagated component error | Interpretation |
|---:|---:|---:|---:|---:|---|
| 32 nK | 0.348846447 | 125,417 | 2.32921 | 0.0467 | strong winding response |
| 40 nK | 0.279077158 | 401,343 | 1.327767 | 0.0195 | near 4/pi |
| 49.8 nK | 0.224158359 | 402,301 | 0.281886 | 0.0093 | winding largely collapsed |

At 40 nK, L = 32 lies slightly above 4/pi while both independent late L = 64 families lie below it. Together with the low/high-temperature bracket, this is BKT-like finite-size behaviour. It is not enough to quote T_BKT: the L = 32 points were not subjected to the same multi-seed protocol, and a transition estimate still requires a denser temperature grid, at least three useful sizes, independent seed families at the fit points, and a finite-size BKT scaling fit.

## Measurement limitations

The O(1) observables such as energy, particle number, and winding were sampled hundreds of thousands to more than one million times per reported block. The O(N) density matrix and correlation arrays were intentionally measured only every 10000 updates and have only roughly one hundred samples in a 64 x 64 block. Their reported errors are not converged, so no correlation-length or algebraic-exponent claim is made from this test.

There is also a specific estimator defect at zero separation. In a canonical 3 x 3, N = 13 regression, the exact identity is `Density_Matrix[0] = N / 9 = 1.444444...`; the earlier seed returned 1.42728 +/- 0.00167. A fresh independent seed returned 1.43226 +/- 0.00158, still about 7.7 reported standard errors low. This confirms that the deficit is not a fluctuation of one seed. Do not use `Density_Matrix[0]` for normalization or an on-site occupation check until that branch of `measure_density_matrix()` is corrected and revalidated. The tested off-site density-matrix values, energies, density correlations, and winding do not share this on-site estimator branch.

The optional `MATSUBARA=ON` direct `G(k=0,tau)` histogram is not saved and restored with the other in-progress Green-function histograms. A checkpoint continuation can therefore lose its partial accumulation. The distributed BKT executables use `MATSUBARA=OFF`, so this does not alter any result tabulated above; a Matsubara-enabled build needs a separate fix and test before use.

MPI also has two practical details:

- The ALPS scheduler checks completion periodically, so short jobs can overshoot the requested sweep count.
- MPI completion is aggregated over ranks, while thermalisation and checkpoints are rank-local. A continuation must use the same number of ranks.

## Lower-severity robustness findings

These do not affect the supplied Rb-87 files, but they matter if new parameter files or source variants are introduced:

- The constructor now rejects fewer than ten sweeps, negative thermalisation, nonpositive or nonfinite beta/E_off/C_worm/dtol_scale, zero measurement intervals, negative nmax, invalid probability groups, and an ensemble/executable mismatch. This closes the malformed-input cases found in the audit; it is not a formal proof that every possible bad INI is diagnosed.
- Two internal link-direction invariants in interaction-passing code are enforced only with C `assert`. Assertions are disabled in the Release build, so a runtime exception would be safer future hardening. No such failure was observed: all consistency tests passed and the independent parser reconstructed every retained large-run vertex and winding counter exactly.
- Release executable bytes are not bit-for-bit reproducible because the PE linker writes a timestamp/checksum. Dependency provenance, build options, imports, and the distributed package hashes were nevertheless verified.

## What is established and what remains

Established by this audit:

- reproducible pinned ALPSCore/MinGW build;
- functioning 24-rank canonical MPI run on the desktop;
- exact N in retained closed configurations;
- removal of the narrow-window and ordered-initial-state mixing problems;
- exact reconstruction of N, event pairing, occupation bounds, vertex count, and winding in every newly audited 64 x 64 rank checkpoint;
- late-block stationarity in the new seed 230926 family after discarding its first drifting production block;
- agreement of two independent strict-window L = 64 base-seed families at 40 nK within conservative rank-level uncertainty;
- quantitative grand-canonical chemical-potential response consistent with beta Var(N) to 3.15 percent over the tested finite mu interval;
- start erasure and x-y winding isotropy in the L = 16, 24, and 32 grand-canonical scan; and
- low/near/high-temperature winding behaviour consistent with a BKT crossover hypothesis, together with an explicit detection of unresolved crossover block drift.

Still required for a publishable transition estimate:

- temperatures spaced more finely around 40 nK;
- several lattice sizes, preferably including 48, 64, and larger if affordable;
- at least two independent seeds at every point used in the finite-size fit;
- at least three plain Resume blocks for comparable new L = 64 starts, followed by enough reset blocks to demonstrate rather than assume stationarity;
- explicit finite-size BKT scaling rather than one universal-jump comparison;
- grid-spacing checks to quantify the continuum-discretisation error;
- denser O(N) measurements if correlation functions are part of the analysis.
