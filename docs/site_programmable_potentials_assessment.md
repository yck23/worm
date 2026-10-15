# Feasibility of programmable on-site potentials

Assessment date: 15 October 2026. Source reviewed: commit `c4ebe8a`, active
`src/` implementation. Branch: `assess/site-programmable-potentials`.

Static, independently specified potentials at every simulation site are highly
feasible. The existing nonuniform Hamiltonian and worm updates already provide
the essential physics. This is mainly input, measurement, persistence, and
validation work, rather than a new Monte Carlo algorithm. The current
nonuniform mode is a developer hook, not yet a reliable user-facing feature.

The assessment assumes the primary target is the Bose-Hubbard model used in
the Rb-87 workflows. XXZ has an analogous site-dependent longitudinal field.
No simulation source was changed for this assessment.

## Meaning of programmable

| Requested capability | Assessment |
| --- | --- |
| One site, barrier, trap, disorder, or arbitrary static spatial pattern | Natural extension of the existing nonuniform mode |
| Scripted sequence of independently equilibrated patterns | Straightforward after static support; each pattern needs equilibration and separate statistics |
| Potential depending on imaginary time | Substantial update-algorithm work; existing constant-rate assumptions no longer suffice |
| Physical-time ramps, quenches, pulses, or feedback dynamics | Requires a separate dynamics solver; Monte Carlo sweeps are not physical time |

This program samples a finite-temperature partition function in imaginary time,
as described in the [implementation paper](https://arxiv.org/abs/2204.12262).
The model API has no time argument (`src/model.hpp:26`), and MOVEWORM uses
constant local energy differences for exponential waiting times
(`src/worm.update.cpp:297`, `:439`). Piecewise imaginary-time profiles would
need additional stopping events and consistent proposal/weight calculations.

## Existing support and concrete gaps

1. **The Hamiltonian is already site-aware.** With `UNIFORM=OFF`,
   `BoseHubbard::site_weight_diag` evaluates
   `U_on[s]*n*(n-1)/2 - mu[s]*n` (`src/model.hpp:107`). `Diag_energy` passes
   the actual site into the model (`src/worm.hpp:126`), as do the integrated
   and measurement energy reconstructions (`src/worm.cpp:960`, `:1026`).
   No new worm move appears necessary for a static diagonal potential.

2. **There is no usable potential-map input yet.** `worm::initialize()`
   hardcodes identical `t=1`, `U=4`, `mu=0.7`, and `V_nn=0` on all sites/bonds
   (`src/worm.cpp:334`). Scalar model input definitions are themselves inside
   `UNISYS` guards (`src/model.hpp:81`). Simply rebuilding with
   `-DUNIFORM=OFF` will not use the Rb-87 INI's desired scalar couplings.

3. **Uniform energy caches must not be used for a varying potential.** The
   default build enables `UNISYS` (`src/CMakeLists.txt:30`); its occupation
   lookup table is evaluated at site zero and shared across all sites
   (`src/worm.cpp:242`). Start with the nonuniform path. A later optimization
   could retain uniform bond tables while making diagonal energies site-aware.

4. **Local density measurement is disabled in the nonuniform path.**
   `Density_Distribution` is registered unconditionally (`src/worm.cpp:208`),
   but populated only inside the `UNISYS`/periodic/single-basis guards in
   `measure_corrfun()` (`src/worm.cpp:611`). Extract per-site occupations into
   an independent measurement using `MyModel->stateval`, including open
   boundaries. Winding registration and accumulation are also guarded
   (`src/worm.cpp:213`, `:478`), although the winding is still computed.
   Restore global winding output for applicable periodic geometries if the
   BKT workflow is to continue with an inhomogeneous potential.

5. **Nonuniform checkpoint restore is broken by code inspection.** Model
   vectors are initially empty (`src/model.hpp:66`) and populated only by
   `initialize()`. Restore skips that function (`src/worm.run.cpp:33`), and
   `load()` reconstructs energies without populating the vectors
   (`src/worm.output.cpp:117`, `:229`). Model arrays are not explicitly saved
   by `save()` (`src/worm.output.cpp:72`). Fix this before using the feature
   for production or parameter scans; indexing empty vectors is undefined
   behavior. This failure was identified by inspection, not a crash test.

6. **Correlations need an explicit definition for an inhomogeneous sample.**
   Density-density output currently uses site zero as its reference
   (`src/worm.cpp:621`); the one-body estimator bins relative displacement
   (`src/worm.cpp:682`). These are not a full site-pair matrix. Initially
   provide local densities and energies, then selected reference sites or
   requested pairs if needed. Existing documentation also flags an unrelated
   on-site `Density_Matrix` estimator defect (`RB87_BKT_AUDIT.md:208`).

## Suggested implementation

Expose a physical site potential `onsite_potential`, distinct from the existing
nearest-neighbor interaction parameter `V_nn`. Use the convention

    E_i(n) = U_i*n*(n-1)/2 + V_i*n - mu_global*n
    mu_effective[i] = mu_global - V_i

Positive `V_i` repels particles and negative `V_i` attracts them. Arbitrary real
static diagonal offsets preserve positive diagonal Boltzmann factors, so they
do not introduce a new sign problem into this sign-positive hopping model.

Keep the physical potential and global chemical potential distinguishable in
storage and energy reporting. Folding them into `mu[i]` is enough for the
sampling weights, but the existing `site_weight_diag(..., false)` removes the
whole chemical-potential term. A physical-energy estimator should retain
`V_i*n_i` when excluding only `-mu_global*N`. Currently reported total energy
includes the chemical-potential contribution.

A practical initial interface would retain scalar `t_hop`, `U_on`, `V_nn`, and
`mu` in both modes, with an optional dense potential file or sparse
`site_id,value` overrides on a zero background. Read and validate the profile
once at initialization. Check sizes, coordinate ranges, duplicate entries,
finite values, and units. A small generator script can produce arbitrary
patterns without adding an expression language to the C++ simulation.

Square-lattice indexing is zero-based `site_id = x + Lx*y`. General indexing
is `basis + N_BASIS*(x + Lx*(y + Ly*z))`
(`src/lattice.hpp:97`, `:177`). Publish this convention and save the site
coordinates with the resolved map. Use the same energy units as `t_hop`,
`U_on`, and `mu`. Site control here means control of simulation grid sites;
matching experimental spatial resolution requires a separate discretization
and potential-calibration check.

Move Hamiltonian setup out of fresh-worldline initialization, so fresh and
restored runs both have valid arrays before any energy access. Save the
resolved map and its lattice metadata inside checkpoint and result HDF5 files;
store the source filename/hash only as provenance. Restore must not depend on
the original file still existing or containing the same values. MPI replicas
must use the same resolved Hamiltonian.

For equilibrium scans, start a fresh run at each map initially. Warm starts
are a later convenience: changing a Hamiltonian requires rebuilding cached
energies, discarding old statistics, and rethermalizing. The existing
`--reset-statistics` flag alone does not provide this transition. In fixed-N
runs, site potentials redistribute particles while total N stays fixed; in
grand-canonical runs, total N can change as well.

Potential storage is O(number of sites): one double per site is 32 KiB for
64x64 or 128 KiB for 128x128. Local update scaling is unchanged, though the
nonuniform path has extra lookup/arithmetic cost. Strong barriers can slow
mixing, and attractive wells require occupation-cutoff (`nmax`) convergence
checks. Neither runtime overhead nor convergence has been benchmarked here.

## Verification performed

- Configured an isolated square-lattice build with `UNIFORM=OFF`, using the
  installed ALPSCore dependency, and compiled the serial `qmc_worm` target.
  Build passed; the compiler reported two existing unused-variable warnings.
- Ran the unmodified nonuniform mode on a periodic 4x4 lattice, beta=1,
  nmax=4, 20 thermalization sweeps, and 200 production measurements, using
  its hardcoded couplings. Periodic and final configuration checks passed.
- The run log reported zero density samples. HDF5 inspection found 200
  `Total_Energy` samples and no `Density_Distribution` or
  `Winding_number_squared` result group, confirming the measurement gap.
- Build configuration, executable, logs, INI, and HDF5 fixtures are under
  the ignored `build/site-potential-assessment/` directory. No packaged
  executable was replaced. The Windows run required `TMPDIR`, as in the
  repository's supplied runners.

This is a compilation/fresh-run smoke check of the existing nonuniform path,
not a validation of a nonconstant potential. No arbitrary-map implementation,
restart test, MPI test, or converged physical benchmark was performed.

## Acceptance checks for implementation

1. Zero map reproduces the uniform solver statistically. A constant map is
   equivalent to shifting mu by minus that constant.
2. At fixed N, a constant map leaves density statistics unchanged and shifts
   the reported energy by `V*N`.
3. Compare a small nonuniform finite-temperature system against exact
   diagonalization, checking every site density and energy; also test the
   zero-hopping limit. Verify the sign and location of a single-site defect.
4. Confirm sum of mean site occupations equals mean total N; test deep-well
   nmax convergence and independent-start mixing with large barriers.
5. Test save/restore with the original map file removed or changed; verify
   identical restored Hamiltonian and valid cached energies. Cover canonical,
   grand-canonical, serial, and MPI variants.
6. Reject malformed profiles and verify site-coordinate mapping, including
   open boundaries. Existing test_mpi builds force UNISYS and therefore do
   not cover this extension (`test_mpi/BoseHubbard_GrandCanonical/CMakeLists.txt:31`).

Estimated engineering effort, assuming the present working build environment:
roughly one day for a basic static-map prototype, and 3-5 working days for a
usable, documented implementation with density output, checkpoint/MPI handling,
and small-system validation. Broad site-pair correlations, performance tuning,
and scientific convergence campaigns are additional work. These are planning
estimates, not measured implementation times.
