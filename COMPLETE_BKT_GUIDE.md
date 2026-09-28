# From the 2D Bose-Hubbard path integral to an audited BKT calculation

This document is a self-contained guide to the physics, Monte Carlo algorithm, code, Windows build, Rb-87 mapping, and numerical workflow in this repository.

Guide status: 28 September 2026, repository HEAD fb73840 plus the audited working-tree changes described here. To read the formatted version in VS Code, press Ctrl+Shift+V. Backticks and dollar signs visible in a plain-text editor are Markdown/LaTeX formatting marks, not characters to type unless they occur inside a command block.

The immediate target is a homogeneous, periodic, two-dimensional Bose gas represented by a 64 x 64 Bose-Hubbard lattice. The supplied parameter set approximates the Rb-87 regime. It is a controlled lattice calculation of BKT-like physics.

The shortest operational instructions remain in [BASIC_RUN_GUIDE.md](BASIC_RUN_GUIDE.md). The numerical evidence is tabulated in [RB87_BKT_AUDIT.md](RB87_BKT_AUDIT.md), and the reproducible compiler/dependency details are in [WINDOWS_BUILD.md](WINDOWS_BUILD.md).

## 0. Orientation: what is fixed and what is measured

The independently replicated 64 x 64 reference workflow does **not** ask the simulation to find a density from a chemical potential. It uses the canonical MPI executable, so the user chooses the exact total particle number before the run. A separate mu-controlled path is described in [GRAND_CANONICAL_RUN_GUIDE.md](GRAND_CANONICAL_RUN_GUIDE.md). That path has passed independent-start and response tests at smaller audited sizes, while a 64 x 64 grand-canonical result still requires its own long convergence campaign.

For the 64 x 64 input:

| Question | Answer |
|---|---|
| Which ensemble is used? | Canonical, with exact fixed N |
| Which input fixes N? | `canonical = 9948` |
| What filling does this give? | 9948 / 4096 = 2.4287109375 particles per site |
| Does `mu = -3.47` determine N? | No |
| What should `Number_of_particles` report? | Exactly 9948 in every measured closed configuration |
| What does each MPI rank do? | Runs a complete independent 64 x 64 Markov chain |

The value 9948 was chosen before simulation. A weak homogeneous-gas critical-density estimate gives a target filling 2.428635525 at 40 nK under the selected grid mapping, and rounding this filling times 64 squared gives 9948. The calculation then studies the Bose-Hubbard system at that chosen density; it does not solve an equation of state to obtain N.

A fresh chain is also concrete. The canonical initializer puts two particles on each of 4096 sites and adds one more particle to 1756 randomly selected sites. Worm moves subsequently redistribute the occupations, while number-conserving hopping and the canonical closure condition keep the sum at 9948 whenever the worldlines are closed.

The worm temporarily opens a worldline so that one endpoint can travel through space and imaginary time. During that auxiliary open excursion, the imaginary-time-averaged number may depart fractionally from 9948. The `canonical_window = 0.5` constraint limits this temporary displacement. It does not describe a measured uncertainty in N: each ordinary Monte Carlo step returns to a closed configuration before the standard energy, number, and winding measurements are recorded.

The chemical-potential term still requires care. At fixed N, changing `mu` adds the same constant $-\Delta\mu N$ to every closed-state energy, so it cannot change their normalized relative probabilities or tune the density. Nevertheless:

- the stored `Potential_Energy` and `Total_Energy` include $-\mu N$ and therefore shift numerically if `mu` changes;
- the open-worm proposal formulae contain `mu`, so it can change sampling efficiency even though the correct closed canonical distribution is unchanged.

For the supplied values, $\mu N=-3.47\times9948=-34519.56$. To remove the chemical-potential term from a reported potential or total energy, use

$$
E_{\text{without }-\mu N}=E_{\text{reported}}+\mu N
=E_{\text{reported}}-34519.56.
$$

This also explains the file lifecycle. An INI file defines a new chain only when `Fresh` is used. Its parameters are then stored in rank-specific checkpoints. `Resume` and `Production` continue those checkpoint parameters; editing the INI does not redefine an existing chain. A changed temperature, size, N, interaction, cutoff, worm window, seed, or rank count therefore needs a new INI, a new output directory, and a new `Fresh` run.

## 1. The physical problem

### 1.1 Bose-Hubbard Hamiltonian

The program samples the finite-temperature partition function of

$$
H = -t\sum_{\langle i j\rangle}
    \left(b_i^\dagger b_j+b_j^\dagger b_i\right)
  + \frac{U}{2}\sum_i n_i(n_i-1)
  + V\sum_{\langle i j\rangle}n_i n_j
  - \mu\sum_i n_i .
$$

Here $b_i^\dagger$ and $b_i$ create and destroy a boson on site $i$, $n_i=b_i^\dagger b_i$, $t$ is the nearest-neighbour hopping, $U$ is the on-site repulsion, $V$ is a nearest-neighbour density interaction, and $\mu$ is the chemical potential. The present calculation uses $V=0$, sets $t=1$ as the energy unit, and imposes periodic boundary conditions in $x$ and $y$.

At integer filling and sufficiently large $U/t$, the same model can exhibit a Mott insulator. That is not the regime studied here: $U/t=0.152$ is weak and the filling is about 2.43 particles per grid site. The lattice is being used primarily as a spatial discretisation of a weakly interacting continuum Bose gas.

In the exact-particle-number executable, every closed worldline configuration has fixed total particle number $N$. The term $-\mu N$ is then an additive constant, so changing `mu` cannot change normalized closed-state probabilities or tune the density. It does shift the numerical `Potential_Energy` and `Total_Energy` outputs by $-\Delta\mu N$, and it can affect open-worm acceptance and mixing efficiency. Section 0 gives the practical interpretation.

### 1.2 Why the transition is BKT rather than ordinary symmetry breaking

In an infinite uniform two-dimensional system with short-range interactions, thermal fluctuations prevent true long-range order of a continuous $U(1)$ phase at nonzero temperature. The low-temperature phase nevertheless has superfluid rigidity and quasi-long-range order:

$$
g_1(r)=\langle b^\dagger(\mathbf r)b(0)\rangle
\sim r^{-\eta(T)}.
$$

Its topological excitations are vortices. An isolated vortex costs an energy proportional to $\rho_s\ln(L/\xi)$, while its positional entropy also grows logarithmically with system size. At low temperature, vortices occur mainly as bound vortex-antivortex pairs. At the Berezinskii-Kosterlitz-Thouless transition these pairs unbind, screening the phase stiffness and causing the superfluid density to jump in the thermodynamic limit.

With $k_B=1$, the Nelson-Kosterlitz universal jump is

$$
\frac{\rho_s(T_{\rm BKT}^{-})}{T_{\rm BKT}}=\frac{2}{\pi}.
$$

At the transition, the algebraic exponent is $\eta=1/4$. On finite lattices there is no singular jump. Winding fluctuations change smoothly and acquire unusually slow logarithmic size corrections. That is why one temperature and one lattice size cannot establish $T_{\rm BKT}$, even if the winding happens to equal the thermodynamic universal-jump value.

### 1.3 Winding numbers and superfluid response

Periodic imaginary-time paths can wrap around the spatial torus. Let $W_x$ and $W_y$ be their integer net winding numbers. A boundary phase twist couples directly to these topological sectors, and differentiating the free energy with respect to that twist gives the helicity modulus. For this code's normalization in two dimensions,

$$
\rho_s = \frac{\langle W_x^2+W_y^2\rangle}{2\beta}.
$$

Combining this with the universal jump gives the thermodynamic critical criterion

$$
W_\Sigma \equiv \langle W_x^2+W_y^2\rangle
=\frac{4}{\pi}=1.273239545\ldots .
$$

The result stored as `Winding_number_squared` contains the two components separately. They must be added. Agreement of the $x$ and $y$ components is also an important isotropy diagnostic.

The leading finite-size behavior at criticality is logarithmic:

$$
W_\Sigma(L,T_{\rm BKT}) = \frac{4}{\pi}
\left[1+\frac{1}{2\ln(L/L_0)}+\cdots\right],
$$

where $L_0$ is nonuniversal. This slow correction motivates simulations at several sizes and a Weber-Minnhagen-type fit rather than a simple crossing with $4/\pi$.

## 2. From the partition function to continuous-time worldlines

### 2.1 Interaction expansion

Split the Hamiltonian as

$$
H=H_0-H_1,
$$

where $H_0$ contains the diagonal $U$, $V$, and $\mu$ terms in the occupation-number basis, and $H_1$ contains the positive hopping operators proportional to $t$. In the interaction representation,

$$
Z=\operatorname{Tr}\left[e^{-\beta H_0}
\mathcal T_\tau\exp\left(\int_0^\beta d\tau\,H_1(\tau)\right)\right].
$$

Expanding the exponential gives

$$
Z=\sum_{k=0}^{\infty}\int_{0<\tau_1<\cdots<\tau_k<\beta}
d\tau_1\cdots d\tau_k
\sum_{\{\alpha\}}
\prod_{p=0}^{k}e^{-E_{\alpha_p}(\tau_{p+1}-\tau_p)}
\prod_{p=1}^{k}
\langle\alpha_p|H_1|\alpha_{p-1}\rangle,
$$

with $\tau_0=0$, $\tau_{k+1}=\beta$, and basis states $\lvert\alpha_p\rangle$ specified by all site occupations. Each off-diagonal operator moves one boson across a bond. Its bosonic matrix element is, for a hop $j\to i$,

$$
t\sqrt{(n_i+1)n_j}.
$$

The sampled configuration is therefore a set of piecewise-constant occupation worldlines. A hopping operator is a kink joining events on neighboring sites. Between kinks, the state is diagonal and contributes an exponential of its diagonal energy integrated over imaginary time.

This is a continuous-time expansion: the $\tau_p$ are real-valued Monte Carlo variables, so there is no Trotter time step and no time-discretisation error. There are still statistical errors, autocorrelation, finite-size error, the site-occupation cutoff `nmax`, and—in this continuum mapping—spatial lattice-spacing error.

For the unfrustrated Bose-Hubbard parameters used here all configuration weights are nonnegative. Reversing hopping signs in a frustrated geometry could reintroduce a sign problem; this implementation and audit do not address that regime.

### 2.2 Why ordinary local updates are inefficient

A closed worldline configuration contributes to $Z$, but local changes constrained to remain closed have difficulty altering particle number, permutation cycles, and spatial winding. Near a superfluid transition, the slow modes of interest are precisely large phase fluctuations and topological winding sectors. A sampler that rarely changes them can report small formal error bars around a biased mean.

The worm algorithm enlarges configuration space to include one open worldline. Its two endpoints act as a creation and an annihilation operator. Schematically,

$$
Z_{\rm MC}=Z+C_{\rm worm}Z_G,
$$

where $Z_G$ is the partition sum containing the two operators of the one-particle Green function and `C_worm` controls the relative statistical weight of the open and closed sectors. The worm head can move through space and imaginary time. When it closes again, it may have changed the global topology of the worldlines. This makes winding sectors accessible with local, reversible moves.

The open sector also has a physical interpretation:

$$
G(i-j,\tau)=
\langle\mathcal T_\tau b_i(\tau)b_j^\dagger(0)\rangle.
$$

The relative endpoint separation samples this Green function. In the present BKT workflow, however, winding is the best-converged transition observable because the expensive full correlation arrays were intentionally measured much less frequently.

### 2.3 The elementary worm moves

One Monte Carlo update in this code is a complete worm excursion: enter the open sector, move the head repeatedly, and return to the closed sector. The source names correspond to the following physical operations.

- `INSERTWORM`: place a head-tail pair on a worldline and enter the Green-function sector.

- `MOVEWORM`: move the head forward or backward in continuous imaginary time while changing the occupation on the segment swept by the head. The proposal distribution is chosen so that much of the diagonal Boltzmann factor cancels in the acceptance ratio.

- `INSERTKINK`: move the head to a neighboring site and insert a hopping event, including the appropriate bosonic square-root matrix element.

- `DELETEKINK`: reverse that operation by removing a compatible hopping event.

- `PASSINTERACTION`: pass an existing interaction vertex in the special local geometries where this is required.

- `GLUEWORM`: reunite compatible head and tail endpoints and return to a closed configuration contributing to $Z$.

Detailed balance is enforced with Metropolis-Hastings acceptance ratios containing the physical weight ratio and the forward/reverse proposal probabilities. Consequently, changing the proposal probabilities affects efficiency and autocorrelation but, if the reverse moves remain accessible, not the stationary physical distribution.

In the closed sector, `p_insertworm=1` is normalized separately. In the open sector the current probabilities

```text
p_moveworm   = 0.4
p_insertkink = 0.3
p_deletekink = 0.2
p_glueworm   = 0.1
```

sum to one. A source-level `sweep` is one update call that begins and ends in the closed sector. An accepted insertion leads to a complete variable-length worm excursion before the call returns; a rejected insertion stays closed. It is not one attempted update per lattice site and not a physical unit of time. This distinction is essential when comparing the present run lengths with the Metropolis-update counts in Beregi's classical-field simulations.

### 2.4 Estimators and sampling cadence

The number of hopping vertices $k$ is an expansion-order estimator. With the source convention $H=H_0-H_1$, the kinetic energy estimator is

$$
E_{\rm kin}=-\frac{\langle k\rangle}{\beta}.
$$

The diagonal potential energy is accumulated from the occupation segments. Particle number and its square, density histograms, winding components, density-density correlations, and the one-body density matrix are also available. The final audit validated the energy and winding estimators, but found that the current on-site density-matrix element is biased; `Density_Matrix[0]` must not be used until its special estimator branch is corrected and revalidated.

The code distinguishes inexpensive $O(1)$ measurements from $O(N_{\rm sites})$ arrays:

- `Nmeasure=1` records energy, particle number, and winding after every completed worm excursion that is eligible for measurement.

- `Nmeasure2=10000` records the full spatial correlations only every ten thousand excursions in the audited file.

The latter setting makes the 64 x 64 test affordable, but it leaves only about one hundred correlation samples per 15-minute block. Those arrays are not adequate for a reliable fit of $g_1(r)\sim r^{-\eta}$; the BKT conclusion must therefore rest primarily on winding, block stability, and finite-size behavior.

## 3. What the original Sadoune-Pollet code does

### 3.1 Provenance and scope

The upstream project is [LodePollet/worm](https://github.com/LodePollet/worm), accompanying Nicolas Sadoune and Lode Pollet, [Efficient and scalable Path Integral Monte Carlo Simulations with worm-type updates for Bose-Hubbard and XXZ models](https://arxiv.org/abs/2204.12262). The current checkout is commit `fb73840`. Its `upstream/main` and the configured fork's `origin/main` point to that same commit.

There is no historical release tag that cleanly labels “the exact paper snapshot.” For a reproducible description, this guide uses two source-history layers:

1. The original/published implementation and its maintenance lineage through commit `4966283`.
2. The later committed measurement and diagnostic work from `2a4c97e` through the current merge `fb73840`.

The Windows, canonical-mixing, runner, and documentation changes in the present worktree are a third layer discussed in Section 5. This distinction matters: it avoids attributing all current functionality either to the paper authors or to the local adaptation.

The project supports Bose-Hubbard and XXZ models on selectable lattices, with serial and MPI front ends. Its main design goal is an efficient continuous-time worm sampler with local data structures, consistency tests, ALPSCore accumulators, HDF5 results, and resumable checkpoints.

### 3.2 Source tree mapped to the physics

The important files are:

| File | Numerical role | Physical meaning |
|---|---|---|
| [`src/model.hpp`](src/model.hpp) | Supplies diagonal energies and off-diagonal matrix elements. | Implements $U n(n-1)/2-\mu n$, nearest-neighbour $V$, and bosonic hopping square roots. |
| [`src/lattice.hpp`](src/lattice.hpp) | Builds sites, bonds, directions, neighbors, and periodic boundaries. | Defines the spatial graph on which worldlines hop and wind. |
| [`src/worm.Element.hpp`](src/worm.Element.hpp) | Defines one event in a site's chronological container. | Stores imaginary time, occupation before/after the event, linked site, event type/color, and neighbor associations. |
| [`src/worm.hpp`](src/worm.hpp) | Holds Markov-chain state and declarations. | Represents the complete worldline configuration, worm endpoints, vertex count, particle number, winding, estimators, and RNG state. |
| [`src/worm.cpp`](src/worm.cpp) | Defines parameters, initializes configurations, precomputes matrix elements, measures, checks, and resets. | Connects user inputs to the ensemble and constructs the starting worldlines. |
| [`src/worm.update.cpp`](src/worm.update.cpp) | Implements insertion, motion, kink, deletion, passing, and gluing transitions. | Samples the extended $Z+CZ_G$ ensemble while preserving detailed balance. |
| [`src/worm.output.cpp`](src/worm.output.cpp) | Prints, serializes, and restores state. | Makes a Markov chain restartable without losing its current topological/configurational sector. |
| [`src/worm.run.cpp`](src/worm.run.cpp) | Serial command-line front end. | Parses a fresh input or restored checkpoint, runs, validates, and writes results. |
| [`src/worm.run_mpi.cpp`](src/worm.run_mpi.cpp) | MPI command-line front end. | Runs independent chains on multiple ranks and combines their statistical accumulators. |
| [`src/CMakeLists.txt`](src/CMakeLists.txt) | Selects compile-time variants. | Chooses lattice, worldline container, canonical window build, and optional Green-function measurements. |

The core representation is not a global fine time grid. Each spatial site owns a chronological container of the events that change its occupation. Events carry iterators/associations to relevant events on neighboring sites, so the code can find local temporal neighborhoods without repeatedly scanning the entire operator string. Compile-time data-structure options include list, AVL-based, and list-stack variants; the audited Windows binaries use `LIST`.

Each site also has a dummy event at the imaginary-time boundary. This closes the periodic time direction and provides a well-defined diagonal time slice for measurements and consistency tests. Periodicity in imaginary time represents the trace in $Z=\operatorname{Tr}e^{-\beta H}$; periodicity in space is separately controlled by `pbcx`, `pbcy`, and `pbcz`.

### 3.3 Life of one Markov-chain step

At a high level, the original run loop performs the following cycle:

1. In a closed configuration, choose a site/time and insert a worm pair.
2. Select head moves according to the open-sector proposal probabilities.
3. Move in imaginary time and across bonds, inserting or deleting hopping vertices as accepted.
4. Update occupations, vertex count, particle-number bookkeeping, and winding incrementally.
5. Glue the endpoints when compatible, producing a new closed configuration.
6. If thermalisation is complete and the sampling cadence is due, measure the closed configuration.
7. Periodically perform deep configuration consistency tests and write checkpoints.

The code stores update-attempt and acceptance counters, making it possible to diagnose proposals that almost never succeed. Its built-in `test_conf()` checks occupation bounds, local chronology and continuity, linked-vertex associations, open-worm flags, and a potential energy recomputed from the full configuration. It does not independently recompute particle number, hopping-vertex count, or winding. For the final audit, a separate checkpoint parser reconstructed those three quantities directly from every worldline event in all 24 retained 64 x 64 ranks; every reconstruction matched the stored value and every rank had N = 9948.

### 3.4 Canonical and grand-canonical sampling

The ordinary executable samples the grand-canonical ensemble. Closed configurations can change their total particle number after a worm crosses the imaginary-time boundary before closing. In a modest system this is a powerful way to sample number fluctuations. In early short 64 x 64, high-filling tests, however, starts at occupations 1, 2, and 3 remained in incompatible particle-number sectors and showed no useful winding. Those runs demonstrated inadequate burn-in, not a failure of grand-canonical physics. The later audit therefore uses closer low/central/high starts, complete discarded burn-in invocations, and separately reset production blocks.

The optional `initial_particle_number` parameter now lets a fresh ordinary run begin near an expected mean by randomly distributing exactly that many particles initially. It is only an initial condition: the ordinary executable places no later constraint on N, and `mu` controls the equilibrium distribution. This removes an avoidable distant start but does not solve or conceal number-sector autocorrelation. Independent low and high starts must converge to the same late mean and variance before a grand-canonical result is accepted.

The canonical build constrains the open worm by its imaginary-time-averaged particle-number displacement. If the target integer is $N_0$, adding or removing one particle on a segment of duration $\Delta\tau$ changes the time-averaged number by roughly

$$
|\overline N-N_0|=\frac{\Delta\tau}{\beta}.
$$

The `canonical_window` bounds this quantity during the open excursion. Any value strictly below one prevents a newly closed configuration from ending in the integer sector $N_0\pm1$, while still permitting temporary open-sector fluctuations. A very narrow window can be physically exact yet numerically disastrous because the worm cannot travel far enough in imaginary time to reorganize worldlines and winding.

This is why `canonical_window=0.5`, rather than the old hard-coded 0.1, is central to the current calculation. The window changes the path through extended configuration space, not the physical weights of the retained exact-$N$ closed configurations.

### 3.5 Serial versus MPI behavior

The MPI executable does not divide one 64 x 64 lattice over 24 ranks. It launches 24 statistically independent Markov chains, each holding a complete lattice. Rank $r$ uses the base seed plus $r$, has its own checkpoint, and contributes its ALPSCore accumulators to one combined result.

This has four operational consequences:

- Memory is replicated per rank rather than spatially distributed.

- Increasing ranks increases the number of independent samples and improves robustness against one chain's slow sector, but it does not accelerate one worm trajectory.

- A resumed calculation must use the same rank count because all rank-specific checkpoints must exist.

- `thermalization` is a rank-local sweep threshold, whereas the MPI scheduler's reported/completion work is aggregated. Wall-clock stopping is checked periodically, so short runs can overshoot nominal sweep targets.

The desktop has 24 physical cores and 32 logical processors. The build uses 32 parallel compilation jobs; production simulation uses 24 MPI ranks to avoid relying on hyperthreads for 24 memory- and branch-heavy independent chains.

### 3.6 Measurements and checkpoint state

The original observable set includes kinetic, potential, and total energy; particle number and its square; local-density distributions; the one-body density matrix; density-density correlations; and squared winding components. ALPSCore supplies log-binning accumulators, error estimates, and integrated-autocorrelation diagnostics.

A checkpoint stores much more than an output mean. It includes the random-number-generator state, sweep counters, update counters, the complete closed operator/worldline configuration, particle number, winding, vertex count, and ALPS measurement state. Checkpoints are written between completed worm excursions; restore points the head and tail iterators at a dummy event rather than restoring a live open worm. Restoring still continues the same closed-state Markov chain, rather than merely launching a new run with similar parameters.

The physical parameters in a restored checkpoint should be treated as immutable. Editing the INI file after a fresh run does not redefine a chain already stored in checkpoint HDF5. For a different $L$, $T$, $U$, $N$, cutoff, window, or seed, start a new output directory and fresh checkpoints.

## 4. Later committed enhancements already present at HEAD

The current upstream/fork HEAD contains post-baseline work that should not be confused with the uncommitted Windows/Rb-87 adaptation.

### 4.1 Matsubara and imaginary-time Green functions

Commit `2a4c97e` introduced $G(\mathbf k,i\omega_n)$, with subsequent commits correcting the fractional imaginary-time estimator and FFT handling, adding direct $G(\mathbf k=0,\tau)$ bins, and placing the feature behind the `MATSUBARA` compile option. Per-step Green-function measurements were also gated so they do not contaminate accumulators during thermalisation.

The packaged binaries used for the present winding audit were built with `MATSUBARA=OFF`. Therefore these observables are part of the source history but not enabled in this executable. Rebuilding with that option changes measurement cost and should be validated separately before production use.

The final source audit found a checkpoint gap in this optional path: `hist_grtau_re` and `hist_grtau_im` are serialized, but the in-progress direct `hist_g0tau` histogram is not. A resumed `MATSUBARA=ON` run can therefore lose that histogram's partial interval. This has no effect on the distributed BKT binaries, but the optional direct tau estimator should be fixed and tested before production use.

### 4.2 Periodic-boundary consistency correction

The history briefly rejected periodic lattices of linear size three because a chronology check could mistake legitimate neighbor associations for corruption. Commit `c3d41ee` fixed the false positive, restored the original test geometry, and the references were checked against exact diagonalization. This was a diagnostic fix, not a change to the physical ensemble.

### 4.3 Resetting statistics after restore

The current command-line front ends recognize `--reset-statistics` manually after restoring a checkpoint. This superseded an earlier attempt to put reset behavior into restored ALPS parameters, which could not work because those parameters are no longer mutable in the required way.

Resetting statistics preserves the current worldline and RNG state but clears measurement accumulators for a new production block. It is exactly what is needed for block-to-block equilibration checks—provided the chain has already burned in. The current routine refuses to reset if the stored sweep count has not passed the configured thermalisation threshold. After a valid reset it returns the counter to that threshold so the next calls are a clean production block. Section 8 explains why passing the threshold is necessary but not sufficient.

### 4.4 Other committed maintenance

The same commit range includes tighter floating-point diagnostics (`dtol`/`dtol_scale` work), dead-code cleanup, and regenerated canonical exact-diagonalization references containing the particle-number field. These changes improve validation and maintainability without changing the target Bose-Hubbard Hamiltonian.

## 5. Changes made for this Windows/Rb-87 workflow

These are the current worktree adaptations created for this calculation. They are not part of the published algorithm, and at the time of this guide they are uncommitted local changes. They preserve unrelated user work and can be inspected with `git diff` and `git status --short`.

### 5.1 Reproducible, verified ALPSCore build

[`scripts/build_windows.ps1`](scripts/build_windows.ps1) turns an environment-specific manual build into one command. It:

1. obtains ALPSCore and checks out the fixed commit `3606edfbabc44e0fa7d05efa2a50c6a8a481340f`;
2. applies the recorded [`patches/alpscore-mingw.patch`](patches/alpscore-mingw.patch);
3. creates fresh Release CMake caches, avoiding stale network-drive mappings;
4. builds and locally installs ALPSCore;
5. builds ordinary and `CWINDOW` square-lattice variants, both serial and MPI;
6. recursively packages the required MinGW and ALPSCore DLLs; and
7. writes `dist/windows-square/SHA256SUMS.txt`.

The compatibility patch uses `unsigned long` for selected ALPS parameter types under Windows. The matching changes in [`src/lattice.hpp`](src/lattice.hpp) and [`src/worm.cpp`](src/worm.cpp) cover lattice sizes, seeds, run/measurement intervals, and optional Matsubara sizes. This resolves the MinGW/Windows serialization mismatch without changing their numerical meaning.

The final fresh-cache audit used GCC 14.2.0, HDF5 1.14.3, Microsoft MPI, and the pinned ALPSCore revision. All 32 packaged payload entries matched their SHA-256 manifest. The package contains:

```text
qmc_worm.exe
qmc_worm_mpi.exe
qmc_worm_canonical.exe
qmc_worm_canonical_mpi.exe
```

Only Windows system/runtime libraries and `msmpi.dll` remain host dependencies. The build is “verified” in the practical sense of source pinning, clean configuration, hash checking, parser/error-path tests, and end-to-end serial/MPI WORM jobs. ALPSCore's entire upstream unit-test suite was not run because this MinGW configuration includes POSIX-only test targets. That limitation is documented rather than hidden.

### 5.2 Better fresh-state initialization

The noncanonical build accepts `initial_occupancy` and validates it against the model's allowed local occupation range. It now also accepts `initial_particle_number`: a non-negative value constructs a randomized near-uniform fresh state with that total, while -1 retains the uniform `initial_occupancy` behavior. This parameter is an unconstrained starting hint, not a thermodynamic control.

More importantly, a fresh exact-$N$ state no longer puts every remainder particle onto consecutive row-major sites. For $N=9948$ and $64^2=4096$ sites, integer division gives a base occupation of 2 and 1756 excess particles. The old initialization placed all 1756 in one contiguous ordering pattern, making a macroscopic artificial density stripe that the Markov chain had to erase.

The revised initializer fills every site with $\lfloor N/N_s\rfloor$, shuffles the site indices reproducibly using that rank's seeded generator, and adds one particle to the first $N\bmod N_s$ shuffled sites. It also rejects target particle numbers incompatible with `nmin`, `nmax`, and lattice size. This does not make the starting state an equilibrium draw, but it removes a large deterministic inhomogeneity while retaining exact $N$.

### 5.3 Configurable and validated canonical window

The original compile-time constant `can_window=0.1` became the runtime parameter `canonical_window`. The constructor now requires it to be finite and satisfy

```text
0 < canonical_window < 1
```

and the selected value is printed in the run header. The strict upper bound guarantees that an integer closed sector cannot differ from the requested exact $N$. The tuned value 0.5 permits much longer open-worm temporal excursions than 0.1 and dramatically reduces trapping without relaxing exactness in closed configurations.

### 5.4 Correct checkpoint restoration and safe no-data output

After checkpoint load, the total potential energy was already recomputed but the measurement estimator `Epot_measure` was not. It is now recomputed as well, preventing a restored chain from carrying a stale potential-energy measurement baseline.

Both serial and MPI front ends now handle a run that reaches its wall-time limit before collecting any post-thermalisation samples. They test the actual Total_Energy measurement count, print that no production measurements were collected, and save parameters without a misleading zero-count NaN results group.

The end-of-run configuration test now reconstructs the boundary-slice particle number and snapshot potential energy independently of their incrementally maintained counters. It compares every boundary occupation with the saved per-site state, checks the grand-canonical number counter, and in a closed canonical configuration rechecks exact N. The executable therefore fails rather than writing a result if these core measurement invariants disagree.

Ensemble choice is also enforced inside the binaries: the ordinary build accepts only canonical = -1, and the CWINDOW build requires non-negative exact N. The runner still performs an earlier user-friendly check. Additional guards reject non-finite or invalid beta, E_off, C_worm, measurement intervals, update probabilities, dtol_scale, and local occupation cutoffs before sampling starts.

### 5.5 Dedicated runners and parameter sweeps

[`scripts/run_rb87_bkt64.ps1`](scripts/run_rb87_bkt64.ps1) defaults to the canonical MPI executable. `-Ensemble GrandCanonical` instead selects the ordinary MPI executable and its separate parameter template/output directory. The runner verifies that canonical mode receives a non-negative `canonical` value and grand-canonical mode receives `canonical = -1`, preventing an accidental ensemble mismatch. It resolves absolute paths, supplies the packaged DLL directory, creates a local temporary directory, restores the caller's PowerShell location, and refuses to overwrite a fresh run. Before continuation it checks every rank-specific checkpoint. Before replacing a result it makes a timestamped archival copy.

Its modes are:

- `Fresh`: start from the edited INI file and create a new chain/checkpoints.

- `Resume`: continue the same chain and existing accumulators.

- `Production`: continue the same chain but pass `--reset-statistics`, producing a statistically separate measurement block.

The dedicated runner deliberately does not override the physical or Monte Carlo values in the INI file. Its command-line controls are the ensemble, mode, process count, parameter-file path, and output-directory path. This realizes the requested “always edit the parameter file” workflow.

[`scripts/run_temperature_sweep.ps1`](scripts/run_temperature_sweep.ps1) is a separate convenience tool. It copies a base INI for each $(L,T)$ point, sets `beta=1/(T/t)`, computes exact particle number from a supplied density, and records the generated jobs in a CSV manifest. Because the generic sweep runner creates one fresh block per point, it is best used for reconnaissance or for generating parameter files with `-DryRun`; it is not by itself a convergence protocol.

The dedicated grand-canonical audit runner creates low, central, and high initial-number families for every size-temperature point, limits the whole batch to the desktop's 24 physical cores, and separates discarded Fresh/Resume histories from numbered Production blocks. New manifests record a SHA-256 hash for every generated INI. Before launching any child, the runner checks all job paths, parameter hashes, process counts, rank checkpoints, and production-index conflicts, so a bad batch fails before only part of the grid starts.

The BKT analyzer validates the Hamiltonian, ensemble, geometry, seed, starting number, and scan coordinates stored inside every HDF5 block. It pools first and second particle-number moments by their measurement counts within a family, gives independent families equal weight, reports between-start scatter and block chronology, checks x-y winding isotropy, and labels a 4/pi crossing as a raw finite-size diagnostic rather than a transition estimate.

### 5.6 Validation actually performed

After packaging, a fresh two-rank canonical smoke run preserved the requested particle number, generated winding measurements, restored all checkpoints, reset statistics, archived the previous output, and returned PowerShell to its caller directory. A deliberately invalid `canonical_window=1.0` failed with the intended validation error.

After the μ-controlled extension, both variants were rebuilt again and all 32 package hashes passed. A two-rank L = 8 grand-canonical smoke test started from `initial_particle_number=155`, changed particle sectors, restored checkpoints, and completed a reset-statistics block. Its reset block had mean N = 241.686 and variance 543.221. The strong shift from the earlier block mean near 192 proves both that N was free and that this short smoke run was not equilibrated. It validates mechanics, not large-system grand-canonical physics.

Scientific validation then compared grand-canonical starts, canonical windows, occupation cutoffs, temperatures, independent rank chains, and successive production blocks. Those outcomes—and the failed choices—are summarized in Section 8 rather than being inferred merely from a successful executable exit code.

## 6. Mapping the Rb-87 regime to the lattice input

### 6.1 Source parameters

The physical inputs come from two primary sources:

- Abel Beregi, [Probing universality of 2D quantum systems with bilayer Bose gases](https://ora.ox.ac.uk/objects/uuid:b2f4f0a1-8576-4528-bbd3-557d273cfbdd), Oxford DPhil thesis (2024). Appendix D.8 uses Rb-87, dimensionless coupling $\widetilde g=0.076$ for 1 kHz axial confinement, $T=40$ nK, and a 0.5 micrometre numerical grid. The thesis also discusses an approximately 32 micrometre uniform region.

- Daniel Steck, [Rubidium 87 D Line Data](https://steck.us/alkalidata/rubidium87numbers.pdf), for $m=1.44316089500\times10^{-25}$ kg.

The thesis's Appendix D calculation is not this worm code: it uses a classical-field Metropolis method. Its quoted 500,000 updates per site for equilibration and samples separated by 1,000 updates per site cannot be copied numerically into “worm sweeps,” because one sweep here is a variable-length open-worldline excursion. The physical scales can be mapped; the algorithmic update count must be calibrated independently by autocorrelation and block stability.

### 6.2 Spatial discretisation

Start with the homogeneous continuum Hamiltonian

$$
H=\int d^2r\left[
\frac{\hbar^2}{2m}|\nabla\psi|^2
+\frac{g_{2D}}{2}\psi^\dagger\psi^\dagger\psi\psi
-\mu_{\rm cont}\psi^\dagger\psi\right].
$$

On a square grid of spacing $a$, use $\psi(\mathbf r_i)\simeq b_i/a$. The finite-difference Laplacian gives

$$
t=\frac{\hbar^2}{2ma^2},\qquad
U=\frac{g_{2D}}{a^2}.
$$

With $\widetilde g=mg_{2D}/\hbar^2$,

$$
\frac Ut=2\widetilde g=0.152.
$$

For $a=0.5$ micrometres and Rb-87,

| Quantity | Value |
|---|---:|
| $64a$ | 32 micrometres |
| $t/h$ | 232.6009778 Hz |
| $t/k_B$ | 11.1630863 nK |
| $U/t$ | 0.152 |
| $T/t$ at 40 nK | 3.58323844 |
| $\beta t=t/(k_BT)$ at 40 nK | 0.279077158 |

The finite-difference kinetic energy also produces a diagonal $+4t n_i$ term on a square lattice. Therefore the lattice convention has $\mu_{\rm BH}=\mu_{\rm cont}-4t$. The supplied `mu=-3.47` is not used to infer the density here and should not be overinterpreted. At exact fixed $N$, $-\mu N$ is a common closed-state energy offset, although it remains visible in reported energies and in open-worm proposal rates.

For the separate grand-canonical audit, the weak-gas critical-chemical-potential
relation with $\widetilde g=0.076$, $\xi_\mu=13.2$, and $T/t=3.58323844$
gives $\mu_{\rm cont}/t\simeq0.4471$, or
$\mu_{\rm BH}/t\simeq-3.5529$. This independently motivates the tested
mu = -3.53 grid. It does not remove the need to measure the lattice equation of
state or to retune mu when following a fixed-density path.

The thermal wavelength at 40 nK is only about 1.87 grid spacings. Thus $a=0.5$ micrometres is a fairly coarse ultraviolet cutoff. The calculation can show the intended finite-temperature topology, but quantitative continuum-experiment comparison requires smaller $a$ at fixed physical size/density and an extrapolation toward $a\to0$.

### 6.3 Density and exact particle number

For a weak homogeneous 2D gas, use the critical phase-space-density estimate

$$
\mathcal D_c\simeq\ln(380/\widetilde g)=8.51719.
$$

Since

$$
\frac{\lambda_T^2}{a^2}=\frac{4\pi}{T/t},
$$

the corresponding site filling at 40 nK is

$$
n_{\rm site}=\frac{\mathcal D_c}{\lambda_T^2/a^2}
=\frac{\mathcal D_c(T/t)}{4\pi}
=2.428635525.
$$

For each finite size, the workflow rounds $n_{\rm site}L^2$ to the nearest integer:

| $L$ | Exact $N$ |
|---:|---:|
| 16 | 622 |
| 24 | 1399 |
| 32 | 2487 |
| 48 | 5596 |
| 64 | 9948 |
| 80 | 15543 |

This choice keeps the density approximately fixed across sizes. It is an informed critical-region target, not a self-consistent determination of the experimental equation of state. A later precision study should test nearby densities as well as grid spacing.

### 6.4 Temperature conversion

The input is dimensionless inverse temperature in units with $t=1$:

$$
\beta t=\frac{t/k_B}{T_{\rm nK}}
=\frac{11.1630863}{T_{\rm nK}},
\qquad
\frac Tt=\frac{T_{\rm nK}}{11.1630863}.
$$

Useful points are:

| Physical $T$ | $T/t$ passed to the sweep tool | `beta` in an INI file |
|---:|---:|---:|
| 32 nK | 2.866590753 | 0.348846447 |
| 34 nK | 3.045752675 | 0.328326068 |
| 36 nK | 3.224914597 | 0.310085731 |
| 37 nK | 3.314495558 | 0.301705036 |
| 38 nK | 3.404076519 | 0.293765429 |
| 39 nK | 3.493657480 | 0.286232982 |
| 40 nK | 3.583238441 | 0.279077158 |
| 41 nK | 3.672819402 | 0.272270398 |
| 42 nK | 3.762400363 | 0.265787769 |
| 44 nK | 3.941562285 | 0.253706507 |
| 46 nK | 4.120724207 | 0.242675789 |
| 49.8 nK | 4.461131859 | 0.224158360 |

Do not confuse the two interfaces: the INI file takes `beta`, while `run_temperature_sweep.ps1 -Temperatures` takes $T/t$ and computes `beta=1/(T/t)`.

## 7. The audited parameter file, line by line

The starting point is [`parameter_files/Rb87_BKT_64.ini`](parameter_files/Rb87_BKT_64.ini). Comments beginning with `#` are ignored. Parameter values in this file are not overridden by the dedicated Rb runner; only `-Processes`, run mode, input path, and output directory come from its command line.

| Parameter | Audited value | Meaning and reason |
|---|---:|---|
| `model` | `BoseHubbard` | Selects the bosonic model and its matrix elements. |
| `runtimelimit` | 900 | Fifteen wall-clock minutes per invocation. This creates manageable convergence blocks and archives. |
| `sweeps` | 4,800,000 | High ceiling; the wall-time limit normally stops the block first. It prevents an unintended short sweep cap. |
| `thermalization` | 50,000 | Rank-local sweep threshold before measurements. On $L=64$, more than one 900 s invocation may be needed to reach it. |
| `Lx`, `Ly`, `Lz` | 64, 64, 1 | A two-dimensional 64 x 64 grid. |
| `pbcx`, `pbcy`, `pbcz` | 1, 1, 0 | A spatial torus in the two active directions, necessary for winding. |
| `t_hop` | 1.0 | Energy unit. All energies and temperatures are dimensionless relative to $t$. |
| `U_on` | 0.152 | $2\widetilde g$ for the selected continuum mapping. |
| `V_nn` | 0.0 | No nearest-neighbour density interaction in the continuum contact model. |
| `mu` | -3.47 | Does not set density at exact $N$. It supplies the reported $-\mu N$ energy offset and enters open-worm proposals, so it may affect sampling efficiency. |
| `beta` | 0.279077157505874 | 40 nK under the selected Rb/grid conversion. |
| `nmax` | 16 | Local Fock-space cutoff. 8 was insufficient; 10, 12, and 16 were compatible within residual chain scatter, so 16 is conservative. |
| `E_off` | 1.0 | Technical offset used in temporal move proposals; it changes sampling efficiency, not the target Hamiltonian. |
| `canonical` | 9948 | Exact total number in every retained closed configuration. |
| `canonical_window` | 0.5 | Tuned open-worm number window: mobile enough to mix, still strictly below 1 and therefore exact when closed. |
| `initial_occupancy` | 2 | Fresh noncanonical start only. The canonical initializer instead constructs exactly 9948 particles, so this value is ignored by this executable. |
| `seed` | 96040 | Base RNG seed; MPI rank $r$ uses this plus $r$. New replica groups need a different base seed. |
| `Ntest` | 1,000,000 | Interval for expensive deep consistency checks. A final validation is also performed. |
| `Nsave` | 100,000,000 | Update-based checkpoint interval; ALPS scheduling/checkpoint behavior also responds to stopping. The large value avoids excessive network-file writes. |
| `Nmeasure` | 1 | Measure cheap energy/number/winding observables after each eligible completed worm excursion. |
| `Nmeasure2` | 10,000 | Measure expensive full spatial arrays sparsely. Adequate for smoke diagnostics, not precision correlation exponents. |
| `C_worm` | 2.0 | Relative statistical weight of open configurations. It tunes time in the Green-function sector. |
| `p_insertworm` | 1.0 | Closed-sector insertion choice; normalized separately. |
| `p_moveworm` | 0.4 | Open-sector temporal-head-move probability. |
| `p_insertkink` | 0.3 | Open-sector spatial/kink insertion probability. |
| `p_deletekink` | 0.2 | Open-sector reverse kink deletion probability. |
| `p_glueworm` | 0.1 | Open-sector closure proposal probability. Open-sector probabilities sum to one. |
| `outputfile` | `Rb87_BKT_64.out.h5` | Combined result HDF5 filename inside the chosen output directory. |
| `checkpoint` | `Rb87_BKT_64.clone.h5` | Rank-0 checkpoint basename; other ranks append `.1`, `.2`, and so on. |

The proposal and technical parameters are not universal physical constants. Their correctness comes from detailed balance; their quality comes from whether they let the chain decorrelate and change winding sectors. If they are retuned, acceptance counters, autocorrelation, and independent replicas must be audited again.

## 8. What converged, what did not, and why

### 8.1 Grand-canonical number-sector equilibration

Early high-filling tests initialized with one, two, or three particles per site remained near incompatible total-particle sectors and did not develop useful winding within their short run times. This was not evidence that the equilibrium state has zero superfluid response. It showed that those chains had not made enough global imaginary-time topology changes to forget their starts.

The follow-up audit changed the protocol rather than the target distribution: it initialized closer to the expected density, compared three deliberately separated starts, discarded whole Fresh and Resume invocations, and retained multiple reset blocks. At L = 16, the late families agreed and a finite-difference change of mean N with mu matched beta times the number variance to 3.15 percent over a finite delta mu of 0.05. A fixed-mu L = 16, 24, and 32 temperature sweep then tested whether normal and superfluid winding regimes could both be reached and exposed the remaining slow number-sector sampling near the crossover. Full tables and reproducible commands are in [GRAND_CANONICAL_RUN_GUIDE.md](GRAND_CANONICAL_RUN_GUIDE.md).

This makes grand-canonical sampling a validated smaller-size workflow, not an automatic certificate for L = 64. Exact-N sampling still removes the slow number-sector degree of freedom and remains the independently replicated 64 x 64 reference. In either ensemble, winding itself can remain slow.

### 8.2 Canonical-window scan

At $L=16$, two independent replicas were run for each candidate window with long thermalisation and measurement periods:

| Window | Kinetic energy in two replicas | $W_x^2+W_y^2$ in two replicas | Interpretation |
|---:|---:|---:|---|
| 0.1 | -2179, -2127 | 2.467, 2.107 | Strong replica disagreement and unphysical trapping; reject. |
| 0.5 | -2019, -2025 | 1.365, 1.384 | Replicas agree, exact closed $N$, autocorrelation roughly 33-39; use. |
| 1.0 | -2022, -2018 | 1.340, 1.332 | Good mobility, but cannot guarantee exact closed $N$; reject for this ensemble. |
| 2.0 | about -2021 | 1.385, 1.380 | Good mobility, but not an exact-$N$ constraint; reject. |

The 0.1 result is a useful warning. A tighter mathematical constraint did not make the result “more canonical” in a useful numerical sense; both 0.1 and 0.5 keep closed states at exact $N$, but 0.1 prevented sufficiently long imaginary-time excursions. It trapped the Markov chain in atypical configurations. The 0.5 value retained exactness while lowering winding autocorrelation by over an order of magnitude relative to the trapped behavior.

### 8.3 Local occupation cutoff

`nmax=8` shifted potential-energy and winding estimates and was rejected. Results at 10, 12, and 16 were broadly compatible within the remaining between-chain scatter. The production value 16 gives headroom in the local Fock basis. A precision calculation should still record the maximum-occupation histogram or repeat selected points at a larger cutoff; “16 worked here” is not a theorem for every $U$, density, and temperature.

### 8.4 The thermalisation/reset caveat

The first audited source routine behind `--reset-statistics` unconditionally set its sweep counter to the configured thermalisation count. In the first $L=64$, 900-second fresh run, rank checkpoints had only roughly 30,000-37,000 sweeps, below `thermalization=50000`. Calling `Production` immediately therefore skipped the remaining nominal burn-in in the bookkeeping, even though the physical chain state had not magically equilibrated.

The current source closes that loophole: Production now fails with a clear instruction to Resume when a restored chain has not passed the threshold. This is a safety guard, not an equilibration oracle. A chain can pass the local sweep threshold while its particle-number or winding sectors are still drifting.

This explains the early production drift. A fresh independent rerun then showed that even two plain Resume blocks could leave the first reset block elevated. For this particular L = 64 point and 900-second block length, the recommended default is therefore:

```text
Fresh -> Resume -> Resume -> Resume -> Production -> Production -> Production -> ...
```

The three plain resumes continue both the physical chain and its honest sweep counter while their data are treated as burn-in. Only afterward should `Production` reset the measurements. This is safer than relying on the nominal threshold alone because winding-sector equilibration can be slower than local energy equilibration.

Do not infer that exactly three resumes are universal. It is an empirical safer starting point based on the observed evolution of this machine, size, and parameter set. The actual criterion is that multiple subsequent reset production blocks and independent seed groups agree. A visibly drifting first production block is burn-in evidence, not data to retain merely because it was labelled Production.

### 8.5 Independently replicated 64 x 64 diagnostic at 40 nK

The historical seed-96040 audit used one 15-minute fresh invocation followed by four 15-minute reset blocks on 24 ranks. It preserved exact $N=9948$, sampled nonzero $x$ and $y$ winding sectors, and produced over one million cheap-observable samples per block:

| Block | Samples | $\langle W_x^2\rangle$ | $\langle W_y^2\rangle$ | $W_\Sigma$ | Winding $\tau$ range | Kinetic energy |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1,224,312 | 0.633763 | 0.649312 | 1.283075 | 44.1-44.4 | -32216.6 |
| 2 | 1,244,273 | 0.611959 | 0.617730 | 1.229689 | 45.2-45.4 | -32203.1 |
| 3 | 1,259,354 | 0.587457 | 0.589747 | 1.177204 | 46.4-48.3 | -32202.3 |
| 4 | 1,231,369 | 0.597359 | 0.585063 | 1.182422 | 45.3-46.8 | -32209.5 |

Blocks 1 and 2 drift monotonically and must be discarded as equilibration-contaminated. The combined means of blocks 3 and 4 agree in energy, winding sum, and approximate $x/y$ symmetry. Their simple late-block mean is

$$
W_\Sigma(L=64,40\,{\rm nK})\simeq1.1798.
$$

For block 4, the quadrature of the two ALPS component errors is about 0.0108. Treating the 24 rank means as the independent units gives a more conservative standard error of about 0.014. The rank winding sums span 1.057 to 1.301, so this family alone could not prove that a rank-wide offset had decorrelated.

A new seed-230926 family was then run from a fresh randomized exact-N state in a separate directory. Fresh and two plain Resume calls were treated as burn-in, followed by three reset blocks:

| Block | Samples | $\langle W_x^2\rangle$ | $\langle W_y^2\rangle$ | $W_\Sigma$ | Winding $\tau$ range | Kinetic energy | Potential energy |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1,206,814 | 0.601002 | 0.600731 | 1.201734 | 44.2-44.4 | -32218.0 | 37059.0 |
| 2 | 1,238,032 | 0.584875 | 0.571814 | 1.156689 | 45.4-45.7 | -32193.2 | 37067.7 |
| 3 | 1,238,948 | 0.575524 | 0.588271 | 1.163795 | 45.4-46.3 | -32212.1 | 37074.6 |

The first-to-second change is 2.4 paired rank-level standard errors, so block 1 is discarded as residual burn-in. Blocks 2 and 3 agree in winding and all energy components within two paired rank-level standard errors. Their count-weighted winding sum is $1.160243$, with a conservative standard error of about 0.01175 from the 24 combined late rank means.

The old seed-96040 late mean is 1.179783. The independent-family difference is 0.01954, about 1.1 combined conservative standard errors, and their late kinetic, potential, and total energies agree closely. The two-family equal-weight diagnostic is therefore approximately

$$
W_\Sigma(L=64,40\,{\rm nK})\simeq1.170,
$$

without a precision error claim from only two base seeds. Both families lie below $4/\pi$. The point is now independently replicated for energy and winding, but one replicated point is not a BKT transition estimate.

An external parser reconstructed all 24 new rank checkpoints after every long stage. It found zero occupation-range, chronology, continuity, event-pairing, particle-number, vertex-count, or winding mismatches. The final stage had 20/24 ranks in nonzero winding sectors. Direct imaginary-time integration gave exact mean filling and a probability of only a few parts in $10^6$ at the `nmax=16` boundary.

The methodological point remains important: neither a small pooled error nor a Production label proves equilibration. The first new production block still drifted after two resumes. Use at least three plain resumes as the default for a comparable start, then decide from chronology, rank scatter, and independent seeds which reset blocks are admissible.

### 8.6 Temperature bracket at L = 32

Strict exact-$N$, `canonical_window=0.5` runs with $N=2487$ gave:

| Temperature | `beta` | Samples | $W_\Sigma$ | Quadrature component error | Physical reading |
|---:|---:|---:|---:|---:|---|
| 32 nK | 0.348846447 | 125,417 | 2.32921 | 0.0467 | Strong superfluid winding response. |
| 40 nK | 0.279077158 | 401,343 | 1.327767 | 0.0195 | Near the universal-jump scale. |
| 49.8 nK | 0.224158359 | 402,301 | 0.281886 | 0.0093 | Winding largely collapsed. |

This is the expected low/near/high-temperature pattern. At 40 nK, $L=32$ is slightly above $4/\pi$ while both replicated late $L=64$ families are below. That pattern is consistent with BKT-like behavior, but it is not a converged finite-size crossing. The L = 32 points lack the same independent-seed treatment, and three temperatures, two highlighted sizes, and one coarse spatial lattice do not determine logarithmic finite-size scaling reliably.

### 8.7 Grand-canonical finite-size audit

The independent fixed-mu audit used $\mu/t=-3.53$, $U/t=0.152$, $n_{\max}=16$, and three separately seeded starting-number families at every point. Fresh and Resume stages were discarded. Five reset blocks per family were retained at $L=16$ and 24, and three per family at $L=32$, giving 156 retained HDF5 blocks.

| $L$ | $W_\Sigma$ at $\beta t=0.24$ | $W_\Sigma$ at $\beta t=0.27$ | $W_\Sigma$ at $\beta t=0.29$ | $W_\Sigma$ at $\beta t=0.32$ | raw $4/\pi$ crossing |
|---:|---:|---:|---:|---:|---:|
| 16 | 0.5593 | 1.1026 | 1.5216 | 2.1454 | 0.27815 |
| 24 | 0.2576 | 0.9237 | 1.3848 | 2.0244 | 0.28516 |
| 32 | 0.1425 | 0.7299 | 1.2841 | 2.0131 | 0.28961 |

The high-temperature values are below the universal-jump scale, the low-temperature values are above it, and the raw crossings move systematically with size. Both winding components agree to within 2 percent at every point. The deliberately imposed low-to-high initial particle-number spans shrink from 307, 691, and 1229 particles to late family-mean spans no larger than 19.2, 42.6, and 67.7. Together these are convincing checks of start erasure, rotational symmetry, and finite-temperature BKT-like response.

The same audit also shows what has not converged. At $L=24$, $\beta t=0.29$, the family-averaged winding rose from 1.1988 in block 2 to about 1.514 in blocks 4 and 5 while particle number rose. The five-block mean is 1.3848, but the blocks-3-to-5 mean is about 1.4859. Several other points have visible number motion, with a minimum estimated effective number sample count near 141. Thus the raw crossings are useful bracketing diagnostics, not fitted transition temperatures. A final result needs longer blocks, more starts, finer temperatures, occupation-cutoff checks, and a logarithmic finite-size fit. The full density table and exact commands are in [GRAND_CANONICAL_RUN_GUIDE.md](GRAND_CANONICAL_RUN_GUIDE.md).

### 8.8 Measurements that did not converge

The full one-body and density correlation arrays were sampled only once per 10,000 worm excursions. The 64 x 64 blocks therefore have roughly one hundred such measurements and effectively unusable errors for a long-distance algebraic fit. No $\eta=1/4$ or correlation-length claim should be extracted from those files.

In addition to sparse sampling, the final exact-identity regression found a specific on-site estimator defect. For a canonical 3 x 3 system with $N=13$, exact number conservation requires

$$
g_1(0)=\langle b_i^\dagger b_i\rangle=\frac{N}{9}=1.444444\ldots,
$$

whereas one seed returned `Density_Matrix[0] = 1.42728 +/- 0.00167`. A fresh independent seed returned `1.43226 +/- 0.00158`, still about 7.7 reported standard errors low. The independent reproduction confirms an estimator defect rather than one unlucky chain. Do not normalize a correlation function by element zero or use that element as a density check. The off-site elements use a different estimator branch and agreed with exact diagonalisation in the small regression, but the sparse 64 x 64 arrays are still inadequate for exponent fitting.

To use correlations as a second BKT estimator, decrease `Nmeasure2` only after benchmarking the extra $O(L^2)$ cost, then collect enough independent blocks to stabilize $g_1(r)$ at large $r$. Winding remains substantially cheaper and is the primary observable in the workflow below.

## 9. Commands in order: one trustworthy 64 x 64 point

This section is deliberately literal. Run the commands in a normal 64-bit PowerShell window. Markdown backticks are only formatting; do not type them.

### 9.1 Open the repository

```powershell
Set-Location "\\alfs1.physics.ox.ac.uk\al\kuo\projects\worm"
```

The build and run scripts resolve their own absolute paths. No MSYS2 activation command is required: the build script names CMake and MSYS2 tools explicitly, and the run script prepends the packaged DLL directory. Microsoft MPI must already be installed because `mpiexec.exe` and `msmpi.dll` are host prerequisites.

### 9.2 Build once after source changes

```powershell
.\scripts\build_windows.ps1 -Jobs 32
```

Use 32 jobs because compilation benefits from all 32 logical processors on the i9-13900. A successful command creates `dist\windows-square` with the four executables, DLL closure, and checksum manifest. Parameter changes do not require rebuilding.

To verify every packaged file against the recorded hashes:

```powershell
Push-Location .\dist\windows-square
try {
    foreach ($line in Get-Content .\SHA256SUMS.txt) {
        if ($line -notmatch '^([0-9A-Fa-f]{64}) \*(.+)$') { throw "Malformed manifest line: $line" }
        $expected = $Matches[1]
        $file = $Matches[2]
        $actual = (Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash
        if ($actual -ne $expected) { throw "Hash mismatch: $file" }
    }
    Write-Host "All package hashes match."
}
finally {
    Pop-Location
}
```

### 9.3 Make a point-specific parameter file

Keep the audited template unchanged and copy it to a name containing size, temperature, and base seed:

```powershell
Copy-Item .\parameter_files\Rb87_BKT_64.ini .\parameter_files\Rb87_L64_T40_seed96040.ini
notepad .\parameter_files\Rb87_L64_T40_seed96040.ini
```

For the first audited point, retain:

```text
Lx = 64
Ly = 64
beta = 0.279077157505874
U_on = 0.152
canonical = 9948
nmax = 16
canonical_window = 0.5
thermalization = 50000
runtimelimit = 900
seed = 96040
```

Save and close Notepad. Use a new output directory whenever you change any physical/Markov-chain parameter or base seed. Never edit parameters between `Fresh`, `Resume`, and `Production` calls for one checkpoint family.

### 9.4 Start the chain

```powershell
.\scripts\run_rb87_bkt64.ps1 -Mode Fresh -Processes 24 -ParameterFile .\parameter_files\Rb87_L64_T40_seed96040.ini -OutputDirectory .\results\Rb87_L64_T40_seed96040
```

Twenty-four ranks use the 24 physical cores. Each rank is one full independent Markov chain; this is not a domain decomposition. `Fresh` refuses to replace an existing output/checkpoint family. If that happens unintentionally, choose a genuinely new output directory rather than deleting evidence from an earlier run.

Treat all output from this invocation as burn-in.

### 9.5 Continue burn-in without resetting the sweep counter

Run three plain continuations:

```powershell
.\scripts\run_rb87_bkt64.ps1 -Mode Resume -Processes 24 -ParameterFile .\parameter_files\Rb87_L64_T40_seed96040.ini -OutputDirectory .\results\Rb87_L64_T40_seed96040
```

```powershell
.\scripts\run_rb87_bkt64.ps1 -Mode Resume -Processes 24 -ParameterFile .\parameter_files\Rb87_L64_T40_seed96040.ini -OutputDirectory .\results\Rb87_L64_T40_seed96040
```

```powershell
.\scripts\run_rb87_bkt64.ps1 -Mode Resume -Processes 24 -ParameterFile .\parameter_files\Rb87_L64_T40_seed96040.ini -OutputDirectory .\results\Rb87_L64_T40_seed96040
```

These commands retain both chain state and accumulated counters. Their result files are not production data. Their purpose is to let all ranks exceed the nominal 50,000-sweep threshold and give winding sectors time to reorganize before any measurement reset. Three is the safer audited starting count for this L = 64 point; block agreement remains the real criterion.

Always continue with the same 24 ranks. The runner checks `Rb87_BKT_64.clone.h5` and all `.1` through `.23` rank checkpoints before launching.

### 9.6 Collect separate production blocks

Run a reset production block:

```powershell
.\scripts\run_rb87_bkt64.ps1 -Mode Production -Processes 24 -ParameterFile .\parameter_files\Rb87_L64_T40_seed96040.ini -OutputDirectory .\results\Rb87_L64_T40_seed96040
Copy-Item .\results\Rb87_L64_T40_seed96040\Rb87_BKT_64.out.h5 .\results\Rb87_L64_T40_seed96040\production_01.out.h5
```

Then collect a second block:

```powershell
.\scripts\run_rb87_bkt64.ps1 -Mode Production -Processes 24 -ParameterFile .\parameter_files\Rb87_L64_T40_seed96040.ini -OutputDirectory .\results\Rb87_L64_T40_seed96040
Copy-Item .\results\Rb87_L64_T40_seed96040\Rb87_BKT_64.out.h5 .\results\Rb87_L64_T40_seed96040\production_02.out.h5
```

Then a third:

```powershell
.\scripts\run_rb87_bkt64.ps1 -Mode Production -Processes 24 -ParameterFile .\parameter_files\Rb87_L64_T40_seed96040.ini -OutputDirectory .\results\Rb87_L64_T40_seed96040
Copy-Item .\results\Rb87_L64_T40_seed96040\Rb87_BKT_64.out.h5 .\results\Rb87_L64_T40_seed96040\production_03.out.h5
```

Each call starts from the previous chain state but clears the statistics. The explicit copies give stable, human-readable block names; the runner also puts timestamped previous results in its `archive` subdirectory. Three blocks are a starting minimum, not an automatic declaration of convergence. If the last two or three differ beyond their errors, continue producing blocks and discard the drifting prefix.

Agreement of continued blocks is necessary but not sufficient. The audit found that an apparently well-sampled first production block could still drift and that rank scatter exceeded the pooled within-chain error. After the first family appears stable, repeat Sections 9.3-9.7 with a different `seed` and a new parameter filename/output directory. Treat agreement between fresh base-seed families, plus acceptable scatter among their ranks, as the convergence requirement. Never point a new seed at the old checkpoint directory.

### 9.7 Read the key observables

Use the HDF5 utility already installed with MSYS2:

```powershell
$h5dump = "C:\msys64\ucrt64\bin\h5dump.exe"
& $h5dump -d "/simulation/results/Winding_number_squared/mean/value" ".\results\Rb87_L64_T40_seed96040\production_03.out.h5"
& $h5dump -d "/simulation/results/Winding_number_squared/mean/error" ".\results\Rb87_L64_T40_seed96040\production_03.out.h5"
& $h5dump -d "/simulation/results/Winding_number_squared/tau" ".\results\Rb87_L64_T40_seed96040\production_03.out.h5"
& $h5dump -d "/simulation/results/Number_of_particles/mean/value" ".\results\Rb87_L64_T40_seed96040\production_03.out.h5"
& $h5dump -d "/simulation/results/Total_Energy/mean/value" ".\results\Rb87_L64_T40_seed96040\production_03.out.h5"
```

Interpret the two winding entries as $\langle W_x^2\rangle$ and $\langle W_y^2\rangle$, and add them. For a quick approximate error on the sum, add the component errors in quadrature; a final analysis should retain their covariance or bootstrap independent blocks/rank groups.

A production point passes the basic audit only if:

1. measured particle number is exactly the requested integer within file representation;
2. successive production blocks agree in total/kinetic/potential energy;
3. $W_x^2$ and $W_y^2$ agree statistically on the square lattice;
4. the winding sum is stable across the last two or three blocks;
5. multiple ranks and a second base-seed group visit winding sectors; and
6. the result remains stable under any relevant `nmax`, window, or run-length check.

### 9.8 Make an independent replica group

Internal MPI ranks are useful independent chains, but a separate base seed and output family provide a stronger end-to-end check. Create a second file:

```powershell
Copy-Item .\parameter_files\Rb87_L64_T40_seed96040.ini .\parameter_files\Rb87_L64_T40_seed196040.ini
notepad .\parameter_files\Rb87_L64_T40_seed196040.ini
```

Change only:

```text
seed = 196040
```

Then repeat Sections 9.4-9.7 with parameter file `Rb87_L64_T40_seed196040.ini` and output directory `results\Rb87_L64_T40_seed196040`. Agreement between the equilibrated seed groups is stronger evidence than many measurements from one chain history.

## 10. Commands in order: locate and resolve the BKT transition

### 10.1 Generate a coarse map cheaply

First generate, but do not run, a low/central/high scan:

```powershell
.\scripts\run_temperature_sweep.ps1 -ParameterFile .\parameter_files\Rb87_BKT_64.ini -Temperatures 2.866590753,3.224914597,3.583238441,3.941562285,4.461131859 -Sizes 16,24,32 -Canonical -Density 2.428635525 -Sweeps 4800000 -Thermalization 50000 -RuntimeLimit 900 -Processes 8 -SweepName rb87_bkt_coarse_seed96040 -DryRun
```

This creates complete per-point INI files and `results\rb87_bkt_coarse_seed96040\sweep_manifest.csv`. The temperatures correspond to 32, 36, 40, 44, and 49.8 nK. Small sizes and eight ranks are suitable for reconnaissance, where the goal is only to locate the rapid winding crossover at modest cost.

To launch those one-block reconnaissance jobs, repeat the command without `-DryRun`. Reusing the same sweep name is safe before any HDF5 point outputs exist; the generated INI files and manifest are refreshed:

```powershell
.\scripts\run_temperature_sweep.ps1 -ParameterFile .\parameter_files\Rb87_BKT_64.ini -Temperatures 2.866590753,3.224914597,3.583238441,3.941562285,4.461131859 -Sizes 16,24,32 -Canonical -Density 2.428635525 -Sweeps 4800000 -Thermalization 50000 -RuntimeLimit 900 -Processes 8 -SweepName rb87_bkt_coarse_seed96040
```

One automatically generated block per point is not converged production. It is allowed to tell you that 32 nK is clearly low, 49.8 nK clearly high, and the interesting interval is near 37-42 nK. Do not feed these preliminary points directly into a precision BKT fit.

### 10.2 Generate the finite-size fine scan

Generate parameter files around the transition:

```powershell
.\scripts\run_temperature_sweep.ps1 -ParameterFile .\parameter_files\Rb87_BKT_64.ini -Temperatures 3.314495558,3.404076519,3.493657480,3.583238441,3.672819402,3.762400363 -Sizes 24,32,48,64 -Canonical -Density 2.428635525 -Sweeps 4800000 -Thermalization 50000 -RuntimeLimit 900 -Processes 24 -SweepName rb87_bkt_fine_seed96040 -DryRun
```

These are 37, 38, 39, 40, 41, and 42 nK. `-DryRun` is intentional: inspect each generated INI, then run each selected point through the dedicated `Fresh -> Resume -> Resume -> Resume -> Production...` protocol, with its own output directory. This is slower operationally but preserves checkpoints and independent production blocks, which a one-shot sweep cannot certify.

List the generated files with:

```powershell
Get-ChildItem .\results\rb87_bkt_fine_seed96040\parameters\*.ini
```

For each selected INI, run the same commands as Section 9 while changing `-ParameterFile` and assigning a unique `-OutputDirectory`. For example, after identifying the generated $L=48$, 40 nK INI path:

```powershell
.\scripts\run_rb87_bkt64.ps1 -Mode Fresh -Processes 24 -ParameterFile ".\results\rb87_bkt_fine_seed96040\parameters\rb87_bkt_fine_seed96040_L048_T3p583238.ini" -OutputDirectory .\results\production\L48_T40_seed96040
```

For a comparable L = 64 start, follow it with at least three `Resume` calls and at least three `Production` calls using those same two path arguments. Smaller sizes may need less burn-in, but decide that from block and seed agreement rather than silently shortening the protocol. Generated job names round the temperature label to six decimals; use `Get-ChildItem` rather than guessing if a filename differs.

Repeat the fine grid with a second base seed, for example `-BaseSeed 196040` and `-SweepName rb87_bkt_fine_seed196040`, then run the selected files into separate output directories.

### 10.3 Why these sizes and temperatures work

The coarse runs already place 32 nK in the large-winding regime and 49.8 nK in the collapsed-winding regime. A 1 nK spacing from 37 to 42 nK resolves the interval in which $W_\Sigma$ moves through the universal-jump scale without spending the full budget far from the transition.

Sizes 24 and 32 are cheap enough to map the curve, while 48 and 64 expose the slow logarithmic BKT size dependence. At least three sizes are required to fit a nonuniversal logarithmic scale $L_0$; four make sensitivity to dropping the smallest size visible. An $L=80$ check is valuable if resources permit, using $N=15543$, but it should be added only after the 24-64 data and chain mixing are under control.

Using density rather than a fixed $N$ prevents the finite-size scan from inadvertently changing the thermodynamic state. Using exact $N$ prevents the grand-canonical number-sector barrier observed in the audit. Periodic boundaries make winding a well-defined bulk response without edges.

## 11. Turning the scan into a BKT estimate

### 11.1 Build a clean data table

For every $(L,T,\text{base seed})$, record only equilibrated reset blocks. A useful analysis table has these columns:

```text
L, T_nK, beta, N, base_seed, block,
sample_count, total_energy, kinetic_energy,
Wx2, Wx2_error, Wx2_tau, Wy2, Wy2_error, Wy2_tau
```

Add derived columns

$$
W_\Sigma=W_x^2+W_y^2,
\qquad
\rho_s=\frac{W_\Sigma}{2\beta}.
$$

Do not average away the block or seed labels at import time. First plot every block in chronological order. Energy usually equilibrates before winding, so a stable energy trace alone is insufficient. Reject only a clearly identified drifting prefix; do not selectively remove isolated high or low values after seeing the desired transition.

Treat reset blocks and independent seed groups as the primary independent units. Millions of per-worm samples have autocorrelation and share one evolving worldline history. A hierarchical bootstrap over seed groups and accepted blocks is more honest than pretending every recorded sweep is independent.

### 11.2 Visual diagnostics before fitting

For each size, plot $W_\Sigma(T)$ with the horizontal line $4/\pi$. The low-temperature side should have larger winding and the high-temperature side should collapse toward zero. Increasing $L$ should sharpen the crossover, but finite-size critical values need not cross exactly at $4/\pi$.

Also plot:

- $W_x^2-W_y^2$, which should fluctuate around zero;

- energy per site versus block number;

- winding autocorrelation estimate versus $L$ and $T$;

- independent-seed differences normalized by their errors; and

- the fraction of runs/blocks that actually visit nonzero winding sectors.

A point with tiny printed error but poor replica agreement is not converged. Increase burn-in, production length, or the number of independent seeds before adding it to the finite-size fit.

### 11.3 Weber-Minnhagen finite-size fit

At each trial temperature, fit the size dependence to

$$
\frac{\pi W_\Sigma(L,T)}{4}
=A(T)\left[1+\frac{1}{2\ln(L/L_0(T))}\right],
$$

using at least three sizes. At the BKT transition, the thermodynamic amplitude satisfies $A(T_{\rm BKT})=1$. Below it, the limiting stiffness is larger; above it, the BKT critical form ceases to describe large sizes as winding collapses.

The practical procedure is:

1. Fit $A$ and $L_0$ independently at each sampled temperature.
2. Inspect residuals and goodness of fit; a numerical optimizer returning parameters is not enough.
3. Interpolate $A(T)$ through one near unity to estimate $T_{\rm BKT}$.
4. Repeat while dropping the smallest size, because subleading corrections are largest there.
5. Bootstrap whole blocks/replica groups and repeat the interpolation to obtain statistical uncertainty.
6. Quote the spread from fit ranges, temperature interpolation, cutoff checks, and lattice spacing as separate systematic uncertainties.

Because BKT corrections are logarithmic, conventional power-law finite-size extrapolation can look smooth while giving the wrong critical temperature. The logarithm is not an optional cosmetic correction; it is a consequence of the marginal renormalization-group flow at the transition.

### 11.4 Optional correlation-function cross-check

With much denser `Nmeasure2` sampling, fit the long-distance one-body density matrix to

$$
g_1(r)\sim r^{-\eta(T)}
$$

below the transition, with $\eta(T_{\rm BKT})=1/4$. On a periodic finite square, use a periodic-distance or conformal finite-size form rather than fitting all raw points to a simple infinite-line power law. Exclude ultraviolet distances of a few grid spacings and test fit-window stability.

This is an independent physical signature, but it is more expensive here than winding. The current audited files do not have enough correlation measurements for it.

### 11.5 What would justify a transition claim

A defensible result would contain all of the following:

- a low/high-temperature bracket and a fine grid through the critical region;

- at least three well-separated sizes, preferably 24, 32, 48, and 64 or larger;

- at least two independent base seeds per important point;

- stable post-burn-in production blocks with exact $N$, isotropic winding, and stable energy;

- a finite-size logarithmic fit whose inferred $T_{\rm BKT}$ is stable to removing the smallest size;

- occupation-cutoff and run-length checks;

- a spatial-grid-spacing systematic before quantitative comparison to continuum Rb-87; and

- clear separation of statistical uncertainty from mapping/model systematics.

The current audit verifies the build, ensemble, Hamiltonian, energy estimator, winding estimator, and retained checkpoint invariants. It independently replicates the late L = 64, 40 nK energy and winding point and demonstrates low/near/high-temperature BKT-like winding. It still lacks the multi-size, multi-temperature, replicated data set and continuum-spacing study required by the list above. It therefore must not be presented as a precision measurement of $T_{\rm BKT}$.

## 12. Common mistakes and recovery

### Rebuilding after every parameter edit

Do not rebuild. The INI file is parsed at runtime. Rebuild only after source, compiler, dependency, or CMake-option changes.

### Starting production immediately after a short fresh run

On $L=64$, the fresh 900-second block did not bring every rank to 50,000 sweeps, and the independent rerun's first production still drifted after two resumes. Use at least three `Resume` burn-in blocks before the first reset `Production` for a comparable start, then judge actual block stability.

### Editing the INI while continuing old checkpoints

Do not do this. A checkpoint contains the established chain and restored parameters. Start a new output directory and `Fresh` run for any changed physics, cutoff, canonical window, lattice, or seed.

### Changing the MPI rank count

Do not change it within one checkpoint family. Rank-specific files are part of the state. Start a fresh family if the process count must change.

### Treating the sweep tool's temperature as beta

The sweep tool takes $T/t$; the INI takes $\beta t$. They are reciprocals. At 40 nK, pass `3.583238441` to `-Temperatures`, but put approximately `0.279077158` in `beta`.

### Comparing one finite size directly to 4/pi

The universal jump is thermodynamic and has logarithmic finite-size corrections. One crossing gives a useful bracket, not $T_{\rm BKT}$.

### Trusting formal error bars while blocks drift

Within-chain errors quantify fluctuations after the chain reaches the sampled region. They do not detect an unforgotten initial stripe, trapped number sector, or slowly drifting winding distribution. Compare chronology, blocks, and seeds.

### Trusting only the combined MPI mean

Each MPI rank is a complete chain. Inspect agreement among rank or independent-seed means as well as the combined result. A persistent rank-wide offset is sampled only once by that rank, even though it contributes many correlated log-binning blocks.

### Using `canonical_window` at or above one

The modified executable rejects it. Such a window can permit the closed integer sector to change by one and is not the strict exact-$N$ calculation documented here.

### Interpreting `mu` as the density control in this run

The exact `canonical` value controls density. At fixed $N$, `mu` adds a common constant to closed-state energies and does not alter their normalized relative probabilities. It does shift the raw potential and total energies, and it can change open-worm efficiency; it is therefore neither a density control nor a meaningless parser placeholder.

### Claiming an algebraic exponent from the present correlation arrays

`Nmeasure2=10000` leaves too few expensive measurements. Reduce it, benchmark cost, and collect independent long blocks before fitting $\eta$.

## 13. Model and verification boundaries

This calculation includes a homogeneous single-layer contact-interacting Bose gas approximated on a periodic square lattice. It does not include Beregi's trap shape, hard/soft box wall, bilayer tunnel coupling, imaging response, disorder, or experimental nonequilibrium protocol.

The 0.5 micrometre spatial grid is coarse relative to the 40 nK thermal wavelength. The inferred density uses a weak-gas critical-density relation rather than a measured equation-of-state fit. These are larger conceptual systematics than the final decimal places printed in `beta`.

The validated executable has optional Matsubara measurements disabled, and the sparse current correlation sampling does not support a Green-function exponent. The current on-site density-matrix estimator is known to fail its exact $g_1(0)=N/N_{\rm sites}$ identity, and the optional direct Matsubara tau histogram is not checkpoint-complete. Neither issue affects the audited energy/winding path. The complete POSIX-oriented ALPSCore test suite was not run under MinGW, although the pinned source, patch application, packaged dependencies, hashes, error guards, exact-diagonalisation regressions, serial/MPI smoke tests, checkpoint continuation, and worldline-invariant reconstruction were checked.

## 14. References

- Nicolas Sadoune and Lode Pollet, [Efficient and scalable Path Integral Monte Carlo Simulations with worm-type updates for Bose-Hubbard and XXZ models](https://arxiv.org/abs/2204.12262).

- [LodePollet/worm source repository](https://github.com/LodePollet/worm).

- N. V. Prokof'ev, B. V. Svistunov, and I. S. Tupitsyn, [Exact, complete, and universal continuous-time worldline Monte Carlo approach to the statistics of discrete quantum systems](https://doi.org/10.1016/S0375-9601(97)00957-2).

- M. Boninsegni, N. V. Prokof'ev, and B. V. Svistunov, [Worm algorithm for continuous-space path integral Monte Carlo simulations](https://doi.org/10.1103/PhysRevE.74.036701).

- E. L. Pollock and D. M. Ceperley, [Path-integral computation of superfluid densities](https://doi.org/10.1103/PhysRevB.36.8343).

- D. R. Nelson and J. M. Kosterlitz, [Universal Jump in the Superfluid Density of Two-Dimensional Superfluids](https://doi.org/10.1103/PhysRevLett.39.1201).

- H. Weber and P. Minnhagen, [Monte Carlo determination of the critical temperature for the two-dimensional XY model](https://doi.org/10.1103/PhysRevB.37.5986).

- Y.-D. Hsieh, Y.-J. Kao, and A. W. Sandvik, [Finite-size scaling method for the Berezinskii-Kosterlitz-Thouless transition](https://arxiv.org/abs/1302.2900).

- Abel Beregi, [Probing universality of 2D quantum systems with bilayer Bose gases](https://ora.ox.ac.uk/objects/uuid:b2f4f0a1-8576-4528-bbd3-557d273cfbdd), Oxford DPhil thesis (2024).

- Daniel A. Steck, [Rubidium 87 D Line Data](https://steck.us/alkalidata/rubidium87numbers.pdf).

## 15. Minimal final checklist

Before using any number in a report, confirm:

```text
[ ] Pinned build completed and the canonical MPI executable is used.
[ ] One immutable INI and one unique output directory define each chain family.
[ ] L, beta, U/t, exact N, nmax, and canonical_window are recorded.
[ ] Fresh was followed by sufficient plain Resume burn-in before reset statistics.
[ ] The MPI rank count stayed fixed for all continuations.
[ ] At least three equilibrated Production blocks and two base seeds agree.
[ ] N is exact, Wx and Wy agree, energy is stable, and winding sectors are visited.
[ ] Several temperatures bracket the collapse of winding.
[ ] Several sizes are fitted with BKT logarithmic finite-size scaling.
[ ] Cutoff, fit-range, and spatial-grid systematics are reported.
```
