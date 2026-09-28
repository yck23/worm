# Grand-canonical Bose-Hubbard and BKT guide

This is the guide to use when **chemical potential is fixed and particle number is measured**.

Audit status: 28 September 2026. The grand-canonical implementation, checkpoint restore, independent-start protocol, and number response have passed. A fixed-mu temperature audit now covers L = 16, 24, and 32 with three independent starts per point and repeated reset production blocks. It resolves normal, crossover, and superfluid winding regimes, while also showing residual slow number-sector motion near the crossover; it is a validated diagnostic workflow, not a precision transition-temperature result. The separately audited 64 x 64 reference calculation remains canonical, because a successful smaller lattice does not prove that an L = 64 grand-canonical chain has equilibrated.

If Markdown punctuation is visible, open the formatted preview in VS Code with Ctrl+Shift+V. The command blocks in this guide use tildes instead of backtick fences so the raw file is less distracting.

## 1. The basic idea

There are two different questions one can ask.

| Question | Ensemble | User supplies | Simulation measures |
|---|---|---|---|
| What happens at exactly N particles? | Canonical | N and temperature | Energy and winding; N stays exact |
| What density is selected by a reservoir? | Grand canonical | Chemical potential mu and temperature | Mean N, fluctuations, energy, and winding |

For a grand-canonical run:

- set canonical = -1;
- run qmc_worm_mpi.exe, selected by the runner option -Ensemble GrandCanonical;
- choose mu and beta in the parameter file;
- read the resulting mean particle number from Number_of_particles.

The line initial_particle_number is **only a starting guess**. It does not fix N and it does not overwrite mu.

The shortest correct statement is:

~~~text
input:  mu, beta, L, t, U, ...
output: mean N and Var(N)
mean density = mean N / (Lx times Ly)
~~~

## 2. Why particle number fluctuates

The Bose-Hubbard Hamiltonian used here is

$$
H_\mu =
-t\sum_{\langle ij\rangle}(b_i^\dagger b_j+b_j^\dagger b_i)
+\frac{U}{2}\sum_i n_i(n_i-1)
+V\sum_{\langle ij\rangle}n_i n_j
-\mu\sum_i n_i .
$$

It is useful to call the hopping-plus-interaction part without the final term H0. Then

$$
H_\mu=H_0-\mu N
$$

and the grand partition function is

$$
\Xi(\mu,T)=\sum_N {\rm Tr}_N\,
\exp[-\beta(H_0-\mu N)] .
$$

Increasing mu makes larger-N sectors statistically more favourable. The actual mean N also depends on beta, U, t, lattice size, boundary conditions, and the local occupation cutoff.

The response identity is

$$
\frac{\partial\langle N\rangle}{\partial\mu}
=\beta\left(\langle N^2\rangle-\langle N\rangle^2\right).
$$

The program measures both N and N squared, so this identity can be checked numerically. Per site, the quantity reported by the helper is

$$
\chi_N/L^2=\frac{\beta\,{\rm Var}(N)}{L^2}.
$$

This is the number susceptibility per site. Some literature defines an isothermal compressibility with an additional density normalization, so state the convention rather than using the word compressibility alone.

## 3. What N means in the worldline algorithm

The finite-temperature path integral represents bosons as occupation-number worldlines wrapped around imaginary time from 0 to beta.

- A closed configuration contributes to the partition function.
- A hopping event is a spatial corner, or kink, in a worldline.
- The number of worldlines crossing a complete imaginary-time slice is the particle number N.
- A worm update temporarily opens a line, moves one endpoint through space and imaginary time, inserts or deletes kinks, and finally closes the line again.

In the ordinary grand-canonical executable, a worm can cross the imaginary-time boundary before closing. That changes the integer number of closed worldlines and therefore moves the Markov chain from N to N plus or minus one. The source updates number_of_particles at that boundary crossing and measures N only after the update has returned to a closed sector.

The open worm is sampling machinery, not a physical configuration whose ordinary energy is retained. One source-level sweep is one complete update call. It may contain a long worm excursion, so a sweep is neither one local proposal nor physical time.

## 4. How BKT physics appears

With periodic x and y boundaries, worldlines can wrap around the spatial torus. The code measures the squared winding components

$$
\langle W_x^2\rangle,\qquad \langle W_y^2\rangle .
$$

For t = 1 in two dimensions, the finite-size helicity modulus is

$$
\rho_s(L,T)=\frac{\langle W_x^2+W_y^2\rangle}{2\beta}.
$$

The thermodynamic Nelson-Kosterlitz jump gives the reference

$$
W_\Sigma \equiv \langle W_x^2+W_y^2\rangle
=\frac{4}{\pi}=1.273239545\ldots .
$$

This reference is a useful diagnostic, not a one-size definition of the transition. A finite lattice has a smooth crossover and logarithmic BKT size corrections. A defensible transition temperature requires several sizes, a fine temperature grid, independent chains, and a finite-size scaling analysis. The [Nelson-Kosterlitz paper](https://doi.org/10.1103/PhysRevLett.39.1201) establishes the universal jump, and [Hsieh, Kao, and Sandvik](https://arxiv.org/abs/1302.2900) discuss finite-size BKT scaling.

The original implementation is the continuous-time worm code described by [Sadoune and Pollet](https://arxiv.org/abs/2204.12262). The Rb-87 physical starting point is motivated by Abel Beregi's [Oxford DPhil thesis](https://ora.ox.ac.uk/objects/uuid:b2f4f0a1-8576-4528-bbd3-557d273cfbdd). The weak-gas critical-density comparison comes from [Prokof'ev, Ruebenacker, and Svistunov](https://arxiv.org/abs/cond-mat/0106075).

## 5. Files and what they do

| File | Purpose |
|---|---|
| [parameter_files/Rb87_BKT_64_grand_canonical.ini](parameter_files/Rb87_BKT_64_grand_canonical.ini) | Editable one-point template |
| [scripts/run_rb87_bkt64.ps1](scripts/run_rb87_bkt64.ps1) | Runs one parameter file and one chain family |
| [scripts/run_gc_bkt_audit.ps1](scripts/run_gc_bkt_audit.ps1) | Generates and runs low/mid/high starts for a size-temperature grid |
| [scripts/check_gc_burnin.py](scripts/check_gc_burnin.py) | Reads current discarded burn-in means and checkpoint endpoint sectors |
| [scripts/analyze_gc_bkt_audit.py](scripts/analyze_gc_bkt_audit.py) | Combines only numbered retained production blocks |
| [scripts/analyze_gc_response.py](scripts/analyze_gc_response.py) | Recomputes the two-mu fluctuation-response check from retained blocks |
| [scripts/summarize_bkt_result.py](scripts/summarize_bkt_result.py) | Explains one HDF5 result |

Editing an INI file does not require rebuilding. Rebuild only after C++ source, compiler, or ALPSCore changes.

## 6. Important parameter-file lines

The current one-point template contains the following physical controls:

~~~text
Lx = 64
Ly = 64
t_hop = 1.0
U_on = 0.152
V_nn = 0.0
mu = -3.53
beta = 0.279077157505874
nmax = 16
canonical = -1
initial_particle_number = 9948
~~~

Their meanings are:

- **Lx and Ly** set the number of spatial grid sites.
- **t_hop = 1** defines the energy unit.
- **U_on** is the on-site repulsion in units of t.
- **mu** is the reservoir chemical potential in units of t. It controls the equilibrium number distribution.
- **beta** is inverse temperature in units of 1/t.
- **nmax** is the largest allowed occupation on one site.
- **canonical = -1** selects the grand-canonical behavior of the ordinary executable.
- **initial_particle_number** constructs only the fresh starting configuration. It never constrains later N.
- **canonical_window** may remain in the shared template, but the ordinary grand-canonical executable ignores it. It neither limits nor fixes N.

The value mu = -3.53 is also physically sensible for the intended 40 nK
critical region. Inserting the weak-gas constant xi_mu = 13.2 and coupling
g-tilde = 0.076 into the Prokof'ev-Svistunov critical-chemical-potential
relation gives mu_cont/t about 0.447. The square-grid shift
mu_BH = mu_cont - 4t then gives mu_BH/t about -3.553. This analytic estimate
motivates the scan; the measured lattice equation of state, cutoff tests, and
finite-size scaling still determine the numerical result.

The main numerical controls are:

~~~text
runtimelimit = 900
sweeps = 4800000
thermalization = 1000
Nmeasure = 1
Nmeasure2 = 1000
C_worm = 2.0
p_moveworm = 0.4
p_insertkink = 0.3
p_deletekink = 0.2
p_glueworm = 0.1
~~~

The value thermalization = 1000 is a technical accumulator threshold, **not a claim that 1000 worm excursions equilibrate the physics**. The workflow discards complete Fresh and Resume invocations as explicit burn-in. This also ensures that all MPI ranks acquire the sparse spatial observables before ALPS attempts to merge them.
The current executable refuses a Production reset if a rank has not genuinely passed this threshold; use Resume first. Passing the threshold is still only a mechanical prerequisite, not evidence of physical equilibration.

The wall-time limit is soft. The stop condition is checked after a complete worm update, so a long near-critical worm can carry the process beyond the nominal number of seconds.

Do not tune the proposal probabilities merely to make a run appear faster. They affect autocorrelation and sector mobility. Any change requires a new convergence audit.

## 7. What the command line does and does not override

For the one-point runner, physics values come from the chosen parameter file. The command line selects:

- the canonical or grand-canonical executable;
- Fresh, Resume, or Production mode;
- MPI process count;
- input path;
- output directory.

It does not silently replace mu, beta, U, lattice size, or the worm probabilities.

For the grid audit runner, switches such as -Mu and -Betas are used once to **generate explicit INI files** under the audit result directory. Those generated files are the record of what ran.

Fresh reads an INI and creates checkpoints. Resume and Production restore the parameters stored in those checkpoints. Therefore:

- edit the parameter file before Fresh;
- never expect an edit to change an existing checkpoint;
- any changed mu, beta, size, seed, or model parameter needs a new output directory and a new Fresh run;
- keep the MPI process count unchanged when restoring rank-specific checkpoints.

## 8. Run one editable parameter file

Open PowerShell and enter:

~~~powershell
Set-Location '\\alfs1.physics.ox.ac.uk\al\kuo\projects\worm'
notepad .\parameter_files\Rb87_BKT_64_grand_canonical.ini
~~~

Save the file, then choose a new output directory:

~~~powershell
$parameterFile = '.\parameter_files\Rb87_BKT_64_grand_canonical.ini'
$outputDirectory = '.\results\my_grand_canonical_run'
~~~

Start a new chain:

~~~powershell
.\scripts\run_rb87_bkt64.ps1 -Ensemble GrandCanonical -Mode Fresh -Processes 24 -ParameterFile $parameterFile -OutputDirectory $outputDirectory
~~~

Continue without resetting measurements:

~~~powershell
.\scripts\run_rb87_bkt64.ps1 -Ensemble GrandCanonical -Mode Resume -Processes 24 -ParameterFile $parameterFile -OutputDirectory $outputDirectory
~~~

Use Resume for discarded burn-in. Inspect the evolving state:

~~~powershell
python .\scripts\summarize_bkt_result.py (Join-Path $outputDirectory 'Rb87_BKT_64_grand_canonical.out.h5')
~~~

Once independent starts and late blocks are stable, reset the accumulators and retain a production block:

~~~powershell
.\scripts\run_rb87_bkt64.ps1 -Ensemble GrandCanonical -Mode Production -Processes 24 -ParameterFile $parameterFile -OutputDirectory $outputDirectory
Copy-Item (Join-Path $outputDirectory 'Rb87_BKT_64_grand_canonical.out.h5') (Join-Path $outputDirectory 'production_01.out.h5')
~~~

Repeat Production with production_02, production_03, and so on. Production continues the worldlines and RNG state but resets measured statistics. It does not create a statistically independent fresh chain.

## 9. Recommended independent-start temperature audit

The audit helper makes three starts for every point:

| Family | Default starting density |
|---|---:|
| Low | 1.8 |
| Mid | 2.4 |
| High | 3.0 |

These values are deliberately separated. Agreement later demonstrates that the result is not inherited from initial_particle_number.

The following command creates and runs a 32 x 32 fixed-mu grid. Twelve jobs times two MPI ranks use all 24 physical desktop cores:

~~~powershell
Set-Location '\\alfs1.physics.ox.ac.uk\al\kuo\projects\worm'
.\scripts\run_gc_bkt_audit.ps1 -Mode Fresh -AuditName my_gc_L32 -Sizes 32 -Betas 0.24,0.27,0.29,0.32 -Mu -3.53 -ProcessesPerChain 2 -RuntimeLimit 600
~~~

For the present Rb-87 grid, where t/kB = 11.1630863 nK, these points mean:

| beta | T/t | Physical temperature |
|---:|---:|---:|
| 0.24 | 4.1667 | 46.51 nK |
| 0.27 | 3.7037 | 41.35 nK |
| 0.29 | 3.4483 | 38.49 nK |
| 0.32 | 3.1250 | 34.88 nK |

Temperature increases when beta decreases. The selected range therefore runs
from the colder expected superfluid side to the hotter expected normal side.

Treat that entire call as burn-in. Inspect it:

~~~powershell
py -3.13 .\scripts\check_gc_burnin.py .\results\my_gc_L32
~~~

Continue all checkpoints, still as discarded burn-in:

~~~powershell
.\scripts\run_gc_bkt_audit.ps1 -Mode Resume -AuditName my_gc_L32
py -3.13 .\scripts\check_gc_burnin.py .\results\my_gc_L32
~~~

Do not demand identical checkpoint endpoints. Each endpoint is one draw from a broad equilibrium number distribution. What matters is that low, middle, and high families occupy overlapping ranges and no longer preserve separate bands tied to their starts.

Now make at least three reset blocks. Near the crossover, five or longer blocks
may be needed because particle-number autocorrelation is slow:

~~~powershell
.\scripts\run_gc_bkt_audit.ps1 -Mode Production -AuditName my_gc_L32 -ProductionIndex 1
.\scripts\run_gc_bkt_audit.ps1 -Mode Production -AuditName my_gc_L32 -ProductionIndex 2
.\scripts\run_gc_bkt_audit.ps1 -Mode Production -AuditName my_gc_L32 -ProductionIndex 3
~~~

Analyze only the numbered blocks:

~~~powershell
py -3.13 .\scripts\analyze_gc_bkt_audit.py .\results\my_gc_L32
~~~

The analysis directory contains:

- audit_blocks.csv: every retained block;
- audit_points.csv: averages built from independent-start family means;
- AUDIT_SUMMARY.md: readable table and raw crossing diagnostic;
- audit_summary.png: winding and density against beta, if matplotlib is available.

The analyzer deliberately does not read the unnumbered result.out.h5 burn-in file.

### Add L = 16 and L = 24 without oversubscribing the desktop

With one MPI process per independent family, two sizes times four beta values
times three starts use the 24 physical cores exactly:

~~~powershell
.\scripts\run_gc_bkt_audit.ps1 -Mode Fresh -AuditName my_gc_L16_L24 -Sizes 16,24 -Betas 0.24,0.27,0.29,0.32 -Mu -3.53 -ProcessesPerChain 1 -RuntimeLimit 300
.\scripts\run_gc_bkt_audit.ps1 -Mode Resume -AuditName my_gc_L16_L24
.\scripts\run_gc_bkt_audit.ps1 -Mode Production -AuditName my_gc_L16_L24 -ProductionIndex 1
.\scripts\run_gc_bkt_audit.ps1 -Mode Production -AuditName my_gc_L16_L24 -ProductionIndex 2
.\scripts\run_gc_bkt_audit.ps1 -Mode Production -AuditName my_gc_L16_L24 -ProductionIndex 3
~~~

Combine these points with the L = 32 audit:

~~~powershell
py -3.13 .\scripts\analyze_gc_bkt_audit.py .\results\my_gc_L16_L24 .\results\my_gc_L32 --output-dir .\results\my_gc_combined_analysis
~~~

Do not run the L = 16/L = 24 and L = 32 batches simultaneously: each command
above is already designed to occupy all 24 physical cores. Inspect burn-in and
add Resume or Production blocks when the convergence diagnostics require them.

## 10. How to decide whether it converged

A successful executable exit is necessary but not sufficient. Check all of the following.

1. **Start erasure.** Low, middle, and high initial families give compatible late mean N.
2. **Block stationarity.** Production 1, 2, and preferably 3 agree without a sustained drift.
3. **Number fluctuations.** Var(N) is positive and stable. A changing N alone is not enough.
4. **Response check.** A nearby-mu finite difference is compatible with beta Var(N).
5. **Winding stability.** Winding agrees across starts and blocks.
6. **Spatial isotropy.** Wx squared and Wy squared agree within sampling scatter on the square lattice.
7. **Cutoff control.** Selected points should be repeated at larger nmax if the occupation tail may approach 16.
8. **Finite-size control.** A transition claim needs several L values.

The pooled ALPS error assumes that each rank's retained history represents the same stationary distribution. If rank or start means are much more dispersed, use the between-family scatter as the honest uncertainty and run longer.

Large autocorrelation estimates mean that millions of recorded excursions may correspond to far fewer independent number or winding samples. Sample count is not effective sample count.

## 11. Fixed mu is not fixed density

At fixed mu, changing beta generally changes mean N. This is a legitimate grand-canonical thermodynamic path:

~~~text
mu fixed
temperature changed
density measured and allowed to change
~~~

If the scientific goal is a constant-density temperature scan, tune mu separately at every beta:

1. choose two or more nearby mu values;
2. run independent-start production at each;
3. interpolate mean density against mu;
4. start a fresh audit at the interpolated mu;
5. verify the target density and repeat for the next beta.

Never change mu and Resume an old checkpoint.

## 12. Energy labels

The historical HDF5 names need interpretation:

- Kinetic_Energy is the hopping energy.
- Potential_Energy is interaction energy minus mu times N.
- Total_Energy is the expectation of H0 minus mu N.

The physical internal energy without the reservoir term is

$$
\langle H_0\rangle =
\langle H_0-\mu N\rangle+\mu\langle N\rangle .
$$

The summary helper prints both forms. It does not attach an error to the reconstructed H0 because the covariance between Total_Energy and N is not stored separately.

## 13. Independent numerical checks already completed

### L = 16 start-family test

At beta = 0.2790771575, three independently seeded families were started well below, near, and well above equilibrium. Five retained blocks per family at mu = -3.58 gave:

| Quantity | Result |
|---|---:|
| Mean N | 517.943 plus or minus 0.885 |
| Mean density | 2.02321 |
| Var(N) | 7114.00 plus or minus 124.12 |
| Winding sum | 0.77072 plus or minus 0.00581 |

The low/mid/high family means for N were 518.35, 519.23, and 516.25. Their agreement shows loss of the starting-number bias.

At mu = -3.53, three retained blocks per family gave:

| Quantity | Result |
|---|---:|
| Mean N | 612.968 plus or minus 2.829 |
| Mean density | 2.39441 |
| Var(N) | 6949.16 plus or minus 560.11 |
| Winding sum | 1.30067 plus or minus 0.01586 |

The low/mid/high family means for N were 609.50, 610.83, and 618.57. The value 2.39441 is about 1.41 percent below the earlier target filling 2.42864, which demonstrates why fixed mu must not be described as fixed density.

Using the two chemical potentials,

$$
\frac{\Delta\langle N\rangle}{\Delta\mu}=1900.50
$$

while beta times the endpoint-average Var(N) is 1962.35. Their ratio is 0.9685, a 3.15 percent difference. The descriptive uncertainties are 59.28 and 80.05 respectively, and the finite difference spans 0.05 in mu, so exact equality is not expected. This is a strong independent check that mu controls the sampled number distribution quantitatively.

These corrected values pool N and N squared with their measurement counts inside each independent family, then give the three families equal weight. The uncertainties are descriptive standard errors across those family means. Reproduce the calculation with:

~~~powershell
py -3.13 .\scripts\analyze_gc_response.py .\results\grand_canonical_bkt_audit_20260924\pilot_L16_mu_m3p58 .\results\grand_canonical_bkt_audit_20260924\tune_L16_mu_m3p53
~~~

### Fixed-mu finite-size temperature audit

The final audit held mu = -3.53, U/t = 0.152, nmax = 16, and periodic square geometry fixed. Every size-temperature point used independently seeded low, middle, and high starting-number families. Fresh and Resume outputs were discarded. L = 16 and 24 retained five reset blocks per family; L = 32 retained three. The final analysis therefore read 156 HDF5 production blocks.

These archived roots were originally prepared with a legacy manifest that predates the ParameterSHA256 column. The analyzer independently checked the physics parameters, geometry, seed, and initial N stored inside every HDF5 block against that manifest. The current runner is stricter: every newly prepared INI is hashed, and Fresh, Resume, or Production refuses a missing or changed hash before launching any child process. Consequently, reproduce or extend the science with a new audit name rather than continuing the legacy roots below.

| L | beta | T/t | density | Wx squared | Wy squared | winding sum | start SEM |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 0.24 | 4.1667 | 2.21016 | 0.27695 | 0.28239 | 0.55934 | 0.01457 |
| 16 | 0.27 | 3.7037 | 2.33528 | 0.55220 | 0.55042 | 1.10262 | 0.02108 |
| 16 | 0.29 | 3.4483 | 2.45261 | 0.75548 | 0.76607 | 1.52155 | 0.01068 |
| 16 | 0.32 | 3.1250 | 2.60664 | 1.06836 | 1.07701 | 2.14538 | 0.02371 |
| 24 | 0.24 | 4.1667 | 2.11527 | 0.12867 | 0.12890 | 0.25758 | 0.00134 |
| 24 | 0.27 | 3.7037 | 2.32154 | 0.46477 | 0.45893 | 0.92369 | 0.03242 |
| 24 | 0.29 | 3.4483 | 2.41865 | 0.68611 | 0.69868 | 1.38479 | 0.02344 |
| 24 | 0.32 | 3.1250 | 2.56123 | 1.01084 | 1.01356 | 2.02439 | 0.03692 |
| 32 | 0.24 | 4.1667 | 2.12315 | 0.07114 | 0.07133 | 0.14247 | 0.00332 |
| 32 | 0.27 | 3.7037 | 2.27497 | 0.36526 | 0.36460 | 0.72986 | 0.01737 |
| 32 | 0.29 | 3.4483 | 2.39352 | 0.64385 | 0.64023 | 1.28408 | 0.03329 |
| 32 | 0.32 | 3.1250 | 2.56956 | 1.01099 | 1.00212 | 2.01312 | 0.01820 |

Start SEM is the descriptive standard error of the three independent-family means, not a precision confidence interval. The deliberately imposed initial N spans were 307, 691, and 1229 particles for L = 16, 24, and 32. Their retained family-mean spans fell to at most 19.2, 42.6, and 67.7 particles respectively. The x-y winding mismatch was below 2 percent at every point. These are strong start-erasure and square-lattice isotropy checks.

All three sizes show winding far below 4/pi at beta = 0.24 and well above it at beta = 0.32. Linear interpolation through the thermodynamic 4/pi line gives raw finite-size diagnostics of beta = 0.27815, 0.28516, and 0.28961 for L = 16, 24, and 32, corresponding to about 40.13, 39.15, and 38.55 nK for the stated grid mapping. The ordered size trend and resolved normal/low-temperature regimes are finite-temperature BKT-like behavior.

They are not precision transition estimates. At L = 24, beta = 0.29, the family-averaged winding moved from 1.1988 in block 2 to 1.5136 and 1.5147 in blocks 4 and 5 while mean N also rose. Pooling all five blocks gives 1.3848, whereas pooling only blocks 3-5 gives about 1.4859. L = 24, beta = 0.27 and several other points also retain visible block scatter, and the minimum estimated effective N sample count is only about 141. The correct conclusion is that the chains are approaching a common distribution and robustly bracket the crossover, but longer blocks and a finer beta grid are required before a finite-size fit.

Reproduce the complete table with:

~~~powershell
py -3.13 .\scripts\analyze_gc_bkt_audit.py .\results\gc_bkt_L16_L24_mu_m3p53_20260924 .\results\gc_bkt_L32_mu_m3p53_20260924 --output-dir .\results\gc_bkt_combined_mu_m3p53_20260928
~~~

On this desktop, use the py -3.13 launcher for these HDF5 analysis scripts. The plain python command inherits an MSYS2 Python path containing an incomplete h5py namespace and is not reliable in the current shell environment.

### What the early failed runs taught us

Short L = 64 runs begun at one, two, and three particles per site remained in widely separated number sectors and showed no useful winding. Those data were not equilibrium results and are not used.

An L = 8 smoke test changed N, restored checkpoints, and reset statistics, but its mean moved strongly between blocks. It validated mechanics, not physics.

One L = 16 Resume call reached the first expensive DensDens_CorrFun sample on only some MPI ranks. ALPS correctly refused to merge an observable present on only part of the communicator. The checkpoints remained valid, and another Resume completed successfully. The audit runner's smaller technical threshold avoids this edge case.

The known on-site Density_Matrix estimator discrepancy from the exact-diagonalisation regression is unrelated to N, energy, number response, or winding. Density_Matrix is not used for the conclusions in this guide.

## 14. What is enough for a BKT claim

The completed coarse L = 16, 24, and 32 scan establishes that the code reaches:

- a high-temperature regime with winding well below 4/pi;
- a low-temperature regime with winding above 4/pi;
- an intermediate crossover region.

That is evidence of finite-temperature BKT-like behavior. It is **not** enough to quote a precision transition temperature.

For a publishable estimate:

1. run at least three useful sizes, for example L = 16, 24, 32, and ideally larger;
2. use a finer beta grid around the size-dependent crossover;
3. retain several blocks from several independent starts at each point;
4. repeat selected points at a larger nmax;
5. fit the logarithmic finite-size form rather than simply solving Wsum = 4/pi;
6. state whether the path is fixed mu or density-matched;
7. propagate both within-chain and between-chain uncertainty.

## 15. Rebuilding

The packaged binaries are already rebuilt against the verified ALPSCore installation. Parameter edits and new sweeps do not require a rebuild.

Only after source or dependency changes:

~~~powershell
Set-Location '\\alfs1.physics.ox.ac.uk\al\kuo\projects\worm'
.\scripts\build_windows.ps1 -Jobs 32
~~~

The exact build and ALPSCore layout are documented in [WINDOWS_BUILD.md](WINDOWS_BUILD.md).
