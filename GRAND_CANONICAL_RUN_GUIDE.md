# Grand-canonical Bose-Hubbard and BKT guide

This is the guide to use when **chemical potential is fixed and particle number is measured**.

Independent audit completed: 2 October 2026, for commit bd0a6d9 against the preceding PR merge fb73840. Source comparison, fresh-run comparisons, small-system finite-temperature checks, checkpoint checks, the verified ALPSCore package, and an independent recalculation of the archived results are complete. The supported usage is grand-canonical number, energy, and winding measurements with explicit convergence checks. The L = 16, 24, and 32 scan shows BKT-like regimes, but residual drift near the crossover prevents a precision transition-temperature claim. The inherited on-site Density_Matrix estimator issue remains; do not use it to measure density or normalize correlations.

For the commands in order, go straight to [Section 8: run one editable parameter file](#8-run-one-editable-parameter-file). It starts with a 32 x 32 example. For an automated temperature scan with independent starting states, continue to Section 9. The separately audited 64 x 64 reference is canonical; its success does not establish convergence of a 64 x 64 grand-canonical run.

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
- **nmax** is the largest allowed occupation on one site. This truncates the model's state space, so it is not just a speed setting; check selected results at a larger value.
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

For this mapping, a 0.5 micrometre grid and the Rb-87 mass give t/kB = 11.1630863 nK, so beta = 11.1630863 / temperature_in_nK. The values U/t = 0.152 and beta about 0.279077 therefore represent the chosen weak-coupling, 40 nK starting point. The starting number 9948 on a 64 x 64 grid comes from an approximate critical-density formula, not a measured particle number in the thesis. It is only an initial guess in this ensemble.

This is a homogeneous periodic-grid model inspired by the thesis, not a reproduction of its trap or bilayer. The quoted parameters occur in Appendix D.9; D.8 describes its classical-field Metropolis method, which differs from this quantum worm algorithm. At 40 nK the thermal wavelength is only about 1.87 grid spacings, and the 1 kHz axial level spacing corresponds to about 48 nK. Quantitative comparison with the experiment therefore needs grid-spacing and axial-excitation checks as well as Monte Carlo convergence.

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

The value thermalization = 1000 is a technical accumulator threshold, **not a claim that 1000 worm excursions equilibrate the physics**. The workflow discards complete Fresh and Resume invocations as explicit burn-in. The modest threshold and Nmeasure2 reduce the risk that only some MPI ranks acquire sparse spatial observables before ALPS tries to merge them; every rank still needs enough completed sweeps.
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

Even run-length and measurement settings are restored from the checkpoint. Editing runtimelimit, sweeps, thermalization, or C_worm in the original INI does not retune a resumed chain. To collect more data with its existing settings, repeat Resume or Production as appropriate. To use changed settings through this runner, start Fresh in another directory.

The latest commit leaves the worm move code and the primary energy, number, and winding formulae unchanged from the preceding PR. Default fresh grand-canonical runs matched in three seeded comparisons. New initial-number controls change the starting state, and checkpoint/validation changes affect restart behaviour. The canonical-window setting can substantially change mixing in canonical runs; it has no effect in this grand-canonical executable. Do not expect identical finite-run results after changing starting states or sampling settings.

## 8. Run one editable parameter file

### 8.1 Open the project; build only if needed

Open PowerShell and enter:

~~~powershell
Set-Location '\\alfs1.physics.ox.ac.uk\al\kuo\projects\worm'
~~~

The verified executable package already exists on this desktop. You do not need to activate a Python or ALPSCore environment to run it: the runner sets the DLL search path. If the package is missing, or you have changed the C++ source or dependency, build it once:

~~~powershell
.\scripts\build_windows.ps1 -Jobs 32
~~~

ALPSCore's source is in build-repro/alps-src and its installed libraries are in build-repro/alps-install. The runnable package is in dist/windows-square. These are different parts of the same build, not separate environments you must activate. See [WINDOWS_BUILD.md](WINDOWS_BUILD.md) for compiler and MPI prerequisites on another computer.

### 8.2 Make your own parameter file

Run this once for a new experiment. If this filename already exists, choose another name so you preserve it:

~~~powershell
$parameterFile = '.\parameter_files\my_gc_L32_T40.ini'
if (Test-Path -LiteralPath $parameterFile) { throw 'Choose a new parameter filename.' }
Copy-Item .\parameter_files\Rb87_BKT_64_grand_canonical.ini $parameterFile
notepad $parameterFile
~~~

In Notepad, replace the existing values of these lines, then save and close. Do not append duplicate keys; leave the other template lines unchanged:

~~~text
Lx = 32
Ly = 32
mu = -3.53
beta = 0.279077157505874
canonical = -1
initial_particle_number = 2487
seed = 140001
outputfile = "my_gc.out.h5"
checkpoint = "my_gc.clone.h5"
~~~

This keeps t = 1, U = 0.152, nmax = 16 and periodic x/y boundaries from the template. You control chemical potential and temperature. The number 2487 is a starting guess, not a target that the simulation holds fixed. The measured mean N need not equal it.

### 8.3 Choose where results go and start

Run the following in the same PowerShell window. These variables are only filename and process-count shortcuts; they do not override physics in the INI:

~~~powershell
$parameterFile = '.\parameter_files\my_gc_L32_T40.ini'
$outputDirectory = '.\results\my_gc_L32_T40'
$processes = 2
$resultFile = Join-Path $outputDirectory 'my_gc.out.h5'
.\scripts\run_rb87_bkt64.ps1 -Ensemble GrandCanonical -Mode Fresh -Processes $processes -ParameterFile $parameterFile -OutputDirectory $outputDirectory
~~~

Keep -Ensemble GrandCanonical on every call: the runner otherwise defaults to canonical mode. Its historical filename contains bkt64, but Lx and Ly in your INI determine the actual size.

Two MPI processes run two separately sampled copies and combine their measurements; they do not split the lattice into two pieces. You can choose up to 24 processes before Fresh to use this desktop's 24 physical cores for one family. Keep that count unchanged for all its restarts. More cores improve sampling throughput but do not make an unequilibrated chain equilibrated.

Each call has the template's 900-second wall-time limit, checked after complete updates, so allow some overrun. The sweep limit can also end a call. Wait for completion before issuing the next command in that same result directory. A source sweep is a worm update, not a second of simulated atomic motion.

### 8.4 Continue the initial, discarded part of the run

Run Resume, inspect the result, and repeat as needed:

~~~powershell
.\scripts\run_rb87_bkt64.ps1 -Ensemble GrandCanonical -Mode Resume -Processes $processes -ParameterFile $parameterFile -OutputDirectory $outputDirectory
py -3.13 .\scripts\summarize_bkt_result.py $resultFile
~~~

Use py -3.13 on this desktop; the plain python command currently finds an unsuitable MSYS2 environment. The analysis needs h5py and numpy, already available in the audited Python installation.

These Fresh/Resume measurements are burn-in: they help you judge how the state is changing but are not retained as equilibrium data. Resume reports accumulated means, which can conceal recent drift. Compare separate later blocks and independently started families too. If a call ends before any measurements exist, the summary cannot read a results group; continue with Resume rather than claiming a zero result.

There is no guaranteed number of Resume calls. The line thermalization = 1000 only enables accumulation after a technical threshold; it does not establish equilibrium. In particular, movement of N alone does not prove that its full distribution has been explored.

### 8.5 Save separate measurement blocks

After burn-in, run this block to continue the same state but start new measurement accumulators:

~~~powershell
$blockNumber = 1
$savedBlock = Join-Path $outputDirectory ('production_{0:D2}.out.h5' -f $blockNumber)
if (Test-Path -LiteralPath $savedBlock) { throw 'Choose an unused block number.' }
.\scripts\run_rb87_bkt64.ps1 -Ensemble GrandCanonical -Mode Production -Processes $processes -ParameterFile $parameterFile -OutputDirectory $outputDirectory
Copy-Item -LiteralPath $resultFile -Destination $savedBlock
py -3.13 .\scripts\summarize_bkt_result.py $savedBlock
~~~

Run the same block with blockNumber = 2, then 3, and continue with unused numbers. Only save a block after the runner reports "Completed Production"; if it fails, do not relabel the previous result as new data. The copy preserves each block because the next invocation replaces my_gc.out.h5. The runner also archives the previous output, but those archive files can contain discarded or overlapping histories; do not pool them indiscriminately.

Read the summary as follows:

- Mean N is the measured mean number, and mean N divided by the number of sites is the mean filling.
- Var(N) describes physical number fluctuations; the error on mean N describes uncertainty in its estimated mean. They are not the same thing.
- The winding sum measures the response to a phase twist. Its comparison with 4/pi is a finite-size diagnostic, not by itself a transition measurement.
- Compare N, energy, and winding between late blocks. A systematic drift means more equilibration or sampling is needed, even when an individual error bar looks small.

Three blocks are a useful first comparison, not a convergence guarantee. Production does not create a new independent chain or discard another automatic warm-up interval. Consecutive blocks can remain correlated. Retain only a stationary portion, with the discarded prefix documented, and require agreement between independent starts.

### 8.6 Restart later, change the physics, or repeat independently

| What you want | What to do |
|---|---|
| Continue burn-in | Same files, directory and process count; use Resume |
| Collect another separately accumulated block | Same setup; use Production and a new saved block number |
| Change mu, beta, size, seed or other INI settings | Edit a new copy; choose a new output directory; use Fresh |
| Test independence from the starting state | Same physical parameters, different seed and initial_particle_number; new directory and Fresh |

For example, at 32 x 32 repeat with starting numbers around 1843, 2487 and 3072 and different seeds. Each family must lose memory of that guess and produce compatible late means. Section 9 automates this idea over a temperature grid.

After closing PowerShell, reopen it, return to the project, and re-enter the four variable assignments in Section 8.3. Then use Resume or Production, not Fresh, for that existing calculation. Do not recreate the INI or delete the checkpoints. The files my_gc.clone.h5 and my_gc.clone.h5.1 are the two rank-specific checkpoints, not the numbered measurement blocks.

To change temperature while keeping mu fixed, copy the INI, calculate beta = 11.1630863 / temperature_in_nK, and enter the resulting number in its beta line before a new Fresh run. Enter a number, not that arithmetic expression. Expect the resulting density to change. For a fixed-density temperature scan instead, you must tune mu separately at every temperature.

## 9. Recommended independent-start temperature audit

This helper is an alternative to manually copying an INI for every point. Its command-line choices generate the parameter files and a manifest, then run the batch. If you want to edit every INI yourself, repeat Section 8 for each size, temperature and independent start instead.

Do not edit the helper's generated INIs: their hashes are recorded in the manifest and changed files are rejected. To change that batch's settings, create a new AuditName. Resume and Production reuse the original manifest, not newly supplied size or temperature choices.

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

Temperature decreases when beta increases. In the table's order, the scan runs
from the hotter expected normal side to the colder expected superfluid side.

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

The Bose-Hubbard hopping-plus-interaction energy without the reservoir term is

$$
\langle H_0\rangle =
\langle H_0-\mu N\rangle+\mu\langle N\rangle .
$$

The summary helper prints both forms. It does not attach an error to the reconstructed H0 because the covariance between Total_Energy and N is not stored separately.

The helper's label "physical internal energy" refers to that lattice H0. If comparing with the continuum finite-difference Hamiltonian, restore its kinetic-energy zero as well:

~~~text
continuum internal energy = Total_Energy + (mu_BH + 4*t)*mean_N
~~~

The extra 4*t per particle comes from the square-grid Laplacian. It changes the energy convention, not the sampled distribution. The winding-sum error printed by the helper is also only a component-quadrature estimate; it omits the covariance between the x and y measurements.

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

The known on-site Density_Matrix estimator discrepancy was confirmed again in the independent finite-temperature check. It is inherited from the preceding PR and is not used for N, energy, number response, or winding. Read density from Number_of_particles divided by the site count, not from Density_Matrix[0]. No claim of corrected correlation-function normalization is made here.

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
