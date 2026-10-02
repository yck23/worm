# Start here: running the 64 x 64 Rb-87 Bose-Hubbard calculation

This is the beginner guide for fixed particle number. It starts with what the simulation controls, then explains the worm update, the files, and the commands. If you want to choose chemical potential and let the simulation find the density, use the step-by-step [grand-canonical instructions](GRAND_CANONICAL_RUN_GUIDE.md#8-run-one-editable-parameter-file) instead.

If this file is open as plain text, Markdown symbols will be visible. In VS Code, press Ctrl+Shift+V to open the formatted preview. Only text inside a PowerShell command block is meant to be typed.

## The most important fact

This workflow uses the **canonical ensemble**. You choose the total particle number before the run.

For the supplied 64 x 64 calculation:

| Quantity | Value |
|---|---:|
| Number of lattice sites | 64 x 64 = 4096 |
| Exact particle number N | 9948 |
| Mean filling N / number of sites | 2.4287109375 |
| Ensemble | Canonical: N is fixed |

The line that fixes the number is:

```text
canonical = 9948
```

The chemical potential does **not** find the particle number in this run. The simulation is explicitly told to keep N = 9948.

If you instead want to fix chemical potential and measure a fluctuating particle number, use the separate [GRAND_CANONICAL_RUN_GUIDE.md](GRAND_CANONICAL_RUN_GUIDE.md). That path has passed independent-start, fluctuation-response, checkpoint, and finite-temperature winding tests at L = 16, 24, and 32. Some crossover points retain slow block-to-block number motion, so the scan is a diagnostic rather than a precision transition estimate. A 64 x 64 grand-canonical result still needs its own long convergence test; success at smaller sizes does not automatically validate it.

## 1. What physical system is being simulated?

The model is a homogeneous two-dimensional Bose-Hubbard model on a square grid with periodic boundaries. Periodic means that leaving the right edge returns through the left edge, and similarly in the vertical direction. Topologically, the grid is a torus.

The important physical parameters are:

| Parameter | Meaning | Supplied value |
|---|---|---:|
| Lx, Ly | Number of sites in each direction | 64, 64 |
| t_hop | Hopping energy; used as the energy unit | 1.0 |
| U_on | Repulsion between particles on one site, in units of t | 0.152 |
| V_nn | Nearest-neighbour density interaction | 0.0 |
| beta | Inverse temperature in units of 1/t | 0.279077157505874 |
| canonical | Exact total particle number | 9948 |
| nmax | Largest allowed occupation of one site | 16 |

The lattice is a spatial discretisation of a weakly interacting Rb-87 gas. It is not a literal model of the experimental trap, bilayer, imaging system, or box walls.

## 2. How was N = 9948 obtained?

N = 9948 is a chosen physical input, not a number discovered by the Monte Carlo calculation.

For the present mapping:

- the grid spacing is 0.5 micrometres;
- the hopping scale is t/kB = 11.1630863 nK;
- 40 nK therefore corresponds to T/t = 3.58323844 and beta = 0.279077158;
- a weak-gas critical phase-space-density estimate gives a target filling of about 2.428635525 particles per site.

The integer particle number is then chosen as

```text
N = round(2.428635525 x 64 x 64) = 9948.
```

Rounding makes the actual filling 9948/4096 = 2.4287109375. This is close to the target but is not an experimental equation-of-state fit. A serious result should also test nearby densities.

When a fresh canonical run is created, the code constructs N = 9948 directly:

1. It puts 2 particles on every one of the 4096 sites. This accounts for 8192 particles.
2. It randomly chooses 1756 sites and puts one additional particle on each.
3. The initial sum is therefore 8192 + 1756 = 9948.
4. Worm updates redistribute those particles, but every retained physical configuration still has total N = 9948.

The output observable **Number_of_particles** should consequently be exactly 9948. It is a useful consistency check, but it is not an estimate with a physical number fluctuation in this canonical calculation.

## 3. What does the chemical potential do here?

The Hamiltonian still contains the term -mu N, and the parameter file contains:

```text
mu = -3.47
```

In a grand-canonical calculation, mu controls the statistical preference for different particle-number sectors, and the simulation measures a fluctuating mean N.

This calculation is different. Since N is fixed, -mu N has the same value for every closed physical configuration. It therefore cancels from relative probabilities among those configurations. Changing mu cannot tune the density in this exact-N run.

There are two qualifications:

1. The reported **Potential_Energy** and **Total_Energy** include the -mu N term. To obtain the interaction-only potential energy or the hopping-plus-interaction total energy, add mu N to the reported value:

```text
energy without the -mu*N term = reported energy + mu*N
                                  = reported energy - 34519.56
```

2. The code also uses mu while proposing temporary open-worm moves. It can therefore affect acceptance and mixing efficiency even though it does not set the canonical density. Keep the audited value unless you deliberately retune and revalidate the sampler.

The energy obtained by adding mu N is the Bose-Hubbard lattice internal energy. For the continuum finite-difference energy convention, add a further 4*t*N to restore the square-grid kinetic-energy zero. This distinction matters when comparing simulation energies with a continuum model or the thesis; it does not alter the run.

If the scientific question instead requires chemical potential to determine the mean particle number, use a grand-canonical executable and measure **Number_of_particles**. That is not the workflow in this guide. Early short 64 x 64 tests started far apart and remained trapped in different number sectors. Later 16 x 16 and 32 x 32 audits showed that the same machinery can converge when the starts, burn-in, and retained blocks are controlled. The exact-N workflow remains the independently replicated 64 x 64 reference.

The executable, parameter template, commands, and required start-independence test are provided in [GRAND_CANONICAL_RUN_GUIDE.md](GRAND_CANONICAL_RUN_GUIDE.md).

## 4. What is a worm update?

The finite-temperature quantum problem is represented by particle occupations evolving around imaginary time from 0 to beta.

- A particle traces a **worldline** through space and imaginary time.
- A hopping event moves one particle between neighbouring sites. In a worldline picture this event is a corner, often called a **kink**.
- A physical contribution to the partition function has closed worldlines: there are no loose ends.

Changing only closed lines locally makes it difficult to rearrange long exchange cycles and spatial winding. The worm algorithm temporarily opens a line by creating two endpoints. One endpoint then moves through space and imaginary time, inserting or deleting hopping events. Finally, the two endpoints join and the worldlines are closed again.

The open configuration is update machinery, not a physical state being reported as an ordinary energy or particle-number sample. In this code one source-level **sweep** is one update call that starts and ends closed. If insertion succeeds, that call contains the complete variable-length worm excursion; a rejected insertion simply remains closed. A sweep is not one attempted move per lattice site and is not physical time.

The setting

```text
canonical_window = 0.5
```

allows the imaginary-time-averaged number to move temporarily while the worm is open, but by less than half a particle from the target. Because a closed sector has an integer particle number, a window strictly below 1 prevents closure at N + 1 or N - 1. Thus **canonical_window = 0.5 does not mean that measured N fluctuates by plus or minus 0.5**. Closed measured configurations remain exactly at N = 9948.

The value 0.5 was selected because a much narrower value trapped the worm and gave poor winding-sector mixing. Do not set it to 1 or larger; the executable rejects that choice.

## 5. What does each file do?

There are five distinct objects:

| Object | Purpose |
|---|---|
| Parameter file | Your editable description of a new run |
| Runner script | Chooses Fresh, Resume, or Production and starts MPI |
| Executable | The compiled Monte Carlo program |
| Checkpoint files | Complete saved Markov-chain states, one per MPI rank |
| Result file | HDF5 measurements combined from all MPI ranks |

The normal flow is:

```text
edited INI file --Fresh--> checkpoints + result
checkpoints      --Resume--> continued checkpoints + burn-in result
checkpoints  --Production--> continued checkpoints + reset production result
```

The starting parameter file is:

```text
parameter_files\Rb87_BKT_64.ini
```

The dedicated runner is:

```text
scripts\run_rb87_bkt64.ps1
```

The commands in this guide use the runner's default exact-N MPI executable. The same runner can select its other executable with -Ensemble GrandCanonical, as explained in the grand-canonical guide. It does not override beta, U, N, lattice size, cutoff, or worm settings.

### A crucial checkpoint rule

**Fresh** reads the complete INI file and saves its parameters into the checkpoints.

**Resume** and **Production** restore the existing checkpoint, including its run-length and measurement settings. Editing the INI after Fresh does not change those settings or the physical parameters inside that chain. The runner reads the INI during continuation to locate the named output/checkpoint files and validate the requested ensemble. It does not use it to retune a saved simulation.

Therefore:

- finish editing before Fresh;
- do not edit that run's INI while continuing it;
- for any changed temperature, size, N, U, cutoff, worm window, seed, or MPI rank count, make a new INI copy, choose a new output directory, and start Fresh.

## 6. Build the program once

Open PowerShell and move to the repository:

```powershell
Set-Location '\\alfs1.physics.ox.ac.uk\al\kuo\projects\worm'
```

Build with all 32 logical processors:

```powershell
.\scripts\build_windows.ps1 -Jobs 32
```

The audited build uses the pinned and verified ALPSCore installation packaged by the build script. The executables and required DLLs are placed in:

```text
dist\windows-square
```

Build again only after changing source code, compiler/dependency installation, or CMake options. Editing an INI file never requires a rebuild.

## 7. Make a parameter file for one run

Keep the supplied file as a template. Make a named copy for each independent chain family:

```powershell
Copy-Item .\parameter_files\Rb87_BKT_64.ini .\parameter_files\Rb87_L64_T40_myrun01.ini
notepad .\parameter_files\Rb87_L64_T40_myrun01.ini
```

For this first run, leave the audited physics and worm settings unchanged. Change **seed** to a new positive integer, for example:

```text
seed = 140001
```

The outputfile and checkpoint lines are only filenames inside the output directory. They can keep their supplied values because every chain family will use a separate directory.

The parameter groups are:

| Group | Parameters | What to do initially |
|---|---|---|
| Physical state | Lx, Ly, beta, t_hop, U_on, V_nn, canonical | Change only when defining a new physical point |
| Hilbert-space approximation | nmax | Keep 16 for this regime |
| Sampler controls | canonical_window, C_worm, move probabilities, E_off | Keep the audited values |
| Run length | runtimelimit, sweeps, thermalization, Nmeasure, Nmeasure2 | Keep the audited values for the first run |
| Independent-chain identity | seed | Use a different value for a genuinely new replica family |
| File names | outputfile, checkpoint | Usually leave unchanged when output directories are unique |

## 8. Run one complete 64 x 64 chain family

Use 24 MPI processes because this desktop has 24 physical cores. Each MPI rank simulates a complete 64 x 64 lattice with a different random-number stream. The lattice is **not** divided into 24 pieces.

Burn-in is the initial part of a Markov chain that is discarded while it forgets the deliberately simple starting occupation pattern. The INI requests 50,000 sweeps before measurements, but a 900-second block may end before every rank reaches that count. Resume continues both the saved worldlines and their counters.

Set the two paths once in the same PowerShell window:

```powershell
$parameterFile = '.\parameter_files\Rb87_L64_T40_myrun01.ini'
$outputDirectory = '.\results\Rb87_L64_T40_myrun01'
```

### Step 1: create a fresh chain

```powershell
.\scripts\run_rb87_bkt64.ps1 -Mode Fresh -Processes 24 -ParameterFile $parameterFile -OutputDirectory $outputDirectory
```

Fresh refuses to overwrite an existing result. If the directory already belongs to another run, choose a new name rather than deleting it.

### Step 2: continue burn-in three times

Run these in order:

```powershell
.\scripts\run_rb87_bkt64.ps1 -Mode Resume -Processes 24 -ParameterFile $parameterFile -OutputDirectory $outputDirectory
.\scripts\run_rb87_bkt64.ps1 -Mode Resume -Processes 24 -ParameterFile $parameterFile -OutputDirectory $outputDirectory
.\scripts\run_rb87_bkt64.ps1 -Mode Resume -Processes 24 -ParameterFile $parameterFile -OutputDirectory $outputDirectory
```

Fresh plus these Resume blocks are burn-in. Do not treat their measurements as final data. Three Resume blocks are the safer starting prescription found by the 64 x 64 audit; convergence still has to be checked from later blocks.

### Step 3: collect at least three production blocks

Production continues the same equilibrated chain state but clears the old measurement accumulators first. After each block, make a simply named copy:

```powershell
.\scripts\run_rb87_bkt64.ps1 -Mode Production -Processes 24 -ParameterFile $parameterFile -OutputDirectory $outputDirectory
Copy-Item (Join-Path $outputDirectory 'Rb87_BKT_64.out.h5') (Join-Path $outputDirectory 'production_01.out.h5')

.\scripts\run_rb87_bkt64.ps1 -Mode Production -Processes 24 -ParameterFile $parameterFile -OutputDirectory $outputDirectory
Copy-Item (Join-Path $outputDirectory 'Rb87_BKT_64.out.h5') (Join-Path $outputDirectory 'production_02.out.h5')

.\scripts\run_rb87_bkt64.ps1 -Mode Production -Processes 24 -ParameterFile $parameterFile -OutputDirectory $outputDirectory
Copy-Item (Join-Path $outputDirectory 'Rb87_BKT_64.out.h5') (Join-Path $outputDirectory 'production_03.out.h5')
```

Each invocation has a 900-second wall-time limit. The runner also archives the previous result automatically before replacing it.

Always use the same MPI process count when resuming a checkpoint family. There is one checkpoint per rank.

## 9. Read the basic results

The result is an HDF5 file rather than a text table. For the third production block:

```powershell
$h5dump = 'C:\msys64\ucrt64\bin\h5dump.exe'
$result = Join-Path $outputDirectory 'production_03.out.h5'

& $h5dump -d '/simulation/results/Number_of_particles/mean/value' $result
& $h5dump -d '/simulation/results/Total_Energy/mean/value' $result
& $h5dump -d '/simulation/results/Winding_number_squared/mean/value' $result
& $h5dump -d '/simulation/results/Winding_number_squared/mean/error' $result
```

Interpret these as follows:

- **Number_of_particles** should equal 9948 exactly.
- **Total_Energy** includes the -mu N offset explained above.
- **Winding_number_squared** returns two entries, one for x and one for y.
- Add the x and y winding entries to obtain Wsum.
- Similar x and y values are an isotropy check.

For this code's normalization, the infinite-system BKT jump corresponds to

```text
Wsum = <Wx^2> + <Wy^2> = 4/pi = 1.273239545.
```

A single 64 x 64 value near 4/pi is not a transition-temperature measurement. Finite systems have smooth crossovers and logarithmic size corrections.

## 10. How do I know the run has converged?

Exact N alone does not demonstrate convergence; N is constrained by construction. Check the quantities that still have to explore configuration space:

1. Compare Total_Energy and both winding components across production_01, production_02, and production_03.
2. Discard an early production block if the later blocks have clearly moved away from it.
3. Continue making Production blocks until the late blocks agree within their uncertainties and show no chronological drift.
4. Repeat the whole Fresh-to-Production sequence with another seed and another output directory.
5. Require independent seed families to agree. Twenty-four MPI ranks improve statistics, but they all begin as one simultaneously launched family and do not replace a deliberate second seed family.

The audit found that winding equilibrates more slowly than energy. A small internal error bar does not rescue a result whose block means continue to drift.

Do not currently use **Density_Matrix[0]** to normalize correlations. An exact-identity regression found that this on-site estimator is biased. The audited energy, exact particle-number, and winding path is not affected by that defect.

## 11. Change temperature or size correctly

For this Rb-87 mapping, convert physical temperature to the INI value with

```text
beta = 11.1630863 / temperature_in_nK.
```

Examples:

| Temperature | beta |
|---:|---:|
| 32 nK | 0.348846446875 |
| 40 nK | 0.279077157500 |
| 49.8 nK | 0.224158359438 |

To keep approximately the same density when changing size, use

```text
canonical = round(2.428635525 x Lx x Ly).
```

For example, use N = 2487 at 32 x 32 and N = 9948 at 64 x 64.

For every changed temperature or size:

1. copy the template to a new INI filename;
2. edit beta, Lx, Ly, canonical, and seed;
3. choose a new output directory;
4. start with Fresh;
5. repeat burn-in, production blocks, and convergence checks.

Do not change mu in an attempt to change density in this canonical workflow.

## 12. Make a temperature-and-size scan

The sweep helper accepts dimensionless T/t, not beta and not nK. This dry run creates the point-specific INI files and a manifest without launching simulations:

```powershell
.\scripts\run_temperature_sweep.ps1 -ParameterFile .\parameter_files\Rb87_BKT_64.ini -Temperatures 2.866590756,3.583238445,4.461131865 -Sizes 32,64 -Canonical -Density 2.428635525 -Sweeps 4800000 -Thermalization 50000 -RuntimeLimit 900 -Processes 24 -SweepName rb87_bkt_check -DryRun
```

Here the temperatures correspond to approximately 32, 40, and 49.8 nK. The helper holds the filling fixed, calculates a separate exact integer N for each size, and writes beta = 1/(T/t) into every generated INI. It does not use mu to obtain that density.

Remove **-DryRun** to perform one fresh block at each point. Those one-block jobs are reconnaissance only. A BKT result requires continuing each important point through burn-in and independent production blocks, then comparing several sizes with finite-size scaling.

## 13. Five common questions

### Do command-line parameters overwrite my INI?

The dedicated Rb runner does not override its physical parameters. Its command line selects the INI path, output directory, MPI process count, and run mode. The separate sweep helper deliberately generates a new INI for each point with its requested size, temperature, run length, seed, and exact N.

### Can I edit the INI and then Resume?

No. Resume restores the parameters saved in the checkpoint. Make a new parameter file and start Fresh in a new directory.

### Does canonical_window = 0.5 mean N fluctuates?

No. It permits a temporary fractional imaginary-time-averaged displacement while a worm is open. Closed measured states still have exactly N = 9948.

### Does a Production command create an independent run?

No. It keeps the same worldline and random-number-generator state and only resets accumulated measurements. A genuinely independent family needs a new seed, a new directory, and Fresh.

### Does one point below 4/pi prove a BKT transition?

No. Recovering the transition requires a temperature scan, several lattice sizes, converged independent runs, and BKT finite-size scaling.

For the derivation, code map, audit evidence, and full BKT analysis workflow, continue with [COMPLETE_BKT_GUIDE.md](COMPLETE_BKT_GUIDE.md). The completed numerical checks are in [RB87_BKT_AUDIT.md](RB87_BKT_AUDIT.md), and build details are in [WINDOWS_BUILD.md](WINDOWS_BUILD.md).
