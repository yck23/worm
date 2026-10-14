# Your first run: a 4 x 4 Bose-Hubbard example

Start here. Steps 1-4 run a small example and print its results. On this desktop the checked run took about 39 seconds. The program is already built, and the example parameter file is already supplied.

## 1. Open a PowerShell terminal in VS Code

Use the menu Terminal > New Terminal. Choose PowerShell from the terminal's profile menu if needed.

Paste each command below into that terminal and press Enter. Wait for it to finish before pasting the next command. Copy only the command inside each box.

## 2. Go to the project folder

~~~powershell
Set-Location '\\alfs1.physics.ox.ac.uk\al\kuo\projects\worm'
~~~

Run all the following commands from this folder. You can leave VS Code's terminal here throughout the example.

## 3. Run the supplied example

~~~powershell
.\scripts\run_rb87_bkt64.ps1 -Ensemble GrandCanonical -Mode Fresh -Processes 1 -ParameterFile .\parameter_files\basic_bose_square.ini -OutputDirectory .\results\first_example
~~~

Wait about 30-60 seconds. When it finishes, look for:

~~~text
# Finished successfully.
Completed Fresh.
~~~

The program prints many other lines; you can ignore them for this first exercise.

The script's name contains rb87_bkt64, but it runs the size written in your parameter file. This example is 4 x 4.

If you have already run this example, Fresh will say the output exists. To continue it, use the Resume command in Section 6. To start it again from scratch, change first_example to another folder name, such as first_example_again, in both the run and read commands.

## 4. Print the results

~~~powershell
py -3.13 .\scripts\summarize_bkt_result.py .\results\first_example\basic_bose_square.out.h5
~~~

You should see lines similar to these. The numerical values and sample count may differ:

~~~text
ensemble: grand canonical (mu fixed, N measured)
lattice: 4 x 4 = 16 sites
beta: 1
mu: 0.7
samples: 372576
<N>: 22.96 ...
<N>/sites: 1.435 ...
Var(N): 4.27 ...
~~~

Mean N is the average number of particles in the entire lattice. Dividing it by 16 gives particles per site. Var(N) measures how much the particle number fluctuates.

You have now run the simulation and read its output. You can stop here for your first test. This short exercise checks the execution procedure; a physics study needs longer sampling and convergence checks.

## What did those commands use?

Your input is [parameter_files/basic_bose_square.ini](parameter_files/basic_bose_square.ini). It is a plain text file containing the physical parameters and run settings.

The wrapper script starts the compiled executable in dist/windows-square and puts this run's files in results/first_example. You do not need to move into the executable's folder or activate an environment.

Inside results/first_example:

- basic_bose_square.out.h5 contains measurements; use the Python command above to read it.
- basic_bose_square.clone.h5 contains the saved simulation state, used to continue later.
- tmp and archive are folders managed by the runner.

The main lines in the supplied INI are:

~~~text
Lx = 4
Ly = 4
t_hop = 1.0
U_on = 4.0
mu = 0.7
beta = 1.0
canonical = -1
initial_particle_number = 16
runtimelimit = 30
~~~

Here t_hop sets the energy unit, U_on sets the on-site repulsion, mu controls particle number, and beta is inverse temperature. canonical = -1 means that N is allowed to fluctuate. initial_particle_number = 16 only chooses the starting state; the measured mean N can differ from 16.

These are simple dimensionless practice parameters. For the Rb-87 calculation, use the separate Rb-87 template after completing this exercise.

## 5. Optional: edit one parameter and run another example

Make your own copy. Use these two commands, one at a time:

~~~powershell
Copy-Item .\parameter_files\basic_bose_square.ini .\parameter_files\my_bose_square.ini
~~~

~~~powershell
notepad .\parameter_files\my_bose_square.ini
~~~

In the editor, find the existing beta line and replace it with:

~~~text
beta = 0.5
~~~

Save and close the editor. This doubles temperature because beta is inverse temperature. Leave the other lines unchanged.

Run the edited file in a new output folder:

~~~powershell
.\scripts\run_rb87_bkt64.ps1 -Ensemble GrandCanonical -Mode Fresh -Processes 1 -ParameterFile .\parameter_files\my_bose_square.ini -OutputDirectory .\results\second_example
~~~

Read that run:

~~~powershell
py -3.13 .\scripts\summarize_bkt_result.py .\results\second_example\basic_bose_square.out.h5
~~~

The result should now say beta: 0.5. Its filename remains basic_bose_square.out.h5 because that is the outputfile line in your INI. The new folder keeps it separate from your first result.

For later experiments, choose a new INI filename and output folder. Edit and save the INI before starting Fresh.

## 6. Optional: continue the original run

Return to the project folder with the command in Section 2 if you have reopened the terminal.

This continues the first example for another short invocation:

~~~powershell
.\scripts\run_rb87_bkt64.ps1 -Ensemble GrandCanonical -Mode Resume -Processes 1 -ParameterFile .\parameter_files\basic_bose_square.ini -OutputDirectory .\results\first_example
~~~

Read it again using the command in Section 4. Resume retains the previous measurements and adds more.

To practise collecting a new measurement block while continuing the saved state:

~~~powershell
.\scripts\run_rb87_bkt64.ps1 -Ensemble GrandCanonical -Mode Production -Processes 1 -ParameterFile .\parameter_files\basic_bose_square.ini -OutputDirectory .\results\first_example
~~~

Read it again with Section 4's command. Production clears the measurement counters and continues the same saved state. It does not make an independent run.

Keep Processes 1 for these restarts. Editing the original INI does not change a saved run's parameters; use a new output folder and Fresh when you want changed settings.

## When you are ready for larger lattices

For chemical-potential-controlled Rb-87 runs, follow [GRAND_CANONICAL_RUN_GUIDE.md, Section 8](GRAND_CANONICAL_RUN_GUIDE.md#8-run-one-editable-parameter-file). Its longer calls and convergence checks are the next stage.

For the original fixed-particle-number 64 x 64 reference, use [COMPLETE_BKT_GUIDE.md, Section 9](COMPLETE_BKT_GUIDE.md#9-commands-in-order-one-trustworthy-64-x-64-point). The detailed physics and algorithm explanation are in that same guide.

If an executable is missing on another installation, the build procedure is in [WINDOWS_BUILD.md](WINDOWS_BUILD.md). The existing verified package on this desktop is sufficient for this exercise.
