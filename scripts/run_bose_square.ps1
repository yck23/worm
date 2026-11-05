<#
.SYNOPSIS
    Run the packaged 2D square-lattice Bose-Hubbard executable.

.EXAMPLE
    .\scripts\run_bose_square.ps1

.EXAMPLE
    .\scripts\run_bose_square.ps1 -Lx 4 -Ly 4 -Sweeps 200 -Thermalization 20 -RuntimeLimit 5

.EXAMPLE
    .\scripts\run_bose_square.ps1 -Processes 4 -ParameterFile .\my_parameters.ini

.EXAMPLE
    .\scripts\run_bose_square.ps1 -Lx 32 -Ly 32 -PotentialFile .\potentials\speckle.dat
    # -PotentialFile automatically selects the programmable-potential package.

.EXAMPLE
    .\scripts\run_bose_square.ps1 -Canonical -ProgrammablePotential -DisablePotential -Lx 32 -Ly 32 -ParameterFile .\parameter_files\BoseHubbard_speckle_canonical.ini
    # A fresh fixed-N run with zero added site potential, keeping the INI map path.
#>
[CmdletBinding()]
param(
    [ValidateRange(3, 10000)] [int]$Lx = 8,
    [ValidateRange(3, 10000)] [int]$Ly = 8,
    [ValidateRange(1, 32)] [int]$Processes = 1,
    [ValidateRange(1, 2147483647)] [int]$RuntimeLimit = 900,
    [ValidateRange(1, 2147483647)] [int]$Sweeps = 100000,
    [ValidateRange(0, 2147483647)] [int]$Thermalization = 100,
    [string]$ParameterFile,
    [string]$JobName,
    [string]$OutputDirectory,
    [switch]$Canonical,
    [switch]$ProgrammablePotential,
    [switch]$DisablePotential,
    [string]$PotentialFile
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
# Supplying a potential file also selects the build that can load it.
if ($PotentialFile) { $ProgrammablePotential = $true }
$packageName = if ($ProgrammablePotential) { 'windows-square-programmable' } else { 'windows-square' }
$packageDir = Join-Path $repoRoot "dist\$packageName"
$serialExe = Join-Path $packageDir 'qmc_worm.exe'
$mpiExe = Join-Path $packageDir 'qmc_worm_mpi.exe'
if ($Canonical) {
    $serialExe = Join-Path $packageDir 'qmc_worm_canonical.exe'
    $mpiExe = Join-Path $packageDir 'qmc_worm_canonical_mpi.exe'
}
$mpiExec = 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe'

if (-not $ParameterFile) { $ParameterFile = Join-Path $packageDir 'BoseHubbard.ini' }
if (-not $OutputDirectory) { $OutputDirectory = Join-Path $repoRoot 'results' }
if (-not $JobName) { $JobName = 'bose_square_{0}x{1}_{2}' -f $Lx, $Ly, (Get-Date -Format 'yyyyMMdd_HHmmss') }
$ParameterFile = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($ParameterFile)
$OutputDirectory = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($OutputDirectory)
if ($PotentialFile) {
    $PotentialFile = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($PotentialFile)
    # The solver checks file existence only when enabled, honoring either
    # the INI flag or -DisablePotential without requiring an unused file.
}

foreach ($required in @($serialExe, $ParameterFile)) {
    if (-not (Test-Path $required)) { throw "Required file not found: $required. Run scripts\build_windows.ps1 with the matching -ProgrammablePotential option first." }
}
if ($Processes -gt 1 -and -not (Test-Path $mpiExec)) { throw "Microsoft MPI launcher not found: $mpiExec" }
if ($Processes -gt 1 -and -not (Test-Path $mpiExe)) { throw "MPI executable not found: $mpiExe" }

New-Item -ItemType Directory -Force $OutputDirectory | Out-Null
$tmpDir = Join-Path $OutputDirectory 'tmp'
New-Item -ItemType Directory -Force $tmpDir | Out-Null
$env:TMPDIR = ($tmpDir -replace '\\','/') + '/'
$env:Path = "$packageDir;$env:Path"

$simulationArgs = @(
    ([IO.Path]::GetFullPath($ParameterFile)),
    "--Lx=$Lx", "--Ly=$Ly", '--Lz=1',
    '--pbcx=1', '--pbcy=1', '--pbcz=0',
    "--runtimelimit=$RuntimeLimit",
    "--sweeps=$Sweeps",
    "--thermalization=$Thermalization",
    "--outputfile=$JobName.out.h5",
    "--checkpoint=$JobName.clone.h5"
)

if ($PotentialFile) { $simulationArgs += "--onsite_potential_file=$PotentialFile" }
if ($DisablePotential) { $simulationArgs += "--onsite_potential_enabled=0" }

Write-Host "Running ${JobName}: ${Lx}x${Ly}, $Processes process(es), $Sweeps sweeps, canonical=$Canonical, programmable=$ProgrammablePotential, disablePotential=$DisablePotential."
$simulationExit = $null
Push-Location -LiteralPath $OutputDirectory
try {
    if ($Processes -gt 1) {
        & $mpiExec -n $Processes $mpiExe @simulationArgs
    } else {
        & $serialExe @simulationArgs
    }
    $simulationExit = $LASTEXITCODE
}
finally {
    Pop-Location
}
if ($simulationExit -ne 0) { throw "Simulation failed with exit code $simulationExit." }

Write-Host "Completed: $(Join-Path $OutputDirectory "$JobName.out.h5")"
Write-Host "Checkpoint base: $(Join-Path $OutputDirectory "$JobName.clone.h5")"
