<#
.SYNOPSIS
    Run or continue an Rb-87 canonical or grand-canonical calculation.
.EXAMPLE
    .\scripts\run_rb87_bkt64.ps1 -Mode Fresh
.EXAMPLE
    .\scripts\run_rb87_bkt64.ps1 -Mode Production
.EXAMPLE
    .\scripts\run_rb87_bkt64.ps1 -Ensemble GrandCanonical -Mode Fresh
#>
[CmdletBinding()]
param(
    [ValidateSet('Fresh', 'Resume', 'Production')]
    [string]$Mode = 'Fresh',
    [ValidateSet('Canonical', 'GrandCanonical')]
    [string]$Ensemble = 'Canonical',
    [ValidateRange(1, 32)]
    [int]$Processes = 24,
    [string]$ParameterFile,
    [string]$OutputDirectory
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$packageDir = Join-Path $repoRoot 'dist\windows-square'
$isCanonical = $Ensemble -eq 'Canonical'
$executableName = if ($isCanonical) { 'qmc_worm_canonical_mpi.exe' } else { 'qmc_worm_mpi.exe' }
$executable = Join-Path $packageDir $executableName
$mpiExec = 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe'
if (-not $ParameterFile) {
    $parameterName = if ($isCanonical) { 'Rb87_BKT_64.ini' } else { 'Rb87_BKT_64_grand_canonical.ini' }
    $ParameterFile = Join-Path $repoRoot "parameter_files\$parameterName"
}
if (-not $OutputDirectory) {
    $outputName = if ($isCanonical) { 'Rb87_BKT_64' } else { 'Rb87_BKT_64_grand_canonical' }
    $OutputDirectory = Join-Path $repoRoot "results\$outputName"
}

$parameterPath = [IO.Path]::GetFullPath($ParameterFile)
$outputPath = [IO.Path]::GetFullPath($OutputDirectory)

function Get-IniFileName {
    param(
        [Parameter(Mandatory)] [string]$Text,
        [Parameter(Mandatory)] [string]$Name
    )

    $pattern = '(?m)^\s*' + [regex]::Escape($Name) + '\s*=\s*(.+?)\s*$'
    $matches = [regex]::Matches($Text, $pattern)
    if ($matches.Count -ne 1) { throw "Parameter '$Name' must occur exactly once in $parameterPath." }
    $value = $matches[0].Groups[1].Value.Trim().Trim('"').Trim("'")
    if (-not $value -or [IO.Path]::GetFileName($value) -ne $value) {
        throw "Parameter '$Name' must be a filename without a directory: $value"
    }
    return $value
}

function Get-IniValue {
    param(
        [Parameter(Mandatory)] [string]$Text,
        [Parameter(Mandatory)] [string]$Name
    )

    $pattern = '(?m)^\s*' + [regex]::Escape($Name) + '\s*=\s*(.+?)\s*$'
    $matches = [regex]::Matches($Text, $pattern)
    if ($matches.Count -ne 1) { throw "Parameter '$Name' must occur exactly once in $parameterPath." }
    return $matches[0].Groups[1].Value.Trim().Trim('"').Trim("'")
}

$parameterText = [IO.File]::ReadAllText($parameterPath)
$checkpointName = Get-IniFileName -Text $parameterText -Name 'checkpoint'
$outputName = Get-IniFileName -Text $parameterText -Name 'outputfile'
$canonicalText = Get-IniValue -Text $parameterText -Name 'canonical'
if ($canonicalText -notmatch '^-?\d+$') { throw "Parameter 'canonical' must be an integer: $canonicalText" }
$canonicalValue = [int]$canonicalText
if ($isCanonical -and $canonicalValue -lt 0) {
    throw "Canonical mode requires canonical >= 0, but the file contains $canonicalValue."
}
if (-not $isCanonical -and $canonicalValue -ne -1) {
    throw "GrandCanonical mode requires canonical = -1, but the file contains $canonicalValue."
}
$checkpointPath = Join-Path $outputPath $checkpointName
$outputFile = Join-Path $outputPath $outputName

foreach ($required in @($executable, $mpiExec, $parameterPath)) {
    if (-not (Test-Path -LiteralPath $required -PathType Leaf)) {
        throw "Required file not found: $required. Run scripts\build_windows.ps1 first."
    }
}

New-Item -ItemType Directory -Force $outputPath | Out-Null
$tmpDir = Join-Path $outputPath 'tmp'
New-Item -ItemType Directory -Force $tmpDir | Out-Null
$env:TMPDIR = ($tmpDir -replace '\\', '/') + '/'
$env:Path = "$packageDir;$env:Path"

if ($Mode -eq 'Fresh') {
    foreach ($existing in @($checkpointPath, $outputFile)) {
        if (Test-Path -LiteralPath $existing) {
            throw "Output already exists: $existing. Use -Mode Resume or -Mode Production, or choose another OutputDirectory."
        }
    }
    $simulationArgs = @($parameterPath)
}
else {
    foreach ($rank in 0..($Processes - 1)) {
        $rankCheckpoint = if ($rank -eq 0) { $checkpointPath } else { "$checkpointPath.$rank" }
        if (-not (Test-Path -LiteralPath $rankCheckpoint -PathType Leaf)) {
            throw "Missing checkpoint for rank ${rank}: $rankCheckpoint. Continue with the same process count used for Fresh."
        }
    }

    if (Test-Path -LiteralPath $outputFile -PathType Leaf) {
        $archiveDir = Join-Path $outputPath 'archive'
        New-Item -ItemType Directory -Force $archiveDir | Out-Null
        $archiveName = '{0}_{1}' -f (Get-Date -Format 'yyyyMMdd_HHmmss_fff'), $outputName
        Copy-Item -LiteralPath $outputFile -Destination (Join-Path $archiveDir $archiveName)
    }

    $simulationArgs = @($checkpointName)
    if ($Mode -eq 'Production') { $simulationArgs += '--reset-statistics' }
}

$simulationExit = $null
Push-Location -LiteralPath $outputPath
try {
    Write-Host "Ensemble=$Ensemble, mode=$Mode, processes=$Processes, executable=$executableName, output=$outputPath"
    & $mpiExec -n $Processes $executable @simulationArgs
    $simulationExit = $LASTEXITCODE
}
finally {
    Pop-Location
}
if ($simulationExit -ne 0) { throw "Simulation failed with exit code $simulationExit." }

Write-Host "Completed $Mode."
Write-Host "Result: $outputFile"
Write-Host "Checkpoint: $checkpointPath"
