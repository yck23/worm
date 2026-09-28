<#
.SYNOPSIS
    Run a temperature and lattice-size sweep for the 2D Bose-Hubbard model.
.DESCRIPTION
    Reads one user-edited Bose-Hubbard INI file, creates a complete parameter
    file for every (size, temperature) point, and invokes run_bose_square.ps1
    sequentially. Temperature is converted to inverse temperature as beta=1/T.
    With -Canonical, the exact particle number is rounded from Density*Lx*Ly
    independently for every lattice size and the canonical executable is used.
.EXAMPLE
    .\scripts\run_temperature_sweep.ps1 -ParameterFile .\parameter_files\BoseHubbard.ini -Temperatures 0.35,0.40,0.45 -Sizes 8,12,16 -Processes 24 -SweepName bkt_coarse
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory)] [double[]]$Temperatures,
    [Parameter(Mandatory)] [int[]]$Sizes,
    [ValidateRange(1, 32)] [int]$Processes = 1,
    [ValidateRange(1, 2147483647)] [int]$RuntimeLimit = 3600,
    [ValidateRange(1, 2147483647)] [int]$Sweeps = 100000,
    [ValidateRange(0, 2147483647)] [int]$Thermalization = 1000,
    [ValidateRange(0, 4294967295)] [uint32]$BaseSeed = 1000,
    [switch]$Canonical,
    [double]$Density = 0.0,
    [string]$ParameterFile,
    [string]$SweepName = 'bkt_temperature_sweep',
    [string]$OutputDirectory,
    [switch]$SkipExisting,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$runner = Join-Path $PSScriptRoot 'run_bose_square.ps1'
if (-not $ParameterFile) { $ParameterFile = Join-Path $repoRoot 'parameter_files\BoseHubbard.ini' }
$parameterPath = [IO.Path]::GetFullPath($ParameterFile)
$safeSweepName = [regex]::Replace($SweepName, '[^A-Za-z0-9_-]', '_')
if (-not $safeSweepName) { throw 'SweepName must contain at least one letter, number, underscore, or hyphen.' }
if (-not $OutputDirectory) { $OutputDirectory = Join-Path $repoRoot "results\$safeSweepName" }
$outputPath = [IO.Path]::GetFullPath($OutputDirectory)
$generatedParameterDir = Join-Path $outputPath 'parameters'
$manifestPath = Join-Path $outputPath 'sweep_manifest.csv'

foreach ($required in @($runner, $parameterPath)) {
    if (-not (Test-Path -LiteralPath $required -PathType Leaf)) { throw "Required file not found: $required" }
}
if ($Temperatures.Count -eq 0) { throw 'Provide at least one temperature.' }
if ($Sizes.Count -eq 0) { throw 'Provide at least one lattice size.' }
foreach ($temperature in $Temperatures) {
    if ([double]::IsNaN($temperature) -or [double]::IsInfinity($temperature) -or $temperature -le 0) {
        throw "Temperature must be finite and greater than zero: $temperature"
    }
}
foreach ($size in $Sizes) {
    if ($size -lt 3) { throw "Periodic square lattices require size 3 or greater: $size" }
}
if ($Canonical -and
    ([double]::IsNaN($Density) -or [double]::IsInfinity($Density) -or $Density -le 0)) {
    throw 'Canonical sweeps require a finite -Density greater than zero.'
}

New-Item -ItemType Directory -Force $outputPath, $generatedParameterDir | Out-Null
$baseParameters = [IO.File]::ReadAllText($parameterPath)
$invariant = [Globalization.CultureInfo]::InvariantCulture

function Set-IniValue {
    param(
        [Parameter(Mandatory)] [string]$Text,
        [Parameter(Mandatory)] [string]$Name,
        [Parameter(Mandatory)] [string]$Value
    )

    $pattern = '(?m)^\s*' + [regex]::Escape($Name) + '\s*=.*$'
    $expression = [regex]::new($pattern)
    $found = $expression.Matches($Text)
    if ($found.Count -gt 1) { throw "Parameter '$Name' occurs more than once in $parameterPath." }
    if ($found.Count -eq 1) { return $expression.Replace($Text, "$Name = $Value", 1) }

    $newline = if ($Text.Contains("`r`n")) { "`r`n" } else { "`n" }
    return $Text.TrimEnd() + $newline + "$Name = $Value" + $newline
}

function Save-Manifest {
    param([object[]]$Rows)
    @($Rows) | Export-Csv -LiteralPath $manifestPath -NoTypeInformation -Encoding UTF8
}

$rows = @()
$runIndex = 0
foreach ($size in $Sizes) {
    foreach ($temperature in $Temperatures) {
        $beta = 1.0 / $temperature
        $particles = $null
        if ($Canonical) {
            $particleValue = [math]::Round(
                $Density * [double]$size * [double]$size,
                0,
                [MidpointRounding]::AwayFromZero
            )
            if ($particleValue -lt 1 -or $particleValue -gt [int]::MaxValue) {
                throw "Invalid canonical particle number for L=$size and density=$Density."
            }
            $particles = [int]$particleValue
        }
        $temperatureText = $temperature.ToString('0.######', $invariant)
        $betaText = $beta.ToString('G17', $invariant)
        $temperatureLabel = $temperatureText.Replace('-', 'm').Replace('.', 'p')
        $jobName = '{0}_L{1:D3}_T{2}' -f $safeSweepName, $size, $temperatureLabel
        $seed64 = [uint64]$BaseSeed + [uint64]$runIndex
        if ($seed64 -gt [uint32]::MaxValue) { throw 'Generated seed exceeds the 32-bit range.' }
        $seed = [uint32]$seed64
        $generatedParameter = Join-Path $generatedParameterDir "$jobName.ini"
        $outputFile = Join-Path $outputPath "$jobName.out.h5"

        $pointParameters = $baseParameters
        $effectiveValues = [ordered]@{
            beta = $betaText
            seed = $seed.ToString($invariant)
            Lx = $size.ToString($invariant)
            Ly = $size.ToString($invariant)
            Lz = '1'
            pbcx = '1'
            pbcy = '1'
            pbcz = '0'
            runtimelimit = $RuntimeLimit.ToString($invariant)
            sweeps = $Sweeps.ToString($invariant)
            thermalization = $Thermalization.ToString($invariant)
        }
        foreach ($entry in $effectiveValues.GetEnumerator()) {
            $pointParameters = Set-IniValue -Text $pointParameters -Name $entry.Key -Value $entry.Value
        }
        if ($Canonical) {
            $pointParameters = Set-IniValue -Text $pointParameters -Name 'canonical' -Value $particles.ToString($invariant)
        }
        [IO.File]::WriteAllText($generatedParameter, $pointParameters, [Text.UTF8Encoding]::new($false))

        $row = [pscustomobject][ordered]@{
            Temperature = $temperatureText
            Beta = $betaText
            Lx = $size
            Ly = $size
            Canonical = [bool]$Canonical
            Density = $(if ($Canonical) { $Density.ToString('G17', $invariant) } else { '' })
            Particles = $(if ($Canonical) { $particles } else { '' })
            Seed = $seed
            Processes = $Processes
            JobName = $jobName
            ParameterFile = $generatedParameter
            OutputFile = $outputFile
            Status = 'planned'
        }
        $rows += $row
        Save-Manifest -Rows $rows

        if ($DryRun) {
            $row.Status = 'generated-only'
            Save-Manifest -Rows $rows
            Write-Host "Generated: L=$size, T=$temperatureText, beta=$betaText"
            $runIndex++
            continue
        }

        if (Test-Path -LiteralPath $outputFile) {
            if ($SkipExisting) {
                $row.Status = 'skipped-existing'
                Save-Manifest -Rows $rows
                Write-Host "Skipping existing point: L=$size, T=$temperatureText"
                $runIndex++
                continue
            }
            throw "Output already exists: $outputFile. Use -SkipExisting or choose a new SweepName."
        }

        Write-Host "Starting sweep point: L=$size, T=$temperatureText, beta=$betaText, seed=$seed"
        $runnerArguments = @(
            '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', $runner,
            '-ParameterFile', $generatedParameter,
            '-Lx', $size, '-Ly', $size,
            '-Processes', $Processes,
            '-RuntimeLimit', $RuntimeLimit,
            '-Sweeps', $Sweeps,
            '-Thermalization', $Thermalization,
            '-JobName', $jobName,
            '-OutputDirectory', $outputPath
        )
        if ($Canonical) { $runnerArguments += '-Canonical' }
        & powershell.exe @runnerArguments
        if ($LASTEXITCODE -ne 0) {
            $row.Status = "failed-$LASTEXITCODE"
            Save-Manifest -Rows $rows
            throw "Sweep point failed: L=$size, T=$temperatureText"
        }

        $row.Status = 'completed'
        Save-Manifest -Rows $rows
        $runIndex++
    }
}

Write-Host "Temperature sweep finished. Manifest: $manifestPath"
