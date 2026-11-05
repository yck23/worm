<#
.SYNOPSIS
    Rebuild and package WORM for Windows/MinGW with square-lattice and MPI support.
.DESCRIPTION
    Pins ALPSCore, applies the repository's MinGW compatibility patch, disables
    ALPSCore tests and POSIX-only signal/stacktrace features, and collects runtime
    DLLs. Programmable potential builds use separate build and package directories.
.PARAMETER AlpsInstall
    Reuse an existing compatible ALPSCore installation without modifying it.
.PARAMETER AlpsRuntimeDirectory
    Directory containing the reused ALPSCore runtime DLLs. Defaults to the sibling
    alps-build directory because the pinned ALPSCore does not install Windows DLLs.
.EXAMPLE
    .\scripts\build_windows.ps1 -ProgrammablePotential -Jobs 4
#>
[CmdletBinding()]
param(
    [ValidateRange(1, 64)] [int]$Jobs = [Environment]::ProcessorCount,
    [switch]$ProgrammablePotential,
    [string]$AlpsInstall,
    [string]$AlpsRuntimeDirectory
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$workRoot = Join-Path $repoRoot 'build-repro'
$alpsSource = Join-Path $workRoot 'alps-src'
$alpsBuild = Join-Path $workRoot 'alps-build'
$reuseAlps = -not [string]::IsNullOrWhiteSpace($AlpsInstall)
if ($reuseAlps) {
    $AlpsInstall = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($AlpsInstall)
    if (-not $AlpsRuntimeDirectory) {
        $AlpsRuntimeDirectory = Join-Path (Split-Path -Parent $AlpsInstall) 'alps-build'
    }
} else {
    if ($AlpsRuntimeDirectory) { throw '-AlpsRuntimeDirectory requires -AlpsInstall.' }
    $AlpsInstall = Join-Path $workRoot 'alps-install'
    $AlpsRuntimeDirectory = $alpsBuild
}
$AlpsRuntimeDirectory = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($AlpsRuntimeDirectory)
$buildName = if ($ProgrammablePotential) { 'worm-square-programmable' } else { 'worm-square' }
$packageName = if ($ProgrammablePotential) { 'windows-square-programmable' } else { 'windows-square' }
$potentialOption = if ($ProgrammablePotential) { 'ON' } else { 'OFF' }
$wormBuildRelative = "build-repro/$buildName"
$wormCanonicalBuildRelative = "build-repro/$buildName-canonical"
$wormBuild = Join-Path $repoRoot $wormBuildRelative
$wormCanonicalBuild = Join-Path $repoRoot $wormCanonicalBuildRelative
$packageDir = Join-Path $repoRoot "dist\$packageName"
$patchFile = Join-Path $repoRoot 'patches\alpscore-mingw.patch'
$alpsCommit = '3606edfbabc44e0fa7d05efa2a50c6a8a481340f'
$cmake = 'C:\Program Files\CMake\bin\cmake.exe'
$objdump = 'C:\msys64\ucrt64\bin\objdump.exe'

foreach ($required in @($cmake, $objdump, $patchFile)) {
    if (-not (Test-Path -LiteralPath $required)) { throw "Required build input not found: $required" }
}
New-Item -ItemType Directory -Force $workRoot | Out-Null

function Invoke-InMappedRepo([string]$Command) {
    $mappedCommand = $Command.Replace('{ROOT}', '!CD!')
    # cmd cannot start with a UNC current directory. Start locally, then let
    # pushd map the repository. Merge native diagnostics inside cmd so captured
    # CMake/compiler warnings do not become terminating PowerShell errors.
    Push-Location -LiteralPath ([IO.Path]::GetTempPath())
    try {
        & cmd.exe /d /v:on /c ('pushd "{0}" && {1} 2>&1' -f $repoRoot, $mappedCommand)
        if ($LASTEXITCODE -ne 0) { throw "Command failed with exit code ${LASTEXITCODE}: $Command" }
    }
    finally { Pop-Location }
}

function ConvertTo-CMakePath([string]$Path) {
    # MinGW tools need a drive path, so use cmd's mapping for paths on this share.
    $rootUri = [Uri]($repoRoot.TrimEnd('\') + '\')
    $pathUri = [Uri]$Path
    if ($rootUri.IsUnc -and $pathUri.IsUnc -and
        $rootUri.Host -eq $pathUri.Host -and $rootUri.Segments[1] -eq $pathUri.Segments[1]) {
        return '{ROOT}/' + [Uri]::UnescapeDataString($rootUri.MakeRelativeUri($pathUri).ToString())
    }
    return $Path.Replace('\', '/')
}

if ($reuseAlps) {
    $alpsConfig = Join-Path $AlpsInstall 'share\ALPSCore\ALPSCoreConfig.cmake'
    if (-not (Test-Path -LiteralPath $alpsConfig)) { throw "ALPSCore installation not found: $alpsConfig" }
    if (-not (Test-Path -LiteralPath $AlpsRuntimeDirectory -PathType Container)) {
        throw "ALPSCore runtime directory not found: $AlpsRuntimeDirectory. Set -AlpsRuntimeDirectory."
    }
    Write-Host "Reusing ALPSCore installation: $AlpsInstall"
} else {
    if (-not (Test-Path -LiteralPath (Join-Path $alpsSource '.git'))) {
        Invoke-InMappedRepo 'git clone https://github.com/ALPSCore/ALPSCore.git build-repro/alps-src'
    }
    & git -C $alpsSource checkout --detach $alpsCommit
    if ($LASTEXITCODE -ne 0) { throw 'Could not check out the pinned ALPSCore commit.' }

    & git -C $alpsSource apply --reverse --check $patchFile 2>$null
    $patchAlreadyApplied = $LASTEXITCODE -eq 0
    if (-not $patchAlreadyApplied) {
        & git -C $alpsSource apply --check $patchFile
        if ($LASTEXITCODE -ne 0) { throw 'ALPSCore patch does not apply cleanly.' }
        & git -C $alpsSource apply $patchFile
        if ($LASTEXITCODE -ne 0) { throw 'Failed to apply ALPSCore patch.' }
    }

    $alpsConfigure = @(
        ('"{0}"' -f $cmake), '--fresh',
        '-S build-repro/alps-src', '-B build-repro/alps-build', '-G "MinGW Makefiles"',
        '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_INSTALL_PREFIX=build-repro/alps-install',
        '-DTesting=OFF', '-DDocumentation=OFF',
        '-DCMAKE_CXX_FLAGS="-DALPS_UTILITY_NO_STACKTRACE -DALPS_NO_SIGNALS"',
        '-DHDF5_ROOT=C:/msys64/ucrt64', '-DBOOST_ROOT=C:/msys64/ucrt64'
    ) -join ' '
    Invoke-InMappedRepo $alpsConfigure
    Invoke-InMappedRepo ('"{0}" --build build-repro/alps-build --parallel {1}' -f $cmake, $Jobs)
    Invoke-InMappedRepo ('"{0}" --install build-repro/alps-build' -f $cmake)
}

$alpsCmakeDir = ConvertTo-CMakePath (Join-Path $AlpsInstall 'share\ALPSCore')
foreach ($canonicalBuild in @($false, $true)) {
    $buildRelative = if ($canonicalBuild) { $wormCanonicalBuildRelative } else { $wormBuildRelative }
    $canonicalOption = if ($canonicalBuild) { 'ON' } else { 'OFF' }
    $wormConfigure = @(
        ('"{0}"' -f $cmake), '--fresh',
        '-S src', ('-B "{0}"' -f $buildRelative), '-G "MinGW Makefiles"',
        '-DCMAKE_BUILD_TYPE=Release', '-DLATTICE=square',
        "-DCWINDOW=$canonicalOption", "-DPROGRAMMABLE_POTENTIAL=$potentialOption",
        ('-DALPSCore_DIR="{0}"' -f $alpsCmakeDir),
        '-DBOOST_ROOT=C:/msys64/ucrt64', '-DHDF5_ROOT=C:/msys64/ucrt64',
        '-DCMAKE_CXX_FLAGS="-DALPS_UTILITY_NO_STACKTRACE -DALPS_NO_SIGNALS"'
    ) -join ' '
    Invoke-InMappedRepo $wormConfigure
    Invoke-InMappedRepo ('"{0}" --build "{1}" --parallel {2}' -f $cmake, $buildRelative, $Jobs)
}

# Replace only this variant's generated package after both builds succeed.
$resolvedPackage = [IO.Path]::GetFullPath($packageDir)
$resolvedDist = [IO.Path]::GetFullPath((Join-Path $repoRoot 'dist')).TrimEnd('\') + '\'
if (-not $resolvedPackage.StartsWith($resolvedDist, [StringComparison]::OrdinalIgnoreCase) -or
    (Split-Path -Leaf $resolvedPackage) -ne $packageName) {
    throw "Refusing to replace package outside the expected dist directory: $resolvedPackage"
}
if (Test-Path -LiteralPath $resolvedPackage) { Remove-Item -LiteralPath $resolvedPackage -Recurse -Force }
New-Item -ItemType Directory -Force $packageDir | Out-Null

Copy-Item -LiteralPath (Join-Path $wormBuild 'qmc_worm.exe') -Destination $packageDir -Force
Copy-Item -LiteralPath (Join-Path $wormBuild 'qmc_worm_mpi.exe') -Destination $packageDir -Force
Copy-Item -LiteralPath (Join-Path $wormCanonicalBuild 'qmc_worm.exe') -Destination (Join-Path $packageDir 'qmc_worm_canonical.exe') -Force
Copy-Item -LiteralPath (Join-Path $wormCanonicalBuild 'qmc_worm_mpi.exe') -Destination (Join-Path $packageDir 'qmc_worm_canonical_mpi.exe') -Force
Copy-Item -LiteralPath (Join-Path $repoRoot 'parameter_files\BoseHubbard.ini') -Destination $packageDir -Force
$alpsDlls = @(Get-ChildItem -LiteralPath $AlpsRuntimeDirectory -Recurse -Filter 'libalps*.dll' -File)
$alpsDlls | Copy-Item -Destination $packageDir -Force

# Recursively collect non-system MinGW DLL dependencies into the package.
$searchDirs = @('C:\msys64\ucrt64\bin') +
    @($alpsDlls | ForEach-Object DirectoryName | Select-Object -Unique)
$queue = [Collections.Generic.Queue[string]]::new()
Get-ChildItem -LiteralPath $packageDir -File | Where-Object Extension -in '.exe','.dll' | ForEach-Object { $queue.Enqueue($_.FullName) }
$inspected = [Collections.Generic.HashSet[string]]::new([StringComparer]::OrdinalIgnoreCase)
while ($queue.Count -gt 0) {
    $binary = $queue.Dequeue()
    if (-not $inspected.Add($binary)) { continue }
    $imports = & $objdump -p $binary 2>$null | ForEach-Object {
        if ($_ -match '^\s*DLL Name:\s*(.+?)\s*$') { $Matches[1] }
    }
    foreach ($dllName in $imports) {
        $destination = Join-Path $packageDir $dllName
        if (Test-Path -LiteralPath $destination) { continue }
        $source = $null
        foreach ($dir in $searchDirs) {
            $candidate = Join-Path $dir $dllName
            if (Test-Path -LiteralPath $candidate) { $source = $candidate; break }
        }
        if ($source) {
            Copy-Item -LiteralPath $source -Destination $destination -Force
            $queue.Enqueue($destination)
        }
    }
}

$manifest = Get-ChildItem -LiteralPath $packageDir -File |
    Where-Object Name -ne 'SHA256SUMS.txt' |
    Sort-Object Name |
    ForEach-Object {
        $hash = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash
        '{0} *{1}' -f $hash, $_.Name
    }
[IO.File]::WriteAllLines((Join-Path $packageDir 'SHA256SUMS.txt'), $manifest)
Write-Host "Build and packaging completed: $packageDir"
$runSwitch = if ($ProgrammablePotential) { ' -ProgrammablePotential' } else { '' }
Write-Host "Run a smoke test with: .\scripts\run_bose_square.ps1$runSwitch -Lx 4 -Ly 4 -Sweeps 200 -Thermalization 20 -RuntimeLimit 5"
