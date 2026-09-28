<#
.SYNOPSIS
    Rebuild and package WORM for Windows/MinGW with square-lattice and MPI support.
.DESCRIPTION
    Pins ALPSCore, applies the repository's MinGW compatibility patch, disables
    ALPSCore tests and POSIX-only signal/stacktrace features, builds with all
    requested workers, and collects runtime DLLs into dist\windows-square.
#>
[CmdletBinding()]
param(
    [ValidateRange(1, 64)] [int]$Jobs = [Environment]::ProcessorCount
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$workRoot = Join-Path $repoRoot 'build-repro'
$alpsSource = Join-Path $workRoot 'alps-src'
$alpsBuild = Join-Path $workRoot 'alps-build'
$alpsInstall = Join-Path $workRoot 'alps-install'
$wormBuild = Join-Path $workRoot 'worm-square'
$wormCanonicalBuild = Join-Path $workRoot 'worm-square-canonical'
$packageDir = Join-Path $repoRoot 'dist\windows-square'
$patchFile = Join-Path $repoRoot 'patches\alpscore-mingw.patch'
$alpsCommit = '3606edfbabc44e0fa7d05efa2a50c6a8a481340f'
$cmake = 'C:\Program Files\CMake\bin\cmake.exe'
$objdump = 'C:\msys64\ucrt64\bin\objdump.exe'

foreach ($required in @($cmake, $objdump, $patchFile)) {
    if (-not (Test-Path $required)) { throw "Required build input not found: $required" }
}
New-Item -ItemType Directory -Force $workRoot | Out-Null

function Invoke-InMappedRepo([string]$Command) {
    $mappedCommand = $Command.Replace('{ROOT}', '!CD!')
    & cmd.exe /d /v:on /c ('pushd "{0}" && {1}' -f $repoRoot, $mappedCommand)
    if ($LASTEXITCODE -ne 0) { throw "Command failed with exit code ${LASTEXITCODE}: $Command" }
}

if (-not (Test-Path (Join-Path $alpsSource '.git'))) {
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
    '"C:\Program Files\CMake\bin\cmake.exe"',
    '--fresh',
    '-S build-repro/alps-src', '-B build-repro/alps-build', '-G "MinGW Makefiles"',
    '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_INSTALL_PREFIX=build-repro/alps-install',
    '-DTesting=OFF', '-DDocumentation=OFF',
    '-DCMAKE_CXX_FLAGS="-DALPS_UTILITY_NO_STACKTRACE -DALPS_NO_SIGNALS"',
    '-DHDF5_ROOT=C:/msys64/ucrt64', '-DBOOST_ROOT=C:/msys64/ucrt64'
) -join ' '
Invoke-InMappedRepo $alpsConfigure
Invoke-InMappedRepo ('"{0}" --build build-repro/alps-build --parallel {1}' -f $cmake, $Jobs)
Invoke-InMappedRepo ('"{0}" --install build-repro/alps-build' -f $cmake)

$wormConfigure = @(
    '"C:\Program Files\CMake\bin\cmake.exe"',
    '--fresh',
    '-S src', '-B build-repro/worm-square', '-G "MinGW Makefiles"',
    '-DCMAKE_BUILD_TYPE=Release', '-DLATTICE=square',
    '-DALPSCore_DIR={ROOT}/build-repro/alps-install/share/ALPSCore',
    '-DBOOST_ROOT=C:/msys64/ucrt64', '-DHDF5_ROOT=C:/msys64/ucrt64',
    '-DCMAKE_CXX_FLAGS="-DALPS_UTILITY_NO_STACKTRACE -DALPS_NO_SIGNALS"'
) -join ' '
Invoke-InMappedRepo $wormConfigure
Invoke-InMappedRepo ('"{0}" --build build-repro/worm-square --parallel {1}' -f $cmake, $Jobs)

$wormCanonicalConfigure = @(
    '"C:\Program Files\CMake\bin\cmake.exe"',
    '--fresh',
    '-S src', '-B build-repro/worm-square-canonical', '-G "MinGW Makefiles"',
    '-DCMAKE_BUILD_TYPE=Release', '-DLATTICE=square', '-DCWINDOW=ON',
    '-DALPSCore_DIR={ROOT}/build-repro/alps-install/share/ALPSCore',
    '-DBOOST_ROOT=C:/msys64/ucrt64', '-DHDF5_ROOT=C:/msys64/ucrt64',
    '-DCMAKE_CXX_FLAGS="-DALPS_UTILITY_NO_STACKTRACE -DALPS_NO_SIGNALS"'
) -join ' '
Invoke-InMappedRepo $wormCanonicalConfigure
Invoke-InMappedRepo ('"{0}" --build build-repro/worm-square-canonical --parallel {1}' -f $cmake, $Jobs)

# Replace only the generated package after both builds have succeeded, so stale
# DLLs and old checksum manifests cannot survive an incremental rebuild.
if (Test-Path $packageDir) { Remove-Item -LiteralPath $packageDir -Recurse -Force }
New-Item -ItemType Directory -Force $packageDir | Out-Null

Copy-Item (Join-Path $wormBuild 'qmc_worm.exe') $packageDir -Force
Copy-Item (Join-Path $wormBuild 'qmc_worm_mpi.exe') $packageDir -Force
Copy-Item (Join-Path $wormCanonicalBuild 'qmc_worm.exe') (Join-Path $packageDir 'qmc_worm_canonical.exe') -Force
Copy-Item (Join-Path $wormCanonicalBuild 'qmc_worm_mpi.exe') (Join-Path $packageDir 'qmc_worm_canonical_mpi.exe') -Force
Copy-Item (Join-Path $repoRoot 'parameter_files\BoseHubbard.ini') $packageDir -Force
Get-ChildItem $alpsBuild -Recurse -Filter 'libalps*.dll' | Copy-Item -Destination $packageDir -Force

# Recursively collect non-system MinGW DLL dependencies into the package.
$searchDirs = @('C:\msys64\ucrt64\bin') +
    @(Get-ChildItem $alpsBuild -Recurse -Filter 'libalps*.dll' | ForEach-Object DirectoryName | Select-Object -Unique)
$queue = [Collections.Generic.Queue[string]]::new()
Get-ChildItem $packageDir -File | Where-Object Extension -in '.exe','.dll' | ForEach-Object { $queue.Enqueue($_.FullName) }
$inspected = [Collections.Generic.HashSet[string]]::new([StringComparer]::OrdinalIgnoreCase)
while ($queue.Count -gt 0) {
    $binary = $queue.Dequeue()
    if (-not $inspected.Add($binary)) { continue }
    $imports = & $objdump -p $binary 2>$null | ForEach-Object {
        if ($_ -match '^\s*DLL Name:\s*(.+?)\s*$') { $Matches[1] }
    }
    foreach ($dllName in $imports) {
        $destination = Join-Path $packageDir $dllName
        if (Test-Path $destination) { continue }
        $source = $null
        foreach ($dir in $searchDirs) {
            $candidate = Join-Path $dir $dllName
            if (Test-Path $candidate) { $source = $candidate; break }
        }
        if ($source) {
            Copy-Item $source $destination -Force
            $queue.Enqueue($destination)
        }
    }
}

$manifest = Get-ChildItem $packageDir -File |
    Where-Object Name -ne 'SHA256SUMS.txt' |
    Sort-Object Name |
    ForEach-Object {
        $hash = (Get-FileHash $_.FullName -Algorithm SHA256).Hash
        '{0} *{1}' -f $hash, $_.Name
    }
[IO.File]::WriteAllLines((Join-Path $packageDir 'SHA256SUMS.txt'), $manifest)
Write-Host "Build and packaging completed: $packageDir"
Write-Host 'Run a smoke test with: .\scripts\run_bose_square.ps1 -Lx 4 -Ly 4 -Sweeps 200 -Thermalization 20 -RuntimeLimit 5'
