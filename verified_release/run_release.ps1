$ErrorActionPreference = 'Stop'
$releaseRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$pythonExe = if ($env:PYTHON_EXE) { $env:PYTHON_EXE } else { 'python' }

$manifest = Import-Csv -LiteralPath (Join-Path $releaseRoot 'INPUT_MANIFEST.csv')
foreach ($row in $manifest) {
    $inputPath = Join-Path $releaseRoot $row.relative_path
    if (-not (Test-Path -LiteralPath $inputPath -PathType Leaf)) {
        throw "Missing required input: $($row.relative_path)"
    }
    $actual = (Get-FileHash -Algorithm SHA256 -LiteralPath $inputPath).Hash
    if ($actual -ne $row.sha256) {
        throw "SHA-256 mismatch: $($row.relative_path)"
    }
}

& $pythonExe (Join-Path $releaseRoot '05_scripts\83_crc_leakage_free_validation.py')
if ($LASTEXITCODE -ne 0) { throw 'CRC leakage-free analysis failed' }
& $pythonExe (Join-Path $releaseRoot '05_scripts\84_crc_singlecell_whole_atlas_pseudobulk.py')
if ($LASTEXITCODE -ne 0) { throw 'CRC whole-atlas single-cell analysis failed' }
& $pythonExe (Join-Path $releaseRoot '05_scripts\85_framework_boundary_benchmark.py')
if ($LASTEXITCODE -ne 0) { throw 'Framework boundary benchmark failed' }

$outputs = Get-ChildItem -LiteralPath (Join-Path $releaseRoot '06_results') -File -Recurse | Sort-Object FullName
$rows = foreach ($file in $outputs) {
    [pscustomobject]@{
        relative_path = $file.FullName.Substring($releaseRoot.Length + 1).Replace('\','/')
        sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $file.FullName).Hash
        bytes = $file.Length
    }
}
$rows | Export-Csv -LiteralPath (Join-Path $releaseRoot 'RUN_SHA256.csv') -NoTypeInformation -Encoding utf8
Write-Host "Release run complete: $($outputs.Count) output files"
