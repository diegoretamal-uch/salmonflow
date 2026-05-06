# SalmonFlow — Windows Quick-start (PowerShell)
# Usage: .\run.ps1 [FastqDir] [ReferencesDir] [OutputDir]
param(
    [string]$FastqDir  = ".\data\input",
    [string]$RefDir    = ".\data\references",
    [string]$OutDir    = ".\data\output"
)

# Resolve to absolute paths (Docker Desktop requires them)
$FastqDir = (Resolve-Path $FastqDir).Path
$RefDir   = (Resolve-Path $RefDir).Path
$OutDir   = (Resolve-Path $OutDir).Path

Write-Host ""
Write-Host "  SalmonFlow" -ForegroundColor Cyan
Write-Host "  FASTQs:     $FastqDir"
Write-Host "  References: $RefDir"
Write-Host "  Output:     $OutDir"
Write-Host ""
Write-Host "  Starting... Open http://localhost:3838" -ForegroundColor Green
Write-Host ""

docker run --rm -p 3838:3838 `
    -v "${FastqDir}:/data/input" `
    -v "${RefDir}:/data/references" `
    -v "${OutDir}:/data/output" `
    salmonflow
