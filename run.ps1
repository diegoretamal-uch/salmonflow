# SalmonFlow — Windows Quick-start (PowerShell)
# Usage: .\run.ps1 [FastqDir] [ReferencesDir] [OutputDir]
#
# Paths default to data\ inside this folder. Any absolute path works.
param(
    [string]$FastqDir = ".\data\input",
    [string]$RefDir   = ".\data\references",
    [string]$OutDir   = ".\data\output"
)

# Create directories if they don't exist
New-Item -ItemType Directory -Force -Path $FastqDir  | Out-Null
New-Item -ItemType Directory -Force -Path $RefDir    | Out-Null
New-Item -ItemType Directory -Force -Path $OutDir    | Out-Null
New-Item -ItemType Directory -Force -Path ".\data\tmp" | Out-Null

# Resolve to absolute paths (required for Docker volume mounts)
$FastqDir = (Resolve-Path $FastqDir).Path
$RefDir   = (Resolve-Path $RefDir).Path
$OutDir   = (Resolve-Path $OutDir).Path
$TmpDir   = (Resolve-Path ".\data\tmp").Path

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
    -v "${TmpDir}:/data/tmp" `
    salmonflow
