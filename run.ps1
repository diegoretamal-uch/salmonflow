# SalmonFlow — Windows Quick-start (PowerShell)
# Usage: .\run.ps1 [FastqDir] [ReferencesDir] [OutputDir]
#
# Paths default to data\ inside this folder. Any absolute path works.
#
# If the iDEP image is present locally, iDEP is started alongside
# SalmonFlow so the Results tab can hand off to it. Set the environment
# variable SALMONFLOW_NO_IDEP=1 to skip that.
param(
    [string]$FastqDir = ".\data\input",
    [string]$RefDir   = ".\data\references",
    [string]$OutDir   = ".\data\output"
)

$IdepImage = if ($env:IDEP_IMAGE) { $env:IDEP_IMAGE } else { "gexijin/idep:latest" }
$IdepPort  = if ($env:IDEP_PORT)  { $env:IDEP_PORT }  else { "3839" }
$Network   = "salmonflow-net"

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

# ── iDEP (optional) ───────────────────────────────────────────
# Started only if the image is already pulled — it is large, and we
# never trigger that download implicitly.
$IdepStatus = "disabled (SALMONFLOW_NO_IDEP=1)"

if ($env:SALMONFLOW_NO_IDEP -ne "1") {
    docker image inspect $IdepImage 2>&1 | Out-Null
    if ($LASTEXITCODE -eq 0) {

        docker network inspect $Network 2>&1 | Out-Null
        if ($LASTEXITCODE -ne 0) { docker network create $Network | Out-Null }

        $running = docker ps -q -f "name=^idep$"
        if ($running) {
            $IdepStatus = "already running on http://localhost:$IdepPort"
        } else {
            docker rm -f idep 2>&1 | Out-Null
            docker run -d --name idep --network $Network `
                -p "${IdepPort}:3838" $IdepImage 2>&1 | Out-Null
            if ($LASTEXITCODE -eq 0) {
                $IdepStatus = "started on http://localhost:$IdepPort"
            } else {
                $IdepStatus = "failed to start (is port $IdepPort already in use?)"
            }
        }
    } else {
        $IdepStatus = "not installed - run: docker pull $IdepImage"
    }
}

# Put SalmonFlow on the shared network only if it exists, so the app
# can reach iDEP by hostname for its status indicator.
$NetArgs = @()
docker network inspect $Network 2>&1 | Out-Null
if ($LASTEXITCODE -eq 0) { $NetArgs = @("--network", $Network) }

Write-Host ""
Write-Host "  SalmonFlow" -ForegroundColor Cyan
Write-Host "  FASTQs:     $FastqDir"
Write-Host "  References: $RefDir"
Write-Host "  Output:     $OutDir"
Write-Host "  iDEP:       $IdepStatus"
Write-Host ""
Write-Host "  Starting... Open http://localhost:3838" -ForegroundColor Green
Write-Host ""

docker run --rm -p 3838:3838 `
    @NetArgs `
    -e "IDEP_PORT=$IdepPort" `
    -v "${FastqDir}:/data/input" `
    -v "${RefDir}:/data/references" `
    -v "${OutDir}:/data/output" `
    -v "${TmpDir}:/data/tmp" `
    salmonflow

# SalmonFlow runs in the foreground; iDEP stays up after it exits so a
# session can be finished in the browser. Stop it with:
#   docker rm -f idep
