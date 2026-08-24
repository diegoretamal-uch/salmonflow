#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
# SalmonFlow — Quick-start script
# Usage: ./run.sh [FASTQ_DIR] [REFERENCES_DIR] [OUTPUT_DIR]
#
# Paths default to data/ inside this folder. Any absolute path works.
#
# If the iDEP image is present locally, iDEP is started alongside
# SalmonFlow so the Results tab can hand off to it. Set
# SALMONFLOW_NO_IDEP=1 to skip that.
# ──────────────────────────────────────────────────────────────
set -euo pipefail

FASTQ_DIR="${1:-$(pwd)/data/input}"
REF_DIR="${2:-$(pwd)/data/references}"
OUT_DIR="${3:-$(pwd)/data/output}"

IDEP_IMAGE="${IDEP_IMAGE:-gexijin/idep:latest}"
IDEP_PORT="${IDEP_PORT:-3839}"
NETWORK="salmonflow-net"

# Create directories if they don't exist
mkdir -p "$FASTQ_DIR" "$REF_DIR" "$OUT_DIR" "$(pwd)/data/tmp"

# ── iDEP (optional) ───────────────────────────────────────────
# Started only if the image is already pulled — it is large, and we
# never trigger that download implicitly.
IDEP_STATUS="disabled (SALMONFLOW_NO_IDEP=1)"

if [ "${SALMONFLOW_NO_IDEP:-0}" != "1" ]; then
  if docker image inspect "$IDEP_IMAGE" >/dev/null 2>&1; then

    docker network inspect "$NETWORK" >/dev/null 2>&1 \
      || docker network create "$NETWORK" >/dev/null

    if [ -n "$(docker ps -q -f name='^idep$')" ]; then
      IDEP_STATUS="already running on http://localhost:${IDEP_PORT}"
    else
      # Remove a stopped leftover so the name is free
      docker rm -f idep >/dev/null 2>&1 || true
      if docker run -d --name idep --network "$NETWORK" \
           -p "${IDEP_PORT}:3838" "$IDEP_IMAGE" >/dev/null 2>&1; then
        IDEP_STATUS="started on http://localhost:${IDEP_PORT}"
      else
        IDEP_STATUS="failed to start (is port ${IDEP_PORT} already in use?)"
      fi
    fi
  else
    IDEP_STATUS="not installed — run: docker pull ${IDEP_IMAGE}"
  fi
fi

# Put SalmonFlow on the shared network only if it exists, so the app
# can reach iDEP by hostname for its status indicator.
NET_ARGS=()
if docker network inspect "$NETWORK" >/dev/null 2>&1; then
  NET_ARGS=(--network "$NETWORK")
fi

echo ""
echo "  SalmonFlow"
echo "  FASTQs:     $FASTQ_DIR"
echo "  References: $REF_DIR"
echo "  Output:     $OUT_DIR"
echo "  iDEP:       $IDEP_STATUS"
echo ""
echo "  Starting... Open http://localhost:3838"
echo ""

docker run --rm -p 3838:3838 \
  "${NET_ARGS[@]}" \
  -e "IDEP_PORT=${IDEP_PORT}" \
  -v "${FASTQ_DIR}:/data/input" \
  -v "${REF_DIR}:/data/references" \
  -v "${OUT_DIR}:/data/output" \
  -v "$(pwd)/data/tmp:/data/tmp" \
  salmonflow

# SalmonFlow runs in the foreground; iDEP stays up after it exits so a
# session can be finished in the browser. Stop it with:
#   docker rm -f idep
