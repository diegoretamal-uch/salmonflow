#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
# SalmonFlow — Quick-start script
# Usage: ./run.sh [FASTQ_DIR] [REFERENCES_DIR] [OUTPUT_DIR]
# ──────────────────────────────────────────────────────────────
set -euo pipefail

FASTQ_DIR="${1:-$(pwd)/data/input}"
REF_DIR="${2:-$(pwd)/data/references}"
OUT_DIR="${3:-$(pwd)/data/output}"

echo "╔═══════════════════════════════════════════╗"
echo "║           🐟  SalmonFlow  🐟              ║"
echo "╠═══════════════════════════════════════════╣"
echo "║  FASTQs:      $FASTQ_DIR"
echo "║  References:  $REF_DIR"
echo "║  Output:      $OUT_DIR"
echo "╚═══════════════════════════════════════════╝"
echo ""
echo "Starting container... Open http://localhost:3838"
echo ""

docker run --rm -p 3838:3838 \
  -v "${FASTQ_DIR}:/data/input" \
  -v "${REF_DIR}:/data/references" \
  -v "${OUT_DIR}:/data/output" \
  salmonflow
