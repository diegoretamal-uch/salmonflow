#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
# SalmonFlow — Quick-start script
# Usage: ./run.sh [FASTQ_DIR] [REFERENCES_DIR] [OUTPUT_DIR]
#
# Paths default to data/ inside this folder. Any absolute path works.
# ──────────────────────────────────────────────────────────────
set -euo pipefail

FASTQ_DIR="${1:-$(pwd)/data/input}"
REF_DIR="${2:-$(pwd)/data/references}"
OUT_DIR="${3:-$(pwd)/data/output}"

# Create directories if they don't exist
mkdir -p "$FASTQ_DIR" "$REF_DIR" "$OUT_DIR" "$(pwd)/data/tmp"

echo ""
echo "  SalmonFlow"
echo "  FASTQs:     $FASTQ_DIR"
echo "  References: $REF_DIR"
echo "  Output:     $OUT_DIR"
echo ""
echo "  Starting... Open http://localhost:3838"
echo ""

docker run --rm -p 3838:3838 \
  -v "${FASTQ_DIR}:/data/input" \
  -v "${REF_DIR}:/data/references" \
  -v "${OUT_DIR}:/data/output" \
  -v "$(pwd)/data/tmp:/data/tmp" \
  salmonflow
