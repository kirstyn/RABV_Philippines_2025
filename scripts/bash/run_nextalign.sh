#!/usr/bin/env bash
set -euo pipefail

# Usage: ./run_nextalign.sh <input_fasta> <reference_fasta>
INPUT_FASTA=$1
REFERENCE_FASTA=$2

# Extract base name of input (without path and extension)
BASENAME=$(basename "$INPUT_FASTA" .fasta)

# Output folder
OUTPUT_DIR=$(dirname "$INPUT_FASTA")/${BASENAME}_nextalign
mkdir -p "$OUTPUT_DIR"

echo ">>> Running Nextalign..."
nextalign run \
  -r "${REFERENCE_FASTA}" \
  -O "${OUTPUT_DIR}" \
  "${INPUT_FASTA}"

# Nextalign outputs default filenames inside OUTPUT_DIR
# Rename files to prepend input basename
for f in "${OUTPUT_DIR}"/*; do
    fname=$(basename "$f")
    mv "$f" "${OUTPUT_DIR}/${BASENAME}.${fname}"
done

echo ">>> Alignment complete!"
echo "Output files written to: ${OUTPUT_DIR}"