#!/usr/bin/env bash
set -euo pipefail

export PATH="/home/aswesh/mamba/envs/ewingella_pan/bin:$PATH"
BASE="/mnt/d/Ewingella_americana_research"

echo "=== 1. Running Panaroo on 41 strict E. americana genomes ==="
rm -rf "$BASE/results/panaroo_strict_41"
panaroo \
    -i $(cat "$BASE/data/gff_panaroo/strict_41_list.txt") \
    -o "$BASE/results/panaroo_strict_41" \
    --clean-mode strict \
    --remove-invalid-genes \
    -a core \
    --aligner mafft \
    --core_threshold 0.95 \
    -t 8

echo "=== 2. Running Panaroo on 45 Ewingella-complex + Y. enterocolitica outgroup ==="
rm -rf "$BASE/results/panaroo_complex_45"
panaroo \
    -i $(cat "$BASE/data/gff_panaroo/complex_45_list.txt") \
    -o "$BASE/results/panaroo_complex_45" \
    --clean-mode strict \
    --remove-invalid-genes \
    -a core \
    --aligner mafft \
    --core_threshold 0.95 \
    -t 8

echo "Panaroo runs complete!"
