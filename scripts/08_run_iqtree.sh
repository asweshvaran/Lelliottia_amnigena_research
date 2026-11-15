#!/usr/bin/env bash
set -euo pipefail

BASE="/mnt/d/Ewingella_americana_research"
mkdir -p "$BASE/results/phylogeny"

echo "=== 1. Running IQ-TREE 2 on 45 Ewingella-complex + Y. enterocolitica outgroup ==="
mkdir -p "$BASE/results/phylogeny/complex_45"
iqtree2 \
    -s "$BASE/results/panaroo_complex_45/core_gene_alignment_filtered.aln" \
    -m GTR+F+G4 \
    -B 1000 \
    -T AUTO \
    --prefix "$BASE/results/phylogeny/complex_45/ewingella_complex" \
    --redo

echo "=== 2. Running IQ-TREE 2 on 41 strict E. americana ==="
mkdir -p "$BASE/results/phylogeny/strict_41"
iqtree2 \
    -s "$BASE/results/panaroo_strict_41/core_gene_alignment_filtered.aln" \
    -m GTR+F+G4 \
    -B 1000 \
    -T AUTO \
    --prefix "$BASE/results/phylogeny/strict_41/ewingella_strict" \
    --redo

echo "IQ-TREE 2 phylogenomics complete!"
