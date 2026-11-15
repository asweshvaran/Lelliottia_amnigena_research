#!/usr/bin/env bash
set -euo pipefail

export PATH="/home/aswesh/mamba/envs/ewingella_pan/bin:$PATH"
BASE="/mnt/d/Ewingella_americana_research"

echo "=== 1. Resuming Panaroo-MSA on 41 strict E. americana ==="
panaroo-msa \
    -o "$BASE/results/panaroo_strict_41" \
    -a core \
    --aligner mafft \
    --core_threshold 0.95 \
    --resume \
    -t 8

echo "=== 2. Running Panaroo on 45 Ewingella-complex + Y. enterocolitica outgroup (500-gene core subset) ==="
rm -rf "$BASE/results/panaroo_complex_45"
panaroo \
    -i $(cat "$BASE/data/gff_panaroo/complex_45_list.txt") \
    -o "$BASE/results/panaroo_complex_45" \
    --clean-mode strict \
    --remove-invalid-genes \
    -a core \
    --aligner mafft \
    --core_threshold 0.95 \
    --core_subset 500 \
    -t 8

echo "=== 3. Running IQ-TREE 2 on 41 strict E. americana ==="
mkdir -p "$BASE/results/phylogeny/strict_41"
iqtree2 \
    -s "$BASE/results/panaroo_strict_41/core_gene_alignment_filtered.aln" \
    -m GTR+F+G4 \
    -B 1000 \
    -T AUTO \
    --prefix "$BASE/results/phylogeny/strict_41/ewingella_strict" \
    --redo

echo "=== 4. Running IQ-TREE 2 on 45 Ewingella-complex + outgroup ==="
mkdir -p "$BASE/results/phylogeny/complex_45"
iqtree2 \
    -s "$BASE/results/panaroo_complex_45/core_gene_alignment_filtered.aln" \
    -m GTR+F+G4 \
    -B 1000 \
    -T AUTO \
    --prefix "$BASE/results/phylogeny/complex_45/ewingella_complex" \
    --redo

echo "ALL PHYLOGENOMICS COMPLETE!"
