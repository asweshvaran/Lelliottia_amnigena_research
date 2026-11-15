#!/usr/bin/env bash
set -euo pipefail

BASE="/mnt/d/Ewingella_americana_research"
mkdir -p "$BASE/results/ani"

# Create list of all 57 E. americana + 1 outgroup
ls "$BASE/data/genomes_fna/"*.fna "$BASE/data/outgroup_fna/"*.fna > "$BASE/results/ani/genome_list.txt"

echo "Running all-vs-all FastANI on $(wc -l < "$BASE/results/ani/genome_list.txt") genomes..."
/home/aswesh/mamba/envs/ewingella_ani/bin/fastANI \
    --ql "$BASE/results/ani/genome_list.txt" \
    --rl "$BASE/results/ani/genome_list.txt" \
    -t 8 \
    -o "$BASE/results/ani/fastani_all_vs_all.tsv"

echo "FastANI complete: $BASE/results/ani/fastani_all_vs_all.tsv"
