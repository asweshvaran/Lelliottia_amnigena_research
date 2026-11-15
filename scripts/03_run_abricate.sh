#!/usr/bin/env bash
set -euo pipefail

BASE="/mnt/d/Ewingella_americana_research"
ABRICATE="/home/aswesh/mamba/envs/ewingella_ani/bin/abricate"
mkdir -p "$BASE/results/abricate"

export PATH="/home/aswesh/mamba/envs/ewingella_ani/bin:$PATH"

echo "Updating ABRicate databases..."
$ABRICATE --setupdb

for db in card ncbi resfinder vfdb plasmidfinder; do
    echo "Screening against $db..."
    $ABRICATE --db "$db" --threads 8 --minid 80 --mincov 60 "$BASE/data/genomes_fna/"*.fna > "$BASE/results/abricate/abricate_${db}.tsv"
    $ABRICATE --summary "$BASE/results/abricate/abricate_${db}.tsv" > "$BASE/results/abricate/summary_${db}.tsv"
done

echo "ABRicate screening complete!"
