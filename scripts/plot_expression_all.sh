#!/usr/bin/env bash
# Both expression heatmaps: what the libraries contain, and what the cell was doing once the
# surviving mitochondrial and structural RNA is set aside. Pixi's task shell rejects loops, so the
# two calls live here rather than inline in the manifest.
set -euo pipefail

python scripts/plot_expression.py \
  --top results/rnaseq/top_expressed.tsv \
  --out docs/figures/top_expressed.png

python scripts/plot_expression.py \
  --top results/rnaseq/top_expressed_nuclear.tsv \
  --out docs/figures/top_expressed_nuclear.png \
  --title "Most expressed nuclear genes, matched RNA-seq (Calu-6)" \
  --note "Mitochondrial and small structural RNA set aside; TPM over union exons, shaded on log10."
