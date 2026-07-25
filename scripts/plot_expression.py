#!/usr/bin/env python
"""The top expressed genes of the matched RNA-seq arm, as a gene by library heatmap.

The quantity is magnitude, so the colour ramp is one hue running light to dark — never a rainbow,
whose bands invent boundaries the data does not have. Values span four orders of magnitude, so the
cells are shaded on log10 TPM while the printed number stays the TPM itself; the shading orders the
matrix and the number is the datum. Every cell carries its value, so the figure is also the table.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# One hue, light to dark. Blues is a single-hue sequential ramp; the darkest end is the largest
# value, which is the direction a reader assumes without being told.
_RAMP = "Blues"
_INK = "#1c1c1c"
_MUTED = "#5a5a5a"


def _library_label(name: str) -> str:
    """`calu6_g418_rep1_rnaseq` reads as `G418 rep1` once the constant parts are dropped."""

    parts = name.replace("_rnaseq", "").split("_")
    condition = parts[1].replace("untreated", "untreated").upper() if len(parts) > 1 else name
    replicate = parts[-1] if parts[-1].startswith("rep") else ""
    return f"{condition}\n{replicate}" if replicate else condition


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", required=True, type=Path, help="top_expressed.tsv")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--title", default="Most expressed genes, matched RNA-seq (Calu-6)", help="Figure title."
    )
    parser.add_argument(
        "--note",
        default="Transcripts per million over union exons; shaded on log10, labelled with TPM.",
        help="The line under the figure.",
    )
    args = parser.parse_args()

    table = pd.read_csv(args.top, sep="\t")
    libraries = [c for c in table.columns if c not in ("gene_id", "gene_symbol", "mean_tpm")]
    values = table[libraries].to_numpy(dtype=float)
    shaded = np.log10(values + 1)

    height = 0.34 * len(table) + 2.2
    figure, axes = plt.subplots(figsize=(1.35 * len(libraries) + 3.6, height))
    image = axes.imshow(shaded, cmap=_RAMP, aspect="auto", vmin=0, vmax=shaded.max())

    axes.set_xticks(range(len(libraries)), [_library_label(c) for c in libraries], fontsize=9)
    axes.set_yticks(range(len(table)), table["gene_symbol"], fontsize=9)
    axes.tick_params(length=0)
    for edge in axes.spines.values():
        edge.set_visible(False)

    # A 2 px surface gap between cells, so adjacent fills read as separate marks.
    axes.set_xticks(np.arange(-0.5, len(libraries), 1), minor=True)
    axes.set_yticks(np.arange(-0.5, len(table), 1), minor=True)
    axes.grid(which="minor", color="white", linewidth=2)
    axes.tick_params(which="minor", length=0)

    # The value in every cell, inked light or dark by the shade under it rather than by a fixed
    # choice that would vanish at one end of the ramp.
    ceiling = shaded.max()
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            value = values[row, column]
            text = f"{value:,.0f}" if value >= 10 else f"{value:.1f}"
            dark = shaded[row, column] > 0.6 * ceiling
            axes.text(
                column,
                row,
                text,
                ha="center",
                va="center",
                fontsize=8,
                color="white" if dark else _INK,
            )

    bar = figure.colorbar(image, ax=axes, fraction=0.025, pad=0.02)
    bar.set_label("log₁₀ (TPM + 1)", fontsize=9, color=_MUTED)
    bar.ax.tick_params(labelsize=8, length=0, colors=_MUTED)
    bar.outline.set_visible(False)

    axes.set_title(
        args.title,
        fontsize=12,
        color=_INK,
        pad=12,
        loc="left",
    )
    figure.text(
        0.01,
        0.005,
        args.note,
        fontsize=8,
        color=_MUTED,
    )
    figure.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
