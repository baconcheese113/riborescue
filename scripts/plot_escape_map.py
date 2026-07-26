"""The cross-modality escape map as a report figure: base-editing reachability, and readthrough
overlaid on it as the continuous axis it is.

Reads results/escape_map.tsv and writes results/figures/escape_map.{png,svg}. Reachability is
geometric — whether a base editor can be placed on the stop — never eligibility.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Blue / orange are the CVD-safe categorical pair; grey is the remainder, not a series.
COLOURS = {
    "base_editable_exact": "#2a78d6",
    "base_editable_alternative": "#eb6834",
    "not_base_editable_under_panel": "#b8b6ad",
}
# Darkened text-tone of each fill, for labels that must clear contrast on the white surface.
INKED = {
    "base_editable_exact": "#1c5cab",
    "base_editable_alternative": "#b8471d",
    "not_base_editable_under_panel": "#6b6960",
}
LABELS = {
    "base_editable_exact": "Base-editable\nexact",
    "base_editable_alternative": "Base-editable\nalternative",
    "not_base_editable_under_panel": "Not editable\nunder panel",
}
# Single-word labels under the stacked segments, so a narrow segment's label does not run into its
# neighbour's.
SHORT = {
    "base_editable_exact": "Exact",
    "base_editable_alternative": "Alternative",
    "not_base_editable_under_panel": "Not editable",
}
ORDER = list(COLOURS)
INK, MUTED, SURFACE = "#1b1f27", "#5c6670", "#ffffff"


def plot(escape: pd.DataFrame, out: Path) -> None:
    scoreable = escape[escape["scoreable"]]
    total, n = len(escape), len(scoreable)
    counts = {c: int((scoreable["reach_class"] == c).sum()) for c in ORDER}

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "text.color": INK,
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
        }
    )
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1.4, 1]})
    fig.suptitle(
        "Cross-modality reachability of pathogenic nonsense variants",
        x=0.04,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )

    # Panel A — the scoreable set as one horizontal stacked bar, 2px white gaps between segments,
    # each segment directly labelled with its count, share, and class name. No separate legend.
    left = 0.0
    for c in ORDER:
        width = counts[c]
        a.barh(0, width, left=left, height=0.62, color=COLOURS[c], edgecolor=SURFACE, linewidth=2)
        mid = left + width / 2
        light = c != "not_base_editable_under_panel"
        a.text(
            mid,
            0.08,
            f"{width:,}",
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
            color="#ffffff" if light else INK,
        )
        a.text(
            mid,
            -0.13,
            f"{width / n * 100:.0f}%",
            ha="center",
            va="center",
            fontsize=9,
            color="#eaf1fb" if light else MUTED,
        )
        a.text(
            mid,
            -0.52,
            SHORT[c],
            ha="center",
            va="top",
            fontsize=9,
            color=INKED[c],
            fontweight="bold",
        )
        left += width
    a.set_xlim(0, n)
    a.set_ylim(-1.15, 0.7)
    a.axis("off")
    a.text(
        0,
        0.52,
        f"{n:,} variants placed on a stop-forming codon",
        color=INK,
        fontsize=11,
        fontweight="bold",
    )
    a.text(
        0,
        -1.02,
        f"Of {total:,} pathogenic nonsense variants, {total - n:,} could not be placed. Every "
        "reachable candidate\nuses ABE7.10 — reverting a stop needs a T→C edit on the opposite "
        "strand, and no stop codon\ncarries a cytosine for a C→T editor to act on.",
        color=MUTED,
        fontsize=8.5,
        linespacing=1.5,
    )

    # Panel B — readthrough stays continuous, shown per base-editing class rather than bucketed.
    b.set_title(
        "Best predicted readthrough, by class",
        loc="left",
        fontsize=11,
        fontweight="bold",
        color=INK,
        pad=8,
    )
    data, positions, colours = [], [], []
    for i, c in enumerate(ORDER):
        vals = (
            pd.to_numeric(
                scoreable.loc[scoreable["reach_class"] == c, "best_readthrough"], errors="coerce"
            ).dropna()
            * 100
        )
        data.append(vals.to_numpy())
        positions.append(i)
        colours.append(c)
    bp = b.boxplot(
        data,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        showfliers=False,
        orientation="horizontal",
        medianprops={"color": INK, "linewidth": 1.5},
        whiskerprops={"color": "#b0b4ac"},
        capprops={"color": "#b0b4ac"},
    )
    for patch, c in zip(bp["boxes"], colours, strict=True):
        patch.set_facecolor(COLOURS[c])
        patch.set_edgecolor(SURFACE)
        patch.set_linewidth(2)
        patch.set_alpha(0.6 if c == "not_base_editable_under_panel" else 0.92)
    b.set_yticks(range(len(ORDER)))
    b.set_yticklabels([LABELS[c].replace("\n", " ") for c in ORDER], fontsize=8.5)
    for tick, c in zip(b.get_yticklabels(), ORDER, strict=True):
        tick.set_color(INKED[c])
    b.set_xlabel("best predicted readthrough (%)", fontsize=9, color=MUTED)
    b.invert_yaxis()
    for spine in ("top", "right", "left"):
        b.spines[spine].set_visible(False)
    b.spines["bottom"].set_color("#d9ddd6")
    b.tick_params(left=False)
    b.grid(axis="x", color="#eceee9", linewidth=1)
    b.set_axisbelow(True)
    b.annotate(
        "Readthrough is a continuous overlay, not a route: a variant not base-editable under\n"
        "this panel may still have some predicted readthrough, and the reverse.",
        xy=(0, 0),
        xytext=(0.0, -0.30),
        textcoords="axes fraction",
        fontsize=8,
        color=MUTED,
        va="top",
        linespacing=1.5,
    )

    fig.text(
        0.04,
        0.015,
        "Geometric reachability under the BE4max + ABE7.10 / SpCas9-NGG panel: "
        "whether an editor can be placed on the stop. Not editing efficiency, off-target "
        "activity, delivery, tissue, splice, or eligibility.",
        fontsize=7.5,
        color=MUTED,
    )
    fig.subplots_adjust(left=0.04, right=0.985, top=0.85, bottom=0.24, wspace=0.42)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".svg"))
    print(f"wrote {out.with_suffix('.png')} and {out.with_suffix('.svg')}")


if __name__ == "__main__":
    frame = pd.read_csv("results/escape_map.tsv", sep="\t")
    plot(frame, Path("results/figures/escape_map"))
