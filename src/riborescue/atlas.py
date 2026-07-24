"""The native-stop safety atlas: where readthrough carries a ribosome past a *normal* stop codon.

Two layers are kept apart on purpose, because one is measured and one is extrapolated.

The empirical anchor here is the readthrough assay turned on its subject, not its control. The
same extension window that measured a drug's effect past premature stops measures, per transcript,
how far a ribosome runs past the transcript's own native stop — so the GSE144140 libraries already
carry, for every transcript deep enough to score, the downstream occupancy G418 induces at a
stop. That is an anchor, not ground truth: it covers only what HEK293T expressed at depth, and its
per-transcript counts are small.

The predicted layer — scoring every MANE canonical stop with the amenability model — lives elsewhere
and is never fused into this one. The model was trained around premature stops; a canonical stop is
out of its distribution, and fusing a measurement with an extrapolation into one safety number
would hide exactly the disagreement worth seeing.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd

from riborescue.readthrough import _sequences

__all__ = [
    "AMINO_ACIDS",
    "native_stop_occupancy",
    "translate_extension",
]

# The standard genetic code, stops as "*". Only what a C-terminal extension needs.
AMINO_ACIDS = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def native_stop_occupancy(
    counts: pd.DataFrame,
    arms: Mapping[str, Sequence[str]],
    lengths: Sequence[int],
    excluded_transcripts: frozenset[str] = frozenset(),
    min_depth: int = 100,
    min_utr3: int = 50,
) -> pd.DataFrame:
    """Per transcript, the in-frame occupancy past its native stop, in each treatment arm.

    `arms` maps an arm name to its libraries; `lengths` is the footprint set the contrast used. The
    occupancy is pooled across an arm's replicates, so a transcript carries one number per arm,
    than one per library. Depth — the coding-frame P-sites the ratio rests on — travels with it, and
    a transcript is marked included only where every arm clears `min_depth`, because an occupancy
    resting on a handful of reads is noise dressed as a measurement.

    The same exclusions the readthrough assay applies are applied here, and for the same reason: a
    3' untranslated region that runs into another gene's coding sequence fills the downstream window
    with that gene's ribosomes, reading as native-stop readthrough that never happened. Those
    transcripts, the genes that read through their stop natively, and windows too short to separate
    the termination peak from what follows are dropped rather than reported as susceptible.
    """

    kept = counts[
        counts["length"].isin(list(lengths))
        & ~counts["transcript"].isin(excluded_transcripts)
        & (counts["l_utr3"] >= min_utr3)
    ]
    per_arm: dict[str, pd.DataFrame] = {}
    for arm, libraries in arms.items():
        rows = kept[kept["sample"].isin(list(libraries))]
        pooled = rows.groupby("transcript")[["extension_frame0", "cds_frame0"]].sum()
        ext0, cds0 = pooled["extension_frame0"], pooled["cds_frame0"]
        # A stop with no downstream reads is perfectly tight — occupancy 0, not missing. Excluding
        # those would keep only the stops that already leak. The relative standard error is a
        # Poisson error on a ratio of counts, and it is what is genuinely undefined at zero.
        rel_se = (1 / ext0.where(ext0 > 0) + 1 / cds0.where(cds0 > 0)) ** 0.5
        per_arm[arm] = pd.DataFrame(
            {
                f"{arm}_occupancy": ext0 / cds0.where(cds0 > 0),
                f"{arm}_depth": cds0,
                f"{arm}_rel_se": rel_se,
            }
        )

    atlas = pd.concat(per_arm.values(), axis=1).reset_index()

    depth_columns = [f"{arm}_depth" for arm in arms]
    atlas[depth_columns] = atlas[depth_columns].fillna(0).astype(int)
    atlas["included"] = (atlas[depth_columns] >= min_depth).all(axis=1)
    return atlas.sort_values("transcript").reset_index(drop=True)


def translate_extension(
    transcripts: Path,
    annotation: pd.DataFrame,
    extensions: pd.DataFrame,
) -> pd.DataFrame:
    """The peptide a readthrough ribosome would add past each native stop, and where it ends.

    The extension is new sequence: residues translated from just past the native stop to the next
    in-frame stop. It is not part of the annotated protein and does not overlap any domain of it —
    what it is, is a novel C-terminal tail whose own sequence is the thing worth reading.
    """

    sequences = _sequences(transcripts)
    windows = extensions.dropna(subset=["extension"]).set_index("transcript")["extension"].to_dict()
    coords = annotation.set_index("transcript")[["l_utr5", "l_cds"]].to_dict("index")

    rows = []
    for transcript, extension in windows.items():
        sequence = sequences.get(str(transcript))
        meta = coords.get(str(transcript))
        if sequence is None or meta is None:
            continue
        # The coding sequence ends with its stop; the extension begins right after it.
        start = int(meta["l_utr5"]) + int(meta["l_cds"])
        window = int(extension) - 3
        codons = [sequence[i : i + 3].upper() for i in range(start, start + window, 3)]
        peptide = "".join(AMINO_ACIDS.get(codon, "X") for codon in codons)
        rows.append(
            {
                "transcript": str(transcript),
                "extension_nt": int(extension),
                "extension_aa": len(peptide),
                "extension_peptide": peptide,
            }
        )
    return pd.DataFrame(rows)
