"""aenmd's NMD-escape verdict, read on the same MANE Select transcript as the rule tier.

aenmd (Klonowski et al., Bioinformatics 2023) is a real, independent implementation of the NMD
escape rules the ensemble computes by hand — the model tier's rule-based member (ADR-0017). It
annotates every Ensembl v105 transcript a variant overlaps; this keeps only the MANE Select
transcript, so its verdict sits beside the rule tier's on the same molecule, and every variant aenmd
did not score — because it filtered the variant out (splice overlap, non-coding on that transcript)
or lacks the transcript — is marked unavailable with the reason rather than silently dropped.

The join key is aenmd's own: ``{chrom}:{pos:09d}|{ref}|{alt}`` in NCBI (no-"chr") naming.
"""

from __future__ import annotations

import pandas as pd

__all__ = [
    "AENMD_ESCAPE_RULES",
    "aenmd_verdicts",
    "mane_ensembl_by_gene",
    "model_agreement",
    "read_aenmd_rules",
]

# The five escape rules aenmd reports, in its own column names. is_ptc is the precondition, not one.
AENMD_ESCAPE_RULES = ("is_last", "is_penultimate", "is_css_proximal", "is_single", "is_407plus")


def read_aenmd_rules(path: str) -> pd.DataFrame:
    """aenmd's per-variant x transcript rule table, Ensembl transcripts stripped of version."""

    rules = pd.read_csv(path, sep="\t")
    rules["transcript"] = rules["transcript"].str.split(".").str[0]
    return rules


def mane_ensembl_by_gene(summary_path: str) -> dict[int, str]:
    """NCBI gene id to its MANE Select Ensembl transcript (version stripped), from the MANE summary.

    The rule tier scores each gene on its MANE Select RefSeq transcript; aenmd reports Ensembl ids.
    MANE Select is a one-to-one RefSeq/Ensembl pairing, so this map is how the two tiers meet on the
    same molecule.
    """

    summary = pd.read_csv(summary_path, sep="\t", dtype=str)
    gene = summary["#NCBI_GeneID"].str.removeprefix("GeneID:").astype(int)
    ensembl = summary["Ensembl_nuc"].str.split(".").str[0]
    select = summary["MANE_status"].str.contains("MANE Select", na=False)
    return dict(zip(gene[select], ensembl[select], strict=True))


def aenmd_verdicts(
    rules: pd.DataFrame,
    nonsense: pd.DataFrame,
    mane_by_gene: dict[int, str],
) -> pd.DataFrame:
    """One row per input variant: aenmd's escape call on its MANE transcript, or why it has none.

    `nonsense` supplies the coordinates the key is built from (variant_id, gene_id, chrom, pos, ref,
    alt). A variant is available only when aenmd annotated its MANE Select transcript; escape is the
    OR of aenmd's five escape rules. Unavailable rows carry a reason: no MANE transcript for the
    gene, a transcript aenmd's Ensembl set lacked, or a variant aenmd filtered on that transcript.
    """

    variants = nonsense.drop_duplicates("variant_id").copy()
    variants["mane_enst"] = variants["gene_id"].astype(int).map(mane_by_gene)
    variants["key"] = (
        variants["chrom"].astype(str)
        + ":"
        + variants["pos"].astype(int).map(lambda p: f"{p:09d}")
        + "|"
        + variants["ref"].astype(str)
        + "|"
        + variants["alt"].astype(str)
    )

    scored = rules.drop_duplicates(["key", "transcript"])
    annotated_transcripts = set(rules["transcript"].unique())
    merged = variants.merge(
        scored, left_on=["key", "mane_enst"], right_on=["key", "transcript"], how="left"
    )

    escape_cols = [f"aenmd_{rule}" for rule in AENMD_ESCAPE_RULES]
    merged = merged.rename(columns={rule: f"aenmd_{rule}" for rule in AENMD_ESCAPE_RULES})
    available = merged["is_ptc"].notna()
    out = pd.DataFrame({"variant_id": merged["variant_id"], "aenmd_available": available})

    # Unavailability is split by cause, because the causes are different claims. A MANE transcript
    # aenmd's older Ensembl set never carried (transcript_absent) is a build-version gap; one it
    # carried but did not annotate this variant on (variant_filtered) is aenmd's own splice/coding
    # filtering. Only the second is a statement about the variant.
    reason = pd.Series("", index=merged.index)
    reason[merged["mane_enst"].isna()] = "no_mane_ensembl"
    absent = merged["mane_enst"].notna() & ~merged["mane_enst"].isin(annotated_transcripts)
    reason[~available & absent] = "transcript_absent"
    reason[~available & merged["mane_enst"].isin(annotated_transcripts)] = "variant_filtered"
    out["reason"] = reason.values

    out["aenmd_is_ptc"] = merged["is_ptc"].fillna(False).astype(bool)
    for col in escape_cols:
        out[col] = merged[col].fillna(False).astype(bool)
    out["aenmd_escape"] = out[escape_cols].any(axis=1) & available
    return out


def model_agreement(rule_predictors: pd.DataFrame, verdicts: pd.DataFrame) -> dict:
    """How the hand-rolled rule tier and the aenmd tool compare, over variants aenmd scored.

    The two implement the same published rules, so they should largely agree; the value is in the
    exceptions, which come from aenmd's own edge-case handling (its splice filtering, transcript
    set, exact boundary conventions) rather than a different idea of NMD. The counts here are the
    concrete check that the ensemble's `full_rules` reproduces a real, independent implementation.
    """

    available = verdicts[verdicts["aenmd_available"].fillna(False).astype(bool)]
    merged = rule_predictors.merge(available, on="variant_id", how="inner")
    total = len(merged)
    summary: dict[str, object] = {
        "scoreable": len(rule_predictors),
        "aenmd_available": len(available),
        "both_available": total,
        "aenmd_available_fraction": round(float(verdicts["aenmd_available"].mean()), 4)
        if len(verdicts)
        else 0.0,
    }
    if not total:
        return summary
    full = merged["escape_full_rules"].astype(bool)
    guideline = merged["escape_guideline"].astype(bool)
    aenmd = merged["aenmd_escape"].astype(bool)
    return summary | {
        "aenmd_escape": int(aenmd.sum()),
        "aenmd_escape_fraction": round(float(aenmd.mean()), 4),
        "full_rules_vs_aenmd_agree": int((full == aenmd).sum()),
        "full_rules_vs_aenmd_agree_fraction": round(float((full == aenmd).mean()), 4),
        "full_escape_aenmd_decay": int((full & ~aenmd).sum()),
        "full_decay_aenmd_escape": int((~full & aenmd).sum()),
        "guideline_vs_aenmd_agree_fraction": round(float((guideline == aenmd).mean()), 4),
    }
