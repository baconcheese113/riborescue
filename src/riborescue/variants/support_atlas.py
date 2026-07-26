"""What the model had actually seen when it scored a variant, term by term.

The published model is factorised, so it composes a context it has never met out of marginals it has
estimated separately. That is legitimate — an unseen complete context is still predictable when
every contributing factor level is supported — and it is exactly why a single "was this context
observed?" flag is the wrong thing to report. Support is a property of *terms*, not of contexts, and
a variant can be well supported on all of them while its exact context was never measured.

So each scored variant carries a count per term: how many library variants shared its stop codon,
its upstream triplet, its downstream triplet, the stop-by-downstream cell the model interacts, and
its complete context. A cell with no observations at all yields a dependent design column, which
`identifiable_columns` drops and the prediction treats as zero effect, falling back to the marginal;
that variant is flagged separately, because its interval is narrower than its evidence.

**None of this is disease evidence.** The library is a synthetic reporter — native sequence around a
premature stop in an EGFP-T2A-mCherry construct, in HEK293T — so a supporting observation is an
assay measurement of a variant sharing model terms, not a clinical observation about a patient's
gene. `context_analogues` names its output accordingly and the naming is load-bearing.
"""

from __future__ import annotations

import pandas as pd

__all__ = [
    "SUPPORT_TERMS",
    "context_analogues",
    "support_atlas",
    "term_support",
]

SUPPORT_TERMS = ("stop_type", "up_123nt", "down_123nt")
"""The factors the published model estimates. Support is counted for each, and for their products
where the model interacts them."""

_KEYS = {
    "stop": ("stop_type",),
    "upstream": ("up_123nt",),
    "downstream": ("down_123nt",),
    "interaction_cell": ("stop_type", "down_123nt"),
    "complete_context": ("stop_type", "up_123nt", "down_123nt"),
}
"""What each support column counts. `interaction_cell` is the only product the model fits, and
`complete_context` is the exact nine nucleotides — reported because a reader will ask, not because
the model needs it."""


def _normalise(frame: pd.DataFrame, columns=SUPPORT_TERMS) -> pd.DataFrame:
    """Triplets are spelled as RNA in the library and as DNA in the scored table, or the reverse."""

    return frame.assign(
        **{
            column: frame[column].astype(str).str.upper().str.replace("U", "T")
            for column in columns
        }
    )


def term_support(library: pd.DataFrame) -> dict[str, pd.Series]:
    """How many measured variants back each level of each term, counted once from the library."""

    normalised = _normalise(library)
    return {name: normalised.groupby(list(keys)).size() for name, keys in _KEYS.items()}


def support_atlas(scored: pd.DataFrame, library: pd.DataFrame) -> pd.DataFrame:
    """Every scored variant, with the measured support behind each term the model used.

    `aliased_interaction` marks the variants whose stop-by-downstream cell the library never
    observed. Their prediction is the marginal, which is a weaker statement than the interval alone
    conveys, and it is the one support failure that changes what the model actually computed rather
    than only how confident a reader should be.
    """

    counts = term_support(library)
    normalised = _normalise(scored)
    atlas = (
        scored[["variant_id"]].copy()
        if "variant_id" in scored.columns
        else pd.DataFrame(index=scored.index)
    )
    for name, keys in _KEYS.items():
        index = (
            pd.MultiIndex.from_frame(normalised[list(keys)])
            if len(keys) > 1
            else pd.Index(normalised[keys[0]])
        )
        atlas[f"support_{name}"] = counts[name].reindex(index).fillna(0).to_numpy().astype(int)
    atlas["aliased_interaction"] = atlas["support_interaction_cell"] == 0
    atlas["measured_exactly"] = atlas["support_complete_context"] > 0
    # The weakest term is what a reader should look at first: a prediction is composed of all of
    # them, so it is no better supported than the thinnest one it rests on.
    atlas["weakest_term"] = atlas[
        ["support_stop", "support_upstream", "support_downstream", "support_interaction_cell"]
    ].min(axis=1)
    return atlas


def context_analogues(
    scored: pd.DataFrame,
    library: pd.DataFrame,
    *,
    label: str,
    limit: int = 10,
) -> pd.DataFrame:
    """The measured reporter variants sharing a scored variant's complete context.

    One row per (scored variant, measured variant) pair, carrying the measurement itself. This is
    the strongest evidence the library can offer about a scored variant, and it is still an assay
    reading: `source` records which measurement it is, and no column of this table is a clinical
    observation about the gene the scored variant sits in.

    `limit` caps how many analogues a variant carries, because a common context has hundreds and a
    reader needs the evidence rather than the census. `support_complete_context` in the atlas is the
    uncapped count, so a truncated list is visible as one.
    """

    keys = list(_KEYS["complete_context"])
    measured = _normalise(library).assign(source=label)
    columns = ["gene", "RT_binomial", "source", *keys]
    left = _normalise(scored)[["variant_id", "gene_symbol", *keys]]
    paired = left.merge(measured[columns], on=keys, how="inner", suffixes=("", "_measured"))
    paired = paired.rename(columns={"gene": "analogue_gene", "RT_binomial": "analogue_readthrough"})
    # Same gene first — a reporter measurement of a different variant in the same gene is the
    # closest thing the library has — then by measured readthrough so the list is stable.
    paired["same_gene"] = paired["analogue_gene"].str.upper() == paired["gene_symbol"].str.upper()
    paired = paired.sort_values(
        ["variant_id", "same_gene", "analogue_readthrough"], ascending=[True, False, False]
    )
    return paired.groupby("variant_id", sort=False).head(limit).reset_index(drop=True)
