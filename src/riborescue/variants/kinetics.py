"""Whether codon occupancy carries structure the published factorisation cannot express.

ADR-0020. The published design already spans every function of the nine nucleotides a scored variant
carries: `up_123nt` and `down_123nt` are 64-level factors, so a per-codon score is a linear
combination of columns the design matrix already holds, and adding it changes nothing at all. Not
approximately nothing — the column is dropped as rank-deficient and the fit is identical, for any
kinetic quantity whatever.

What a codon score can be is a *structured* hypothesis about which function of that sequence
matters. Two places remain where it is not already spanned. Across the up/down boundary,
`k_up:k_down` is one column no combination of the baseline's reproduces — the rank-1 version of the
saturated 4,096-level interaction ADR-0006 records as untestable naively. And as a replacement for a
64-level factor, where a rate shares strength across codons and a factor knows nothing about a level
it did not see.

The four models below are fitted together and reported together. The saturated interaction is here
as the capacity comparator: it has the information K1 uses and none of its structure, so K1's gain
is read against it rather than against zero.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from riborescue.core.sequences import SYNONYMOUS, as_dna
from riborescue.variants.evaluation import SEED, BootstrapCI, ShuffleKind
from riborescue.variants.readthrough_model import FORMULA, cross_validate, fit_round, r_squared

__all__ = [
    "BASELINE",
    "KINETIC_COLUMNS",
    "MODELS",
    "SUPPORT_STRATA",
    "Criteria",
    "attach",
    "codon_scores",
    "head_to_head",
    "improvement",
    "permute_scores",
    "shuffle",
    "stratified_gain",
    "support",
]

BASELINE = "B"
"""The model every gain is measured against."""

KINETIC_COLUMNS = ("k_up", "k_down")
"""The occupancy of the codon in the P site at termination, and of the first codon past the stop."""

MODELS: dict[str, str] = {
    "B": FORMULA,
    "K1": f"{FORMULA} + stop_type:k_up + k_up:k_down + stop_type:k_up:k_down",
    "K2": "RT_binomial ~ 0 + stop_type + down_123nt + k_up + stop_type:down_123nt",
    "S": f"{FORMULA} + up_123nt:down_123nt",
}
"""The four fits, named before any was run.

`B` is the parity baseline, refit unchanged in every round. `K1` adds seven columns and is the
primary claim. `K2` replaces the upstream factor with its rate, which asks about generalisation
rather than fit. `S` is the capacity comparator and is not expected to be believed.
"""


def codon_scores(table: pd.DataFrame, site: str = "a") -> pd.Series:
    """One occupancy per codon, from a codon table, keyed on the DNA spelling.

    The table carries both ribosomal sites; which one the features are built from is a declared
    choice rather than whichever row came first.
    """

    rows = table.loc[table["site"] == site] if "site" in table.columns else table
    if rows.empty:
        raise ValueError(f"the codon table holds no rows for the {site} site")
    scores = rows.set_index(rows["codon"].map(as_dna))["occupancy"]
    if scores.index.has_duplicates:
        raise ValueError("the codon table names a codon twice")
    return scores.rename("occupancy")


def attach(features: pd.DataFrame, scores: pd.Series) -> pd.DataFrame:
    """Add the kinetic columns to a feature table, refusing a triplet the table cannot score.

    A missing score would otherwise propagate as a silent NaN into a design matrix, where it removes
    a row from one model and not from another and makes the comparison between them meaningless.
    """

    attached = features.copy()
    for column, triplet in zip(KINETIC_COLUMNS, ("up_123nt", "down_123nt"), strict=True):
        codons = features[triplet].map(as_dna)
        attached[column] = codons.map(scores).astype(float)
        if (missing := sorted(set(codons[attached[column].isna()]))) != []:
            raise ValueError(f"the codon table scores no {', '.join(missing)} for {triplet}")
    return attached


def permute_scores(scores: pd.Series, kind: ShuffleKind, seed: int) -> pd.Series:
    """Permute the codon-to-score assignment, leaving every sequence column untouched.

    `global_` permutes the whole table, breaking any dependence on the measured scores at all.
    `context_matched` permutes within synonymous families, which preserves the amino acid, the
    chemistry and the sequence context and destroys only the quantity the kinetic hypothesis names.
    It is the decisive one, and the single-codon families are invariant under it by construction.
    """

    generator = np.random.default_rng(seed)
    if kind is ShuffleKind.global_:
        return pd.Series(
            generator.permutation(scores.to_numpy()), index=scores.index, name=scores.name
        )
    if kind is ShuffleKind.context_matched:
        permuted = scores.copy()
        for family in SYNONYMOUS.values():
            present = [codon for codon in family if codon in scores.index]
            permuted.loc[present] = generator.permutation(scores.loc[present].to_numpy())
        return permuted
    raise ValueError(f"{kind} permutes variants, not the codon table; use `shuffle`")


def shuffle(
    features: pd.DataFrame, scores: pd.Series, kind: ShuffleKind | None, seed: int
) -> pd.DataFrame:
    """A feature table with its kinetic columns broken in one named way, or not at all.

    `within_gene` permutes the kinetic columns among the variants of a gene, so gene identity and
    every sequence column survive and only which variant carries which kinetics does not. The two
    columns move together as a pair, because permuting them apart would break their joint
    distribution as well and the control would test two things at once.
    """

    if kind is None:
        return attach(features, scores)
    if kind is not ShuffleKind.within_gene:
        return attach(features, permute_scores(scores, kind, seed))

    attached = attach(features, scores)
    generator = np.random.default_rng(seed)
    columns = list(KINETIC_COLUMNS)
    for _, rows in attached.groupby("gene", sort=True):
        if len(rows) > 1:
            order = generator.permutation(len(rows))
            attached.loc[rows.index, columns] = rows[columns].to_numpy()[order]
    return attached


def head_to_head(
    features: pd.DataFrame,
    rounds: pd.DataFrame,
    scores: pd.Series,
    *,
    models: dict[str, str] | None = None,
    kind: ShuffleKind | None = None,
    seed: int = SEED,
) -> pd.DataFrame:
    """Every model over every round of one split, one row each.

    All four fit the same rows in the same order, so a difference between them is the formula and
    nothing else. `kind` breaks the kinetic columns in one named way; the baseline and the saturated
    comparator do not read those columns, so they are identical under every shuffle — which is what
    makes a shuffled gain readable as a gain.
    """

    prepared = shuffle(features, scores, kind, seed)
    return pd.DataFrame(
        {
            "model": name,
            "round": prediction.round,
            "r2": prediction.r2,
            "held_out": len(prediction.predictions),
        }
        for name, formula in (models or MODELS).items()
        for prediction in cross_validate(prepared, rounds, formula)
    )


def improvement(scored: pd.DataFrame, baseline: str = BASELINE) -> pd.DataFrame:
    """Each model's gain over the baseline, round by round.

    Paired within round, because the rounds hold out different rows and an unpaired difference of
    means would carry that variation instead of the difference between the models.
    """

    wide = scored.pivot(index="round", columns="model", values="r2")
    if baseline not in wide.columns:
        raise ValueError(f"no {baseline!r} rows to measure against")
    gains = wide.drop(columns=baseline).sub(wide[baseline], axis=0)
    return gains.reset_index().melt(id_vars="round", var_name="model", value_name="gain")


SUPPORT_STRATA = ("unsupported", "thin", "supported")
"""How much training data a held-out variant's interaction cell had.

`unsupported` is a `(stop_type, down_123nt)` cell the training fold never observed. Its column is
dependent, `identifiable_columns` drops it, and the prediction falls back to the marginal — so a
gain concentrated there is a statement about what the baseline cannot fit rather than about what
kinetics knows. The rest split at the median of the observed cells.
"""


def support(features: pd.DataFrame, train_rows: pd.Series) -> pd.Series:
    """Which support stratum each row falls in, given what its training fold contained.

    Computed from the training fold and applied to every row, so the stratum a held-out variant sits
    in is a property of that round rather than of the library as a whole.
    """

    train = features[features.index.isin(train_rows)]
    observed = train.groupby(["stop_type", "down_123nt"]).size()
    counts = pd.Series(
        [
            int(observed.get((stop, down), 0))
            for stop, down in zip(features["stop_type"], features["down_123nt"], strict=True)
        ],
        index=features.index,
    )
    median = float(counts[counts > 0].median())
    values = counts.to_numpy()
    strata = np.where(values == 0, "unsupported", np.where(values < median, "thin", "supported"))
    return pd.Series(strata, index=features.index, name="support")


def stratified_gain(
    features: pd.DataFrame,
    rounds: pd.DataFrame,
    scores: pd.Series,
    *,
    models: dict[str, str] | None = None,
    kind: ShuffleKind | None = None,
    seed: int = SEED,
) -> pd.DataFrame:
    """Each model's held-out performance within each support stratum, round by round.

    A stratum with fewer than three held-out rows in a round has no stable correlation and is
    dropped from that round rather than reported as a number built on two points.
    """

    prepared = shuffle(features, scores, kind, seed)
    rows: list[dict[str, object]] = []
    for name, formula in (models or MODELS).items():
        for round_ in sorted(int(r) for r in rounds["round"].unique()):
            train = rounds.loc[rounds["round"] == round_, "row"]
            prediction = fit_round(prepared, train, round_, formula)
            strata = support(prepared, train).reindex(prediction.predictions.index)
            for stratum in SUPPORT_STRATA:
                inside = strata == stratum
                if int(inside.sum()) < 3:
                    continue
                rows.append(
                    {
                        "model": name,
                        "round": round_,
                        "support": stratum,
                        "r2": r_squared(
                            prediction.predictions[inside], prediction.observed[inside]
                        ),
                        "held_out": int(inside.sum()),
                    }
                )
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class Criteria:
    """ADR-0020's passing rule, evaluated. Every reporting path goes through this.

    The rule is a conjunction and is written as one, so that a claim cannot be made from two of
    three conditions by a reader who looked at the first two. `supported` is the only thing entitled
    to be described as kinetics carrying transferable information.
    """

    grouped_gain: BootstrapCI
    shuffles: dict[str, BootstrapCI]
    confined_to_unsupported: bool

    @property
    def gain_excludes_zero(self) -> bool:
        """Condition 1: the gain under the grouped-by-gene split is positive and clears zero."""

        return self.grouped_gain.point > 0 and not self.grouped_gain.includes_zero()

    @property
    def shuffles_collapse(self) -> bool:
        """Condition 2: every shuffle control's interval includes zero."""

        return bool(self.shuffles) and all(ci.includes_zero() for ci in self.shuffles.values())

    @property
    def survives_support(self) -> bool:
        """Condition 3: the gain is not confined to cells the baseline could not fit."""

        return not self.confined_to_unsupported

    @property
    def supported(self) -> bool:
        return self.gain_excludes_zero and self.shuffles_collapse and self.survives_support

    def failures(self) -> tuple[str, ...]:
        """Which conditions did not hold, named. Empty when the claim is supported."""

        return tuple(
            name
            for name, held in (
                ("gain_excludes_zero", self.gain_excludes_zero),
                ("shuffles_collapse", self.shuffles_collapse),
                ("survives_support", self.survives_support),
            )
            if not held
        )
