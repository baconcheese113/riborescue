"""The front door: decide whether readthrough applies to a variant before any analysis runs.

Most submitted variants are not supported nonsense mutations, so triage is a first-class step that
classifies the variant and says, in plain language, why readthrough does or does not apply.
"""

from dataclasses import dataclass

import pandas as pd

from riborescue.contracts import Consequence, TriageClass

__all__ = ["TriageResult", "classify", "classify_table"]


@dataclass(frozen=True)
class TriageResult:
    triage_class: TriageClass
    applies: bool
    reason: str


_BY_CONSEQUENCE: dict[Consequence, TriageResult] = {
    Consequence.stop_gained: TriageResult(
        TriageClass.supported_nonsense,
        True,
        "a premature stop codon truncates the protein; "
        "readthrough may restore the full-length product",
    ),
    Consequence.missense_variant: TriageResult(
        TriageClass.missense,
        False,
        "an amino-acid substitution, not a premature stop; readthrough has nothing to read through",
    ),
    Consequence.frameshift_variant: TriageResult(
        TriageClass.frameshift,
        False,
        "a shifted reading frame changes the downstream sequence; readthrough cannot restore it",
    ),
    Consequence.splice_acceptor_variant: TriageResult(
        TriageClass.splice,
        False,
        "a splice-site change alters the transcript itself, which readthrough does not address",
    ),
    Consequence.splice_donor_variant: TriageResult(
        TriageClass.splice,
        False,
        "a splice-site change alters the transcript itself, which readthrough does not address",
    ),
    Consequence.splice_region_variant: TriageResult(
        TriageClass.splice,
        False,
        "a splice-region change may alter the transcript, which readthrough does not address",
    ),
    Consequence.synonymous_variant: TriageResult(
        TriageClass.synonymous,
        False,
        "a silent change leaves the protein unchanged; there is nothing to rescue",
    ),
    Consequence.stop_lost: TriageResult(
        TriageClass.stop_loss,
        False,
        "the native stop is lost, extending the protein; readthrough would not help",
    ),
}


def classify(consequence: Consequence, *, transcript_supported: bool = True) -> TriageResult:
    if not transcript_supported:
        return TriageResult(
            TriageClass.unsupported_transcript,
            False,
            "the variant's transcript is outside the supported set, so no analysis is offered",
        )
    result = _BY_CONSEQUENCE.get(consequence)
    if result is None:
        return TriageResult(
            TriageClass.unsupported,
            False,
            f"{consequence.value} is outside the variant classes RiboRescue supports",
        )
    return result


def classify_table(variants: pd.DataFrame) -> pd.DataFrame:
    """Triage a table of variants, appending the verdict, its class and its reason to each row."""

    verdicts = [
        classify(Consequence(row.consequence), transcript_supported=bool(row.transcript_supported))
        for row in variants.itertuples()
    ]
    return variants.assign(
        triage_class=[v.triage_class.value for v in verdicts],
        applies=[v.applies for v in verdicts],
        reason=[v.reason for v in verdicts],
    )
