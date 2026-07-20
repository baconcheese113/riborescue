"""The sequence context around a premature stop, in the form the readthrough model expects.

Placing a genomic variant on a transcript is where a pipeline goes quietly wrong, so every context
is earned rather than assumed. Three facts must hold before a variant is scored: its position falls
in an exon of its gene's MANE Select transcript, the reference base ClinVar reports matches the
transcript's own base there, and substituting the alternate base turns that codon into a stop. A
variant failing any of them carries a reason and no context, never a context computed anyway.

Features are written as the readthrough assay wrote them — RNA, lower case, the triplet before the
stop and the triplet after it — so a ClinVar variant and a measured library variant are described
identically and the model cannot tell them apart by formatting.
"""

from dataclasses import dataclass
from enum import StrEnum
from typing import cast

import pandas as pd
from Bio.Seq import Seq

from riborescue.contracts import REPORTER_DOWNSTREAM_NT
from riborescue.transcripts import STOP_CODONS, TranscriptModel

__all__ = [
    "ContextFailure",
    "ProteinDisagreement",
    "PtcContext",
    "context_for",
    "contexts_for",
    "disagreements_with_protein",
]


def _int(value: object) -> int:
    """Pandas hands back a scalar of its own; the callers here all want a plain integer."""

    return int(cast(int, value))


_TRANSCRIBE = str.maketrans("ACGT", "acgu")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_WINDOW = min(REPORTER_DOWNSTREAM_NT)


class ContextFailure(StrEnum):
    """Why a variant has no scoreable context. Each is a fact about the data, not a tolerance."""

    no_mane_transcript = "no_mane_transcript"
    outside_transcript = "outside_transcript"
    reference_mismatch = "reference_mismatch"
    not_a_premature_stop = "not_a_premature_stop"
    truncated_context = "truncated_context"


@dataclass(frozen=True)
class ProteinDisagreement:
    """A variant whose reference codon does not translate to the residue the protein carries."""

    transcript_id: str
    protein_position: int
    codon: str
    translated: str
    in_protein: str


@dataclass(frozen=True)
class PtcContext:
    """A premature stop codon and the sequence either side of it."""

    transcript_id: str
    protein_position: int
    stop_type: str
    reference_codon: str
    up_123nt: str
    down_123nt: str
    upstream: str
    downstream: str


def _transcribe(sequence: str) -> str:
    return sequence.translate(_TRANSCRIBE)


def context_for(
    model: TranscriptModel, position: int, ref: str, alt: str
) -> PtcContext | ContextFailure:
    """Place a substitution on a transcript and read the context around the stop it creates."""

    offset = model.offset_of(position)
    if offset is None:
        return ContextFailure.outside_transcript

    # A VCF states alleles on the forward strand; a reverse-strand transcript reads the
    # complements.
    expected, substituted = (
        (ref, alt)
        if model.strand == "+"
        else (
            ref.translate(_COMPLEMENT),
            alt.translate(_COMPLEMENT),
        )
    )
    if model.sequence[offset] != expected:
        return ContextFailure.reference_mismatch

    covering = model.codon_covering(position)
    if covering is None:
        return ContextFailure.outside_transcript
    codon_start, codon = covering
    mutated = codon[: offset - codon_start] + substituted + codon[offset - codon_start + 1 :]
    if mutated not in STOP_CODONS or codon in STOP_CODONS:
        return ContextFailure.not_a_premature_stop

    upstream = model.sequence[codon_start - _WINDOW : codon_start]
    downstream = model.sequence[codon_start + 3 : codon_start + 3 + _WINDOW]
    if len(upstream) < _WINDOW or len(downstream) < _WINDOW:
        return ContextFailure.truncated_context

    coding = model.coding_offset
    assert coding is not None
    return PtcContext(
        transcript_id=model.transcript_id,
        protein_position=(codon_start - coding) // 3 + 1,
        stop_type=_transcribe(mutated),
        reference_codon=codon,
        up_123nt=_transcribe(upstream[-3:]),
        down_123nt=_transcribe(downstream[:3]),
        upstream=_transcribe(upstream),
        downstream=_transcribe(downstream),
    )


def contexts_for(variants: pd.DataFrame, models: dict[int, TranscriptModel]) -> pd.DataFrame:
    """Read the context of every variant, keeping the failures and their reasons alongside."""

    rows = []
    for variant in variants.itertuples():
        model = models.get(_int(variant.gene_id))
        found = (
            ContextFailure.no_mane_transcript
            if model is None
            else context_for(model, _int(variant.pos), str(variant.ref), str(variant.alt))
        )
        row: dict[str, object] = {
            "variant_id": variant.variant_id,
            "gene_id": _int(variant.gene_id),
            "gene_symbol": variant.gene_symbol,
            "review_stars": variant.review_stars,
        }
        if isinstance(found, ContextFailure):
            rows.append(row | {"scoreable": False, "reason": found.value})
            continue
        rows.append(
            row
            | {
                "scoreable": True,
                "reason": "",
                "transcript_id": found.transcript_id,
                "protein_position": found.protein_position,
                "stop_type": found.stop_type,
                "up_123nt": found.up_123nt,
                "down_123nt": found.down_123nt,
                "upstream": found.upstream,
                "downstream": found.downstream,
            }
        )
    return pd.DataFrame(rows)


def disagreements_with_protein(
    variants: pd.DataFrame,
    models: dict[int, TranscriptModel],
    proteins: dict[str, str],
) -> list[ProteinDisagreement]:
    """Check every placement against the reference protein and return the ones that disagree.

    A codon read at the right position in the right frame on the right strand translates to the
    residue the protein already has there. Checking that against an independently distributed
    protein sequence catches an error in any of the three at once, which no amount of internal
    consistency would reveal.
    """

    found: list[ProteinDisagreement] = []
    for variant in variants.itertuples():
        model = models.get(_int(variant.gene_id))
        if model is None:
            continue
        context = context_for(model, _int(variant.pos), str(variant.ref), str(variant.alt))
        protein = proteins.get(model.protein_id)
        if isinstance(context, ContextFailure) or protein is None:
            continue
        if context.protein_position > len(protein):
            continue
        translated = str(Seq(context.reference_codon).translate())
        carried = protein[context.protein_position - 1]
        if translated != carried:
            found.append(
                ProteinDisagreement(
                    transcript_id=model.transcript_id,
                    protein_position=context.protein_position,
                    codon=context.reference_codon,
                    translated=translated,
                    in_protein=carried,
                )
            )
    return found
