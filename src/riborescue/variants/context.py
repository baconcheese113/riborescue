"""The sequence context around a premature stop, in the form the readthrough model expects.

Three facts must hold before a variant is scored: the position falls in an exon of its gene's MANE
Select transcript, the reference base matches the transcript's own base there, and substituting the
alternate turns that codon into a stop. Failing any of them yields a reason and no context.

Features are written as the assay wrote them — RNA, lower case, the triplet either side — so a
ClinVar variant and a measured one are described identically.
"""

from dataclasses import dataclass
from enum import StrEnum
from itertools import pairwise
from typing import cast

import pandas as pd
from Bio.Seq import Seq

from riborescue.core.contracts import REPORTER_DOWNSTREAM_NT
from riborescue.variants.transcripts import STOP_CODONS, TranscriptModel

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
    original_aa: str
    up_123nt: str
    down_123nt: str
    upstream: str
    downstream: str
    exon_count: int
    in_last_exon: bool
    nt_to_last_junction: int | None
    nt_from_start: int
    ptc_exon_length: int


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
    junctions = model.junction_offsets
    last_junction = junctions[-1] if junctions else None
    # The exon boundaries in transcript coordinates bracket the stop; the exon it lands in is the
    # one interval that contains its codon start, and its length drives the long-exon NMD rule.
    boundaries = (0, *junctions, len(model.sequence))
    ptc_exon_length = next(
        high - low for low, high in pairwise(boundaries) if low <= codon_start < high
    )
    return PtcContext(
        transcript_id=model.transcript_id,
        protein_position=(codon_start - coding) // 3 + 1,
        stop_type=_transcribe(mutated),
        reference_codon=codon,
        original_aa=str(Seq(codon).translate()),
        up_123nt=_transcribe(upstream[-3:]),
        down_123nt=_transcribe(downstream[:3]),
        upstream=_transcribe(upstream),
        downstream=_transcribe(downstream),
        exon_count=len(model.exons),
        in_last_exon=last_junction is None or codon_start >= last_junction,
        nt_to_last_junction=None if last_junction is None else last_junction - codon_start,
        nt_from_start=codon_start - coding,
        ptc_exon_length=ptc_exon_length,
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
                "original_aa": found.original_aa,
                "up_123nt": found.up_123nt,
                "down_123nt": found.down_123nt,
                "upstream": found.upstream,
                "downstream": found.downstream,
                "exon_count": found.exon_count,
                "in_last_exon": found.in_last_exon,
                "nt_to_last_junction": found.nt_to_last_junction,
                "nt_from_start": found.nt_from_start,
                "ptc_exon_length": found.ptc_exon_length,
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
