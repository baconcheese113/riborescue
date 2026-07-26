"""Data contracts — the vocabulary every RiboRescue table is written and validated against.

The enumerations here are the closed sets a column may hold: which stop codons exist, which
consequences are recognised, what triage can conclude, and which evaluation protocols are named. A
table crossing a pipeline boundary is checked against them in `tables.py`, which is where validation
happens; these are what it validates *to*.

The reporter constants are the geometry of the measured library, and every feature window is read
against them.
"""

from enum import StrEnum

from riborescue._version import CONTRACTS_VERSION

__all__ = [
    "CONTRACTS_VERSION",
    "REPORTER_DOWNSTREAM_NT",
    "REPORTER_UPSTREAM_NT",
    "Consequence",
    "EvalConfig",
    "StopCodon",
    "TriageClass",
]


class StopCodon(StrEnum):
    uaa = "UAA"
    uag = "UAG"
    uga = "UGA"


class EvalConfig(StrEnum):
    published_random_cv = "published_random_cv"
    grouped_by_gene = "grouped_by_gene"
    grouped_by_sequence_cluster = "grouped_by_sequence_cluster"
    external_transfer_test = "external_transfer_test"


class Consequence(StrEnum):
    stop_gained = "stop_gained"
    missense_variant = "missense_variant"
    frameshift_variant = "frameshift_variant"
    splice_acceptor_variant = "splice_acceptor_variant"
    splice_donor_variant = "splice_donor_variant"
    splice_region_variant = "splice_region_variant"
    synonymous_variant = "synonymous_variant"
    stop_lost = "stop_lost"
    inframe_deletion = "inframe_deletion"
    inframe_insertion = "inframe_insertion"


class TriageClass(StrEnum):
    supported_nonsense = "supported_nonsense"
    missense = "missense"
    frameshift = "frameshift"
    splice = "splice"
    synonymous = "synonymous"
    stop_loss = "stop_loss"
    unsupported_transcript = "unsupported_transcript"
    unsupported = "unsupported"


REPORTER_UPSTREAM_NT = 72
"""Nucleotides between the start of the reporter context and the premature stop codon.

Every measured variant in the library places the stop at the same offset, so upstream context is a
constant rather than a per-variant quantity.
"""

REPORTER_DOWNSTREAM_NT = (72, 75)
"""Nucleotides following the premature stop codon, one value per oligo design.

The library ships two designs — 147 nt and 150 nt — that differ only downstream. A window reaching
further than the shorter one cannot be filled for every variant.
"""
