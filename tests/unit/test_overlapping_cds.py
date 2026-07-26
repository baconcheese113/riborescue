"""The neighbouring-gene filter, on an annotation small enough to reason about.

This filter was silently absent for the whole project: it searched for a feature type called
`three_prime_utr`, GENCODE emits `UTR` and leaves the end to be worked out from strand and coding
extent, and so it returned an empty set every time. Nothing noticed, because the other filters were
excluding most of the same transcripts for other reasons.

The production comparison in `tests/controls/test_result_drift.py` would catch a repeat, but only
where the pipeline has been run — a 1.5 GB annotation and a count table are not in the repository,
so it skips on a fresh clone. This does not. Six transcripts in a string, every case the real
annotation contains, and it fails on a fresh clone and in continuous integration alike.
"""

import gzip
from pathlib import Path

import pytest

from riborescue.riboseq.readthrough_assay import overlapping_downstream_cds

# GENCODE columns, and GENCODE's undifferentiated `UTR` rather than a 3'-specific feature type —
# using the real spelling is the whole point, since assuming the other one is the bug.
_ROWS = [
    # A plus-strand transcript whose 3' UTR runs into gene B's coding sequence. Excluded.
    ("chr1", "CDS", 1000, 2000, "+", "GA", "TA"),
    ("chr1", "UTR", 2001, 3000, "+", "GA", "TA"),
    ("chr1", "CDS", 2500, 3500, "+", "GB", "TB"),
    # A plus-strand transcript whose *5'* UTR runs into gene D. Kept: a ribosome reading gene D
    # lands before this transcript's start codon, nowhere near its readthrough window.
    ("chr1", "CDS", 5000, 6000, "+", "GC", "TC"),
    ("chr1", "UTR", 4000, 4999, "+", "GC", "TC"),
    ("chr1", "CDS", 4200, 4800, "+", "GD", "TD"),
    # A minus-strand transcript whose 3' UTR is *below* its coding sequence, running into gene F.
    # Excluded, and the case a plus-strand-only rule would miss.
    ("chr1", "CDS", 10000, 11000, "-", "GE", "TE"),
    ("chr1", "UTR", 9000, 9999, "-", "GE", "TE"),
    ("chr1", "CDS", 9500, 9800, "-", "GF", "TF"),
    # A transcript whose 3' UTR overlaps another isoform of its own gene. Kept: that is ordinary
    # isoform structure, not a neighbouring gene's ribosomes.
    ("chr1", "CDS", 20000, 21000, "+", "GG", "TG"),
    ("chr1", "UTR", 21001, 22000, "+", "GG", "TG"),
    ("chr1", "CDS", 21500, 21800, "+", "GG", "TG2"),
    # A transcript with a neighbour that does not reach it. Kept.
    ("chr1", "CDS", 30000, 31000, "+", "GH", "TH"),
    ("chr1", "UTR", 31001, 32000, "+", "GH", "TH"),
    ("chr1", "CDS", 33000, 34000, "+", "GI", "TI"),
]

_EXCLUDED = {"TA", "TE"}
"""The two whose readthrough window a different gene's ribosomes occupy."""


@pytest.fixture(scope="module")
def annotation(tmp_path_factory) -> Path:
    path = tmp_path_factory.mktemp("gtf") / "tiny.gtf.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("#!genome-build TINY\n")
        for chrom, feature, start, end, strand, gene, transcript in _ROWS:
            attributes = f'gene_id "{gene}"; transcript_id "{transcript}";'
            handle.write(
                f"{chrom}\tTEST\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n"
            )
    return path


def test_the_filter_excludes_something_at_all(annotation):
    # The failure this file exists for: an empty set every time, from looking for a feature type
    # the annotation does not use.
    assert overlapping_downstream_cds(annotation) != frozenset()


def test_exactly_the_contaminated_transcripts_are_excluded(annotation):
    assert set(overlapping_downstream_cds(annotation)) == _EXCLUDED


def test_a_five_prime_utr_is_not_mistaken_for_a_three_prime_one(annotation):
    # TC's UTR overlaps gene D, but it is upstream of the coding sequence and irrelevant here.
    assert "TC" not in overlapping_downstream_cds(annotation)


def test_the_minus_strand_case_is_found(annotation):
    # On the minus strand the 3' UTR has the lower coordinates, so a rule written for the plus
    # strand alone silently keeps every contaminated minus-strand transcript.
    assert "TE" in overlapping_downstream_cds(annotation)


def test_an_overlap_with_the_transcripts_own_gene_does_not_count(annotation):
    assert "TG" not in overlapping_downstream_cds(annotation)


def test_a_neighbour_that_does_not_reach_the_window_does_not_count(annotation):
    assert "TH" not in overlapping_downstream_cds(annotation)


def test_an_annotation_with_no_coding_sequence_yields_nothing(tmp_path):
    path = tmp_path / "empty.gtf.gz"
    with gzip.open(path, "wt") as handle:
        handle.write('chr1\tTEST\tUTR\t1\t100\t.\t+\t.\tgene_id "G"; transcript_id "T";\n')
    assert overlapping_downstream_cds(path) == frozenset()
