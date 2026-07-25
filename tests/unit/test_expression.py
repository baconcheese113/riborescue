"""TPM, composition and ranking over a hand-worked count table.

Three genes of known length and count, so every TPM below is arithmetic a reader can redo: a gene
twice as long with twice the counts has the same abundance, which is the whole point of dividing by
length before normalising to the library.
"""

import pandas as pd
import pytest

from riborescue.riboseq.expression import (
    composition,
    gene_class,
    gene_symbols,
    library_depth,
    read_counts,
    top_expressed,
    tpm,
)


@pytest.fixture
def counts() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
            "length": [1000, 2000, 1000],
            "lib_a": [100, 200, 700],
            "lib_b": [300, 100, 100],
        }
    )


_SYMBOLS = {"ENSG1": "ACTB", "ENSG2": "MT-CO1", "ENSG3": "RN7SL1"}


def test_a_gene_twice_as_long_with_twice_the_counts_has_the_same_abundance(counts):
    # ENSG1: 100 reads / 1 kb = 100 RPK. ENSG2: 200 / 2 kb = 100 RPK. Equal abundance.
    values = tpm(counts).set_index("gene_id")["lib_a"]
    assert values["ENSG1"] == pytest.approx(values["ENSG2"])


def test_every_library_sums_to_a_million(counts):
    values = tpm(counts)
    assert values["lib_a"].sum() == pytest.approx(1e6)
    assert values["lib_b"].sum() == pytest.approx(1e6)


def test_depth_reports_the_assigned_reads_behind_each_library(counts):
    depth = library_depth(counts)
    assert depth["lib_a"] == 1000
    assert depth["lib_b"] == 500


@pytest.mark.parametrize(
    ("symbol", "expected"),
    [
        ("MT-CO1", "mitochondrial"),
        ("MT-RNR2", "mitochondrial"),
        ("RN7SL1", "structural"),
        ("RNU6-1", "structural"),
        ("RPPH1", "structural"),
        ("ACTB", "nuclear"),
        ("MTOR", "nuclear"),  # starts with MT but is not mitochondrial — the hyphen is the tell
    ],
)
def test_gene_class_separates_the_surviving_contaminants_from_message(symbol, expected):
    assert gene_class(symbol) == expected


def test_composition_says_what_share_of_each_library_is_not_message(counts):
    shares = composition(tpm(counts), _SYMBOLS).set_index("gene_class")
    # lib_a: RPK 100 / 100 / 700 → nuclear 11.1%, mitochondrial 11.1%, structural 77.8%.
    assert shares.loc["structural", "lib_a"] == pytest.approx(77.8, abs=0.1)
    assert shares.loc["nuclear", "lib_a"] == pytest.approx(11.1, abs=0.1)
    assert shares[["lib_a", "lib_b"]].sum().tolist() == pytest.approx([100.0, 100.0])


def test_excluding_the_contaminants_changes_which_gene_leads(counts):
    values = tpm(counts)
    everything = top_expressed(values, _SYMBOLS, top=3)
    assert everything["gene_symbol"].iloc[0] == "RN7SL1"  # the structural RNA dominates lib_a

    nuclear = top_expressed(values, _SYMBOLS, top=3, exclude=("mitochondrial", "structural"))
    assert list(nuclear["gene_symbol"]) == ["ACTB"]  # only one gene is ordinary message


def test_counts_are_read_with_the_library_named_by_its_sample(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text(
        "# Program:featureCounts\n"
        "Geneid\tChr\tStart\tEnd\tStrand\tLength\tresults/x/lib_a.sorted.bam\n"
        "ENSG1\tchr1\t1\t2\t+\t1000\t100\n"
    )
    table = read_counts(path)
    assert list(table.columns) == ["gene_id", "length", "lib_a"]


def test_symbols_come_from_gene_records_only(tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        "#comment\n"
        'chr1\tH\tgene\t1\t9\t.\t+\t.\tgene_id "ENSG1.2"; gene_name "ACTB";\n'
        'chr1\tH\texon\t1\t9\t.\t+\t.\tgene_id "ENSG9.1"; gene_name "NOPE";\n'
    )
    assert gene_symbols(gtf) == {"ENSG1.2": "ACTB"}
