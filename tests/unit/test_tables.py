from pathlib import Path

import pandas as pd
import pytest
from pandera.errors import SchemaError, SchemaErrors

from riborescue.core.tables import (
    PathogenicNonsense,
    ReadthroughLabels,
    TriageInput,
    TriageOutput,
    read_table,
    write_table,
)
from riborescue.variants.triage import classify_table

DATA = Path(__file__).parents[2] / "pipeline/tests/data"


def test_the_committed_label_fixture_validates():
    labels = ReadthroughLabels.validate(read_table(DATA / "labels.tsv"))
    assert len(labels) == 3
    assert bool(labels["censored"].sum()) is True


def test_a_missing_replicate_sd_stays_missing():
    labels = ReadthroughLabels.validate(read_table(DATA / "labels.tsv"))
    assert labels["replicate_sd"].isna().sum() == 1


def test_the_committed_variant_fixture_triages_into_a_valid_output_table():
    triaged = TriageOutput.validate(
        classify_table(TriageInput.validate(read_table(DATA / "variants.tsv")))
    )
    assert triaged["applies"].sum() == 2
    assert set(triaged.columns) >= {"triage_class", "applies", "reason"}


def test_an_efficiency_outside_the_unit_interval_is_refused():
    labels = read_table(DATA / "labels.tsv")
    labels.loc[0, "readthrough_efficiency"] = 1.4
    with pytest.raises((SchemaError, SchemaErrors)):
        ReadthroughLabels.validate(labels, lazy=True)


def test_an_unexpected_column_is_refused():
    variants = read_table(DATA / "variants.tsv").assign(score_1=0.0)
    with pytest.raises((SchemaError, SchemaErrors)):
        TriageInput.validate(variants, lazy=True)


def test_an_unknown_consequence_is_refused():
    variants = pd.DataFrame(
        {"variant_id": ["X:c.1A>T"], "consequence": ["start_lost"], "transcript_supported": [True]}
    )
    with pytest.raises((SchemaError, SchemaErrors)):
        TriageInput.validate(variants, lazy=True)


def test_a_variant_table_survives_a_write_and_read(tmp_path: Path):
    """Chromosomes mix digits and letters, so a column type inferred piecemeal would break here."""

    variants = pd.DataFrame(
        {
            "variant_id": [f"NC_00000{i}:g.100{i}C>T" for i in range(1, 4)],
            "allele_id": [1, 2, 3],
            "chrom": ["7", "22", "X"],
            "pos": [100, 200, 300],
            "ref": ["C", "G", "A"],
            "alt": ["T", "A", "T"],
            "gene_symbol": ["CFTR", "IDUA", "DMD"],
            "gene_id": [1080, 3425, 1756],
            "clinical_significance": ["Pathogenic"] * 3,
            "review_status": ["practice_guideline"] * 3,
            "review_stars": [4, 4, 4],
            "conditions": ["Cystic_fibrosis", "MPS_I", "Muscular_dystrophy"],
        }
    )
    written = tmp_path / "variants.tsv"
    write_table(PathogenicNonsense.validate(variants), written)
    assert PathogenicNonsense.validate(read_table(written)) is not None
