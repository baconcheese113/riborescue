from pathlib import Path

import pandas as pd
import pytest
from pandera.errors import SchemaError, SchemaErrors

from riborescue.tables import ReadthroughLabels, TriageInput, TriageOutput, read_table
from riborescue.triage import classify_table

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
