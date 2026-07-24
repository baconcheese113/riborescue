import pandas as pd

from riborescue.variants.diseases import normalize_conditions, parse_conditions


def test_a_fully_cross_referenced_condition_is_mapped():
    [condition] = parse_conditions(
        "Cystic_fibrosis", "MedGen:C0010674,OMIM:219700,Orphanet:586,MONDO:MONDO:0009061"
    )
    assert condition.name == "Cystic_fibrosis"
    assert condition.medgen == "C0010674"
    assert condition.omim == "219700"
    assert condition.orphanet == "586"
    assert condition.mondo == "MONDO:0009061"
    assert condition.mapping_status == "mapped"
    assert condition.reason == ""


def test_a_medgen_only_condition_is_labelled_not_dropped():
    [condition] = parse_conditions("Muscular_dystrophy", "MedGen:C0026850")
    assert condition.mapping_status == "medgen_only"
    assert "no OMIM or Orphanet" in condition.reason


def test_a_placeholder_condition_is_recognised_as_not_a_disease():
    [not_provided] = parse_conditions("not_provided", "MedGen:CN517202")
    [not_specified] = parse_conditions("not_specified", "MedGen:CN169374")
    assert not_provided.mapping_status == "placeholder"
    assert not_specified.mapping_status == "placeholder"
    assert "not provided" in not_provided.reason


def test_the_plain_cui_not_provided_placeholder_is_also_caught():
    # C3661900 is a regular MedGen CUI, not a CN* id, but it still means "not provided".
    [condition] = parse_conditions("not_provided", "MedGen:C3661900")
    assert condition.mapping_status == "placeholder"


def test_a_condition_with_no_medgen_is_unmapped():
    [condition] = parse_conditions("Some_rare_thing", "OMIM:123456")
    assert condition.medgen == ""
    assert condition.mapping_status == "unmapped"


def test_parallel_lists_stay_index_aligned():
    conditions = parse_conditions(
        "Cystic_fibrosis|not_specified",
        "MedGen:C0010674,OMIM:219700,Orphanet:586|MedGen:CN169374",
    )
    assert [c.name for c in conditions] == ["Cystic_fibrosis", "not_specified"]
    assert conditions[0].mapping_status == "mapped"
    assert conditions[1].mapping_status == "placeholder"


def test_repeated_omim_ids_are_all_kept():
    [condition] = parse_conditions(
        "Retinitis_pigmentosa", "MedGen:C0035334,OMIM:268000,OMIM:PS268000"
    )
    assert condition.omim == "268000;PS268000"


def test_more_names_than_xrefs_leaves_the_extra_unmapped():
    conditions = parse_conditions("A|B", "MedGen:C0000001,OMIM:1")
    assert conditions[1].name == "B"
    assert conditions[1].mapping_status == "unmapped"


def test_normalize_preserves_one_variant_to_many_conditions():
    variants = pd.DataFrame(
        [
            {
                "variant_id": "V1",
                "gene_symbol": "CFTR",
                "conditions": "Cystic_fibrosis|Bronchiectasis",
                "condition_xrefs": "MedGen:C0010674,OMIM:219700,Orphanet:586|MedGen:C0006267",
            }
        ]
    )
    rows = normalize_conditions(variants)
    assert list(rows["condition_name"]) == ["Cystic_fibrosis", "Bronchiectasis"]
    assert list(rows["medgen"]) == ["C0010674", "C0006267"]
    assert set(rows["variant_id"]) == {"V1"}


def test_a_variant_with_no_conditions_contributes_no_rows():
    variants = pd.DataFrame(
        [{"variant_id": "V1", "gene_symbol": "X", "conditions": "", "condition_xrefs": ""}]
    )
    assert normalize_conditions(variants).empty


def test_the_normalized_columns_are_stable_even_when_empty():
    variants = pd.DataFrame(
        [{"variant_id": "V1", "gene_symbol": "X", "conditions": "", "condition_xrefs": ""}]
    )
    rows = normalize_conditions(variants)
    assert list(rows.columns) == [
        "variant_id",
        "gene_symbol",
        "condition_name",
        "medgen",
        "omim",
        "orphanet",
        "mondo",
        "mesh",
        "mapping_status",
        "reason",
    ]
