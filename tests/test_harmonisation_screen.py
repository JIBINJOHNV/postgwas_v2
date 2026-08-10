"""Tests for dataset-level user-facing harmonisation output."""

from postgwas.modules.harmonisation.service import (
    _study_decisions_block,
    _validated_column_mapping_block,
)


def test_study_decisions_are_shown_as_aligned_user_facing_values():
    text = _study_decisions_block({
        "effect_type": "odds_ratio",
        "se_scale": "log_odds",
        "pvalue_type": "raw",
        "eaf_is_maf": False,
        "strand": "forward",
    })

    assert text.splitlines() == [
        "    🧠  Study-wide decisions",
        "        🔬  Effect type      : odds_ratio",
        "        🔬  SE scale         : log_odds",
        "        🔬  P-value type     : raw",
        "        🧬  Frequency type   : effect allele frequency",
        "        🧬  Strand consensus : forward",
    ]


def test_maf_like_study_decision_is_displayed_as_provisional():
    text = _study_decisions_block({
        "effect_type": "beta",
        "pvalue_type": "raw",
        "eaf_is_maf": True,
        "strand": "forward",
    })

    assert "MAF-like; reference confirmation required" in text
    assert "SE scale         : not applicable to beta" in text


def test_validated_column_mapping_uses_public_sample_sheet_fields():
    columns = {
        "gwas_outputname": "ADHD2022", "sumstat_file": "/data/adhd.tsv.gz",
        "trait_type": "case_control", "delimiter": "whitespace",
        "chr_col": "CHR", "pos_col": "BP", "chr_pos_col": "NA",
        "snp_id_col": "SNP", "ea_col": "A1", "oa_col": "A2",
        "eaf_col": "FRQ_U", "beta_or_col": "OR", "se_col": "SE",
        "imp_z_col": "NA", "pval_col": "P", "ncontrol_col": "Nco",
        "ncase_col": "Nca", "imp_info_col": "INFO",
        "info_source_detail": "internal column INFO", "fixed_info": None,
        "declared_effect_type": "odds_ratio",
        "declared_pvalue_type": "raw",
        "ncontrol": "NA", "ncase": "NA",
        "infofile": "NA", "infocolumn": "NA",
        "eaffile": "NA", "eafcolumn": "NA",
    }
    decisions = {
        "effect_type": "odds_ratio", "pvalue_type": "raw", "eaf_is_maf": False,
    }

    text = _validated_column_mapping_block(columns, decisions)

    assert "    🔬  Validated sample-sheet values" in text
    assert "        🔹  Dataset" in text
    assert "        🧬  Variant columns" in text
    assert "        🔬  Association columns" in text
    assert "        🧮  Sample size" in text
    assert "        🔬  Quality and external sources" in text
    for field in (
        "dataset_id", "input_file", "trait_type", "chromosome_column",
        "position_column", "chromosome_position_column", "variant_id_column",
        "effect_allele_column", "other_allele_column",
        "effect_allele_frequency_column", "effect_column", "effect_type",
        "standard_error_column", "z_score_column", "p_value_column",
        "p_value_type", "control_count_column", "case_count_column",
        "control_count", "case_count", "imputation_info_column",
        "resolved_info_source", "fixed_info",
        "external_info_file", "external_info_column", "external_eaf_file",
        "external_eaf_column", "delimiter",
    ):
        assert field in text
    assert "dataset_id                     : ADHD2022" in text
    assert "effect_column                  : OR" in text
    assert "effect_type                    : odds_ratio" in text
    assert "z_score_column                 : Not provided — will not be used" in text
    assert "external_info_file             : Not provided — will not be used" in text
    assert "          🔹  dataset_id" in text
    assert "│" not in text
