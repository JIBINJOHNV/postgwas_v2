"""Structured read-stage and final field-completeness reporting."""

from pathlib import Path

import polars as pl

from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.summary_statistics_io import (
    _input_field_completeness,
    finalise_field_completeness_report,
    read_summary_statistics,
    write_field_completeness_report,
)


OUTPUT_LAYOUT = {
    "rejected_directory": "rejected",
    "input_reject": "rejected/{dataset_id}_input.tsv",
    "duplicates": "{dataset_id}_duplicates.tsv",
}


def _field(report, name):
    return next(field for field in report["fields"] if field["field"] == name)


def _study(tmp_path: Path) -> Path:
    source = tmp_path / "study.tsv"
    source.write_text(
        "SNP\tCHR\tPOS\tEA\tOA\tBETA\tZ\tP\tSE\tEAF\n"
        "beta_row\t1\t100\tA\tG\t0.1\t.\t0.05\t.\t0.2\n"
        "z_row\t1\t101\tC\tT\t.\t2.0\t0.04\t.\t.\n"
        "missing_effect\t1\t102\tG\tA\t.\t.\t0.03\t.\t.\n",
        encoding="utf-8",
    )
    return source


def _mapping():
    return {
        "gwas_outputname": "study",
        "chr_col": "CHR",
        "pos_col": "POS",
        "snp_id_col": "SNP",
        "ea_col": "EA",
        "oa_col": "OA",
        "beta_or_col": "BETA",
        "imp_z_col": "Z",
        "pval_col": "P",
        "se_col": "SE",
        "eaf_col": None,
        "eaffile": "/references/study_eaf.tsv.gz",
        "eafcolumn": "EAF",
        "imp_info_col": None,
        "infofile": "/references/study_info.tsv.gz",
        "infocolumn": "INFO",
        "ncontrol_col": None,
        "ncontrol": 1000,
    }


def test_report_distinguishes_individual_missingness_and_beta_or_z_gate(tmp_path):
    result = read_summary_statistics(
        sumstat_file=str(_study(tmp_path)),
        output_dir=str(tmp_path),
        sample_column_dict=_mapping(),
        output_layout=OUTPUT_LAYOUT,
        table_delimiter="\t",
        policies=default_policies(),
        input_line_count=3,
    )

    retained = result[0]
    removed_missing = result[-2]
    report = result[-1]

    assert retained["SNP"].to_list() == ["beta_row", "z_row"]
    assert removed_missing == 1
    assert report["input_measurement_stage"] == (
        "post_parse_before_numeric_normalisation"
    )
    assert report["rows_removed_at_read_mandatory_gate"] == 1
    assert report["post_parse_missing_by_column"]["BETA"] == 2

    beta = _field(report, "beta")
    zscore = _field(report, "zscore")
    standard_error = _field(report, "se")
    eaf = _field(report, "eaf")
    info = _field(report, "info")

    assert beta["post_parse_missing"] == 2
    assert zscore["post_parse_missing"] == 2
    assert beta["read_requirement"] == "alternative"
    assert zscore["read_requirement"] == "alternative"
    assert standard_error["post_parse_missing"] == 3
    assert standard_error["read_requirement"] == "downstream_validated"
    assert standard_error["recovery_step"] == "chromosome step 08 or 06"
    assert "not row-wise filled" in standard_error["lifecycle_note"]
    assert eaf["source_kind"] == "configured_alternative"
    assert eaf["source_config_keys"] == ["eaffile", "eafcolumn"]
    assert info["source_config_keys"] == ["infofile", "infocolumn"]

    alternative = report["read_alternative_groups"][0]
    assert alternative["fields"] == ["beta", "zscore"]
    assert alternative["post_parse_missing_all"] == 1
    assert alternative["post_parse_missing_all_fraction"] == 1 / 3


def test_selected_internal_eaf_is_not_reported_as_external_row_backfill(tmp_path):
    mapping = _mapping()
    mapping["eaf_col"] = "EAF"
    report = read_summary_statistics(
        sumstat_file=str(_study(tmp_path)),
        output_dir=str(tmp_path),
        sample_column_dict=mapping,
        output_layout=OUTPUT_LAYOUT,
        table_delimiter="\t",
        policies=default_policies(),
        input_line_count=3,
    )[-1]

    eaf = _field(report, "eaf")
    assert eaf["source_kind"] == "study_column"
    assert eaf["source_config_keys"] == ["eaf_col"]
    assert eaf["post_parse_missing"] == 2
    assert eaf["read_requirement"] == "downstream_validated"
    assert "does not row-wise backfill" in eaf["lifecycle_note"]


def test_combined_coordinates_are_reported_as_configured_alternative():
    mapping = _mapping()
    mapping.update({
        "chr_col": None,
        "pos_col": None,
        "chr_pos_col": "CHRPOS",
    })
    report = _input_field_completeness(
        pl.DataFrame({
            "CHRPOS": ["1:100"],
            "EA": ["A"],
            "OA": ["G"],
            "BETA": [0.1],
            "Z": [None],
            "P": [0.05],
            "SE": [None],
        }),
        mapping,
        default_policies(),
    )

    chromosome = _field(report, "chr")
    position = _field(report, "pos")
    assert chromosome["source_config_keys"] == ["chr_pos_col"]
    assert position["source_config_keys"] == ["chr_pos_col"]
    assert chromosome["input_column"] == "CHRPOS"
    assert position["post_parse_missing"] == 0
    assert chromosome["read_requirement"] == (
        "required_via_configured_alternative"
    )
    assert position["read_action"] == (
        "validate_configured_alternative_during_normalisation"
    )


def test_absent_z_does_not_borrow_effect_missingness():
    mapping = _mapping()
    mapping["imp_z_col"] = None
    report = _input_field_completeness(pl.DataFrame({"BETA": [0.1, None]}), mapping, default_policies())
    assert _field(report, "zscore")["post_parse_missing"] is None
    assert _field(report, "zscore")["input_column"] is None
    assert _field(report, "beta")["post_parse_missing"] == 1


def test_conversion_evidence_does_not_change_values_and_marks_deferred_info():
    from postgwas.modules.harmonisation.summary_statistics_io import normalise_summary_statistics_values

    frame = pl.DataFrame({"P": ["6.91E-08", "bad", None], "INFO": ["0.8,0.9", "bad", None]})
    kwargs = dict(chr_col=None, numeric_columns=["P", "INFO"], compound_numeric_columns=["INFO"])
    expected = normalise_summary_statistics_values(frame, **kwargs)
    counts = {}
    observed = normalise_summary_statistics_values(frame, conversion_counts=counts, **kwargs)
    assert observed.equals(expected)
    assert observed["P"].to_list() == [6.91e-8, None, None]
    assert counts["P"] == {"status": "assessed", "rows": 3, "non_missing": 2, "invalid": 1}
    assert counts["INFO"]["status"] == "compound_values_deferred"


def test_screen_has_separate_fields_precise_small_percentages_and_resolved_actions():
    from postgwas.modules.harmonisation.service import _field_completeness_block

    mapping = _mapping()
    mapping.update(imp_z_col=None, imp_info_col="INFO")
    report = _input_field_completeness(
        pl.DataFrame({"BETA": ["bad"], "SNP": [None], "INFO": [None]}),
        mapping, default_policies(),
    )
    _field(report, "snp").update(input_rows=6774224, post_parse_missing=7)
    report["numeric_conversion_by_column"] = {"BETA": {"status": "assessed", "rows": 1, "non_missing": 1, "invalid": 1}}
    rendered = " ".join(_field_completeness_block(report).split())
    for text in ["Input completeness and numeric conversion", "Missing in input", "Invalid numeric",
                 "7 / 6,774,224 (<0.01%)", "set to missing", "Calculate from harmonised BETA and SE",
                 "chromosome_position_effect_other", "Keep variants with missing INFO"]:
        assert text in rendered
    for text in ["read role", "final gate", "configured beta_or_col"]:
        assert text not in rendered


def test_reader_records_conversion_denominator_after_initial_filtering(tmp_path):
    source = tmp_path / "numeric.tsv"
    source.write_text("CHR\tPOS\tEA\tOA\tBETA\tP\n1\t100\tA\tG\t0.1\tbad\n1\t101\tA\tG\t0.2\t6.91E-08\n1\t102\tA\tG\t0.3\tNA\n")
    mapping = _mapping()
    mapping.update(snp_id_col=None, imp_z_col=None, se_col=None)
    report = read_summary_statistics(
        sumstat_file=str(source), output_dir=str(tmp_path), sample_column_dict=mapping,
        output_layout=OUTPUT_LAYOUT, table_delimiter="\t", policies=default_policies(), input_line_count=3,
    )[-1]
    assert _field(report, "pval")["post_parse_missing"] == 1
    assert report["numeric_conversion_by_column"]["P"] == {"status": "assessed", "rows": 2, "non_missing": 2, "invalid": 1}
    path = tmp_path / "completeness.tsv"
    write_field_completeness_report(path, report, "\t")
    table = pl.read_csv(path, separator="\t")
    assert '"invalid": 1' in table.filter(pl.col("field") == "pval").item(0, "numeric_conversion")


def test_final_report_aggregates_post_recovery_gate_without_double_counting(tmp_path):
    report = read_summary_statistics(
        sumstat_file=str(_study(tmp_path)),
        output_dir=str(tmp_path),
        sample_column_dict=_mapping(),
        output_layout=OUTPUT_LAYOUT,
        table_delimiter="\t",
        policies=default_policies(),
        input_line_count=3,
    )[-1]
    required = report["final_check"]["required_fields"]
    summaries = {
        "1": {
            "stage_qc": {
                "final_completeness_qc": {
                    "initial_variants": 100,
                    "final_variants": 97,
                    "missing_by_field": {
                        field: (2 if field == "eaf" else 1 if field == "se" else 0)
                        for field in required
                    },
                }
            }
        },
        "2": {
            "stage_qc": {
                "final_completeness_qc": {
                    "initial_variants": 50,
                    "final_variants": 49,
                    "missing_by_field": {
                        field: (1 if field == "eaf" else 0)
                        for field in required
                    },
                }
            }
        },
    }

    final = finalise_field_completeness_report(report, summaries, ["1", "2"])

    assert final["final_check"]["status"] == "complete"
    assert final["final_check"]["rows_entering_gate"] == 150
    assert final["final_check"]["rows_retained_after_gate"] == 146
    assert final["final_check"]["measurement"] == (
        "sequential rule attribution in final_check.require order"
    )
    assert _field(final, "eaf")["final_missing"] == 3
    assert _field(final, "se")["final_missing"] == 1
    assert _field(final, "info")["final_status"] == "not_required"

    path = tmp_path / "study_field_completeness.tsv"
    write_field_completeness_report(path, final, delimiter="\t")
    table = pl.read_csv(path, separator="\t")
    eaf_row = table.filter(
        (pl.col("record_type") == "field") & (pl.col("field") == "eaf")
    )
    alternative_row = table.filter(
        pl.col("record_type") == "read_alternative_group"
    )

    assert eaf_row.item(0, "final_missing") == 3
    assert eaf_row.item(0, "final_order") == 3
    assert eaf_row.item(0, "source_config_keys") == '["eaffile", "eafcolumn"]'
    assert eaf_row.item(0, "read_gate_rows_rejected_total") == 1
    assert eaf_row.item(0, "final_expected_chromosomes") == '["1", "2"]'
    assert eaf_row.item(0, "final_assessed_chromosomes") == '["1", "2"]'
    assert eaf_row.item(0, "final_report_status") == "complete"
    assert "does not row-wise backfill" in eaf_row.item(0, "lifecycle_note")
    assert alternative_row.item(0, "field") == '["beta", "zscore"]'
    assert alternative_row.item(0, "post_parse_missing") == 1


def test_final_report_marks_partial_chromosome_and_keep_count_scope(tmp_path):
    report = read_summary_statistics(
        sumstat_file=str(_study(tmp_path)),
        output_dir=str(tmp_path),
        sample_column_dict=_mapping(),
        output_layout=OUTPUT_LAYOUT,
        table_delimiter="\t",
        policies=default_policies(),
        input_line_count=3,
    )[-1]
    report["final_check"]["on_missing"] = "keep"
    required = report["final_check"]["required_fields"]
    summaries = {
        "1": {
            "stage_qc": {
                "final_completeness_qc": {
                    "initial_variants": 10,
                    "final_variants": 10,
                    "missing_by_field": {
                        field: (1 if field in {"eaf", "se"} else 0)
                        for field in required
                    },
                }
            }
        },
        "2": {"stage_qc": {"error": "chromosome failed before final gate"}},
    }

    final = finalise_field_completeness_report(
        report,
        summaries,
        ["1"],
        expected_chromosomes=["1", "2"],
    )

    assert final["final_check"]["status"] == "partial"
    assert final["final_check"]["unassessed_chromosomes"] == ["2"]
    assert final["final_check"]["measurement"] == (
        "independent per-field matches; rows may appear in multiple counts"
    )
    assert _field(final, "eaf")["final_missing"] == 1
    assert _field(final, "eaf")["final_status"] == (
        "partial_chromosome_coverage"
    )
