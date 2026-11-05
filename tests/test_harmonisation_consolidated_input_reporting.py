"""Consolidated input and inference presentation retains measured evidence."""

from copy import deepcopy
from pathlib import Path

import polars as pl
import pytest

from postgwas.modules.harmonisation import service
from postgwas.modules.harmonisation.study_properties import resolve_study_properties


def test_consolidated_input_preserves_distinct_denominators_and_one_removal_ledger():
    report = {
        "input_rows": 11749765,
        "rows_removed_at_read_mandatory_gate": 315449,
        "fields": [
            {"field": "chr", "input_column": "CHR", "input_rows": 11749765,
             "post_parse_missing": 315449},
            {"field": "beta", "input_column": "BETA", "input_rows": 11749765,
             "post_parse_missing": 0},
            {"field": "eaf", "input_column": None, "source_kind": "configured_alternative",
             "source_values": {"eaffile": "/references/supplied.tsv.gz", "eafcolumn": "EUR"}},
        ],
        "numeric_conversion_by_column": {
            "BETA": {"status": "assessed", "invalid": 0, "non_missing": 11434316},
        },
    }
    original = deepcopy(report)
    columns = {"gwas_outputname": "fetal", "sumstat_file": "/input/fetal.tsv.gz",
               "trait_type": "quantitative", "delimiter": "tab"}
    ledger = dict(input_variants=11749765, variants_read=11749765,
        missing_required=315449, invalid_coordinates=0, unsupported_chromosomes=0,
        non_standard_alleles=0, duplicate_variants=0, ready_variants=11434316,
        ready_snps=11434316, ready_indels_or_other=0)
    rendered = service._field_completeness_block(
        report, validation_summary=ledger, sample_columns=columns)
    text = " ".join(rendered.split())
    assert rendered.count("Input validation") == 1
    assert "Input completeness and numeric conversion" not in rendered
    assert "Input validation summary" not in rendered
    assert "Validated sample-sheet values" not in rendered
    assert "Initial removals" not in rendered
    assert "315,449 / 11,749,765 (2.68%)" in text
    assert "0 / 11,434,316 (0.00%) of non-missing values assessed" in text
    assert "Missing read-stage mandatory values" in text and "Ready for harmonisation" in text
    assert rendered.count("/references/supplied.tsv.gz") == 1
    assert "/input/fetal.tsv.gz" in text and "fetal" in text
    assert report == original


def test_consolidated_input_keeps_unequal_read_counts_and_user_exclusion():
    text = " ".join(service._field_completeness_block(
        {"input_rows": 8, "fields": []},
        validation_summary=dict(input_variants=12, variants_read=10,
            excluded_chromosomes=2, missing_required=1, invalid_coordinates=1,
            unsupported_chromosomes=0, non_standard_alleles=0, duplicate_variants=0,
            ready_variants=6, ready_snps=5, ready_indels_or_other=1),
        sample_columns={"gwas_outputname": "counts"}).split())
    assert "Input variants : 12" in text
    assert "Variants read : 10" in text
    assert "User-excluded chromosomes : 2" in text
    assert "Variants assessed : 8" in text
    assert "Ready for harmonisation : 6" in text
    assert "Ready indels / other variants : 1" in text


def test_input_sources_reuse_distinct_file_references_without_reading(tmp_path, monkeypatch):
    from pathlib import Path
    from postgwas.config import load_configuration
    from postgwas.core.input_validation import InputValidationSession, record_file_validation
    from postgwas.core.validation_reporting import FileValidationDisplay

    paths = [tmp_path / "one" / "same.tsv", tmp_path / "two" / "same.tsv"]
    report = {"input_rows": 1, "fields": [
        {"field": "eaf", "input_column": None, "source_values": {"eaffile": str(paths[0]), "eafcolumn": "EUR"}},
        {"field": "info", "input_column": None, "source_values": {"infofile": str(paths[1]), "infocolumn": "INFO"}},
    ]}
    ledger = dict(input_variants=1, variants_read=1, missing_required=0,
        invalid_coordinates=0, unsupported_chromosomes=0, non_standard_alleles=0,
        duplicate_variants=0, ready_variants=1, ready_snps=1, ready_indels_or_other=0)
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        for path in paths:
            record_file_validation(path, "reference", checks=("schema",), status="passed")
        assert display.file_reference(str(paths[0])) == "File 1"
        assert display.file_reference(str(paths[1])) == "File 2"
        monkeypatch.setattr(Path, "open", lambda *a, **k: pytest.fail("input renderer reopened a source"))
        text = service._field_completeness_block(report, validation_summary=ledger,
            sample_columns={"gwas_outputname": "references"})
    assert "File 1" in text and "File 2" in text
    assert str(paths[0]) not in text and str(paths[1]) not in text


@pytest.mark.parametrize("emit", [True, False])
def test_orchestration_can_defer_inference_cards_without_changing_decisions(capsys, emit):
    frame = pl.DataFrame({"BETA": [-.2, .2, -.5, .5], "P": [.05, .2, .4, .8]})
    mapping = {"beta_or_col": "BETA", "pval_col": "P"}
    expected = resolve_study_properties(frame, mapping)
    capsys.readouterr()
    actual = resolve_study_properties(frame, mapping, report_inference_warnings=emit)
    output = capsys.readouterr().out
    assert actual == expected
    assert ("Effect type automatically inferred" in output) is emit
    assert ("P-value type automatically inferred" in output) is emit
    assert actual["effect_type_inference_warning"]
    assert actual["pvalue_type_inference_warning"]
    assert frame["BETA"].to_list() == [-.2, .2, -.5, .5]


def test_deferred_inference_preserves_canonical_warning_severity(tmp_path, capsys):
    from postgwas.core.pipeline_logging import PipelineLogger
    from postgwas.modules.harmonisation.policies import default_policies

    frame = pl.DataFrame({"BETA": [-.2, .2, -.5, .5], "P": [.05, .2, .4, .8]})
    policies = default_policies()
    with PipelineLogger("inferred", "dataset", str(tmp_path), policies=policies, screen_level="ERROR") as logger:
        decisions = resolve_study_properties(frame,
            {"beta_or_col": "BETA", "pval_col": "P"}, logger=logger, policies=policies,
            report_inference_warnings=False)
        assert logger.warning_count == 2
    visible = capsys.readouterr().out
    assert "Effect type automatically inferred" not in visible
    assert "P-value type automatically inferred" not in visible
    log = Path(logger.log_path).read_text()
    assert sum("WARNING" in line for line in log.splitlines()) == 2
    assert "effect-type inference" in log and "P-value-type inference" in log
    assert decisions["effect_type_inference_warning"] and decisions["pvalue_type_inference_warning"]
