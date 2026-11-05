"""AF reports keep source approval separate from original input semantics."""

from copy import deepcopy
from html import unescape

import pytest

from postgwas.core.pipeline_logging import _chromosome_analysis_lines
from postgwas.modules.harmonisation.html_report import render_dataset_report
from postgwas.modules.harmonisation.study_reporting import render_study_decisions


def selected():
    return {
        "eaf_column": "effect_allele_frequency", "eaf_selected_source": "replacement_reference",
        "eaf_input_interpretation": "maf", "eaf_is_maf_initial": True,
        "eaf_effective_frequency_screen": {"rows": 2000, "usable": 2000,
            "unusable": 0, "at_or_below_0_5": 2000, "maf_like": True},
        "eaf_reference_validation": {"status": "replacement_selected", "study_variant_count": 2000,
            "comparable_variants": 1800, "rule_1_passed": False, "rule_2_passed": True,
            "direct_eaf_correlation": -.590259, "mean_eaf_error": .557096,
            "supplied_to_reference_maf_correlation": .99363, "mean_maf_error": .0121395},
        "eaf_replacement": {"file": "/reference/replacement.tsv", "column": "EUR"},
        "eaf_replacement_validation": {"status": "passed", "contract": {"row_bound": True},
            "study_variant_count": 2000, "comparable_variants": 1700,
            "correlation": .998799, "mean_absolute_difference": .000181164,
            "effective_frequency_screen": {"rows": 2000, "usable": 1750, "unusable": 250,
                "at_or_below_0_5": 350, "maf_like": False}},
    }


def reports(decisions):
    return [" ".join(render_study_decisions(decisions).split()),
        " ".join(unescape(render_dataset_report({}, {"dataset": {"study_decisions": decisions}})).split())]


def test_approved_replacement_has_readable_shared_summary_and_preserves_signed_evidence():
    decisions = selected()
    before = deepcopy(decisions)
    for report in reports(decisions):
        assert "Replacement approved" in report
        assert "-0.590259" in report and "0.99363" in report
        assert "0.998799" in report
        assert "1,700 / 2,000 (85.00%)" in report
        assert "250" in report and "final removal count" in report
        assert "Original column" in report and "Preserved unchanged" in report
        assert "Initial selected AF A" not in report and "Replacement B" not in report
        assert "provisionally" not in report.lower()
    assert decisions == before


@pytest.mark.parametrize("status,contract,want", [
    ("passed", None, "approval not recorded"),
    ("invalid", None, "Invalid"),
    ("inconclusive", None, "Inconclusive"),
    ("pending", None, "pending"),
])
def test_incomplete_or_failed_replacement_never_claims_approval(status, contract, want):
    decisions = selected()
    decisions["eaf_replacement_validation"].update(status=status, contract=contract,
        failure_reason="replacement_inputs_changed", failure_message="source changed [EUR]")
    for report in reports(decisions):
        assert "Replacement approved" not in report
        assert want in report
        assert "source changed [EUR]" in report


def test_reference_inventory_is_one_summary_instead_of_one_line_per_chromosome():
    decisions = selected()
    resources = {str(i): {"default_eaf_file": f"/panel/chr{i}.tsv", "default_eaf_reference_column": "EUR"}
                 for i in range(1, 23)}
    report = render_study_decisions(decisions, resource_maps=resources)
    assert "22 chromosome files" in report
    assert "Comparison reference" not in report
    report = render_dataset_report({}, {"dataset": {"study_decisions": decisions,
        "resource_preflight": {"resources_by_chromosome": resources}}})
    assert "22 chromosome files" in report


def test_complete_input_does_not_repeat_zero_missing_and_identical_usable_count():
    report = reports(selected())[0]
    assert "Input usable frequencies" not in report
    assert "Input missing or invalid frequencies" not in report


def test_chromosome_uses_approved_replacement_without_original_maf_warning():
    decision = selected()
    qc = {"decision_source": "replacement_reference", "eaf_selected_source": "replacement_reference",
        "study_decision_eaf_is_maf": True, "eaf_provenance": "reference_imputed",
        "eaf_replacement_validation": decision["eaf_replacement_validation"],
        "eaf_replacement_source_file": "/reference/replacement.tsv",
        "eaf_replacement_source_column": "EUR", "final_eaf_col": "EUR_"}
    report = "\n".join(_chromosome_analysis_lines({"stage_qc": {"eaf_qc": qc}}, 120))
    assert "Approved replacement EAF" in report
    assert "whole-study decision reused" in report
    assert "MAF suspected / unresolved" not in report
    assert "Study column EUR_" not in report
    assert "inverted using" not in report


@pytest.mark.parametrize("status", ["failed", "invalid", "inconclusive"])
def test_native_failure_cannot_look_like_pending_or_retained_input(status):
    decisions = {"eaf_column": "FREQ", "eaf_reference_validation": {
        "status": status, "failure_reason": "inputs_changed"}}
    for report in reports(decisions):
        assert "stopped before chromosome analysis" in report
        assert "inputs changed" in report
        assert "assessment pending" not in report


def test_unapproved_worker_assignment_remains_unauthorised():
    qc = {"eaf_selected_source": None, "decision_source": "replacement_reference",
        "eaf_replacement_validation": {"status": "passed"}}
    report = "\n".join(_chromosome_analysis_lines({"stage_qc": {"eaf_qc": qc}}, 120))
    assert "approval not recorded" in report
    assert "Assignment not authorised" in report
    assert "Approved replacement EAF" not in report


@pytest.mark.parametrize("evaluated,want", [(None, "Outcome not recorded"), (False, "Not evaluated")])
def test_missing_or_unmeasured_comparisons_do_not_imply_failed_agreement(evaluated, want):
    decisions = selected()
    evidence = decisions["eaf_reference_validation"]
    if evaluated is None:
        evidence.pop("rule_1_passed")
        evidence.pop("rule_2_passed")
    else:
        evidence.update(whole_study_rules_evaluated=False, rule_1_passed=False, rule_2_passed=False,
            direct_eaf_correlation=None, mean_eaf_error=None, supplied_to_reference_maf_correlation=None, mean_maf_error=None)
    for report in reports(decisions):
        assert want in report
        assert "Failed; correlation -0.590259" not in report
