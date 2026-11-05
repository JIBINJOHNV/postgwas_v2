"""Saved GWAS-row screens remain distinct from source-panel descriptions."""

from copy import deepcopy
import html
import json

import pytest

from postgwas.modules.harmonisation.html_report import render_dataset_report
from postgwas.modules.harmonisation.study_reporting import render_study_decisions


def _screen(*, rows=210, usable=200, low=190, maf_like=False):
    return {
        "scope": "merged_gwas_rows", "allele_direction": "original_gwas_effect_allele",
        "column_present": True, "rows": rows, "usable": usable, "unusable": rows - usable,
        "at_or_below_0_5": low, "low_fraction_of_usable": low / usable if usable else 0.0,
        "decision_cutoff": .95, "maf_like": maf_like,
    }


def _decisions():
    return {
        "eaf_column": None, "eaf_is_maf_source": "no_study_frequency_column",
        "eaf_external_validation_required": True, "eaf_selected_source": "external_eaf",
        "eaf_external_validation": {
            "status": "passed", "contract": {"row_bound": True}, "study_variant_count": 210,
            "decision_reason": "external_frequency_not_predominantly_minor",
            "raw_frequency_screen": {**_screen(rows=1000, usable=1000, low=999, maf_like=True),
                "scope": "source_panel_rows", "allele_direction": "source_counted_allele"},
            "effective_frequency_screen": _screen(),
            "source": {"file": "/panel/initial.tsv", "column": "EUR"},
        },
    }


def _reports(decisions):
    return (
        " ".join(render_study_decisions(decisions).split()),
        " ".join(html.unescape(render_dataset_report({}, {"dataset": {"study_decisions": decisions}})).split()),
    )


def test_external_reports_screen_the_saved_gwas_overlap_and_skip_default_diagnostic():
    decisions = _decisions()
    before = deepcopy(decisions)
    for report in _reports(decisions):
        assert "Screen scope" in report
        assert "merged GWAS rows" in report
        assert "original listed GWAS effect allele" in report
        assert "190 / 200 (95.00%)" in report
        assert "strictly above 95.00%" in report
        assert "999 / 1,000" not in report
        assert "Reference comparison" in report
        assert "Not required; input was not MAF-like" in report
        assert "Supplied external frequency; no native GWAS frequency column mapped" in report
    assert decisions == before


@pytest.mark.parametrize("usable,low,reason,result", [
    (200, 191, "replacement_frequency_maf_like", "MAF-like; replacement rejected"),
    (0, 0, "no_usable_replacement_frequencies", "No usable frequencies; replacement rejected"),
])
def test_replacement_screen_failure_is_terminal_without_fabricated_consistency_metrics(usable, low, reason, result):
    decisions = _decisions()
    decisions.update(eaf_selected_source="replacement_reference")
    decisions["eaf_external_validation"].update(status="replacement_selected", contract=None)
    decisions["eaf_replacement_validation"] = {
        "status": "invalid", "decision_reason": reason, "study_variant_count": 210,
        "source": {"file": "/panel/replacement.tsv", "column": "EUR"},
        "effective_frequency_screen": _screen(usable=usable, low=low, maf_like=bool(usable)),
    }
    for report in _reports(decisions):
        assert "Replacement distribution" in report
        assert result in report
        assert reason.replace("_", " ") in report
        assert "No further default comparison, replacement or fallback" in report
        assert "accepted for assignment" not in report
        assert "Replacement correlation" not in report
        assert "/panel/replacement.tsv" in report


@pytest.mark.parametrize("maf_source", ["no_study_frequency_column", "dataset_reference_maf_support"])
def test_external_to_replacement_reports_preserve_absent_native_and_final_adoption(maf_source):
    decisions = _decisions()
    decisions.update(eaf_selected_source="replacement_reference", eaf_is_maf_source=maf_source)
    decisions["eaf_external_validation"].update(status="replacement_selected", contract=None)
    decisions["eaf_replacement_validation"] = {
        "status": "passed", "contract": {"row_bound": True},
        "source": {"file": "/panel/replacement.tsv", "column": "EUR"},
        "effective_frequency_screen": _screen(usable=200, low=100),
    }
    for report in _reports(decisions):
        assert "Replacement approved" in report
        assert "external EAF accepted for assignment" not in report
        assert "Original study MAF remains unchanged" not in report
        assert "Original column" in report and "Preserved unchanged" in report
        assert "no native GWAS frequency column mapped" in report


def test_native_screen_records_original_effect_allele_scope():
    decisions = {
        "eaf_column": "FREQ", "eaf_is_maf_initial": False,
        "eaf_evidence": _screen(), "eaf_reference_validation": {"status": "not_required"},
    }
    for report in _reports(decisions):
        assert "Screen scope" in report
        assert "original listed GWAS effect allele" in report
        assert "190 / 200 (95.00%)" in report


def test_terminal_native_screen_has_one_owner_for_counts_fraction_and_cutoff():
    decisions = {"eaf_column": "FREQ", "eaf_is_maf_initial": False,
        "eaf_evidence": {**_screen(), "input_missing": 7, "unparseable": 1,
            "non_finite": 1, "out_of_range": 1, "min": .1, "max": .8, "median": .3},
        "eaf_maf_decision_cutoff": .95, "eaf_reference_validation": {"status": "not_required"}}
    report = _reports(decisions)[0]
    assert report.count("190 / 200 (95.00%)") == 1
    assert report.count("strictly above 95.00%") == 1
    assert report.count("Input usable frequencies : 200") == 1
    assert report.count("Input missing or invalid frequencies : 10") == 1
    assert "Usable frequencies : 200" not in report
    assert "Missing or invalid : 10" not in report
    assert "Missing input: 7" in report
    assert "Detailed evidence" in report


def test_native_screen_retains_its_saved_cutoff_without_a_second_legacy_field():
    decisions = {"eaf_column": "FREQ", "eaf_is_maf_initial": False,
        "eaf_evidence": {**_screen(), "decision_cutoff": .912},
        "eaf_reference_validation": {"status": "not_required"}}
    report = _reports(decisions)[0]
    assert report.count("strictly above 91.20%") == 1
    assert "strictly above not available" not in report


def test_native_initial_and_replacement_screens_keep_independent_labels_and_denominators():
    decisions = {"eaf_column": "FREQ", "eaf_is_maf_initial": True,
        "eaf_evidence": _screen(low=199, maf_like=True), "eaf_selected_source": "replacement_reference",
        "eaf_reference_validation": {"status": "replacement_selected"},
        "eaf_replacement_validation": {"status": "passed", "contract": {"row_bound": True},
            "source": {"file": "/panel/replacement.tsv", "column": "EUR"},
            "effective_frequency_screen": _screen(usable=180, low=120)}}
    report = _reports(decisions)[0]
    assert report.count("199 / 200 (99.50%)") == 1
    assert report.count("120 / 180 (66.67%)") == 1
    assert "Input usable frequencies : 200" in report
    assert "Unresolved at preparation : 30" in report


@pytest.mark.parametrize("replacement", [False, True])
def test_external_source_report_summarises_files_without_dumping_identity_maps(replacement):
    decisions = _decisions()
    source = {"column": "EUR", "mapping": {"chr": "CHROM", "pos": "POS", "a1": "REF", "a2": "ALT"},
        "declared_files": {str(c): "/panel/declared_chr%s.tsv" % c for c in range(1, 23)},
        "resolved_files": {str(c): "/cache/resolved_chr%s.tsv" % c for c in range(1, 23)},
        "identities": {str(c): {"path": "/cache/resolved_chr%s.tsv" % c, "size": c} for c in range(1, 23)}}
    key = "eaf_replacement_validation" if replacement else "eaf_external_validation"
    if replacement:
        decisions.update(eaf_selected_source="replacement_reference")
        decisions["eaf_external_validation"].update(status="replacement_selected", contract=None)
        decisions[key] = {"status": "passed", "contract": {"row_bound": True}}
    decisions[key]["source"] = source
    before = deepcopy(decisions)
    report = _reports(decisions)[0]
    assert "column EUR" in report
    assert "22 declared files" in report
    assert "source identities" in report
    assert "allele mappings" in report
    assert "/cache/resolved_chr" not in report
    assert "/panel/declared_chr" not in report
    assert decisions == before


def test_html_describes_saved_semantic_approval_and_later_output_validation():
    from postgwas.modules.harmonisation.html_report import _decision_section
    report = _decision_section({}, {"dataset": {"study_decisions": {"eaf_column": "AF"}}})
    assert "whole-study frequency source approval" in report
    assert "chromosome and output QC" in report
    assert "subsequent chromosome-level reference validation" not in report


@pytest.mark.parametrize("external", [False, True])
def test_zero_usable_initial_screen_is_unavailable_rather_than_acceptable_non_maf(external):
    if external:
        decisions = _decisions()
        decisions["eaf_external_validation"].update(status="inconclusive", contract=None,
            decision_reason="no_usable_external_frequencies", effective_frequency_screen=_screen(usable=0, low=0))
    else:
        decisions = {"eaf_column": "AF", "eaf_evidence": _screen(usable=0, low=0),
            "eaf_reference_validation": {"status": "failed", "decision_reason": "no_usable_study_frequencies"}}
    report = _reports(decisions)[0]
    assert "Input distribution : No usable frequencies" in report
    assert "accepted for assignment" not in report
    assert "EAF retained" not in report


def _replacement_vcf_values(tmp_path, *, external=True, initial_source=None, include_initial=True):
    from postgwas.config import load_configuration
    from postgwas.modules.harmonisation.policies import default_policies
    from postgwas.modules.harmonisation.vcf_provenance import build_harmonisation_vcf_provenance
    initial = str(tmp_path / "initial_panel.tsv")
    replacement = str(tmp_path / "replacement_panel.tsv")
    native_column = None if external else "GWAS_MAF"
    decisions = {"effect_type": "beta", "pvalue_type": "raw", "eaf_column": native_column,
        "frequency_type": "effect_allele_frequency", "eaf_selected_source": "replacement_reference",
        "eaf_input_interpretation": "maf", "eaf_is_maf_source": "replacement_reference",
        "eaf_reference_validation": {"status": "replacement_selected"},
        "eaf_replacement": {"file": replacement, "column": "B_AF", "resolved_files": {"1": replacement}}}
    if external and include_initial:
        decisions["eaf_initial_frequency"] = {"source_kind": "external_eaf",
            "column": "__postgwas_initial_frequency", "source": initial_source if initial_source is not None else {
                "column": "SOURCE_AF", "declared_files": {"1": initial}, "resolved_files": {"1": initial}}}
    config = load_configuration().modules.harmonisation.vcf_processing.model_dump()
    values = build_harmonisation_vcf_provenance(
        sample={"gwas_outputname": "study", "sumstat_file": str(tmp_path / "study.tsv"),
            "output_root": str(tmp_path / "output"), "beta_or_col": "BETA", "se_col": "SE",
            "imp_z_col": "Z", "pval_col": "P", "eaf_col": native_column,
            "eaffile": replacement, "eafcolumn": "B_AF", "info_source": "fixed_cli", "fixed_info": .99},
        study_decisions=decisions, resource_preflight={"require_default_eaf": False,
            "resources_by_chromosome": {"1": {"user_eaf_file": replacement}}},
        completed_chromosomes=["1"], policies=default_policies(),
        resource_directory=str(tmp_path), output_directory=str(tmp_path / "output"),
        run_manifest=str(tmp_path / "manifest.json"), vcf_config=config)
    return values, config


def test_external_a_to_b_vcf_headers_keep_both_sources_and_explicit_native_absence(tmp_path):
    from postgwas.modules.harmonisation.vcf_processing import _render_provenance_headers
    values, config = _replacement_vcf_values(tmp_path)
    assert values["af_input_column"] == "SOURCE_AF"
    assert values["af_native_input_column"] == "not_provided"
    assert values["af_initial_source"] == "external user-provided reference proxy; not study-measured AF"
    assert values["af_initial_reference_files"] == str((tmp_path / "initial_panel.tsv").resolve())
    assert values["af_reference_column"] == "B_AF"
    assert values["af_reference_files"] == str((tmp_path / "replacement_panel.tsv").resolve())
    assert values["af_input_type"] == "minor_allele_frequency"
    assert "initial external AF" in values["af_type_decision"]
    assert "original study MAF" not in values["af_harmonisation"].lower()
    assert "initial external AF A preserved" in values["af_swapped_formula"]
    assert "__postgwas_initial_frequency" not in json.dumps(values)
    assert "initial_panel.tsv" in values["resolved_resource_files"]
    headers = _render_provenance_headers(values, vcf_config=config, input_build="GRCh37",
        output_build="GRCh38", merge_group="grch38")
    observed = {line.split("=", 1)[0][2:]: json.loads(line.split("=", 1)[1]) for line in headers}
    for key in ("af_input_column", "af_native_input_column", "af_initial_source", "af_initial_reference_files", "af_reference_column"):
        assert observed[config["provenance"]["headers"][key]] == values[key]


def test_native_a_replacement_vcf_provenance_retains_the_existing_study_contract(tmp_path):
    values, _ = _replacement_vcf_values(tmp_path, external=False)
    assert values["af_input_column"] == values["af_native_input_column"] == "GWAS_MAF"
    assert values["af_initial_source"] == "study column"
    assert values["af_initial_reference_files"] == "not_provided"
    assert "Original study MAF preserved" in values["af_harmonisation"]
    assert "original study MAF is preserved" in values["af_swapped_formula"]


@pytest.mark.parametrize("source", [
    {"column": None, "resolved_files": {"1": "/panel/initial.tsv"}},
    {"column": "SOURCE_AF", "resolved_files": {}},
])
def test_external_a_vcf_provenance_rejects_incomplete_original_source_instead_of_using_b(tmp_path, source):
    with pytest.raises(ValueError, match="initial external AF"):
        _replacement_vcf_values(tmp_path, initial_source=source)


def test_external_a_vcf_provenance_does_not_use_b_when_initial_descriptor_is_missing(tmp_path):
    with pytest.raises(ValueError, match="initial external AF"):
        _replacement_vcf_values(tmp_path, include_initial=False)
