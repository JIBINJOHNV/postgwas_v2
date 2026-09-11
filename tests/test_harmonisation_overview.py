"""Scientist-facing opening summary must preserve evidence and denominators."""

from copy import deepcopy
from html import unescape
import re

import pytest

from postgwas.modules.harmonisation.html_report import (
    render_dataset_report,
    render_run_report,
)


_HEADINGS = (
    "1.1 Study and run",
    "1.2 Input validation",
    "1.3 Study-wide decisions",
    "1.4 Strand and allele-frequency harmonisation",
    "1.5 Sample size and statistical handling",
    "1.6 Export, VCF creation, and liftover",
    "1.7 Reference-frequency QC",
    "1.8 Final VCF QC",
    "1.9 Input–VCF concordance",
    "1.10 Final takeaway",
)


def _opening(html: str, dataset_id: str = "overview_study") -> str:
    """Limit assertions to section 1, not matching evidence further below."""
    marker = f'id="dataset-{dataset_id}-section-1"'
    return html.split(marker, 1)[1].split('<section class="report-section"', 1)[0]


def _block(html: str, number: int) -> str:
    heading = _HEADINGS[number - 1]
    section = html.split(f"<h4>{heading}</h4>", 1)[1]
    if number < len(_HEADINGS):
        section = section.split(f"<h4>{_HEADINGS[number]}</h4>", 1)[0]
    return section


def _text(html: str) -> str:
    return " ".join(unescape(re.sub(r"<[^>]*>", " ", html)).split())


def _row(html: str, label: str) -> str:
    match = re.search(r'<th scope="row">' + re.escape(label) + r"</th><td>(.*?)</td>", html)
    assert match is not None, f"Missing overview row: {label}"
    return _text(match[1])


def _metric(html: str, label: str) -> str:
    match = re.search(
        r'<div class="metric-label">' + re.escape(label)
        + r'</div><div class="metric-value">(.*?)</div>', html,
    )
    assert match is not None, f"Missing overview metric: {label}"
    return _text(match[1])


def _population_row(html: str, label: str) -> list[str]:
    for row in re.findall(r"<tr>(.*?)</tr>", html):
        cells = [_text(cell) for cell in re.findall(r"<td>(.*?)</td>", row)]
        if len(cells) > 1 and cells[0] == label:
            return cells
    raise AssertionError(f"Missing population-comparison row: {label}")


def _evidence(dataset_id: str = "overview_study") -> tuple[dict, dict]:
    """Small persisted-evidence fixture; no scientific data or tool execution."""
    record = {
        "dataset_id": dataset_id, "status": "OK",
        "input_file": f"/study/{dataset_id}.tsv.gz",
        "input_variants": 100, "parsed_variants": 100,
        "input_ready_variants": 95, "input_ready_snps": 90,
        "input_ready_indels_or_other_variants": 5,
        "harmonised_variants": 94, "rejected_variants": 6,
        "row_accounting_balanced": True, "row_accounting_complete": True,
        "genome_build": "GRCh37", "genome_build_source": "detector",
        "effect_type": "odds_ratio", "effect_type_source": "detector",
        "effect_scale": "log-odds BETA", "p_value_type": "raw",
        "p_value_type_source": "sample_sheet", "frequency_type": "effect_allele_frequency",
        "frequency_type_source": "study_level_statistic", "strand_consensus": "forward",
        "strand_status": "complete", "strand_metadata_complete": True,
        "strand_chromosomes_completed": 1, "strand_chromosomes_expected": 1,
        "strand_reference_panel": "alignment_panel", "strand_reference_population": "EUR",
        "strand_variants_evaluated": 95, "strand_matched": 94,
        "strand_forward": 90, "strand_forward_swapped": 4,
        "strand_reverse_complement": 0, "strand_reverse_complement_swapped": 0,
        "strand_reference_unmatched": 1, "strand_removed_total": 1,
        "final_vcf_variants": 94, "final_vcf_snps": 89,
        "final_vcf_indels_or_other_variants": 5,
        "qc_passed_variants": 88, "qc_passed_snps": 84, "qc_failed_snps": 5,
        "qc_active_rule_count": 2, "qc_passed_snp_percent": 100 * 84 / 89,
        "comparison_af_panel": "comparison_panel", "comparison_af_population": "EUR",
        "closest_population": "EUR", "af_comparable_variants": 90,
        "af_concordant_variants": 89, "af_mismatched_variants": 1,
        "af_concordance_percent": 100 * 89 / 90,
        "af_mismatch_percent": 100 / 90, "af_difference_cutoff": 0.17,
        "study_af_missing": 0, "reference_af_missing": 4,
        "qc_passed_low_neff_variants": 2, "qc_passed_low_neff_percent": 100 * 2 / 88,
        "neff_minimum_threshold": 160, "qc_passed_missing_or_invalid_neff": 0,
        "qc_passed_missing_imputation_score": 0,
    }
    manifest = {
        "sample_id": dataset_id, "status": "OK", "qc_genome_build": "GRCh37",
        "sumstat_file": record["input_file"], "output_folder": f"/results/{dataset_id}",
        "elapsed_seconds": 60,
        "config": {"declared_effect_type": "auto", "declared_pvalue_type": "raw",
                   "beta_or_col": "OR_input", "se_col": "SE_input", "pval_col": "P_input",
                   "eaf_col": "frequency_input", "imp_info_col": "quality_input",
                   "info_source_detail": "internal column quality_input"},
        "policies": {"build.mode": "auto", "effect.type": "auto", "pvalue.type": "auto",
                     "info.on_missing": "keep", "strand.af_discordance_action": "warn"},
        "report_context": {"configuration": {"harmonisation": {
            "vcf_processing": {"target_builds": {"GRCh37": "GRCh38"}},
        }}},
        "pre_vcf_summary": {
            "total_variant_infile": 100, "total_variant_read": 100,
            "total_variant_removed_missing_values": 2,
            "total_variant_removed_null_coords": 0,
            "total_variant_removed_non_standard_alleles": 0,
            "total_variant_removed_duplicates": 3,
            "total_variant_remaining_for_harmonisation": 95,
            "total_variant_ready_snps": 90, "total_variant_ready_indels_or_other": 5,
            "total_variant_in_vcf_input": 94,
            "total_variant_removed_chromosome_harmonisation": 1,
            "total_variant_rejected_end_to_end": 6,
            "row_accounting_balanced": True, "row_accounting_complete": True,
        },
        "dataset": {
            "completed": ["1"], "failed": [], "chromosomes": {"1": {"status": "ok"}},
            "genome_build": {"inferred_build": "GRCh37", "mode": "auto", "forced": False,
                             "matches": {"GRCh37": 80, "GRCh38": 1},
                             "testable_variants": 81, "input_variants": 100},
            "study_decisions": {
                "effect_type": "odds_ratio", "effect_type_source": "detector",
                "pvalue_type": "raw", "pvalue_type_source": "sample_sheet",
                "se_scale": "log_odds", "se_scale_source": "dataset_z_pvalue_cross_check",
                "frequency_type": "effect_allele_frequency", "eaf_is_maf_initial": False,
                "eaf_provenance": "study_supplied", "eaf_provenance_counts": {"study_supplied": 1},
                "strand": "forward", "strand_consensus": {"dominant_fraction": 1.0},
                "info_score_type": "standard_info", "info_score_type_source": "automatic_full_dataset",
            },
            "chromosome_summaries": {"1": {
                "status": "ok", "rows_in": 95, "rows_out": 94, "rejected": 1, "balanced": True,
                "steps": [{"number": 16, "total": 16, "status": "ok", "rows_out": 94,
                           "extra": {"liftover_accounting": {"input": 94, "rejected": 1,
                                     "swap_excluded": 2, "swap_policy": "exclude", "final": 91},
                                     "target_vcf": f"/results/{dataset_id}_GRCh38.vcf.gz"}}],
                "stage_qc": {
                    "beta_or_oddsratio_qc": {"conversion": "OR_to_Beta_log_transform_applied",
                                             "se_rescaled_by_or": False},
                    "eaf_qc": {"variants_removed": 1, "final_variants": 94,
                               "eaf_provenance": "study_supplied",
                               "strand_non_palindromic_af_comparable": 90,
                               "strand_non_palindromic_af_discordant": 3,
                               "strand_af_discordance_action": "warn"},
                    "sample_size_qc": {"Neff_status": "calculated_from_case_control",
                                        "ncase_source": "column:case_input", "ncontrol_source": "column:control_input",
                                        "sample_size_stats": {"Neff": {"min": 100, "max": 240}}},
                    "pval_detection_and_conversion_qc": {"variants_with_zero_pvalues_replaced": 0},
                    "info_qc": {"source": "study_column:quality_input", "info_score_type": "standard_info",
                                "rescaled_within_tolerance": 2, "rejected_out_of_range": 0,
                                "rejected_missing": 0},
                },
            }},
        },
        "vcf_provenance": {
            "af_source": "study column", "af_input_column": "frequency_input",
            "af_harmonisation": "Frequency oriented to the final VCF ALT allele",
            "info_source": "study column", "info_input_column": "quality_input",
            "info_interpretation": "Study-provided imputation-quality score",
            "effect_harmonisation": "Converted supplied odds ratios to log-odds BETA",
            "effect_formula": "BETA = ln(OR)", "se_input_scale": "log-odds",
            "se_harmonisation": "Supplied SE retained without scale conversion",
            "sample_size_source": "case and control counts supplied by column or fixed value",
            "sample_size_formula": "Neff = 4 / (1/Ncase + 1/Ncontrol)",
        },
        "merged_vcfs": {"merge_status": "OK", "grch37": f"/results/{dataset_id}_GRCh37.vcf.gz",
                        "grch38": f"/results/{dataset_id}_GRCh38.vcf.gz"},
        "qc_assessment": {"raw_variants": 94, "qc_passed_variants": 88, "excluded": 6,
                          "retained_fraction": 88 / 94},
        "population_frequency_qc": {"status": "complete", "closest_population": "EUR",
                                    "comparable_variants": 85},
    }
    return record, manifest


def test_ten_overview_blocks_are_inside_section_one_in_requested_order():
    html = render_dataset_report(*_evidence())
    opening = _opening(html)
    assert re.findall(r"<h4>(1\.\d+ [^<]*)</h4>", opening) == list(_HEADINGS)
    assert "2. Inputs and initial checks" not in opening
    assert "Chromosome-wise processing" in html


def test_run_report_includes_ten_blocks_for_every_dataset_in_sample_order():
    first, first_manifest = _evidence("first")
    second, second_manifest = _evidence("second")
    second["status"] = "PREFLIGHT_FAILED"
    second["failure_reason"] = "Missing study reference"
    html = render_run_report([first, second], {"first": first_manifest, "second": second_manifest})
    assert html.index('id="dataset-first"') < html.index('id="dataset-second"')
    for dataset in ("first", "second"):
        opening = _opening(html, dataset)
        assert re.findall(r"<h4>(1\.\d+ [^<]*)</h4>", opening) == list(_HEADINGS)
    assert "Missing study reference" in _opening(html, "second")


def test_overview_preserves_input_objects_and_does_not_open_scientific_files(monkeypatch):
    record, manifest = _evidence()
    original = deepcopy((record, manifest))

    def unexpected_open(*args, **kwargs):
        raise AssertionError("Rendering must not open a scientific input or output")

    monkeypatch.setattr("builtins.open", unexpected_open)
    monkeypatch.setattr("pathlib.Path.open", unexpected_open)
    render_dataset_report(record, manifest)
    assert (record, manifest) == original


@pytest.mark.parametrize("requested,expected", [(None, "Not recorded"), (False, "Not requested"), (True, "Not recorded")])
def test_missing_concordance_result_is_not_an_invented_pass_or_request_state(requested, expected):
    record, manifest = _evidence()
    if requested is not None:
        manifest["report_context"] = {"concordance_requested": requested}
    block = _text(_block(_opening(render_dataset_report(record, manifest)), 9))
    assert expected in block
    assert "not a passing" in block.lower() or "not a pass" in block.lower()
    if requested is None:
        assert "Not requested" not in block


def test_overview_distinguishes_recorded_zero_from_missing_input_removal_count():
    record, manifest = _evidence()
    block = _block(_opening(render_dataset_report(record, manifest)), 2)
    assert _row(block, "Invalid coordinates — removed") == "0 / 100 (0.0000% of rows read)"
    assert _row(block, "Unsupported chromosomes — removed") == "Not recorded / 100"
    assert "95 / 100" in _row(block, "Ready for harmonisation")
    assert "90 / 95" in _row(block, "Ready SNPs")
    assert "5 / 95" in _row(block, "Ready indels / other variants")


def test_rare_variants_do_not_round_to_zero_or_full_coverage_percentages():
    record = {"dataset_id": "overview_study", "input_variants": 10_000_000,
              "parsed_variants": 10_000_000, "input_ready_variants": 10_000_000,
              "input_ready_snps": 9_999_999, "input_ready_indels_or_other_variants": 1}
    block = _block(_opening(render_dataset_report(record, {})), 2)
    assert _row(block, "Ready SNPs") == "9,999,999 / 10,000,000 (>99.9999% of ready variants)"
    assert _row(block, "Ready indels / other variants") == "1 / 10,000,000 (<0.0001% of ready variants)"
    assert _row(block, "Ready for harmonisation") == "10,000,000 / 10,000,000 (100.0000% of rows read)"


def test_overview_build_support_uses_testable_denominator_not_all_input_rows():
    opening = _opening(render_dataset_report(*_evidence()))
    block = _text(_block(opening, 3))
    assert "80 / 81" in block
    assert "testable" in block.lower()
    assert "80 / 100" not in block
    assert "automatically inferred" in block and "declared in sample sheet" in block
    assert "odds_ratio" in block and "raw" in block


def test_overview_distinguishes_configured_build_from_inference():
    record, manifest = _evidence()
    record["genome_build_source"] = "policy"
    manifest["dataset"]["genome_build"]["forced"] = True
    block = _block(_opening(render_dataset_report(record, manifest)), 3)
    assert _row(block, "Genome build") == "GRCh37 · configured in YAML"
    assert "automatically inferred" in _row(block, "Effect type")
    assert "declared in sample sheet" in _row(block, "P-value type")


def test_overview_vcf_qc_keeps_build_snp_denominator_and_unfiltered_output_explicit():
    opening = _opening(render_dataset_report(*_evidence()))
    block = _text(_block(opening, 8))
    assert "GRCh37" in block and "84 / 89" in block
    assert "88 / 94" in block
    assert "SNP" in block and "indel" in block.lower()
    assert "unfiltered" in block.lower()
    assert "overlap" in block.lower()
    assert "GRCh38" not in block


def test_overview_liftover_separates_plugin_rejections_from_policy_exclusions():
    block = _text(_block(_opening(render_dataset_report(*_evidence())), 6))
    assert "91" in block and "GRCh38" in block
    assert "rejection" in block.lower()
    assert "successfully lifted" in block.lower() and "swap" in block.lower()
    assert "merged" in block.lower() and "chromosome" in block.lower()
    assert "not" in block.lower()


def test_overview_keeps_eaf_and_info_sources_separate_from_reference_validation():
    opening = _opening(render_dataset_report(*_evidence()))
    strand = _text(_block(opening, 4))
    statistical = _text(_block(opening, 5))
    reference = _text(_block(opening, 7))
    assert "frequency_input" in strand and "alignment_panel" in strand
    assert "quality_input" in statistical
    assert "study_supplied" in strand and "internal column quality_input" in statistical
    assert "comparison_panel" in reference and "0.17" in reference
    assert "ancestry" in reference.lower() and "not" in reference.lower()
    assert "warn" in strand.lower() and "changed" in strand.lower()


def test_population_correlation_and_selected_reference_af_keep_distinct_denominators():
    record, manifest = _evidence()
    manifest["population_frequency_qc"]["populations"] = {
        "EUR": {"comparable_variants": 85, "pearson_correlation": 0.91},
    }
    block = _block(_opening(render_dataset_report(record, manifest)), 7)
    assert _row(block, "Population-correlation comparable variants") == "85"
    assert _row(block, "Reference-population AF correlation") == "0.91"
    assert _row(block, "Concordant AF") == "89 / 90 (98.8889% of comparable variants)"


def test_headline_detector_results_do_not_replace_different_applied_declarations():
    record, manifest = _evidence()
    record.update(
        genome_build_detected="GRCh37", effect_type_detected="beta",
        effect_type_source="sample_sheet", p_value_type_detected="raw",
        p_value_type="neglog10",
    )
    manifest["config"].update(declared_effect_type="odds_ratio", declared_pvalue_type="neglog10")
    manifest["dataset"]["study_decisions"].update(
        effect_type_detected="beta", effect_type_source="sample_sheet",
        pvalue_type="neglog10", pvalue_type_detected="raw",
    )
    opening = _opening(render_dataset_report(record, manifest))
    assert _metric(opening, "Inferred genome build") == "GRCh37"
    assert _metric(opening, "Detected effect type") == "beta"
    assert _metric(opening, "Detected P-value type") == "raw"
    assert _metric(opening, "Most similar reference population") == "EUR"
    assert "Closest AF population" not in opening
    decisions = _block(opening, 3)
    assert _row(decisions, "Effect type") == "odds_ratio · declared in sample sheet"
    assert _row(decisions, "P-value type") == "neglog10 · declared in sample sheet"
    assert "not proof" in _text(opening)
    assert opening.index('class="metric-label">Inferred genome build') < opening.index("<h4>1.1")


@pytest.mark.parametrize("explicit_null", [False, True])
def test_missing_detector_results_do_not_fall_back_to_applied_or_forced_settings(explicit_null):
    record, manifest = _evidence()
    record["genome_build_source"] = "policy"
    manifest["dataset"]["genome_build"].update(forced=True, mode="GRCh37")
    if explicit_null:
        record.update(genome_build_detected=None, effect_type_detected=None, p_value_type_detected=None)
    opening = _opening(render_dataset_report(record, manifest))
    for label in ("Inferred genome build", "Detected effect type", "Detected P-value type"):
        assert _metric(opening, label) == "Not recorded"
    assert _row(_block(opening, 3), "Genome build") == "GRCh37 · configured in YAML"
    assert "odds_ratio" in _row(_block(opening, 3), "Effect type")
    assert "raw" in _row(_block(opening, 3), "P-value type")


def test_all_population_comparisons_are_in_opening_with_separate_missingness_denominators():
    record, manifest = _evidence()
    metrics = {
        "AFR": (0.72, 0.12, 0.8, 4), "EAS": (0.81, 0.09, 0.7, 2),
        "EUR": (0.98, 0.02, 0.9, 0), "SAS": (0.87, 0.07, 0.6, 3),
    }
    population = manifest["population_frequency_qc"]
    population["populations"] = {
        label: {"comparable_variants": 85, "pearson_correlation": correlation,
                "mean_absolute_difference": difference,
                "inverted_mean_absolute_difference": inverted}
        for label, (correlation, difference, inverted, _) in metrics.items()
    }
    population["fields"] = {
        label: {"total_records": 94, "missing": missing, "missing_fraction": missing / 94,
                "invalid": 0, "invalid_fraction": 0.0}
        for label, (_, _, _, missing) in metrics.items()
    }
    html = render_dataset_report(record, manifest)
    reference = _block(_opening(html), 7)
    detail = html.split('id="dataset-overview_study-section-6"', 1)[1].split(
        '<section class="report-section"', 1,
    )[0]
    for label, (correlation, difference, inverted, missing) in metrics.items():
        common = [label, "85", str(correlation), str(difference)]
        assert _population_row(reference, label) == [
            *common, f"{missing} / 94 ({100 * missing / 94:.4f}% of VCF records)",
        ]
        assert _population_row(detail, label) == [*common, str(inverted)]
    assert "Inverted mean absolute AF difference" not in reference
    assert "Inverted mean absolute AF difference" in detail
    assert "Missing percentage" in detail and "Invalid percentage" in detail
    assert "Closest population" not in html
    assert "Closest reference population" not in html
    assert "Most similar reference population" in detail
    assert "it is not confirmation of ancestry" in detail
    assert "common eligible records" in _text(reference)
    assert "missingness uses all its records" in _text(reference)
    assert "not confirmed ancestry" in _text(reference)
    assert "not independent evidence" in _text(reference)
    assert _row(reference, "Concordant AF") == "89 / 90 (98.8889% of comparable variants)"
    del population["fields"]["SAS"]["missing"]
    reference = _block(_opening(render_dataset_report(record, manifest)), 7)
    assert _population_row(reference, "SAS")[-1] == "Not recorded / 94"


@pytest.mark.parametrize("status", [None, "failed", "inconclusive"])
def test_absent_population_comparisons_do_not_become_a_closest_population_or_pass(status):
    record, manifest = _evidence()
    record["closest_population"] = None
    manifest["population_frequency_qc"] = {} if status is None else {"status": status}
    if status == "inconclusive":
        manifest["population_frequency_qc"].update(
            comparable_variants=0, closest_population=None,
            decision_reason="insufficient comparable population-frequency evidence",
            populations={"EUR": {"comparable_variants": 0, "pearson_correlation": None,
                                  "mean_absolute_difference": None}},
            fields={"EUR": {"total_records": 94, "missing": 94}},
        )
    opening = _opening(render_dataset_report(record, manifest))
    reference = _block(opening, 7)
    assert _metric(opening, "Most similar reference population") == "Not recorded"
    assert _row(reference, "Most similar reference population") == "Not recorded"
    assert _row(reference, "Population comparison status") == (status or "Not recorded")
    if status == "inconclusive":
        assert _population_row(reference, "EUR") == [
            "EUR", "0", "Not recorded", "Not recorded", "94 / 94 (100.0000% of VCF records)",
        ]
        assert _row(reference, "Population-correlation comparable variants") == "0"
        assert "insufficient comparable population-frequency evidence" in reference
    else:
        assert "Population-frequency comparisons not recorded" in _text(reference)
        assert "not a passing similarity result" in _text(reference)


def test_arbitrary_population_labels_and_headline_values_are_html_escaped():
    record, manifest = _evidence()
    unsafe = '<img src=x onerror="alert(1)">'
    record.update(closest_population=unsafe, effect_type_detected=unsafe)
    labels = ("Panel-27", unsafe)
    manifest["population_frequency_qc"].update(
        closest_population=unsafe,
        populations={label: {"comparable_variants": 85, "pearson_correlation": 0.91,
                             "mean_absolute_difference": 0.03} for label in labels},
    )
    opening = _opening(render_dataset_report(record, manifest))
    reference = _block(opening, 7)
    for label in labels:
        assert _population_row(reference, label)[:4] == [label, "85", "0.91", "0.03"]
    assert _metric(opening, "Most similar reference population") == unsafe
    assert _metric(opening, "Detected effect type") == unsafe
    assert "<img" not in opening
    assert "&lt;img" in reference and "&quot;alert(1)&quot;" in reference


def test_run_report_keeps_headline_detections_and_population_metrics_per_dataset():
    first, first_manifest = _evidence("first")
    second, second_manifest = _evidence("second")
    entries = ((first, first_manifest, "GRCh37", "beta", "raw", "reference_one", 0.91),
               (second, second_manifest, "GRCh38", "odds_ratio", "neglog10", "reference_two", 0.87))
    for record, manifest, build, effect, pvalue, population, correlation in entries:
        record.update(genome_build_detected=build, effect_type_detected=effect,
                      p_value_type_detected=pvalue, closest_population=population)
        manifest["population_frequency_qc"].update(
            closest_population=population,
            populations={population: {"comparable_variants": 85, "pearson_correlation": correlation,
                                      "mean_absolute_difference": 0.03}},
        )
    html = render_run_report([first, second], {"first": first_manifest, "second": second_manifest})
    for record, _, build, effect, pvalue, population, correlation in entries:
        opening = _opening(html, record["dataset_id"])
        assert _metric(opening, "Inferred genome build") == build
        assert _metric(opening, "Detected effect type") == effect
        assert _metric(opening, "Detected P-value type") == pvalue
        assert _metric(opening, "Most similar reference population") == population
        assert _population_row(_block(opening, 7), population)[:4] == [population, "85", str(correlation), "0.03"]
        other_population = "reference_two" if population == "reference_one" else "reference_one"
        assert other_population not in opening


@pytest.mark.parametrize("kind", ["external", "fixed"])
def test_overview_distinguishes_external_info_proxy_from_user_assigned_constant(kind):
    record, manifest = _evidence()
    manifest["config"]["imp_info_col"] = None
    manifest["config"]["info_source"] = "external" if kind == "external" else "fixed_cli"
    if kind == "external":
        manifest["config"]["info_source_detail"] = "external proxy quality_ref from /reference/quality.tsv"
    else:
        manifest["config"]["info_source_detail"] = "user-assigned fixed value"
        manifest["config"]["fixed_info"] = 0.83
    block = _block(_opening(render_dataset_report(record, manifest)), 5)
    if kind == "external":
        assert "external proxy quality_ref" in _row(block, "INFO source")
        assert _row(block, "Fixed INFO assigned by user") == "Not used"
    else:
        assert _row(block, "INFO source") == "user-assigned fixed value"
        assert _row(block, "Fixed INFO assigned by user") == "0.83"


def test_known_internal_info_does_not_report_an_unused_fixed_value_as_missing_evidence():
    record, manifest = _evidence()
    manifest["config"]["info_source"] = "internal"
    manifest["config"]["fixed_info"] = None
    block = _block(_opening(render_dataset_report(record, manifest)), 5)
    assert _row(block, "INFO source") == "internal column quality_input"
    assert _row(block, "Fixed INFO assigned by user") == "Not used"


def test_absent_info_source_does_not_invent_internal_external_or_fixed_provenance():
    opening = _opening(render_dataset_report({"dataset_id": "overview_study"}, {}))
    block = _block(opening, 5)
    assert _row(block, "INFO source") == "Not recorded"
    assert _row(block, "Fixed INFO assigned by user") == "Not recorded"
    assert _row(block, "INFO rounding corrections") == "Not recorded"


def test_external_eaf_uses_external_column_and_preserves_reference_template():
    record, manifest = _evidence()
    manifest["config"].update(eaf_col=None, eafcolumn="EUR_reference",
                              eaffile="/reference/panel_chr{chromosome}.tsv.gz")
    manifest["dataset"]["study_decisions"]["eaf_provenance"] = "reference_imputed"
    block = _block(_opening(render_dataset_report(record, manifest)), 4)
    assert _row(block, "EAF source") == "reference_imputed · EUR_reference"
    assert _row(block, "External EAF column") == "EUR_reference"
    assert "panel_chr{chromosome}.tsv.gz" in _row(block, "External EAF file / template")


def test_partial_chromosome_counts_are_labelled_at_the_count_not_only_in_status():
    record, manifest = _evidence()
    record.update(status="PARTIAL", strand_chromosomes_expected=2)
    manifest["dataset"]["failed"] = ["2"]
    manifest["dataset"]["chromosomes"]["2"] = {"status": "failed"}
    opening = _opening(render_dataset_report(record, manifest))
    statistic = _row(_block(opening, 5), "INFO rounding corrections")
    assert statistic == "2 — partial evidence from 1 / 2 chromosomes"
    assert "91 — partial evidence from 1 / 2 chromosomes" in _text(_block(opening, 6))
    assert "Incomplete dataset" in _block(opening, 10)


def test_concordance_failure_retains_variant_union_and_statistic_specific_denominators():
    record, manifest = _evidence()
    manifest["concordance_validation"] = {
        "status": "failed",
        "summary": {"variant_types": {"snps": {
            "input_unique_variants": 5, "vcf_unique_variants": 4,
            "exact_matched_variants": 3, "variant_union": 6,
            "input_only_variants": 2, "vcf_only_variants": 1,
        }}},
        "metrics": {"effect": {"concordant": 1, "checked": 2, "mismatches": 1,
                                "unavailable_in_input": 0, "missing_in_vcf": 0,
                                "excluded_orientation": 1}},
    }
    block = _text(_block(_opening(render_dataset_report(record, manifest)), 9))
    assert "Failed" in block
    assert "3 / 6 (50.0000% of input/VCF variant union)" in block
    assert "1 / 2 (50.0000% of checked values)" in block
    assert "Orientation excluded" in block
    assert "do not by themselves prove a statistical error" in block


@pytest.mark.parametrize("status", ["PARTIAL", "FAILED", "INTERRUPTED", "PREFLIGHT_FAILED", "RUNNING", "NOT_RUN"])
def test_missing_or_incomplete_run_evidence_does_not_become_success(status):
    record = {"dataset_id": "overview_study", "status": status}
    opening = _opening(render_dataset_report(record, {}))
    assert re.findall(r"<h4>(1\.\d+ [^<]*)</h4>", opening) == list(_HEADINGS)
    assert "Not recorded" in opening
    assert "100.0000%" not in opening
    assert "ready for downstream analysis" not in _text(opening).lower()
    takeaway = _text(_block(opening, 10))
    assert "No headline alert" not in takeaway
    assert "Incomplete dataset" in takeaway if status == "PARTIAL" else "Completion not established" in takeaway


def test_overview_escapes_untrusted_sources_and_failure_messages():
    record, manifest = _evidence()
    unsafe = '<img src=x onerror="alert(1)">'
    record["failure_reason"] = unsafe
    manifest["vcf_provenance"]["af_input_column"] = unsafe
    manifest["config"]["eaf_col"] = unsafe
    opening = _opening(render_dataset_report(record, manifest))
    assert "<img" not in opening
    assert "&lt;img" in opening
    assert "&quot;alert(1)&quot;" in opening
    assert "&lt;img" in _block(opening, 4)
