"""Regression coverage for the scientific fine-mapping HTML report."""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
from html.parser import HTMLParser
from pathlib import Path

import pandas as pd
import pytest

from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.modules.fine_mapping.output_layout import resolve_output_paths
from postgwas.modules.fine_mapping.preflight import (
    FineMappingInputEvidence,
    FineMappingPreflight,
    FineMappingReferenceEvidence,
    FineMappingToolEvidence,
)
from postgwas.modules.fine_mapping.reporting import (
    FineMappingReportError,
    build_fine_mapping_report,
    render_fine_mapping_html,
    write_fine_mapping_html_report,
)


def _preflight(tmp_path: Path) -> FineMappingPreflight:
    summary = tmp_path / "study_susie.tsv"
    loci = tmp_path / "loci.tsv"
    reference = tmp_path / "reference"
    plink = tmp_path / "plink"
    for path in (summary, loci, plink):
        path.write_text("fixture\n", encoding="utf-8")
    for suffix in (".bed", ".bim", ".fam"):
        reference.with_suffix(suffix).write_text("fixture\n", encoding="utf-8")
    return FineMappingPreflight(
        engine="susie",
        genome_build="GRCh37",
        inputs=FineMappingInputEvidence(
            summary_statistics_file=summary,
            summary_statistics_rows=100,
            summary_statistics_columns=("SNP", "CHR", "BP", "REF", "ALT", "EZ", "NEF", "LP"),
            summary_statistics_chromosomes=("1",),
            locus_summary_variants=3,
            locus_file=loci,
            input_loci=2,
            eligible_loci=2,
            loci_below_lp_threshold=0,
            loci_excluded_by_mhc=0,
            locus_chromosomes=("1",),
        ),
        reference=FineMappingReferenceEvidence(
            prefix=str(reference),
            bed_file=reference.with_suffix(".bed"),
            bim_file=reference.with_suffix(".bim"),
            fam_file=reference.with_suffix(".fam"),
            variants=500,
            samples=250,
            chromosomes=("1",),
            bed_bytes=1000,
            locus_variant_id_matches=3,
            locus_variant_ids_requested=3,
            coordinate_concordant_variant_matches=3,
            loci_with_reference_variants=2,
            loci_requested=2,
        ),
        tools=(FineMappingToolEvidence("PLINK", str(plink), "PLINK v1.90"),),
        r_runtime={"version": "R 4.4; susieR available"},
    )


def _report_fixture(tmp_path: Path, *, empty: bool = False):
    module = load_configuration().modules.fine_mapping
    output = tmp_path / "fine_mapping"
    paths = resolve_output_paths(
        output, module.output_layout.model_dump(), "study",
    )
    paths["combined_results_directory"].mkdir(parents=True)
    paths["quality_control_directory"].mkdir(parents=True)
    final_set = paths["downstream_flames_directory"] / "primary__00001_set.txt"
    if not empty:
        final_set.parent.mkdir(parents=True)
        final_set.write_text(
            "index\tcred1\tprob1\n"
            "1\t1:100:A_G\t0.70\n"
            "2\t1:120:C_T\t0.26\n",
            encoding="utf-8",
        )
        final_rows = [{
            "Filename": "source.txt",
            "GenomicLocus": "chr1:90-130",
            "genome_build": "GRCh37",
            "credible_set": "L1",
            "component_index": 1,
            "target_coverage": 0.95,
            "achieved_component_coverage": 0.951,
            "global_pip_sum": 0.96,
            "credible_set_membership_definition": "susie_component_alpha_configured_coverage",
            "probability_definition": "susie_model_wide_pip",
            "n_variants": 2,
            "component_log10bf": 4.2,
            "min_abs_ld": 0.72,
            "analysis_round": "primary",
            "final_selection_reason": "non_overlapping_primary_result",
            "overlap_group_id": "",
            "final_genomic_locus": "chr1:90-130",
            "source_primary_genomic_loci": "chr1:90-130",
            "source_primary_credible_set_files": "source.txt",
            "source_primary_status_files": "locus.tsv",
            "source_primary_manifest_files": "manifest.tsv",
            "source_primary_warning_reasons": "",
            "source_primary_failure_reasons": "",
            "source_result_file": "manifest.tsv",
            "final_credible_set_file": str(final_set),
            "warning_reason": "nef_relative_range_exceeds_threshold",
            "failure_reason": "",
        }]
    else:
        final_rows = []
    provenance_columns = [
        "analysis_round", "final_selection_reason", "overlap_group_id",
        "final_genomic_locus", "source_primary_genomic_loci",
        "source_primary_credible_set_files", "source_primary_status_files",
        "source_primary_manifest_files", "source_primary_warning_reasons",
        "source_primary_failure_reasons", "source_result_file",
        "final_credible_set_file", "warning_reason", "failure_reason",
    ]
    final_frame = pd.DataFrame(final_rows)
    if final_frame.empty:
        final_frame = pd.DataFrame(columns=provenance_columns)
    final_path = paths["combined_results_directory"] / module.overlap_resolution.final_combined_filename
    final_frame.to_csv(final_path, sep="\t", index=False)

    primary_index = paths["primary_flames_directory"] / module.overlap_resolution.index_filename
    if not empty:
        primary_index.parent.mkdir(parents=True)
        pd.DataFrame([{
            "Filename": str(final_set),
            "GenomicLocus": "chr1:90-130",
            "Annotfiles": "unused_annotation.txt",
        }]).to_csv(primary_index, sep="\t", index=False)

    locus_status = paths["susie_qc_file"]
    pd.DataFrame([
        {
            "genomic_locus": "chr1:90-130",
            "stage": "final",
            "converged": True,
            "warning_reason": "nef_relative_range_exceeds_threshold",
            "failure_reason": "",
            "input_variant_count": 3,
            "selected_nef": 1000,
            "nef_relative_range": 0.10,
            "ld_z_mismatch_lambda": 0.002,
            "ld_min_eigenvalue": 0.01,
        },
        {
            "genomic_locus": "chr1:500-700",
            "stage": "failed",
            "converged": False,
            "warning_reason": "",
            "failure_reason": "severe_ld_z_mismatch",
            "input_variant_count": 20,
        },
    ]).to_csv(locus_status, sep="\t", index=False)
    preflight_path = paths["preflight_validation_file"]
    pd.DataFrame([
        {
            "category": "input",
            "check": "summary_statistics",
            "status": "passed",
            "value": 100,
            "path": str(tmp_path / "study_susie.tsv"),
            "detail": "required columns and row values validated",
        }
    ]).to_csv(preflight_path, sep="\t", index=False)
    overlap_path = paths["quality_control_directory"] / module.overlap_resolution.plan_filename
    pd.DataFrame([
        {
            "overlap_group_id": "",
            "primary_genomic_locus": "chr1:90-130",
            "chromosome": 1,
            "primary_start": 90,
            "primary_end": 130,
            "joint_genomic_locus": "",
            "resolution_status": "primary_retained",
            "resolution_reason": "no_primary_result_overlap",
            "primary_warning_reason": "nef_relative_range_exceeds_threshold",
            "joint_warning_reason": "",
            "joint_failure_reason": "",
        }
    ]).to_csv(overlap_path, sep="\t", index=False)
    overlap_summary = paths["quality_control_directory"] / module.overlap_resolution.summary_filename
    pd.DataFrame([{"status": "success"}]).to_csv(
        overlap_summary, sep="\t", index=False,
    )
    recovery = paths["quality_control_directory"] / module.engines.susie.recovery_audit_filename
    pd.DataFrame([
        {
            "genomic_locus": "chr1:90-130",
            "sequence": 1,
            "stage": "primary_fit",
            "action": "fit",
            "status": "success",
            "reason": "primary_fit_converged",
            "ld_status": "valid",
        }
    ]).to_csv(recovery, sep="\t", index=False)
    args = argparse.Namespace(
        dataset_id="study",
        finemap_method="susie",
        genome_build="GRCh37",
        fine_mapping_html_report=module.html_report.model_dump(),
        resolved_fine_mapping_configuration=module.model_dump(mode="json"),
        fine_mapping_validation=module.validation.model_dump(),
    )
    result = {
        "status": "completed_no_credible_sets" if empty else "success",
        "output_dir": str(output),
        "flames_input": None if empty else str(final_set.parent),
        "flames_index": None if empty else str(primary_index),
        "n_credible_sets": 0 if empty else 1,
        "locus_status": str(locus_status),
        "preflight_validation": str(preflight_path),
        "overlap_resolution": str(overlap_path),
        "overlap_resolution_summary": str(overlap_summary),
        "final_combined_credible_sets": str(final_path),
        "recovery_audit": str(recovery),
        "n_attempted": 2,
        "n_successful": 1,
        "n_failed": 1,
        "n_warnings": 1,
        "warning_reason_counts": {"nef_relative_range_exceeds_threshold": 1},
        "failure_reason_counts": {"severe_ld_z_mismatch": 1},
        "n_overlap_groups": 0,
        "n_joint_rerun_successful": 0,
        "n_joint_rerun_failed": 0,
        "n_joint_rerun_without_credible_sets": 0,
        "n_final_credible_sets": 0 if empty else 1,
    }
    return args, _preflight(tmp_path), result, paths


def test_report_contains_complete_results_qc_provenance_and_interpretation(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)

    destination = write_fine_mapping_html_report(
        args, preflight, result, paths,
    )

    assert destination == paths["html_report_file"]
    html = destination.read_text(encoding="utf-8")
    assert "Fine-mapping report" in html
    assert "SuSiE-RSS" in html
    assert "Credible sets and variant details" in html
    assert '<section id="variants">' not in html
    assert "Complete credible-set details" not in html
    assert "Initial fine-mapping: locus details" in html
    assert "Additional analysis: overlapping loci" in html
    assert "Input and resource validation" in html
    assert "SuSiE fitting and recovery audit" in html
    assert "model-wide posterior inclusion probability" in html
    assert "not by itself establish biological causality" in html
    assert "1:100:A_G" in html
    assert "nef_relative_range_exceeds_threshold" in html
    assert "severe_ld_z_mismatch" in html
    assert "authoritative" not in html.lower()
    assert "Analysis completed — review warnings or exclusions" in html
    assert "Warnings and failed or skipped loci remain part of the interpretation" in html
    assert "Reference coordinate matches" in html
    assert "Loci with reference variants" in html
    assert 'class="data-browser"' in html
    assert 'data-page-size="50"' in html
    assert "Search all rows" in html
    assert "study_fine_mapping_report.html" in html
    assert "PLINK BED" in html
    assert "PLINK BIM" in html
    assert "PLINK FAM" in html
    assert "Diagnostic plots" not in html


def test_report_preserves_valid_completed_outcome_without_credible_sets(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path, empty=True)

    destination = write_fine_mapping_html_report(
        args, preflight, result, paths,
    )

    html = destination.read_text(encoding="utf-8")
    assert "complete · no credible sets" in html
    assert "No credible set was retained" in html
    assert "not evidence that the locus contains no causal variant" in html


def test_finemap_report_keeps_snp_and_model_probabilities_distinct(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    args.finemap_method = "finemap"
    args.resolved_fine_mapping_configuration = deepcopy(
        args.resolved_fine_mapping_configuration
    )
    args.resolved_fine_mapping_configuration["engine"] = "finemap"
    preflight = replace(
        preflight,
        engine="finemap",
        tools=(
            FineMappingToolEvidence("PLINK 2", "/tools/plink2", "PLINK v2"),
            FineMappingToolEvidence("FINEMAP", "/tools/finemap", "FINEMAP 1.4.2"),
        ),
        r_runtime=None,
    )
    final = pd.read_csv(result["final_combined_credible_sets"], sep="\t")
    final["engine"] = "FINEMAP"
    final["model_k"] = 2
    final["model_posterior_probability"] = 0.81
    final["achieved_coverage"] = 0.96
    final["credible_set_membership_definition"] = (
        "finemap_0.95_coverage_credible_set"
    )
    final["probability_definition"] = "finemap_snp_posterior_in_selected_model"
    final.to_csv(result["final_combined_credible_sets"], sep="\t", index=False)
    paths["finemap_model_summary_file"].parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {
            "Locus": "chr1:90-130",
            "File_Name": "chr1:90-130.cred2",
            "K_Causal": 2,
            "Posterior_Prob": 0.81,
            "Selected": True,
            "Validation_Status": "valid",
            "Failure_Reason": "",
            "Failure_Detail": "",
        }
    ]).to_csv(paths["finemap_model_summary_file"], index=False)

    destination = write_fine_mapping_html_report(
        args, preflight, result, paths,
    )

    html = destination.read_text(encoding="utf-8")
    assert "SNP posterior probability within the selected FINEMAP model" in html
    assert "causal-count model probability are different quantities" in html
    assert "FINEMAP causal-count model selection" in html
    assert "Posterior Prob" in html
    assert "0.81" in html


def test_report_rejects_invalid_authoritative_posterior_probability(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    final = pd.read_csv(result["final_combined_credible_sets"], sep="\t")
    credible_file = Path(final.loc[0, "final_credible_set_file"])
    credible_file.write_text(
        "index\tcred1\tprob1\n1\t1:100:A_G\t1.2\n",
        encoding="utf-8",
    )

    with pytest.raises(FineMappingReportError, match=r"outside \[0, 1\]"):
        build_fine_mapping_report(args, preflight, result, paths)


def test_report_rejects_configured_columns_that_match_no_audit_fields(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    args.fine_mapping_html_report["susie_recovery_columns"] = ["not_a_column"]

    with pytest.raises(FineMappingReportError, match="do not match"):
        build_fine_mapping_report(args, preflight, result, paths)


def test_html_report_configuration_is_schema_validated(tmp_path):
    invalid = tmp_path / "fine_mapping.yaml"
    invalid.write_text(
        "html_report:\n  page_size: 0\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="page_size"):
        from postgwas.config import load_run_configuration_for_module

        load_run_configuration_for_module("fine_mapping", invalid)


def test_initial_counts_remain_distinct_when_loci_are_reanalysed_together(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    # Two initially analysed loci produce three sets; their joint analysis
    # produces one final set. Final counts must never replace initial counts.
    index = pd.read_csv(result["flames_index"], sep="\t")
    extra = pd.DataFrame([
        {"Filename": "initial_set_2.txt", "GenomicLocus": "chr1:90-130", "Annotfiles": "a.txt"},
        {"Filename": "initial_set_3.txt", "GenomicLocus": "chr1:120-180", "Annotfiles": "b.txt"},
    ])
    pd.concat([index, extra], ignore_index=True).to_csv(result["flames_index"], sep="\t", index=False)
    result.update(n_credible_sets=3, n_attempted=3, n_successful=2, n_failed=1,
                  n_overlap_groups=1, n_joint_rerun_successful=1)
    report = build_fine_mapping_report(args, preflight, result, paths)
    summary = report["summary"]
    assert summary["initial_credible_sets"] == 3
    assert summary["initial_loci_with_credible_sets"] == 2
    assert summary["initial_loci_without_credible_sets"] == 0
    assert summary["final_credible_sets"] == 1
    html = render_fine_mapping_html(report)
    overview = html.split('<section id="overview">')[1].split('</section>')[0]
    assert 'Credible sets found</div><div class="metric-value">3</div>' in overview
    assert 'Loci with credible sets</div><div class="metric-value">2</div>' in overview
    assert "Final credible sets" not in overview
    assert "Set variant memberships" not in overview
    assert "Reference ID matches" not in overview
    assert "Highest posterior value" not in overview
    assert html.index('<section id="loci">') < html.index('<section id="overlap">')
    assert html.index('<section id="overlap">') < html.index('<section id="credible-sets">')
    assert "analysed again as one larger region" in html
    assert "without adding another flank" in html
    assert "their initial sets are not restored" in html
    assert "the new sets replace the initial sets" in html
    assert "Groups without a valid new set are excluded" in html
    assert "The results shown below are unchanged" not in html


def test_no_overlap_explains_why_no_additional_fine_mapping_was_needed(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    html = render_fine_mapping_html(build_fine_mapping_report(args, preflight, result, paths))
    assert "No additional fine-mapping was needed" in html
    assert "initial credible sets were kept unchanged" in html
    assert "Overlapping regions can contain the same association signal" in html
    assert "Groups producing credible sets" not in html
    assert "The initial analysis found 1 credible set." in html
    assert "The results shown below are unchanged from the initial analysis" in html
    assert "A longer bar means stronger support for that variant within this model" in html
    assert "These are the sets retained after the overlap check" not in html


def test_initial_no_set_count_excludes_failed_loci(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path, empty=True)
    report = build_fine_mapping_report(args, preflight, result, paths)
    assert report["summary"]["initial_loci_without_credible_sets"] == 1
    assert report["summary"]["initial_loci_with_credible_sets"] == 0
    assert report["summary"]["primary_loci_failed"] == 1
    html = render_fine_mapping_html(report)
    assert "No credible sets are available in the final results" in html
    assert "The results shown below are unchanged" not in html


@pytest.mark.parametrize("problem", ["missing_count", "missing_index", "wrong_count", "too_many_loci", "duplicate_file", "missing_locus"])
def test_report_rejects_unverifiable_initial_results(tmp_path, problem):
    args, preflight, result, paths = _report_fixture(tmp_path)
    if problem == "missing_count":
        result.pop("n_credible_sets")
    elif problem == "missing_index":
        Path(result["flames_index"]).unlink()
    elif problem == "wrong_count":
        result["n_credible_sets"] = 2
    elif problem == "too_many_loci":
        result["n_successful"] = 0
    else:
        index = pd.read_csv(result["flames_index"], sep="\t")
        if problem == "duplicate_file":
            index = pd.concat([index, index], ignore_index=True)
            result["n_credible_sets"] = 2
        else:
            index.loc[0, "GenomicLocus"] = None
        index.to_csv(result["flames_index"], sep="\t", index=False)
    with pytest.raises(FineMappingReportError, match="Initial"):
        build_fine_mapping_report(args, preflight, result, paths)


@pytest.mark.parametrize("status", ["completed_with_excluded_overlap_groups", "completed_no_credible_sets", "failed_no_authoritative_credible_sets"])
def test_incomplete_outcomes_never_receive_a_clean_completion_label(tmp_path, status):
    args, preflight, result, paths = _report_fixture(tmp_path, empty=True)
    result.update(status=status, n_overlap_groups=1, n_joint_rerun_failed=1)
    html = render_fine_mapping_html(build_fine_mapping_report(args, preflight, result, paths))
    assert 'badge good' not in html
    assert "authoritative" not in html.lower()
    assert "Groups failed or excluded" in html


def test_warning_explanations_use_affected_loci_and_resolved_thresholds(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    qc = pd.read_csv(result["locus_status"], sep="\t")
    combined = "nef_relative_range_exceeds_threshold;ld_z_mismatch_warning"
    qc.loc[0, ["warning_reason", "nef_min", "nef_max", "nef_relative_range", "selected_nef", "ld_z_mismatch_lambda"]] = [combined, 700, 1000, 0.3, 1000, 0.025]
    # An unrelated failed locus must not inflate the observed warning ranges.
    qc.loc[1, ["nef_relative_range", "ld_z_mismatch_lambda"]] = [9, 0.99]
    qc.to_csv(result["locus_status"], sep="\t", index=False)
    result["warning_reason_counts"] = {combined: 1}
    config = args.resolved_fine_mapping_configuration
    config["sample_size"]["relative_range_warning_threshold"] = 0.2
    config["engines"]["susie"]["ld_validation"].update(
        mismatch_warning_threshold=0.02, mismatch_failure_threshold=0.15,
    )
    report = build_fine_mapping_report(args, preflight, result, paths)
    warnings = report["warning_explanations"]
    assert len(warnings) == 2
    assert [row["loci"] for row in warnings] == [1, 1]
    assert [row["Value"] for row in warnings[0]["measurements"]] == ["30%", "20%", "700–1,000", "1,000"]
    assert [row["Value"] for row in warnings[1]["measurements"]] == ["0.025", 0.02, 0.15]
    html = render_fine_mapping_html(report)
    assert "Effective sample size differs between variants" in html
    assert "GWAS statistics and reference LD do not fully agree" in html
    assert "maximum NEF − minimum NEF" in html
    assert "not guarantees of accuracy" in html
    assert "does not identify the cause" in html
    assert "genomic inflation λ" in html
    assert "affected-locus counts should not be added together" in html
    assert combined in html
    assert "https://pmc.ncbi.nlm.nih.gov/articles/PMC9337707/" in html
    assert "https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html" in html


@pytest.mark.parametrize("value", [None, "invalid", float("inf")])
def test_missing_warning_diagnostics_are_not_reported_as_zero(tmp_path, value):
    args, preflight, result, paths = _report_fixture(tmp_path)
    qc = pd.read_csv(result["locus_status"], sep="\t")
    qc.loc[0, "nef_relative_range"] = value
    qc.to_csv(result["locus_status"], sep="\t", index=False)
    report = build_fine_mapping_report(args, preflight, result, paths)
    assert report["warning_explanations"][0]["measurements"][0]["Value"] == "Not fully recorded in locus QC"


def test_unknown_warning_codes_remain_visible_without_invented_explanations(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    code = "unrecognised_warning_<check>"
    qc = pd.read_csv(result["locus_status"], sep="\t")
    qc["warning_reason"] = code
    qc.to_csv(result["locus_status"], sep="\t", index=False)
    result["warning_reason_counts"] = {code: len(qc)}
    report = build_fine_mapping_report(args, preflight, result, paths)
    assert report["warning_explanations"] == []
    html = render_fine_mapping_html(report)
    assert "unrecognised_warning_&lt;check&gt;" in html
    assert "All recorded warning and failure codes" in html


def test_clean_loci_do_not_receive_scientific_warning_explanations(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    qc = pd.read_csv(result["locus_status"], sep="\t")
    qc["warning_reason"] = None
    qc.to_csv(result["locus_status"], sep="\t", index=False)
    result["warning_reason_counts"] = {}
    assert build_fine_mapping_report(args, preflight, result, paths)["warning_explanations"] == []


class _DisclosureParser(HTMLParser):
    """Inspect nesting and literal content without requiring a browser package."""

    def __init__(self, html):
        super().__init__()
        self.nodes = []
        self.stack = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        if tag == "details":
            node = {"attrs": dict(attrs), "text": "", "parent": self.stack[-1] if self.stack else None}
            self.nodes.append(node)
            self.stack.append(node)

    def handle_endtag(self, tag):
        if tag == "details":
            self.stack.pop()

    def handle_data(self, data):
        for node in self.stack:
            node["text"] += data


def test_combined_browser_preserves_each_set_and_all_its_members(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    final = pd.read_csv(result["final_combined_credible_sets"], sep="\t")
    second = final.iloc[0].copy()
    second_file = Path(final.loc[0, "final_credible_set_file"]).with_name("second.txt")
    second_file.write_text("index\tcred1\tprob1\n1\t1:100:A_G\t1\n2\t2:200:C_T\t0\n")
    second["final_credible_set_file"] = str(second_file)
    second["final_genomic_locus"] = "chr2:190-210"
    second["analysis_round"] = "joint"
    # Same component ID and shared variant must stay attached to their own set.
    pd.concat([final, second.to_frame().T], ignore_index=True).to_csv(
        result["final_combined_credible_sets"], sep="\t", index=False,
    )
    result["n_final_credible_sets"] = 2
    report = build_fine_mapping_report(args, preflight, result, paths)
    html = render_fine_mapping_html(report)
    nodes = _DisclosureParser(html).nodes
    cards = [n for n in nodes if n["attrs"].get("class") == "cs-card"]
    variants = [n for n in nodes if n["attrs"].get("class") == "variant-detail"]
    assert len(cards) == 2
    assert len(variants) == 4
    assert all("open" not in n["attrs"] and "hidden" not in n["attrs"] for n in cards + variants)
    assert [v["parent"] is cards[0] for v in variants] == [True, True, False, False]
    for card, row in zip(cards, report["credible_sets"]):
        assert row["final_genomic_locus"] in card["text"]
        assert row["final_credible_set_file"] in card["text"]
        assert "nef_relative_range_exceeds_threshold" in card["text"]
    assert "Initial fine-mapping" in cards[0]["text"]
    assert "Additional fine-mapping" in cards[1]["text"]
    assert "1:120:C_T" not in cards[1]["text"]
    assert "2:200:C_T" not in cards[0]["text"]
    assert all("Source file" in v["text"] and "Rank" in v["text"] for v in variants)
    assert 'style="width:0%"' in html
    assert 'style="width:100%"' in html
    assert "variant-table-data" not in html
    assert "credible-set-table-data" not in html
    assert "Search sets and variants" in html


def test_variant_details_preserve_low_probabilities_and_escape_identifiers(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    final = pd.read_csv(result["final_combined_credible_sets"], sep="\t")
    source = Path(final.loc[0, "final_credible_set_file"])
    source.write_text("index\tcred1\tprob1\n1\t<script>alert(1)</script>\t0.96\n2\tlow_variant\t0.0000001\n")
    html = render_fine_mapping_html(build_fine_mapping_report(args, preflight, result, paths))
    variants = [n for n in _DisclosureParser(html).nodes if n["attrs"].get("class") == "variant-detail"]
    assert len(variants) == 2
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html
    assert "low_variant" in variants[1]["text"]
    assert "1e-07" in variants[1]["text"]


def test_all_members_are_displayed_without_a_top_variant_limit(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    final = pd.read_csv(result["final_combined_credible_sets"], sep="\t")
    count = 60
    source = Path(final.loc[0, "final_credible_set_file"])
    source.write_text("index\tcred1\tprob1\n" + "".join(
        f"{i}\tvariant_{i}\t{1 / count}\n" for i in range(1, count + 1)
    ))
    final["n_variants"] = count
    final.to_csv(result["final_combined_credible_sets"], sep="\t", index=False)
    report = build_fine_mapping_report(args, preflight, result, paths)
    html = render_fine_mapping_html(report)
    variants = [n for n in _DisclosureParser(html).nodes if n["attrs"].get("class") == "variant-detail"]
    assert len(variants) == count
    assert "variant_60" in variants[-1]["text"]
    assert "plotted_variants_per_credible_set" not in args.fine_mapping_html_report


def test_retired_plot_limit_is_rejected_in_custom_configuration(tmp_path):
    from postgwas.config import load_run_configuration_for_module

    config = tmp_path / "old_report.yaml"
    config.write_text("html_report:\n  plotted_variants_per_credible_set: 12\n")
    with pytest.raises(ConfigurationError, match="plotted_variants_per_credible_set"):
        load_run_configuration_for_module("fine_mapping", config)


def test_report_links_recorded_upstream_outputs_once_and_marks_missing_files(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    clumping = tmp_path / "clumping & QC.html"
    formatter = tmp_path / "formatter.html"
    missing = tmp_path / "missing.html"
    for path in (clumping, formatter):
        path.write_text("<!doctype html><title>Earlier result</title>")
    upstream = {
        "ld_clump": {"html_report": str(clumping)},
        "formatter": {"susie": {"html_report": str(formatter)},
                      "finemap": {"html_report": str(formatter)}},
        "another_module": [{"html_report": str(missing)}],
    }
    destination = write_fine_mapping_html_report(
        args, preflight, result, paths, upstream_results=upstream,
    )
    html = destination.read_text()
    section = html.split('<section id="upstream-reports">')[1].split('</section>')[0]
    assert "LD clumping" in section
    assert "Summary-statistics formatting" in section
    assert section.count('formatter.html</a>') == 1
    assert "clumping &amp; QC.html" in section
    assert 'href="../../clumping &amp; QC.html"' in section
    assert "Recorded report file not found" in section
    assert "missing.html</a>" not in section
    assert '<a href="#upstream-reports">Clumping and input details</a>' in html
    assert upstream["formatter"]["finemap"]["html_report"] == str(formatter)


def test_direct_report_does_not_guess_earlier_reports_from_nearby_files(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    (tmp_path / "unrelated_clumping_report.html").write_text("unrelated run")
    report = build_fine_mapping_report(args, preflight, result, paths)
    assert report["upstream_reports"] == []
    html = render_fine_mapping_html(report)
    assert 'href="#upstream-reports"' not in html
    assert '<section id="upstream-reports">' not in html


def test_upstream_report_collection_rejects_invalid_paths_and_accepts_absent_reports():
    from postgwas.core.io.reports import collect_html_reports

    assert collect_html_reports(None) == []
    assert collect_html_reports({"module": {"html_report": None, "log_file": "log.txt"}}) == []
    with pytest.raises(ValueError, match="html_report must be a file path"):
        collect_html_reports({"module": {"html_report": 42}})


def test_clumping_details_are_embedded_open_before_fine_mapping_loci(tmp_path):
    args, preflight, result, paths = _report_fixture(tmp_path)
    source = tmp_path / "clumping.html"
    source.write_text('''<html><head><style>body{display:none}</style></head><body>
    <main><section><h2>Genomic Risk Loci</h2><p>Boundary is not LD membership.</p>
    <div class="result-tools"><label>Search</label><input/><span>1 row</span></div>
    <table id="source-table"><thead><tr><th>Locus</th><th>P value</th></tr></thead>
    <tbody><tr><td>chr1:90-130</td><td>1e-30</td></tr></tbody></table>
    <script>alert("not copied")</script></section>
    <section><h2>Lead SNPs</h2><p onclick="bad()">rs123 &amp; rs456</p>
    <a href="javascript:bad()">Source name</a></section></main></body></html>''')
    report = build_fine_mapping_report(
        args, preflight, result, paths,
        upstream_results={"ld_clump": {"html_report": str(source)}},
    )
    html = render_fine_mapping_html(report)
    assert '<details class="upstream-detail" open>' in html
    assert '<h2>Genomic Risk Loci</h2>' in html
    assert '<h2>Lead SNPs</h2>' in html
    assert '<td>1e-30</td>' in html
    assert 'Boundary is not LD membership.' in html
    assert 'rs123 &amp; rs456' in html
    assert 'Source name' in html
    assert 'javascript:bad()' not in html
    assert 'onclick=' not in html
    assert 'alert("not copied")' not in html
    assert 'id="source-table"' not in html
    assert 'class="result-tools"' not in html
    assert html.index('id="upstream-reports"') < html.index('<section id="loci">')


def test_locus_rows_are_readable_without_running_javascript(tmp_path):
    import json
    import re

    args, preflight, result, paths = _report_fixture(tmp_path)
    html = render_fine_mapping_html(build_fine_mapping_report(args, preflight, result, paths))
    section = html.split('<section id="loci">')[1].split('</section>')[0]
    body = re.search(r'<tbody>(.*?)</tbody>', section, re.S).group(1)
    payload = json.loads(re.search(r'id="locus-table-data">(.*?)</script>', section, re.S).group(1))
    assert body.count('<tr>') == len(payload) == 2
    assert 'chr1:90-130' in body
    assert 'chr1:500-700' in body
    assert 'severe_ld_z_mismatch' in body
    assert 'class="table-tools" hidden' in section


def test_unavailable_static_report_content_is_explicit(tmp_path):
    from postgwas.core.io.reports import read_static_report_content

    source = tmp_path / "nonstandard.html"
    source.write_text('<body><p>No main region</p></body>')
    assert read_static_report_content(source) == ""
