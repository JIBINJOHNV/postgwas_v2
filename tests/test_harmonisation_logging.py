"""Contract tests for the harmonisation audit log."""

import io
import inspect
from pathlib import Path

import pytest

from postgwas.core.pipeline_logging import (
    PipelineLogger,
    print_chromosome_summary,
)
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.service import (
    CHR_STEP_TOTAL,
    process_one_chromosome,
)


def _read(logger: PipelineLogger) -> str:
    logger.close()
    return Path(logger.log_path).read_text(encoding="utf-8")


def test_run_and_validation_logs_use_the_canonical_formatter_at_an_exact_path(tmp_path):
    path = tmp_path / "run_metadata" / "postgwas.log"
    logger = PipelineLogger(
        "postgwas", "run", str(path.parent),
        level="DEBUG", screen_level="ERROR", log_path=str(path),
    )
    logger.error("preflight failed")

    text = _read(logger)
    assert Path(logger.log_path) == path
    assert "] RUN    FAILED   preflight failed" in text
    assert " | ERROR | " not in text


def test_chromosome_worker_requires_prevalidated_split_files():
    source = inspect.getsource(process_one_chromosome)

    assert CHR_STEP_TOTAL == 16
    assert "fix_chr_pos_allele_column" not in source
    assert source.index("harmonise_effect_estimates(") < source.index(
        "harmonise_allele_frequency("
    )


def test_log_records_values_and_outcomes_without_repeating_documentation(tmp_path):
    logger = PipelineLogger(
        "study", "1", str(tmp_path), policies=default_policies(),
        level="INFO", screen_level="WARNING",
    )

    with logger.step(
        1, 2, "Validate p-values", "validate_pvalues", rows_in=10,
        policy_keys=["pvalue.clip_high"],
    ) as step:
        step.qc(
            "missing p-values", "This explanation must not be logged when zero rows are affected.",
            10, 10, reason="pvalue_missing",
        )
        step.qc(
            "out-of-range p-values", "Clip the affected values to the configured boundaries.",
            10, 10, changed=2, warn=True,
        )
        step.decide("p-value representation", "all sampled values are in [0, 1]", "raw")
        step.set_rows(10)

    text = _read(logger)
    assert "STEP     1/2 name='Validate p-values' function=validate_pvalues" in text
    assert "INPUT    variants rows=10" in text
    assert "PARAM    pvalue.clip_high=1.0" in text
    assert "PASS     missing p-values affected=0 rows=10" in text
    assert "This explanation must not be logged" not in text
    assert "WARNING  out-of-range p-values affected=2 removed=0 changed=2" in text
    assert "ACTION   Clip the affected values" in text
    assert "RESULT   out-of-range p-values rows_in=10 rows_out=10 removed=0 changed=2" in text
    assert "OBSERVED p-value representation evidence='all sampled values are in [0, 1]'" in text
    assert "DECIDE   p-value representation value=raw" in text
    assert "STATUS   1/2 status=OK rows_in=10 rows_out=10 removed=0" in text
    assert "Settings in use:" not in text
    assert "What it does:" not in text


def test_removed_rows_use_an_explicit_action_instead_of_a_reason_paragraph(tmp_path):
    logger = PipelineLogger("study", "3", str(tmp_path), policies=default_policies())
    with logger.step(1, 1, "Effect QC", "effect_qc", rows_in=10) as step:
        step.qc(
            "beta zero", "The effect size is exactly zero.", 10, 8,
            reason="beta_zero",
        )
        step.set_rows(8)

    text = _read(logger)
    assert "ACTION   remove_variants count=2 reason=beta_zero" in text
    assert "The effect size is exactly zero" not in text


def test_matched_but_retained_rows_are_not_reported_as_changed(tmp_path):
    logger = PipelineLogger("study", "4", str(tmp_path), policies=default_policies())
    with logger.step(1, 1, "Effect QC", "effect_qc", rows_in=10) as step:
        step.qc(
            "beta zero", "Zero is a valid null effect.", 10, 10,
            reason="beta_zero", matched=3, outcome="kept",
        )
        step.set_rows(10)

    text = _read(logger)
    assert "OBSERVED beta zero affected=3 removed=0 changed=0 matched=3" in text
    assert "ACTION   keep_variants count=3 reason=beta_zero" in text
    assert "RESULT   beta zero rows_in=10 rows_out=10 removed=0 changed=0 matched=3 outcome=kept" in text


def test_chromosome_screen_block_shows_concise_outcomes_without_duplicate_qc():
    stream = io.StringIO()
    summary = {
        "steps": [
            {
                "number": 1, "title": "Load the chromosome file",
                "rows_in": None, "rows_out": 112557, "removed": 0,
                "elapsed": 0.2, "status": "ok",
            },
            {
                "number": 12, "title": "SNP identifier",
                "rows_in": 112557, "rows_out": 112557, "removed": 0,
                "elapsed": 0.4, "status": "ok",
                "extra": {
                    "source": "study_column:SNP",
                    "identifiers_rewritten": 0,
                    "identifiers_filled": 1,
                },
            },
            {
                "number": 16, "title": "Annotate and lift over",
                "rows_in": 112557, "rows_out": 112557, "removed": 0,
                "elapsed": 0.8, "status": "ok",
                "extra": {"target_vcf": "/results/study_chr8_GRCh38_CSQ.vcf.gz"},
            },
        ],
        "qc_actions": [
            {
                "step": "16/16 liftover", "check": "LIFTED",
                "before": 112557, "after": 112251, "removed": 306,
                "changed": 0, "matched": 0, "outcome": "removed",
            },
        ],
        "stage_qc": {
            "eaf_qc": {
                "decision_source": "internal_eaf_study_decision",
                "final_eaf_col": "FRQ_U_186843",
                "study_decision_eaf_is_maf": False,
                "internal_af_initial_non_null_rows": 112557,
                "internal_af_outside_unit_interval_count": 0,
                "variants_removed": 0,
                "strand_orientation": {
                    "status": "success",
                    "reference_file": "/references/GRCh37_1000G_freq_chr8.tsv.gz",
                    "reference_population_column": "EUR",
                    "study_strand_consensus": "forward",
                    "actions": {
                        "forward": 10,
                        "forward_swapped": 110000,
                        "reverse_complement": 20,
                        "reverse_complement_swapped": 30,
                    },
                    "reference_unmatched": 111,
                    "palindromic_ambiguous": 2256,
                    "reference_ambiguous": 1,
                },
            },
            "sample_size_qc": {
                "ncase_source": "column:Nca",
                "ncontrol_source": "column:Nco",
                "Neff_status": "calculated_from_case_control",
                "min_value": 0,
                "removed_by_reason": {},
            },
            "effect_from_z_qc": {
                "has_beta": True,
                "has_se": True,
                "effect_column": "OR",
                "se_column": "SE",
                "variants_removed_total": 0,
            },
            "beta_or_oddsratio_qc": {
                "initial_variants": 112557,
                "effect_type": "odds_ratio",
                "effect_decision_source": "study-level decision",
                "conversion": "OR_to_Beta_log_transform_applied",
                "or_non_positive_count": 0,
                "or_non_positive_action": "null",
                "se_input_scale": "log_odds",
                "se_rescaled_by_or": False,
                "final_total": 112557,
            },
            "pval_detection_and_conversion_qc": {
                "pvalue_column": "P",
                "detected_scale": "raw",
                "pvalue_decision_source": "study-level decision",
                "initial_total_variants_with_pvalues": 112557,
                "after_filter_variants_with_pvalues": 112557,
                "variants_with_pvalues_clipped_low": 0,
                "variants_with_pvalues_clipped_high": 0,
                "pvalue_out_of_range_action": "clip",
                "pvalue_clip_low": 1e-300,
                "pvalue_clip_high": 0.999999,
            },
            "se_from_beta_pvalue_qc": {
                "initial_variants": 112557,
                "status": "SE already present",
                "total_variants": 112557,
                "se_undefined": 0,
            },
            "z_from_beta_se_qc": {
                "z_source": "calculated_from_beta_se",
                "variants_with_zero_beta": 261,
                "beta_zero_action": "keep",
                "variants_removed_invalid_beta_se": 0,
            },
            "effect_statistics_validation_qc": {
                "removed_total": 0,
                "removed_by_reason": {},
                "z_pval_concordance": {"action": "off"},
            },
            "info_qc": {
                "source": "study_column:INFO",
                "info_column": "INFO",
                "info_clip_min": 0,
                "info_clip_max": 1,
                "info_clip_tolerance": 1.05,
                "info_out_of_range_action": "clip",
                "info_on_missing_action": "keep",
                "rescaled_within_tolerance": 8,
                "out_of_range": 0,
                "missing_info_after": 0,
                "rejected_out_of_range": 0,
                "rejected_missing": 0,
            },
            "final_completeness_qc": {
                "required_fields": ["chr", "pos", "eaf", "beta", "se", "zscore", "pval"],
                "on_missing": "reject",
                "missing_by_field": {},
                "removed_total": 0,
            },
        },
        "rows_in": 112557,
        "rows_out": 112557,
        "rejected": 0,
        "reject_counts": {},
        "warnings": 0,
        "errors": 0,
        "log_path": "/results/logs/study_chr8.log",
    }

    text = print_chromosome_summary(
        "8", "study", 1, "OK", 168.0, "", summary,
        stream=stream, width=128,
    )
    flat = " ".join(text.split())

    assert "    🧬  Chromosome 8" in text
    assert "      ✅  Status               : completed in 2 m 48 s (attempt 1)" in text
    assert "🧮  Variants read" in text
    assert ": 112,557" in text
    assert "Allele frequency" in text
    assert "Strand orientation" in text
    assert "reference GRCh37_1000G_freq_chr8.tsv.gz, AF column EUR" in flat
    assert "study consensus forward; 110,060 matched" in flat
    assert "forward 10, forward-swapped 110,000" in flat
    assert "reverse-complement 20, reverse-complement-swapped 30" in flat
    assert "111 unmatched, 2,256 palindromic ambiguous and 1 reference ambiguous removed" in flat
    assert "internal column FRQ_U_186843 interpreted as EAF" in flat
    assert "110,030 allele-swapped frequencies inverted" in flat
    assert "cases from column Nca; controls from column Nco" in flat
    assert "effect from column OR and SE from column SE" in flat
    assert "positive OR converted with BETA = log(OR)" in flat
    assert "column P; raw scale" in flat
    assert "accepted range 1e-300 to 0.999999" in flat
    assert "supplied SE used directly" in flat
    assert "calculated as BETA / SE; 261 exact zero effects kept" in flat
    assert "study column INFO" in flat
    assert "8 corrected" in flat
    assert "study column SNP" in flat
    assert "1 missing IDs filled" in flat
    assert "7 required fields checked" in flat
    assert "🔹  Final result" in text
    assert "112,557 harmonised; none removed." in text
    assert "🔻  VCF and liftover" in text
    assert "112,251 of 112,557 retained after liftover (306 not lifted)" in flat
    assert "🔹  Messages" in text
    assert "no warnings or errors." in text
    assert "study_chr8_GRCh38_CSQ.vcf.gz" in text
    assert "study_chr8.log (detailed log)" in text
    assert "/results/" not in text
    assert "Load the chromosome file" not in text
    assert "derive_z_score_from_effect_and_standard_error" not in text
    assert "harmonize_sample_sizes" not in text
    assert "│" not in text
    assert stream.getvalue().startswith("\n    🧬  Chromosome 8")
    assert stream.getvalue().endswith("\n\n")


def test_chromosome_screen_reports_external_frequency_orientation_counts():
    summary = {
        "stage_qc": {
            "eaf_qc": {
                "decision_source": "external_eaf",
                "external_eaf_file": "reference_chr3.tsv.gz",
                "external_eaf_column": "EUR",
                "external_merge_direct_match_rows": 90,
                "external_merge_flip_match_rows": 7,
                "external_merge_unmatched_rows": 3,
                "variants_removed": 3,
            },
        },
        "warnings": 0,
        "errors": 0,
        "log_path": "/results/logs/study_chr3.log",
    }

    text = print_chromosome_summary(
        "3", "study", 1, "OK", 1.0, "", summary,
        stream=io.StringIO(), width=128,
    )
    flat = " ".join(text.split())

    assert "external file reference_chr3.tsv.gz, column EUR" in flat
    assert "90 direct matches, 7 allele-swapped matches, 3 unmatched; 3 removed" in flat


def test_chromosome_screen_labels_z_only_effect_as_standardized():
    summary = {
        "stage_qc": {
            "effect_from_z_qc": {
                "has_beta": False,
                "has_se": False,
                "has_zscore": True,
                "z_column": "Z",
                "beta_computed": True,
                "se_computed": True,
                "effect_estimate_scale": "standardized_effect_estimate",
                "effect_estimate_citation": "Zhu et al. 2016; PMID 27019110",
                "variants_removed_total": 0,
            },
        },
        "warnings": 2,
        "errors": 0,
        "log_path": "/results/logs/study_chr3.log",
    }

    text = print_chromosome_summary(
        "3", "study", 1, "OK", 1.0, "", summary,
        stream=io.StringIO(), width=128,
    )
    flat = " ".join(text.split())

    assert "standardized BETA estimate and SE derived from Z, EAF and NEFF" in flat
    assert "Zhu et al. 2016; PMID 27019110" in flat


def test_chromosome_screen_reports_absent_strand_actions_as_zero():
    summary = {
        "stage_qc": {
            "eaf_qc": {
                "strand_orientation": {
                    "status": "success",
                    "reference_file": "reference_chr5.tsv.gz",
                    "reference_population_column": "EUR",
                    "study_strand_consensus": "forward",
                    "actions": {"forward_swapped": 24227},
                    "reference_unmatched": 478,
                    "palindromic_ambiguous": 546,
                    "reference_ambiguous": 0,
                },
            },
        },
        "warnings": 0,
        "errors": 0,
    }

    text = print_chromosome_summary(
        "5", "study", 1, "OK", 1.0, "", summary,
        stream=io.StringIO(), width=128,
    )
    flat = " ".join(text.split())

    assert "24,227 matched" in flat
    assert "forward 0, forward-swapped 24,227" in flat
    assert "reverse-complement 0, reverse-complement-swapped 0" in flat
    assert "unknown" not in flat


def test_failed_chromosome_screen_never_claims_there_were_no_errors():
    text = print_chromosome_summary(
        "5", "study", 2, "FAILED", 4.0,
        "Reference frequency file was not found.",
        {"log_path": "/results/logs/study_chr5.log"},
        stream=io.StringIO(), width=128,
    )

    assert "🧬  Chromosome 5" in text
    assert "❌  Status" in text
    assert "failed in 4.0 s (attempt 2)" in text
    assert "❌  Messages" in text
    assert "1 error; see the log" in text
    assert "Reference frequency file was not found." in text
    assert "no warnings or errors" not in text


def test_failure_is_actionable_on_screen_and_complete_in_file(tmp_path):
    logger = PipelineLogger(
        "study", "2", str(tmp_path), policies=default_policies(),
        level="INFO", screen_level="WARNING",
    )

    with pytest.raises(FileNotFoundError):
        with logger.step(1, 1, "Read input", "read_input"):
            raise FileNotFoundError("missing.tsv")

    screen = logger.screen_text()
    text = _read(logger)
    assert "FAILED" in screen
    assert "ACTION" in screen
    assert "Traceback" not in screen
    assert "status=FAILED" in text
    assert "FileNotFoundError: missing.tsv" in text
    assert "Traceback" in text
