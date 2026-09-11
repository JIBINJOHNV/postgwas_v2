"""Presentation-only harmonisation HTML report tests."""

from pathlib import Path

from postgwas.modules.harmonisation.html_report import (
    render_dataset_report,
    render_run_report,
    write_dataset_report,
)


def _record(dataset_id: str = "study_one") -> dict:
    return {
        "dataset_id": dataset_id,
        "input_file": "/inputs/%s.tsv.gz" % dataset_id,
        "status": "OK",
        "input_variants": 100,
        "input_ready_variants": 95,
        "input_ready_snps": 90,
        "input_ready_indels_or_other_variants": 5,
        "genome_build": "GRCh38",
        "genome_build_source": "detector",
        "effect_type": "beta",
        "effect_type_source": "sample_sheet",
        "effect_scale": "BETA retained on the supplied additive scale",
        "p_value_type": "raw",
        "p_value_type_source": "sample_sheet",
        "frequency_type": "effect_allele_frequency",
        "frequency_type_source": "study_level_statistic",
        "strand_consensus": "forward",
        "strand_status": "complete",
        "strand_metadata_complete": True,
        "strand_reference_panel": "ALFA",
        "strand_reference_population": "EUR",
        "strand_reference_files": '["/reference/chr1.tsv.gz"]',
        "strand_chromosomes_expected": 1,
        "strand_chromosomes_completed": 1,
        "strand_variants_evaluated": 95,
        "strand_matched": 94,
        "strand_forward": 90,
        "strand_forward_swapped": 4,
        "strand_reverse_complement": 0,
        "strand_reverse_complement_swapped": 0,
        "strand_reference_unmatched_detected": 1,
        "strand_reference_unmatched_retained": 0,
        "strand_reference_unmatched": 1,
        "strand_palindromic_frequency_conflict": 0,
        "strand_palindromic_orientation_unavailable": 0,
        "strand_palindromic_frequency_discordant": 0,
        "strand_palindromic_ambiguous": 0,
        "strand_reference_ambiguous": 0,
        "strand_removed_total": 1,
        "strand_accounting_balanced": True,
        "harmonised_variants": 94,
        "rejected_variants": 6,
        "final_vcf_variants": 94,
        "final_vcf_snps": 89,
        "final_vcf_indels_or_other_variants": 5,
        "qc_active_rule_count": 6,
        "qc_passed_variants": 88,
        "qc_passed_snps": 84,
        "qc_failed_snps": 5,
        "qc_passed_snp_percent": 94.382,
        "qc_failed_snp_percent": 5.618,
        "closest_population": "EUR",
        "comparison_af_panel": "ALFA",
        "comparison_af_population": "EUR",
        "af_comparable_variants": 90,
        "af_concordant_variants": 89,
        "af_concordance_percent": 98.8889,
        "af_mismatched_variants": 1,
        "af_mismatch_percent": 1.1111,
        "study_af_missing": 0,
        "reference_af_missing": 4,
        "af_difference_cutoff": 0.2,
        "pre_vcf_invalid_effect_statistic_removals": 0,
        "qc_passed_missing_or_invalid_neff": 0,
        "qc_passed_missing_imputation_score": 0,
        "neff_reference_quantile": 0.9,
        "neff_reference_value": 100000,
        "neff_minimum_fraction_of_reference": 2 / 3,
        "neff_minimum_threshold": 66666.6667,
        "raw_low_neff_variants": 5,
        "raw_low_neff_percent": 5.3191,
        "qc_passed_low_neff_variants": 2,
        "qc_passed_low_neff_percent": 2.2727,
        "qc_passed_neff_upper_outliers": 0,
        "manifest": "/outputs/%s_manifest.json" % dataset_id,
        "screen_report": "/outputs/%s_screen.txt" % dataset_id,
        "html_report": "/outputs/%s_report.html" % dataset_id,
    }


def _manifest(dataset_id: str = "study_one") -> dict:
    return {
        "sample_id": dataset_id,
        "status": "OK",
        "sumstat_file": "/inputs/%s.tsv.gz" % dataset_id,
        "output_folder": "/outputs/%s" % dataset_id,
        "started": "2026-08-26 10:00:00",
        "completed_at": "2026-08-26T10:10:00-04:00",
        "elapsed_seconds": 600.0,
        "pre_vcf_summary": {
            "total_variant_infile": 100,
            "total_variant_read": 100,
            "total_variant_removed_missing_values": 2,
            "total_variant_removed_null_coords": 0,
            "total_variant_removed_unsupported_chromosomes": 0,
            "total_variant_removed_non_standard_alleles": 0,
            "total_variant_removed_duplicates": 3,
            "total_variant_remaining_for_harmonisation": 95,
            "total_variant_ready_snps": 90,
            "total_variant_ready_indels_or_other": 5,
            "total_variant_removed_chromosome_harmonisation": 1,
            "total_variant_rejected_end_to_end": 6,
            "total_variant_in_vcf_input": 94,
            "row_accounting_balanced": True,
            "row_accounting_complete": True,
        },
        "dataset": {
            "genome_build": {
                "inferred_build": "GRCh38",
                "matches": {"GRCh38": 80},
                "testable_variants": 81,
            },
            "study_decisions": {
                "effect_type": "beta",
                "pvalue_type": "raw",
                "frequency_type": "effect_allele_frequency",
                "strand_consensus": {"dominant_fraction": 0.99},
            },
            "chromosomes": {"1": {"status": "ok", "attempts": 1}},
            "chromosome_summaries": {
                "1": {
                    "status": "ok",
                    "attempt": 1,
                    "rows_in": 95,
                    "rows_out": 94,
                    "rejected": 1,
                    "balanced": True,
                    "elapsed": 10.5,
                    "warnings": 0,
                    "errors": 0,
                    "stage_qc": {
                        "eaf_qc": {
                            "strand_orientation": {
                                "actions": {"forward": 90, "forward_swapped": 4},
                                "reference_unmatched_detected": 1,
                                "reference_unmatched_retained": 0,
                                "reference_unmatched": 1,
                                "palindromic_ambiguous": 0,
                                "palindromic_frequency_conflict": 0,
                                "palindromic_orientation_unavailable": 0,
                                "palindromic_frequency_discordant": 0,
                                "reference_ambiguous": 0,
                                "reference_file": "/reference/chr1.tsv.gz",
                            }
                        }
                    },
                }
            },
            "rejected_variants_file": "/outputs/rejected.tsv.gz",
            "reject_reason_table": "/outputs/reject_reasons.tsv",
        },
        "qc_assessment": {
            "definition": "all active rules",
            "reports": {
                "json": "/outputs/vcf_qc.json",
                "summary": "/outputs/vcf_qc.tsv",
            },
        },
        "population_frequency_qc": {
            "status": "complete",
            "report": "/outputs/population_qc.json",
        },
        "merged_vcfs": {
            "grch37": "/outputs/study_GRCh37.vcf.gz",
            "grch38": "/outputs/study_GRCh38.vcf.gz",
            "merge_status": "OK",
        },
    }


def test_dataset_report_reuses_recorded_input_vcf_and_chromosome_evidence(tmp_path):
    destination = tmp_path / "study_report.html"
    written = write_dataset_report(destination, _record(), _manifest())

    text = written.read_text(encoding="utf-8")
    assert text.startswith("<!doctype html>")
    assert "Input and harmonisation accounting" in text
    assert "VCF content and virtual QC" in text
    assert "Chromosome-wise processing" in text
    assert "100" in text
    assert "94" in text
    assert "ALFA" in text
    assert "QC-passed low-Neff variants" in text
    assert "2.2727%" in text
    assert "QC-passed Neff upper-tail diagnostic" in text
    assert "Reference unmatched detected" in text
    assert "Reference unmatched retained by opt-in policy" in text
    assert "Unmatched retained" in text
    assert "Unmatched removed" in text
    assert "/reference/chr1.tsv.gz" in text
    assert "does not reread summary statistics or VCFs" in text


def test_run_report_contains_full_sections_for_every_selected_dataset():
    records = [_record("study_one"), _record("study_two")]
    manifests = {
        "study_one": _manifest("study_one"),
        "study_two": _manifest("study_two"),
    }

    text = render_run_report(records, manifests)

    assert "2 selected datasets" in text
    assert 'id="dataset-study_one"' in text
    assert 'id="dataset-study_two"' in text
    assert text.count("Input and harmonisation accounting") == 2
    assert text.count("VCF content and virtual QC") == 2
    assert text.count("Chromosome-wise processing") == 2


def test_report_escapes_user_controlled_text():
    record = _record()
    record["dataset_id"] = '<script>alert("x")</script>'
    record["failure_reason"] = "bad <input>"

    text = render_dataset_report(record, {})

    assert "<script>" not in text
    assert "&lt;script&gt;" in text
    assert "bad &lt;input&gt;" in text
