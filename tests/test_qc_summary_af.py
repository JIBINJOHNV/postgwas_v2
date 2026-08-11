"""Scientific contract tests for the single-pass merged-VCF QC assessment."""

import json
from pathlib import Path
import shutil
from unittest.mock import patch

import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.policies import PolicyError, default_policies
from postgwas.modules.harmonisation.qc_reporting import (
    harmonisation_qc_summary_lines,
    harmonisation_qc_takeaway_lines,
    summarise_strand_orientation,
)
from postgwas.modules.qc_summary.assessment import (
    assess_variant_table,
    extract_vcf_assessment_table,
    run_vcf_qc_assessment,
)


def _assessment_table(path: Path) -> Path:
    path.write_text(
        "CHROM\tPOS\tREF\tALT\tINFO_AF\tINFO_EXTERNAL_AF\t"
        "FORMAT_AF\tFORMAT_SI\tFORMAT_LP\tFORMAT_NEF\n"
        "1\t100\tA\tG\t0.2\t0.21\t0.2\t0.9\t8\t100\n"
        "1\t101\tC\tT\t0.3\t.\t0.3\t0.9\t8\t100\n"
        "1\t102\tA\tC\t0.001\t0.001\t0.001\t0.9\t8\t100\n"
        "1\t103\tG\tA\t0.2\t0.2\t0.2\t0.5\t8\t100\n"
        "1\t104\tA\tAT\t0.2\t0.2\t0.2\t0.9\t8\t100\n"
        "1\t105\tA\tT\t0.5\t0.5\t0.5\t0.9\t8\t100\n"
        "6\t30000000\tA\tG\t0.2\t0.2\t0.2\t0.9\t8\t100\n"
        "1\t106\tA\tAT\t0.001\t0.5\t0.001\t0.5\t8\t1000\n",
        encoding="utf-8",
    )
    return path


def _harmonisation_config():
    module = load_configuration().modules.harmonisation
    return module.vcf_processing.model_dump(), dict(module.output_layout.root)


def test_assessment_applies_every_rule_to_raw_data_and_combines_failures(tmp_path):
    vcf_config, _ = _harmonisation_config()
    assessment = assess_variant_table(
        _assessment_table(tmp_path / "raw_vcf.tsv"),
        default_policies(),
        "EUR",
        vcf_config["qc_fields"],
        delimiter=vcf_config["table_delimiter"],
        null_values=vcf_config["table_null_values"],
    )

    assert assessment["raw"]["num_records"] == 8
    assert assessment["raw"]["num_snps"] == 6
    assert assessment["raw"]["num_non_snps"] == 2
    assert assessment["raw"]["ts_tv_ratio"] == 2.0
    assert assessment["raw"]["external_af_missing"] == 1
    assert assessment["raw"]["af_comparable"] == 7
    assert assessment["raw"]["af_difference_above_cutoff"] == 1
    assert assessment["raw"]["effective_sample_size_available"] == 8
    assert assessment["raw"]["effective_sample_size_mean"] == pytest.approx(212.5)
    assert assessment["raw"]["effective_sample_size_above_outlier_threshold"] == 0

    rules = {rule["key"]: rule for rule in assessment["rules"]}
    assert rules["allele_frequency"]["failed_raw"] == 2
    assert rules["imputation_quality"]["failed_raw"] == 2
    assert rules["external_af_concordance"]["failed_raw"] == 2
    assert rules["variant_type"]["failed_raw"] == 2
    assert rules["palindromic_variants"]["failed_raw"] == 1
    assert rules["mhc_region"]["failed_raw"] == 1

    assert assessment["qc_passed"]["num_records"] == 1
    assert assessment["qc_passed"]["effective_sample_size_available"] == 1
    assert assessment["qc_passed"]["effective_sample_size_mean"] == 100.0
    assert assessment["qc_passed"]["effective_sample_size_standard_deviation"] is None
    assert assessment["excluded_total"] == 7
    assert assessment["rule_match_total"] == 10
    assert assessment["active_rule_count"] == 6
    assert assessment["accounting_balanced"] is True
    assert assessment["overlap_variants"] == 1
    assert assessment["extra_rule_matches"] == 3
    assert "no filtered VCF is created" in assessment["definition"]


def test_sample_size_outliers_are_reported_before_and_after_combined_qc(tmp_path):
    vcf_config, _ = _harmonisation_config()
    policies = default_policies().with_overrides({
        "qc.sample_size_outlier_standard_deviations": 1.5,
    })
    assessment = assess_variant_table(
        _assessment_table(tmp_path / "raw_vcf.tsv"),
        policies,
        "EUR",
        vcf_config["qc_fields"],
        delimiter=vcf_config["table_delimiter"],
        null_values=vcf_config["table_null_values"],
    )

    assert assessment["sample_size_outlier_standard_deviations"] == 1.5
    assert assessment["raw"]["effective_sample_size_outlier_threshold"] == pytest.approx(
        212.5 + 1.5 * 318.1980515339464
    )
    assert assessment["raw"]["effective_sample_size_above_outlier_threshold"] == 1
    assert assessment["qc_passed"]["effective_sample_size_above_outlier_threshold"] == 0


def test_unusable_sample_sizes_are_reported_but_not_used_in_distribution(tmp_path):
    table = _assessment_table(tmp_path / "raw_vcf.tsv")
    contents = table.read_text(encoding="utf-8")
    contents = contents.replace("\t100\n", "\t.\n", 1).replace("\t100\n", "\t0\n", 1)
    table.write_text(contents, encoding="utf-8")
    vcf_config, _ = _harmonisation_config()

    assessment = assess_variant_table(
        table,
        default_policies(),
        "EUR",
        vcf_config["qc_fields"],
        delimiter=vcf_config["table_delimiter"],
        null_values=vcf_config["table_null_values"],
    )

    assert assessment["raw"]["effective_sample_size_available"] == 6
    assert assessment["raw"]["effective_sample_size_missing_or_invalid"] == 2


def test_bcftools_extracts_all_qc_fields_in_one_data_query(tmp_path):
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("bcftools is not installed")
    vcf = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    table = tmp_path / "assessment.tsv"
    vcf_config, _ = _harmonisation_config()

    extract_vcf_assessment_table(
        vcf_path=vcf,
        table_path=table,
        dataset_id="study",
        external_af_name="EUR",
        vcf_fields=vcf_config["qc_fields"],
        table_delimiter=vcf_config["table_delimiter"],
        io_buffer_bytes=vcf_config["io_buffer_bytes"],
        bcftools_bin=bcftools,
    )

    rows = table.read_text(encoding="utf-8").splitlines()
    assert len(rows) == 3
    assert rows[0].split("\t") == [
        "CHROM", "POS", "REF", "ALT", "INFO_AF",
        "INFO_EXTERNAL_AF", "FORMAT_AF", "FORMAT_SI", "FORMAT_LP", "FORMAT_NEF",
    ]
    assert rows[1].split("\t") == [
        "1", "100", "A", "G", "0.2", "0.21", "0.2", "0.9", "8", "1000",
    ]
    assert rows[2].split("\t")[5] == "."


def test_qc_vcf_query_fields_are_configuration_driven(tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"placeholder")
    table = tmp_path / "assessment.tsv"
    vcf_config, _ = _harmonisation_config()
    fields = dict(vcf_config["qc_fields"])
    fields.update({
        "study_info_af": "%INFO/STUDY_AF",
        "external_info_af": "%INFO/REF_{external_af}",
        "study_format_af": "[%EAF]",
        "imputation_format": "[%INFO_SCORE]",
        "log_pvalue_format": "[%MLOG10P]",
        "effective_sample_size_format": "[%EFFECTIVE_N]",
    })
    with patch(
        "postgwas.modules.qc_summary.assessment.extract_vcf_table",
        return_value=str(table),
    ) as extract:
        extract_vcf_assessment_table(
            vcf_path=vcf,
            table_path=table,
            dataset_id="study",
            external_af_name="EUR",
            vcf_fields=fields,
            table_delimiter=vcf_config["table_delimiter"],
            io_buffer_bytes=vcf_config["io_buffer_bytes"],
            bcftools_bin="configured-bcftools",
        )
    columns = extract.call_args.args[3]
    assert columns["INFO_AF"] == "%INFO/STUDY_AF"
    assert columns["INFO_EXTERNAL_AF"] == "%INFO/REF_EUR"
    assert columns["FORMAT_AF"] == "[%EAF]"
    assert columns["FORMAT_NEF"] == "[%EFFECTIVE_N]"


def test_sample_size_outlier_threshold_is_validated_and_defaults_to_five():
    assert default_policies().get(
        "qc.sample_size_outlier_standard_deviations"
    ) == 5.0
    with pytest.raises(PolicyError, match="qc.sample_size_outlier_standard_deviations"):
        default_policies().with_overrides({
            "qc.sample_size_outlier_standard_deviations": 0,
        })


def test_configuration_rejects_an_invalid_info_interval_before_analysis():
    with pytest.raises(PolicyError) as error:
        default_policies().with_overrides({
            "filter": {"info_cutoff": 0.9, "info_max": 0.8},
        })
    assert "filter.info_cutoff" in str(error.value)
    assert "filter.info_max" in str(error.value)


def test_run_writes_reports_but_never_writes_a_filtered_vcf(tmp_path):
    raw_vcf = tmp_path / "study_GRCh37_merged.vcf.gz"
    raw_vcf.write_bytes(b"placeholder")
    vcf_config, output_layout = _harmonisation_config()
    output_layout.update({
        "qc_assessment_summary": "reports/{dataset_id}/{build}/metrics.tsv",
        "qc_filter_rules": "reports/{dataset_id}/{build}/rules.tsv",
        "qc_assessment_json": "reports/{dataset_id}/{build}/assessment.json",
        "qc_assessment_temporary": "reports/{dataset_id}/{build}/.working_",
    })

    def fake_extract(**kwargs):
        _assessment_table(Path(kwargs["table_path"]))
        return str(kwargs["table_path"])

    with patch(
        "postgwas.modules.qc_summary.assessment.extract_vcf_assessment_table",
        side_effect=fake_extract,
    ):
        result = run_vcf_qc_assessment(
            vcf_path=raw_vcf,
            output_directory=tmp_path,
            dataset_id="study",
            genome_build="GRCh37",
            external_af_name="EUR",
            vcf_fields=vcf_config["qc_fields"],
            output_layout=output_layout,
            table_delimiter=vcf_config["table_delimiter"],
            table_null_values=vcf_config["table_null_values"],
            table_null_output=vcf_config["table_null_output"],
            temporary_table_suffix=vcf_config["temporary_table_suffix"],
            io_buffer_bytes=vcf_config["io_buffer_bytes"],
            policies=default_policies(),
            bcftools_bin="configured-bcftools",
        )

    assert Path(result["reports"]["summary"]).is_file()
    assert Path(result["reports"]["rules"]).is_file()
    assert Path(result["reports"]["json"]).is_file()
    assert Path(result["reports"]["summary"]) == (
        tmp_path / "reports" / "study" / "GRCh37" / "metrics.tsv"
    )
    assert not list(tmp_path.rglob("*filtered*.vcf*"))
    assert not list(tmp_path.rglob(".working_*.tsv"))
    saved = json.loads(Path(result["reports"]["json"]).read_text(encoding="utf-8"))
    assert saved["qc_passed"]["num_records"] == 1
    summary = Path(result["reports"]["summary"]).read_text(encoding="utf-8")
    assert "raw\teffective_sample_size_mean\t212.5" in summary
    assert "qc_passed\teffective_sample_size_mean\t100.0" in summary


def test_screen_report_has_raw_rule_and_combined_final_sections(tmp_path):
    vcf_config, _ = _harmonisation_config()
    assessment = assess_variant_table(
        _assessment_table(tmp_path / "raw.tsv"), default_policies(), "EUR",
        vcf_config["qc_fields"],
        delimiter=vcf_config["table_delimiter"],
        null_values=vcf_config["table_null_values"],
    )
    assessment["raw_vcf"] = str(tmp_path / "study_GRCh37_merged.vcf.gz")
    assessment["reports"] = {
        "summary": str(tmp_path / "study_qc_assessment.tsv"),
        "rules": str(tmp_path / "study_qc_filter_rules.tsv"),
    }
    lines = harmonisation_qc_summary_lines(
        {
            "total_variant_infile": 10,
            "total_variant_read": 10,
            "total_variant_remaining_for_harmonisation": 10,
            "total_variant_removed_palindromic_ambiguous": 1,
            "total_variant_removed_reference_unmatched": 1,
            "total_variant_removed_reference_ambiguous": 0,
            "total_variant_removed_chromosome_harmonisation": 2,
            "total_variant_passed_chromosome_harmonisation": 8,
            "total_variant_in_vcf_input": 8,
        },
        assessment,
    )
    output = "\n".join(lines)
    compact_output = " ".join(output.split())

    assert "1. Before VCF creation" in output
    assert "2. Raw merged VCF" in output
    assert "Effective sample-size distribution" in output
    assert "Above upper outlier threshold" in output
    assert "3. QC conditions assessed on the raw VCF" in output
    assert "Rule 1 · Study allele frequency" in output
    assert "4. Final QC-passed assessment" in output
    assert "FORMAT/NEF values above threshold" in output
    assert "no QC-filtered VCF is created" in compact_output
    assert "tested against every active condition" in compact_output
    assert "all 6 active conditions are applied together" in compact_output
    assert "10 total = 7 unique failed variants + 3 additional overlapping matches" in compact_output
    assert "8 raw = 1 QC-passed + 7 failed combined QC" in compact_output
    assert "Passed initial input validation" in output
    assert "Palindromic ambiguous removed" in output
    assert "Reference unmatched removed" in output
    assert "Reference ambiguous removed" in output
    assert "All chromosome-stage removals" in output
    assert "Passed chromosome harmonisation" in output
    assert "10 initial = 8 sent to GWAS-to-VCF + 2 removed" in compact_output

    # Every field at the same hierarchy depth uses one colon column. Remove
    # the warning symbol's variation selector before measuring terminal text.
    aligned_fields = []
    for indentation in (6, 8, 12):
        fields = [
            line.replace("\ufe0f", "")
            for line in lines
            if len(line) - len(line.lstrip(" ")) == indentation and " : " in line
        ]
        assert fields
        assert len({line.index(" : ") for line in fields}) == 1
        aligned_fields.extend(fields)
    assert len({line.index(" : ") for line in aligned_fields}) == 1

    assert next(line for line in lines if "Frequency and imputation fields" in line).startswith(
        " " * 10
    )
    assert next(line for line in lines if "Rule 1 ·" in line).startswith(" " * 10)
    rule_detail = next(
        line for line in lines
        if "Missing FORMAT/AF" in line and "raw VCF" in line
    )
    assert len(rule_detail) - len(rule_detail.lstrip(" ")) == 12
    assert next(line for line in lines if "Final virtual subset metrics" in line).startswith(
        " " * 10
    )


def test_final_takeaways_reuse_existing_counts_and_show_precise_af_percentages(tmp_path):
    vcf_config, _ = _harmonisation_config()
    assessment = assess_variant_table(
        _assessment_table(tmp_path / "raw.tsv"), default_policies(), "EUR",
        vcf_config["qc_fields"],
        delimiter=vcf_config["table_delimiter"],
        null_values=vcf_config["table_null_values"],
    )
    lines = harmonisation_qc_takeaway_lines(
        {
            "total_variant_infile": 10,
            "total_variant_read": 10,
            "total_variant_in_vcf_input": 8,
            "total_variant_with_invalid_beta_se": 0,
        },
        assessment,
        genome_build_info={
            "inferred_build": "GRCh37",
            "matches": {"GRCh37": 99},
            "testable_variants": 100,
        },
        study_decisions={
            "effect_type": "odds_ratio",
            "effect_type_detected": "odds_ratio",
            "effect_type_source": "detector",
            "pvalue_type": "raw",
            "pvalue_type_detected": "raw",
            "pvalue_type_source": "detector",
            "frequency_type": "effect_allele_frequency",
            "strand": "forward",
            "strand_consensus": {"dominant_fraction": 1.0},
        },
        population_frequency={
            "closest_population": "EUR",
            "selected_population_check": {
                "normalized_population": "EUR",
                "status": "match",
            },
            "external_file_checks": [],
        },
        dataset_summary={
            "completed": ["1"],
            "failed": [],
            "chromosomes": {"1": {"status": "ok"}},
            "chromosome_summaries": {
                "1": {
                    "stage_qc": {
                        "eaf_qc": {
                            "strand_orientation": {
                                "status": "success",
                                "initial_variants": 10,
                                "final_variants": 8,
                                "reference_file": "/reference/panel_chr1.tsv.gz",
                                "reference_population_column": "EUR",
                                "study_strand_consensus": "forward",
                                "actions": {
                                    "forward": 7,
                                    "forward_swapped": 1,
                                },
                                "reference_unmatched": 1,
                                "palindromic_ambiguous": 1,
                                "reference_ambiguous": 0,
                            }
                        }
                    }
                }
            },
            "reconciliation": {"balanced": True},
            "rejected_variants_rows": 2,
            "rejected_variants_file": str(tmp_path / "rejected.tsv.gz"),
        },
        merge_summary={
            "merge_status": "OK",
            "required_merge_failures": [],
        },
        primary_outputs={
            "GRCh37": "input.vcf.gz",
            "GRCh38": "target.vcf.gz",
            "gwas2vcf": "raw.vcf.gz",
        },
        final_status="OK",
        reference_source="1000G",
        reference_population="EUR",
    )
    output = " ".join("\n".join(lines).split())

    assert "Final QC takeaways" in output
    assert "10 / 10 variants" in output
    assert "8 variants used → 8 variants in final VCF" in output
    assert "6 SNPs · 2 indels / other variants" in output
    assert "1 / 6 SNPs passed all 6 active rules (16.67%)" in output
    assert "5 SNPs failed ≥1 active rule (83.33%)" in output
    assert "6 / 7 concordant (85.7143%)" in output
    assert "1 mismatched (14.2857%)" in output
    assert (
        "0 missing study AF · 1 missing reference AF · "
        "|AF difference| ≤ 0.2"
    ) in output
    assert "2 variants removed with reasons recorded" in output
    assert "GRCh37 · automatically detected" in output
    assert "1000G · EUR · 1 chromosome reference file" in output
    assert "8 total · forward 7 · forward-swapped 1" in output
    assert "10 evaluated = 8 retained + 2 removed" in output
    assert "odds_ratio · automatically detected" in output
    assert "raw · automatically detected" in output
    assert "7 merged-VCF variants failed ≥1 of 6 active rules" in output
    assert "Variant and chromosome counts reconciled" in output
    assert "1 / 1 completed" in output
    assert "Final merged VCF" in output
    assert "Indexed and validated" in output
    assert "3 primary outputs completed · dataset status OK" in output

    card_headings = [
        "Variant flow",
        "Reference alignment",
        "Allele frequency",
        "Statistical quality",
        "Audit trail",
        "VCF integrity",
    ]
    for heading in card_headings:
        line = next(line for line in lines if heading in line)
        assert len(line) - len(line.lstrip(" ")) == 8
    assert lines.count("") == len(card_headings)
    detail_lines = [
        line.replace("\ufe0f", "")
        for line in lines
        if len(line) - len(line.lstrip(" ")) == 12 and " : " in line
    ]
    assert detail_lines
    assert len({line.index(" : ") for line in detail_lines}) == 1
    assert next(line for line in lines if "Input / read" in line).startswith(" " * 12)


def test_strand_summary_combines_completed_chromosomes_without_rescanning_rows():
    def chromosome(initial, final, actions, unmatched, palindromic, ambiguous, name):
        return {
            "stage_qc": {
                "eaf_qc": {
                    "strand_orientation": {
                        "status": "success",
                        "initial_variants": initial,
                        "final_variants": final,
                        "reference_file": "/reference/%s.tsv.gz" % name,
                        "reference_population_column": "EUR",
                        "study_strand_consensus": "forward",
                        "actions": actions,
                        "reference_unmatched": unmatched,
                        "palindromic_ambiguous": palindromic,
                        "reference_ambiguous": ambiguous,
                    }
                }
            }
        }

    summary = summarise_strand_orientation({
        "completed": ["1", "2"],
        "failed": [],
        "chromosomes": {"1": {"status": "ok"}, "2": {"status": "ok"}},
        "chromosome_summaries": {
            "1": chromosome(
                10, 8,
                {"forward": 5, "forward_swapped": 2, "reverse_complement": 1},
                1, 0, 1, "panel_chr1",
            ),
            "2": chromosome(
                20, 17,
                {
                    "forward": 15,
                    "forward_swapped": 1,
                    "reverse_complement_swapped": 1,
                },
                1, 1, 1, "panel_chr2",
            ),
        },
    })

    assert summary == {
        "status": "success",
        "chromosomes_expected": 2,
        "chromosomes_completed": 2,
        "chromosomes_summarized": 2,
        "complete_chromosome_coverage": True,
        "metadata_complete": True,
        "reference_population": "EUR",
        "reference_file_count": 2,
        "reference_files": [
            "/reference/panel_chr1.tsv.gz",
            "/reference/panel_chr2.tsv.gz",
        ],
        "study_consensus": "forward",
        "variants_evaluated": 30,
        "variants_retained": 25,
        "variants_matched": 25,
        "forward": 20,
        "forward_swapped": 3,
        "reverse_complement": 1,
        "reverse_complement_swapped": 1,
        "disabled": 0,
        "reference_unmatched": 2,
        "palindromic_ambiguous": 1,
        "reference_ambiguous": 2,
        "removed_total": 5,
        "accounting_balanced": True,
    }


def test_strand_summary_marks_partial_disabled_and_malformed_evidence():
    disabled = {
        "stage_qc": {
            "eaf_qc": {
                "strand_orientation": {
                    "status": "disabled",
                    "initial_variants": 4,
                    "final_variants": 4,
                    "actions": {"disabled": 4},
                }
            }
        }
    }
    result = summarise_strand_orientation({
        "completed": ["1"],
        "failed": [],
        "chromosomes": {"1": {"status": "ok"}},
        "chromosome_summaries": {"1": disabled},
    })
    assert result["status"] == "disabled"
    assert result["removed_total"] == 0
    assert result["accounting_balanced"] is True

    malformed = {
        "stage_qc": {
            "eaf_qc": {
                "strand_orientation": {
                    "status": "success",
                    "initial_variants": 4.5,
                    "final_variants": 4,
                    "actions": {"forward": 4},
                    "reference_unmatched": 0,
                    "palindromic_ambiguous": 0,
                    "reference_ambiguous": 0,
                }
            }
        }
    }
    result = summarise_strand_orientation({
        "completed": ["1"],
        "failed": ["2"],
        "chromosomes": {
            "1": {"status": "ok"},
            "2": {"status": "failed"},
        },
        "chromosome_summaries": {"1": malformed},
    })
    assert result["status"] == "partial"
    assert result["complete_chromosome_coverage"] is False
    assert result["variants_evaluated"] is None
    assert result["accounting_balanced"] is False
