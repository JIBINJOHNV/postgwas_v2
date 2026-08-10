"""Scientific and failure-path tests for input-versus-VCF validation."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import patch

import polars as pl

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.concordance.analysis import (
    compare_input_to_vcf,
    prepare_input_table,
)
from postgwas.modules.harmonisation.concordance.errors import (
    ConcordanceValidationError,
)
from postgwas.modules.harmonisation.concordance.cli import get_validation_parser
from postgwas.modules.harmonisation.concordance.service import (
    _manifest_context,
    _vcf_build,
    run_concordance_validation,
)
from postgwas.modules.harmonisation.policies import load_policies


def _row(**updates):
    values = {
        "dataset_id": "study",
        "input_file": Path("input.tsv"),
        "delimiter": "tab",
        "chromosome_column": "CHR",
        "position_column": "BP",
        "chromosome_position_column": None,
        "effect_allele_column": "EA",
        "other_allele_column": "OA",
        "effect_allele_frequency_column": "EAF",
        "effect_column": "BETA",
        "standard_error_column": "SE",
        "z_score_column": "Z",
        "p_value_column": "P",
        "effect_type": "beta",
        "p_value_type": "raw",
        "external_eaf_file": None,
        "external_eaf_column": None,
    }
    values.update(updates)
    return SimpleNamespace(**values)


def _settings():
    return load_configuration().modules.harmonisation.concordance_validation


def _policies():
    return load_policies(load_configuration().modules.harmonisation.policies)


class ConcordanceAnalysisTests(unittest.TestCase):
    def test_standalone_cli_has_the_requested_contract(self):
        parser = get_validation_parser()
        args = parser.parse_args([
            "--sample-sheet", "studies.csv",
            "--dataset-id", "study1",
            "--vcf", "study1_GRCh37_merged.vcf.gz",
            "--output-directory", "results",
        ])
        self.assertEqual(args.dataset_id, "study1")
        self.assertEqual(args.vcf, "study1_GRCh37_merged.vcf.gz")
        help_text = parser.format_help()
        self.assertIn("effect estimates, allele frequencies, and Z scores", help_text)
        self.assertIn("Indels are compared exactly as written", help_text)

    def test_top_level_validate_flag_dispatches_to_standalone_command(self):
        from postgwas.__main__ import main

        with patch(
            "postgwas.modules.harmonisation.concordance.cli.main", return_value=7,
        ) as validation_main, patch(
            "sys.argv", ["postgwas", "--validate", "--help"],
        ):
            self.assertEqual(main(), 7)
            validation_main.assert_called_once_with()
            self.assertEqual(sys.argv, ["postgwas --validate", "--help"])

    def test_vcf_contig_assembly_wins_over_liftover_command_history(self):
        header = (
            "##contig=<ID=1,length=248956422,assembly=GRCh38>\n"
            "##bcftools_liftoverCommand=input_GRCh37.vcf output_GRCh38.vcf\n"
        )
        self.assertEqual(_vcf_build(header, ["GRCh37", "GRCh38"]), "GRCh38")

    def test_vcf_build_detection_uses_configured_build_names(self):
        header = "##contig=<ID=1,length=10,assembly=CustomBuild>\n"
        self.assertEqual(
            _vcf_build(header, ["ReferenceA", "CustomBuild"]),
            "CustomBuild",
        )

    def test_manifest_uses_reference_confirmed_eaf_not_initial_maf_suspicion(self):
        chromosome_summary = {
            "status": "ok",
            "stage_qc": {"eaf_qc": {"maf_reference_decision": "eaf"}},
        }
        build, decisions = _manifest_context({
            "dataset": {
                "genome_build": {"inferred_build": "GRCh37"},
                "study_decisions": {
                    "eaf_is_maf": True,
                    "eaf_is_maf_source": "study_level_statistic",
                },
                "completed": ["1", "2"],
                "chromosome_summaries": {
                    "1": chromosome_summary,
                    "2": chromosome_summary,
                },
            }
        })

        self.assertEqual(build, "GRCh37")
        self.assertTrue(decisions["eaf_is_maf_initial"])
        self.assertFalse(decisions["eaf_is_maf"])
        self.assertEqual(decisions["frequency_type"], "effect_allele_frequency")
        self.assertEqual(decisions["eaf_reference_decisions"]["eaf"], 2)

    def test_manifest_leaves_mixed_frequency_evidence_unresolved(self):
        _, decisions = _manifest_context({
            "dataset": {
                "study_decisions": {"eaf_is_maf": True},
                "completed": ["1", "2"],
                "chromosome_summaries": {
                    "1": {
                        "status": "ok",
                        "stage_qc": {
                            "eaf_qc": {"maf_reference_decision": "eaf"}
                        },
                    },
                    "2": {
                        "status": "ok",
                        "stage_qc": {
                            "eaf_qc": {"maf_reference_decision": "inconclusive"}
                        },
                    },
                },
            }
        })

        self.assertIsNone(decisions["eaf_is_maf"])
        self.assertEqual(decisions["frequency_type"], "unresolved")

    def test_manifest_without_frequency_provenance_does_not_assume_eaf(self):
        _, decisions = _manifest_context({})

        self.assertIsNone(decisions["eaf_is_maf"])
        self.assertEqual(decisions["frequency_type"], "unresolved")
        self.assertEqual(
            decisions["eaf_is_maf_source"], "frequency_decision_missing"
        )

    def test_explicit_absent_internal_frequency_uses_external_eaf_contract(self):
        _, decisions = _manifest_context({
            "dataset": {
                "study_decisions": {
                    "eaf_is_maf": None,
                    "eaf_is_maf_source": "no_frequency_column",
                }
            }
        })

        self.assertFalse(decisions["eaf_is_maf"])
        self.assertEqual(decisions["frequency_type"], "effect_allele_frequency")
        self.assertEqual(
            decisions["eaf_is_maf_source"], "external_effect_allele_frequency"
        )

    def test_matches_all_supported_orientations_without_normalizing_indels(self):
        input_frame = pl.DataFrame({
            "CHR": [1, 1, 1, 1, 1, 1],
            "BP": [101, 102, 103, 104, 105, 106],
            "EA": ["A", "T", "A", "G", "A", "AT"],
            "OA": ["G", "C", "C", "T", "T", "A"],
            "EAF": [0.2, 0.7, 0.1, 0.4, 0.2, 0.3],
            "BETA": [0.2, -0.3, 0.4, 0.5, 0.6, 0.7],
            "SE": [0.1] * 6,
            "Z": [2.0, -3.0, 4.0, 5.0, 6.0, 7.0],
            "P": [0.05] * 6,
        })
        vcf_frame = pl.DataFrame({
            "CHROM": [1] * 6,
            "POS": [101, 102, 103, 104, 105, 106],
            "ID": ["v1", "v2", "v3", "v4", "v5", "v6"],
            "REF": ["G", "T", "G", "C", "T", "A"],
            "ALT": ["A", "C", "T", "A", "A", "AT"],
            "ES": [0.2, 0.3, 0.4, -0.5, 0.6, 0.7],
            "SE": [0.1] * 6,
            "EZ": [2.0, 3.0, 4.0, -5.0, 6.0, 7.0],
            "AF": [0.2, 0.3, 0.1, 0.6, 0.2, 0.3],
        })

        result = compare_input_to_vcf(
            input_frame,
            vcf_frame,
            _row(),
            effect_type="beta",
            p_value_type="raw",
            eaf_is_maf=False,
            settings=_settings(),
            policies=_policies(),
        )

        self.assertEqual(result.status, "PASS")
        self.assertEqual(result.summary["matched_variants"], 6)
        self.assertEqual(result.summary["variant_types"]["snps"]["exact_matched_variants"], 5)
        self.assertEqual(result.summary["variant_types"]["indels"]["exact_matched_variants"], 1)
        self.assertFalse(
            result.summary["variant_types"]["indels"][
                "representation_adjustment_applied"
            ]
        )
        self.assertEqual(result.summary["palindromic_variants_excluded_from_orientation"], 1)
        self.assertEqual(result.metric_summary["effect"]["concordant"], 5)
        self.assertEqual(result.metric_summary["allele_frequency"]["concordant"], 5)
        self.assertEqual(result.metric_summary["z_score"]["concordant"], 5)
        self.assertEqual(
            set(result.matched.get_column("match_type")),
            {
                "direct",
                "allele_swapped",
                "strand_complement",
                "strand_complement_swapped",
                "palindromic_as_listed",
            },
        )
        indel = result.matched.filter(pl.col("input_pos") == 106).row(0, named=True)
        self.assertEqual((indel["vcf_ref"], indel["vcf_alt"]), ("A", "AT"))

    def test_stratifies_snps_and_possible_normalized_indel_representations(self):
        input_frame = pl.DataFrame({
            "CHR": [1, 1], "BP": [101, 200], "EA": ["A", "AA"],
            "OA": ["G", "A"], "EAF": [0.2, 0.3], "BETA": [0.2, 0.4],
            "SE": [0.1, 0.1], "Z": [2.0, 4.0], "P": [0.05, 0.01],
        })
        vcf_frame = pl.DataFrame({
            "CHROM": [1, 1], "POS": [101, 199], "ID": ["snp", "indel"],
            "REF": ["G", "AA"], "ALT": ["A", "A"],
            "ES": [0.2, 0.4], "SE": [0.1, 0.1], "EZ": [2.0, 4.0],
            "AF": [0.2, 0.3],
        })

        result = compare_input_to_vcf(
            input_frame, vcf_frame, _row(), effect_type="beta",
            p_value_type="raw", eaf_is_maf=False, settings=_settings(),
            policies=_policies(),
        )

        snps = result.summary["variant_types"]["snps"]
        indels = result.summary["variant_types"]["indels"]
        self.assertEqual(snps["status"], "PASS")
        self.assertEqual(snps["exact_matched_variants"], 1)
        self.assertEqual(indels["status"], "WARNING")
        self.assertEqual(indels["exact_matched_variants"], 0)
        self.assertEqual(indels["maximum_possible_representation_pairs"], 1)
        self.assertTrue(indels["representation_adjustment_applied"])
        self.assertEqual(result.status, "WARNING")
        self.assertEqual(result.input_only.filter(pl.col("input_is_snp").not_()).height, 1)
        self.assertEqual(result.vcf_only.filter(pl.col("vcf_is_snp").not_()).height, 1)

        strict = _settings().model_copy(
            update={"indel_representation_action": "fail"}
        )
        strict_result = compare_input_to_vcf(
            input_frame, vcf_frame, _row(), effect_type="beta",
            p_value_type="raw", eaf_is_maf=False, settings=strict,
            policies=_policies(),
        )
        self.assertEqual(strict_result.status, "FAIL")
        self.assertFalse(
            strict_result.summary["variant_types"]["indels"][
                "representation_adjustment_applied"
            ]
        )

    def test_does_not_assume_indel_normalization_when_snp_concordance_fails(self):
        input_frame = pl.DataFrame({
            "CHR": [1, 1], "BP": [101, 200], "EA": ["A", "AA"],
            "OA": ["G", "A"], "EAF": [0.2, 0.3], "BETA": [0.2, 0.4],
            "SE": [0.1, 0.1], "Z": [2.0, 4.0], "P": [0.05, 0.01],
        })
        vcf_frame = pl.DataFrame({
            "CHROM": [1, 1], "POS": [102, 199], "ID": ["snp", "indel"],
            "REF": ["G", "AA"], "ALT": ["A", "A"],
            "ES": [0.2, 0.4], "SE": [0.1, 0.1], "EZ": [2.0, 4.0],
            "AF": [0.2, 0.3],
        })

        result = compare_input_to_vcf(
            input_frame, vcf_frame, _row(), effect_type="beta",
            p_value_type="raw", eaf_is_maf=False, settings=_settings(),
            policies=_policies(),
        )

        self.assertEqual(result.summary["variant_types"]["snps"]["status"], "FAIL")
        self.assertFalse(
            result.summary["variant_types"]["indels"][
                "representation_adjustment_applied"
            ]
        )
        self.assertEqual(result.status, "FAIL")

    def test_does_not_pair_unmatched_indels_from_different_chromosomes(self):
        input_frame = pl.DataFrame({
            "CHR": [1, 1], "BP": [101, 200], "EA": ["A", "AA"],
            "OA": ["G", "A"], "EAF": [0.2, 0.3], "BETA": [0.2, 0.4],
            "SE": [0.1, 0.1], "Z": [2.0, 4.0], "P": [0.05, 0.01],
        })
        vcf_frame = pl.DataFrame({
            "CHROM": [1, 2], "POS": [101, 199], "ID": ["snp", "indel"],
            "REF": ["G", "AA"], "ALT": ["A", "A"],
            "ES": [0.2, 0.4], "SE": [0.1, 0.1], "EZ": [2.0, 4.0],
            "AF": [0.2, 0.3],
        })

        result = compare_input_to_vcf(
            input_frame, vcf_frame, _row(), effect_type="beta",
            p_value_type="raw", eaf_is_maf=False, settings=_settings(),
            policies=_policies(),
        )

        indels = result.summary["variant_types"]["indels"]
        self.assertEqual(indels["maximum_possible_representation_pairs"], 0)
        self.assertFalse(indels["representation_adjustment_applied"])
        self.assertEqual(result.status, "FAIL")

    def test_calculates_z_from_effect_and_two_sided_p_value(self):
        frame = pl.DataFrame({
            "CHR": [1], "BP": [101], "EA": ["A"], "OA": ["G"],
            "EAF": [0.2], "BETA": [1.0], "P": [0.04550026389635842],
        })
        prepared = prepare_input_table(
            frame,
            _row(standard_error_column=None, z_score_column=None),
            effect_type="beta",
            p_value_type="raw",
            eaf_is_maf=False,
            settings=_settings(),
            policies=_policies(),
        )
        self.assertAlmostEqual(prepared.item(0, "input_z"), 2.0, places=10)
        self.assertEqual(prepared.item(0, "input_z_source"), "calculated_effect_and_p")

    def test_equal_length_multibase_variant_is_not_classified_as_an_indel(self):
        input_frame = pl.DataFrame({
            "CHR": [1], "BP": [101], "EA": ["AC"], "OA": ["GT"],
            "EAF": [0.2], "BETA": [0.2], "SE": [0.1], "Z": [2.0],
            "P": [0.05],
        })
        vcf_frame = pl.DataFrame({
            "CHROM": [1], "POS": [101], "ID": ["mnv"], "REF": ["GT"],
            "ALT": ["AC"], "ES": [0.2], "SE": [0.1], "EZ": [2.0],
            "AF": [0.2],
        })

        result = compare_input_to_vcf(
            input_frame, vcf_frame, _row(), effect_type="beta",
            p_value_type="raw", eaf_is_maf=False, settings=_settings(),
            policies=_policies(),
        )

        self.assertEqual(result.summary["variant_types"]["indels"]["input_unique_variants"], 0)
        self.assertEqual(
            result.summary["variant_types"]["other_variants"][
                "exact_matched_variants"
            ],
            1,
        )

    def test_raw_p_value_of_one_produces_zero_z_without_artificial_clipping(self):
        frame = pl.DataFrame({
            "CHR": [1], "BP": [101], "EA": ["A"], "OA": ["G"],
            "EAF": [0.2], "BETA": [0.0], "P": [1.0],
        })
        prepared = prepare_input_table(
            frame,
            _row(standard_error_column=None, z_score_column=None),
            effect_type="beta",
            p_value_type="raw",
            eaf_is_maf=False,
            settings=_settings(),
            policies=_policies(),
        )
        self.assertEqual(prepared.item(0, "input_z"), 0.0)

    def test_external_frequency_is_oriented_to_the_input_effect_allele(self):
        input_frame = pl.DataFrame({
            "CHR": [1, 1], "BP": [101, 102], "EA": ["A", "T"],
            "OA": ["G", "C"], "BETA": [0.2, -0.3], "SE": [0.1, 0.1],
            "Z": [2.0, -3.0], "P": [0.05, 0.05],
        })
        external = pl.DataFrame({
            "CHROM": [1, 1], "POS": [101, 102], "ALT": ["A", "C"],
            "REF": ["G", "T"], "FREQ": [0.2, 0.3],
        })
        vcf_frame = pl.DataFrame({
            "CHROM": [1, 1], "POS": [101, 102], "ID": ["v1", "v2"],
            "REF": ["G", "T"], "ALT": ["A", "C"],
            "ES": [0.2, 0.3], "SE": [0.1, 0.1], "EZ": [2.0, 3.0],
            "AF": [0.2, 0.3],
        })
        row = _row(
            effect_allele_frequency_column=None,
            external_eaf_file=Path("frequencies.tsv"),
            external_eaf_column="FREQ",
        )
        config = load_configuration()

        result = compare_input_to_vcf(
            input_frame,
            vcf_frame,
            row,
            effect_type="beta",
            p_value_type="raw",
            eaf_is_maf=False,
            settings=config.modules.harmonisation.concordance_validation,
            policies=load_policies(config.modules.harmonisation.policies),
            external_eaf_frame=external,
            external_eaf_mapping=config.modules.harmonisation.external_eaf_mapping,
        )

        self.assertEqual(result.status, "PASS")
        self.assertEqual(result.metric_summary["allele_frequency"]["concordant"], 2)

    def test_reference_confirmed_eaf_is_compared_without_maf_folding(self):
        input_frame = pl.DataFrame({
            "CHR": [1], "BP": [101], "EA": ["A"], "OA": ["G"],
            "EAF": [0.2], "BETA": [0.2], "SE": [0.1], "Z": [2.0],
            "P": [0.05],
        })
        vcf_frame = pl.DataFrame({
            "CHROM": [1], "POS": [101], "ID": ["v1"], "REF": ["G"],
            "ALT": ["A"], "ES": [0.2], "SE": [0.1], "EZ": [2.0],
            "AF": [0.8],
        })

        result = compare_input_to_vcf(
            input_frame, vcf_frame, _row(), effect_type="beta",
            p_value_type="raw", eaf_is_maf=False, settings=_settings(),
            policies=_policies(),
        )

        self.assertEqual(result.metric_summary["allele_frequency"]["mismatches"], 1)
        self.assertEqual(result.status, "FAIL")

    def test_unresolved_frequency_type_skips_only_frequency_concordance(self):
        input_frame = pl.DataFrame({
            "CHR": [1], "BP": [101], "EA": ["A"], "OA": ["G"],
            "EAF": [0.2], "BETA": [0.2], "SE": [0.1], "Z": [2.0],
            "P": [0.05],
        })
        vcf_frame = pl.DataFrame({
            "CHROM": [1], "POS": [101], "ID": ["v1"], "REF": ["G"],
            "ALT": ["A"], "ES": [0.2], "SE": [0.1], "EZ": [2.0],
            "AF": [0.8],
        })

        result = compare_input_to_vcf(
            input_frame, vcf_frame, _row(), effect_type="beta",
            p_value_type="raw", eaf_is_maf=None, settings=_settings(),
            policies=_policies(),
        )

        self.assertEqual(result.metric_summary["effect"]["concordant"], 1)
        self.assertEqual(result.metric_summary["z_score"]["concordant"], 1)
        self.assertEqual(result.metric_summary["allele_frequency"]["checked"], 0)
        self.assertEqual(
            result.metric_summary["allele_frequency"]["frequency_type"],
            "unresolved",
        )
        self.assertEqual(result.status, "WARNING")

    def test_missing_input_effect_is_reported_as_warning_not_false_pass(self):
        input_frame = pl.DataFrame({
            "CHR": [1], "BP": [101], "EA": ["A"], "OA": ["G"],
            "EAF": [0.2], "Z": [2.0], "P": [0.05],
        })
        vcf_frame = pl.DataFrame({
            "CHROM": [1], "POS": [101], "ID": ["v1"], "REF": ["G"],
            "ALT": ["A"], "ES": [0.2], "SE": [0.1], "EZ": [2.0],
            "AF": [0.2],
        })
        result = compare_input_to_vcf(
            input_frame,
            vcf_frame,
            _row(effect_column=None, standard_error_column=None),
            effect_type="beta",
            p_value_type="raw",
            eaf_is_maf=False,
            settings=_settings(),
            policies=_policies(),
        )
        self.assertEqual(result.status, "WARNING")
        self.assertEqual(result.metric_summary["effect"]["checked"], 0)
        self.assertEqual(result.metric_summary["effect"]["unavailable_in_input"], 1)
        self.assertEqual(result.metric_summary["z_score"]["concordant"], 1)

    def test_uses_harmonisation_duplicate_action_instead_of_first_input_row(self):
        input_frame = pl.DataFrame({
            "CHR": [1, 1], "BP": [101, 101], "EA": ["A", "A"],
            "OA": ["G", "G"], "EAF": [0.2, 0.2], "BETA": [0.2, 0.2],
            "SE": [None, 0.1], "Z": [2.0, 2.0], "P": [0.05, 0.05],
        })
        duplicate_report = input_frame.with_columns(
            pl.Series("duplicate_input_row", [1, 2]),
            pl.Series("duplicate_action", ["removed", "kept"])
        )
        vcf_frame = pl.DataFrame({
            "CHROM": [1], "POS": [101], "ID": ["v1"], "REF": ["G"],
            "ALT": ["A"], "ES": [0.2], "SE": [0.1], "EZ": [2.0],
            "AF": [0.2],
        })

        result = compare_input_to_vcf(
            input_frame,
            vcf_frame,
            _row(),
            effect_type="beta",
            p_value_type="raw",
            eaf_is_maf=False,
            settings=_settings(),
            policies=_policies(),
            duplicate_report_frame=duplicate_report,
        )

        self.assertEqual(result.matched.item(0, "input_row"), 2)
        self.assertEqual(result.metric_summary["effect"]["concordant"], 1)
        self.assertEqual(result.summary["input_duplicate_rows"], 1)

    def test_duplicate_input_row_preserves_sample_size_selection(self):
        input_frame = pl.DataFrame({
            "CHR": [1, 1], "BP": [101, 101], "EA": ["A", "A"],
            "OA": ["G", "G"], "EAF": [0.2, 0.2], "BETA": [0.2, 0.2],
            "SE": [0.1, 0.1], "Z": [2.0, 2.0], "P": [0.05, 0.05],
            "N": [100, 1000],
        })
        duplicate_report = input_frame.with_columns(
            pl.Series("duplicate_input_row", [1, 2]),
            pl.Series("duplicate_action", ["removed", "kept"]),
        )
        vcf_frame = pl.DataFrame({
            "CHROM": [1], "POS": [101], "ID": ["v1"], "REF": ["G"],
            "ALT": ["A"], "ES": [0.2], "SE": [0.1], "EZ": [2.0],
            "AF": [0.2],
        })

        result = compare_input_to_vcf(
            input_frame, vcf_frame, _row(), effect_type="beta",
            p_value_type="raw", eaf_is_maf=False, settings=_settings(),
            policies=_policies(), duplicate_report_frame=duplicate_report,
        )

        self.assertEqual(result.matched.item(0, "input_row"), 2)

    def test_combined_coordinate_parser_matches_harmonisation(self):
        frame = pl.DataFrame({
            "VARIANT": [
                "1:100000", "1_100000", "1:100000:C:T",
                "1:100000.0", "1:100000-100001",
            ],
            "EA": ["A"] * 5, "OA": ["G"] * 5, "EAF": [0.2] * 5,
            "BETA": [0.2] * 5, "SE": [0.1] * 5, "Z": [2.0] * 5,
            "P": [0.05] * 5,
        })
        prepared = prepare_input_table(
            frame,
            _row(
                chromosome_column=None,
                position_column=None,
                chromosome_position_column="VARIANT",
            ),
            effect_type="beta", p_value_type="raw", eaf_is_maf=False,
            settings=_settings(), policies=_policies(),
        )

        self.assertEqual(prepared["input_chrom"].to_list(), ["1"] * 5)
        self.assertEqual(prepared["input_pos"].to_list(), [100000] * 5)
        self.assertEqual(prepared["valid_input_variant"].to_list(), [True] * 5)

    def test_duplicate_vcf_records_fail_at_the_configured_default(self):
        input_frame = pl.DataFrame({
            "CHR": [1], "BP": [101], "EA": ["A"], "OA": ["G"],
            "EAF": [0.2], "BETA": [0.2], "SE": [0.1], "Z": [2.0],
            "P": [0.05],
        })
        vcf_frame = pl.DataFrame({
            "CHROM": [1, 1], "POS": [101, 101], "ID": ["v1", "v1-copy"],
            "REF": ["G", "G"], "ALT": ["A", "A"], "ES": [0.2, 0.2],
            "SE": [0.1, 0.1], "EZ": [2.0, 2.0], "AF": [0.2, 0.2],
        })
        result = compare_input_to_vcf(
            input_frame,
            vcf_frame,
            _row(),
            effect_type="beta",
            p_value_type="raw",
            eaf_is_maf=False,
            settings=_settings(),
            policies=_policies(),
        )
        self.assertEqual(result.status, "FAIL")
        self.assertEqual(result.summary["vcf_duplicate_records"], 1)
        self.assertEqual(result.vcf_duplicates.height, 1)

    def test_service_streams_and_combines_chromosome_partitions(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "input.tsv"
            vcf_path = root / "study.vcf.gz"
            manifest_path = root / "manifest.json"
            input_frame = pl.DataFrame({
                "CHR": [2, 1, 2],
                "BP": [202, 101, 203],
                "EA": ["A", "C", "T"],
                "OA": ["G", "T", "C"],
                "EAF": [0.2, 0.3, 0.4],
                "BETA": [0.2, -0.3, 0.4],
                "SE": [0.1, 0.1, 0.1],
                "Z": [2.0, -3.0, 4.0],
                "P": [0.05, 0.01, 0.001],
            })
            input_frame.write_csv(input_path, separator="\t")
            vcf_path.write_bytes(b"placeholder")
            manifest_path.write_text(
                '{"dataset":{"genome_build":{"inferred_build":"GRCh37"},'
                '"study_decisions":{"effect_type":"beta","pvalue_type":"raw",'
                '"eaf_is_maf":false,"eaf_is_maf_source":"study_level_statistic"}}}',
                encoding="utf-8",
            )
            extracted = pl.DataFrame({
                "CHROM": [1, 2, 2],
                "POS": [101, 202, 203],
                "ID": ["v1", "v2", "v3"],
                "REF": ["T", "G", "C"],
                "ALT": ["C", "A", "T"],
                "ES": [-0.3, 0.2, 0.4],
                "SE": [0.1, 0.1, 0.1],
                "EZ": [-3.0, 2.0, 4.0],
                "AF": [0.3, 0.2, 0.4],
            })

            def write_extracted_table(
                _vcf, table, _dataset, _columns, _bcftools, **kwargs,
            ):
                extracted.write_csv(table, separator=kwargs["delimiter"])
                return str(table)

            configuration = load_configuration()
            settings = _settings().model_copy(update={"write_all_matches": True})
            expected = compare_input_to_vcf(
                input_frame,
                extracted,
                _row(input_file=input_path),
                effect_type="beta",
                p_value_type="raw",
                eaf_is_maf=False,
                settings=settings,
                policies=_policies(),
            )
            with patch(
                "postgwas.modules.harmonisation.concordance.service.run_checked_command",
                return_value="##contig=<ID=1,assembly=GRCh37>\n",
            ), patch(
                "postgwas.modules.harmonisation.concordance.service.extract_vcf_table",
                side_effect=write_extracted_table,
            ), patch("polars.read_csv", side_effect=AssertionError("eager CSV read")):
                result = run_concordance_validation(
                    row=_row(input_file=input_path),
                    vcf_path=vcf_path,
                    output_root=root / "results",
                    settings=settings,
                    policies=_policies(),
                    threads=1,
                    bcftools="bcftools",
                    output_layout=dict(
                        configuration.modules.harmonisation.output_layout.root
                    ),
                    vcf_config=(
                        configuration.modules.harmonisation.vcf_processing.model_dump()
                    ),
                    file_log_level=configuration.logging.file_level,
                    screen_log_level=configuration.logging.console_level,
                    run_manifest=manifest_path,
                    show_screen=False,
                )

            self.assertEqual(result["status"], "PASS")
            self.assertEqual(result["summary"], expected.summary)
            self.assertEqual(result["metrics"], expected.metric_summary)
            self.assertEqual(result["summary"]["matched_variants"], 3)
            self.assertEqual(result["metrics"]["effect"]["concordant"], 3)
            self.assertAlmostEqual(result["metrics"]["effect"]["pearson"], 1.0)
            self.assertAlmostEqual(result["metrics"]["effect"]["spearman"], 1.0)
            matches = pl.read_csv(
                result["reports"]["all_matches"], separator="\t",
            )
            self.assertEqual(matches["input_row"].to_list(), [1, 2, 3])
            log = Path(result["reports"]["log"]).read_text(encoding="utf-8")
            self.assertIn("2 canonical chromosome partition(s)", log)
            concordance_directory = Path(result["reports"]["summary"]).parent
            self.assertFalse(any(
                path.is_dir()
                for path in concordance_directory.glob(".study_vcf_*")
            ))

    def test_failure_still_writes_a_complete_log_and_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            row = _row(input_file=root / "input.tsv")
            configuration = load_configuration()
            with self.assertRaisesRegex(ConcordanceValidationError, "does not exist"):
                run_concordance_validation(
                    row=row,
                    vcf_path=root / "missing.vcf.gz",
                    output_root=root / "results",
                    settings=_settings(),
                    policies=_policies(),
                    threads=1,
                    bcftools="bcftools",
                    output_layout=dict(
                        configuration.modules.harmonisation.output_layout.root
                    ),
                    vcf_config=configuration.modules.harmonisation.vcf_processing.model_dump(),
                    file_log_level=configuration.logging.file_level,
                    screen_log_level=configuration.logging.console_level,
                    show_screen=False,
                )
            base = (
                root / "results" / "study" / "harmonisation"
                / "00_harmonised_sumstat"
            )
            log = (base / "logs" / "study_concordance_validation.log").read_text()
            summary = (
                base / "qc_summary" / "concordance"
                / "study_concordance_summary.tsv"
            ).read_text()
            self.assertIn("Validation failed", log)
            self.assertIn("Concordance validation finished", log)
            self.assertIn("\tstatus\tFAIL", summary)


if __name__ == "__main__":
    unittest.main()
