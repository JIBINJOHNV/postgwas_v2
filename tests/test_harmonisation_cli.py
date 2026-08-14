"""CLI/config boundary tests; no scientific work is started."""

import argparse
import csv
import io
import json
import tempfile
import unittest
from argparse import Namespace
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import yaml

from postgwas.core.errors import ConfigurationError
from postgwas.core.paths import configured_output_path
from postgwas.cli.compute import get_compute_parser, memory_limit, resolve_compute_args
from postgwas.modules.harmonisation.cli import (
    RunSummaryError,
    ScreenReportError,
    _run_validated_rows,
    get_harmonisation_parser,
    run_harmonisation,
)
from postgwas.modules.harmonisation.service import PipelineError
from postgwas.config import load_configuration
from postgwas.modules.harmonisation.sample_sheet import load_harmonisation_sample_sheet


def _successful_engine_result(**kwargs):
    """Write the manifest promised by the mocked successful engine."""
    engine_input = kwargs["sample_column_dict"]
    output_layout = kwargs["default_cfg"]["output_layout"]
    manifest_path = configured_output_path(
        engine_input["output_folder"],
        output_layout["run_manifest"],
        dataset_id=engine_input["gwas_outputname"],
    )
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    assessment_path = manifest_path.parent / "qc_assessment.json"
    assessment_path.write_text(
        json.dumps({
            "raw": {
                "num_records": 0,
                "num_snps": 0,
                "num_non_snps": 0,
                "af_comparable": 0,
                "af_difference_above_cutoff": 0,
                "study_af_missing": 0,
                "external_af_missing": 0,
            },
            "qc_passed": {
                "num_records": 0,
                "num_snps": 0,
                "effective_sample_size_missing_or_invalid": 0,
                "format_si_missing": 0,
                "effective_sample_size_above_outlier_threshold": 0,
            },
            "active_rule_count": 0,
            "af_difference_cutoff": 0.2,
        }),
        encoding="utf-8",
    )
    manifest_path.write_text(
        json.dumps({
            "status": "OK",
            "pre_vcf_summary": {
                "total_variant_with_invalid_beta_se": 0,
            },
            "qc_assessment": {
                "raw_variants": 0,
                "qc_passed_variants": 0,
                "reports": {"json": str(assessment_path)},
            },
            "dataset": {
                "genome_build": {
                    "inferred_build": "GRCh37",
                    "forced": False,
                },
                "study_decisions": {
                    "effect_type": "beta",
                    "effect_type_detected": "beta",
                    "effect_type_source": "detector",
                    "pvalue_type": "raw",
                    "pvalue_type_detected": "raw",
                    "pvalue_type_source": "detector",
                },
                "completed": [],
                "failed": [],
                "chromosomes": {},
                "chromosome_summaries": {},
            },
        }),
        encoding="utf-8",
    )
    return {"status": "OK", "manifest": str(manifest_path)}


class HarmonisationCLITests(unittest.TestCase):
    def test_configurable_cli_actions_do_not_own_defaults(self):
        parser = get_harmonisation_parser()
        configurable = {
            "sample_sheet",
            "dataset_id",
            "run_config",
            "resource_directory",
            "output_directory",
            "comparison_af_source",
            "comparison_af_column",
            "threads",
            "memory_gb",
            "seed",
            "validate",
            "fixed_info",
            "zero_p_se_action",
            "show_screen",
        }
        actions = {action.dest: action for action in parser._actions}
        for destination in configurable:
            self.assertEqual(actions[destination].default, argparse.SUPPRESS)

    def test_help_displays_internal_defaults(self):
        help_text = get_harmonisation_parser().format_help()
        normalized_help = " ".join(help_text.split())
        self.assertIn("Default: ALFA", help_text)
        self.assertIn("Default: EUR", help_text)
        self.assertIn("PostGWAS chooses automatically", help_text)
        self.assertIn("available physical memory", normalized_help)
        self.assertIn("REQUIRED", help_text)
        self.assertIn("Options: ALFA, 1000G", help_text)
        self.assertIn("GRCh37_ALFA_freq_chr[1..22,X,Y].vcf.gz", help_text)
        self.assertIn("GRCh37_1000G_freq_chr[1..22,X,Y].vcf.gz", help_text)
        self.assertIn("Disabled by default", help_text)
        self.assertIn("--zero-p-se-action", help_text)
        self.assertIn("Default: fail", help_text)
        self.assertIn("--fixed-info", help_text)
        self.assertIn("between 0 and 1", help_text)
        self.assertIn("--show-screen", help_text)
        self.assertIn("--hide-screen", help_text)
        self.assertIn("Default: true", help_text)

    def test_screen_display_cli_flags_are_mutually_exclusive(self):
        parser = get_harmonisation_parser()
        self.assertFalse(hasattr(parser.parse_args([]), "show_screen"))
        self.assertTrue(parser.parse_args(["--show-screen"]).show_screen)
        self.assertFalse(parser.parse_args(["--hide-screen"]).show_screen)
        with self.assertRaises(SystemExit):
            parser.parse_args(["--show-screen", "--hide-screen"])

    def test_comparison_af_source_accepts_only_supported_panels(self):
        parser = get_harmonisation_parser()
        self.assertEqual(
            parser.parse_args(["--comparison-af-source", "ALFA"]).comparison_af_source,
            "ALFA",
        )
        with self.assertRaises(SystemExit):
            parser.parse_args(["--comparison-af-source", "unsupported"])

    def test_fixed_info_is_validated_by_the_cli_parser(self):
        parser = get_harmonisation_parser()
        for supplied, expected in (("0", 0.0), ("0.99", 0.99), ("1", 1.0)):
            with self.subTest(supplied=supplied):
                self.assertEqual(
                    parser.parse_args(["--fixed-info", supplied]).fixed_info,
                    expected,
                )
        for invalid in ("-0.01", "1.01", "nan", "inf", "not-a-number"):
            with self.subTest(invalid=invalid), self.assertRaises(SystemExit):
                parser.parse_args(["--fixed-info", invalid])

    def test_legacy_flags_are_not_part_of_the_parser(self):
        options = {
            option
            for action in get_harmonisation_parser()._actions
            for option in action.option_strings
        }
        self.assertFalse(
            {"--config", "--defaults", "--parameter-config", "--nthreads", "--max-mem"}
            & options
        )

    def test_compute_parser_is_shared_and_config_driven(self):
        parser = get_compute_parser()
        actions = {action.dest: action for action in parser._actions}
        self.assertEqual(actions["threads"].default, argparse.SUPPRESS)
        self.assertEqual(actions["memory_gb"].default, argparse.SUPPRESS)
        self.assertEqual(actions["seed"].default, argparse.SUPPRESS)
        args = resolve_compute_args(Namespace())
        config = load_configuration()
        self.assertEqual(args.threads, config.execution.threads)
        self.assertEqual(args.memory_gb, config.execution.memory_gb)
        self.assertEqual(args.seed, config.execution.random_seed)
        self.assertEqual(memory_limit(12.5), "12.5G")

    def test_sample_sheet_is_required_in_new_mode(self):
        with self.assertRaisesRegex(ConfigurationError, "sample-sheet is required"):
            run_harmonisation(Namespace(output_directory="results"))

    def test_explicit_cli_values_are_passed_to_configuration_loader(self):
        args = Namespace(
            sample_sheet="sheet.csv",
            resource_directory="resources",
            output_directory="output",
            comparison_af_source="ALFA",
            comparison_af_column="AFR",
            threads=4,
            memory_gb=20,
            fixed_info=0.99,
            zero_p_se_action="approximate",
            show_screen=False,
        )
        fake_config = unittest.mock.Mock()
        fake_config.run.output_directory = "output"
        fake_config.logging.filename = "postgwas.log"
        with patch(
            "postgwas.modules.harmonisation.cli.load_configuration",
            return_value=fake_config,
        ) as load, patch(
            "postgwas.modules.harmonisation.cli.require_binaries",
        ), patch(
            "postgwas.modules.harmonisation.cli.load_harmonisation_sample_sheet",
            return_value=[],
        ), patch(
            "postgwas.modules.harmonisation.cli._run_validated_rows",
            return_value={},
        ):
            run_harmonisation(args)
        overrides = load.call_args.kwargs["cli_overrides"]
        self.assertEqual(overrides["resources.root"], "resources")
        self.assertEqual(overrides["execution.threads"], 4)
        self.assertEqual(
            overrides["modules.harmonisation.comparison_af.column"], "AFR"
        )
        self.assertEqual(
            overrides[
                "modules.harmonisation.policies.pvalue.zero_missing_se"
            ],
            "approximate",
        )
        self.assertEqual(
            overrides["modules.harmonisation.fixed_info.value"], 0.99
        )
        self.assertFalse(
            overrides["logging.show_screen"]
        )

    def test_run_configuration_cannot_activate_fixed_info_without_cli(self):
        fixture = (
            Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        )
        rows = load_harmonisation_sample_sheet(fixture)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            config = load_configuration(cli_overrides={
                "resources.root": str(resources),
                "run.output_directory": str(root / "output"),
                "modules.harmonisation.fixed_info.value": 0.99,
            })
            with self.assertRaisesRegex(ConfigurationError, "--fixed-info"):
                _run_validated_rows(Namespace(), config, rows)

    def test_validated_rows_write_metadata_before_engine_execution(self):
        fixture = (
            Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        )
        rows = load_harmonisation_sample_sheet(fixture)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            output = root / "output"
            config = load_configuration(
                cli_overrides={
                    "resources.root": str(resources),
                    "run.output_directory": str(output),
                    "execution.retries": 0,
                }
            )
            with patch(
                "postgwas.modules.harmonisation.cli.run_harmonisation_pipeline",
                side_effect=_successful_engine_result,
            ) as engine:
                result = _run_validated_rows(Namespace(), config, rows)
            metadata = output / rows[0].dataset_id / "run_metadata"
            self.assertTrue((metadata / "resolved_config.yaml").is_file())
            resolved = yaml.safe_load(
                (metadata / "resolved_config.yaml").read_text(encoding="utf-8")
            )
            self.assertEqual(
                list(resolved["modules"]), ["harmonisation", "qc_summary"]
            )
            self.assertEqual(
                set(resolved["resources"]), {"root", "executables"}
            )
            self.assertTrue((metadata / "sample_sheet_row.json").is_file())
            self.assertTrue((metadata / "command.txt").is_file())
            self.assertTrue((metadata / config.logging.filename).is_file())
            run_log = output / "run_metadata" / config.logging.filename
            self.assertIn(
                "Run completed successfully: 1 dataset(s).",
                run_log.read_text(encoding="utf-8"),
            )
            self.assertEqual(result[rows[0].dataset_id]["status"], "OK")
            defaults = engine.call_args.kwargs["default_cfg"]
            self.assertNotIn("default_eaf", defaults)
            self.assertNotIn("default_info", defaults)
            self.assertIn("resource_layout", defaults)
            engine_input = engine.call_args.kwargs["sample_column_dict"]
            self.assertIn("resource_folder", engine_input)
            self.assertNotIn("resourse_folder", engine_input)

    def test_each_selected_dataset_is_announced_with_position_and_input(self):
        fixture = (
            Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        )
        base = load_harmonisation_sample_sheet(fixture)[0]
        rows = [
            base.model_copy(update={
                "dataset_id": "study_one",
                "input_file": "/data/study_one.tsv.gz",
            }),
            base.model_copy(update={
                "dataset_id": "study_two",
                "input_file": "/data/study_two.tsv.gz",
            }),
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            output = root / "output"
            config = load_configuration(cli_overrides={
                "resources.root": str(resources),
                "run.output_directory": str(output),
                "execution.retries": 0,
            })
            screen = io.StringIO()
            def engine(**kwargs):
                print(
                    "ENGINE OUTPUT %s"
                    % kwargs["sample_column_dict"]["gwas_outputname"]
                )
                return _successful_engine_result(**kwargs)

            with patch(
                "postgwas.modules.harmonisation.cli.run_harmonisation_pipeline",
                side_effect=engine,
            ), redirect_stdout(screen):
                result = _run_validated_rows(Namespace(), config, rows)

            text = screen.getvalue()
            first = "Starting dataset 1 of 2 — study_one"
            second = "Starting dataset 2 of 2 — study_two"
            self.assertLess(text.index(first), text.index(second))
            self.assertIn("Summary statistics   : /data/study_one.tsv.gz", text)
            self.assertIn("Summary statistics   : /data/study_two.tsv.gz", text)
            self.assertEqual(set(result), {"study_one", "study_two"})
            run_log = (
                output / "run_metadata" / config.logging.filename
            ).read_text(encoding="utf-8")
            self.assertIn(
                "Starting dataset 1 of 2: dataset_id=study_one", run_log,
            )
            self.assertIn(
                "Starting dataset 2 of 2: dataset_id=study_two", run_log,
            )
            with (output / "run_metadata" / "harmonisation_run_summary.csv").open(
                newline="", encoding="utf-8",
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            self.assertEqual(
                [(item["dataset_id"], item["status"]) for item in summary_rows],
                [("study_one", "OK"), ("study_two", "OK")],
            )
            for row, summary_row in zip(rows, summary_rows):
                report_path = (
                    output
                    / row.dataset_id
                    / "harmonisation"
                    / (row.dataset_id + "_screen_report.txt")
                ).resolve()
                self.assertEqual(summary_row["screen_report"], str(report_path))
                self.assertEqual(
                    result[row.dataset_id]["screen_report"], str(report_path)
                )
                report = report_path.read_text(encoding="utf-8")
                self.assertIn("Preparing PostGWAS harmonisation", report)
                self.assertIn("Starting dataset", report)
                self.assertIn("ENGINE OUTPUT %s" % row.dataset_id, report)
                self.assertIn("Harmonisation run summary", report)
                self.assertIn("Harmonisation complete — 2 dataset(s)", report)
                other = "study_two" if row.dataset_id == "study_one" else "study_one"
                self.assertNotIn("ENGINE OUTPUT %s" % other, report)
                self.assertNotIn("— %s" % other, report)

    def test_hide_screen_suppresses_stdout_but_still_writes_complete_report(self):
        fixture = (
            Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        )
        rows = load_harmonisation_sample_sheet(fixture)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            output = root / "output"
            config = load_configuration(cli_overrides={
                "resources.root": str(resources),
                "run.output_directory": str(output),
                "execution.retries": 0,
                "logging.show_screen": False,
            })

            def engine(**kwargs):
                print("ENGINE OUTPUT WHILE HIDDEN")
                return _successful_engine_result(**kwargs)

            screen = io.StringIO()
            with patch(
                "postgwas.modules.harmonisation.cli.run_harmonisation_pipeline",
                side_effect=engine,
            ), redirect_stdout(screen):
                result = _run_validated_rows(Namespace(), config, rows)

            self.assertEqual(screen.getvalue(), "")
            report_path = Path(result[rows[0].dataset_id]["screen_report"])
            report = report_path.read_text(encoding="utf-8")
            self.assertIn("Preparing PostGWAS harmonisation", report)
            self.assertIn("ENGINE OUTPUT WHILE HIDDEN", report)
            self.assertIn("Harmonisation run summary", report)
            self.assertIn("Harmonisation complete — 1 dataset(s)", report)

    def test_screen_report_initialization_failure_stops_before_analysis(self):
        fixture = (
            Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        )
        rows = load_harmonisation_sample_sheet(fixture)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            config = load_configuration(cli_overrides={
                "resources.root": str(resources),
                "run.output_directory": str(root / "output"),
                "execution.retries": 0,
            })
            with patch.object(Path, "write_text", side_effect=OSError("read only")):
                with self.assertRaisesRegex(ScreenReportError, "screen report"):
                    _run_validated_rows(Namespace(), config, rows)

    def test_run_summary_combines_study_and_chromosome_evidence(self):
        fixture = (
            Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        )
        rows = load_harmonisation_sample_sheet(fixture)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            output = root / "output"
            config = load_configuration(cli_overrides={
                "resources.root": str(resources),
                "run.output_directory": str(output),
                "execution.retries": 0,
            })

            def engine(**kwargs):
                result = _successful_engine_result(**kwargs)
                manifest_path = Path(result["manifest"])
                document = json.loads(manifest_path.read_text(encoding="utf-8"))
                document.update({
                    "pre_vcf_summary": {
                        "total_variant_infile": 30,
                        "total_variant_remaining_for_harmonisation": 30,
                        "total_variant_ready_snps": 28,
                        "total_variant_ready_indels_or_other": 2,
                        "total_variant_passed_chromosome_harmonisation": 25,
                        "total_variant_with_invalid_beta_se": 3,
                    },
                    "population_frequency_qc": {"closest_population": "EUR"},
                })
                document["qc_assessment"].update({
                    "raw_variants": 25,
                    "qc_passed_variants": 24,
                })
                Path(document["qc_assessment"]["reports"]["json"]).write_text(
                    json.dumps({
                        "raw": {
                            "num_records": 25,
                            "num_snps": 23,
                            "num_non_snps": 2,
                            "af_comparable": 20,
                            "af_difference_above_cutoff": 2,
                            "study_af_missing": 1,
                            "external_af_missing": 5,
                        },
                        "qc_passed": {
                            "num_records": 24,
                            "num_snps": 22,
                            "effective_sample_size_missing_or_invalid": 0,
                            "format_si_missing": 0,
                            "effective_sample_size_above_outlier_threshold": 1,
                        },
                        "active_rule_count": 6,
                        "af_difference_cutoff": 0.2,
                    }),
                    encoding="utf-8",
                )
                document["dataset"].update({
                    "genome_build": {
                        "inferred_build": "GRCh38",
                        "forced": False,
                    },
                    "study_decisions": {
                        "effect_type": "odds_ratio",
                        "effect_type_detected": "odds_ratio",
                        "effect_type_source": "detector",
                        "pvalue_type": "neglog10",
                        "pvalue_type_detected": "neglog10",
                        "pvalue_type_source": "detector",
                        "frequency_type": "effect_allele_frequency",
                        "eaf_is_maf_source": "reference_confirmed_eaf",
                        "strand": "forward",
                    },
                    "completed": ["1", "2"],
                    "failed": [],
                    "chromosomes": {
                        "1": {"status": "ok"},
                        "2": {"status": "ok"},
                    },
                    "reconciliation": {"rows_exported": 25},
                    "chromosome_summaries": {
                        chromosome: {
                            "stage_qc": {
                                "eaf_qc": {
                                    "strand_orientation": orientation
                                }
                            }
                        }
                        for chromosome, orientation in {
                            "1": {
                                "status": "success",
                                "initial_variants": 10,
                                "final_variants": 8,
                                "reference_file": "/reference/panel_chr1.tsv.gz",
                                "reference_population_column": "EUR",
                                "study_strand_consensus": "forward",
                                "actions": {"forward": 7, "forward_swapped": 1},
                                "reference_unmatched": 1,
                                "palindromic_ambiguous": 1,
                                "reference_ambiguous": 0,
                            },
                            "2": {
                                "status": "success",
                                "initial_variants": 20,
                                "final_variants": 17,
                                "reference_file": "/reference/panel_chr2.tsv.gz",
                                "reference_population_column": "EUR",
                                "study_strand_consensus": "forward",
                                "actions": {
                                    "forward": 15,
                                    "forward_swapped": 1,
                                    "reverse_complement_swapped": 1,
                                },
                                "reference_unmatched": 1,
                                "palindromic_ambiguous": 1,
                                "reference_ambiguous": 1,
                            },
                        }.items()
                    },
                })
                manifest_path.write_text(json.dumps(document), encoding="utf-8")
                return result

            with patch(
                "postgwas.modules.harmonisation.cli.run_harmonisation_pipeline",
                side_effect=engine,
            ):
                _run_validated_rows(Namespace(), config, rows)

            with (output / "run_metadata" / "harmonisation_run_summary.csv").open(
                newline="", encoding="utf-8",
            ) as handle:
                record = next(csv.DictReader(handle))
            self.assertEqual(record["genome_build"], "GRCh38")
            self.assertEqual(record["genome_build_detected"], "GRCh38")
            self.assertEqual(record["effect_type"], "odds_ratio")
            self.assertEqual(record["effect_type_automatically_detected"], "true")
            self.assertEqual(
                record["effect_scale"],
                "odds ratio harmonised to log-odds BETA",
            )
            self.assertEqual(record["p_value_type"], "neglog10")
            self.assertEqual(record["p_value_type_automatically_detected"], "true")
            self.assertEqual(record["strand_chromosomes_summarized"], "2")
            self.assertEqual(record["strand_metadata_complete"], "true")
            self.assertEqual(record["strand_variants_evaluated"], "30")
            self.assertEqual(record["strand_matched"], "25")
            self.assertEqual(record["strand_forward"], "22")
            self.assertEqual(record["strand_forward_swapped"], "2")
            self.assertEqual(record["strand_reverse_complement_swapped"], "1")
            self.assertEqual(record["strand_removed_total"], "5")
            self.assertEqual(record["strand_accounting_balanced"], "true")
            self.assertEqual(record["strand_reference_file_count"], "2")
            self.assertEqual(
                record["strand_reference_files"],
                '["/reference/panel_chr1.tsv.gz","/reference/panel_chr2.tsv.gz"]',
            )
            self.assertEqual(record["strand_reference_panel"], "ALFA")
            self.assertEqual(record["strand_reference_population"], "EUR")
            self.assertEqual(record["closest_population"], "EUR")
            self.assertEqual(record["comparison_af_panel"], "ALFA")
            self.assertEqual(record["comparison_af_population"], "EUR")
            self.assertEqual(record["input_variants"], "30")
            self.assertEqual(record["input_ready_variants"], "30")
            self.assertEqual(record["input_ready_snps"], "28")
            self.assertEqual(
                record["input_ready_indels_or_other_variants"], "2",
            )
            self.assertEqual(record["harmonised_variants"], "25")
            self.assertEqual(record["final_vcf_variants"], "25")
            self.assertEqual(record["qc_passed_variants"], "24")
            self.assertEqual(record["final_vcf_snps"], "23")
            self.assertEqual(record["final_vcf_indels_or_other_variants"], "2")
            self.assertEqual(record["qc_active_rule_count"], "6")
            self.assertEqual(record["qc_passed_snps"], "22")
            self.assertEqual(record["qc_failed_snps"], "1")
            self.assertAlmostEqual(float(record["qc_passed_snp_percent"]), 95.6521739)
            self.assertAlmostEqual(float(record["qc_failed_snp_percent"]), 4.3478261)
            self.assertEqual(record["af_comparable_variants"], "20")
            self.assertEqual(record["af_concordant_variants"], "18")
            self.assertEqual(float(record["af_concordance_percent"]), 90.0)
            self.assertEqual(record["af_mismatched_variants"], "2")
            self.assertEqual(float(record["af_mismatch_percent"]), 10.0)
            self.assertEqual(record["study_af_missing"], "1")
            self.assertEqual(record["reference_af_missing"], "5")
            self.assertEqual(float(record["af_difference_cutoff"]), 0.2)
            self.assertEqual(
                record["pre_vcf_invalid_effect_statistic_removals"], "3",
            )
            self.assertEqual(record["qc_passed_missing_or_invalid_neff"], "0")
            self.assertEqual(record["qc_passed_missing_imputation_score"], "0")
            self.assertEqual(record["qc_passed_neff_upper_outliers"], "1")
            fields = list(record)
            self.assertLess(fields.index("input_ready_snps"), fields.index("genome_build"))
            self.assertLess(fields.index("genome_build"), fields.index("effect_type"))
            self.assertLess(fields.index("effect_type"), fields.index("strand_status"))
            self.assertLess(fields.index("strand_status"), fields.index("harmonised_variants"))
            self.assertLess(
                fields.index("harmonised_variants"),
                fields.index("closest_population"),
            )
            self.assertLess(fields.index("closest_population"), fields.index("final_vcf_variants"))

    def test_completed_dataset_without_qc_json_fails_required_run_summary(self):
        fixture = (
            Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        )
        rows = load_harmonisation_sample_sheet(fixture)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            output = root / "output"
            config = load_configuration(cli_overrides={
                "resources.root": str(resources),
                "run.output_directory": str(output),
                "execution.retries": 0,
            })

            def incomplete_engine(**kwargs):
                result = _successful_engine_result(**kwargs)
                manifest_path = Path(result["manifest"])
                document = json.loads(manifest_path.read_text(encoding="utf-8"))
                document["qc_assessment"].pop("reports")
                manifest_path.write_text(json.dumps(document), encoding="utf-8")
                return result

            with patch(
                "postgwas.modules.harmonisation.cli.run_harmonisation_pipeline",
                side_effect=incomplete_engine,
            ):
                with self.assertRaisesRegex(
                    RunSummaryError, "required QC assessment JSON",
                ):
                    _run_validated_rows(Namespace(), config, rows)

    def test_partial_dataset_is_not_reported_as_successful(self):
        fixture = (
            Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        )
        rows = load_harmonisation_sample_sheet(fixture)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            output = root / "output"
            config = load_configuration(cli_overrides={
                "resources.root": str(resources),
                "run.output_directory": str(output),
                "execution.retries": 0,
            })

            def partial_engine(**kwargs):
                result = _successful_engine_result(**kwargs)
                manifest_path = Path(result["manifest"])
                document = json.loads(manifest_path.read_text(encoding="utf-8"))
                document["status"] = "PARTIAL"
                manifest_path.write_text(json.dumps(document), encoding="utf-8")
                result["status"] = "PARTIAL"
                return result

            screen = io.StringIO()
            with patch(
                "postgwas.modules.harmonisation.cli.run_harmonisation_pipeline",
                side_effect=partial_engine,
            ), redirect_stdout(screen):
                result = _run_validated_rows(Namespace(), config, rows)

            self.assertEqual(result[rows[0].dataset_id]["status"], "PARTIAL")
            self.assertIn("completed with status PARTIAL", screen.getvalue())
            self.assertNotIn("Dataset completed successfully", screen.getvalue())
            with (output / "run_metadata" / "harmonisation_run_summary.csv").open(
                newline="", encoding="utf-8",
            ) as handle:
                record = next(csv.DictReader(handle))
            self.assertEqual(record["status"], "PARTIAL")
            run_log = (
                output / "run_metadata" / config.logging.filename
            ).read_text(encoding="utf-8")
            self.assertIn("Run completed with 0 OK and 1 partial", run_log)
            self.assertNotIn("Run completed successfully", run_log)

    def test_integrated_validation_uses_same_build_merged_vcf(self):
        fixture = Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        rows = load_harmonisation_sample_sheet(fixture)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            output = root / "output"
            manifest = root / "run_manifest.json"
            assessment = root / "qc_assessment.json"
            assessment.write_text(
                json.dumps({
                    "raw": {
                        "num_records": 0,
                        "num_snps": 0,
                        "num_non_snps": 0,
                        "af_comparable": 0,
                        "af_difference_above_cutoff": 0,
                        "study_af_missing": 0,
                        "external_af_missing": 0,
                    },
                    "qc_passed": {
                        "num_records": 0,
                        "num_snps": 0,
                        "effective_sample_size_missing_or_invalid": 0,
                        "format_si_missing": 0,
                        "effective_sample_size_above_outlier_threshold": 0,
                    },
                    "active_rule_count": 0,
                    "af_difference_cutoff": 0.2,
                }),
                encoding="utf-8",
            )
            manifest.write_text(
                json.dumps({
                    "dataset": {
                        "genome_build": {"inferred_build": "GRCh37"},
                    },
                    "qc_assessment": {
                        "reports": {"json": str(assessment)},
                    },
                }),
                encoding="utf-8",
            )
            config = load_configuration(cli_overrides={
                "resources.root": str(resources),
                "run.output_directory": str(output),
                "execution.retries": 0,
                "modules.harmonisation.concordance_validation.enabled": True,
            })
            harmonisation_result = {
                "GRCh37": str(root / "study_GRCh37_merged.vcf.gz"),
                "GRCh38": str(root / "study_GRCh38_merged.vcf.gz"),
                "manifest": str(manifest),
                "status": "OK",
            }
            validation_result = {
                "status": "PASS", "reports": {"summary": str(root / "summary.tsv")}
            }
            with patch(
                "postgwas.modules.harmonisation.cli.run_harmonisation_pipeline",
                return_value=harmonisation_result,
            ), patch(
                "postgwas.modules.harmonisation.concordance.service.run_concordance_validation",
                return_value=validation_result,
            ) as validation:
                result = _run_validated_rows(Namespace(), config, rows)

            call = validation.call_args.kwargs
            self.assertEqual(call["vcf_path"], harmonisation_result["GRCh37"])
            self.assertEqual(call["run_manifest"], harmonisation_result["manifest"])
            self.assertIsNotNone(call["external_eaf_mapping"])
            self.assertEqual(
                result[rows[0].dataset_id]["concordance_validation"]["status"],
                "PASS",
            )

    def test_engine_failure_finishes_both_logs_with_actionable_context(self):
        fixture = Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        rows = load_harmonisation_sample_sheet(fixture)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            output = root / "output"
            config = load_configuration(
                cli_overrides={
                    "resources.root": str(resources),
                    "run.output_directory": str(output),
                    "execution.retries": 0,
                }
            )
            with patch(
                "postgwas.modules.harmonisation.cli.run_harmonisation_pipeline",
                side_effect=RuntimeError("adapter reported a scientific failure"),
            ):
                with self.assertRaises(PipelineError):
                    _run_validated_rows(Namespace(), config, rows)

            dataset_log = (
                output / rows[0].dataset_id / "run_metadata" / config.logging.filename
            ).read_text(encoding="utf-8")
            run_log = (
                output / "run_metadata" / config.logging.filename
            ).read_text(encoding="utf-8")
            self.assertIn("Dataset failed permanently", dataset_log)
            self.assertIn("adapter reported a scientific failure", dataset_log)
            self.assertIn("Traceback", dataset_log)
            self.assertIn("Run failed: all 1 dataset(s) failed", run_log)
            with (output / "run_metadata" / "harmonisation_run_summary.csv").open(
                newline="", encoding="utf-8",
            ) as handle:
                summary = next(csv.DictReader(handle))
            self.assertEqual(summary["status"], "FAILED")
            self.assertIn(
                "adapter reported a scientific failure", summary["failure_reason"],
            )

    def test_global_preflight_failure_records_every_selected_dataset(self):
        fixture = Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            output = root / "output"
            args = Namespace(
                sample_sheet=str(fixture),
                resource_directory=str(resources),
                output_directory=str(output),
            )
            with patch(
                "postgwas.modules.harmonisation.cli.require_binaries",
                side_effect=RuntimeError("bcftools plugin unavailable"),
            ):
                with self.assertRaisesRegex(RuntimeError, "plugin unavailable"):
                    run_harmonisation(args)

            with (output / "run_metadata" / "harmonisation_run_summary.csv").open(
                newline="", encoding="utf-8",
            ) as handle:
                records = list(csv.DictReader(handle))
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0]["status"], "PREFLIGHT_FAILED")
            self.assertIn("bcftools plugin unavailable", records[0]["failure_reason"])
            report_path = Path(records[0]["screen_report"])
            self.assertEqual(
                report_path,
                (
                    output
                    / records[0]["dataset_id"]
                    / "harmonisation"
                    / (records[0]["dataset_id"] + "_screen_report.txt")
                ).resolve(),
            )
            self.assertIn(
                "bcftools plugin unavailable",
                report_path.read_text(encoding="utf-8"),
            )

    def test_user_interruption_is_recorded_before_propagation(self):
        fixture = Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        rows = load_harmonisation_sample_sheet(fixture)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            output = root / "output"
            config = load_configuration(
                cli_overrides={
                    "resources.root": str(resources),
                    "run.output_directory": str(output),
                }
            )
            with patch(
                "postgwas.modules.harmonisation.cli.run_harmonisation_pipeline",
                side_effect=KeyboardInterrupt,
            ):
                with self.assertRaises(KeyboardInterrupt):
                    _run_validated_rows(Namespace(), config, rows)

            run_log = (
                output / "run_metadata" / config.logging.filename
            ).read_text(encoding="utf-8")
            self.assertIn("Run interrupted by the user", run_log)
            with (output / "run_metadata" / "harmonisation_run_summary.csv").open(
                newline="", encoding="utf-8",
            ) as handle:
                summary = next(csv.DictReader(handle))
            self.assertEqual(summary["status"], "INTERRUPTED")


if __name__ == "__main__":
    unittest.main()
