"""CLI/config boundary tests; no scientific work is started."""

import argparse
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

import yaml

from postgwas.core.errors import ConfigurationError
from postgwas.cli.compute import get_compute_parser, memory_limit, resolve_compute_args
from postgwas.modules.harmonisation.cli import (
    _run_validated_rows,
    get_harmonisation_parser,
    run_harmonisation,
)
from postgwas.modules.harmonisation.service import PipelineError
from postgwas.config import load_configuration
from postgwas.modules.harmonisation.sample_sheet import load_harmonisation_sample_sheet


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
        }
        actions = {action.dest: action for action in parser._actions}
        for destination in configurable:
            self.assertEqual(actions[destination].default, argparse.SUPPRESS)

    def test_help_displays_internal_defaults(self):
        help_text = get_harmonisation_parser().format_help()
        normalized_help = " ".join(help_text.split())
        self.assertIn("Default: 1000G", help_text)
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
                return_value={"status": "OK"},
            ) as engine:
                result = _run_validated_rows(Namespace(), config, rows)
            metadata = output / rows[0].dataset_id / "run_metadata"
            self.assertTrue((metadata / "resolved_config.yaml").is_file())
            resolved = yaml.safe_load(
                (metadata / "resolved_config.yaml").read_text(encoding="utf-8")
            )
            self.assertEqual(list(resolved["modules"]), ["harmonisation"])
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

    def test_integrated_validation_uses_same_build_merged_vcf(self):
        fixture = Path(__file__).parent / "data" / "harmonisation" / "manifest_v2.csv"
        rows = load_harmonisation_sample_sheet(fixture)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources"
            resources.mkdir()
            output = root / "output"
            manifest = root / "run_manifest.json"
            manifest.write_text(
                '{"dataset":{"genome_build":{"inferred_build":"GRCh37"}}}',
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


if __name__ == "__main__":
    unittest.main()
