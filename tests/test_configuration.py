"""Tests for the canonical layered configuration system."""

from argparse import Namespace
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import yaml

from postgwas.config import (
    load_configuration,
    load_module_configuration,
    load_run_configuration_for_module,
    resolved_configuration_values,
    write_resolved_configuration,
)
from postgwas.config.exporter import (
    render_module_configuration,
    render_pipeline_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.errors import ConfigurationError


class ConfigurationTests(unittest.TestCase):
    def test_fine_mapping_example_configurations_are_reloadable(self):
        example_directory = (
            Path(__file__).parents[1] / "examples" / "fine_mapping"
        )
        direct = load_configuration(example_directory / "direct_susie.yaml")
        pipeline = load_configuration(example_directory / "pipeline_susie.yaml")

        self.assertEqual(direct.modules.fine_mapping.engine, "susie")
        self.assertEqual(direct.modules.fine_mapping.locus_type, "range")
        self.assertEqual(direct.modules.fine_mapping.locus_window_kb, 0)
        self.assertEqual(
            pipeline.pipeline.modules,
            ["ld_clumping", "formatting", "fine_mapping"],
        )
        self.assertEqual(
            pipeline.modules.ld_clumping.methods,
            ["standard"],
        )
        self.assertTrue(pipeline.modules.fine_mapping.enabled)
        self.assertEqual(pipeline.modules.fine_mapping.genome_build.value, "GRCh37")

    def test_pipeline_resume_is_enabled_by_default(self):
        configuration = load_configuration()
        self.assertTrue(configuration.run.resume)
        self.assertFalse(configuration.run.overwrite)
        policy = configuration.run.resume_policy
        self.assertEqual(policy.checkpoint_validation, "sha256")
        self.assertEqual(
            policy.checkpoint_directory,
            Path("run_metadata/checkpoints"),
        )
        self.assertEqual(policy.direct_manifest, "{command}_direct.yaml")
        self.assertEqual(
            policy.pipeline_stage_manifest,
            "{stage_number}_{module}.yaml",
        )
        self.assertEqual(policy.audit_log, "checkpoint_events.log")
        self.assertEqual(
            policy.partial_results,
            "resume_validated_stages",
        )
        self.assertEqual(
            policy.changed_parameters,
            "warn_and_restart",
        )
        self.assertEqual(
            policy.changed_inputs,
            "warn_and_restart",
        )
        self.assertEqual(
            policy.unvalidated_outputs,
            "warn_and_restart",
        )
        self.assertTrue(configuration.logging.show_screen)
        self.assertTrue(configuration.logging.show_progress)
        self.assertEqual(configuration.logging.progress_refresh_seconds, 1.0)
        self.assertEqual(configuration.logging.terminal_label_width, 42)
        self.assertEqual(
            configuration.logging.screen_log_file,
            "run_metadata/screen.log",
        )

    def test_resume_policy_paths_and_manifest_patterns_are_schema_validated(self):
        invalid = (
            ("run.resume_policy.checkpoint_directory", "../checkpoints"),
            ("run.resume_policy.direct_manifest", "direct.yaml"),
            ("run.resume_policy.pipeline_stage_manifest", "{module}.yaml"),
            ("run.resume_policy.audit_log", "/tmp/checkpoints.log"),
        )
        for key, value in invalid:
            with self.subTest(key=key), self.assertRaises(ConfigurationError):
                load_configuration(cli_overrides={key: value})

    def test_shared_runtime_cli_values_map_to_canonical_run_configuration(self):
        overrides = explicit_overrides(
            Namespace(resume=False, overwrite=True),
            {},
        )

        self.assertEqual(
            overrides,
            {"run.resume": False, "run.overwrite": True},
        )
        configuration = load_configuration(cli_overrides=overrides)
        self.assertFalse(configuration.run.resume)
        self.assertTrue(configuration.run.overwrite)

    def test_screen_display_can_be_disabled_from_canonical_configuration(self):
        configuration = load_configuration(
            cli_overrides={"logging.show_screen": False},
        )
        self.assertFalse(configuration.logging.show_screen)

    def test_screen_log_file_must_remain_inside_the_output_directory(self):
        with self.assertRaisesRegex(ConfigurationError, "screen_log_file"):
            load_configuration(
                cli_overrides={"logging.screen_log_file": "../screen.log"},
            )

    def test_terminal_progress_can_be_disabled_from_canonical_configuration(self):
        configuration = load_configuration(
            cli_overrides={"logging.show_progress": False},
        )
        self.assertFalse(configuration.logging.show_progress)

    def test_progress_refresh_interval_is_schema_validated(self):
        with self.assertRaisesRegex(ConfigurationError, "progress_refresh_seconds"):
            load_configuration(
                cli_overrides={"logging.progress_refresh_seconds": 0},
            )

    def test_terminal_label_width_is_schema_validated(self):
        with self.assertRaisesRegex(ConfigurationError, "terminal_label_width"):
            load_configuration(
                cli_overrides={"logging.terminal_label_width": 10},
            )

    def test_module_override_without_config_is_qualified_to_module_namespace(self):
        module = load_module_configuration(
            "gcta_gene", cli_overrides={"method": "fastbat_set"},
        )
        self.assertEqual(module.method, "fastbat_set")

    def test_imputation_configuration_exposes_only_the_implemented_engine(self):
        module = load_configuration().modules.imputation

        self.assertEqual(module.engine, "pred_ld")
        self.assertEqual(list(module.engines.model_dump()), ["pred_ld"])
        with self.assertRaisesRegex(ConfigurationError, "engine"):
            load_configuration(
                cli_overrides={"modules.imputation.engine": "raiss"},
            )

    def test_module_loader_preserves_explicit_root_override(self):
        module = load_module_configuration(
            "gcta_gene",
            cli_overrides={
                "method": "mbat_combo",
                "execution.threads": 3,
            },
        )
        self.assertEqual(module.method, "mbat_combo")

    def test_module_only_run_loader_preserves_screen_display_override(self):
        with tempfile.TemporaryDirectory() as directory:
            module_file = Path(directory) / "gcta_gene.yaml"
            module_file.write_text("method: fastbat_set\n", encoding="utf-8")
            configuration = load_run_configuration_for_module(
                "gcta_gene",
                module_file,
                module_overrides={"logging.show_screen": False},
            )

        self.assertFalse(configuration.logging.show_screen)
        self.assertEqual(configuration.modules.gcta_gene.method, "fastbat_set")

    def test_harmonisation_export_styles_share_values_and_analysis_order(self):
        rendered = {
            style: render_module_configuration("harmonisation", style=style)
            for style in ("full", "minimal", "values")
        }
        parsed = {style: yaml.safe_load(text) for style, text in rendered.items()}
        self.assertEqual(parsed["full"], parsed["minimal"])
        self.assertEqual(parsed["minimal"], parsed["values"])
        self.assertFalse(any(
            line.lstrip().startswith("#")
            for line in rendered["values"].splitlines()
        ))
        self.assertGreater(rendered["full"].count("#"), rendered["minimal"].count("#"))
        policy_text = rendered["values"].split("policies:\n", 1)[1]
        self.assertLess(policy_text.index("  input:"), policy_text.index("  chromosome:"))
        self.assertLess(
            policy_text.index("  chromosome:"), policy_text.index("  execution:")
        )
        chromosome_line = next(
            line
            for line in rendered["values"].splitlines()
            if line.strip().startswith("allowed_after_split:")
        )
        self.assertIn("[", chromosome_line)
        self.assertIn("22", chromosome_line)

    def test_exported_harmonisation_values_are_reloadable(self):
        with tempfile.TemporaryDirectory() as directory:
            exported = Path(directory) / "harmonisation.yaml"
            exported.write_text(
                render_module_configuration("harmonisation", style="values"),
                encoding="utf-8",
            )
            config = load_module_configuration("harmonisation", exported)
        self.assertEqual(config.comparison_af.source, "ALFA")
        self.assertEqual(config.default_eaf.source, "ALFA")
        self.assertNotIn("filter", config.policies)
        self.assertEqual(
            load_configuration().modules.qc_summary.rules.maf_min, 0.01
        )

    def test_formatting_export_styles_share_values_and_keep_chromosomes_on_one_line(self):
        rendered = {
            style: render_module_configuration("formatting", style=style)
            for style in ("full", "minimal", "values")
        }
        parsed = {style: yaml.safe_load(text) for style, text in rendered.items()}
        self.assertEqual(parsed["full"], parsed["minimal"])
        self.assertEqual(parsed["minimal"], parsed["values"])
        self.assertNotIn("#", rendered["values"])
        self.assertGreater(rendered["full"].count("#"), rendered["minimal"].count("#"))
        chromosome_line = next(
            line for line in rendered["values"].splitlines()
            if line.startswith("chromosomes:")
        )
        self.assertIn("'22'", chromosome_line)
        self.assertEqual(
            parsed["values"]["ldsc_sample_prevalence"],
            {"aggregation": "median"},
        )

    def test_formatting_export_can_select_only_ldsc_and_remains_reloadable(self):
        rendered = render_module_configuration(
            "formatting", style="values", formats=["ldsc"],
        )
        document = yaml.safe_load(rendered)

        self.assertEqual(document["formats"], ["ldsc"])
        self.assertEqual(list(document["exports"]), ["ldsc"])
        self.assertIn("ldsc_reference", document)
        self.assertIn("ldsc_sample_prevalence", document)
        self.assertIn("study_design", document)
        self.assertEqual(
            document["variant_identifiers"]["target_duplicate_policies"],
            {},
        )
        self.assertNotIn("mixer", document)
        self.assertNotIn("chromosomes", document)
        self.assertNotIn("custom_output", document)
        self.assertNotIn("format_order", document)
        self.assertNotIn("module_formats", document)
        self.assertNotIn("resolved_config", document)

        with tempfile.TemporaryDirectory() as directory:
            exported = Path(directory) / "formatting.yaml"
            exported.write_text(rendered, encoding="utf-8")
            reloaded = load_module_configuration("formatting", exported)

        self.assertEqual(reloaded.formats, ["ldsc"])
        self.assertEqual(
            reloaded.exports["ldsc"].output_file,
            "{dataset_id}_ldsc_input.tsv",
        )

    def test_formatting_export_normalizes_multiple_targets_to_configured_order(self):
        document = yaml.safe_load(render_module_configuration(
            "formatting", style="values", formats=["ldsc", "magma"],
        ))

        self.assertEqual(document["formats"], ["magma", "ldsc"])
        self.assertEqual(list(document["exports"]), ["magma", "ldsc"])
        self.assertEqual(
            document["variant_identifiers"]["target_duplicate_policies"],
            {},
        )

    def test_format_selection_is_rejected_for_other_configuration_modules(self):
        with self.assertRaisesRegex(
            ConfigurationError,
            "only with --module formatting",
        ):
            render_module_configuration("mixer", formats=["ldsc"])

    def test_format_selection_rejects_empty_unknown_and_duplicate_targets(self):
        invalid_selections = (
            ([], "requires at least one"),
            (["unknown"], "Unknown formatter target"),
            (["ldsc", "ldsc"], "must not repeat target"),
        )
        for formats, message in invalid_selections:
            with self.subTest(formats=formats), self.assertRaisesRegex(
                ConfigurationError, message,
            ):
                render_module_configuration("formatting", formats=formats)

    def test_config_cli_exports_only_selected_formatter_target(self):
        from postgwas.__main__ import main

        with tempfile.TemporaryDirectory() as directory:
            output_file = Path(directory) / "formatting.yaml"
            with patch(
                "sys.argv",
                [
                    "postgwas",
                    "config",
                    "export",
                    "--module",
                    "formatting",
                    "--format",
                    "ldsc",
                    "--style",
                    "minimal",
                    "--output",
                    str(output_file),
                ],
            ):
                self.assertEqual(main(), 0)
            document = yaml.safe_load(output_file.read_text(encoding="utf-8"))

        self.assertEqual(document["formats"], ["ldsc"])
        self.assertEqual(list(document["exports"]), ["ldsc"])

    def test_config_cli_rejects_format_selection_for_pipeline_export(self):
        from postgwas.__main__ import main

        error = StringIO()
        with patch(
            "sys.argv",
            [
                "postgwas",
                "config",
                "export",
                "--pipeline",
                "finemap",
                "--format",
                "ldsc",
            ],
        ), redirect_stderr(error):
            self.assertEqual(main(), 2)

        self.assertIn(
            "--format can be used only with --module formatting",
            error.getvalue(),
        )

    def test_pipeline_export_contains_only_selected_modules(self):
        document = yaml.safe_load(
            render_pipeline_configuration(["finemap"], style="values")
        )
        self.assertEqual(
            list(document["modules"]),
            ["ld_annotation", "ld_clumping", "formatting", "fine_mapping"],
        )
        self.assertEqual(
            document["pipeline"]["modules"],
            ["ld_annotation", "ld_clumping", "formatting", "fine_mapping"],
        )
        self.assertTrue(all(module["enabled"] for module in document["modules"].values()))

    def test_standard_finemap_export_omits_unused_ld_annotation(self):
        config_file = (
            Path(__file__).parents[1]
            / "examples"
            / "fine_mapping"
            / "pipeline_susie.yaml"
        )
        document = yaml.safe_load(
            render_pipeline_configuration(
                ["finemap"], config_file=config_file, style="values",
            )
        )

        self.assertEqual(
            document["pipeline"]["modules"],
            ["ld_clumping", "formatting", "fine_mapping"],
        )
        self.assertNotIn("ld_annotation", document["modules"])

    def test_flames_pipeline_export_retains_documented_magmacovar_defaults(self):
        rendered = render_pipeline_configuration(["flames"], style="values")
        document = yaml.safe_load(rendered)
        magmacovar = document["modules"]["magmacovar"]

        self.assertEqual(magmacovar["model"], [])
        self.assertEqual(magmacovar["direction"], "two-sided")
        with tempfile.TemporaryDirectory() as directory:
            exported = Path(directory) / "flames_pipeline.yaml"
            exported.write_text(rendered, encoding="utf-8")
            reloaded = load_configuration(exported)

        self.assertEqual(
            reloaded.modules.magmacovar.model,
            [],
        )
        self.assertEqual(reloaded.modules.magmacovar.direction, "two-sided")

    def test_mixer_pipeline_export_selects_its_formatter_contract(self):
        document = yaml.safe_load(
            render_pipeline_configuration(["mixer"], style="values")
        )
        self.assertEqual(list(document["modules"]), ["formatting", "mixer"])
        self.assertEqual(document["modules"]["formatting"]["formats"], ["mixer"])

    def test_gcta_gene_pipeline_export_selects_its_formatter_contract(self):
        document = yaml.safe_load(
            render_pipeline_configuration(["gcta_gene"], style="values")
        )
        self.assertEqual(list(document["modules"]), ["formatting", "gcta_gene"])
        self.assertEqual(
            document["modules"]["formatting"]["formats"], ["gcta_gene"]
        )

    def test_gcta_cojo_pipeline_export_reuses_shared_gcta_formatter_contract(self):
        document = yaml.safe_load(
            render_pipeline_configuration(["gcta_cojo"], style="values")
        )
        self.assertEqual(list(document["modules"]), ["formatting", "gcta_cojo"])
        self.assertEqual(
            document["modules"]["formatting"]["formats"], ["gcta_gene"]
        )

    def test_mixer_export_styles_have_identical_ordered_values(self):
        rendered = {
            style: render_module_configuration("mixer", style=style)
            for style in ("full", "minimal", "values")
        }
        parsed = {style: yaml.safe_load(text) for style, text in rendered.items()}
        self.assertEqual(parsed["full"], parsed["minimal"])
        self.assertEqual(parsed["minimal"], parsed["values"])
        self.assertNotIn("#", rendered["values"])
        self.assertGreater(rendered["full"].count("#"), rendered["minimal"].count("#"))
        self.assertLess(
            rendered["values"].index("bim_file_pattern:"),
            rendered["values"].index("ld_file_pattern:"),
        )

    def test_magma_export_styles_are_annotated_reloadable_and_ordered(self):
        rendered = {
            style: render_module_configuration("magma", style=style)
            for style in ("full", "minimal", "values")
        }
        parsed = {style: yaml.safe_load(text) for style, text in rendered.items()}
        self.assertEqual(parsed["full"], parsed["minimal"])
        self.assertEqual(parsed["minimal"], parsed["values"])
        self.assertFalse(
            any(
                line.lstrip().startswith("#")
                for line in rendered["values"].splitlines()
            )
        )
        self.assertGreater(rendered["full"].count("#"), rendered["minimal"].count("#"))
        self.assertLess(
            rendered["values"].index("input:"),
            rendered["values"].index("snp_harmonisation:"),
        )

    def test_mixer_chromosome_range_is_configured_without_a_code_maximum(self):
        with tempfile.TemporaryDirectory() as directory:
            module_file = Path(directory) / "mixer.yaml"
            module_file.write_text("chromosomes: ['23']\n", encoding="utf-8")
            mixer = load_module_configuration("mixer", module_file)

        self.assertEqual(mixer.chromosomes, ["23"])

    def test_global_config_shorthand_writes_yaml_only(self):
        from postgwas.__main__ import main

        output = StringIO()
        with patch(
            "sys.argv",
            [
                "postgwas",
                "--config",
                "--module",
                "harmonisation",
                "--style",
                "values",
            ],
        ), redirect_stdout(output):
            self.assertEqual(main(), 0)
        rendered = yaml.safe_load(output.getvalue())
        self.assertEqual(rendered["comparison_af"]["source"], "ALFA")
        self.assertEqual(rendered["default_eaf"]["source"], "ALFA")

    def test_pipeline_shorthand_expands_finemap_dependencies(self):
        from postgwas.__main__ import main

        output = StringIO()
        with patch(
            "sys.argv",
            ["postgwas", "--config", "--pipeline", "finemap", "--style", "values"],
        ), redirect_stdout(output):
            self.assertEqual(main(), 0)
        document = yaml.safe_load(output.getvalue())
        self.assertEqual(
            document["pipeline"]["modules"],
            ["ld_annotation", "ld_clumping", "formatting", "fine_mapping"],
        )

    def test_config_cli_does_not_advertise_profiles(self):
        from postgwas.config.cli import build_export_parser, build_parser

        self.assertNotIn("--profile", build_parser().format_help())
        help_text = build_export_parser().format_help()
        self.assertIn("--module MODULE", help_text)
        self.assertIn("--format FORMAT [FORMAT ...]", help_text)
        self.assertIn("Available modules:\n", help_text)
        self.assertIn("--run-config PATH", help_text)
        self.assertNotIn("postgwas config export --module", help_text)

    def test_packaged_defaults_validate(self):
        config = load_configuration()
        self.assertEqual(config.config_version, 1)
        self.assertEqual(config.modules.filtering.maf_min, 0.01)
        self.assertEqual(config.modules.fine_mapping.engine, "susie")
        self.assertEqual(
            config.modules.ld_annotation.bed_filename_template,
            "{genome_build}_{population}_ldetect.bed.gz",
        )
        self.assertEqual(
            config.modules.ld_annotation.info_field_template,
            "{population}_LDblock",
        )
        self.assertEqual(
            config.modules.ld_annotation.output_filename_template,
            "{dataset_id}_ldblock.vcf.gz",
        )
        self.assertEqual(
            config.modules.ld_annotation.summary_filename_template,
            "{dataset_id}_ldblock_summary.csv",
        )
        self.assertEqual(
            config.modules.ld_annotation.html_report_filename_template,
            "{dataset_id}_ldblock_report.html",
        )
        self.assertEqual(
            config.modules.ld_annotation.canonical_log_filename_template,
            "{dataset_id}_ldblock.log",
        )
        self.assertIsNone(config.modules.ld_annotation.inputs.vcf)
        self.assertIsNone(config.modules.ld_annotation.inputs.ld_region_dir)
        self.assertIsNone(config.modules.ld_annotation.inputs.dataset_id)
        self.assertIsNone(config.modules.ld_annotation.output_directory)
        self.assertNotIn(
            "ld_blocks",
            type(config.resources.populations["EUR"]).model_fields,
        )
        with self.assertRaises(ConfigurationError):
            load_configuration(
                cli_overrides={
                    "resources.populations.EUR.ld_blocks": "blocks.bed.gz",
                }
            )

    def test_ld_annotation_bed_filename_template_is_schema_validated(self):
        invalid = (
            "ldetect.bed.gz",
            "{genome_build}_{population}_{unexpected}.bed.gz",
            "../{genome_build}_{population}_ldetect.bed.gz",
            "{genome_build}_{population}_ldetect.txt.gz",
        )
        for template in invalid:
            with self.subTest(template=template), self.assertRaises(
                ConfigurationError
            ):
                load_configuration(
                    cli_overrides={
                        "modules.ld_annotation.bed_filename_template": template,
                    },
                )

    def test_ld_annotation_output_and_info_templates_are_schema_validated(self):
        invalid = {
            "modules.ld_annotation.info_field_template": (
                "LDblock",
                "{population}_{unexpected}",
                "{population} LDblock",
            ),
            "modules.ld_annotation.info_description_template": (
                "LD block for {population}",
                'LD "block" for {population} {genome_build}',
                "LD block for {population} {genome_build} {unexpected}",
            ),
            "modules.ld_annotation.output_filename_template": (
                "ldblock.vcf.gz",
                "../{dataset_id}_ldblock.vcf.gz",
                "{dataset_id}_ldblock.vcf",
            ),
            "modules.ld_annotation.summary_filename_template": (
                "ldblock_summary.csv",
                "../{dataset_id}_ldblock_summary.csv",
                "{dataset_id}_ldblock_summary.tsv",
            ),
            "modules.ld_annotation.html_report_filename_template": (
                "ldblock_report.html",
                "../{dataset_id}_ldblock_report.html",
                "{dataset_id}_ldblock_report.htm",
            ),
            "modules.ld_annotation.canonical_log_filename_template": (
                "ldblock.log",
                "../{dataset_id}_ldblock.log",
                "{dataset_id}_ldblock.txt",
            ),
        }
        for key, values in invalid.items():
            for value in values:
                with self.subTest(key=key, value=value), self.assertRaises(
                    ConfigurationError
                ):
                    load_configuration(cli_overrides={key: value})

    def test_annotation_and_clumping_info_fields_match_in_one_pipeline(self):
        with self.assertRaisesRegex(
            ConfigurationError,
            "info_field_template.*vcf_fields.ld_block",
        ):
            load_configuration(
                cli_overrides={
                    "pipeline.modules": ["ld_annotation", "ld_clumping"],
                    "modules.ld_annotation.info_field_template": (
                        "{population}_CUSTOM_BLOCK"
                    ),
                },
            )

    def test_precedence_is_cli_user_profile_defaults(self):
        with tempfile.TemporaryDirectory() as directory:
            config_file = Path(directory) / "run.yaml"
            config_file.write_text(
                "modules:\n  filtering:\n    maf_min: 0.02\n",
                encoding="utf-8",
            )
            config = load_configuration(
                config_file,
                profile="minimal",
                cli_overrides={"modules.filtering.maf_min": 0.03},
            )
        self.assertTrue(config.modules.filtering.enabled)
        self.assertEqual(config.modules.filtering.maf_min, 0.03)

    def test_unknown_keys_are_rejected_with_complete_path(self):
        with tempfile.TemporaryDirectory() as directory:
            config_file = Path(directory) / "invalid.yaml"
            config_file.write_text(
                "modules:\n  filtering:\n    maf_mim: 0.01\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                ConfigurationError, r"modules\.filtering\.maf_mim"
            ):
                load_configuration(config_file)

    def test_scientific_ranges_are_validated(self):
        with tempfile.TemporaryDirectory() as directory:
            config_file = Path(directory) / "invalid.yaml"
            config_file.write_text(
                "modules:\n  ld_clumping:\n    clump_r2: 1.5\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ConfigurationError, "less than or equal to 1"):
                load_configuration(config_file)

    def test_dataset_id_cannot_escape_configured_output_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            config_file = Path(directory) / "invalid.yaml"
            config_file.write_text(
                "run:\n  dataset_id: ../outside\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ConfigurationError, "run.dataset_id"):
                load_configuration(config_file)

    def test_module_only_file_uses_the_same_model_and_defaults(self):
        with tempfile.TemporaryDirectory() as directory:
            config_file = Path(directory) / "filtering.yaml"
            config_file.write_text("enabled: true\ninfo_min: 0.9\n", encoding="utf-8")
            config = load_module_configuration("sumstat_filter", config_file)
        self.assertTrue(config.enabled)
        self.assertEqual(config.info_min, 0.9)
        self.assertEqual(config.maf_min, 0.01)

    def test_relative_includes_are_layered(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "base.yaml").write_text(
                "execution:\n  threads: 4\n", encoding="utf-8"
            )
            (root / "run.yaml").write_text(
                "include: base.yaml\nexecution:\n  memory_gb: 16\n",
                encoding="utf-8",
            )
            config = load_configuration(root / "run.yaml")
        self.assertEqual(config.execution.threads, 4)
        self.assertEqual(config.execution.memory_gb, 16)

    def test_include_cycles_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "a.yaml").write_text("include: b.yaml\n", encoding="utf-8")
            (root / "b.yaml").write_text("include: a.yaml\n", encoding="utf-8")
            with self.assertRaisesRegex(ConfigurationError, "include cycle"):
                load_configuration(root / "a.yaml")

    def test_resolved_configuration_is_written(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "metadata" / "resolved_config.yaml"
            write_resolved_configuration(load_configuration(), output)
            written = yaml.safe_load(output.read_text(encoding="utf-8"))
        self.assertEqual(written["modules"]["filtering"]["maf_min"], 0.01)
        self.assertEqual(written["config_version"], 1)

    def test_scoped_resolved_configuration_contains_only_used_sections(self):
        config = load_configuration(
            cli_overrides={
                "resources.executables.magma": "/opt/magma",
                "modules.magma.input.gene_location_file": "/data/genes.loc",
            }
        )
        written = resolved_configuration_values(
            config,
            modules=("magma",),
            resource_paths=("executables.magma",),
        )

        self.assertEqual(
            list(written),
            ["config_version", "run", "execution", "logging", "resources", "modules"],
        )
        self.assertEqual(list(written["modules"]), ["magma"])
        self.assertEqual(
            written["modules"]["magma"],
            config.modules.magma.model_dump(mode="json"),
        )
        self.assertEqual(
            written["resources"], {"executables": {"magma": "/opt/magma"}}
        )
        self.assertNotIn("pipeline", written)

    def test_scoped_resolved_configuration_is_reloadable(self):
        config = load_configuration(
            cli_overrides={
                "run.dataset_id": "SCOPED",
                "execution.threads": 3,
                "resources.executables.magma": "/opt/magma",
                "modules.magma.gene_window_upstream_kb": 42,
            }
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "resolved_config.yaml"
            write_resolved_configuration(
                config,
                output,
                modules=("magma",),
                resource_paths=("executables.magma",),
            )
            reloaded = load_configuration(output)

        self.assertEqual(reloaded.run.dataset_id, "SCOPED")
        self.assertEqual(reloaded.execution.threads, 3)
        self.assertEqual(reloaded.resources.executables.magma, "/opt/magma")
        self.assertEqual(reloaded.modules.magma.gene_window_upstream_kb, 42)
        self.assertEqual(reloaded.modules.magma.enabled, config.modules.magma.enabled)

    def test_scoped_resolved_configuration_can_select_module_fields(self):
        config = load_configuration(
            cli_overrides={"modules.formatting.formats": ["magma"]},
        )
        written = resolved_configuration_values(
            config,
            modules=("formatting",),
            module_paths={
                "formatting": ("formats", "exports.magma"),
            },
        )

        formatting = written["modules"]["formatting"]
        self.assertEqual(list(formatting), ["formats", "exports"])
        self.assertEqual(list(formatting["exports"]), ["magma"])
        self.assertNotIn("mixer", formatting)


if __name__ == "__main__":
    unittest.main()
