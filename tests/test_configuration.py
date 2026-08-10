"""Tests for the canonical layered configuration system."""

import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import yaml

from postgwas.config import (
    load_configuration,
    load_module_configuration,
    resolved_configuration_values,
    write_resolved_configuration,
)
from postgwas.config.exporter import (
    render_module_configuration,
    render_pipeline_configuration,
)
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
            ["ld_annotation", "formatting", "ld_clumping", "fine_mapping"],
        )
        self.assertTrue(pipeline.modules.fine_mapping.enabled)
        self.assertEqual(pipeline.modules.fine_mapping.genome_build.value, "GRCh37")

    def test_pipeline_resume_is_enabled_by_default(self):
        configuration = load_configuration()
        self.assertTrue(configuration.run.resume)
        self.assertFalse(configuration.run.overwrite)
        self.assertTrue(configuration.logging.show_progress)
        self.assertEqual(configuration.logging.terminal_label_width, 42)

    def test_terminal_progress_can_be_disabled_from_canonical_configuration(self):
        configuration = load_configuration(
            cli_overrides={"logging.show_progress": False},
        )
        self.assertFalse(configuration.logging.show_progress)

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

    def test_module_loader_preserves_explicit_root_override(self):
        module = load_module_configuration(
            "gcta_gene",
            cli_overrides={
                "method": "mbat_combo",
                "execution.threads": 3,
            },
        )
        self.assertEqual(module.method, "mbat_combo")

    def test_harmonisation_export_styles_share_values_and_analysis_order(self):
        rendered = {
            style: render_module_configuration("harmonisation", style=style)
            for style in ("full", "minimal", "values")
        }
        parsed = {style: yaml.safe_load(text) for style, text in rendered.items()}
        self.assertEqual(parsed["full"], parsed["minimal"])
        self.assertEqual(parsed["minimal"], parsed["values"])
        self.assertNotIn("#", rendered["values"])
        self.assertGreater(rendered["full"].count("#"), rendered["minimal"].count("#"))
        policy_text = rendered["values"].split("policies:\n", 1)[1]
        self.assertLess(policy_text.index("  input:"), policy_text.index("  chromosome:"))
        self.assertLess(policy_text.index("  chromosome:"), policy_text.index("  filter:"))
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
        self.assertEqual(config.comparison_af.source, "1000G")
        self.assertEqual(config.policies["filter"]["maf_cutoff"], 0.01)

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

    def test_pipeline_export_contains_only_selected_modules(self):
        document = yaml.safe_load(
            render_pipeline_configuration(["finemap"], style="values")
        )
        self.assertEqual(
            list(document["modules"]),
            ["ld_annotation", "formatting", "ld_clumping", "fine_mapping"],
        )
        self.assertEqual(
            document["pipeline"]["modules"],
            ["ld_annotation", "formatting", "ld_clumping", "fine_mapping"],
        )
        self.assertTrue(all(module["enabled"] for module in document["modules"].values()))

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
        self.assertEqual(yaml.safe_load(output.getvalue())["comparison_af"]["source"], "1000G")

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
            ["ld_annotation", "formatting", "ld_clumping", "fine_mapping"],
        )

    def test_config_cli_does_not_advertise_profiles(self):
        from postgwas.config.cli import build_export_parser, build_parser

        self.assertNotIn("--profile", build_parser().format_help())
        help_text = build_export_parser().format_help()
        self.assertIn("--module MODULE", help_text)
        self.assertIn("Available modules:\n", help_text)
        self.assertIn("--run-config PATH", help_text)
        self.assertNotIn("postgwas config export --module", help_text)

    def test_packaged_defaults_validate(self):
        config = load_configuration()
        self.assertEqual(config.config_version, 1)
        self.assertEqual(config.modules.filtering.maf_min, 0.01)
        self.assertEqual(config.modules.fine_mapping.engine, "susie")

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
