"""Configuration behavior specific to harmonisation."""

import io
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml
import polars as pl

from postgwas.config import load_configuration
from postgwas.config.exporter import render_module_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.modules.harmonisation.cli import _engine_defaults
from postgwas.modules.harmonisation.policies import (
    PolicyError,
    default_policies,
    registry_path,
)
from postgwas.modules.harmonisation.p_values import (
    CLIPPED_HIGH_COLUMN,
    harmonise_p_values,
)
from postgwas.modules.harmonisation.service import _resolve_policies, build_resource_map
from postgwas.modules.harmonisation.genome_build import infer_genome_build
from postgwas.modules.harmonisation.study_properties import (
    StudyPropertyError,
    resolve_study_properties,
)
from postgwas.modules.harmonisation.shared.runtime import emit_high_visibility_warning
from postgwas.modules.harmonisation.summary_statistics_io import (
    load_summary_statistics_table,
)
from postgwas.modules.harmonisation.input_validation import validate_header


class HarmonisationConfigTests(unittest.TestCase):
    def test_quantitative_trait_help_matches_enforced_runtime_behavior(self):
        registry = yaml.safe_load(Path(registry_path()).read_text(encoding="utf-8"))
        help_text = registry["policies"]["sample_size"]["trait_type"]["help"]

        self.assertIn("any case-count field stops the run", help_text)
        self.assertIn("represents total sample size", help_text)
        self.assertNotIn("only prints a warning", help_text)
        self.assertNotIn("not a guard", help_text)

    def test_dataset_trait_type_reaches_sample_size_policy(self):
        defaults = _engine_defaults(load_configuration())

        self.assertEqual(
            _resolve_policies(defaults, {"trait_type": "case_control"}).get(
                "sample_size.trait_type"
            ),
            "binary",
        )
        self.assertEqual(
            _resolve_policies(defaults, {"trait_type": "quantitative"}).get(
                "sample_size.trait_type"
            ),
            "quantitative",
        )
        self.assertEqual(
            _resolve_policies(defaults, {"trait_type": "auto"}).get(
                "sample_size.trait_type"
            ),
            "auto",
        )

    def test_declared_trait_type_overrides_generic_sample_size_policy(self):
        defaults = _engine_defaults(load_configuration())
        policies = _resolve_policies(
            defaults,
            {
                "trait_type": "case_control",
                "policies": {"sample_size": {"trait_type": "quantitative"}},
            },
        )

        self.assertEqual(policies.get("sample_size.trait_type"), "binary")

    def test_sample_sheet_scientific_declarations_override_generic_policies(self):
        defaults = _engine_defaults(load_configuration())
        policies = _resolve_policies(
            defaults,
            {
                "declared_effect_type": "beta",
                "declared_pvalue_type": "negln",
                "delimiter": "comma",
                "policies": {
                    "effect": {"type": "odds_ratio"},
                    "pvalue": {"type": "raw"},
                    "input": {"delimiter": "tab"},
                },
            },
        )

        self.assertEqual(policies.get("effect.type"), "beta")
        self.assertEqual(policies.get("pvalue.type"), "negln")
        self.assertEqual(policies.get("input.delimiter"), "comma")

    def test_auto_sample_sheet_values_preserve_yaml_policies(self):
        policies = _resolve_policies(
            {
                "policies": {
                    "effect": {"type": "beta"},
                    "pvalue": {"type": "raw"},
                    "input": {"delimiter": "tab"},
                }
            },
            {
                "declared_effect_type": "auto",
                "declared_pvalue_type": "auto",
                "delimiter": "auto",
            },
        )

        self.assertEqual(policies.get("effect.type"), "beta")
        self.assertEqual(policies.get("pvalue.type"), "raw")
        self.assertEqual(policies.get("input.delimiter"), "tab")

    def test_declared_effect_is_cross_checked_but_remains_authoritative(self):
        decisions = resolve_study_properties(
            pl.DataFrame({
                "BETA": [0.90, 0.98, 1.02, 1.10],
                "P": [0.5, 0.4, 0.3, 0.2],
                "EAF": [0.6, 0.7, 0.8, 0.9],
            }),
            {
                "beta_or_col": "BETA",
                "pval_col": "P",
                "eaf_col": "EAF",
                "declared_effect_type": "beta",
                "declared_pvalue_type": "raw",
            },
            policies=default_policies().with_overrides({"effect.type": "beta"}),
        )

        self.assertEqual(decisions["effect_type"], "beta")
        self.assertEqual(decisions["effect_type_source"], "sample_sheet")
        self.assertEqual(decisions["effect_type_detected"], "odds_ratio")
        self.assertFalse(decisions["effect_type_matches_declaration"])
        self.assertEqual(len(decisions["declaration_mismatches"]), 1)

    def test_declared_type_mismatch_can_fail_before_chromosome_work(self):
        policies = default_policies().with_overrides({
            "effect.type": "beta",
            "validation.declaration_mismatch_action": "fail",
        })
        with self.assertRaisesRegex(StudyPropertyError, "DECLARATION MISMATCH"):
            resolve_study_properties(
                pl.DataFrame({
                    "BETA": [0.90, 1.00, 1.10],
                    "P": [0.5, 0.4, 0.3],
                    "EAF": [0.6, 0.7, 0.8],
                }),
                {
                    "beta_or_col": "BETA",
                    "pval_col": "P",
                    "eaf_col": "EAF",
                    "declared_effect_type": "beta",
                    "declared_pvalue_type": "raw",
                },
                policies=policies,
            )

    def test_declared_pvalue_is_also_cross_checked(self):
        decisions = resolve_study_properties(
            pl.DataFrame({
                "BETA": [-0.2, -0.1, 0.1, 0.2, -0.3, 0.3, -0.4, 0.4, -0.5, 0.5],
                "P": [0.1, 0.2, 0.3, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 2.0],
                "EAF": [0.6] * 10,
            }),
            {
                "beta_or_col": "BETA",
                "pval_col": "P",
                "eaf_col": "EAF",
                "declared_effect_type": "beta",
                "declared_pvalue_type": "raw",
            },
            policies=default_policies().with_overrides({
                "effect.type": "beta",
                "pvalue.type": "raw",
            }),
        )

        self.assertEqual(decisions["pvalue_type"], "raw")
        self.assertEqual(decisions["pvalue_type_detected"], "neglog10")
        self.assertFalse(decisions["pvalue_type_matches_declaration"])

    def test_declared_delimiter_is_independently_detected(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            source = Path(directory) / "study.csv"
            source.write_text("CHR,POS,P\n1,10,0.5\n", encoding="utf-8")
            output = io.StringIO()
            original = sys.stdout
            try:
                sys.stdout = output
                valid, _problems = validate_header(
                    {"delimiter": "tab", "chr_col": "CHR", "pos_col": "POS", "pval_col": "P"},
                    str(source),
                    policies=default_policies().with_overrides({"input.delimiter": "tab"}),
                )
            finally:
                sys.stdout = original

        self.assertFalse(valid)
        self.assertIn("DECLARATION MISMATCH", output.getvalue())
        self.assertIn("automatic detection found ','", output.getvalue())

    def test_declaration_mismatch_default_and_terminal_colour(self):
        self.assertEqual(
            default_policies().get("validation.declaration_mismatch_action"),
            "warn",
        )

        class InteractiveBuffer(io.StringIO):
            def isatty(self):
                return True

        stream = InteractiveBuffer()
        original = sys.stdout
        try:
            sys.stdout = stream
            emit_high_visibility_warning(None, "declaration differs")
        finally:
            sys.stdout = original
        self.assertIn("\033[1;91m", stream.getvalue())
        self.assertIn("declaration differs", stream.getvalue())

    def test_raw_pvalue_equal_to_one_is_valid_and_unchanged(self):
        frame, qc, _ = harmonise_p_values(
            chromosome="1",
            df=pl.DataFrame({"P": [1.0, 0.5]}),
            sample_column_dict={"pval_col": "P"},
            decision="raw",
            policies=default_policies(),
        )
        self.assertEqual(default_policies().get("pvalue.clip_high"), 1.0)
        self.assertEqual(frame["LP"].to_list(), [1.0, 0.5])
        self.assertEqual(frame[CLIPPED_HIGH_COLUMN].to_list(), [False, False])
        self.assertEqual(qc["variants_with_pvalues_clipped_high"], 0)

    def test_converter_summary_location_is_owned_by_canonical_output_layout(self):
        config = load_configuration().modules.harmonisation
        self.assertEqual(
            config.output_layout.root["dataset_directory"],
            "{dataset_id}/harmonisation",
        )
        self.assertNotIn("analysis_directory", config.output_layout.root)
        self.assertEqual(
            config.output_layout.root["gwas2vcf_summary"],
            "qc_summary/{dataset_id}_gwas2vcf_summary.tsv",
        )
        self.assertEqual(
            config.output_layout.root["chromosome_source_snapshot"],
            "rejected/{dataset_id}_chr{chromosome}_source.parquet",
        )
        self.assertEqual(
            config.output_layout.root["dataset_reject"],
            "rejected/{dataset_id}_rejected_variants.tsv.gz",
        )
        self.assertEqual(
            config.output_layout.root["duplicates"],
            "qc_summary/{dataset_id}_duplicates.tsv",
        )
        self.assertEqual(
            config.output_layout.root["qc_summary"],
            "qc_summary/{dataset_id}_QC_summary.txt",
        )

    def test_export_has_no_duplicate_module_and_policy_settings(self):
        exported = yaml.safe_load(
            render_module_configuration("harmonisation", style="values")
        )
        duplicate_module_keys = {
            "genome_build",
            "effect_type",
            "p_value_type",
            "delimiter",
            "frequency",
            "output_schema",
            "drop_mitochondrial",
            "minimum_position",
            "info_min",
            "maf_min",
            "include_indels",
            "fail_on_chromosome_error",
        }
        self.assertFalse(duplicate_module_keys & set(exported))
        self.assertFalse({
            "schema_inference_rows", "null_values", "p_value_tail",
            "p_value_floor", "p_value_ceiling",
        } & set(exported["concordance_validation"]))
        self.assertIn("standard_error", exported["concordance_validation"])
        self.assertEqual(
            exported["concordance_validation"]["palindromic_action"],
            "compare_resolved",
        )
        self.assertNotIn(
            "indel_representation_action", exported["concordance_validation"],
        )
        self.assertEqual(
            set(exported["concordance_validation"]["failure"]),
            {
                "maximum_value_mismatch_fraction",
                "maximum_vcf_duplicate_records",
                "maximum_invalid_vcf_records",
            },
        )
        self.assertEqual(exported["policies"]["build"]["mode"], "auto")
        self.assertEqual(
            exported["policies"]["strand"][
                "palindromic_af_discordance_action"
            ],
            "reject",
        )
        self.assertEqual(
            exported["policies"]["external_reference"],
            {
                "exact_duplicate_action": "keep_one",
                "non_identical_duplicate_action": "discard_all",
            },
        )
        self.assertNotIn("deduplicate_reference", exported["policies"]["eaf"])
        self.assertNotIn("deduplicate_reference", exported["policies"]["info"])
        self.assertEqual(
            exported["policies"]["eaf"]["palindromic_handling"],
            "ignore",
        )
        self.assertNotIn(
            "palindromic_resolve_tolerance", exported["policies"]["eaf"]
        )
        self.assertEqual(
            list(exported["policies"]["build"]),
            [
                "mode", "confidence_ratio", "min_match_count",
                "min_match_fraction", "deduplicate_reference",
            ],
        )
        self.assertEqual(exported["policies"]["filter"]["maf_cutoff"], 0.01)
        self.assertEqual(
            list(exported["policies"]["filter"]),
            [
                "lp_cutoff", "lp_missing", "maf_cutoff", "af_missing",
                "info_cutoff", "info_max", "info_missing", "af_diff_cutoff",
                "include_indels", "exclude_palindromic",
                "palindromic_af_lower", "palindromic_af_upper", "remove_mhc",
                "mhc_chrom", "mhc_start", "mhc_end",
            ],
        )
        self.assertTrue(
            exported["policies"]["execution"]["fail_dataset_on_chr_error"]
        )

    def test_module_defaults_and_policy_registry_share_one_canonical_file(self):
        project = Path(__file__).parents[1]
        canonical = (
            project
            / "src"
            / "postgwas"
            / "config"
            / "defaults"
            / "modules"
            / "harmonisation.yaml"
        )
        old_location = (
            project
            / "src"
            / "postgwas"
            / "modules"
            / "harmonisation"
            / "policy_defaults.yaml"
        )
        document = yaml.safe_load(canonical.read_text(encoding="utf-8"))

        self.assertEqual(Path(registry_path()).resolve(), canonical.resolve())
        self.assertIn("module", document)
        self.assertIn("policies", document)
        self.assertFalse(old_location.exists())
        self.assertEqual(default_policies().get("filter.maf_cutoff"), 0.01)
        self.assertEqual(default_policies().get("sample_size.cases_only"), "fail")
        self.assertEqual(default_policies().get("sample_size.min_value"), 1)

    def test_mandatory_fields_are_the_only_input_missingness_policy(self):
        policies = default_policies()
        self.assertIn("columns.mandatory", policies)
        self.assertNotIn("columns.missing_rule", policies)
        self.assertNotIn("columns.min_non_null_fraction", policies)

    def test_z_score_division_floor_is_configuration_driven(self):
        self.assertEqual(
            default_policies().get("validation.se_division_floor"),
            1.0e-12,
        )

    def test_hardware_defaults_are_resolved_from_internal_rules(self):
        with patch("postgwas.config.loader.detected_logical_threads", return_value=16), patch(
            "postgwas.config.loader.detected_physical_memory_gb", return_value=64
        ):
            config = load_configuration()
        self.assertEqual(config.execution.threads, 15)
        self.assertEqual(config.execution.memory_gb, 57)

    def test_explicit_execution_values_override_automatic_values(self):
        config = load_configuration(
            cli_overrides={"execution.threads": 3, "execution.memory_gb": 12}
        )
        self.assertEqual(config.execution.threads, 3)
        self.assertEqual(config.execution.memory_gb, 12)

    def test_no_default_eaf_or_info_source_reaches_engine(self):
        defaults = _engine_defaults(load_configuration())
        for key in (
            "default_eaf",
            "default_info",
            "default_eaf_column",
            "default_info_column",
        ):
            self.assertNotIn(key, defaults)
        fixed = load_configuration().modules.harmonisation.fixed_info
        self.assertIsNone(fixed.value)
        self.assertEqual(fixed.column, "__postgwas_fixed_info")

    def test_comparison_af_comes_from_configuration(self):
        config = load_configuration(
            cli_overrides={
                "modules.harmonisation.comparison_af.source": "1000G",
                "modules.harmonisation.comparison_af.column": "AFR",
            }
        )
        defaults = _engine_defaults(config)
        self.assertEqual(defaults["default_comparison_af_file"], "1000G")
        self.assertEqual(defaults["default_comparison_af_column"], "AFR")

    def test_population_frequency_qc_uses_exactly_five_configured_info_fields(self):
        config = load_configuration().modules.harmonisation.population_frequency_qc
        defaults = _engine_defaults(load_configuration())

        self.assertEqual(config.study_field, "%INFO/AF")
        self.assertEqual(
            config.population_fields,
            {
                "AFR": "%INFO/AFR",
                "EAS": "%INFO/EAS",
                "EUR": "%INFO/EUR",
                "SAS": "%INFO/SAS",
            },
        )
        self.assertEqual(
            defaults["population_frequency_qc"], config.model_dump()
        )

    def test_default_eaf_table_contract_comes_from_configuration(self):
        config = load_configuration().modules.harmonisation
        expected_sources = ["ALFA", "wgs_ukb", "panukb", "1000G", "fingen"]

        self.assertEqual(config.comparison_af.source, "ALFA")
        self.assertEqual(config.default_eaf.source, "ALFA")
        self.assertEqual(config.default_eaf.column, "EUR")
        self.assertEqual(config.default_eaf.available_sources, expected_sources)
        self.assertEqual(
            set(config.default_eaf.resource_examples), set(expected_sources)
        )
        self.assertEqual(config.default_eaf_mapping.chromosome, "CHROM")
        self.assertEqual(config.default_eaf_mapping.position, "POS")
        self.assertEqual(config.default_eaf_mapping.effect_allele, "ALT")
        self.assertEqual(config.default_eaf_mapping.other_allele, "REF")
        self.assertEqual(config.default_eaf_mapping.delimiter, "tab")
        self.assertEqual(
            _engine_defaults(load_configuration())["default_eaf_colmap"],
            {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"},
        )
        self.assertEqual(
            _engine_defaults(load_configuration())["default_eaf_reference_source"],
            "ALFA",
        )
        self.assertEqual(
            _engine_defaults(load_configuration())["default_eaf_reference_column"],
            "EUR",
        )
        self.assertEqual(
            config.resource_layout.default_eaf,
            "{build}/default_af/tab_files/"
            "{build}_{source}_freq_chr{chromosome}.tsv.gz",
        )

    def test_default_eaf_source_must_be_declared_in_canonical_registry(self):
        for source in ("ALFA", "wgs_ukb", "panukb", "1000G", "fingen"):
            config = load_configuration(cli_overrides={
                "modules.harmonisation.default_eaf.source": source,
            })
            self.assertEqual(config.modules.harmonisation.default_eaf.source, source)

        with self.assertRaisesRegex(
            ConfigurationError, "source must be one of available_sources"
        ):
            load_configuration(cli_overrides={
                "modules.harmonisation.default_eaf.source": "unknown_panel",
            })

    def test_missing_vcf_id_format_uses_bcftools_escaped_separators(self):
        config = load_configuration().modules.harmonisation.vcf_processing

        self.assertEqual(
            config.missing_id_format,
            "+%CHROM\\_%POS\\_%REF\\_%ALT",
        )

        with self.assertRaisesRegex(ConfigurationError, "must start with"):
            load_configuration(cli_overrides={
                "modules.harmonisation.vcf_processing.missing_id_format":
                    "%CHROM_%POS",
            })

    def test_eaf_cutoff_has_one_canonical_policy_location(self):
        config = load_configuration(
            cli_overrides={
                "modules.harmonisation.policies.eaf.maf_decision_cutoff": 0.9
            }
        )
        defaults = _engine_defaults(config)
        self.assertEqual(defaults["policies"]["eaf"]["maf_decision_cutoff"], 0.9)
        self.assertNotIn("maf_eaf_decision_cutoff", defaults)

    def test_maf_reference_confirmation_has_conservative_canonical_defaults(self):
        defaults = default_policies()

        self.assertEqual(defaults.get("eaf.reference_minor_fraction_cutoff"), 0.95)
        self.assertEqual(defaults.get("eaf.maf_reference_error_margin"), 0.02)
        self.assertNotIn("eaf.on_maf_check_inconclusive", defaults)

    @patch(
        "postgwas.modules.harmonisation.service.resolve_resource_file",
        side_effect=lambda input_file, **_: input_file,
    )
    @patch("postgwas.modules.harmonisation.service.validate_path")
    def test_resource_map_validates_default_eaf_only_when_required(
        self, validate_path_mock, _find_resource_mock
    ):
        validated_paths = []
        validate_path_mock.return_value = validated_paths.append

        arguments = dict(
            chromosome="1",
            grch_version="GRCh37",
            target_build="GRCh38",
            resource_folder="/reference",
            user_eaf_file="external_eaf.tsv",
            default_eaf_reference_source="1000G",
            default_comparison_af_file="ALFA",
            resource_layout=_engine_defaults(load_configuration())["resource_layout"],
            user_info_file="external_info.tsv",
            dbsnp="dbSNP157",
            user_eaf_column=None,
            default_eaf_reference_column="EUR",
            default_comparison_af_column=None,
            user_info_column=None,
        )
        resources = build_resource_map(**arguments, require_default_eaf=True)

        self.assertIn("default_eaf_file", resources)
        self.assertNotIn("default_info_file", resources)
        self.assertEqual(len(validated_paths), 7)
        self.assertIn("GRCh37_1000G_freq_chr1.tsv.gz", validated_paths[0])
        self.assertTrue(any(
            "GRCh37_ALFA_freq_chr1.vcf.gz" in path for path in validated_paths
        ))

        validated_paths.clear()
        build_resource_map(**arguments, require_default_eaf=False)
        self.assertEqual(len(validated_paths), 6)
        self.assertFalse(any(path.endswith(".tsv.gz") for path in validated_paths))

    def test_resource_layout_can_be_overridden_in_configuration(self):
        config = load_configuration(
            cli_overrides={
                "modules.harmonisation.resource_layout.fasta":
                    "genomes/{build}/{chromosome}.fa"
            }
        )
        self.assertEqual(
            _engine_defaults(config)["resource_layout"]["fasta"],
            "genomes/{build}/{chromosome}.fa",
        )

    def test_all_reference_table_mappings_reach_the_engine(self):
        config = load_configuration(
            cli_overrides={
                "modules.harmonisation.external_info_mapping.chromosome": "INFO_CHR",
                "modules.harmonisation.external_eaf_mapping.delimiter": "comma",
                "modules.harmonisation.build_check_mapping.position": "BUILD_POS",
            }
        )
        defaults = _engine_defaults(config)

        self.assertEqual(defaults["external_info_colmap"]["chr"], "INFO_CHR")
        self.assertEqual(defaults["external_eaf_colmap"]["delimiter"], "comma")
        self.assertEqual(defaults["build_check_colmap"]["pos"], "BUILD_POS")

    def test_build_inference_uses_configured_build_names_columns_and_delimiter(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.csv"
            second = root / "second.csv"
            first.write_text(
                "C,P,R,A\n1,100,A,G\n1,200,C,T\n", encoding="utf-8"
            )
            second.write_text(
                "C,P,R,A\n1,101,A,G\n1,201,C,T\n", encoding="utf-8"
            )
            result = infer_genome_build(
                pl.DataFrame({
                    "chr": ["1", "1"],
                    "pos": [100, 200],
                    "ea": ["A", "T"],
                    "oa": ["G", "C"],
                }),
                {"ReferenceA": str(first), "ReferenceB": str(second)},
                {
                    "chr": "C", "pos": "P", "a1": "R", "a2": "A",
                    "delimiter": "comma",
                },
                {"chr_col": "chr", "pos_col": "pos", "ea_col": "ea", "oa_col": "oa"},
                policies=default_policies().with_overrides({
                    "build.min_match_count": 1,
                }),
            )

        self.assertEqual(result["inferred_build"], "ReferenceA")
        self.assertEqual(result["matches"], {"ReferenceA": 2, "ReferenceB": 0})

    def test_build_inference_log_identifies_every_reference_file_and_mapping(self):
        from tempfile import TemporaryDirectory

        from postgwas.core.pipeline_logging import PipelineLogger

        with TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.csv"
            second = root / "second.csv"
            first.write_text("C,P,R,A\n1,100,G,A\n", encoding="utf-8")
            second.write_text("C,P,R,A\n1,101,G,A\n", encoding="utf-8")
            logger = PipelineLogger(
                "study", "dataset", directory, policies=default_policies(),
                level="INFO", screen_level="ERROR",
            )
            infer_genome_build(
                pl.DataFrame({
                    "chr": ["1"], "pos": [100], "ea": ["A"], "oa": ["G"],
                }),
                {"ReferenceA": str(first), "ReferenceB": str(second)},
                {"chr": "C", "pos": "P", "a1": "R", "a2": "A", "delimiter": "comma"},
                {"chr_col": "chr", "pos_col": "pos", "ea_col": "ea", "oa_col": "oa"},
                logger=logger, policies=default_policies(),
            )
            logger.close()
            text = Path(logger.log_path).read_text(encoding="utf-8")

        self.assertIn("INPUT    Genome-build reference build=ReferenceA", text)
        self.assertIn("file=%s" % first, text)
        self.assertIn("file=%s" % second, text)
        self.assertIn("reference_allele_column=R", text)
        self.assertIn("alternate_allele_column=A", text)
        self.assertIn("coordinate-testable variants 1", text)
        self.assertIn("ReferenceA coordinate hits 1", text)
        self.assertIn("minimum allele matches", text)
        self.assertIn("100; minimum match fraction 80.00%", text)

    def test_build_inference_counts_each_study_row_once_per_reference(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.csv"
            second = root / "second.csv"
            first.write_text(
                "C,P,R,A\n1,100,A,T\n1,100,A,T\n1,100,T,A\n",
                encoding="utf-8",
            )
            second.write_text("C,P,R,A\n1,101,A,T\n", encoding="utf-8")
            policies = default_policies().with_overrides(
                {"build.deduplicate_reference": False}
            )
            result = infer_genome_build(
                pl.DataFrame({
                    "chr": ["1"],
                    "pos": [100],
                    "ea": ["A"],
                    "oa": ["T"],
                }),
                {"ReferenceA": str(first), "ReferenceB": str(second)},
                {
                    "chr": "C", "pos": "P", "a1": "R", "a2": "A",
                    "delimiter": "comma",
                },
                {"chr_col": "chr", "pos_col": "pos", "ea_col": "ea", "oa_col": "oa"},
                policies=policies,
            )

        self.assertEqual(result["matches"], {"ReferenceA": 1, "ReferenceB": 0})
        self.assertEqual(result["match_fraction"], 1.0)

    def test_build_inference_accepts_reverse_complement_and_reports_strand_evidence(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.csv"
            second = root / "second.csv"
            first.write_text("C,P,R,A\n1,100,G,A\n1,200,T,C\n", encoding="utf-8")
            second.write_text("C,P,R,A\n1,101,G,A\n1,201,T,C\n", encoding="utf-8")
            result = infer_genome_build(
                pl.DataFrame({
                    "chr": ["1", "1"], "pos": [100, 200],
                    "ea": ["T", "G"], "oa": ["C", "A"],
                }),
                {"ReferenceA": str(first), "ReferenceB": str(second)},
                {"chr": "C", "pos": "P", "a1": "R", "a2": "A", "delimiter": "comma"},
                {"chr_col": "chr", "pos_col": "pos", "ea_col": "ea", "oa_col": "oa"},
                policies=default_policies().with_overrides({
                    "build.min_match_count": 1,
                }),
            )

        self.assertEqual(result["inferred_build"], "ReferenceA")
        self.assertEqual(result["matches"], {"ReferenceA": 2, "ReferenceB": 0})
        self.assertEqual(result["strand_evidence"]["ReferenceA"]["reverse"], 2)

    def test_build_defaults_require_substantial_evidence_and_deduplicate(self):
        policies = default_policies()

        self.assertEqual(policies.get("build.min_match_count"), 100)
        self.assertEqual(policies.get("build.min_match_fraction"), 0.8)
        self.assertEqual(policies.get("build.confidence_ratio"), 0.9)
        self.assertTrue(policies.get("build.deduplicate_reference"))
        with self.assertRaisesRegex(PolicyError, "build.min_match_count"):
            policies.with_overrides({"build.min_match_count": 0})

    def test_build_fraction_uses_only_exact_coordinate_testable_rows(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.csv"
            second = root / "second.csv"
            first.write_text("C,P,R,A\n1,100,A,G\n", encoding="utf-8")
            second.write_text("C,P,R,A\n1,101,A,G\n", encoding="utf-8")
            result = infer_genome_build(
                pl.DataFrame({
                    "chr": ["1", "1", "1"],
                    "pos": [100, 999, 1000],
                    "ea": ["A", "A", "A"],
                    "oa": ["G", "G", "G"],
                }),
                {"ReferenceA": str(first), "ReferenceB": str(second)},
                {"chr": "C", "pos": "P", "a1": "R", "a2": "A", "delimiter": "comma"},
                {"chr_col": "chr", "pos_col": "pos", "ea_col": "ea", "oa_col": "oa"},
                policies=default_policies().with_overrides({
                    "build.min_match_count": 1,
                }),
            )

        self.assertEqual(result["inferred_build"], "ReferenceA")
        self.assertEqual(result["input_variants"], 3)
        self.assertEqual(result["testable_variants"], 1)
        self.assertEqual(result["untestable_variants"], 2)
        self.assertEqual(result["coordinate_hits"], {"ReferenceA": 1, "ReferenceB": 0})
        self.assertEqual(result["percentages"], {"ReferenceA": 100.0, "ReferenceB": 0.0})

    def test_build_inference_rejects_evidence_below_absolute_default(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.csv"
            second = root / "second.csv"
            first.write_text(
                "C,P,R,A\n1,100,A,G\n1,101,A,G\n1,102,A,G\n",
                encoding="utf-8",
            )
            second.write_text("C,P,R,A\n1,200,A,G\n", encoding="utf-8")
            result = infer_genome_build(
                pl.DataFrame({
                    "chr": ["1", "1", "1"],
                    "pos": [100, 101, 102],
                    "ea": ["A", "A", "A"],
                    "oa": ["G", "G", "G"],
                }),
                {"ReferenceA": str(first), "ReferenceB": str(second)},
                {"chr": "C", "pos": "P", "a1": "R", "a2": "A", "delimiter": "comma"},
                {"chr_col": "chr", "pos_col": "pos", "ea_col": "ea", "oa_col": "oa"},
                policies=default_policies(),
            )

        self.assertEqual(result["inferred_build"], "Ambiguous")
        self.assertIn("build.min_match_count", result["ambiguous_reason"])
        self.assertIn("only 3", result["ambiguous_reason"])

    def test_build_inference_accepts_exact_absolute_default_boundary(self):
        from tempfile import TemporaryDirectory

        positions = list(range(100, 200))
        with TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.csv"
            second = root / "second.csv"
            first.write_text(
                "C,P,R,A\n" + "".join(
                    "1,%d,A,G\n" % position for position in positions
                ),
                encoding="utf-8",
            )
            second.write_text("C,P,R,A\n1,1000,A,G\n", encoding="utf-8")
            result = infer_genome_build(
                pl.DataFrame({
                    "chr": ["1"] * len(positions),
                    "pos": positions,
                    "ea": ["A"] * len(positions),
                    "oa": ["G"] * len(positions),
                }),
                {"ReferenceA": str(first), "ReferenceB": str(second)},
                {"chr": "C", "pos": "P", "a1": "R", "a2": "A", "delimiter": "comma"},
                {"chr_col": "chr", "pos_col": "pos", "ea_col": "ea", "oa_col": "oa"},
                policies=default_policies(),
            )

        self.assertEqual(result["inferred_build"], "ReferenceA")
        self.assertEqual(result["matches"]["ReferenceA"], 100)
        self.assertEqual(result["match_fraction"], 1.0)

    def test_build_inference_rejects_low_fraction_of_testable_rows(self):
        from tempfile import TemporaryDirectory

        positions = list(range(100, 110))
        reference_rows = [
            "1,%d,A,%s\n" % (position, "G" if index < 5 else "C")
            for index, position in enumerate(positions)
        ]
        with TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.csv"
            second = root / "second.csv"
            first.write_text("C,P,R,A\n" + "".join(reference_rows), encoding="utf-8")
            second.write_text("C,P,R,A\n1,1000,A,G\n", encoding="utf-8")
            result = infer_genome_build(
                pl.DataFrame({
                    "chr": ["1"] * len(positions),
                    "pos": positions,
                    "ea": ["A"] * len(positions),
                    "oa": ["G"] * len(positions),
                }),
                {"ReferenceA": str(first), "ReferenceB": str(second)},
                {"chr": "C", "pos": "P", "a1": "R", "a2": "A", "delimiter": "comma"},
                {"chr_col": "chr", "pos_col": "pos", "ea_col": "ea", "oa_col": "oa"},
                policies=default_policies().with_overrides({
                    "build.min_match_count": 1,
                }),
            )

        self.assertEqual(result["inferred_build"], "Ambiguous")
        self.assertEqual(result["testable_variants"], 10)
        self.assertEqual(result["matches"]["ReferenceA"], 5)
        self.assertEqual(result["match_fraction"], 0.5)
        self.assertIn("build.min_match_fraction", result["ambiguous_reason"])

    def test_declared_build_bypasses_automatic_evidence_thresholds(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            root = Path(directory)
            grch37 = root / "grch37.csv"
            grch38 = root / "grch38.csv"
            grch37.write_text("C,P,R,A\n1,100,A,G\n", encoding="utf-8")
            grch38.write_text("C,P,R,A\n1,101,A,G\n", encoding="utf-8")
            result = infer_genome_build(
                pl.DataFrame({
                    "chr": ["1"], "pos": [100], "ea": ["A"], "oa": ["G"],
                }),
                {"GRCh37": str(grch37), "GRCh38": str(grch38)},
                {"chr": "C", "pos": "P", "a1": "R", "a2": "A", "delimiter": "comma"},
                {"chr_col": "chr", "pos_col": "pos", "ea_col": "ea", "oa_col": "oa"},
                policies=default_policies().with_overrides({
                    "build.mode": "GRCh37",
                }),
            )

        self.assertEqual(result["inferred_build"], "GRCh37")
        self.assertTrue(result["forced"])
        self.assertEqual(result["matches"]["GRCh37"], 1)
        self.assertIsNone(result["ambiguous_reason"])

    def test_duplicate_header_failure_has_no_hidden_count_threshold(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            source = Path(directory) / "duplicate.csv"
            source.write_text("CHR,CHR\n1,2\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "repeated column names"):
                load_summary_statistics_table(
                    str(source), {"chr_col": "CHR"}, policies=default_policies()
                )

    def test_duplicate_header_warning_uses_first_column_when_configured(self):
        from tempfile import TemporaryDirectory

        policies = default_policies().with_overrides(
            {"validation.on_duplicate_header": "warn"}
        )
        with TemporaryDirectory() as directory:
            source = Path(directory) / "duplicate.csv"
            source.write_text("CHR,CHR\n1,2\n", encoding="utf-8")
            frame, _ = load_summary_statistics_table(
                str(source), {"chr_col": "CHR"}, policies=policies
            )

        self.assertEqual(frame["CHR"].to_list(), ["1"])
        self.assertIn("CHR_duplicated_0", frame.columns)


if __name__ == "__main__":
    unittest.main()
