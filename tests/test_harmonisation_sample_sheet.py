"""Validation of the portable harmonisation sample-sheet fixture."""

import unittest
from pathlib import Path
from unittest.mock import patch

from postgwas.core.errors import ConfigurationError
from postgwas.config import load_configuration
from postgwas.modules.harmonisation.sample_sheet import (
    load_harmonisation_sample_sheet,
    to_harmonisation_input,
)
from postgwas.modules.harmonisation.input_validation import validate_config
from postgwas.modules.harmonisation.policies import load_policies


FIXTURE_DIR = Path(__file__).parent / "data" / "harmonisation"
MANIFEST = FIXTURE_DIR / "manifest_v2.csv"


class HarmonisationSampleSheetTests(unittest.TestCase):
    @staticmethod
    def _write_manifest_with_values(path, **updates):
        lines = MANIFEST.read_text(encoding="utf-8").splitlines()
        header = lines[0].split(",")
        values = lines[1].split(",")
        for name, value in updates.items():
            values[header.index(name)] = value
        path.write_text(
            ",".join(header) + "\n" + ",".join(values) + "\n",
            encoding="utf-8",
        )

    def test_adhd_sample_sheet_uses_the_exact_fixture_path(self):
        rows = load_harmonisation_sample_sheet(MANIFEST)
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row.dataset_id, "ADHD2022_iPSYCH_deCODE_PGC")
        self.assertEqual(row.trait_type, "case_control")
        self.assertEqual(row.effect_type, "odds_ratio")
        self.assertTrue(row.input_file.is_file())
        self.assertEqual(
            row.input_file,
            Path(
                "/Users/JJOHN41/Documents/developing_software/postgwas_v2/"
                "tests/data/harmonisation/adhd2022_full_gwas.meta.gz"
            ),
        )

    def test_preflight_checks_metadata_without_reading_gwas_content(self):
        original_open = Path.open
        fixture = (FIXTURE_DIR / "adhd2022_full_gwas.meta.gz").resolve()

        def guarded_open(path, *args, **kwargs):
            if path.resolve() == fixture:
                raise AssertionError("sample-sheet preflight opened the GWAS file")
            return original_open(path, *args, **kwargs)

        with patch.object(Path, "open", guarded_open):
            row = load_harmonisation_sample_sheet(MANIFEST)[0]
        self.assertGreater(row.input_file.stat().st_size, 0)

    def test_engine_input_retains_sample_sheet_metadata_for_reporting(self):
        row = load_harmonisation_sample_sheet(MANIFEST)[0]
        values = to_harmonisation_input(
            row,
            resource_directory=FIXTURE_DIR,
            output_directory=FIXTURE_DIR,
            output_layout=dict(
                load_configuration().modules.harmonisation.output_layout.root
            ),
        )

        self.assertEqual(values["trait_type"], "case_control")
        self.assertEqual(values["declared_effect_type"], "odds_ratio")
        self.assertEqual(values["declared_pvalue_type"], "raw")
        self.assertEqual(values["delimiter"], "whitespace")
        self.assertIsNone(values["chr_pos_col"])
        self.assertIsNone(values["imp_z_col"])
        self.assertIsNone(values["infofile"])
        self.assertIsNone(values["provided_external_info_file"])
        self.assertIsNone(values["provided_external_eaf_file"])
        self.assertNotIn("NA", values.values())
        self.assertEqual(
            Path(values["output_folder"]),
            FIXTURE_DIR / row.dataset_id / "harmonisation",
        )
        self.assertEqual(Path(values["output_root"]), FIXTURE_DIR)

    def test_engine_output_root_uses_the_configured_dataset_layout(self):
        row = load_harmonisation_sample_sheet(MANIFEST)[0]
        values = to_harmonisation_input(
            row,
            resource_directory=FIXTURE_DIR,
            output_directory=FIXTURE_DIR,
            output_layout={"dataset_directory": "datasets/{dataset_id}/custom"},
        )
        self.assertEqual(
            Path(values["output_folder"]),
            FIXTURE_DIR / "datasets" / row.dataset_id / "custom",
        )

    def test_cli_fixed_info_is_used_only_when_the_sample_sheet_has_no_source(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            temporary = Path(directory) / "manifest.csv"
            self._write_manifest_with_values(
                temporary,
                imputation_info_column="NA",
                external_info_file="NA",
                external_info_column="NA",
            )
            row = load_harmonisation_sample_sheet(temporary)[0]
            with self.assertRaisesRegex(ConfigurationError, "--fixed-info"):
                to_harmonisation_input(
                    row,
                    resource_directory=FIXTURE_DIR,
                    output_directory=FIXTURE_DIR,
                    output_layout={"dataset_directory": "{dataset_id}"},
                )
            values = to_harmonisation_input(
                row,
                resource_directory=FIXTURE_DIR,
                output_directory=FIXTURE_DIR,
                output_layout={"dataset_directory": "{dataset_id}"},
                fixed_info=0.99,
                fixed_info_column="__postgwas_fixed_info",
            )

        self.assertEqual(values["fixed_info"], 0.99)
        self.assertEqual(values["fixed_info_column"], "__postgwas_fixed_info")
        self.assertEqual(values["info_source"], "fixed_cli")
        self.assertIsNone(values["imp_info_col"])
        self.assertIsNone(values["infofile"])
        self.assertIsNone(values["infocolumn"])

    def test_internal_info_has_priority_over_external_and_fixed_sources(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            temporary = Path(directory) / "manifest.csv"
            self._write_manifest_with_values(
                temporary,
                external_info_file="does_not_exist.tsv",
                external_info_column="REFERENCE_INFO",
            )
            row = load_harmonisation_sample_sheet(temporary)[0]
            values = to_harmonisation_input(
                row,
                resource_directory=FIXTURE_DIR,
                output_directory=FIXTURE_DIR,
                output_layout={"dataset_directory": "{dataset_id}"},
                fixed_info=0.99,
                fixed_info_column="__postgwas_fixed_info",
            )

        self.assertEqual(values["imp_info_col"], "INFO")
        self.assertIsNone(values["infofile"])
        self.assertIsNone(values["infocolumn"])
        self.assertIsNone(values["fixed_info"])
        self.assertEqual(values["info_source"], "internal")
        self.assertTrue(
            values["provided_external_info_file"].endswith("does_not_exist.tsv")
        )
        self.assertIn("internal column INFO", values["info_source_detail"])
        self.assertIn("has priority", row.normalisation_warnings[0])

    def test_external_info_has_priority_over_the_cli_fixed_fallback(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            temporary = Path(directory) / "manifest.csv"
            self._write_manifest_with_values(
                temporary,
                imputation_info_column="NA",
                external_info_file=str(MANIFEST),
                external_info_column="REFERENCE_INFO",
            )
            row = load_harmonisation_sample_sheet(temporary)[0]
            values = to_harmonisation_input(
                row,
                resource_directory=FIXTURE_DIR,
                output_directory=FIXTURE_DIR,
                output_layout={"dataset_directory": "{dataset_id}"},
                fixed_info=0.99,
                fixed_info_column="__postgwas_fixed_info",
            )

        self.assertEqual(values["info_source"], "external")
        self.assertEqual(values["infofile"], str(MANIFEST.resolve()))
        self.assertEqual(values["infocolumn"], "REFERENCE_INFO")
        self.assertIsNone(values["fixed_info"])

    def test_missing_eaf_error_names_dataset_and_fields_to_complete(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            temporary = Path(directory) / "manifest.csv"
            self._write_manifest_with_values(
                temporary,
                effect_allele_frequency_column="NA",
                external_eaf_file="NA",
                external_eaf_column="NA",
            )
            with self.assertRaises(ConfigurationError) as caught:
                load_harmonisation_sample_sheet(temporary)

        message = str(caught.exception)
        self.assertIn(
            "sample-sheet row 2 for dataset 'ADHD2022_iPSYCH_deCODE_PGC'",
            message,
        )
        self.assertIn("ADHD2022_iPSYCH_deCODE_PGC", message)
        self.assertIn("has no allele-frequency source", message)
        self.assertIn("effect_allele_frequency_column", message)
        self.assertIn("external_eaf_file", message)
        self.assertIn("external_eaf_column", message)
        self.assertNotIn("pydantic.dev", message)

    def test_partial_and_conflicting_eaf_errors_explain_the_correction(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            partial = Path(directory) / "partial.csv"
            self._write_manifest_with_values(
                partial,
                effect_allele_frequency_column="NA",
                external_eaf_file=str(MANIFEST),
                external_eaf_column="NA",
            )
            with self.assertRaisesRegex(
                ConfigurationError,
                "incomplete external EAF source: external_eaf_column is missing",
            ):
                load_harmonisation_sample_sheet(partial)

            conflicting = Path(directory) / "conflicting.csv"
            self._write_manifest_with_values(
                conflicting,
                external_eaf_file=str(MANIFEST),
                external_eaf_column="EUR",
            )
            with self.assertRaisesRegex(
                ConfigurationError,
                "provides two allele-frequency sources",
            ):
                load_harmonisation_sample_sheet(conflicting)

    def test_external_chromosome_templates_are_accepted_and_prefixes_explained(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            root = Path(directory)
            frequency = root / "panel_chr1.tsv.gz"
            frequency.write_text("CHROM\tPOS\tREF\tALT\tEUR\n", encoding="utf-8")
            template = root / "panel_chr{chromosome}.tsv.gz"
            valid = root / "valid.csv"
            self._write_manifest_with_values(
                valid,
                effect_allele_frequency_column="NA",
                external_eaf_file=str(template),
                external_eaf_column="EUR",
            )

            row = load_harmonisation_sample_sheet(valid)[0]
            self.assertEqual(row.external_eaf_file, template.resolve(strict=False))

            prefix = root / "panel"
            invalid = root / "invalid.csv"
            self._write_manifest_with_values(
                invalid,
                effect_allele_frequency_column="NA",
                external_eaf_file=str(prefix),
                external_eaf_column="EUR",
            )
            with self.assertRaises(ConfigurationError) as caught:
                load_harmonisation_sample_sheet(invalid)

        message = str(caught.exception)
        self.assertIn("filename prefix, not a file", message)
        self.assertIn("panel_chr1.tsv.gz", message)
        self.assertIn("{chromosome}", message)

    def test_engine_config_validation_accepts_external_chromosome_template(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            root = Path(directory)
            frequency = root / "panel_chr1.tsv"
            frequency.write_text(
                "CHROM\tPOS\tREF\tALT\tEUR\n1\t10\tA\tG\t0.2\n",
                encoding="utf-8",
            )
            template = root / "panel_chr{chromosome}.tsv"
            manifest = root / "manifest.csv"
            self._write_manifest_with_values(
                manifest,
                effect_allele_frequency_column="NA",
                external_eaf_file=str(template),
                external_eaf_column="EUR",
            )
            row = load_harmonisation_sample_sheet(manifest)[0]
            config = load_configuration()
            values = to_harmonisation_input(
                row,
                resource_directory=FIXTURE_DIR,
                output_directory=root,
                output_layout={"dataset_directory": "{dataset_id}"},
            )

            ok, problems = validate_config(
                values,
                policies=load_policies(config.modules.harmonisation.policies),
                external_reference_delimiters={"eaffile": "auto"},
            )

        self.assertTrue(ok, [str(problem) for problem in problems])

    def test_external_reference_delimiter_is_independent_of_study_delimiter(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            root = Path(directory)
            frequency = root / "panel_chr1.tsv"
            frequency.write_text(
                "CHROM\tPOS\tREF\tALT\tEUR\n1\t10\tA\tG\t0.2\n",
                encoding="utf-8",
            )
            template = root / "panel_chr{chromosome}.tsv"
            manifest = root / "manifest.csv"
            self._write_manifest_with_values(
                manifest,
                delimiter="space",
                effect_allele_frequency_column="NA",
                external_eaf_file=str(template),
                external_eaf_column="EUR",
            )
            row = load_harmonisation_sample_sheet(manifest)[0]
            config = load_configuration()
            values = to_harmonisation_input(
                row,
                resource_directory=FIXTURE_DIR,
                output_directory=root,
                output_layout={"dataset_directory": "{dataset_id}"},
            )
            policies = load_policies(
                config.modules.harmonisation.policies
            ).with_overrides({"input.delimiter": "space"})

            ok, problems = validate_config(
                values,
                policies=policies,
                external_reference_delimiters={"eaffile": "auto"},
            )

        self.assertTrue(ok, [str(problem) for problem in problems])

    def test_quantitative_trait_rejects_case_counts(self):
        text = MANIFEST.read_text(encoding="utf-8")
        text = text.replace("case_control", "quantitative", 1)
        with self.subTest("schema validation"):
            from tempfile import TemporaryDirectory

            with TemporaryDirectory() as directory:
                temporary = Path(directory) / "manifest.csv"
                temporary.write_text(text, encoding="utf-8")
                with self.assertRaisesRegex(ConfigurationError, "must not provide case-count"):
                    load_harmonisation_sample_sheet(temporary)

    def test_invalid_inferable_values_warn_and_become_auto(self):
        text = MANIFEST.read_text(encoding="utf-8")
        text = text.replace("case_control", "unknown_trait", 1)
        text = text.replace("odds_ratio", "unknown_effect", 1)
        text = text.replace(",raw,", ",unknown_pvalue,", 1)
        text = text.replace(",whitespace\n", ",unknown_delimiter\n", 1)
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            temporary = Path(directory) / "manifest.csv"
            temporary.write_text(text, encoding="utf-8")
            row = load_harmonisation_sample_sheet(temporary)[0]
        self.assertEqual(row.trait_type, "auto")
        self.assertEqual(row.effect_type, "auto")
        self.assertEqual(row.p_value_type, "auto")
        self.assertEqual(row.delimiter, "auto")
        self.assertEqual(len(row.normalisation_warnings), 4)

    def test_natural_log_pvalue_declaration_is_not_supported(self):
        text = MANIFEST.read_text(encoding="utf-8").replace(
            ",raw,", ",negln,", 1
        )
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            temporary = Path(directory) / "manifest.csv"
            temporary.write_text(text, encoding="utf-8")
            row = load_harmonisation_sample_sheet(temporary)[0]

        self.assertEqual(row.p_value_type, "auto")
        self.assertTrue(
            any("p_value_type" in warning for warning in row.normalisation_warnings)
        )

    def test_duplicate_dataset_ids_are_case_insensitive(self):
        lines = MANIFEST.read_text(encoding="utf-8").splitlines()
        second = lines[1].replace(
            "ADHD2022_iPSYCH_deCODE_PGC",
            "adhd2022_ipsych_decode_pgc",
            1,
        ).replace(str(FIXTURE_DIR / "adhd2022_full_gwas.meta.gz"), str(MANIFEST), 1)
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            temporary = Path(directory) / "manifest.csv"
            temporary.write_text("\n".join([lines[0], lines[1], second]) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ConfigurationError, "Duplicate dataset_id"):
                load_harmonisation_sample_sheet(temporary)

    def test_duplicate_normalized_input_files_are_rejected(self):
        lines = MANIFEST.read_text(encoding="utf-8").splitlines()
        second = lines[1].replace("ADHD2022_iPSYCH_deCODE_PGC", "ADHD_SECOND", 1)
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            temporary = Path(directory) / "manifest.csv"
            temporary.write_text("\n".join([lines[0], lines[1], second]) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ConfigurationError, "Duplicate input_file"):
                load_harmonisation_sample_sheet(temporary)

    def test_header_names_are_normalized(self):
        text = MANIFEST.read_text(encoding="utf-8").replace(
            "dataset_id", " Dataset-ID ", 1
        )
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as directory:
            temporary = Path(directory) / "manifest.csv"
            temporary.write_text(text, encoding="utf-8")
            row = load_harmonisation_sample_sheet(temporary)[0]
        self.assertEqual(row.dataset_id, "ADHD2022_iPSYCH_deCODE_PGC")


if __name__ == "__main__":
    unittest.main()
