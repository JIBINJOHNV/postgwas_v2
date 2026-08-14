"""Scientific and orchestration regression tests for MAGMA."""

import argparse
from argparse import Namespace
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import polars as pl
import pytest
import yaml

from postgwas.config import (
    load_configuration,
    load_module_configuration,
    resolved_configuration_values,
)
from postgwas.core.statistics import adjust_p_values
from postgwas.modules.magma.cli import build_parser
from postgwas.modules.magma.errors import MagmaError
from postgwas.modules.magma.analysis import (
    _batch_plan,
    _significance_outcome,
    annotate_gene_set_results,
    correct_gene_p_values,
    correct_gene_set_p_values,
    parse_gene_set_file,
    prepare_magma_variant_inputs,
    resolve_gene_set_identifiers,
)
from postgwas.modules.magma.annotations import (
    map_regulatory_elements_to_genes,
    merge_gene_annotations,
    validate_gene_annotation,
)


class RecordingLogger:
    def __init__(self):
        self.records = []

    def record(self, marker, subject, **values):
        self.records.append((marker, subject, values))


def _reference_files(tmp_path: Path) -> Path:
    prefix = tmp_path / "reference"
    prefix.with_suffix(".bim").write_text(
        "1 rs1 0 100 A G\n"
        "1 rs2 0 200 C T\n"
        "1 rs3 0 300 A C\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".bed").write_bytes(b"BED")
    prefix.with_suffix(".fam").write_text("family sample 0 0 0 -9\n", encoding="utf-8")
    return prefix


def _magma_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    reference = _reference_files(tmp_path)
    p_values = tmp_path / "p_values.tsv"
    p_values.write_text(
        "SNP\tP\tN_COL\n"
        "rs1\t0.05\t1000\n"
        "rs2\t0.01\t1000\n"
        "rs2\t0.001\t1000\n"
        "rs3\t0.02\t1000\n"
        "absent\t0.03\t1000\n",
        encoding="utf-8",
    )
    locations = tmp_path / "locations.tsv"
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\n"
        "rs1\t1\t100\tA\tG\n"
        "rs2\t1\t200\tC\tT\n"
        "rs3\t1\t300\tC\tA\n"
        "absent\t1\t400\tG\tT\n",
        encoding="utf-8",
    )
    genes = tmp_path / "genes.loc"
    genes.write_text("1 1 50 350 +\n", encoding="utf-8")
    return locations, p_values, reference, genes


def test_cli_has_no_independent_defaults_or_analysis_imports():
    parser = build_parser()
    for destination in (
        "snp_location_file", "p_value_file", "magma_ld_reference",
        "gene_location_file", "gene_set_file",
        "window_upstream", "window_downstream", "gene_model",
        "sample_size_column", "minimum_snp_overlap",
        "minimum_gene_id_overlap", "resolve_variants_to_reference",
        "magma_mapping", "primary_magma_mapping",
        "magma",
        "run_config", "resume", "overwrite", "dataset_id",
        "output_directory", "threads", "memory_gb", "seed",
    ):
        action = next(item for item in parser._actions if item.dest == destination)
        assert action.default == argparse.SUPPRESS
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; from postgwas.modules.magma.cli import build_parser; "
            "build_parser(); "
            "assert 'postgwas.modules.magma.service' not in sys.modules; "
            "assert 'statsmodels' not in sys.modules",
        ],
        check=True,
    )


def test_gene_model_validation_matches_magma_p_value_models(tmp_path):
    valid = tmp_path / "valid.yaml"
    valid.write_text("gene_model: snp-wise=top,0.1\n", encoding="utf-8")
    assert load_module_configuration("magma", valid).gene_model == "snp-wise=top,0.1"

    invalid = tmp_path / "invalid.yaml"
    invalid.write_text("gene_model: snp-wise=all\n", encoding="utf-8")
    with pytest.raises(Exception, match="supported with SNP p-value input"):
        load_module_configuration("magma", invalid)


def test_bim_extension_must_be_a_required_reference_companion(tmp_path):
    invalid = tmp_path / "invalid.yaml"
    invalid.write_text("input:\n  bim_extension: .map\n", encoding="utf-8")
    with pytest.raises(Exception, match="bim_extension"):
        load_module_configuration("magma", invalid)


def test_prepared_magma_table_delimiter_is_restricted_to_whitespace(tmp_path):
    invalid = tmp_path / "invalid.yaml"
    invalid.write_text(
        "input:\n  output_table_delimiter: ','\n",
        encoding="utf-8",
    )
    with pytest.raises(Exception, match="output_table_delimiter"):
        load_module_configuration("magma", invalid)


def test_configuration_precedence_is_canonical_yaml_then_explicit_cli(tmp_path):
    from postgwas.modules.magma.service import _resolved_configuration

    config_file = tmp_path / "magma.yaml"
    config_file.write_text("gene_window_upstream_kb: 50\n", encoding="utf-8")
    configuration = _resolved_configuration(
        Namespace(
            run_config=str(config_file),
            window_upstream=25,
            resolve_variants_to_reference=True,
        )
    )


def test_resolved_magma_metadata_omits_unselected_mapping_definitions():
    from postgwas.modules.magma.service import _resolved_magma_metadata_paths

    configuration = load_configuration()
    module = configuration.modules.magma
    module.mapping.definitions["unused"] = module.mapping.definitions[
        "positional"
    ].model_copy(update={"display_name": "Unused mapping"})

    observed = resolved_configuration_values(
        configuration,
        modules=("magma",),
        resource_paths=("executables.magma",),
        module_paths={"magma": _resolved_magma_metadata_paths(module)},
    )

    assert list(observed["modules"]["magma"]["mapping"]["definitions"]) == [
        "positional"
    ]
    assert configuration.modules.magma.gene_window_upstream_kb == 35
    assert configuration.modules.magma.gene_window_downstream_kb == 10
    assert (
        configuration.modules.magma.snp_harmonisation
        .resolve_variants_to_reference
        is False
    )


def test_reference_intersection_is_disabled_by_default(tmp_path, monkeypatch):
    locations, p_values, reference, _ = _magma_inputs(tmp_path)
    module = load_configuration().modules.magma
    assert module.snp_harmonisation.resolve_variants_to_reference is False
    monkeypatch.setattr(
        "postgwas.modules.magma.analysis.read_reference_bim_matches",
        lambda *args, **kwargs: pytest.fail("disabled intersection read the BIM"),
    )

    result = prepare_magma_variant_inputs(
        p_values,
        locations,
        reference,
        tmp_path / "p_out.tsv",
        tmp_path / "loc_out.tsv",
        module,
        RecordingLogger(),
    )

    assert pd.read_csv(tmp_path / "p_out.tsv", sep="\t")["SNP"].tolist() == [
        "rs1", "rs2", "rs3", "absent",
    ]
    assert result["qc"] == {
        "input_rows": 5,
        "input_unique_variants": 4,
        "reference_intersection_enabled": False,
        "reference_variant_count": None,
        "reference_id_match_rows": None,
        "reference_unique_id_matches": None,
        "reference_compatible_rows": None,
        "reference_unique_compatible_variants": None,
        "reference_coordinate_mismatch_rows": None,
        "reference_allele_mismatch_rows": None,
        "not_in_reference_rows": None,
        "overlap_fraction": None,
        "duplicates_resolved_by_lowest_p": 1,
        "retained_rows": 4,
    }


def test_default_mode_accepts_and_validates_existing_bim_ids(tmp_path):
    reference = _reference_files(tmp_path)
    p_values = tmp_path / "p.tsv"
    p_values.write_text(
        "SNP\tP\tN_COL\nrs1\t0.1\t100\nrs2\t0.2\t100\n",
        encoding="utf-8",
    )
    locations = tmp_path / "loc.tsv"
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\n"
        "rs1\t1\t100\tA\tG\n"
        "rs2\t1\t200\tC\tT\n",
        encoding="utf-8",
    )

    result = prepare_magma_variant_inputs(
        p_values,
        locations,
        reference,
        tmp_path / "p_out.tsv",
        tmp_path / "loc_out.tsv",
        load_configuration().modules.magma,
        RecordingLogger(),
    )

    assert pd.read_csv(tmp_path / "p_out.tsv", sep="\t")["SNP"].tolist() == [
        "rs1", "rs2",
    ]
    assert result["qc"]["reference_intersection_enabled"] is False
    assert result["qc"]["reference_id_match_rows"] is None
    assert result["qc"]["reference_unique_id_matches"] is None
    assert result["qc"]["reference_compatible_rows"] is None
    assert result["qc"]["reference_unique_compatible_variants"] is None
    assert result["qc"]["not_in_reference_rows"] is None


def test_reference_intersection_is_allele_aware_and_keeps_lowest_duplicate_p(
    tmp_path,
):
    locations, p_values, reference, _ = _magma_inputs(tmp_path)
    p_output = tmp_path / "output" / "p_values.tsv"
    location_output = tmp_path / "output" / "locations.tsv"
    logger = RecordingLogger()

    module = load_configuration().modules.magma
    module.snp_harmonisation.resolve_variants_to_reference = True
    result = prepare_magma_variant_inputs(
        p_values,
        locations,
        reference,
        p_output,
        location_output,
        module,
        logger,
    )

    observed = pd.read_csv(p_output, sep="\t")
    assert observed["SNP"].tolist() == ["rs1", "rs2", "rs3"]
    assert observed.loc[observed["SNP"] == "rs2", "P"].item() == pytest.approx(0.001)
    assert result["qc"] == {
        "input_rows": 5,
        "input_unique_variants": 4,
        "reference_intersection_enabled": True,
        "reference_variant_count": 3,
        "reference_id_match_rows": 4,
        "reference_unique_id_matches": 3,
        "reference_compatible_rows": 4,
        "reference_unique_compatible_variants": 3,
        "reference_coordinate_mismatch_rows": 0,
        "reference_allele_mismatch_rows": 0,
        "not_in_reference_rows": 1,
        "overlap_fraction": pytest.approx(0.75),
        "duplicates_resolved_by_lowest_p": 1,
        "retained_rows": 3,
    }
    first_location = location_output.read_text(encoding="utf-8").splitlines()[0]
    assert first_location == "rs1\t1\t100"
    assert not first_location.startswith("SNP")


def test_reference_intersection_enforces_configured_unique_overlap(tmp_path):
    locations, p_values, reference, _ = _magma_inputs(tmp_path)
    module = load_configuration().modules.magma
    module.snp_harmonisation.resolve_variants_to_reference = True
    module.snp_harmonisation.minimum_overlap_fraction = 0.8

    with pytest.raises(
        MagmaError,
        match=r"Only 3/4 unique formatter variants \(75.00%\).*80.00%",
    ):
        prepare_magma_variant_inputs(
            p_values,
            locations,
            reference,
            tmp_path / "p_out.tsv",
            tmp_path / "loc_out.tsv",
            module,
            RecordingLogger(),
        )


def test_reference_intersection_accepts_formatter_created_unique_ids(tmp_path):
    reference = tmp_path / "reference"
    reference.with_suffix(".bim").write_text(
        "1 1_100_A_G 0 100 G A\n"
        "1 1_200_C_T 0 200 T C\n",
        encoding="utf-8",
    )
    reference.with_suffix(".bed").write_bytes(b"BED")
    reference.with_suffix(".fam").write_text(
        "family sample 0 0 0 -9\n", encoding="utf-8",
    )
    p_values = tmp_path / "p.tsv"
    p_values.write_text(
        "SNP\tP\tN_COL\n"
        "1_100_A_G\t0.01\t1000\n"
        "1_200_C_T\t0.02\t1000\n",
        encoding="utf-8",
    )
    locations = tmp_path / "locations.tsv"
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\n"
        "1_100_A_G\t1\t100\tA\tG\n"
        "1_200_C_T\t1\t200\tC\tT\n",
        encoding="utf-8",
    )

    module = load_configuration().modules.magma
    module.snp_harmonisation.resolve_variants_to_reference = True
    result = prepare_magma_variant_inputs(
        p_values,
        locations,
        reference,
        tmp_path / "p_out.tsv",
        tmp_path / "loc_out.tsv",
        module,
        RecordingLogger(),
    )

    assert pd.read_csv(tmp_path / "p_out.tsv", sep="\t")["SNP"].tolist() == [
        "1_100_A_G", "1_200_C_T",
    ]
    assert result["qc"]["reference_id_match_rows"] == 2
    assert result["qc"]["reference_unique_id_matches"] == 2
    assert result["qc"]["reference_compatible_rows"] == 2
    assert result["qc"]["reference_unique_compatible_variants"] == 2


def test_coordinate_identifier_construction_rejects_fractional_positions(tmp_path):
    reference = _reference_files(tmp_path)
    p_values = tmp_path / "p.tsv"
    p_values.write_text("SNP\tP\tN_COL\nunknown\t0.1\t100\n", encoding="utf-8")
    locations = tmp_path / "loc.tsv"
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\nunknown\t1\t300.5\tA\tC\n",
        encoding="utf-8",
    )

    with pytest.raises(MagmaError, match="non-integer positions"):
        prepare_magma_variant_inputs(
            p_values,
            locations,
            reference,
            tmp_path / "p_out.tsv",
            tmp_path / "loc_out.tsv",
            load_configuration().modules.magma,
            RecordingLogger(),
        )


def test_variant_preparation_rejects_configured_invalid_chromosome_label(tmp_path):
    reference = _reference_files(tmp_path)
    p_values = tmp_path / "p.tsv"
    p_values.write_text("SNP\tP\tN_COL\nrs1\t0.1\t100\n", encoding="utf-8")
    locations = tmp_path / "loc.tsv"
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\nrs1\t.\t100\tA\tG\n",
        encoding="utf-8",
    )

    with pytest.raises(MagmaError, match="invalid chromosome labels"):
        prepare_magma_variant_inputs(
            p_values,
            locations,
            reference,
            tmp_path / "p_out.tsv",
            tmp_path / "loc_out.tsv",
            load_configuration().modules.magma,
            RecordingLogger(),
        )


def test_snp_location_always_requires_alleles(tmp_path):
    reference = _reference_files(tmp_path)
    p_values = tmp_path / "p.tsv"
    p_values.write_text("SNP\tP\tN_COL\nunknown\t0.1\t100\n", encoding="utf-8")
    locations = tmp_path / "loc.tsv"
    locations.write_text("SNP\tCHR\tBP\nunknown\t1\t300\n", encoding="utf-8")

    with pytest.raises(MagmaError, match="location table: ALT, REF"):
        prepare_magma_variant_inputs(
            p_values,
            locations,
            reference,
            tmp_path / "p_out.tsv",
            tmp_path / "loc_out.tsv",
            load_configuration().modules.magma,
            RecordingLogger(),
        )


def test_every_p_value_row_requires_location_and_allele_evidence(tmp_path):
    reference = _reference_files(tmp_path)
    p_values = tmp_path / "p.tsv"
    p_values.write_text("SNP\tP\tN_COL\nrs1\t0.1\t100\n", encoding="utf-8")
    locations = tmp_path / "loc.tsv"
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\nrs2\t1\t200\tC\tT\n",
        encoding="utf-8",
    )

    with pytest.raises(MagmaError, match="no matching SNP-location record"):
        prepare_magma_variant_inputs(
            p_values,
            locations,
            reference,
            tmp_path / "p_out.tsv",
            tmp_path / "loc_out.tsv",
            load_configuration().modules.magma,
            RecordingLogger(),
        )


def test_duplicate_location_identifier_cannot_hide_conflicting_alleles(tmp_path):
    reference = _reference_files(tmp_path)
    p_values = tmp_path / "p.tsv"
    p_values.write_text("SNP\tP\tN_COL\nrs1\t0.1\t100\n", encoding="utf-8")
    locations = tmp_path / "loc.tsv"
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\n"
        "rs1\t1\t100\tA\tG\n"
        "rs1\t1\t100\tC\tT\n",
        encoding="utf-8",
    )

    with pytest.raises(MagmaError, match="conflicting coordinates or allele pairs"):
        prepare_magma_variant_inputs(
            p_values,
            locations,
            reference,
            tmp_path / "p_out.tsv",
            tmp_path / "loc_out.tsv",
            load_configuration().modules.magma,
            RecordingLogger(),
        )


def test_reference_intersection_excludes_id_match_with_incompatible_alleles(tmp_path):
    reference = _reference_files(tmp_path)
    p_values = tmp_path / "p.tsv"
    p_values.write_text(
        "SNP\tP\tN_COL\nrs1\t0.1\t100\nrs2\t0.2\t100\n",
        encoding="utf-8",
    )
    locations = tmp_path / "loc.tsv"
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\n"
        "rs1\t1\t100\tC\tT\n"
        "rs2\t1\t200\tC\tT\n",
        encoding="utf-8",
    )

    module = load_configuration().modules.magma
    module.snp_harmonisation.resolve_variants_to_reference = True
    result = prepare_magma_variant_inputs(
        p_values,
        locations,
        reference,
        tmp_path / "p_out.tsv",
        tmp_path / "loc_out.tsv",
        module,
        RecordingLogger(),
    )

    assert pd.read_csv(tmp_path / "p_out.tsv", sep="\t")["SNP"].tolist() == ["rs2"]
    assert result["qc"]["reference_unique_id_matches"] == 2
    assert result["qc"]["reference_unique_compatible_variants"] == 1
    assert result["qc"]["reference_coordinate_mismatch_rows"] == 0
    assert result["qc"]["reference_allele_mismatch_rows"] == 1


def test_reference_intersection_excludes_id_match_at_incompatible_coordinate(tmp_path):
    reference = _reference_files(tmp_path)
    p_values = tmp_path / "p.tsv"
    p_values.write_text(
        "SNP\tP\tN_COL\nrs1\t0.1\t100\nrs2\t0.2\t100\n",
        encoding="utf-8",
    )
    locations = tmp_path / "loc.tsv"
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\n"
        "rs1\t1\t101\tA\tG\n"
        "rs2\t1\t200\tC\tT\n",
        encoding="utf-8",
    )

    module = load_configuration().modules.magma
    module.snp_harmonisation.resolve_variants_to_reference = True
    result = prepare_magma_variant_inputs(
        p_values,
        locations,
        reference,
        tmp_path / "p_out.tsv",
        tmp_path / "loc_out.tsv",
        module,
        RecordingLogger(),
    )

    assert pd.read_csv(tmp_path / "p_out.tsv", sep="\t")["SNP"].tolist() == ["rs2"]
    assert result["qc"]["reference_unique_id_matches"] == 2
    assert result["qc"]["reference_unique_compatible_variants"] == 1
    assert result["qc"]["reference_coordinate_mismatch_rows"] == 1
    assert result["qc"]["reference_allele_mismatch_rows"] == 0


@pytest.mark.parametrize(
    ("method", "expected"),
    [
        ("bonferroni", [0.04, 0.16, 0.12, 0.008]),
        ("holm", [0.03, 0.06, 0.06, 0.008]),
        ("fdr_bh", [0.02, 0.04, 0.04, 0.008]),
    ],
)
def test_multiple_testing_corrections(method, expected):
    values = np.array([0.01, 0.04, 0.03, 0.002])
    np.testing.assert_allclose(adjust_p_values(values, method), expected)


def test_sidak_correction_is_numerically_stable():
    values = np.array([1e-300, 0.01])
    expected = -np.expm1(2 * np.log1p(-values))
    np.testing.assert_allclose(adjust_p_values(values, "sidak"), expected)


def test_gene_results_add_configured_bonferroni_and_fdr_without_reordering(tmp_path):
    module = load_module_configuration("magma")
    genes = tmp_path / "study.genes.out"
    genes.write_text(
        "# MAGMA gene analysis\n"
        "GENE CHR START STOP NSNPS NPARAM N ZSTAT P\n"
        "1 1 10 20 3 1 1000 2.0 0.01\n"
        "2 1 30 40 2 1 1000 0.0 1.0\n"
        "3 1 50 60 4 1 1000 8.0 0.0\n",
        encoding="utf-8",
    )
    output = tmp_path / "corrected.tsv"
    logger = RecordingLogger()

    corrected = correct_gene_p_values(genes, output, module, logger)
    observed = pd.read_csv(output, sep="\t")

    assert corrected.height == 3
    assert observed["GENE"].tolist() == [1, 2, 3]
    assert observed.columns.tolist()[-2:] == [
        "P_bonferroni_corr",
        "P_fdr_bh_corr",
    ]
    np.testing.assert_allclose(
        observed["P_bonferroni_corr"],
        [0.03, 1.0, 0.0],
    )
    np.testing.assert_allclose(
        observed["P_fdr_bh_corr"],
        [0.015, 1.0, 0.0],
    )
    correction_records = [
        values
        for marker, subject, values in logger.records
        if marker == "PARAM" and subject == "multiple_testing_family"
    ]
    assert correction_records == [
        {"family": "all_genes", "method": "bonferroni", "tests": 3},
        {"family": "all_genes", "method": "fdr_bh", "tests": 3},
    ]


def test_significance_summary_uses_configured_threshold_and_corrected_columns(
    tmp_path,
):
    config_file = tmp_path / "magma.yaml"
    config_file.write_text(
        "multiple_testing:\n  reporting_significance_threshold: 0.01\n",
        encoding="utf-8",
    )
    module = load_module_configuration("magma", config_file)
    frame = pl.DataFrame(
        {
            "P": [0.005, 0.01, 0.02],
            "P_bonferroni_corr": [0.015, 0.03, 0.06],
            "P_fdr_bh_corr": [0.01, 0.015, 0.02],
        }
    )

    message, counts, fields = _significance_outcome(
        frame,
        module.result_schema.gene_p_value_column,
        module.multiple_testing.gene_methods,
        module,
        "gene",
        "genes",
    )

    assert message == (
        "3 genes tested at the configured p ≤ 0.01 threshold: "
        "2 nominally significant; 0 significant after global Bonferroni "
        "correction; 1 significant after global BH-FDR correction."
    )
    assert counts == {
        "tested": 3,
        "significance_threshold": 0.01,
        "nominal_significant": 2,
        "adjusted_significant": {"bonferroni": 0, "fdr_bh": 1},
    }
    assert fields == [
        ("count", "Genes tested", 3),
        ("analysis", "Reporting threshold", "p ≤ 0.01"),
        ("info", "Nominally significant", 2),
        ("analysis", "Significant after global Bonferroni", 0),
        ("analysis", "Significant after global BH-FDR", 1),
    ]


def test_gene_result_correction_schema_is_configuration_driven(tmp_path):
    config_file = tmp_path / "magma.yaml"
    config_file.write_text(
        "multiple_testing:\n"
        "  gene_methods: [fdr_bh]\n"
        "result_schema:\n"
        "  global_correction_column_pattern: 'ADJUSTED_{method}'\n"
        "  report_delimiter: ','\n",
        encoding="utf-8",
    )
    module = load_module_configuration("magma", config_file)
    genes = tmp_path / "study.genes.out"
    genes.write_text("GENE P\n1 0.01\n2 0.04\n", encoding="utf-8")
    output = tmp_path / "corrected.csv"

    correct_gene_p_values(genes, output, module, RecordingLogger())
    observed = pd.read_csv(output)

    assert observed.columns.tolist() == ["GENE", "P", "ADJUSTED_fdr_bh"]
    np.testing.assert_allclose(observed["ADJUSTED_fdr_bh"], [0.02, 0.04])


@pytest.mark.parametrize("invalid_p", ["NA", "1.01", "-0.01"])
def test_gene_result_correction_rejects_invalid_or_missing_p_values(
    tmp_path, invalid_p,
):
    genes = tmp_path / "study.genes.out"
    genes.write_text("GENE P\n1 %s\n" % invalid_p, encoding="utf-8")

    with pytest.raises(MagmaError, match="invalid values"):
        correct_gene_p_values(
            genes,
            tmp_path / "corrected.tsv",
            load_module_configuration("magma"),
            RecordingLogger(),
        )
    assert not (tmp_path / "corrected.tsv").exists()


def test_gene_set_report_schema_and_corrections_are_configuration_driven(tmp_path):
    config_file = tmp_path / "magma.yaml"
    config_file.write_text(
        "multiple_testing:\n"
        "  global_methods: [bonferroni]\n"
        "  families:\n"
        "    custom:\n"
        "      pattern: '^GO_'\n"
        "      methods: [fdr_bh]\n"
        "result_schema:\n"
        "  global_correction_column_pattern: 'ADJUSTED_{method}'\n"
        "  family_correction_column_pattern: '{family}_ADJUSTED_{method}'\n"
        "  report_dataset_column: study\n"
        "  report_gene_set_description_column: set_description\n"
        "  report_input_genes_column: source_genes\n"
        "  report_common_genes_column: common_ids\n"
        "  report_common_gene_p_values_column: common_p_values\n"
        "  report_total_genes_column: set_size\n"
        "  report_common_gene_count_column: tested_size\n"
        "  report_delimiter: ','\n"
        "  report_null_value: MISSING\n",
        encoding="utf-8",
    )
    module = load_module_configuration("magma", config_file)
    gene_sets = tmp_path / "sets.gmt"
    gene_sets.write_text("GO_SET\tdescription\t1\t3\n", encoding="utf-8")
    raw_gene_sets = tmp_path / "study.gsa.out"
    raw_gene_sets.write_text("VARIABLE P\nGO_SET 0.01\n", encoding="utf-8")
    genes = tmp_path / "study.genes.out"
    genes.write_text("GENE P\n1 0.02\n2 0.20\n", encoding="utf-8")
    logger = RecordingLogger()

    corrected = correct_gene_set_p_values(
        raw_gene_sets, tmp_path / "corrected.csv", module, logger,
    )
    parsed_gene_sets = parse_gene_set_file(gene_sets, module, logger)
    output = annotate_gene_set_results(
        genes,
        parsed_gene_sets,
        corrected,
        tmp_path / "annotated.csv",
        "STUDY",
        module,
        logger,
    )
    observed = pd.read_csv(output)

    assert "ADJUSTED_bonferroni" in observed.columns
    assert "custom_ADJUSTED_fdr_bh" in observed.columns
    assert {"VARIABLE", "FULL_NAME"} <= set(observed.columns)
    assert observed.loc[0, "VARIABLE"] == "GO_SET"
    assert observed.loc[0, "FULL_NAME"] == "GO_SET"
    assert observed.loc[0, "set_description"] == "description"
    assert observed.loc[0, "source_input_genes"] == "1,3"
    assert observed.loc[0, "source_genes"] == "1,3"
    assert not {"V1", "V2", "V3"} & set(observed.columns)
    assert observed.loc[0, "common_ids"] == 1
    assert observed.loc[0, "common_p_values"] == pytest.approx(0.02)
    assert observed.loc[0, "set_size"] == 2
    assert observed.loc[0, "tested_size"] == 1
    assert observed.loc[0, "study"] == "STUDY"


def test_native_magma_gene_sets_have_configured_columns(tmp_path):
    module = load_module_configuration("magma")
    gene_sets = tmp_path / "sets.txt"
    gene_sets.write_text("GOBP_SET\t1 2 3\n", encoding="utf-8")

    observed = parse_gene_set_file(gene_sets, module, RecordingLogger())

    assert observed.columns == [
        "FULL_NAME", "gene_set_description", "source_input_genes", "input_genes",
    ]
    assert observed.to_dicts() == [{
        "FULL_NAME": "GOBP_SET",
        "gene_set_description": None,
        "source_input_genes": "1,2,3",
        "input_genes": "1,2,3",
    }]


def test_gene_set_ids_use_primary_location_ids_without_translation(tmp_path):
    module = load_module_configuration("magma")
    locations = tmp_path / "genes.loc"
    locations.write_text(
        "79501 1 69091 70008 + OR4F5\n"
        "729759 1 367659 368597 + OR4F29\n",
        encoding="utf-8",
    )
    gene_sets = pl.DataFrame(
        {
            "FULL_NAME": ["SET"],
            "gene_set_description": [None],
            "input_genes": ["79501,729759"],
        }
    )

    resolved, reference_ids, summary = resolve_gene_set_identifiers(
        locations, gene_sets, 0.5, module,
    )

    assert resolved.equals(gene_sets)
    assert reference_ids == {"79501", "729759"}
    assert summary["identifier_source"] == "primary_gene_id"
    assert summary["translated_unique_ids"] == 0


def test_gene_set_ids_translate_sixth_location_column_to_primary_ids(tmp_path):
    module = load_module_configuration("magma")
    locations = tmp_path / "genes.loc"
    locations.write_text(
        "ENSG1 1 10 20 + GENE1\n"
        "ENSG2 1 30 40 - GENE2\n"
        "ENSG3 2 50 60 + GENE2\n",
        encoding="utf-8",
    )
    gene_sets = pl.DataFrame(
        {
            "FULL_NAME": ["SET"],
            "gene_set_description": [None],
            "input_genes": ["GENE1,GENE2,UNMAPPED1,UNMAPPED2,UNMAPPED3"],
        }
    )

    resolved, reference_ids, summary = resolve_gene_set_identifiers(
        locations, gene_sets, 0.5, module,
    )

    assert reference_ids == {"ENSG1", "ENSG2", "ENSG3"}
    assert resolved["input_genes"].to_list() == [
        "ENSG1,ENSG2,ENSG3,UNMAPPED1,UNMAPPED2,UNMAPPED3"
    ]
    assert summary["identifier_source"] == "alternate_gene_id"
    assert summary["translated_unique_ids"] == 2
    assert summary["one_to_many_identifiers"] == 1
    assert summary["match_fraction"] == 1.0


@pytest.mark.parametrize(
    ("record", "message"),
    [
        ("ENSG1 1 start 20 + GENE1\n", "integer start and end"),
        ("ENSG1 1 10 20 ? GENE1\n", "invalid gene ID, chromosome"),
        ("ENSG1 1 10 20 + GENE1 EXTRA\n", "must contain five columns"),
    ],
)
def test_gene_location_identifier_resolution_validates_six_column_contract(
    tmp_path, record, message,
):
    module = load_module_configuration("magma")
    locations = tmp_path / "genes.loc"
    locations.write_text(record, encoding="utf-8")
    gene_sets = pl.DataFrame(
        {
            "FULL_NAME": ["SET"],
            "gene_set_description": [None],
            "input_genes": ["GENE1"],
        }
    )

    with pytest.raises(MagmaError, match=message):
        resolve_gene_set_identifiers(locations, gene_sets, 0.5, module)


def test_gene_set_identifier_resolution_fails_when_neither_column_matches(tmp_path):
    module = load_module_configuration("magma")
    locations = tmp_path / "genes.loc"
    locations.write_text("ENSG1 1 10 20 + GENE1\n", encoding="utf-8")
    gene_sets = pl.DataFrame(
        {
            "FULL_NAME": ["SET"],
            "gene_set_description": [None],
            "input_genes": ["ABSENT1,ABSENT2"],
        }
    )

    with pytest.raises(MagmaError, match="neither primary gene-location column 1"):
        resolve_gene_set_identifiers(locations, gene_sets, 0.5, module)


def test_external_annotation_validation_uses_exact_bim_identifier_overlap(tmp_path):
    reference = _reference_files(tmp_path)
    annotation = tmp_path / "external.genes.annot"
    annotation.write_text(
        "ENSG1 1:50:150 rs1 absent .\nENSG2 1:175:225 rs2\n",
        encoding="utf-8",
    )
    module = load_module_configuration("magma")
    module.annotation_validation.minimum_bim_variant_overlap_fraction = 0.60

    observed = validate_gene_annotation(
        annotation, reference, module, RecordingLogger(),
    )

    assert observed["genes"] == 2
    assert observed["matched_annotation_variants"] == 2
    assert observed["bim_overlap_fraction"] == pytest.approx(2 / 3)
    assert observed["excluded_placeholder_assignments"] == 1


def test_nmagma_annotation_union_deduplicates_gene_variant_memberships(tmp_path):
    first = tmp_path / "first.genes.annot"
    second = tmp_path / "second.genes.annot"
    first.write_text("GENE1 1:10:20 rs1 rs2 .\n", encoding="utf-8")
    second.write_text(
        "GENE1 1:10:20 rs2 rs3\nGENE2 1:30:40 rs4\n",
        encoding="utf-8",
    )
    locations = tmp_path / "protein_coding_gene.loc"
    locations.write_text(
        "GENE1 1 10 20 +\nGENE2 1 30 40 -\n",
        encoding="utf-8",
    )
    output = tmp_path / "merged.genes.annot"

    result = merge_gene_annotations(
        [first, second],
        locations,
        output,
        load_module_configuration("magma"),
        RecordingLogger(),
    )

    assert output.read_text(encoding="utf-8").splitlines() == [
        "GENE1\t1:10:20\trs1\trs2\trs3",
        "GENE2\t1:30:40\trs4",
    ]
    assert result["unique_gene_variant_assignments"] == 4
    assert result["excluded_placeholder_assignments"] == 1
    assert result["coordinate_disagreement_records"] == 0


def test_nmagma_annotation_union_uses_canonical_gene_coordinates(tmp_path):
    first = tmp_path / "first.genes.annot"
    second = tmp_path / "second.genes.annot"
    first.write_text("GENE1 1:10:20 rs1\n", encoding="utf-8")
    second.write_text("GENE1 1:11:20 rs2\n", encoding="utf-8")
    locations = tmp_path / "protein_coding_gene.loc"
    locations.write_text("GENE1 1 10 20 +\n", encoding="utf-8")
    output = tmp_path / "merged.genes.annot"

    result = merge_gene_annotations(
        [first, second],
        locations,
        output,
        load_module_configuration("magma"),
        RecordingLogger(),
    )

    assert output.read_text(encoding="utf-8") == "GENE1\t1:10:20\trs1\trs2\n"
    assert result["coordinate_disagreement_records"] == 1


def test_nmagma_annotation_union_excludes_genes_absent_from_location_reference(
    tmp_path,
):
    component = tmp_path / "component.genes.annot"
    component.write_text(
        "GENE1 1:10:20 rs1\nNONCODING NA rs2 rs3\n",
        encoding="utf-8",
    )
    locations = tmp_path / "protein_coding_gene.loc"
    locations.write_text("GENE1 1 10 20 +\n", encoding="utf-8")
    output = tmp_path / "merged.genes.annot"

    result = merge_gene_annotations(
        [component],
        locations,
        output,
        load_module_configuration("magma"),
        RecordingLogger(),
    )

    assert output.read_text(encoding="utf-8") == "GENE1\t1:10:20\trs1\n"
    assert result["genes_absent_from_location_reference"] == 1
    assert result["assignments_excluded_outside_location_reference"] == 2
    assert (
        result["invalid_component_coordinate_records"]
        == 1
    )


def test_nmagma_annotation_union_replaces_invalid_canonical_gene_coordinate(
    tmp_path,
):
    component = tmp_path / "component.genes.annot"
    component.write_text("GENE1 NA rs1\n", encoding="utf-8")
    locations = tmp_path / "protein_coding_gene.loc"
    locations.write_text("GENE1 1 10 20 +\n", encoding="utf-8")

    output = tmp_path / "merged.genes.annot"
    result = merge_gene_annotations(
        [component],
        locations,
        output,
        load_module_configuration("magma"),
        RecordingLogger(),
    )

    assert output.read_text(encoding="utf-8") == "GENE1\t1:10:20\trs1\n"
    assert result["coordinate_disagreement_records"] == 1
    assert result["invalid_component_coordinate_records"] == 1


def test_nmagma_configuration_requires_explicit_positional_windows(tmp_path):
    config = tmp_path / "invalid.yaml"
    config.write_text(
        "mapping:\n"
        "  definitions:\n"
        "    positional:\n"
        "      method: n_magma\n"
        "      gene_location_file: genes.loc\n"
        "      gene_annotation_files: [component.genes.annot]\n",
        encoding="utf-8",
    )

    with pytest.raises(Exception, match="requires explicit annotation windows"):
        load_module_configuration("magma", config)


def test_chrom_magma_uses_deterministic_minimum_element_p_without_correction(
    tmp_path,
):
    mapping = tmp_path / "elements.tsv"
    mapping.write_text(
        "enhancer_b 1 20 30 GENE1 2.0\n"
        "enhancer_a 1 10 15 GENE1 1.0\n"
        "enhancer_c 1 40 50 GENE2 3.0\n",
        encoding="utf-8",
    )
    locations = tmp_path / "elements.loc"
    locations.write_text(
        "enhancer_a 1 10 15\n"
        "enhancer_b 1 20 30\n"
        "enhancer_c 1 40 50\n",
        encoding="utf-8",
    )
    element_results = pd.DataFrame(
        {
            "GENE": ["enhancer_b", "enhancer_a", "enhancer_c"],
            "P": [0.01, 0.01, 0.20],
            "ZSTAT": [2.5, 2.5, 0.8],
        }
    )
    output = tmp_path / "chrom_magma.tsv"

    observed = map_regulatory_elements_to_genes(
        element_results,
        locations,
        mapping,
        output,
        load_module_configuration("magma"),
        RecordingLogger(),
    )

    assert observed["GENE"].to_list() == ["GENE1", "GENE2"]
    assert observed["selected_regulatory_element"].to_list() == [
        "enhancer_a", "enhancer_c",
    ]
    assert observed["mapped_regulatory_element_count"].to_list() == [2, 1]
    assert not any("corr" in column.lower() for column in observed.columns)


def test_chrom_magma_rejects_disagreeing_element_coordinates(tmp_path):
    mapping = tmp_path / "elements.tsv"
    mapping.write_text(
        "enhancer_a 1 10 15 GENE1 1.0\n",
        encoding="utf-8",
    )
    locations = tmp_path / "elements.loc"
    locations.write_text(
        "enhancer_a 1 11 15\n",
        encoding="utf-8",
    )

    with pytest.raises(MagmaError, match="disagree for 1 tested elements"):
        map_regulatory_elements_to_genes(
            pd.DataFrame({"GENE": ["enhancer_a"], "P": [0.01]}),
            locations,
            mapping,
            tmp_path / "chrom_magma.tsv",
            load_module_configuration("magma"),
            RecordingLogger(),
        )


def test_chrom_magma_configuration_rejects_calibrated_gene_p_value_label(tmp_path):
    config = tmp_path / "invalid.yaml"
    config.write_text(
        "mapping:\n"
        "  definitions:\n"
        "    positional:\n"
        "      method: chrom_magma\n"
        "      result_statistic_type: calibrated_gene_p_value\n"
        "      regulatory_element_location_file: elements.loc\n"
        "      element_to_gene_file: elements.tsv\n"
        "      annotation_window_upstream_kb: 0\n"
        "      annotation_window_downstream_kb: 0\n",
        encoding="utf-8",
    )

    with pytest.raises(Exception, match="minimum regulatory-element p-value"):
        load_module_configuration("magma", config)


def test_downstream_modules_reject_chrom_magma_ranking_as_gene_association():
    from postgwas.pipeline.runners import _validated_magma_gene_result

    context = {
        "magma": {
            "primary_mapping": "chrom_magma_ovarian",
            "mapping_analyses": {
                "chrom_magma_ovarian": {
                    "result_statistic_type": "minimum_regulatory_element_p_value",
                }
            },
        }
    }

    with pytest.raises(ValueError, match="requires calibrated MAGMA gene results"):
        _validated_magma_gene_result(context, "PoPS")


def test_batching_obeys_configured_threads_memory_and_minimum_gene_count(tmp_path):
    annotation = tmp_path / "study.genes.annot"
    annotation.write_text("\n".join("gene%d" % value for value in range(4000)), encoding="utf-8")
    configuration = load_configuration(
        cli_overrides={"execution.threads": 8, "execution.memory_gb": 32},
    )
    assert _batch_plan(annotation, configuration) == (4000, 2, 2)


def test_batched_gene_analysis_uses_configured_batch_and_native_prefixes(
    tmp_path, monkeypatch,
):
    from postgwas.modules.magma.analysis import _run_gene_associations

    configuration = load_configuration()
    module = configuration.modules.magma
    annotation = tmp_path / "annotation.genes.annot"
    annotation.write_text("gene1\n", encoding="utf-8")
    paths = {
        "harmonised_p_values": tmp_path / "01_inputs" / "p_values.tsv",
        "gene_batch_prefix": (
            tmp_path / "02_intermediates" / "positional" / "batches" / "study"
        ),
        "gene_prefix": (
            tmp_path
            / "02_intermediates"
            / "positional"
            / "native_outputs"
            / "study"
        ),
        "genes_raw": tmp_path / "native.genes.raw",
        "genes_out": tmp_path / "native.genes.out",
    }
    commands = []

    monkeypatch.setattr(
        "postgwas.modules.magma.analysis._batch_plan",
        lambda *_arguments: (4000, 2, 2),
    )
    monkeypatch.setattr(
        "postgwas.modules.magma.analysis._run_command",
        lambda command, purpose, *_arguments: commands.append((command, purpose)),
    )

    _run_gene_associations(
        "magma",
        "reference",
        annotation,
        paths,
        module.mapping.definitions["positional"],
        configuration,
        RecordingLogger(),
    )

    batch_commands = [
        command
        for command, purpose in commands
        if " batch " in purpose and "batch merge" not in purpose
    ]
    assert len(batch_commands) == 2
    assert {
        command[command.index("--out") + 1] for command in batch_commands
    } == {str(paths["gene_batch_prefix"])}
    merge_command = next(command for command, purpose in commands if "batch merge" in purpose)
    assert merge_command[merge_command.index("--merge") + 1] == str(
        paths["gene_batch_prefix"]
    )
    assert merge_command[merge_command.index("--out") + 1] == str(
        paths["gene_prefix"]
    )


@pytest.mark.parametrize("with_gene_sets", [False, True])
def test_service_publishes_only_a_complete_mocked_run(
    tmp_path, monkeypatch, capsys, with_gene_sets,
):
    from postgwas.modules.magma.service import run_magma_direct

    locations, p_values, reference, genes = _magma_inputs(tmp_path)
    output = tmp_path / "results"
    args = Namespace(
        snp_location_file=str(locations),
        p_value_file=str(p_values),
        magma_ld_reference=str(reference),
        gene_location_file=str(genes),
        magma=sys.executable,
        dataset_id="STUDY",
        output_directory=str(output),
        overwrite=True,
        resolve_variants_to_reference=True,
        threads=1,
        memory_gb=16,
    )
    if with_gene_sets:
        genes.write_text(
            "".join(
                "G%d 1 %d %d +\n" % (value, value * 100, value * 100 + 50)
                for value in range(1, 11)
            ),
            encoding="utf-8",
        )
        gene_sets = tmp_path / "sets.gmt"
        gene_sets.write_text(
            "SET\tdescription\tG1\tG2\tG3\tG4\tG5\tG6\n",
            encoding="utf-8",
        )
        args.gene_set_file = str(gene_sets)

    commands = []

    def fake_command(arguments, purpose, **kwargs):
        commands.append(list(arguments))
        if purpose == "Read MAGMA version":
            return "MAGMA version: v1.10 (custom)"
        for path in kwargs.get("expected_outputs", ()):
            destination = Path(path)
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.name.endswith(".genes.out"):
                rows = (
                    "G1 1 100 150 1 1 1000 3.0 0.001\n"
                    "G7 1 700 750 1 1 1000 3.0 0.001\n"
                    "G8 1 800 850 1 1 1000 3.0 0.001\n"
                    "G9 1 900 950 1 1 1000 3.0 0.001\n"
                    if with_gene_sets
                    else "1 1 50 350 3 1 1000 1.96 0.05\n"
                )
                content = "GENE CHR START STOP NSNPS NPARAM N ZSTAT P\n" + rows
            elif destination.name.endswith(".gsa.out"):
                content = "VARIABLE P\nSET 0.05\n"
            elif destination.name.endswith(".genes.annot"):
                content = (
                    "G1 1:100:150 rs1\n"
                    "G7 1:700:750 rs2\n"
                    "G8 1:800:850 rs3\n"
                    "G9 1:900:950 rs3\n"
                    if with_gene_sets
                    else "1 1:50:350 rs1 rs2 rs3\n"
                )
            else:
                content = "gene1\n"
            destination.write_text(content, encoding="utf-8")
        return ""

    monkeypatch.setattr(
        "postgwas.modules.magma.analysis.run_checked_command", fake_command,
    )
    result = run_magma_direct(args)

    assert result["variant_preparation"]["qc"]["retained_rows"] == 3
    assert Path(result["magma_genes_raw"]).is_file()
    assert Path(result["magma_genes_out"]).is_file()
    corrected = Path(result["magma_genes_corrected"])
    assert corrected.is_file()
    assert corrected.name == "STUDY_magma_genes_corrected.tsv"
    assert pd.read_csv(corrected, sep="\t").columns.tolist()[-2:] == [
        "P_bonferroni_corr",
        "P_fdr_bh_corr",
    ]
    if with_gene_sets:
        assert Path(result["magma_pathway"]).is_file()
        assert result["gene_id_validation"]["overlapping_unique_ids"] == 6
        assert result["gene_id_validation"]["reference_unique_ids"] == 10
        assert result["tested_gene_coverage"]["overlapping_unique_ids"] == 1
        assert result["tested_gene_coverage"]["reference_unique_ids"] == 4
    assert not (output / ".partial" / "STUDY_magma").exists()
    log = output / "05_logs" / "STUDY_magma.log"
    log_text = log.read_text(encoding="utf-8")
    assert "magma_run status=COMPLETED" in log_text
    resolved = yaml.safe_load(
        (output / "00_run_metadata" / "resolved_config.yaml").read_text(
            encoding="utf-8"
        )
    )
    completion = yaml.safe_load(
        (output / "00_run_metadata" / "STUDY_magma_completion.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert completion["status"] == "COMPLETED"
    assert completion["scientific_validation"] == "passed_by_stage"
    assert completion["outputs"]
    assert list(resolved["modules"]) == ["magma"]
    assert resolved["resources"] == {
        "executables": {"magma": sys.executable}
    }
    version_command = next(command for command in commands if "--version" in command)
    assert version_command[-1] == "--version"
    gene_command = next(command for command in commands if "--gene-model" in command)
    assert gene_command[gene_command.index("--seed") + 1] == "10"
    screen = capsys.readouterr().out
    compact_screen = " ".join(screen.split())
    compact_log = " ".join(log_text.split())
    assert "MAGMA analysis progress" in screen
    if with_gene_sets:
        assert "Completed 3/8 · Positional MAGMA · validate competitive gene sets" in screen
        assert "Completed 8/8 · Positional MAGMA · annotate gene-set results" in screen
        assert "All 8 stages completed" in screen
        assert "Gene sets tested : 1" in compact_screen
        assert "Tested gene-set results : 1" in compact_screen
        assert "Reference genes : 10" in compact_screen
        assert "Overlapping identifiers : 6/6 (100.00%)" in compact_screen
        assert "Study-tested genes : 4" in compact_screen
        assert "Gene-set genes tested : 1/6" in compact_screen
    else:
        assert "Completed 4/4 · Positional MAGMA · correct p-values" in screen
        assert "All 4 stages completed" in screen
    assert "Formatter rows : 5" in compact_screen
    assert "Prepared variants : 3" in compact_screen
    assert "BIM intersection : applied" in compact_screen
    assert "Unique BIM ID matches : 3" in compact_screen
    assert "Compatible unique variants : 3" in compact_screen
    assert "Rows absent from BIM : 1" in compact_screen
    assert "Rows with coordinate mismatch : 0" in compact_screen
    assert "Rows with allele-pair mismatch : 0" in compact_screen
    assert "Duplicate rows consolidated : 1 (lowest-p-value rule)" in compact_screen
    assert "Genes tested : %d" % (4 if with_gene_sets else 1) in compact_screen
    assert "Reporting threshold : p ≤ 0.05" in compact_screen
    expected_significant = 4 if with_gene_sets else 1
    assert "Nominally significant : %d" % expected_significant in compact_screen
    assert (
        "Significant after global Bonferroni : %d" % expected_significant
        in compact_screen
    )
    assert (
        "Significant after global BH-FDR : %d" % expected_significant
        in compact_screen
    )
    progress_field_lines = [
        line
        for line in screen.splitlines()
        if any(
            label in line
            for label in (
                "Input rows",
                "Formatter rows",
                "Prepared variants",
                "Genes tested",
                "Reporting threshold",
            )
        )
        and " : " in line
    ]
    assert len(progress_field_lines) >= 4
    assert len({line.index(" : ") for line in progress_field_lines}) == 1
    assert "stage_outcome" in log_text
    assert "reference_unique_compatible_variants=3" in compact_log
    assert "reference_coordinate_mismatch_rows=0" in compact_log
    assert "reference_allele_mismatch_rows=0" in compact_log

    args.overwrite = False
    args.resume = True
    resumed = run_magma_direct(args)
    assert resumed["resumed"] is True
    assert resumed["resume_mode"] == "validated_checkpoint"


def test_existing_output_detection_excludes_external_mapping_annotations(tmp_path):
    from postgwas.config.models.modules.magma import MagmaMappingDefinition
    from postgwas.modules.magma.analysis import resolve_magma_output_paths
    from postgwas.modules.magma.service import _existing_primary_outputs

    module = load_module_configuration("magma")
    positional = module.mapping.definitions["positional"]
    external_annotation = tmp_path / "resources" / "Brain_Cortex.genes.annot"
    external_annotation.parent.mkdir()
    external_annotation.write_text("1 1:10:20 rs1\n", encoding="utf-8")
    external = MagmaMappingDefinition.model_validate({
        **positional.model_dump(),
        "method": "emagma",
        "display_name": "eMAGMA · brain cortex",
        "context": "brain cortex",
        "source_name": "test eMAGMA reference",
        "source_version": "test",
        "source_url": "https://example.org/emagma",
        "gene_annotation_file": str(external_annotation),
    })
    module = module.model_copy(update={
        "mapping": module.mapping.model_copy(update={
            "selected": ["positional", "emagma_brain_cortex"],
            "definitions": {
                "positional": positional,
                "emagma_brain_cortex": external,
            },
        }),
    })
    output = tmp_path / "run"

    assert _existing_primary_outputs(output, "STUDY", module) == []

    generated = resolve_magma_output_paths(
        output, "STUDY", module, "emagma_brain_cortex",
    )["genes_out"]
    generated.parent.mkdir(parents=True)
    generated.write_text("GENE P\n1 0.1\n", encoding="utf-8")

    assert _existing_primary_outputs(output, "STUDY", module) == [generated]


def test_magma_output_layout_separates_inputs_intermediates_and_results(tmp_path):
    from postgwas.modules.magma.analysis import resolve_magma_output_paths

    module = load_module_configuration("magma")
    paths = resolve_magma_output_paths(
        tmp_path, "STUDY", module, "positional",
    )

    assert module.output_layout.completion_manifest.startswith("00_run_metadata/")
    assert paths["harmonised_p_values"].relative_to(tmp_path).parts[0] == "01_inputs"
    assert paths["annotation_prefix"].relative_to(tmp_path).parts[:3] == (
        "02_intermediates", "positional", "annotations",
    )
    assert paths["gene_batch_prefix"].relative_to(tmp_path).parts[:3] == (
        "02_intermediates", "positional", "batches",
    )
    assert paths["genes_out"].relative_to(tmp_path).parts[:3] == (
        "02_intermediates", "positional", "native_outputs",
    )
    assert paths["corrected_genes"].relative_to(tmp_path).parts[:2] == (
        "03_results", "positional",
    )
    assert paths["mapping_comparison"].relative_to(tmp_path).parts[0] == (
        "04_comparisons"
    )


def test_service_failure_keeps_isolated_partial_output_and_finalizes_log(
    tmp_path, monkeypatch,
):
    from postgwas.modules.magma.service import run_magma_direct

    locations, p_values, reference, genes = _magma_inputs(tmp_path)
    output = tmp_path / "failed"
    args = Namespace(
        snp_location_file=str(locations),
        p_value_file=str(p_values),
        magma_ld_reference=str(reference),
        gene_location_file=str(genes),
        magma=sys.executable,
        dataset_id="STUDY",
        output_directory=str(output),
        overwrite=True,
        resolve_variants_to_reference=True,
    )

    def fail_after_version(arguments, purpose, **kwargs):
        if purpose == "Read MAGMA version":
            return "MAGMA version: v1.10 (custom)"
        raise MagmaError("simulated MAGMA failure")

    monkeypatch.setattr(
        "postgwas.modules.magma.analysis.run_checked_command", fail_after_version,
    )
    with pytest.raises(MagmaError, match="simulated MAGMA failure"):
        run_magma_direct(args)

    staging = output / ".partial" / "STUDY_magma"
    assert staging.is_dir()
    assert not list(output.glob("02_intermediates/**/*.genes.out"))
    log = output / "05_logs" / "STUDY_magma.log"
    assert "FAILED" in log.read_text(encoding="utf-8")
    partial_checkpoint = yaml.safe_load(
        (
            output
            / "00_run_metadata"
            / "STUDY_magma_completion.yaml"
        ).read_text(encoding="utf-8")
    )
    assert partial_checkpoint["status"] == "PARTIAL"
    assert partial_checkpoint["scientific_validation"] == "incomplete"


def test_preflight_failure_does_not_create_partial_run_marker(tmp_path):
    from postgwas.modules.magma.service import run_magma_direct

    locations, p_values, _, genes = _magma_inputs(tmp_path)
    output = tmp_path / "preflight_failure"
    args = Namespace(
        snp_location_file=str(locations),
        p_value_file=str(p_values),
        magma_ld_reference=str(tmp_path / "missing_reference"),
        gene_location_file=str(genes),
        magma=sys.executable,
        dataset_id="STUDY",
        output_directory=str(output),
        overwrite=False,
        resume=True,
    )

    with pytest.raises(MagmaError, match="LD reference is incomplete or empty"):
        run_magma_direct(args)

    assert not (output / ".partial" / "STUDY_magma").exists()
    assert "FAILED" in (
        output / "05_logs" / "STUDY_magma.log"
    ).read_text(encoding="utf-8")


def test_empty_partial_tree_is_restarted_without_overwrite(tmp_path, monkeypatch):
    from postgwas.modules.magma.service import run_magma_direct

    locations, p_values, reference, genes = _magma_inputs(tmp_path)
    output = tmp_path / "empty_partial"
    staging = output / ".partial" / "STUDY_magma"
    stale_empty_directory = staging / "unused"
    stale_empty_directory.mkdir(parents=True)
    args = Namespace(
        snp_location_file=str(locations),
        p_value_file=str(p_values),
        magma_ld_reference=str(reference),
        gene_location_file=str(genes),
        magma=sys.executable,
        dataset_id="STUDY",
        output_directory=str(output),
        overwrite=False,
        resume=True,
        resolve_variants_to_reference=True,
    )

    def fail_after_preflight(arguments, purpose, **kwargs):
        if purpose == "Read MAGMA version":
            return "MAGMA version: v1.10 (custom)"
        raise MagmaError("simulated analysis failure after staging recovery")

    monkeypatch.setattr(
        "postgwas.modules.magma.analysis.run_checked_command",
        fail_after_preflight,
    )
    with pytest.raises(MagmaError, match="after staging recovery"):
        run_magma_direct(args)

    assert not stale_empty_directory.exists()
    log_text = (output / "05_logs" / "STUDY_magma.log").read_text(
        encoding="utf-8"
    )
    assert "magma_empty_staging_removed" in log_text
    assert "reason=no_partial_files" in log_text


def test_nonempty_partial_tree_remains_protected_without_overwrite(tmp_path):
    from postgwas.modules.magma.service import run_magma_direct

    output = tmp_path / "nonempty_partial"
    staging = output / ".partial" / "STUDY_magma"
    partial_output = staging / "01_inputs" / "prepared.tsv"
    partial_output.parent.mkdir(parents=True)
    partial_output.write_text("SNP\tP\nrs1\t0.05\n", encoding="utf-8")

    with pytest.raises(MagmaError, match="isolated incomplete MAGMA run"):
        run_magma_direct(
            Namespace(
                dataset_id="STUDY",
                output_directory=str(output),
                overwrite=False,
                resume=True,
            )
        )

    assert partial_output.read_text(encoding="utf-8") == "SNP\tP\nrs1\t0.05\n"


def test_configuration_failure_writes_canonical_log(tmp_path):
    from postgwas.modules.magma.service import run_magma_direct

    output = tmp_path / "invalid"
    args = Namespace(
        dataset_id="STUDY",
        output_directory=str(output),
        gene_model="snp-wise=all",
    )
    with pytest.raises(Exception, match="supported with SNP p-value input"):
        run_magma_direct(args)

    log = output / "05_logs" / "STUDY_magma.log"
    assert log.is_file()
    assert "MAGMA configuration failed" in log.read_text(encoding="utf-8")
