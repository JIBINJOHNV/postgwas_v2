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
from postgwas.core.contracts import RunContext
from postgwas.core.errors import ConfigurationError
from postgwas.core.preflight import PipelinePreflightEvidence
from postgwas.core.statistics import adjust_p_values
from postgwas.modules.magma.cli import build_parser
from postgwas.modules.magma.errors import MagmaError
from postgwas.modules.magma.reporting import (
    MAGMA_PIPELINE_STAGES,
    _ranked_association_lines,
    magma_variant_input_outcome_fields,
    render_magma_results_screen,
    write_magma_pipeline_html,
)
from postgwas.modules.magma.analysis import (
    MagmaReferencePreflight,
    _annotation_window_argument,
    _batch_plan,
    _prepare_gene_set_plans,
    _read_table,
    _significance_outcome,
    annotate_gene_set_results,
    correct_gene_p_values,
    correct_gene_set_p_values,
    parse_gene_set_file,
    prepare_magma_variant_inputs,
    resolve_gene_set_identifiers,
    validate_magma_gene_set_references,
)
from postgwas.modules.magma.service import (
    MagmaPipelineResources,
    _magma_pipeline_execution_configuration,
    preflight_magma_pipeline,
)
from postgwas.pipeline.registry import REGISTRY
from preflight_support import pipeline_input_vcf_evidence
from postgwas.modules.magma.annotations import (
    map_regulatory_elements_to_genes,
    merge_gene_annotations,
    prepare_scoped_gene_annotation,
    read_gene_locations,
    validate_gene_annotation,
    write_pathway_compatible_gene_locations,
)


class RecordingLogger:
    def __init__(self):
        self.records = []

    def record(self, marker, subject, **values):
        self.records.append((marker, subject, values))

    def warn(self, message, indent=0):
        self.records.append(("WARNING", "message", {"message": message}))


def test_pipeline_magma_preflight_validates_resources_before_formatter(tmp_path):
    prefix = tmp_path / "reference"
    Path(str(prefix) + ".bed").write_bytes(bytes((0x6C, 0x1B, 0x01, 0x00)))
    Path(str(prefix) + ".bim").write_text("1 rs1 0 100 A G\n", encoding="utf-8")
    Path(str(prefix) + ".fam").write_text("F1 I1 0 0 0 -9\n", encoding="utf-8")
    gene_locations = tmp_path / "genes.loc"
    gene_locations.write_text(
        "1 1 10 20 + GENE1\n",
        encoding="utf-8",
    )
    magma = tmp_path / "magma"
    magma.write_text("#!/bin/sh\necho 'MAGMA v1.10'\n", encoding="utf-8")
    magma.chmod(0o755)
    evidence = preflight_magma_pipeline(
        Namespace(
            dataset_id="STUDY",
            genome_build="GRCh37",
            magma_ld_reference=str(prefix),
            gene_location_file=str(gene_locations),
            magma=str(magma),
            _pipeline_requested_modules=("magma",),
        ),
        preflight_evidence=pipeline_input_vcf_evidence(),
    )

    assert REGISTRY.get("magma").preflight.endswith(":preflight_magma_pipeline")
    assert isinstance(evidence, PipelinePreflightEvidence)
    assert isinstance(evidence.resources, MagmaPipelineResources)
    assert evidence.resources.reference.version == "1.10"
    assert evidence.resources.primary_gene_location_metrics["genes"] == 1
    assert evidence.resources.log_events
    assert "formatter-created" in evidence.deferred_checks[0]


def test_pipeline_magma_preflight_rejects_vcf_build_mismatch(tmp_path):
    input_evidence = pipeline_input_vcf_evidence()
    input_evidence["input_vcf"]["harmonised"]["genome_build"] = "GRCh38"

    with pytest.raises(MagmaError, match="declares genome build GRCh38"):
        preflight_magma_pipeline(
            Namespace(
                dataset_id="STUDY",
                genome_build="GRCh37",
                _pipeline_requested_modules=("magma",),
            ),
            preflight_evidence=input_evidence,
        )


def test_pipeline_magma_configuration_keeps_validated_values_and_stage_path(
    tmp_path,
):
    configuration = load_configuration(cli_overrides={
        "run.output_directory": str(tmp_path / "pipeline-root"),
        "modules.magma.input.ld_reference_prefix": str(tmp_path / "reference"),
    })
    resources = MagmaPipelineResources(
        configuration=configuration,
        reference=MagmaReferencePreflight(
            ld_reference_prefix=str(tmp_path / "reference"),
            executable=sys.executable,
            version="1.10",
            gene_set_plans={},
            analysis_scope={},
        ),
        include_gene_sets=True,
        primary_gene_location_metrics=None,
        file_identities=(),
        log_events=(),
    )
    stage = tmp_path / "04_magma"
    snp_locations = tmp_path / "formatter" / "STUDY_magma_snp_loc.tsv"
    p_values = tmp_path / "formatter" / "STUDY_magma_p_values.tsv"
    derived = _magma_pipeline_execution_configuration(
        Namespace(
            output_directory=str(stage),
            snp_location_file=str(snp_locations),
            p_value_file=str(p_values),
        ),
        resources,
    )

    assert derived.run.output_directory == stage.resolve()
    assert configuration.run.output_directory == tmp_path / "pipeline-root"
    assert derived.modules.magma.input.snp_location_file == str(snp_locations)
    assert derived.modules.magma.input.p_value_file == str(p_values)
    assert configuration.modules.magma.input.snp_location_file is None
    assert configuration.modules.magma.input.p_value_file is None
    assert (
        derived.modules.magma.input.ld_reference_prefix
        == configuration.modules.magma.input.ld_reference_prefix
    )


def test_gene_only_dependency_does_not_read_configured_pathway_file(tmp_path):
    missing_pathway = tmp_path / "configured-but-unused.gmt"
    configuration = load_configuration(cli_overrides={
        "modules.magma.input.gene_set_file": str(missing_pathway),
    })
    logger = RecordingLogger()

    plans = validate_magma_gene_set_references(
        configuration, logger, include_gene_sets=False,
    )

    assert plans == {
        "positional": {"status": "not_requested", "reason": None},
    }
    assert not missing_pathway.exists()
    assert logger.records == []


def test_gene_only_dependency_reports_omit_pathway_sections(tmp_path):
    configuration = load_configuration()
    module = configuration.modules.magma
    analysis = {
        "display_name": "Positional MAGMA",
        "result_statistic_interpretation": "Gene association",
        "gene_significance": {
            "tested": 1,
            "nominal_significant": 0,
            "adjusted_significant": {"bonferroni": 0, "fdr_bh": 0},
        },
        "gene_scope": {"input_units": 1, "retained_units": 1, "excluded_units": 0},
        "gene_set_analysis": {"status": "not_requested", "reason": None},
        "magma_gene_results": str(tmp_path / "STUDY_genes.tsv"),
    }
    result = {
        "variant_preparation": {
            "qc": {
                "input_rows": 3,
                "input_unique_variants": 3,
                "reference_unique_id_matches": 3,
                "not_in_reference_rows": 0,
                "reference_intersection_enabled": True,
                "retained_rows": 3,
                "duplicate_rows_removed": 0,
                "excluded_chromosome_rows": 0,
                "excluded_mhc_rows": 0,
            },
        },
        "mapping_analyses": {"positional": analysis},
        "primary_mapping": "positional",
        "magma_mapping_comparison": str(tmp_path / "mapping.tsv"),
        "magma_pipeline_summary_html": str(tmp_path / "report.html"),
        "magma_pipeline_summary_csv": str(tmp_path / "summary.csv"),
        "magma_version": "1.10",
    }
    records = [{
        "step": 1,
        "stage": "Validate gene results",
        "status": "completed",
        "summary": "Gene-only validation passed.",
        "details": "",
        "output": "",
    }]

    screen = render_magma_results_screen(
        result,
        configuration,
        dataset_id="STUDY",
        output_directory=tmp_path,
        log_file=tmp_path / "magma.log",
        label_width=38,
        include_pathway_results=False,
    )
    report = write_magma_pipeline_html(
        records,
        tmp_path / "report.html",
        dataset_id="STUDY",
        schema=module.pipeline_summary_schema,
        result=result,
        configuration=configuration,
        include_pathway_results=False,
    ).read_text(encoding="utf-8")

    assert "pathway" not in screen.lower()
    assert "pathway" not in report.lower()
    assert "MAGMA gene association" in report


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


def test_magma_pipeline_groups_variant_input_preparation_and_reports_provenance(
    tmp_path,
):
    configuration = load_configuration()
    fields = magma_variant_input_outcome_fields(
        configuration,
        {"rows_in": 100},
        {"variant_id_type": "unique"},
        {
            "reference_unique_id_matches": 90,
            "input_unique_variants": 100,
            "not_in_reference_rows": 10,
            "reference_intersection_enabled": False,
            "excluded_chromosomes": ["Y", "MT"],
            "excluded_chromosome_rows": 0,
            "mhc_region": {
                "chromosome": "6", "start": 28477797, "end": 33448354,
            },
            "excluded_mhc_rows": 0,
            "duplicate_policy": "lowest_p",
            "duplicate_groups_detected": 0,
            "duplicate_rows_removed": 0,
            "retained_rows": 100,
        },
        {
            "snp_loc_file": tmp_path / "01_inputs" / "study.locations.tsv",
            "pval_file": tmp_path / "01_inputs" / "study.associations.tsv",
        },
    )
    values = [(field[1], field[2]) for field in fields if len(field) == 3]

    assert len(MAGMA_PIPELINE_STAGES) == 12
    assert MAGMA_PIPELINE_STAGES[5] == "Prepare variant-level inputs for MAGMA"
    assert ("SNP", "CHROM + POS + REF + ALT → CHROM_POS_REF_ALT") in values
    assert ("P", "FORMAT/LP → raw P = 10⁻ᴸᴾ") in values
    assert ("N_COL", "FORMAT/SS → per-variant total N") in values
    assert ("Present in BIM", "90 (90.00%)") in values
    assert ("BIM intersection", "not applied") in values
    assert ("Unmatched variants", "retained in MAGMA input") in values
    assert ("Chromosome policy", "exclude Y and MT") in values
    assert ("MHC SNP policy", "exclude") in values
    assert ("Retention after analysis", "retained") in values


def test_cli_has_no_independent_defaults_or_analysis_imports():
    parser = build_parser()
    for destination in (
        "snp_location_file", "p_value_file", "magma_ld_reference",
        "gene_location_file", "gene_set_file",
        "magma_positional_gene_id_type", "magma_positional_source_name",
        "magma_positional_source_version", "magma_positional_source_url",
        "magma_positional_context",
        "window_upstream", "window_downstream", "gene_model",
        "sample_size_column", "minimum_snp_overlap",
        "minimum_gene_id_overlap", "gene_set_identifier_mismatch",
        "alternate_gene_id_duplicate_policy",
        "gene_location_alternate_id_type",
        "duplicate_policy",
        "resolve_variants_to_reference",
        "magma_mapping", "primary_magma_mapping",
        "mhc_policy", "mhc_chrom", "mhc_start", "mhc_end",
        "exclude_chromosomes",
        "magma_memory_per_worker_gb",
        "run_config", "resume", "overwrite", "dataset_id",
        "output_directory", "threads", "memory_gb", "seed",
    ):
        action = next(item for item in parser._actions if item.dest == destination)
        assert action.default == argparse.SUPPRESS
    assert "magma" not in {action.dest for action in parser._actions}
    duplicate_action = next(
        item for item in parser._actions if item.dest == "duplicate_policy"
    )
    assert tuple(duplicate_action.choices) == ("err", "lowest_p", "remove")
    rendered_help = " ".join(parser.format_help().split())
    assert "--duplicate-policy {err,lowest_p,remove}" in rendered_help
    assert "Default: lowest_p" in rendered_help
    assert "Available options: err, lowest_p, remove" in rendered_help
    assert "--mhc-policy {include,exclude_snps,exclude_genes,exclude_both}" in rendered_help
    assert "Default: exclude_both" in rendered_help
    assert "Default: Y MT" in rendered_help
    assert "--gene-set-identifier-mismatch {skip,error}" in rendered_help
    assert "Default: skip" in rendered_help
    assert "--magma-memory-per-worker-gb GB" in rendered_help
    assert "Default: 16 GB" in rendered_help
    memory_action = next(
        item for item in parser._actions
        if item.dest == "magma_memory_per_worker_gb"
    )
    assert memory_action.type("12.5") == 12.5
    for invalid in ("0", "-1", "inf", "nan"):
        with pytest.raises(argparse.ArgumentTypeError):
            memory_action.type(invalid)
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
    config_file.write_text(
        "gene_window_upstream_kb: 50\n"
        "snp_harmonisation:\n"
        "  duplicate_policy: err\n",
        encoding="utf-8",
    )
    configuration = _resolved_configuration(
        Namespace(
            run_config=str(config_file),
            genome_build="GRCh37",
            window_upstream=25,
            duplicate_policy="remove",
            gene_set_identifier_mismatch="error",
            alternate_gene_id_duplicate_policy="error",
            gene_location_alternate_id_type="ensembl",
            magma_positional_gene_id_type="ensembl",
            magma_positional_source_name="PoPS gene annotation",
            magma_positional_source_version="gene_annot_jun10",
            magma_positional_source_url="https://github.com/FinucaneLab/pops",
            magma_positional_context="Full PoPS GRCh37 gene universe",
            resolve_variants_to_reference=True,
            mhc_policy="exclude_genes",
            mhc_chrom="chr6",
            mhc_start=100,
            mhc_end=200,
            exclude_chromosomes=["MT"],
            magma_memory_per_worker_gb=12,
        )
    )

    assert configuration.modules.magma.genome_build.value == "GRCh37"
    assert configuration.modules.magma.gene_window_upstream_kb == 25
    assert configuration.modules.magma.snp_harmonisation.duplicate_policy == "remove"
    assert configuration.modules.magma.gene_sets.identifier_mismatch_action == "error"
    assert (
        configuration.modules.magma.gene_sets.alternate_id_duplicate_policy
        == "error"
    )
    assert configuration.modules.magma.input.alternate_gene_id_type == "ensembl"
    positional = configuration.modules.magma.mapping.definitions["positional"]
    assert positional.gene_id_type == "ensembl"
    assert positional.source_name == "PoPS gene annotation"
    assert positional.source_version == "gene_annot_jun10"
    assert positional.source_url == "https://github.com/FinucaneLab/pops"
    assert positional.context == "Full PoPS GRCh37 gene universe"
    assert (
        configuration.modules.magma.snp_harmonisation
        .resolve_variants_to_reference
        is True
    )
    assert configuration.modules.magma.mhc.policy == "exclude_genes"
    assert configuration.modules.magma.mhc.region_override.model_dump() == {
        "chromosome": "chr6", "start": 100, "end": 200,
    }
    assert configuration.modules.magma.chromosomes.exclude == ["MT"]
    assert configuration.modules.magma.batching.memory_per_process_gb == 12


def test_duplicate_policy_is_schema_validated_and_defaults_to_lowest_p(tmp_path):
    module = load_module_configuration("magma")
    assert module.snp_harmonisation.duplicate_policy == "lowest_p"

    for policy in ("err", "lowest_p", "remove"):
        configured = tmp_path / (policy + ".yaml")
        configured.write_text(
            "snp_harmonisation:\n  duplicate_policy: %s\n" % policy,
            encoding="utf-8",
        )
        assert (
            load_module_configuration("magma", configured)
            .snp_harmonisation.duplicate_policy
            == policy
        )

    invalid = tmp_path / "invalid.yaml"
    invalid.write_text(
        "snp_harmonisation:\n  duplicate_policy: first\n",
        encoding="utf-8",
    )
    with pytest.raises(Exception, match="duplicate_policy"):
        load_module_configuration("magma", invalid)


def test_default_magma_scope_excludes_mhc_y_and_mt_but_retains_x():
    configuration = load_configuration()
    module = configuration.modules.magma
    assert module.mhc.policy == "exclude_both"
    assert module.chromosomes.exclude == ["Y", "MT"]
    assert "X" not in module.chromosomes.exclude
    assert configuration.resources.genomes["GRCh37"].regions["mhc"].model_dump() == {
        "chromosome": "6",
        "start": 28477797,
        "end": 33448354,
    }


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
    assert configuration.modules.magma.html_report.page_size == 50
    assert "P_fdr_bh_corr" in configuration.modules.magma.html_report.gene_columns
    assert configuration.modules.magma.screen_summary.top_gene_rows == 10
    assert configuration.modules.magma.screen_summary.top_pathway_rows == 10
    assert (
        configuration.modules.magma.screen_summary.p_value_significant_digits
        == 3
    )


def test_magma_html_report_configuration_is_schema_validated(tmp_path):
    invalid = tmp_path / "invalid_magma.yaml"
    invalid.write_text("html_report:\n  page_size: 0\n", encoding="utf-8")

    with pytest.raises(ConfigurationError, match="page_size"):
        load_module_configuration("magma", invalid)


def test_magma_screen_summary_configuration_is_schema_validated(tmp_path):
    invalid = tmp_path / "invalid_magma.yaml"
    invalid.write_text(
        "screen_summary:\n  top_gene_rows: 0\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="top_gene_rows"):
        load_module_configuration("magma", invalid)

    invalid.write_text(
        "screen_summary:\n  p_value_significant_digits: 1\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="p_value_significant_digits"):
        load_module_configuration("magma", invalid)


def test_magma_ranked_screen_values_use_configured_compact_precision():
    rendered = "\n".join(_ranked_association_lines(
        [{
            "name": "GENE1",
            "p_value": 2.5863e-13,
            "adjusted": {
                "bonferroni": 0.0191309,
                "fdr_bh": 0.00871407,
            },
        }],
        heading="Top gene association",
        method_labels={"bonferroni": "Bonferroni", "fdr_bh": "BH-FDR"},
        p_value_significant_digits=3,
        indent=2,
    ))

    assert "P = 2.59e-13" in rendered
    assert "Bonferroni = 0.0191" in rendered
    assert "BH-FDR = 0.00871" in rendered
    assert "adjusted P" not in rendered


def test_bim_overlap_is_assessed_without_default_reference_filtering(tmp_path):
    locations, p_values, reference, _ = _magma_inputs(tmp_path)
    module = load_configuration().modules.magma
    assert module.snp_harmonisation.resolve_variants_to_reference is False
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
        "reference_variant_count": 3,
        "reference_id_match_rows": 4,
        "reference_unique_id_matches": 3,
        "not_in_reference_rows": 1,
        "overlap_fraction": pytest.approx(0.75),
        "duplicate_policy": "lowest_p",
        "duplicate_groups_detected": 1,
        "duplicate_rows_detected": 2,
        "duplicate_rows_removed": 1,
        "duplicate_groups_resolved_by_lowest_p": 1,
        "duplicate_groups_removed": 0,
        "mhc_policy": "include",
        "mhc_region": None,
        "mhc_region_source": None,
        "excluded_chromosomes": [],
        "excluded_chromosome_rows": 0,
        "excluded_chromosome_counts": {},
        "excluded_mhc_rows": 0,
        "excluded_scope_rows": 0,
        "retained_rows": 4,
    }


def test_magma_scope_filters_mhc_y_and_mt_variants_but_retains_x(tmp_path):
    module = load_configuration().modules.magma
    p_values = tmp_path / "p.tsv"
    locations = tmp_path / "loc.tsv"
    p_values.write_text(
        "SNP\tP\tN_COL\n"
        "autosome\t0.1\t100\n"
        "mhc\t0.2\t100\n"
        "x\t0.3\t100\n"
        "y\t0.4\t100\n"
        "mt\t0.5\t100\n",
        encoding="utf-8",
    )
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\n"
        "autosome\t1\t100\tA\tG\n"
        "mhc\t6\t30000000\tA\tC\n"
        "x\tX\t200\tC\tT\n"
        "y\tY\t300\tG\tA\n"
        "mt\tMT\t400\tT\tC\n",
        encoding="utf-8",
    )
    excluded = tmp_path / "excluded.tsv"
    reference = tmp_path / "reference"
    reference.with_suffix(".bim").write_text(
        "1 autosome 0 100 A G\n6 mhc 0 30000000 A C\n"
        "X x 0 200 C T\nY y 0 300 G A\nMT mt 0 400 T C\n",
        encoding="utf-8",
    )
    result = prepare_magma_variant_inputs(
        p_values,
        locations,
        reference,
        tmp_path / "p_out.tsv",
        tmp_path / "loc_out.tsv",
        module,
        RecordingLogger(),
        analysis_scope={
            "mhc_policy": "exclude_both",
            "exclude_mhc_snps": True,
            "exclude_mhc_genes": True,
            "mhc_region": {"chromosome": "6", "start": 28477797, "end": 33448354},
            "mhc_source": "genome_resources",
            "exclude_chromosomes": ["Y", "MT"],
        },
        excluded_output=excluded,
    )

    assert pd.read_csv(tmp_path / "p_out.tsv", sep="\t")["SNP"].tolist() == [
        "autosome", "x",
    ]
    audit = pd.read_csv(excluded, sep="\t")
    assert dict(zip(audit["SNP"], audit["exclusion_reason"])) == {
        "mhc": "mhc_region",
        "y": "excluded_chromosome",
        "mt": "excluded_chromosome",
    }
    assert result["qc"]["excluded_chromosome_counts"] == {"Y": 1, "MT": 1}
    assert result["qc"]["excluded_mhc_rows"] == 1


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
    assert result["qc"]["reference_id_match_rows"] == 2
    assert result["qc"]["reference_unique_id_matches"] == 2
    assert result["qc"]["not_in_reference_rows"] == 0


def test_reference_intersection_is_identifier_based_and_keeps_lowest_duplicate_p(
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
        "not_in_reference_rows": 1,
        "overlap_fraction": pytest.approx(0.75),
        "duplicate_policy": "lowest_p",
        "duplicate_groups_detected": 1,
        "duplicate_rows_detected": 2,
        "duplicate_rows_removed": 1,
        "duplicate_groups_resolved_by_lowest_p": 1,
        "duplicate_groups_removed": 0,
        "mhc_policy": "include",
        "mhc_region": None,
        "mhc_region_source": None,
        "excluded_chromosomes": [],
        "excluded_chromosome_rows": 0,
        "excluded_chromosome_counts": {},
        "excluded_mhc_rows": 0,
        "excluded_scope_rows": 0,
        "retained_rows": 3,
    }
    first_location = location_output.read_text(encoding="utf-8").splitlines()[0]
    assert first_location == "rs1\t1\t100"
    assert not first_location.startswith("SNP")


def test_duplicate_policy_err_stops_before_writing_prepared_inputs(tmp_path):
    locations, p_values, reference, _ = _magma_inputs(tmp_path)
    module = load_module_configuration("magma")
    module.snp_harmonisation.duplicate_policy = "err"
    p_output = tmp_path / "p_out.tsv"
    location_output = tmp_path / "loc_out.tsv"

    with pytest.raises(
        MagmaError,
        match=r"1 duplicated SNP IDs across 2 rows.*rs2",
    ):
        prepare_magma_variant_inputs(
            p_values,
            locations,
            reference,
            p_output,
            location_output,
            module,
            RecordingLogger(),
        )

    assert not p_output.exists()
    assert not location_output.exists()


def test_duplicate_policy_remove_excludes_complete_duplicate_groups(tmp_path):
    locations, p_values, reference, _ = _magma_inputs(tmp_path)
    module = load_module_configuration("magma")
    module.snp_harmonisation.duplicate_policy = "remove"

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
        "rs1", "rs3", "absent",
    ]
    assert result["qc"]["duplicate_policy"] == "remove"
    assert result["qc"]["duplicate_groups_detected"] == 1
    assert result["qc"]["duplicate_rows_detected"] == 2
    assert result["qc"]["duplicate_rows_removed"] == 2
    assert result["qc"]["duplicate_groups_resolved_by_lowest_p"] == 0
    assert result["qc"]["duplicate_groups_removed"] == 1


def test_duplicate_policy_remove_rejects_an_empty_result(tmp_path):
    reference = _reference_files(tmp_path)
    p_values = tmp_path / "p.tsv"
    p_values.write_text(
        "SNP\tP\tN_COL\nrs1\t0.1\t100\nrs1\t0.2\t100\n",
        encoding="utf-8",
    )
    locations = tmp_path / "loc.tsv"
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\nrs1\t1\t100\tA\tG\n",
        encoding="utf-8",
    )
    module = load_module_configuration("magma")
    module.snp_harmonisation.duplicate_policy = "remove"

    with pytest.raises(MagmaError, match="removed every eligible variant"):
        prepare_magma_variant_inputs(
            p_values,
            locations,
            reference,
            tmp_path / "p_out.tsv",
            tmp_path / "loc_out.tsv",
            module,
            RecordingLogger(),
        )


def test_lowest_p_duplicate_ties_preserve_original_input_order(tmp_path):
    reference = _reference_files(tmp_path)
    p_values = tmp_path / "p.tsv"
    p_values.write_text(
        "SNP\tP\tN_COL\n"
        "rs1\t0.1\t100\n"
        "rs1\t0.1\t200\n"
        "rs2\t0.2\t300\n",
        encoding="utf-8",
    )
    locations = tmp_path / "loc.tsv"
    locations.write_text(
        "SNP\tCHR\tBP\tREF\tALT\n"
        "rs1\t1\t100\tA\tG\n"
        "rs2\t1\t200\tC\tT\n",
        encoding="utf-8",
    )

    prepare_magma_variant_inputs(
        p_values,
        locations,
        reference,
        tmp_path / "p_out.tsv",
        tmp_path / "loc_out.tsv",
        load_module_configuration("magma"),
        RecordingLogger(),
    )

    observed = pd.read_csv(tmp_path / "p_out.tsv", sep="\t")
    assert observed.loc[observed["SNP"] == "rs1", "N_COL"].item() == 100


def test_reference_intersection_enforces_configured_unique_overlap(tmp_path):
    locations, p_values, reference, _ = _magma_inputs(tmp_path)
    module = load_configuration().modules.magma
    module.snp_harmonisation.resolve_variants_to_reference = True
    module.snp_harmonisation.minimum_overlap_fraction = 0.8

    with pytest.raises(
        MagmaError,
        match=r"Only 3/4 unique GWAS variant identifiers \(75.00%\).*80.00%",
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


@pytest.mark.parametrize("duplicate_policy", ["err", "lowest_p", "remove"])
def test_duplicate_location_identifier_cannot_hide_conflicting_alleles(
    tmp_path, duplicate_policy,
):
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

    module = load_module_configuration("magma")
    module.snp_harmonisation.duplicate_policy = duplicate_policy

    with pytest.raises(MagmaError, match="conflicting coordinates or allele pairs"):
        prepare_magma_variant_inputs(
            p_values,
            locations,
            reference,
            tmp_path / "p_out.tsv",
            tmp_path / "loc_out.tsv",
            module,
            RecordingLogger(),
        )


def test_reference_intersection_does_not_compare_alleles(tmp_path):
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

    assert pd.read_csv(tmp_path / "p_out.tsv", sep="\t")["SNP"].tolist() == [
        "rs1", "rs2",
    ]
    assert result["qc"]["reference_unique_id_matches"] == 2


def test_reference_intersection_does_not_compare_coordinates(tmp_path):
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

    assert pd.read_csv(tmp_path / "p_out.tsv", sep="\t")["SNP"].tolist() == [
        "rs1", "rs2",
    ]
    assert result["qc"]["reference_unique_id_matches"] == 2


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
    observed = adjust_p_values(values, "sidak")
    np.testing.assert_allclose(observed, expected)
    assert observed[0] > 0
    assert observed[0] == pytest.approx(2e-300, rel=1e-12)


def test_multiple_testing_corrections_handle_boundaries_and_ties():
    boundary_values = np.array([0.0, 1.0])
    for method in ("bonferroni", "sidak", "holm", "fdr_bh"):
        np.testing.assert_allclose(
            adjust_p_values(boundary_values, method),
            [0.0, 1.0],
        )

    tied_values = np.array([0.01, 0.01, 0.5])
    np.testing.assert_allclose(
        adjust_p_values(tied_values, "holm"),
        [0.03, 0.03, 0.5],
    )
    np.testing.assert_allclose(
        adjust_p_values(tied_values, "fdr_bh"),
        [0.015, 0.015, 0.5],
    )


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


def test_gene_results_join_configured_gene_location_metadata(tmp_path):
    module = load_module_configuration("magma")
    genes = tmp_path / "study.genes.out"
    genes.write_text(
        "GENE CHR START STOP NSNPS NPARAM N ZSTAT P\n"
        "ENSG1 1 65 225 3 1 1000 2.0 0.01\n"
        "ENSG2 1 285 440 2 1 1000 0.0 1.0\n",
        encoding="utf-8",
    )
    locations = tmp_path / "genes.loc"
    locations.write_text(
        "ENSG1 1 100 200 + GENE1\n"
        "ENSG2 1 300 400 - GENE2\n",
        encoding="utf-8",
    )
    output = tmp_path / "annotated.tsv"

    correct_gene_p_values(
        genes,
        output,
        module,
        RecordingLogger(),
        gene_location_file=locations,
    )
    observed = pd.read_csv(output, sep="\t")

    assert observed["GENE"].tolist() == ["ENSG1", "ENSG2"]
    assert observed["GENE_REFERENCE_CHR"].tolist() == [1, 1]
    assert observed["GENE_REFERENCE_START"].tolist() == [100, 300]
    assert observed["GENE_REFERENCE_END"].tolist() == [200, 400]
    assert observed["GENE_REFERENCE_STRAND"].tolist() == ["+", "-"]
    assert observed["GENE_SYMBOL"].tolist() == ["GENE1", "GENE2"]
    assert observed.columns.tolist()[-2:] == [
        "P_bonferroni_corr",
        "P_fdr_bh_corr",
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


def test_gene_set_significance_summary_declares_only_the_primary_correction():
    module = load_module_configuration("magma")
    frame = pl.DataFrame(
        {
            "P": [0.001, 0.03],
            "P_bonferroni_corr": [0.002, 0.06],
            "P_sidak_corr": [0.001999, 0.0591],
            "P_holm_corr": [0.002, 0.03],
            "P_fdr_bh_corr": [0.002, 0.06],
        }
    )

    message, counts, fields = _significance_outcome(
        frame,
        module.result_schema.gene_set_p_value_column,
        module.multiple_testing.global_methods,
        module,
        "gene set",
        "gene sets",
        primary_method=module.multiple_testing.primary_method,
    )

    assert "1 significant by the primary global BH-FDR correction" in message
    assert "Bonferroni correction" not in message
    assert counts["primary_correction"] == "fdr_bh"
    assert counts["primary_significant"] == 1
    assert ("analysis", "Primary correction", "Global BH-FDR") in fields
    assert (
        "analysis", "Significant by primary correction", 1,
    ) in fields
    assert not any(
        label.startswith("Significant after global")
        for _kind, label, _value in fields
    )


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
        "  primary_method: bonferroni\n"
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


def test_primary_gene_set_method_must_be_a_configured_global_method(tmp_path):
    config_file = tmp_path / "magma.yaml"
    config_file.write_text(
        "multiple_testing:\n"
        "  primary_method: sidak\n"
        "  global_methods: [bonferroni, fdr_bh]\n",
        encoding="utf-8",
    )

    with pytest.raises(Exception, match="primary_method must occur in global_methods"):
        load_module_configuration("magma", config_file)


def test_magma_table_reader_honours_a_configured_regex_separator(tmp_path):
    source = tmp_path / "mixed.txt"
    source.write_text("A|B;C\n1|2;3\n", encoding="utf-8")
    logger = RecordingLogger()

    observed = _read_table(source, "mixed MAGMA table", r"[|;]", logger)

    assert observed.columns.tolist() == ["A", "B", "C"]
    assert observed.iloc[0].tolist() == [1, 2, 3]
    parser_record = next(
        values
        for marker, subject, values in logger.records
        if marker == "PARAM" and subject == "table_parser"
    )
    assert parser_record["delimiter_pattern"] == r"[|;]"
    assert parser_record["pandas_engine"] == "python"


def test_magma_delimiter_pattern_is_validated_before_analysis(tmp_path):
    config_file = tmp_path / "magma.yaml"
    config_file.write_text(
        "input:\n  table_delimiter_pattern: '['\n",
        encoding="utf-8",
    )

    with pytest.raises(Exception, match="valid regular expression"):
        load_module_configuration("magma", config_file)

    empty_match = tmp_path / "empty_match.yaml"
    empty_match.write_text(
        "result_schema:\n  table_delimiter_pattern: '\\s*'\n",
        encoding="utf-8",
    )
    with pytest.raises(Exception, match="must not match empty text"):
        load_module_configuration("magma", empty_match)


def test_primary_gene_set_result_columns_cannot_collide_with_corrections(tmp_path):
    config_file = tmp_path / "magma.yaml"
    config_file.write_text(
        "result_schema:\n"
        "  primary_adjusted_p_value_column: P_fdr_bh_corr\n",
        encoding="utf-8",
    )

    with pytest.raises(
        Exception, match="configured correction result columns must be unique",
    ):
        load_module_configuration("magma", config_file)


def test_default_gene_set_results_receive_all_four_global_corrections(tmp_path):
    module = load_module_configuration("magma")
    raw_gene_sets = tmp_path / "study.gsa.out"
    raw_gene_sets.write_text(
        "VARIABLE P FULL_NAME\n"
        "SET_A 0.01 GOBP_ALPHA\n"
        "SET_B 0.20 GOBP_BETA\n",
        encoding="utf-8",
    )

    logger = RecordingLogger()

    corrected = correct_gene_set_p_values(
        raw_gene_sets,
        tmp_path / "corrected.tsv",
        module,
        logger,
    )
    observed = pd.read_csv(tmp_path / "corrected.tsv", sep="\t")

    assert module.multiple_testing.primary_method == "fdr_bh"
    assert module.multiple_testing.global_methods == [
        "bonferroni", "sidak", "holm", "fdr_bh",
    ]
    assert {
        "P_bonferroni_corr", "P_sidak_corr", "P_holm_corr", "P_fdr_bh_corr",
    } <= set(corrected.columns)
    correction_columns = [
        column for column in observed.columns if column.endswith("_corr")
    ]
    assert len(correction_columns) == 12
    assert observed.loc[0, correction_columns].notna().sum() == 8
    assert observed["primary_correction_method"].tolist() == ["fdr_bh", "fdr_bh"]
    np.testing.assert_allclose(
        observed["primary_adjusted_p"], observed["P_fdr_bh_corr"],
    )
    assert observed["primary_significant"].tolist() == [True, False]
    primary_record = next(
        values
        for marker, subject, values in logger.records
        if marker == "RESULT" and subject == "primary_gene_set_correction"
    )
    assert primary_record == {
        "family": "all_gene_sets",
        "method": "fdr_bh",
        "method_label": "BH-FDR",
        "threshold": 0.05,
        "tests": 2,
        "significant": 1,
        "source_column": "P_fdr_bh_corr",
        "adjusted_p_value_column": "primary_adjusted_p",
        "significance_column": "primary_significant",
    }


def test_scoped_annotation_and_gene_sets_preserve_sources_and_report_exclusions(
    tmp_path,
):
    module = load_configuration().modules.magma
    source = tmp_path / "source.genes.annot"
    source_text = (
        "AUTO 1:100:200 rs1\n"
        "MHC 6:30000000:30000100 rs2\n"
        "XGENE X:100:200 rs3\n"
        "YGENE Y:100:200 rs4\n"
        "MTGENE MT:100:200 rs5\n"
    )
    source.write_text(source_text, encoding="utf-8")
    scoped = tmp_path / "scoped.genes.annot"
    excluded = tmp_path / "excluded.tsv"
    summary = tmp_path / "summary.tsv"
    result = prepare_scoped_gene_annotation(
        source,
        scoped,
        excluded,
        summary,
        module,
        {
            "mhc_policy": "exclude_both",
            "exclude_mhc_snps": True,
            "exclude_mhc_genes": True,
            "mhc_region": {"chromosome": "6", "start": 28477797, "end": 33448354},
            "mhc_source": "genome_resources",
            "exclude_chromosomes": ["Y", "MT"],
        },
        RecordingLogger(),
    )

    assert source.read_text(encoding="utf-8") == source_text
    assert scoped.read_text(encoding="utf-8").splitlines() == [
        "AUTO 1:100:200 rs1", "XGENE X:100:200 rs3",
    ]
    assert result["excluded_by_reason"] == {
        "excluded_chromosome": 2,
        "mhc_region": 1,
    }


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

    observed_with_metadata, metadata = parse_gene_set_file(
        gene_sets,
        module,
        RecordingLogger(),
        return_metadata=True,
    )
    assert observed_with_metadata.equals(observed)
    assert metadata == {
        "configured_format": "auto",
        "detected_format": "magma",
        "gene_sets": 1,
    }


def test_gmt_parser_keeps_late_descriptions_as_strings(tmp_path):
    """A late description must not depend on Polars schema inference length."""
    module = load_module_configuration("magma")
    gene_sets = tmp_path / "sets.gmt"
    gene_sets.write_text(
        "".join(
            "SET%d\t%s\tGENE%d\n"
            % (index, "" if index < 120 else "late description", index)
            for index in range(121)
        ),
        encoding="utf-8",
    )

    observed = parse_gene_set_file(gene_sets, module, RecordingLogger())

    assert observed.height == 121
    assert observed["gene_set_description"].dtype == pl.String
    assert observed["gene_set_description"][-1] == "late description"


def test_incompatible_gene_sets_can_be_configured_to_fail(tmp_path):
    locations = tmp_path / "genes.loc"
    locations.write_text("ENSG000001 1 100 200 + GENE1\n", encoding="utf-8")
    gene_sets = tmp_path / "sets.gmt"
    gene_sets.write_text("SET\tdescription\tOTHER1\tOTHER2\n", encoding="utf-8")
    config = tmp_path / "magma.yaml"
    config.write_text(
        "input:\n"
        "  gene_location_file: %s\n"
        "gene_sets:\n"
        "  identifier_mismatch_action: error\n" % locations,
        encoding="utf-8",
    )
    module = load_module_configuration("magma", config)

    with pytest.raises(MagmaError, match="complete MAGMA run was stopped"):
        _prepare_gene_set_plans(
            {"positional": gene_sets}, module, RecordingLogger(),
        )


def test_headered_membership_file_is_not_rewritten_for_magma(tmp_path):
    locations = tmp_path / "genes.loc"
    locations.write_text("GENE1 1 100 200 + ALIAS1\n", encoding="utf-8")
    gene_sets = tmp_path / "memberships.tsv"
    gene_sets.write_text(
        "pathway\tgene\nSET1\tGENE1\n",
        encoding="utf-8",
    )
    config = tmp_path / "magma.yaml"
    config.write_text(
        "input:\n"
        "  gene_location_file: %s\n"
        "gene_sets:\n"
        "  input_format: membership\n"
        "  membership_has_header: true\n" % locations,
        encoding="utf-8",
    )
    module = load_module_configuration("magma", config)

    with pytest.raises(
        MagmaError,
        match="cannot contain a header.*does not rewrite pathway inputs",
    ):
        _prepare_gene_set_plans(
            {"positional": gene_sets}, module, RecordingLogger(),
        )


def test_positional_pathway_is_skipped_when_no_alternate_ids_are_available(
    tmp_path,
):
    locations = tmp_path / "genes.loc"
    locations.write_text("ENSG000001 1 100 200 +\n", encoding="utf-8")
    gene_sets = tmp_path / "sets.gmt"
    gene_sets.write_text("SET\tdescription\tGENE1\n", encoding="utf-8")
    config = tmp_path / "magma.yaml"
    config.write_text(
        "input:\n"
        "  gene_location_file: %s\n" % locations,
        encoding="utf-8",
    )
    module = load_module_configuration("magma", config)

    plan = _prepare_gene_set_plans(
        {"positional": gene_sets}, module, RecordingLogger(),
    )["positional"]

    assert plan["status"] == "skipped"
    assert "absent or contains no usable identifiers" in plan["reason"]
    assert "Gene-association analysis will continue" in plan["reason"]


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
            "input_genes": ["79501,729759,OR4F5"],
        }
    )

    resolved, reference_ids, summary = resolve_gene_set_identifiers(
        locations, gene_sets, 0.5, module,
    )

    assert resolved.equals(gene_sets)
    assert reference_ids == {"79501", "729759"}
    assert summary["identifier_source"] == "primary_gene_id"
    assert summary["pathway_identifiers_modified"] is False
    assert summary["location_rows"] == 2
    assert summary["location_primary_unique_ids"] == 2
    assert summary["location_alternate_unique_ids"] == 2
    assert summary["location_ambiguous_alternate_ids"] == 0
    assert summary["direct_primary_matches"] == 2
    assert summary["primary_comparison_unique_ids"] == 2
    assert summary["primary_match_fraction"] == 1.0
    assert summary["alternate_candidate_matches"] == 1
    assert summary["alternate_ids_matched"] == 0
    assert summary["direct_primary_input_fraction"] == pytest.approx(2 / 3)
    assert summary["alternate_input_fraction"] == pytest.approx(1 / 3)
    assert summary["alternate_reference_fraction"] == pytest.approx(1 / 2)
    assert summary["unmatched_input_ids"] == 1
    assert summary["input_resolved_ids"] == 2
    assert summary["input_resolution_fraction"] == pytest.approx(2 / 3)


def test_gene_set_ids_are_not_rewritten_for_non_positional_mapping(tmp_path):
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

    with pytest.raises(MagmaError, match="Pathway identifiers are never rewritten"):
        resolve_gene_set_identifiers(locations, gene_sets, 0.5, module)

    assert gene_sets["input_genes"].to_list() == [
        "GENE1,GENE2,UNMAPPED1,UNMAPPED2,UNMAPPED3"
    ]


def test_positional_gene_sets_remain_unchanged_when_column_six_is_selected(
    tmp_path,
):
    module = load_module_configuration("magma")
    locations = tmp_path / "genes.loc"
    locations.write_text(
        "ENSG1 1 10 20 + GENE1\n"
        "ENSG2 1 30 40 - GENE2\n"
        "ENSG3 2 50 80 + GENE2\n",
        encoding="utf-8",
    )
    gene_sets = pl.DataFrame(
        {
            "FULL_NAME": ["SET"],
            "gene_set_description": [None],
            "input_genes": ["GENE1,GENE2,UNMAPPED"],
        }
    )

    selected, reference_ids, summary = resolve_gene_set_identifiers(
        locations,
        gene_sets,
        0.5,
        module,
        allow_alternate_reference_ids=True,
    )

    assert selected.equals(gene_sets)
    assert reference_ids == {"GENE1", "GENE2"}
    assert summary["identifier_source"] == "alternate_gene_id"
    assert summary["pathway_identifiers_modified"] is False

    destination, metrics = write_pathway_compatible_gene_locations(
        locations,
        tmp_path / "pathway_compatible.loc",
        module,
        RecordingLogger(),
    )
    assert destination.read_text(encoding="utf-8").splitlines() == [
        "GENE1\t1\t10\t20\t+\tENSG1",
        "GENE2\t2\t50\t80\t+\tENSG3",
    ]
    assert metrics == {
        "input_rows": 3,
        "missing_alternate_id_rows": 0,
        "unique_alternate_ids": 2,
        "ambiguous_alternate_ids": 1,
        "duplicate_rows_removed": 1,
        "duplicate_policy": "longest_interval",
        "retained_rows": 2,
    }


def test_pathway_compatible_location_duplicate_policy_can_stop(tmp_path):
    config = tmp_path / "magma.yaml"
    config.write_text(
        "gene_sets:\n"
        "  alternate_id_duplicate_policy: error\n",
        encoding="utf-8",
    )
    module = load_module_configuration("magma", config)
    locations = tmp_path / "genes.loc"
    locations.write_text(
        "ENSG1 1 10 20 + GENE1\n"
        "ENSG2 1 30 40 - GENE1\n",
        encoding="utf-8",
    )

    with pytest.raises(MagmaError, match="assigned to multiple intervals"):
        write_pathway_compatible_gene_locations(
            locations,
            tmp_path / "pathway_compatible.loc",
            module,
            RecordingLogger(),
        )


def test_pathway_compatible_location_from_headered_source_is_headerless(tmp_path):
    config = tmp_path / "magma.yaml"
    config.write_text(
        "input:\n"
        "  gene_location_has_header: true\n",
        encoding="utf-8",
    )
    module = load_module_configuration("magma", config)
    locations = tmp_path / "genes.loc"
    locations.write_text(
        "GENE CHR START END STRAND SYMBOL\n"
        "ENSG1 1 10 20 + GENE1\n"
        "ENSG2 1 30 40 - GENE2\n",
        encoding="utf-8",
    )

    destination, _ = write_pathway_compatible_gene_locations(
        locations,
        tmp_path / "pathway_compatible.loc",
        module,
        RecordingLogger(),
    )

    assert read_gene_locations(
        destination, module, has_header=False,
    ) == {
        "GENE1": ("1", 10, 20, "+", "ENSG1"),
        "GENE2": ("1", 30, 40, "-", "ENSG2"),
    }


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
        resolve_gene_set_identifiers(
            locations,
            gene_sets,
            0.5,
            module,
            allow_alternate_reference_ids=True,
        )


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


def test_annotation_window_argument_is_upstream_then_downstream():
    assert _annotation_window_argument(35, 10) == "window=35,10"
    assert _annotation_window_argument(7, 3) == "window=7,3"


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
        tuple(command[command.index("--batch") + 1:])
        for command in batch_commands
    } == {("1", "2"), ("2", "2")}
    assert {
        command[command.index("--out") + 1] for command in batch_commands
    } == {str(paths["gene_batch_prefix"])}
    common_arguments = [
        command[:command.index("--out")] for command in batch_commands
    ]
    assert common_arguments[0] == common_arguments[1]
    for command in batch_commands:
        assert "duplicate=error" in command
        assert command[command.index("--seed") + 1] == "10"
    merge_command = next(command for command, purpose in commands if "batch merge" in purpose)
    assert merge_command[merge_command.index("--merge") + 1] == str(
        paths["gene_batch_prefix"]
    )
    assert merge_command[merge_command.index("--out") + 1] == str(
        paths["gene_prefix"]
    )


@pytest.mark.parametrize("gene_set_case", ["none", "compatible", "incompatible"])
def test_service_publishes_only_a_complete_mocked_run(
    tmp_path, monkeypatch, capsys, gene_set_case,
):
    from postgwas.modules.magma.service import run_magma_direct

    locations, p_values, reference, genes = _magma_inputs(tmp_path)
    output = tmp_path / "results"
    run_config = tmp_path / "magma.yaml"
    run_config.write_text("html_report:\n  page_size: 2\n", encoding="utf-8")
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
        vcf=str(tmp_path / "study.vcf.gz"),
        run_config=str(run_config),
    )
    context = RunContext({
        "formatter": {
            "magma": {
                "rows_in": 5,
                "rows_out": 5,
                "variant_id_type": "rsid",
                "columns": {
                    "snp_location": ["SNP", "CHR", "BP", "REF", "ALT"],
                    "p_values": ["SNP", "P", "N_COL"],
                },
                "input_vcf_metadata": {
                    "genome_build": "GRCh37",
                    "postgwas_dataset_id": "STUDY",
                },
            }
        }
    })
    with_gene_sets = gene_set_case == "compatible"
    if gene_set_case != "none":
        genes.write_text(
            "".join(
                "G%d 1 %d %d + SYMBOL%d\n"
                % (value, value * 100, value * 100 + 50, value)
                for value in range(1, 11)
            ),
            encoding="utf-8",
        )
        gene_sets = tmp_path / "sets.gmt"
        gene_sets.write_text(
            (
                "SET\tdescription\tSYMBOL1\tSYMBOL2\tSYMBOL3\tSYMBOL4"
                "\tSYMBOL5\tSYMBOL6\n"
                if with_gene_sets
                else "SET\tdescription\tUNMATCHED1\tUNMATCHED2\n"
            ),
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
                if with_gene_sets:
                    rows = (
                        "SYMBOL1 1 100 150 1 1 1000 3.0 0.001\n"
                        "SYMBOL7 1 700 750 1 1 1000 3.0 0.001\n"
                        "SYMBOL8 1 800 850 1 1 1000 3.0 0.001\n"
                        "SYMBOL9 1 900 950 1 1 1000 3.0 0.001\n"
                    )
                elif gene_set_case != "none":
                    rows = (
                        "G1 1 100 150 1 1 1000 3.0 0.001\n"
                        "G7 1 700 750 1 1 1000 3.0 0.001\n"
                        "G8 1 800 850 1 1 1000 3.0 0.001\n"
                        "G9 1 900 950 1 1 1000 3.0 0.001\n"
                    )
                else:
                    rows = "1 1 50 350 3 1 1000 1.96 0.05\n"
                content = "GENE CHR START STOP NSNPS NPARAM N ZSTAT P\n" + rows
            elif destination.name.endswith(".gsa.out"):
                content = "VARIABLE P\nSET 0.05\n"
            elif destination.name.endswith(".genes.annot"):
                if with_gene_sets:
                    content = (
                        "SYMBOL1 1:100:150 rs1\n"
                        "SYMBOL7 1:700:750 rs2\n"
                        "SYMBOL8 1:800:850 rs3\n"
                        "SYMBOL9 1:900:950 rs3\n"
                    )
                elif gene_set_case != "none":
                    content = (
                        "G1 1:100:150 rs1\n"
                        "G7 1:700:750 rs2\n"
                        "G8 1:800:850 rs3\n"
                        "G9 1:900:950 rs3\n"
                    )
                else:
                    content = "1 1:50:350 rs1 rs2 rs3\n"
            else:
                content = "gene1\n"
            destination.write_text(content, encoding="utf-8")
        return ""

    monkeypatch.setattr(
        "postgwas.modules.magma.analysis.run_checked_command", fake_command,
    )
    result = run_magma_direct(args, context)

    assert result["variant_preparation"]["qc"]["retained_rows"] == 3
    assert Path(result["magma_genes_raw"]).is_file()
    assert Path(result["magma_genes_out"]).is_file()
    corrected = Path(result["magma_genes_annotated"])
    assert corrected.is_file()
    assert corrected.name == "STUDY_magma_genes_annotated.tsv"
    assert pd.read_csv(corrected, sep="\t").columns.tolist()[-2:] == [
        "P_bonferroni_corr",
        "P_fdr_bh_corr",
    ]
    if with_gene_sets:
        prepared = result["mapping_analyses"]["positional"]
        derived_locations = Path(
            prepared["pathway_compatible_gene_location"]
        )
        assert derived_locations.is_file()
        assert derived_locations.read_text(encoding="utf-8").splitlines()[0] == (
            "SYMBOL1\t1\t100\t150\t+\tG1"
        )
        assert gene_sets.read_text(encoding="utf-8").startswith(
            "SET\tdescription\tSYMBOL1\tSYMBOL2"
        )
        assert prepared["pathway_file"] == str(gene_sets.resolve())
        gene_set_commands = [
            command for command in commands if "--set-annot" in command
        ]
        assert len(gene_set_commands) == 1
        set_index = gene_set_commands[0].index("--set-annot")
        assert gene_set_commands[0][set_index + 1] == str(gene_sets.resolve())
        assert Path(result["magma_pathway"]).is_file()
        pathway = pd.read_csv(result["magma_pathway"], sep="\t")
        assert {
            "untested_genes",
            "untested_gene_count",
            "primary_correction_method",
            "primary_adjusted_p",
            "primary_significant",
        } <= set(pathway.columns)
        assert pathway["untested_gene_count"].tolist() == [5]
        assert pathway["primary_correction_method"].tolist() == ["fdr_bh"]
        assert pathway["primary_significant"].tolist() == [True]
        assert result["gene_id_validation"]["overlapping_unique_ids"] == 6
        assert result["gene_id_validation"]["reference_unique_ids"] == 10
        assert result["tested_gene_coverage"]["overlapping_unique_ids"] == 1
        assert result["tested_gene_coverage"]["reference_unique_ids"] == 4
        assert result["tested_gene_coverage"][
            "gene_sets_with_exactly_one_reference_gene"
        ] == 1
        assert result["gene_set_analysis"] == {
            "status": "completed",
            "reason": None,
        }
    elif gene_set_case == "incompatible":
        assert result["gene_set_analysis"]["status"] == "skipped"
        assert "below the configured 50.00% minimum" in result[
            "gene_set_analysis"
        ]["reason"]
        assert "magma_pathway" not in result
    assert not (output / ".partial" / "STUDY_magma").exists()
    comparison = pd.read_csv(result["magma_mapping_comparison"], sep="\t")
    assert comparison["mhc_policy"].unique().tolist() == ["exclude_both"]
    assert comparison["mhc_region"].unique().tolist() == ["6:28477797-33448354"]
    assert comparison["excluded_chromosomes"].unique().tolist() == ["Y,MT"]
    expected_gene_set_status = (
        "completed" if with_gene_sets
        else "skipped" if gene_set_case == "incompatible"
        else "not_requested"
    )
    assert comparison["gene_set_analysis_status"].unique().tolist() == [
        expected_gene_set_status
    ]
    log = output / "05_logs" / "STUDY_magma.log"
    log_text = log.read_text(encoding="utf-8")
    assert "magma_run status=COMPLETED" in log_text
    assert "magma_gene_identifier_preflight" in log_text
    assert "random_seed=10" in log_text
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
        "executables": {"magma": sys.executable},
        "genomes": {
            "GRCh37": {
                "regions": {
                    "mhc": {
                        "chromosome": "6",
                        "start": 28477797,
                        "end": 33448354,
                    }
                }
            }
        },
    }
    version_command = next(command for command in commands if "--version" in command)
    assert version_command[-1] == "--version"
    annotation_command = next(command for command in commands if "--annotate" in command)
    assert annotation_command[annotation_command.index("--annotate") + 1] == (
        "window=35,10"
    )
    gene_command = next(command for command in commands if "--gene-model" in command)
    assert gene_command[gene_command.index("--seed") + 1] == "10"
    assert "duplicate=error" in gene_command
    screen = capsys.readouterr().out
    compact_screen = " ".join(screen.split())
    compact_log = " ".join(log_text.split())
    assert "MAGMA analysis progress" in screen
    if with_gene_sets:
        gene_set_command = next(
            command for command in commands if "--gene-results" in command
        )
        assert gene_set_command[gene_set_command.index("--seed") + 1] == "10"
        assert gene_set_command[gene_set_command.index("--gene-results") + 1] == (
            gene_command[gene_command.index("--out") + 1] + ".genes.raw"
        )
        assert ".batch" not in gene_set_command[
            gene_set_command.index("--gene-results") + 1
        ]
        assert (
            "Completed 6/6 · Positional MAGMA · annotate pathway results and "
            "adjust pathway p-values"
            in compact_screen
        )
        assert "All 6 stages completed" in screen
        assert "Gene sets tested : 1" in compact_screen
        assert "Primary correction : Global BH-FDR" in compact_screen
        assert "Significant by primary correction : 1" in compact_screen
        assert "Nominally significant (descriptive) : 1" in compact_screen
        assert (
            "Other global corrections in results : Bonferroni, Šidák, Holm"
            in compact_screen
        )
        assert "Annotated pathway results : 1" in compact_screen
        assert "MAGMA pre-analysis validation" not in compact_screen
        assert "Gene-identifier comparison · Positional MAGMA" not in compact_screen
        assert "Consistency between the two MAGMA input files" not in compact_screen
        assert "Prepared pathway file" not in compact_screen
        assert (
            "Pathway file handling : passed directly to MAGMA without modification"
            in compact_screen
        )
        assert "Pathway records submitted : 1" in compact_screen
        assert "Pathways submitted to MAGMA : 1" in compact_screen
        assert "Pathways without a valid MAGMA result : 0" in compact_screen
        assert "Study-tested genes : 4" in compact_screen
        assert "Gene-set genes tested : 1/6" in compact_screen
    else:
        assert (
            "Completed 4/4 · Positional MAGMA · annotate gene results and adjust "
            "p-values for multiple testing" in compact_screen
        )
        assert "All 4 stages completed" in screen
        if gene_set_case == "incompatible":
            assert "Competitive pathway analysis" in screen
            assert "not performed because gene identifiers were incompatible" in screen
    assert "Formatter rows : 5" in compact_screen
    assert "Prepared variants : 3" in compact_screen
    assert "BIM intersection : applied" in compact_screen
    assert "BIM reference variants : 3" in compact_screen
    assert "GWAS identifiers found in BIM : 3/4 (75.00%)" in compact_screen
    assert "Rows absent from BIM : 1" in compact_screen
    assert "Rows with coordinate mismatch" not in compact_screen
    assert "Rows with allele-pair mismatch" not in compact_screen
    assert "Duplicate policy : lowest_p" in compact_screen
    assert "Duplicate ID groups detected : 1" in compact_screen
    assert "Rows removed after keeping lowest P : 1" in compact_screen
    assert "Genes tested : %d" % (
        4 if gene_set_case != "none" else 1
    ) in compact_screen
    assert "Reporting threshold : p ≤ 0.05" in compact_screen
    expected_significant = 4 if gene_set_case != "none" else 1
    assert "Nominally significant : %d" % expected_significant in compact_screen
    assert (
        "Significant after global Bonferroni : %d" % expected_significant
        in compact_screen
    )
    assert (
        "Significant after global BH-FDR : %d" % expected_significant
        in compact_screen
    )
    progress_screen = screen.split("MAGMA analysis results", 1)[0]
    progress_field_lines = [
        line
        for line in progress_screen.splitlines()
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
    assert "reference_unique_id_matches=3" in compact_log
    assert "reference_coordinate_mismatch_rows" not in compact_log
    assert "reference_allele_mismatch_rows" not in compact_log
    assert "MAGMA pipeline execution summary" not in screen
    assert "MAGMA analysis results" in screen
    assert "Gene association" in screen
    assert "Genes tested : %d" % expected_significant in compact_screen
    assert "Bonferroni significant · adjusted P ≤ 0.05" in compact_screen
    assert "BH-FDR significant · adjusted P ≤ 0.05" in compact_screen
    assert "Top 10 gene associations · ranked by unadjusted P" in compact_screen
    assert "Bonferroni =" in compact_screen
    assert "BH-FDR =" in compact_screen
    assert "Bonferroni-adjusted P =" not in compact_screen
    assert "BH-FDR-adjusted P =" not in compact_screen
    assert "HTML report : STUDY_magma_report.html" in compact_screen
    assert "Results folder : results" in compact_screen
    if with_gene_sets:
        assert "Competitive pathway analysis" in screen
        assert "Pathways supplied : 1" in compact_screen
        assert "Pathways without a tested gene : 0" in compact_screen
        assert "Pathways with only one tested gene : 1" in compact_screen
        assert "Pathways tested : 1" in compact_screen
        assert "Top 10 pathways · ranked by unadjusted P" in compact_screen
        assert "SET · P = 0.05" in compact_screen
        assert "Pathway results : STUDY_magma_gene_sets_annotated.tsv" in compact_screen
    summary_csv = Path(result["magma_pipeline_summary_csv"])
    summary_html = Path(result["magma_pipeline_summary_html"])
    assert summary_csv.is_file()
    assert summary_html.is_file()
    summary = pd.read_csv(summary_csv)
    assert summary["step"].tolist() == list(range(1, 13))
    assert summary["stage"].tolist() == [
        "Validate the input GWAS-VCF",
        "Validate the PLINK LD-reference files and determine the BIM variant-ID format",
        "Validate the gene-location file",
        "Validate the original pathway GMT file",
        "Compare pathway identifiers with gene-location columns and select the "
        "MAGMA gene ID",
        "Prepare variant-level inputs for MAGMA",
        "Create the MAGMA SNP-to-gene annotation",
        "Calculate MAGMA gene-association statistics",
        "Annotate gene results and adjust gene p-values",
        "Run competitive pathway analysis with the original pathway file",
        "Annotate pathway results and adjust pathway p-values",
        "Validate and publish all outputs",
    ]
    assert summary.loc[0, "output"] == args.vcf
    assert (
        "GWAS-VCF structural validation passed; 5 variant records were read."
        == summary.loc[0, "summary"]
    )
    assert "total_variants=5" in summary.loc[0, "details"]
    assert "structural_validation=passed" in summary.loc[0, "details"]
    assert "not independently verifiable" in summary.loc[1, "details"]
    assert "primary_column_1_unique_IDs=" in summary.loc[2, "details"]
    assert "file_structure_validation=passed" in summary.loc[2, "details"]
    if gene_set_case != "none":
        assert "unique_gene_identifiers=" in summary.loc[3, "details"]
        assert "pathway_identifiers_modified=false" in summary.loc[3, "details"]
    assert summary.loc[5, "stage"] == "Prepare variant-level inputs for MAGMA"
    assert "3 variants were retained" in summary.loc[5, "summary"]
    assert "Present_in_BIM=3 (75.00%)" in summary.loc[5, "details"]
    assert "Retention_after_analysis=retained" in summary.loc[5, "details"]
    assert summary.loc[10, "status"] == (
        "completed" if with_gene_sets else
        "skipped" if gene_set_case == "incompatible" else "not requested"
    )
    html = summary_html.read_text(encoding="utf-8")
    assert "MAGMA gene and pathway association" in html
    assert "Results at a glance" in html
    assert "Quality control and analysis coverage" in html
    assert "Most strongly associated genes" in html
    assert "Complete workflow and validation summary" in html
    assert "Variants represented in PLINK BIM" in html
    assert "Annotated genes without valid MAGMA results" in html
    assert "Association does not establish causality" in html
    assert "STUDY_magma_genes_annotated.tsv" in html
    assert 'class="paginated-results"' in html
    assert 'data-page-size="2"' in html
    assert "Search all rows" in html
    assert "Previous" in html
    assert "Next" in html
    assert "All %d gene results are included" % expected_significant in html
    if expected_significant > 2:
        assert "SYMBOL9" in html or "G9" in html
    if with_gene_sets:
        assert "Competitive pathway results" in html
        assert "All 1 pathway results are included" in html
        assert "STUDY_magma_gene_sets_annotated.tsv" in html
        assert "SYMBOL1" in html
        assert "SET" in html
    assert "<th>step</th><th>stage</th><th>status</th>" in html
    assert "Prepare variant-level inputs for MAGMA" in html
    assert "Validate and publish all outputs" in html
    if with_gene_sets:
        assert "primary_correction=fdr_bh" in compact_log
        assert "primary_significant=1" in compact_log
        assert "pandas_engine=c" in compact_log

    args.overwrite = False
    args.resume = True
    resumed = run_magma_direct(args, context)
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
