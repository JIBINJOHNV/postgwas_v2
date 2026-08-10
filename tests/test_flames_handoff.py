"""Regression tests for the fine-mapping to FLAMES interchange contract."""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

import postgwas.modules.flames.service as flames_service
from postgwas.config import load_module_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.modules.flames.errors import FlamesError
from postgwas.modules.flames.cli import build_parser
from postgwas.modules.flames.service import (
    _build_commands,
    _validate_scores,
    preflight_flames_pipeline,
    validate_fine_mapping_index,
)


def _module():
    return load_module_configuration("flames")


@pytest.mark.parametrize(
    "credible_set_name",
    (
        "STUDY_chr1_100_200_CS_L1.txt",
        "chr1_100_200.cred1_L1.txt",
    ),
)
def test_handoff_consumes_only_files_named_by_the_generated_index(
    tmp_path, credible_set_name
):
    interchange = tmp_path / "flames_input"
    interchange.mkdir()
    credible_set = interchange / credible_set_name
    credible_set.write_text(
        "index\tcred1\tprob1\n"
        "1\t1:110:A_G\t0.7\n"
        "2\t1:120:C_T\t0.3\n",
        encoding="utf-8",
    )
    (interchange / "genomic_loci.tsv").write_text(
        "GenomicLocus\tchr\tstart\tend\n1\t1\t100\t200\n",
        encoding="utf-8",
    )
    (interchange / "engine_FLAMES_manifest.tsv").write_text(
        "Filename\ttarget_coverage\nmetadata\t0.95\n",
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {
                "Filename": credible_set_name,
                "GenomicLocus": "1",
                "Annotfiles": str(tmp_path / "annotated.txt"),
            }
        ]
    ).to_csv(interchange / "indexfile.txt", sep="\t", index=False)

    result = validate_fine_mapping_index(interchange, _module())

    assert result["index_path"] == interchange.resolve() / "indexfile.txt"
    assert result["credible_metrics"][0]["pip_mass"] == pytest.approx(1.0)


def test_handoff_rejects_metadata_accidentally_listed_as_a_credible_set(tmp_path):
    interchange = tmp_path / "flames_input"
    interchange.mkdir()
    metadata = interchange / "genomic_loci.tsv"
    metadata.write_text(
        "GenomicLocus\tchr\tstart\tend\n1\t1\t100\t200\n",
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {
                "Filename": metadata.name,
                "GenomicLocus": "1",
                "Annotfiles": str(tmp_path / "annotated.txt"),
            }
        ]
    ).to_csv(interchange / "indexfile.txt", sep="\t", index=False)

    with pytest.raises(FlamesError, match="missing required columns"):
        validate_fine_mapping_index(interchange, _module())


def test_handoff_rejects_multiple_credible_sets_for_one_locus(tmp_path):
    interchange = tmp_path / "flames_input"
    interchange.mkdir()
    for number in (1, 2):
        (interchange / ("locus_%d.txt" % number)).write_text(
            "index cred1 prob1\n1 1:110:A_G 1.0\n",
            encoding="utf-8",
        )
    pd.DataFrame(
        [
            {
                "Filename": "locus_1.txt",
                "GenomicLocus": "1",
                "Annotfiles": str(tmp_path / "annotated_1.txt"),
            },
            {
                "Filename": "locus_2.txt",
                "GenomicLocus": "1",
                "Annotfiles": str(tmp_path / "annotated_2.txt"),
            },
        ]
    ).to_csv(interchange / "indexfile.txt", sep="\t", index=False)

    with pytest.raises(FlamesError, match="one credible set per locus"):
        validate_fine_mapping_index(interchange, _module())


@pytest.mark.parametrize("probability", ("missing", "-0.1", "1.1", "inf"))
def test_handoff_rejects_invalid_posterior_probabilities(tmp_path, probability):
    interchange = tmp_path / "flames_input"
    interchange.mkdir()
    credible_set = interchange / "locus_CS_1.txt"
    credible_set.write_text(
        "index cred1 prob1\n1 1:110:A_G %s\n" % probability,
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {
                "Filename": str(credible_set),
                "GenomicLocus": "1",
                "Annotfiles": str(tmp_path / "annotated.txt"),
            }
        ]
    ).to_csv(interchange / "indexfile.txt", sep="\t", index=False)

    with pytest.raises(FlamesError, match=r"finite values in \[0, 1\]"):
        validate_fine_mapping_index(interchange, _module())


def test_handoff_requires_the_engine_generated_index(tmp_path):
    interchange = tmp_path / "flames_input"
    interchange.mkdir()
    (interchange / "locus_CS_1.txt").write_text(
        "index cred1 prob1\n1 1:110:A_G 1.0\n",
        encoding="utf-8",
    )

    with pytest.raises(FlamesError, match="does not exist or is empty"):
        validate_fine_mapping_index(interchange, _module())


def test_handoff_rejects_variant_identifiers_flames_cannot_parse(tmp_path):
    interchange = tmp_path / "flames_input"
    interchange.mkdir()
    credible_set = interchange / "locus_CS_1.txt"
    credible_set.write_text(
        "index cred1 prob1\n1 rs123 1.0\n",
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {
                "Filename": credible_set.name,
                "GenomicLocus": "1",
                "Annotfiles": str(tmp_path / "annotated.txt"),
            }
        ]
    ).to_csv(interchange / "indexfile.txt", sep="\t", index=False)

    with pytest.raises(FlamesError, match="identifier contract"):
        validate_fine_mapping_index(interchange, _module())


def test_handoff_column_contract_cannot_diverge_from_fine_mapping(tmp_path):
    config_file = tmp_path / "flames.yaml"
    config_file.write_text(
        "input_schema:\n  credible_set_variant_column: variant\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="credible_set_variant_column"):
        load_module_configuration("flames", config_file)


def test_flames_commands_use_one_generated_index_and_no_shell_defaults(tmp_path):
    configuration = SimpleNamespace(modules=SimpleNamespace(flames=_module()))
    index = tmp_path / "indexfile.txt"
    resources = {
        "python": "/configured/python",
        "script": tmp_path / "FLAMES.py",
        "annotation_directory": tmp_path / "annotations",
        "pops": tmp_path / "pops.tsv",
        "magma": tmp_path / "magma.genes.out",
        "magma_covariate": tmp_path / "magma.gsa.out",
        "model_directory": tmp_path / "model",
        "vep_command": None,
        "vep_cache": None,
        "cadd_file": None,
        "tabix": None,
    }

    annotation, scoring = _build_commands(
        configuration, resources, index, tmp_path, "scores",
    )

    assert str(index) in annotation
    assert str(index) in scoring
    assert "--filter" in annotation
    assert annotation[annotation.index("--filter") + 1] == "750000"
    assert "--distance" in scoring
    assert "--weight" not in scoring
    assert "--cmd_vep" not in annotation
    assert "--CADD_file" not in annotation


def test_flames_cli_has_no_independent_analysis_defaults():
    args = build_parser().parse_args([])

    for name in (
        "credible_sets_directory",
        "magma_gene_results_file",
        "magma_covariate_results_file",
        "pops_scores_file",
        "annotation_resource_directory",
        "flames_genome_build",
        "flames_model_directory",
        "flames_vep_mode",
        "flames_cadd_mode",
        "dry_run",
    ):
        assert not hasattr(args, name)


@pytest.mark.parametrize(
    ("model", "direction"),
    (
        ([], "two-sided"),
        (["condition-hide=Average"], "greater"),
    ),
)
def test_flames_pipeline_preflight_accepts_documented_upstream_models(
    monkeypatch, model, direction,
):
    monkeypatch.setattr(
        flames_service,
        "_validate_upstream_resources",
        lambda configuration: {},
    )
    args = argparse.Namespace(
        covariate_model=model,
        covariate_direction=direction,
    )

    preflight_flames_pipeline(args)


def test_removed_flames_magmacovar_contract_is_rejected_as_unknown_config(tmp_path):
    config_file = tmp_path / "flames.yaml"
    config_file.write_text(
        "magmacovar_contract:\n"
        "  model: [condition-hide=Average]\n"
        "  direction: greater\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="magmacovar_contract"):
        load_module_configuration("flames", config_file)


def test_published_model_window_cannot_be_overridden(tmp_path):
    config_file = tmp_path / "flames.yaml"
    config_file.write_text("locus_window_bp: 1000000\n", encoding="utf-8")

    with pytest.raises(ConfigurationError, match="locus_window_bp"):
        load_module_configuration("flames", config_file)


def test_score_validation_accepts_no_prioritized_genes(tmp_path):
    module = _module()
    schema = module.result_schema
    raw_path = tmp_path / "scores.raw"
    prediction_path = tmp_path / "scores.pred"
    raw = pd.DataFrame(
        [
            {
                schema.filename_column: "annotation.tsv",
                schema.locus_column: "1",
                schema.symbol_column: "GENE1",
                schema.gene_column: "ENSG000001",
                schema.xgboost_score_column: 0.2,
                schema.pops_score_column: -0.1,
                schema.raw_score_column: 0.17,
                schema.scaled_score_column: 1.0,
                schema.precision_column: 0.4,
                schema.highest_column: 1,
                schema.causal_column: 0,
            }
        ]
    )
    raw.to_csv(raw_path, sep=module.input_schema.table_delimiter, index=False)
    pd.DataFrame(
        columns=[
            schema.locus_column,
            schema.filename_column,
            schema.symbol_column,
            schema.gene_column,
            schema.scaled_score_column,
            schema.raw_score_column,
            schema.precision_column,
        ]
    ).to_csv(
        prediction_path,
        sep=module.input_schema.table_delimiter,
        index=False,
    )

    metrics = _validate_scores(raw_path, prediction_path, module)

    assert metrics["prioritized_genes"] == 0


def test_score_validation_rejects_incorrect_raw_score_normalization(tmp_path):
    module = _module()
    schema = module.result_schema
    rows = []
    for gene, raw_score in (("ENSG000001", 0.2), ("ENSG000002", 0.8)):
        rows.append(
            {
                schema.filename_column: "annotation.tsv",
                schema.locus_column: "1",
                schema.symbol_column: gene,
                schema.gene_column: gene,
                schema.xgboost_score_column: raw_score,
                schema.pops_score_column: 0.0,
                schema.raw_score_column: raw_score,
                schema.scaled_score_column: 0.5,
                schema.precision_column: 0.4,
                schema.highest_column: int(raw_score == 0.8),
                schema.causal_column: 0,
            }
        )
    raw_path = tmp_path / "scores.raw"
    prediction_path = tmp_path / "scores.pred"
    pd.DataFrame(rows).to_csv(
        raw_path, sep=module.input_schema.table_delimiter, index=False
    )
    pd.DataFrame(
        columns=[
            schema.locus_column,
            schema.filename_column,
            schema.symbol_column,
            schema.gene_column,
            schema.scaled_score_column,
            schema.raw_score_column,
            schema.precision_column,
        ]
    ).to_csv(
        prediction_path,
        sep=module.input_schema.table_delimiter,
        index=False,
    )

    with pytest.raises(FlamesError, match="divided by the locus total"):
        _validate_scores(raw_path, prediction_path, module)


def test_dry_run_validates_both_commands_without_analysis(tmp_path, monkeypatch):
    module = _module()
    files = {}
    for name in (
        "credible.txt",
        "indexfile.txt",
        "magma.txt",
        "magma_covariate.txt",
        "pops.txt",
        "FLAMES.py",
        "model.sav",
        "features.txt",
        "annotation_resource.txt",
    ):
        path = tmp_path / name
        path.write_text("validated fixture\n", encoding="utf-8")
        files[name] = path
    index = pd.DataFrame(
        [{"Filename": "credible.txt", "GenomicLocus": "1", "Annotfiles": "unused"}]
    )
    resources = {
        "python": "/configured/python",
        "script": files["FLAMES.py"],
        "model_directory": tmp_path,
        "model": files["model.sav"],
        "features": files["features.txt"],
        "annotation_directory": tmp_path,
        "annotation_files": {"resource": files["annotation_resource.txt"]},
        "vep_command": None,
        "vep_cache": None,
        "cadd_file": None,
        "cadd_index": None,
        "tabix": None,
        "magma": files["magma.txt"],
        "magma_covariate": files["magma_covariate.txt"],
        "pops": files["pops.txt"],
        "magma_genes": {"ENSG000001"},
        "pops_genes": {"ENSG000001"},
        "input_metrics": {"credible_sets": 1},
        "interchange": {
            "index_path": files["indexfile.txt"],
            "index": index,
            "credible_paths": [files["credible.txt"]],
            "credible_metrics": [{"pip_mass": 1.0}],
        },
    }
    commands = []

    monkeypatch.setattr(
        flames_service,
        "_validate_upstream_resources",
        lambda configuration, logger=None: resources,
    )
    monkeypatch.setattr(
        flames_service,
        "_validate_scientific_inputs",
        lambda configured_module, validated_resources: validated_resources,
    )

    def record_command(command, label, **kwargs):
        commands.append((command, label, kwargs))
        return ""

    monkeypatch.setattr(flames_service, "run_checked_command", record_command)
    args = argparse.Namespace(
        dataset_id="TEST",
        output_directory=str(tmp_path / "output"),
        dry_run=True,
    )

    result = flames_service.run_flames_direct(args)

    assert result["status"] == "dry_run"
    assert [label for _, label, _ in commands] == [
        "FLAMES annotation",
        "FLAMES scoring",
    ]
    assert all(kwargs["dry_run"] is True for _, _, kwargs in commands)
    assert not list((tmp_path / "output").glob("results/*"))


@pytest.mark.parametrize("probabilities", ((0.4, 0.4), (0.7, 0.4)))
def test_handoff_rejects_invalid_cumulative_probability_mass(tmp_path, probabilities):
    interchange = tmp_path / "flames_input"
    interchange.mkdir()
    credible_set = interchange / "locus.txt"
    credible_set.write_text(
        "index cred1 prob1\n1 1:110:A_G %s\n2 1:120:C_T %s\n"
        % probabilities,
        encoding="utf-8",
    )
    pd.DataFrame([{
        "Filename": credible_set.name,
        "GenomicLocus": "1",
        "Annotfiles": str(tmp_path / "annotated.txt"),
    }]).to_csv(interchange / "indexfile.txt", sep="\t", index=False)

    with pytest.raises(FlamesError, match="cumulative PIP"):
        validate_fine_mapping_index(interchange, _module())
