"""Regression tests for the fine-mapping to FLAMES interchange contract."""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

import postgwas.modules.flames.service as flames_service
from postgwas.config import load_configuration, load_module_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.validation_reporting import FileValidationDisplay
from postgwas.modules.flames.errors import FlamesError
from postgwas.modules.flames.cli import build_parser
from postgwas.modules.flames.service import (
    _build_commands,
    _validate_scores,
    preflight_flames_pipeline,
    validate_fine_mapping_index,
)
from postgwas.pipeline.registry import REGISTRY
from preflight_support import pipeline_input_vcf_evidence


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


def test_handoff_accepts_multiple_credible_sets_for_one_locus(tmp_path):
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

    result = validate_fine_mapping_index(interchange, _module())

    assert len(result["credible_paths"]) == 2
    assert result["index"]["GenomicLocus"].tolist() == ["1", "1"]


def test_handoff_rejects_duplicate_credible_set_files(tmp_path):
    interchange = tmp_path / "flames_input"
    interchange.mkdir()
    credible_set = interchange / "locus.txt"
    credible_set.write_text(
        "index cred1 prob1\n1 1:110:A_G 1.0\n",
        encoding="utf-8",
    )
    pd.DataFrame([
        {
            "Filename": credible_set.name,
            "GenomicLocus": "1",
            "Annotfiles": str(tmp_path / "annotated_1.txt"),
        },
        {
            "Filename": credible_set.name,
            "GenomicLocus": "2",
            "Annotfiles": str(tmp_path / "annotated_2.txt"),
        },
    ]).to_csv(interchange / "indexfile.txt", sep="\t", index=False)

    with pytest.raises(FlamesError, match="duplicate credible-set files"):
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
        configuration, resources, index, tmp_path, "scores", tmp_path / "api.jsonl",
    )

    assert str(index) in annotation
    assert str(index) in scoring
    assert "--filter" in annotation
    assert annotation[annotation.index("--filter") + 1] == "750000"
    assert "--distance" in scoring
    assert "--weight" not in scoring
    assert "--cmd_vep" not in annotation
    assert "--CADD_file" not in annotation
    assert annotation[annotation.index("--annotation-api-log") + 1] == str(tmp_path / "api.jsonl")
    assert "--annotation-api-settings" in annotation


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
        "vep_cache_genome_build",
        "flames_cadd_mode",
        "cadd_genome_build",
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
        annotation_resource_directory="reference/FLAMES/Annotation_data",
    )

    preflight_flames_pipeline(
        args, preflight_evidence=pipeline_input_vcf_evidence(),
    )


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


def _annotation_inputs(tmp_path, module, gene_ids):
    features = tmp_path / "features.txt"
    features.write_text("MAGMA_Z\nPoPS_Score\n", encoding="utf-8")
    annotation = tmp_path / "annotation.tsv"
    pd.DataFrame([
        {
            module.result_schema.gene_column: gene,
            module.result_schema.symbol_column: "GENE%d" % (number + 1),
            "MAGMA_Z": 2.0 if number == 0 else 0.0,
            "PoPS_Score": 0.3 if number == 0 else 0.0,
        }
        for number, gene in enumerate(gene_ids)
    ]).to_csv(annotation, sep=module.input_schema.table_delimiter, index=False)
    resources = {
        "features": features,
        "magma_genes": {gene_ids[0]},
        "pops_genes": {gene_ids[0]},
    }
    return annotation, resources


@pytest.mark.parametrize("missing_gene", ("ENSG000002", "ENSG000002.7"))
@pytest.mark.parametrize("missing_sources", (("magma", "pops"), ("magma",), ("pops",)))
def test_annotation_validation_audits_genes_without_magma_or_pops_evidence(
    tmp_path, missing_gene, missing_sources,
):
    module = _module()
    annotation, resources = _annotation_inputs(
        tmp_path, module, ["ENSG000001", missing_gene],
    )
    for source in ("magma", "pops"):
        if source not in missing_sources:
            resources["%s_genes" % source].add(missing_gene)

    metrics = flames_service._validate_annotations(
        [annotation], resources, module,
    )

    assert metrics["annotated_genes"] == 2
    for source in ("magma", "pops"):
        expected = [missing_gene] if source in missing_sources else []
        assert metrics["gene_rows_without_%s_evidence" % source] == len(expected)
        assert metrics["genes_without_%s_evidence" % source] == expected


@pytest.mark.parametrize(
    "invalid_gene",
    ("", "   ", None, float("nan"), "NOT_AN_ENSEMBL_ID"),
    ids=("blank", "whitespace", "none", "nan", "malformed"),
)
def test_annotation_validation_rejects_invalid_gene_ids(tmp_path, invalid_gene):
    module = _module()
    annotation, resources = _annotation_inputs(
        tmp_path, module, ["ENSG000001", invalid_gene],
    )

    with pytest.raises(FlamesError, match="gene identifiers|non-empty genes"):
        flames_service._validate_annotations([annotation], resources, module)


def test_annotation_validation_rejects_duplicate_valid_gene_ids(tmp_path):
    module = _module()
    annotation, resources = _annotation_inputs(
        tmp_path, module, ["ENSG000001", "ENSG000001"],
    )

    with pytest.raises(FlamesError, match="unique"):
        flames_service._validate_annotations([annotation], resources, module)


def test_annotation_validation_counts_missing_gene_rows_across_loci(tmp_path):
    module = _module()
    annotation, resources = _annotation_inputs(
        tmp_path, module, ["ENSG000001", "ENSG000002.7"],
    )
    second_annotation = tmp_path / "second_locus.tsv"
    second_annotation.write_bytes(annotation.read_bytes())

    metrics = flames_service._validate_annotations(
        [annotation, second_annotation], resources, module,
    )

    assert metrics["annotation_files"] == 2
    assert metrics["annotated_genes"] == 4
    for source in ("magma", "pops"):
        assert metrics["gene_rows_without_%s_evidence" % source] == 2
        assert metrics["genes_without_%s_evidence" % source] == ["ENSG000002.7"]


@pytest.mark.parametrize(
    ("missing_gene", "accepted"),
    (("ENSG000002.99", True), ("ENSG000002.98", False)),
)
def test_annotation_validation_uses_configured_gene_pattern(
    tmp_path, missing_gene, accepted,
):
    config_file = tmp_path / "flames.yaml"
    config_file.write_text(
        "input_schema:\n  ensembl_gene_pattern: '^ENSG[0-9]+[.]99$'\n",
        encoding="utf-8",
    )
    module = load_module_configuration("flames", config_file)
    annotation, resources = _annotation_inputs(
        tmp_path, module, ["ENSG000001.99", missing_gene],
    )

    if not accepted:
        with pytest.raises(FlamesError, match="Ensembl gene identifiers"):
            flames_service._validate_annotations([annotation], resources, module)
        return

    metrics = flames_service._validate_annotations([annotation], resources, module)

    assert metrics["gene_rows_without_magma_evidence"] == 1
    assert metrics["gene_rows_without_pops_evidence"] == 1
    assert metrics["genes_without_magma_evidence"] == [missing_gene]
    assert metrics["genes_without_pops_evidence"] == [missing_gene]


@pytest.mark.parametrize("missing_gene", (None, float("nan")))
def test_annotation_validation_rejects_null_ids_with_permissive_pattern(
    tmp_path, missing_gene,
):
    config_file = tmp_path / "flames.yaml"
    config_file.write_text(
        "input_schema:\n  ensembl_gene_pattern: '.*'\n", encoding="utf-8",
    )
    module = load_module_configuration("flames", config_file)
    annotation, resources = _annotation_inputs(
        tmp_path, module, ["ENSG000001", missing_gene],
    )

    with pytest.raises(FlamesError, match="non-empty"):
        flames_service._validate_annotations([annotation], resources, module)


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


def test_dry_run_validates_both_commands_without_analysis(
    tmp_path, monkeypatch, capsys,
):
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
        "annotation_directories": {"ENSG": tmp_path},
        "annotation_inventory": {"resource": files["annotation_resource.txt"]},
        "annotation_files": {"resource": files["annotation_resource.txt"]},
        "vep_command": None,
        "vep_cache": None,
        "cadd_file": None,
        "cadd_index": None,
        "tabix": None,
        "feature_names": ["validated fixture"],
        "magma": files["magma.txt"],
        "magma_covariate": files["magma_covariate.txt"],
        "pops": files["pops.txt"],
        "magma_genes": {"ENSG000001"},
        "pops_genes": {"ENSG000001"},
        "input_metrics": {
            "credible_sets": 1,
            "credible_set_variants": 1,
            "magma_genes": 1,
            "pops_genes": 1,
            "shared_magma_pops_genes": 1,
            "magma_covariates": 1,
        },
        "interchange": {
            "index_path": files["indexfile.txt"],
            "index": index,
            "credible_paths": [files["credible.txt"]],
            "root": tmp_path,
            "credible_metrics": [
                {"pip_mass": 1.0, "chromosome": 1, "variants": 1}
            ],
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
        lambda configured_module, validated_resources: resources,
    )

    def record_command(command, label, **kwargs):
        commands.append((command, label, kwargs))
        return ""

    monkeypatch.setattr(flames_service, "run_checked_command", record_command)
    run_config = tmp_path / "run.yaml"
    run_config.write_text(
        "config_version: 1\nlogging:\n  show_progress: false\n",
        encoding="utf-8",
    )
    args = argparse.Namespace(
        run_config=str(run_config),
        dataset_id="TEST",
        output_directory=str(tmp_path / "output"),
        dry_run=True,
        credible_sets_directory=str(tmp_path),
        magma_gene_results_file=str(files["magma.txt"]),
        magma_covariate_results_file=str(files["magma_covariate.txt"]),
        pops_scores_file=str(files["pops.txt"]),
        annotation_resource_directory=str(tmp_path),
    )

    result = flames_service.run_flames_direct(args)

    assert result["status"] == "dry_run"
    assert [label for _, label, _ in commands] == [
        "FLAMES annotation",
        "FLAMES scoring",
    ]
    assert all(kwargs["dry_run"] is True for _, _, kwargs in commands)
    assert not list((tmp_path / "output").glob("results/*"))
    screen = " ".join(capsys.readouterr().out.split())
    assert "FLAMES analysis progress" in screen
    assert "All 5 stages completed" in screen
    for label in (
        "Indexed loci / credible sets",
        "MAGMA gene-result file",
        "MAGMAcovar result file",
        "PoPS score file",
        "Required annotation directories",
        "Annotation files inventoried",
        "Configured Python imports",
    ):
        assert label in screen
    assert (
        "live API selected; remote content and availability are not pinned"
        in screen
    )


def test_missing_direct_inputs_are_reported_together_before_runtime(
    tmp_path, monkeypatch, capsys,
):
    monkeypatch.setattr(
        flames_service,
        "_validate_upstream_resources",
        lambda *args, **kwargs: pytest.fail(
            "resource and runtime validation must follow required arguments"
        ),
    )
    args = argparse.Namespace(
        dataset_id="STUDY",
        output_directory=str(tmp_path / "output"),
    )

    with pytest.raises(ConfigurationError) as captured:
        flames_service.run_flames_direct(args)

    message = str(captured.value)
    required = {
        "--credible-sets-directory": "modules.flames.credible_sets_directory",
        "--magma-gene-results-file": "modules.flames.magma_gene_results_file",
        "--magma-covariate-results-file": (
            "modules.flames.magma_covariate_results_file"
        ),
        "--pops-scores-file": "modules.flames.pops_scores_file",
        "--flames-annotation-directory": (
            "modules.flames.annotation_resource_directory"
        ),
    }
    for option, configuration_path in required.items():
        assert (
            "Required argument not provided: %s. Provide %s VALUE or set %s "
            "in the run configuration."
            % (option, option, configuration_path)
        ) in message
    service_log = tmp_path / "output" / "logs" / "STUDY_flames.log"
    assert "Required argument not provided" in service_log.read_text(
        encoding="utf-8"
    )
    screen = " ".join(capsys.readouterr().out.split())
    assert "Failed 1/5" in screen
    assert "All 5 stages completed" not in screen
    assert "100%" not in screen


def test_flames_direct_and_pipeline_help_mark_only_external_requirements():
    direct = " ".join(build_parser().format_help().split())
    for invocation in (
        "--credible-sets-directory PATH",
        "--magma-gene-results-file PATH",
        "--magma-covariate-results-file PATH",
        "--pops-scores-file PATH",
        "--flames-annotation-directory PATH",
    ):
        assert "%s Required:" % invocation in direct

    required = {
        option.dest: option.config_path
        for option in REGISTRY.get("flames").required_options
    }
    assert required["annotation_resource_directory"] == (
        "modules.flames.annotation_resource_directory"
    )
    for generated in (
        "credible_sets_directory",
        "magma_gene_results_file",
        "magma_covariate_results_file",
        "pops_scores_file",
    ):
        assert generated not in required


def test_annotation_preflight_validates_and_inventories_every_resource_file(
    tmp_path, monkeypatch, capsys,
):
    base = _module()
    annotation_root = tmp_path / "annotation"
    for relative in base.upstream.required_annotation_directories:
        directory = annotation_root / relative
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "resource.dat").write_text("reference\n", encoding="utf-8")
    for pattern in base.upstream.required_annotation_file_patterns:
        path = annotation_root / pattern.format(genome_build="GRCH37")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("required reference\n", encoding="utf-8")
    model_directory = tmp_path / "model"
    model_directory.mkdir()
    (model_directory / base.upstream.model_file).write_bytes(b"model")
    (model_directory / base.upstream.feature_file).write_text(
        "feature_one\nfeature_two\n", encoding="utf-8"
    )
    module = base.model_copy(update={
        "annotation_resource_directory": str(annotation_root),
        "model_directory": str(model_directory),
    })
    configuration = SimpleNamespace(
        modules=SimpleNamespace(flames=module),
        resources=SimpleNamespace(
            executables=SimpleNamespace(python="python", tabix="tabix")
        ),
        execution=SimpleNamespace(timeout_seconds=60),
    )
    monkeypatch.setattr(
        flames_service,
        "_validate_runtime",
        lambda configuration, module, logger=None: "/validated/python",
    )

    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        resources = flames_service._validate_upstream_resources(configuration)
        display.flush()

    expected = {
        str(path.resolve().relative_to(annotation_root.resolve()))
        for path in annotation_root.rglob("*")
        if path.is_file()
    }
    assert set(resources["annotation_inventory"]) == expected
    assert resources["feature_names"] == ["feature_one", "feature_two"]
    screen = capsys.readouterr().out
    assert "FLAMES annotation resources — AVAILABLE — availability only" in screen
    assert "Inventoried annotation files" in screen
    assert not any(
        Path(path).name in screen
        for path in resources["annotation_inventory"].values()
    )

    empty = annotation_root / base.upstream.required_annotation_directories[0] / "empty.dat"
    empty.touch()
    with pytest.raises(FlamesError, match="does not exist or is empty"):
        flames_service._validate_upstream_resources(configuration)

    mismatched_module = module.model_copy(update={
        "vep_mode": "local",
        "vep_command": "vep",
        "vep_cache": str(tmp_path),
        "vep_cache_genome_build": type(module.genome_build)("GRCh38"),
    })
    mismatched_configuration = SimpleNamespace(
        modules=SimpleNamespace(flames=mismatched_module),
        resources=configuration.resources,
        execution=configuration.execution,
    )
    empty.unlink()
    with pytest.raises(FlamesError, match="does not match FLAMES genome build"):
        flames_service._validate_upstream_resources(mismatched_configuration)


def test_local_annotation_requirements_use_shared_actionable_errors():
    module = _module().model_copy(update={
        "annotation_resource_directory": "annotation",
        "vep_mode": "local",
        "vep_command": None,
        "vep_cache": None,
        "cadd_mode": "local",
        "cadd_file": None,
    })

    with pytest.raises(ConfigurationError) as captured:
        flames_service._require_flames_arguments(module, pipeline=True)

    message = str(captured.value)
    for option, path in (
        ("--vep-command", "modules.flames.vep_command"),
        ("--vep-cache", "modules.flames.vep_cache"),
        (
            "--vep-cache-genome-build",
            "modules.flames.vep_cache_genome_build",
        ),
        ("--cadd-file", "modules.flames.cadd_file"),
        ("--cadd-genome-build", "modules.flames.cadd_genome_build"),
    ):
        assert (
            "Required argument not provided: %s. Provide %s VALUE or set %s "
            "in the run configuration."
            % (option, option, path)
        ) in message


def test_handoff_rejects_insufficient_cumulative_probability_mass(tmp_path):
    interchange = tmp_path / "flames_input"
    interchange.mkdir()
    credible_set = interchange / "locus.txt"
    credible_set.write_text(
        "index cred1 prob1\n1 1:110:A_G %s\n2 1:120:C_T %s\n"
        % (0.4, 0.4),
        encoding="utf-8",
    )
    pd.DataFrame([{
        "Filename": credible_set.name,
        "GenomicLocus": "1",
        "Annotfiles": str(tmp_path / "annotated.txt"),
    }]).to_csv(interchange / "indexfile.txt", sep="\t", index=False)

    with pytest.raises(FlamesError, match="cumulative PIP"):
        validate_fine_mapping_index(interchange, _module())


def test_handoff_accepts_marginal_pip_sum_above_one(tmp_path):
    interchange = tmp_path / "flames_input"
    interchange.mkdir()
    credible_set = interchange / "locus.txt"
    credible_set.write_text(
        "index cred1 prob1\n1 1:110:A_G 0.7\n2 1:120:C_T 0.4\n",
        encoding="utf-8",
    )
    pd.DataFrame([{
        "Filename": credible_set.name,
        "GenomicLocus": "1",
        "Annotfiles": str(tmp_path / "annotated.txt"),
    }]).to_csv(interchange / "indexfile.txt", sep="\t", index=False)

    result = validate_fine_mapping_index(interchange, _module())

    assert result["credible_metrics"][0]["pip_mass"] == pytest.approx(1.1)
