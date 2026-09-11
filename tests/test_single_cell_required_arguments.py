"""Fail closed on absent method inputs before upstream pipeline resource work."""

from argparse import Namespace
from pathlib import Path
import os
import subprocess
import sys

import pytest
import yaml

from postgwas.config import load_configuration
from postgwas.core.errors import MissingRequiredArgumentsError
from postgwas.core.required_arguments import require_resolved_arguments
from postgwas.modules.single_cell import service
from postgwas.modules.single_cell.cli import build_parser
from postgwas.modules.single_cell.cli import main as single_cell_main
from postgwas.modules.single_cell.methods.registry import METHOD_REGISTRY
from postgwas.pipeline import cli as pipeline_cli


@pytest.mark.parametrize(("tools", "pipeline", "expected"), (
    (["magma_celltype"], False, {"single-cell-covariates", "magma-gene-results-file"}),
    (["magma_celltype"], True, {"single-cell-covariates"}),
    (["scdrs"], False, {"scdrs-h5ad-file", "scdrs-gene-set-file"}),
    (["scdrs"], True, {"scdrs-h5ad-file", "scdrs-gene-id-map"}),
    (["ldsc_celltype"], False, {
        "ldsc-celltype-ldcts-file", "ldsc-celltype-baseline-prefix",
        "ldsc-celltype-weights-prefix", "ldsc-celltype-sumstats-file",
    }),
    (["ldsc_celltype"], True, {
        "ldsc-celltype-ldcts-file", "ldsc-celltype-baseline-prefix",
        "ldsc-celltype-weights-prefix", "ldsc-celltype-merge-alleles-file",
    }),
))
def test_selected_mode_reports_exact_required_inputs(tools, pipeline, expected):
    configuration = load_configuration()
    configuration.modules.single_cell.tools = tools
    requirements = service.single_cell_required_arguments(configuration, pipeline=pipeline)
    assert {item.option for item in requirements} == {"--" + name for name in expected}
    with pytest.raises(MissingRequiredArgumentsError) as captured:
        require_resolved_arguments(requirements)
    for item in requirements:
        assert (
            f"Required argument not provided: {item.option}. Provide {item.option} "
            f"VALUE or set {item.configuration_path} in the run configuration."
        ) in str(captured.value)


@pytest.mark.parametrize("pipeline", (False, True))
def test_scdrs_exact_mapping_does_not_require_crosswalk(pipeline):
    configuration = load_configuration()
    module = configuration.modules.single_cell
    module.tools = ["scdrs"]
    module.scdrs.magma_gene_set.source = "magma"
    module.scdrs.magma_gene_set.target_identifier_type = "entrez"
    module.scdrs.magma_gene_set.mapping_mode = "exact"
    requirements = service.single_cell_required_arguments(configuration, pipeline=pipeline)
    flags = {item.option for item in requirements}
    assert "--scdrs-gene-id-map" not in flags
    assert ("--scdrs-magma-gene-results-file" in flags) is not pipeline
    assert "--scdrs-gene-set-file" not in flags


def test_direct_formatter_ldsc_requires_merge_list():
    configuration = load_configuration()
    configuration.modules.single_cell.tools = ["ldsc_celltype"]
    configuration.modules.single_cell.ldsc_celltype.input.source = "formatter"
    requirements = service.single_cell_required_arguments(configuration)
    assert "--ldsc-celltype-merge-alleles-file" in {item.option for item in requirements}


def test_cli_override_and_yaml_values_satisfy_pipeline_requirements():
    configuration = load_configuration()
    configuration.modules.single_cell.tools = ["scdrs"]
    configuration.modules.single_cell.scdrs.input.gene_identifier_map_file = Path("map.tsv")
    args = Namespace(scdrs_h5ad_file="atlas.h5ad", vcf="study.vcf.gz")
    pipeline_cli._require_pipeline_arguments(args, ("single_cell",), configuration)
    configuration.modules.single_cell.scdrs.input.gene_identifier_map_file = None
    with pytest.raises(MissingRequiredArgumentsError, match="--scdrs-gene-id-map"):
        pipeline_cli._require_pipeline_arguments(args, ("single_cell",), configuration)


def test_pipeline_aggregates_common_and_all_method_missing_inputs():
    configuration = load_configuration()
    configuration.modules.single_cell.tools = ["magma_celltype", "scdrs", "ldsc_celltype"]
    with pytest.raises(MissingRequiredArgumentsError) as captured:
        pipeline_cli._require_pipeline_arguments(Namespace(), ("single_cell",), configuration)
    message = str(captured.value)
    assert message.count("Required argument not provided:") == 8
    assert "--vcf" in message
    for generated in ("--magma-gene-results-file", "--scdrs-gene-set-file",
                      "--scdrs-magma-gene-results-file", "--ldsc-celltype-sumstats-file"):
        assert generated not in message


def test_direct_missing_values_log_to_resolved_yaml_before_any_preflight(tmp_path, monkeypatch):
    output = tmp_path / "yaml-output"
    config = tmp_path / "run.yaml"
    config.write_text(yaml.safe_dump({
        "run": {"output_directory": str(output), "dataset_id": "YAML_STUDY"},
        "modules": {"single_cell": {"tools": ["scdrs", "ldsc_celltype"]}},
    }))
    for method in METHOD_REGISTRY.values():
        monkeypatch.setattr(method, "preflight_direct", lambda *a: pytest.fail("premature preflight"))
    with pytest.raises(MissingRequiredArgumentsError) as captured:
        service.run_single_cell_direct(Namespace(run_config=str(config)))
    log = output / "logs" / "YAML_STUDY_single_cell.log"
    assert log.is_file()
    log_text = log.read_text()
    for flag in ("--scdrs-h5ad-file", "--scdrs-gene-set-file", "--ldsc-celltype-ldcts-file",
                 "--ldsc-celltype-baseline-prefix", "--ldsc-celltype-weights-prefix",
                 "--ldsc-celltype-sumstats-file"):
        assert flag in str(captured.value)
        assert flag in log_text
    assert not (output / "results").exists()
    assert not (output / "engines").exists()


@pytest.mark.parametrize("tool", ("scdrs", "ldsc_celltype"))
def test_public_pipeline_missing_values_stop_before_vcf_and_resource_scans(tmp_path, tool):
    output = tmp_path / tool
    # Pass the parser's existence check, but any record/header scan must fail.
    vcf = tmp_path / "unreadable_as_vcf.vcf.gz"
    vcf.write_bytes(b"not a GWAS VCF\n")
    result = subprocess.run([
        sys.executable, "-B", "-m", "postgwas", "pipeline",
        "--modules", "single_cell", "--tools", tool,
        "--vcf", str(vcf), "--dataset-id", "STUDY",
        "--output-directory", str(output), "--hide-screen",
    ], capture_output=True, text=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    assert result.returncode == 2
    logs = list(output.rglob("*screen*.log"))
    assert logs
    combined = "\n".join(path.read_text() for path in logs)
    flag = "--scdrs-h5ad-file" if tool == "scdrs" else "--ldsc-celltype-ldcts-file"
    assert "Required argument not provided: " + flag in " ".join(combined.split())
    assert "modules.single_cell." in combined
    assert "Started 1/1" not in combined
    assert "Completed 1/" not in combined
    assert not any(path.name.endswith("completion.yaml") for path in output.rglob("*"))


def test_provided_nonexistent_inputs_are_not_reported_as_missing(tmp_path):
    args = Namespace(tools=["scdrs"], output_directory=str(tmp_path / "output"),
                     scdrs_h5ad_file=str(tmp_path / "absent.h5ad"),
                     scdrs_gene_set_file=str(tmp_path / "absent.gs"))
    with pytest.raises(service.SingleCellError) as captured:
        service.run_single_cell_direct(args)
    assert "Required argument not provided" not in str(captured.value)
    assert "absent.h5ad" in str(captured.value)


@pytest.mark.parametrize(("arguments", "required", "optional"), (
    (["--tools", "scdrs"], "scdrs_gene_set_file", "scdrs_magma_gene_results_file"),
    (["--tools", "scdrs", "--scdrs-gene-set-source", "magma"],
     "scdrs_gene_id_map", "scdrs_gene_set_file"),
    (["--tools", "ldsc_celltype"], "ldsc_celltype_sumstats_file", "ldsc_celltype_merge_alleles_file"),
    (["--tools", "ldsc_celltype", "--ldsc-celltype-sumstats-source", "formatter"],
     "ldsc_celltype_merge_alleles_file", "scdrs_h5ad_file"),
))
def test_direct_help_marks_only_requirements_for_selected_mode(arguments, required, optional):
    parser = build_parser(arguments)
    formatter = parser._get_formatter()
    actions = {action.dest: action for action in parser._actions}
    assert "[bold red]Required:" in formatter._get_help_string(actions[required])
    assert "[bold red]Required:" not in formatter._get_help_string(actions[optional])


def test_direct_help_resolves_yaml_selection_and_cli_override(tmp_path):
    config = tmp_path / "settings.yaml"
    config.write_text(yaml.safe_dump({"modules": {"single_cell": {"tools": ["scdrs"]}}}))
    parser = build_parser(["--run-config", str(config)])
    actions = {action.dest: action for action in parser._actions}
    assert getattr(actions["scdrs_h5ad_file"], "_postgwas_required_help", False)
    parser = build_parser(["--run-config", str(config), "--tools", "ldsc_celltype"])
    actions = {action.dest: action for action in parser._actions}
    assert not getattr(actions["scdrs_h5ad_file"], "_postgwas_required_help", False)
    assert getattr(actions["ldsc_celltype_ldcts_file"], "_postgwas_required_help", False)


@pytest.mark.parametrize("tool", ("magma_celltype", "scdrs", "ldsc_celltype"))
def test_pipeline_help_uses_same_required_metadata_with_yaml_tools(tmp_path, monkeypatch, tool):
    config = tmp_path / "settings.yaml"
    config.write_text(yaml.safe_dump({"modules": {"single_cell": {"tools": [tool]}}}))
    captured = []
    monkeypatch.setattr(pipeline_cli, "print_full_pipeline_help", lambda modules, parser: captured.append(parser))
    monkeypatch.setattr(sys, "argv", ["pipeline", "--modules", "single_cell", "--run-config", str(config), "--help"])
    with pytest.raises(SystemExit) as stopped:
        pipeline_cli.main()
    assert stopped.value.code == 0
    parser = captured[0]
    actions = {action.dest: action for action in parser._actions}
    configuration = load_configuration(config)
    for requirement in service.single_cell_required_arguments(configuration, pipeline=True):
        action = actions[requirement.option[2:].replace("-", "_")]
        assert "[bold red]Required:" in parser._get_formatter()._get_help_string(action)
    for generated in ("magma_gene_results_file", "scdrs_gene_set_file",
                      "scdrs_magma_gene_results_file", "ldsc_celltype_sumstats_file"):
        assert actions[generated].help == "==SUPPRESS=="


def test_direct_cli_malformed_yaml_returns_failure_and_public_log(tmp_path):
    config = tmp_path / "malformed.yaml"
    config.write_text("modules: [\n")
    output = tmp_path / "output"
    arguments = ["--tools", "scdrs", "--run-config", str(config),
                 "--output-directory", str(output), "--hide-screen"]
    assert single_cell_main(arguments) == 1
    result = subprocess.run([
        sys.executable, "-B", "-m", "postgwas", "single_cell", *arguments,
    ], capture_output=True, text=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    assert result.returncode == 1
    logs = list(output.rglob("*screen*.log"))
    assert logs and "Invalid YAML" in logs[0].read_text()
    assert not (output / "results").exists()


@pytest.mark.parametrize(("arguments", "code"), ((["--help"], 0), (["--not-an-option"], 2)))
def test_direct_help_and_parser_errors_keep_system_exit(arguments, code):
    with pytest.raises(SystemExit) as stopped:
        single_cell_main(arguments)
    assert stopped.value.code == code
