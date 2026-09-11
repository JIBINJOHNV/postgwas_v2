"""Independent model controls for joint PoPS and K-POPS pipelines."""

import argparse
import sys

import pytest

from postgwas.modules.kpops import service as kpops_service
from postgwas.modules.kpops.cli import build_parser
from postgwas.modules.pops import service as pops_service
from postgwas.pipeline import cli as pipeline_cli


def _pipeline_parser(monkeypatch, modules, width=120):
    captured = {}
    monkeypatch.setenv("COLUMNS", str(width))
    monkeypatch.setattr(
        pipeline_cli, "print_full_pipeline_help",
        lambda _modules, parser, **kwargs: captured.update(parser=parser),
    )
    monkeypatch.setattr(sys, "argv", ["postgwas pipeline", "--modules", *modules, "--help"])
    with pytest.raises(SystemExit) as stopped:
        pipeline_cli.main()
    assert stopped.value.code == 0
    return captured["parser"]


@pytest.mark.parametrize("modules", [("pops", "kpops"), ("kpops", "pops"), ("kpops",)])
@pytest.mark.parametrize("width", [80, 120])
def test_kpops_pipeline_namespaces_conflicting_options(monkeypatch, modules, width):
    parser = _pipeline_parser(monkeypatch, modules, width)
    actions = {action.dest: action for action in parser._actions}
    for field in ("gene_universe_policy", "training_chromosomes", "use_magma_covariates"):
        action = actions["kpops_" + field]
        assert action.default == argparse.SUPPRESS
        assert action.help is not argparse.SUPPRESS
    help_text = " ".join(parser.format_help().split())
    assert "--kpops-gene-universe-policy POLICY" in help_text
    assert "--kpops-training-chromosomes CHROM" in help_text
    assert "--no-kpops-use-magma-covariates" in help_text
    assert "Default: intersect" in help_text
    assert "Default: loco" in help_text
    if "pops" in modules:
        assert "--gene-universe-policy POLICY" in help_text
        assert "Default: strict" in help_text
    else:
        assert "gene_universe_policy" not in actions


def test_joint_pipeline_cli_overrides_do_not_cross_module_boundaries(monkeypatch, tmp_path):
    parser = _pipeline_parser(monkeypatch, ("pops", "kpops"))
    run_config = tmp_path / "joint.yaml"
    run_config.write_text(
        "modules:\n  pops:\n    gene_universe_policy: strict\n"
        "    training_chromosomes: ['1']\n    use_magma_covariates: true\n"
        "  kpops:\n    gene_universe_policy: intersect\n"
        "    training_chromosomes: [loco]\n    use_magma_covariates: true\n",
    )
    args = parser.parse_args([
        "--modules", "pops", "kpops", "--run-config", str(run_config),
        "--gene-universe-policy", "intersect", "--training-chromosomes", "2",
        "--ignore-magma-covariates", "--kpops-gene-universe-policy", "strict",
        "--kpops-training-chromosomes", "all", "--kpops-use-magma-covariates",
    ])
    args._pipeline_requested_modules = ("pops", "kpops")
    pops = pops_service._resolved_configuration(args).modules.pops
    kpops = kpops_service._resolved_configuration(args, pipeline=True).modules.kpops
    assert (pops.gene_universe_policy, pops.training_chromosomes, pops.use_magma_covariates) == ("intersect", ["2"], False)
    assert (kpops.gene_universe_policy, kpops.training_chromosomes, kpops.use_magma_covariates) == ("strict", ["all"], True)

    args = argparse.Namespace(
        run_config=str(run_config), _pipeline_requested_modules=("pops", "kpops"),
        gene_universe_policy="intersect", training_chromosomes=["2"], use_magma_covariates=False,
    )
    kpops = kpops_service._resolved_configuration(args, pipeline=True).modules.kpops
    assert (kpops.gene_universe_policy, kpops.training_chromosomes, kpops.use_magma_covariates) == ("intersect", ["loco"], True)
    args = argparse.Namespace(run_config=str(run_config), _pipeline_requested_modules=("pops", "kpops"))
    pops = pops_service._resolved_configuration(args).modules.pops
    assert (pops.gene_universe_policy, pops.training_chromosomes, pops.use_magma_covariates) == ("strict", ["1"], True)


def test_direct_kpops_flags_and_invalid_pipeline_values(monkeypatch):
    direct = build_parser().parse_args([
        "--gene-universe-policy", "strict", "--training-chromosomes", "all", "--no-use-magma-covariates",
    ])
    module = kpops_service._resolved_configuration(direct).modules.kpops
    assert (module.gene_universe_policy, module.training_chromosomes, module.use_magma_covariates) == ("strict", ["all"], False)
    parser = _pipeline_parser(monkeypatch, ("kpops",))
    for argv in (["--kpops-gene-universe-policy", "invalid"], ["--kpops-training-chromosomes"], ["--gene-universe-policy", "strict"]):
        with pytest.raises(ValueError):
            parser.parse_args(argv)
