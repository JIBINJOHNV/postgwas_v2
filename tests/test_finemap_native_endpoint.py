"""Protect the unfiltered p-value endpoint and zero-exit native failures."""

import shutil
import subprocess

import pytest

from postgwas.config import load_configuration
from postgwas.modules.fine_mapping.engines.finemap import runner
from postgwas.modules.fine_mapping.engines.finemap.merge_results import (
    InvalidModelProbabilityError,
    parse_cred_header,
)


def configuration(**overrides):
    module = load_configuration().modules.fine_mapping
    settings = module.engines.finemap.model_dump()
    settings["prob_cred_set"] = module.credible_set_coverage
    return {**settings, **overrides}


@pytest.mark.parametrize("threshold", [1.0, 0.5, 1e-8])
def test_native_pvalue_endpoint_preserves_configured_policy(tmp_path, monkeypatch, threshold):
    commands = []

    def run(command, *_args, **_kwargs):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(runner, "_run_logged", run)
    master = tmp_path / "study.master"
    runner.run_finemap_binary(
        master, configuration(pvalue_snps=threshold), finemap_binary="configured-finemap",
        timeout_seconds=10, termination_grace_seconds=1,
    )
    command = commands[0]
    if threshold == 1:
        assert "--pvalue-snps" not in command
        assert "native unfiltered endpoint" in master.with_suffix(".finemap.log").read_text()
    else:
        assert float(command[command.index("--pvalue-snps") + 1]) == threshold


@pytest.mark.parametrize("stream", ["stdout", "stderr"])
def test_native_error_diagnostic_fails_even_with_exit_zero(tmp_path, monkeypatch, stream):
    completed = subprocess.CompletedProcess([], 0, stdout="", stderr="")
    setattr(completed, stream, "Error : invalid native argument!\n")
    monkeypatch.setattr(runner, "_run_logged", lambda *_args, **_kwargs: completed)
    with pytest.raises(RuntimeError, match="FINEMAP reported an error: Error : invalid native argument"):
        runner.run_finemap_binary(
            tmp_path / "study.master", configuration(), finemap_binary="configured-finemap",
            timeout_seconds=10, termination_grace_seconds=1,
        )
    assert not list(tmp_path.glob("*.cred*"))


@pytest.mark.skipif(shutil.which("finemap") is None, reason="native FINEMAP is unavailable")
def test_native_unfiltered_endpoint_retains_zero_effect_snp(tmp_path):
    (tmp_path / "input.z").write_text(
        "rsid chromosome position allele1 allele2 maf beta se\n"
        "rs1 1 100 A G 0.2 0.4 0.1\n"
        "rs2 1 200 C T 0.3 0.2 0.1\n"
        "rs3 1 300 A C 0.4 0.0 0.1\n"
    )
    (tmp_path / "input.ld").write_text("1 0.2 0\n0.2 1 0.1\n0 0.1 1\n")
    master = tmp_path / "study.master"
    master.write_text(
        "z;ld;snp;config;cred;log;n_samples\n"
        "input.z;input.ld;study.snp;study.config;study.cred;study.log;10000\n"
    )
    succeeded, _ = runner.run_finemap_binary(
        master, configuration(n_causal_snps=2), finemap_binary=shutil.which("finemap"),
        timeout_seconds=30, termination_grace_seconds=1,
    )
    assert succeeded
    assert (tmp_path / "study.config").stat().st_size > 0
    credible_files = list(tmp_path.glob("study.cred*"))
    assert credible_files
    assert all(0 <= parse_cred_header(path) <= 1 for path in credible_files)
    rows = (tmp_path / "study.snp").read_text().splitlines()
    assert len(rows) == 4
    assert "rs3" in "\n".join(rows)


@pytest.mark.parametrize("model_size,probability", [(1, 0), (1, 1), (2, 0.5), (3, 1e-8)])
def test_native_model_header_count_and_probability(tmp_path, model_size, probability):
    path = tmp_path / f"study.cred{model_size}"
    path.write_text(f"# Post-Pr(# of causal SNPs is {model_size}) = {probability}\n")
    assert parse_cred_header(path) == probability


@pytest.mark.parametrize("header", [
    "# Post-Pr(# of causal SNPs is 2) = 0.5",
    "# Post-Pr(# of causal SNPs is 0) = 0.5",
    "# Post-Pr(# of causal SNPs is -1) = 0.5",
    "# Post-Pr(# of causal SNPs is 1.0) = 0.5",
    "# Post-Pr(# of causal SNPs is one) = 0.5",
    "# Post-Pr(# of causal SNPs is 1) = NaN",
    "# Post-Pr(# of causal SNPs is 1) = Inf",
    "# Post-Pr(# of causal SNPs is 1) = 1e309",
    "# Post-Pr(# of causal SNPs is 1) = -0.1",
    "# Post-Pr(# of causal SNPs is 1) = 1.1",
    "# Post-Pr(# of causal SNPs is 1) = 0.5 extra",
    "# Post-Pr(# of causal SNPs is 1) = 0.5\n# Post-Pr = 0.5",
    "# Post-Pr(# of causal SNPs is 1) = 0.5\n# Post-Pr(# of causal SNPs is 1) = 0.5",
])
def test_invalid_native_model_headers_are_rejected(tmp_path, header):
    path = tmp_path / "study.cred1"
    path.write_text(header + "\n")
    with pytest.raises(InvalidModelProbabilityError):
        parse_cred_header(path)


def test_native_model_count_requires_a_matching_credible_filename(tmp_path):
    path = tmp_path / "study.txt"
    path.write_text("# Post-Pr(# of causal SNPs is 1) = 0.5\n")
    with pytest.raises(InvalidModelProbabilityError, match="does not match"):
        parse_cred_header(path)
