"""Focused scientific and execution regressions for MAGMAcovar."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pytest

from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.modules.magmacovar.errors import MagmaCovarError
from postgwas.modules.magmacovar.main import (
    build_magma_covariate_command,
    validate_magma_covariate_inputs,
    validate_magma_covariate_output,
)
from postgwas.modules.magmacovar.service import (
    preflight_magmacovar_pipeline,
    run_magma_covar_direct,
)


def _gene_results(path: Path, genes: int = 10) -> Path:
    rows = ["# VERSION = 110", "# COVAR = NSAMP MAC"]
    rows.extend(
        "%d 1 %d %d 20 5 1000 4.5 %.4f"
        % (gene, gene * 100, gene * 100 + 50, gene / 10)
        for gene in range(1, genes + 1)
    )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _covariates(path: Path, genes: int = 10, *, first_value: str = "0.1") -> Path:
    rows = ["GENE Property Average"]
    rows.extend(
        "%d %s %.3f"
        % (gene, first_value if gene == 1 else gene / 10, gene / 20)
        for gene in range(1, genes + 1)
    )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _fake_magma(path: Path, *, create_results: bool = True, exit_code: int = 0) -> Path:
    results_block = "" if not create_results else """
results.write_text(
    "TOTAL_GENES = 10\\n"
    "TEST_DIRECTION = two-sided (covar)\\n"
    "VARIABLE TYPE NGENES BETA BETA_STD SE P\\n"
    "Property COVAR 10 0.2 0.1 0.05 0.01\\n"
    "Average COVAR 10 -0.1 -0.05 0.04 0.02\\n",
    encoding="utf-8",
)
"""
    path.write_text(
        "#!%s\n" % sys.executable
        + "import pathlib\n"
        + "import sys\n"
        + "arguments = sys.argv[1:]\n"
        + "prefix = pathlib.Path(arguments[arguments.index('--out') + 1])\n"
        + "prefix.parent.mkdir(parents=True, exist_ok=True)\n"
        + "results = pathlib.Path(str(prefix) + '.gsa.out')\n"
        + "native_log = pathlib.Path(str(prefix) + '.log')\n"
        + "native_log.write_text(' '.join(arguments) + '\\n', encoding='utf-8')\n"
        + results_block
        + "raise SystemExit(%d)\n" % exit_code,
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _arguments(tmp_path: Path, magma: Path) -> argparse.Namespace:
    return argparse.Namespace(
        magma_gene_results_file=str(_gene_results(tmp_path / "study.genes.raw")),
        covariates=str(_covariates(tmp_path / "covariates.tsv")),
        magma=str(magma),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
        resume=True,
    )


def test_packaged_configuration_declares_the_scientific_defaults():
    module = load_configuration().modules.magmacovar

    assert module.model == []
    assert module.model_use_cases["marginal"].model == []
    assert module.model_use_cases["marginal"].direction == "two-sided"
    assert module.model_use_cases["tissue_specificity"].model == [
        "condition-hide=Average"
    ]
    assert module.model_use_cases["tissue_specificity"].direction == "greater"
    assert tuple(module.model_placeholders) == ("PROPERTY", "INTERNAL", "PATH")
    assert tuple(module.model_modifiers) == (
        "analyse",
        "correct",
        "condition",
        "condition-hide",
        "condition-residualize",
        "joint",
        "joint-pairs",
    )
    for specification in module.model_modifiers.values():
        assert specification.published_example
        assert specification.example.startswith("--covariate-model ")
    assert module.direction == "two-sided"
    assert module.minimum_genes == 10
    assert module.input.missing_values == "median"
    assert module.input.maximum_missing_fraction == 0.05
    assert module.input.missing_genes == "drop"


def test_command_uses_documented_gene_covar_and_direction_modifiers(tmp_path):
    command = build_magma_covariate_command(
        "magma",
        tmp_path / "study.genes.raw",
        tmp_path / "covariates.tsv",
        tmp_path / "study",
        model=["condition-hide=Average"],
        direction="greater",
        missing_values="median",
        maximum_missing_fraction=0.05,
        missing_genes="drop",
    )

    assert command == [
        "magma",
        "--gene-results", str(tmp_path / "study.genes.raw"),
        "--gene-covar", str(tmp_path / "covariates.tsv"),
        "missing-values=median", "max-miss=0.05",
        "--model", "condition-hide=Average", "direction-covar=greater",
        "--out", str(tmp_path / "study"),
    ]


def test_configuration_rejects_direction_word_that_magma_does_not_accept():
    with pytest.raises(ConfigurationError, match="modules.magmacovar.direction"):
        load_configuration(cli_overrides={"modules.magmacovar.direction": "less"})


@pytest.mark.parametrize(
    "model",
    (
        ["interaction-pairs"],
        ["direction-covar=greater"],
        ["condition-hide"],
        ["joint-pairs=all"],
    ),
)
def test_configuration_rejects_unavailable_or_malformed_model_modifiers(model):
    with pytest.raises(ConfigurationError, match="modules.magmacovar"):
        load_configuration(cli_overrides={"modules.magmacovar.model": model})


def test_input_validation_requires_gene_overlap_and_numeric_properties(tmp_path):
    raw = _gene_results(tmp_path / "study.genes.raw")
    covariates = _covariates(tmp_path / "covariates.tsv", genes=9)

    with pytest.raises(MagmaCovarError, match="Only 9 gene IDs overlap"):
        validate_magma_covariate_inputs(
            raw,
            covariates,
            minimum_genes=10,
            maximum_missing_fraction=0.05,
            missing_genes="drop",
        )

    invalid = _covariates(
        tmp_path / "invalid_covariates.tsv", first_value="not-a-number",
    )
    with pytest.raises(MagmaCovarError, match="neither numeric nor"):
        validate_magma_covariate_inputs(
            raw,
            invalid,
            minimum_genes=10,
            maximum_missing_fraction=0.05,
            missing_genes="drop",
        )


def test_input_validation_rejects_properties_magma_would_silently_discard(tmp_path):
    raw = _gene_results(tmp_path / "study.genes.raw")
    covariates = _covariates(
        tmp_path / "covariates.tsv", first_value="NA",
    )

    with pytest.raises(MagmaCovarError, match="above the configured MAGMA maximum"):
        validate_magma_covariate_inputs(
            raw,
            covariates,
            minimum_genes=10,
            maximum_missing_fraction=0.05,
            missing_genes="drop",
        )


def test_service_runs_fake_magma_and_publishes_only_after_validation(tmp_path):
    magma = _fake_magma(tmp_path / "magma")
    args = _arguments(tmp_path, magma)

    result = run_magma_covar_direct(args)

    output = tmp_path / "results"
    assert result == str(output / "study.gsa.out")
    assert (output / "study.gsa.out").is_file()
    assert (output / "study.log").is_file()
    assert (output / "logs" / "study_magmacovar.log").is_file()
    assert (output / "logs" / "study_magmacovar_resolved_config.yaml").is_file()
    assert (output / "logs" / "study_magmacovar_completion.yaml").is_file()
    assert not (output / ".magmacovar_staging" / "study").exists()

    native_log = (output / "study.log").read_text(encoding="utf-8")
    assert "--gene-covar" in native_log
    assert "missing-values=median" in native_log
    assert "max-miss=0.05" in native_log
    assert "--model direction-covar=two-sided" in native_log
    assert validate_magma_covariate_output(
        output / "study.gsa.out", minimum_genes=10,
    )["tested_properties"] == 2


def test_success_exit_without_required_output_is_a_failure_and_stays_isolated(tmp_path):
    magma = _fake_magma(tmp_path / "magma", create_results=False)
    args = _arguments(tmp_path, magma)

    with pytest.raises(MagmaCovarError, match="expected output is missing or empty"):
        run_magma_covar_direct(args)

    output = tmp_path / "results"
    assert not (output / "study.gsa.out").exists()
    assert (output / ".magmacovar_staging" / "study").exists()
    canonical_log = output / "logs" / "study_magmacovar.log"
    assert "FAILED" in canonical_log.read_text(encoding="utf-8")


def test_resume_revalidates_completion_without_rerunning_magma(tmp_path):
    magma = _fake_magma(tmp_path / "magma")
    args = _arguments(tmp_path, magma)
    first = run_magma_covar_direct(args)
    native_log = Path(first).parent / "study.log"
    before = native_log.read_text(encoding="utf-8")

    second = run_magma_covar_direct(args)

    assert second == first
    assert native_log.read_text(encoding="utf-8") == before


def test_input_sharing_the_output_prefix_is_never_treated_as_owned_output(tmp_path):
    magma = _fake_magma(tmp_path / "magma")
    gene_results = _gene_results(tmp_path / "study.genes.raw")
    args = argparse.Namespace(
        magma_gene_results_file=str(gene_results),
        covariates=str(_covariates(tmp_path / "covariates.tsv")),
        magma=str(magma),
        dataset_id="study",
        output_directory=str(tmp_path),
        resume=True,
    )

    run_magma_covar_direct(args)

    assert gene_results.is_file()
    assert (tmp_path / "study.gsa.out").is_file()


def test_overwrite_preserves_unrecorded_files_with_the_same_prefix(tmp_path):
    magma = _fake_magma(tmp_path / "magma")
    args = _arguments(tmp_path, magma)
    run_magma_covar_direct(args)
    notes = tmp_path / "results" / "study.notes"
    notes.write_text("user-owned\n", encoding="utf-8")
    args.overwrite = True

    run_magma_covar_direct(args)

    assert notes.read_text(encoding="utf-8") == "user-owned\n"


def test_overwrite_refuses_to_remove_unowned_staging_files(tmp_path):
    magma = _fake_magma(tmp_path / "magma")
    args = _arguments(tmp_path, magma)
    unowned = (
        tmp_path / "results" / ".magmacovar_staging" / "study" / "study.notes"
    )
    unowned.parent.mkdir(parents=True)
    unowned.write_text("user-owned\n", encoding="utf-8")
    args.overwrite = True

    with pytest.raises(MagmaCovarError, match="not owned by this output prefix"):
        run_magma_covar_direct(args)

    assert unowned.read_text(encoding="utf-8") == "user-owned\n"


def test_pipeline_preflight_rejects_bad_covariates_before_upstream_work(tmp_path):
    magma = _fake_magma(tmp_path / "magma")
    invalid = _covariates(
        tmp_path / "invalid_covariates.tsv", first_value="not-a-number",
    )
    args = argparse.Namespace(
        covariates=str(invalid),
        magma=str(magma),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    with pytest.raises(MagmaCovarError, match="neither numeric nor"):
        preflight_magmacovar_pipeline(args)

    preflight_log = tmp_path / "results" / "logs" / "study_magmacovar.log"
    assert "configuration failed" in preflight_log.read_text(encoding="utf-8")
