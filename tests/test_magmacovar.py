"""Focused scientific and execution regressions for MAGMAcovar."""

from __future__ import annotations

import argparse
import csv
from io import StringIO
from pathlib import Path
import sys
from unittest.mock import Mock

import pytest
from rich.console import Console
import yaml

from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import PipelineStageController
from postgwas.modules.magma.reporting import (
    MAGMA_GENE_ONLY_STAGE_KEYS,
    MAGMA_STAGE_TITLES,
)
from postgwas.modules.magmacovar.errors import MagmaCovarError
from postgwas.modules.magmacovar.main import (
    build_magma_covariate_command,
    validate_corrected_magma_covariate_output,
    validate_magma_covariate_inputs,
    validate_magma_covariate_output,
    write_corrected_magma_covariate_results,
)
from postgwas.modules.magmacovar.service import (
    preflight_magmacovar_pipeline,
    resolve_magmacovar_configuration,
    run_magma_covar_direct,
)
from postgwas.modules.magmacovar.stages import (
    MAGMACOVAR_STAGES,
    magmacovar_pipeline_progress_plan,
)
from postgwas.pipeline.planner import build_pipeline_plan
from postgwas.pipeline.registry import REGISTRY
from preflight_support import pipeline_input_vcf_evidence


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


def _fake_magma(
    path: Path,
    *,
    create_results: bool = True,
    create_unexpected_output: bool = False,
    exit_code: int = 0,
) -> Path:
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
    unexpected_block = "" if not create_unexpected_output else (
        "pathlib.Path(str(prefix) + '.unexpected').write_text("
        "'unexpected\\n', encoding='utf-8')\n"
    )
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
        + unexpected_block
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
    assert module.multiple_testing.methods == [
        "bonferroni", "sidak", "holm", "fdr_bh",
    ]
    assert module.multiple_testing.primary_method == "bonferroni"
    assert module.multiple_testing.significance_threshold == 0.05
    assert module.result_schema.adjusted_p_value_column_pattern == "P_{method}_corr"
    assert module.reporting.highlight_method == "fdr_bh"
    assert module.reporting.top_property_count == 5
    assert module.reporting.p_value_significant_digits == 4
    assert module.reporting.effect_significant_digits == 4
    assert module.output_layout.corrected_results_file == (
        "{dataset_id}_magmacovar_corrected.tsv"
    )


def test_pipeline_plan_preserves_magma_validation_before_magmacovar_stages():
    plan = magmacovar_pipeline_progress_plan(argparse.Namespace())
    execution = build_pipeline_plan(["magmacovar"])

    assert plan is not None
    assert plan["kind"] == "magmacovar"
    assert plan["label"] == "MAGMAcovar pipeline execution progress"
    upstream_stages = tuple(
        MAGMA_STAGE_TITLES[key] for key in MAGMA_GENE_ONLY_STAGE_KEYS
    )
    assert plan["stages"] == (*upstream_stages, *MAGMACOVAR_STAGES)
    assert len(plan["stages"]) == 12
    assert all("pathway" not in stage.lower() for stage in plan["stages"])
    assert plan["modules"] == {
        "formatter": (1, 4),
        "magma": (4, 8),
        "magmacovar": (9, 12),
    }
    assert plan["magma_stage_numbers"] == {
        key: number for number, key in enumerate(MAGMA_GENE_ONLY_STAGE_KEYS, 1)
    }
    assert plan["stage_numbers"] == {
        "gene_results": 9,
        "covariates": 10,
        "analysis": 11,
        "results": 12,
    }
    assert REGISTRY.get("magmacovar").pipeline_progress_factory == (
        "postgwas.modules.magmacovar.stages:"
        "magmacovar_pipeline_progress_plan"
    )
    assert execution.steps == ("formatter", "magma", "magmacovar")


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


def test_configuration_rejects_invalid_multiple_testing_contracts():
    with pytest.raises(ConfigurationError, match="must occur in methods"):
        load_configuration(cli_overrides={
            "modules.magmacovar.multiple_testing.methods": ["holm"],
            "modules.magmacovar.multiple_testing.primary_method": "bonferroni",
        })
    with pytest.raises(ConfigurationError, match="must contain unique"):
        load_configuration(cli_overrides={
            "modules.magmacovar.multiple_testing.methods": ["holm", "holm"],
        })
    with pytest.raises(ConfigurationError, match="columns must be unique"):
        load_configuration(cli_overrides={
            "modules.magmacovar.result_schema.primary_adjusted_p_value_column": (
                "P_bonferroni_corr"
            ),
        })


@pytest.mark.parametrize(
    ("key", "value"),
    (
        ("top_property_count", 0),
        ("p_value_significant_digits", 11),
        ("effect_significant_digits", 0),
    ),
)
def test_configuration_rejects_invalid_main_findings_settings(key, value):
    with pytest.raises(ConfigurationError, match=rf"reporting\.{key}"):
        load_configuration(cli_overrides={
            f"modules.magmacovar.reporting.{key}": value,
        })


def test_findings_highlight_method_must_be_a_selected_correction():
    with pytest.raises(
        ConfigurationError,
        match="reporting.highlight_method must occur",
    ):
        load_configuration(cli_overrides={
            "modules.magmacovar.multiple_testing.methods": ["bonferroni"],
        })


def test_cli_missing_data_policies_override_the_same_yaml_fields():
    configuration = resolve_magmacovar_configuration(argparse.Namespace(
        covariate_missing_values="mean",
        covariate_max_miss=0.12,
        covariate_missing_genes="fill",
    ))

    module = configuration.modules.magmacovar
    assert module.input.missing_values == "mean"
    assert module.input.maximum_missing_fraction == pytest.approx(0.12)
    assert module.input.missing_genes == "fill"


def test_cli_max_miss_remains_bounded_by_the_canonical_schema():
    with pytest.raises(
        ConfigurationError,
        match="modules.magmacovar.input.maximum_missing_fraction",
    ):
        resolve_magmacovar_configuration(argparse.Namespace(
            covariate_max_miss=0.21,
        ))


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

    with pytest.raises(
        MagmaCovarError,
        match=r"Gene property 'Property'.*above the configured MAGMA maximum",
    ) as captured:
        validate_magma_covariate_inputs(
            raw,
            covariates,
            minimum_genes=10,
            maximum_missing_fraction=0.05,
            missing_genes="drop",
        )
    message = str(captured.value)
    assert "1/10 missing values (0.100)" in message
    assert "missing-genes=drop" in message
    assert "gene IDs overlapping both input files as the denominator" in message
    assert "--covariate-max-miss FRACTION" in message
    assert "modules.magmacovar.input.maximum_missing_fraction" in message
    assert "Changing the missing-values policy does not bypass max-miss" in message


def test_missingness_denominator_follows_missing_genes_policy(tmp_path):
    raw = _gene_results(tmp_path / "study.genes.raw")
    covariates = _covariates(
        tmp_path / "covariates.tsv", genes=9, first_value="NA",
    )

    dropped = validate_magma_covariate_inputs(
        raw,
        covariates,
        minimum_genes=9,
        maximum_missing_fraction=0.15,
        missing_genes="drop",
    )
    property_missingness = dropped["property_missingness"][0]
    assert property_missingness["missing_genes"] == 1
    assert property_missingness["missing_fraction"] == pytest.approx(1 / 9)

    with pytest.raises(
        MagmaCovarError,
        match=r"Gene property 'Property' has 2/10 missing values \(0\.200\)",
    ) as captured:
        validate_magma_covariate_inputs(
            raw,
            covariates,
            minimum_genes=9,
            maximum_missing_fraction=0.15,
            missing_genes="fill",
        )
    message = str(captured.value)
    assert "2/10 missing values (0.200)" in message
    assert "all eligible .genes.raw genes as the denominator" in message
    assert "modules.magmacovar.input.missing_genes to drop" in message
    assert "--covariate-missing-genes drop" in message


def test_missingness_equal_to_maximum_is_accepted_with_fill_policy(tmp_path):
    summary = validate_magma_covariate_inputs(
        _gene_results(tmp_path / "study.genes.raw"),
        _covariates(tmp_path / "covariates.tsv", genes=9),
        minimum_genes=9,
        maximum_missing_fraction=0.1,
        missing_genes="fill",
    )

    property_missingness = summary["property_missingness"][0]
    assert property_missingness["missing_genes"] == 1
    assert property_missingness["missing_fraction"] == pytest.approx(0.1)


def test_corrections_cover_one_complete_gene_property_family(tmp_path):
    configuration = load_configuration()
    module = configuration.modules.magmacovar
    raw = tmp_path / "study.gsa.out"
    raw.write_text(
        "VARIABLE TYPE NGENES BETA BETA_STD SE P\n"
        "Property COVAR 10 0.2 0.1 0.05 0.01\n"
        "Average COVAR 10 -0.1 -0.05 0.04 0.02\n"
        "ignored SET 10 0.5 0.4 0.03 0.001\n",
        encoding="utf-8",
    )
    corrected = tmp_path / "study_magmacovar_corrected.tsv"

    summary = write_corrected_magma_covariate_results(
        raw, corrected, module=module, logger=Mock(),
    )

    with corrected.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["VARIABLE"] for row in rows] == ["Property", "Average"]
    assert [float(row["P_bonferroni_corr"]) for row in rows] == [0.02, 0.04]
    assert [float(row["P_sidak_corr"]) for row in rows] == pytest.approx(
        [0.0199, 0.0396]
    )
    assert [float(row["P_holm_corr"]) for row in rows] == [0.02, 0.02]
    assert [float(row["P_fdr_bh_corr"]) for row in rows] == [0.02, 0.02]
    assert {row["primary_correction_method"] for row in rows} == {"bonferroni"}
    assert {row["primary_significant"] for row in rows} == {"True"}
    assert summary["tested_properties"] == 2
    assert summary["primary_significant_properties"] == 2
    assert summary["significant_properties_by_method"] == {
        "bonferroni": 2,
        "sidak": 2,
        "holm": 2,
        "fdr_bh": 2,
    }
    assert [finding["property"] for finding in summary["top_properties"]] == [
        "Property",
        "Average",
    ]
    assert summary["top_properties"][0]["standardized_beta"] == pytest.approx(
        0.1
    )
    assert summary["top_properties"][0]["adjusted_p_values"][
        "fdr_bh"
    ] == pytest.approx(0.02)
    assert validate_corrected_magma_covariate_output(
        raw, corrected, module=module,
    ) == summary

    text = corrected.read_text(encoding="utf-8")
    corrected.write_text(
        text.replace("\t0.02\t", "\t0.03\t", 1), encoding="utf-8",
    )
    with pytest.raises(MagmaCovarError, match="does not match the native results"):
        validate_corrected_magma_covariate_output(
            raw, corrected, module=module,
        )


def test_full_name_is_used_for_truncated_native_magma_property_names(tmp_path):
    module = load_configuration().modules.magmacovar
    raw = tmp_path / "study.gsa.out"
    raw.write_text(
        "VARIABLE TYPE NGENES BETA BETA_STD SE P FULL_NAME\n"
        "Brain_Anterior_cingulate_cor... COVAR 10 0.2 0.1 0.05 0.01 "
        "Brain_Anterior_cingulate_cortex_BA24\n",
        encoding="utf-8",
    )
    corrected = tmp_path / "corrected.tsv"

    summary = write_corrected_magma_covariate_results(
        raw, corrected, module=module, logger=Mock(),
    )

    with corrected.open("r", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["VARIABLE"] == "Brain_Anterior_cingulate_cortex_BA24"
    assert summary["top_properties"][0]["property"] == (
        "Brain_Anterior_cingulate_cortex_BA24"
    )


def test_top_property_count_is_controlled_by_canonical_configuration(tmp_path):
    module = load_configuration(cli_overrides={
        "modules.magmacovar.reporting.top_property_count": 1,
    }).modules.magmacovar
    raw = tmp_path / "study.gsa.out"
    raw.write_text(
        "VARIABLE TYPE NGENES BETA BETA_STD SE P\n"
        "Second COVAR 10 0.2 0.1 0.05 0.02\n"
        "First COVAR 10 0.3 0.2 0.04 0.01\n",
        encoding="utf-8",
    )

    summary = write_corrected_magma_covariate_results(
        raw, tmp_path / "corrected.tsv", module=module, logger=Mock(),
    )

    assert [finding["property"] for finding in summary["top_properties"]] == [
        "First"
    ]


def test_primary_bonferroni_boundary_is_inclusive_for_461_properties(tmp_path):
    module = load_configuration().modules.magmacovar
    raw = tmp_path / "study.gsa.out"
    rows = ["VARIABLE TYPE NGENES BETA BETA_STD SE P"]
    rows.extend(
        "Cell%d COVAR 10 0.2 0.1 0.05 %.17g"
        % (index, 0.05 / 461 if index == 1 else 1.0)
        for index in range(1, 462)
    )
    raw.write_text("\n".join(rows) + "\n", encoding="utf-8")

    summary = write_corrected_magma_covariate_results(
        raw,
        tmp_path / "corrected.tsv",
        module=module,
        logger=Mock(),
    )

    assert summary["tested_properties"] == 461
    assert summary["primary_significant_properties"] == 1


def test_service_runs_fake_magma_and_publishes_only_after_validation(
    tmp_path, capsys,
):
    magma = _fake_magma(tmp_path / "magma")
    args = _arguments(tmp_path, magma)

    result = run_magma_covar_direct(args)

    output = tmp_path / "results"
    assert result["raw_results"] == str(output / "study.gsa.out")
    assert result["corrected_results"] == str(
        output / "study_magmacovar_corrected.tsv"
    )
    assert result["primary_correction_method"] == "bonferroni"
    assert result["primary_significant_properties"] == 2
    assert result["significant_properties_by_method"]["fdr_bh"] == 2
    assert [finding["property"] for finding in result["top_properties"]] == [
        "Property",
        "Average",
    ]
    assert (output / "study.gsa.out").is_file()
    assert (output / "study_magmacovar_corrected.tsv").is_file()
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
    completion = yaml.safe_load(
        (output / "logs" / "study_magmacovar_completion.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert set(completion["outputs"]) == {
        "raw_results", "corrected_results", "native_log",
    }
    assert completion["metrics"]["gene_results_file"].endswith(
        "study.genes.raw"
    )
    assert completion["metrics"]["covariates_file"].endswith(
        "covariates.tsv"
    )
    assert completion["metrics"]["gene_results_genes"] == 10
    assert completion["metrics"]["covariate_genes"] == 10
    assert completion["metrics"]["overlapping_genes"] == 10
    assert completion["metrics"]["properties"] == 2
    assert completion["metrics"]["significant_properties_by_method"][
        "fdr_bh"
    ] == 2
    assert completion["metrics"]["top_properties"][0]["property"] == (
        "Property"
    )
    screen = capsys.readouterr().out
    compact_screen = " ".join(screen.split())
    assert "MAGMAcovar analysis progress" in screen
    assert "All 4 stages completed" in screen
    for number, stage in enumerate(MAGMACOVAR_STAGES, 1):
        assert "Completed %d/4 · %s" % (number, stage) in compact_screen
    assert "Eligible gene identifiers : 10" in compact_screen
    assert "Gene IDs shared with MAGMA results : 10/10 (100.00%)" in (
        compact_screen
    )
    assert "File-structure validation : passed" in compact_screen
    assert "MAGMAcovar main findings" in screen
    assert "Genes per tested property : 10–10" in compact_screen
    assert "Gene-property model : marginal" in compact_screen
    assert "Property-test direction : two-sided" in compact_screen
    assert "Bonferroni-significant properties : 2 / 2" in compact_screen
    assert "BH-FDR-significant properties : 2 / 2" in compact_screen
    assert "Top-property ranking : lowest raw MAGMA P; showing 2" in (
        compact_screen
    )
    assert (
        "Top property 1 : Property (BETA_STD=0.1; P=0.01; BH-FDR=0.02)"
        in compact_screen
    )
    assert "Native MAGMA result" in screen
    assert "FLAMES handoff" not in screen
    assert "not the corrected TSV" not in screen
    assert "Corrected results for interpretation" in screen
    assert "Validate and publish all MAGMAcovar outputs" not in screen
    assert "Published MAGMAcovar outputs" not in screen
    assert "Completion checkpoint" not in screen
    canonical_log = (
        output / "logs" / "study_magmacovar.log"
    ).read_text(encoding="utf-8")
    assert canonical_log.count("stage_outcome") == 4
    assert "publish_magmacovar_outputs" not in canonical_log
    assert "Published MAGMAcovar outputs" not in canonical_log
    for stage in MAGMACOVAR_STAGES:
        assert stage in canonical_log


def test_main_findings_include_flames_handoff_only_when_flames_is_requested(
    tmp_path, capsys,
):
    args = _arguments(tmp_path, _fake_magma(tmp_path / "magma"))
    args._pipeline_requested_modules = ("magmacovar", "flames")

    run_magma_covar_direct(args)

    screen = capsys.readouterr().out
    assert "FLAMES handoff" in screen
    assert "native .gsa.out with raw P values" in screen


def test_detailed_pipeline_progress_continues_from_magma_through_publication(
    tmp_path,
):
    magma = _fake_magma(tmp_path / "magma")
    args = _arguments(tmp_path, magma)
    plan = magmacovar_pipeline_progress_plan(args)
    assert plan is not None
    stream = StringIO()
    controller = PipelineStageController(
        plan["label"],
        plan["stages"],
        console=Console(
            file=stream,
            force_terminal=True,
            color_system=None,
            width=120,
        ),
        outcome_label_width=38,
    )
    for number in range(1, len(MAGMA_GENE_ONLY_STAGE_KEYS) + 1):
        controller.start(number)
        controller.complete(number, outcome="validated upstream MAGMA stage")
    args._pipeline_progress_plan = plan
    args._pipeline_stage_progress = controller

    result = run_magma_covar_direct(args)

    assert result["tested_properties"] == 2
    assert controller.completed == 11
    assert controller.current == 12
    completion = args._pipeline_stage_completion
    controller.complete(
        12,
        outcome=completion.get("outcome"),
        outcome_fields=completion.get("outcome_fields"),
    )
    controller.close()
    compact_screen = " ".join(stream.getvalue().split())
    assert "MAGMAcovar pipeline execution progress" in compact_screen
    assert "Completed 12/12 · Validate and adjust the MAGMA gene-property results" in (
        compact_screen
    )
    assert "All 12 stages completed" in compact_screen
    assert "pathway" not in compact_screen.lower()
    assert "Eligible gene identifiers : 10" in compact_screen
    assert "Covariate-table genes : 10" in compact_screen
    assert "Gene properties : 2" in compact_screen
    assert "Gene IDs shared with MAGMA results : 10/10 (100.00%)" in (
        compact_screen
    )
    assert "Native COVAR-result validation : passed" in compact_screen
    assert "Published MAGMAcovar outputs" not in compact_screen
    assert "Completion checkpoint" not in compact_screen

    resumed_stream = StringIO()
    resumed_controller = PipelineStageController(
        plan["label"],
        plan["stages"],
        console=Console(
            file=resumed_stream,
            force_terminal=True,
            color_system=None,
            width=120,
        ),
        outcome_label_width=38,
    )
    for number in range(1, len(MAGMA_GENE_ONLY_STAGE_KEYS) + 1):
        resumed_controller.start(number)
        resumed_controller.complete(
            number, outcome="resumed validated upstream MAGMA stage",
        )
    args._pipeline_stage_progress = resumed_controller

    resumed = run_magma_covar_direct(args)

    assert resumed == result
    assert resumed_controller.completed == 11
    assert resumed_controller.current == 12
    resumed_completion = args._pipeline_stage_completion
    resumed_controller.complete(
        12,
        outcome=resumed_completion.get("outcome"),
        outcome_fields=resumed_completion.get("outcome_fields"),
    )
    resumed_controller.close()
    compact_resumed_screen = " ".join(resumed_stream.getvalue().split())
    assert "Input fingerprint : validated and unchanged" in compact_resumed_screen
    assert "Native execution : reused from validated checkpoint" in (
        compact_resumed_screen
    )
    assert "Published outputs" not in compact_resumed_screen
    assert "Completion checkpoint" not in compact_resumed_screen
    assert "All 12 stages completed" in compact_resumed_screen
    assert "pathway" not in compact_resumed_screen.lower()


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


def test_internal_publication_failure_keeps_final_visible_stage_incomplete(
    tmp_path, capsys,
):
    magma = _fake_magma(
        tmp_path / "magma", create_unexpected_output=True,
    )
    args = _arguments(tmp_path, magma)

    with pytest.raises(MagmaCovarError, match="unconfigured staged artifacts"):
        run_magma_covar_direct(args)

    compact_screen = " ".join(capsys.readouterr().out.split())
    assert "Failed 4/4 · Validate and adjust the MAGMA gene-property results" in (
        compact_screen
    )
    assert "All 4 stages completed" not in compact_screen
    assert not (tmp_path / "results" / "study.gsa.out").exists()
    assert not (
        tmp_path / "results" / "logs" / "study_magmacovar_completion.yaml"
    ).exists()


def test_detailed_pipeline_failure_names_the_active_gene_property_stage(tmp_path):
    magma = _fake_magma(tmp_path / "magma", exit_code=3)
    args = _arguments(tmp_path, magma)
    plan = magmacovar_pipeline_progress_plan(args)
    assert plan is not None
    stream = StringIO()
    controller = PipelineStageController(
        plan["label"],
        plan["stages"],
        console=Console(
            file=stream,
            force_terminal=True,
            color_system=None,
            width=120,
        ),
    )
    for number in range(1, len(MAGMA_GENE_ONLY_STAGE_KEYS) + 1):
        controller.start(number)
        controller.complete(number)
    args._pipeline_progress_plan = plan
    args._pipeline_stage_progress = controller

    with pytest.raises(MagmaCovarError):
        run_magma_covar_direct(args)
    controller.fail_active()
    controller.close()

    compact_screen = " ".join(stream.getvalue().split())
    assert controller.completed == 10
    assert "Failed 11/12 · Run MAGMA gene-property analysis" in compact_screen
    assert "All 12 stages completed" not in compact_screen
    assert "100%" not in compact_screen


def test_resume_revalidates_completion_without_rerunning_magma(tmp_path, capsys):
    magma = _fake_magma(tmp_path / "magma")
    args = _arguments(tmp_path, magma)
    first = run_magma_covar_direct(args)
    native_log = Path(first["native_log"])
    before = native_log.read_text(encoding="utf-8")

    second = run_magma_covar_direct(args)

    assert second == first
    assert native_log.read_text(encoding="utf-8") == before
    compact_screen = " ".join(capsys.readouterr().out.split())
    assert "Input fingerprint : validated and unchanged" in compact_screen
    assert "Native execution : reused from validated checkpoint" in compact_screen
    assert "All 4 stages completed" in compact_screen
    assert "Completion checkpoint" not in compact_screen


def test_resume_restarts_an_older_two_artifact_output_contract(tmp_path):
    magma = _fake_magma(tmp_path / "magma")
    args = _arguments(tmp_path, magma)
    first = run_magma_covar_direct(args)
    completion_path = (
        tmp_path / "results" / "logs" / "study_magmacovar_completion.yaml"
    )
    completion = yaml.safe_load(completion_path.read_text(encoding="utf-8"))
    completion["outputs"].pop("corrected_results")
    completion_path.write_text(
        yaml.safe_dump(completion, sort_keys=False), encoding="utf-8",
    )
    Path(first["corrected_results"]).unlink()

    restarted = run_magma_covar_direct(args)

    assert Path(restarted["raw_results"]).is_file()
    assert Path(restarted["corrected_results"]).is_file()
    canonical_log = (
        tmp_path / "results" / "logs" / "study_magmacovar.log"
    ).read_text(encoding="utf-8")
    assert "reason=changed_output_contract" in canonical_log


def test_resume_restarts_when_every_declared_output_is_missing(tmp_path):
    magma = _fake_magma(tmp_path / "magma")
    args = _arguments(tmp_path, magma)
    first_result = run_magma_covar_direct(args)
    first = Path(first_result["raw_results"])
    native_log = first.parent / "study.log"
    corrected = Path(first_result["corrected_results"])
    completion = first.parent / "logs" / "study_magmacovar_completion.yaml"
    first.unlink()
    native_log.unlink()

    restarted = Path(run_magma_covar_direct(args)["raw_results"])

    assert restarted.is_file()
    assert native_log.is_file()
    assert corrected.is_file()
    assert completion.is_file()
    canonical_log = first.parent / "logs" / "study_magmacovar.log"
    assert "reason=incomplete_outputs" in canonical_log.read_text(
        encoding="utf-8"
    )


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
        preflight_magmacovar_pipeline(
            args, preflight_evidence=pipeline_input_vcf_evidence(),
        )

    preflight_log = tmp_path / "results" / "logs" / "study_magmacovar.log"
    assert "configuration failed" in preflight_log.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    ("destination", "flag", "value"),
    (
        ("gene_set_file", "--gene-set-file", "pathways.gmt"),
        ("minimum_gene_id_overlap", "--minimum-gene-id-overlap", 0.5),
        (
            "gene_set_identifier_mismatch",
            "--gene-set-identifier-mismatch",
            "error",
        ),
        (
            "alternate_gene_id_duplicate_policy",
            "--alternate-gene-id-duplicate-policy",
            "error",
        ),
        (
            "gene_location_alternate_id_type",
            "--gene-location-alternate-id-type",
            "symbol",
        ),
    ),
)
def test_pipeline_preflight_rejects_explicit_pathway_arguments(
    tmp_path, destination, flag, value,
):
    args = argparse.Namespace(
        modules=["magmacovar"],
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )
    setattr(args, destination, value)

    with pytest.raises(
        MagmaCovarError,
        match=rf"{flag}.*does not run pathway analysis",
    ):
        preflight_magmacovar_pipeline(
            args, preflight_evidence=pipeline_input_vcf_evidence(),
        )

    preflight_log = tmp_path / "results" / "logs" / "study_magmacovar.log"
    assert flag in preflight_log.read_text(encoding="utf-8")
