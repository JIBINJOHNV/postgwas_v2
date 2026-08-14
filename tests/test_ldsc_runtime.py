"""Regression tests for canonical CBIIT LDSC integration."""

from argparse import Namespace
import argparse
from pathlib import Path
import sys

import pytest
from pydantic import ValidationError
import yaml

from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.modules.ldsc.cli import build_parser, main as ldsc_main
from postgwas.modules.ldsc.ldsc_runner import (
    LDSCError,
    build_h2_command,
    build_munge_sumstats_command,
    extract_ldsc_metrics,
    run_ldsc,
)
from postgwas.modules.ldsc.service import (
    resolve_ldsc_configuration,
    run_ldsc_direct,
)


def _write_inputs(root: Path) -> tuple[Path, Path, Path, Path]:
    sumstats = root / "input.tsv"
    merge = root / "w_hm3.snplist"
    reference = root / "reference"
    weights = root / "weights"
    sumstats.write_text(
        "SNP\tA1\tA2\tZ\tP\tN\tFRQ\tINFO\n"
        "rs1\tA\tG\t1\t0.3\t1000\t0.2\t0.99\n",
        encoding="utf-8",
    )
    merge.write_text("SNP\tA1\tA2\nrs1\tA\tG\n", encoding="utf-8")
    reference.mkdir()
    weights.mkdir()
    module = load_configuration().modules.ldsc
    for chromosome in range(1, module.chromosomes + 1):
        for directory, suffixes in (
            (
                reference,
                (
                    module.reference_layout.ld_score_suffix,
                    module.reference_layout.m_suffix,
                    module.reference_layout.m_5_50_suffix,
                ),
            ),
            (weights, (module.reference_layout.ld_score_suffix,)),
        ):
            for suffix in suffixes:
                (directory / (str(chromosome) + suffix)).write_text(
                    "reference\n", encoding="utf-8",
                )
    return sumstats, merge, reference, weights


def _service_args(root: Path, **values) -> Namespace:
    sumstats, merge, reference, weights = _write_inputs(root)
    supplied = {
        "ldsc_input": str(sumstats),
        "merge_alleles": str(merge),
        "ref_ld_chr": str(reference),
        "w_ld_chr": str(weights),
        "output_directory": str(root / "output"),
        "dataset_id": "STUDY",
        "ldsc_executable": sys.executable,
        "munge_sumstats_executable": sys.executable,
    }
    supplied.update(values)
    return Namespace(**supplied)


def _fake_ldsc_result(**kwargs):
    prefix = Path(kwargs["output_prefix"])
    module = kwargs["configuration"]
    layout = module.output_layout
    munged = Path(str(prefix) + layout.munged_sumstats_suffix)
    munge_log = Path(str(prefix) + layout.upstream_log_suffix)
    observed = Path(
        str(prefix) + layout.observed_prefix_suffix + layout.upstream_log_suffix
    )
    paths = [munged, munge_log, observed]
    liability = None
    liability_outputs = []
    if module.population_prevalence is not None:
        liability = Path(
            str(prefix)
            + layout.liability_prefix_suffix
            + layout.upstream_log_suffix
        )
        paths.append(liability)
        liability_outputs = [str(liability)]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("validated\n", encoding="utf-8")
    return {
        "munged_sumstats": str(munged),
        "munge_log": str(munge_log),
        "h2_observed": str(observed),
        "observed_outputs": [str(observed)],
        "observed_metrics": {
            "h2": "0.20 (0.03)", "intercept": "1.01", "ratio": "0.1",
        },
        "h2_liability": None if liability is None else str(liability),
        "liability_outputs": liability_outputs,
        "liability_metrics": (
            None
            if liability is None
            else {
                "h2": "0.35 (0.05)",
                "intercept": "1.02 (0.01)",
                "ratio": "0.2 (0.1)",
            }
        ),
    }


def test_packaged_defaults_match_pinned_upstream_ldsc():
    module = load_configuration().modules.ldsc
    assert module.minimum_info == 0.9
    assert module.minimum_maf == 0.01
    assert module.minimum_n is None
    assert module.chunksize == 5_000_000
    assert module.keep_maf is False
    assert module.intercept is None
    assert module.two_step is None
    assert module.chisq_max is None
    assert module.n_blocks == 200
    assert module.use_m_5_50 is True
    assert module.print_covariance is False
    assert module.print_delete_values is False
    assert module.sample_prevalence is None
    assert module.population_prevalence is None
    assert module.population is None
    assert module.genome_build is None
    assert module.chromosomes == 22


def test_ldsc_configurable_cli_actions_do_not_own_defaults():
    parser = build_parser()
    destinations = (
        "ldsc_population", "ldsc_genome_build", "ldsc_minimum_info",
        "ldsc_minimum_maf", "ldsc_minimum_n",
        "ldsc_chunksize", "ldsc_keep_maf", "ldsc_intercept",
        "ldsc_two_step", "ldsc_chisq_max", "ldsc_n_blocks",
        "ldsc_use_m_5_50", "ldsc_print_covariance",
        "ldsc_print_delete_values", "samp_prev", "pop_prev",
        "ldsc_executable", "munge_sumstats_executable", "run_config",
        "overwrite",
    )
    for destination in destinations:
        action = next(item for item in parser._actions if item.dest == destination)
        assert action.default == argparse.SUPPRESS

    for destination in ("ldsc_population", "ldsc_genome_build"):
        action = next(item for item in parser._actions if item.dest == destination)
        assert "unset" in action.help


def test_use_m_5_50_help_displays_canonical_default_in_green():
    parser = build_parser()
    action = next(
        item for item in parser._actions
        if item.dest == "ldsc_use_m_5_50"
    )

    assert "[bold green]Default:[/bold green]" in action.help
    assert "[green]true[/green]" in action.help
    assert "Default: true" in parser.format_help()


def test_direct_cli_requires_every_runtime_input_and_complete_examples(tmp_path):
    sumstats, merge, reference, weights = _write_inputs(tmp_path)
    parser = build_parser()
    for destination in (
        "ldsc_input", "merge_alleles", "ref_ld_chr", "w_ld_chr",
        "dataset_id", "output_directory",
    ):
        action = next(item for item in parser._actions if item.dest == destination)
        assert action.required is True

    with pytest.raises(SystemExit):
        parser.parse_args([
            "--ldsc-input", str(sumstats),
            "--merge-alleles", str(merge),
            "--ref-ld-chr", str(reference),
            "--w-ld-chr", str(weights),
        ])

    assert parser.format_help().count(
        "--merge-alleles reference/w_hm3.snplist"
    ) == 2


def test_direct_cli_renders_existing_outputs_as_an_actionable_error(
    tmp_path, capsys,
):
    sumstats, merge, reference, weights = _write_inputs(tmp_path)
    output = tmp_path / "output"
    output.mkdir()
    existing = output / "STUDY.sumstats.gz"
    existing.write_text("previous validated result\n", encoding="utf-8")

    exit_code = ldsc_main([
        "--ldsc-input", str(sumstats),
        "--merge-alleles", str(merge),
        "--ref-ld-chr", str(reference),
        "--w-ld-chr", str(weights),
        "--dataset-id", "STUDY",
        "--output-directory", str(output),
        "--ldsc-executable", sys.executable,
        "--munge-sumstats-executable", sys.executable,
    ])

    captured = capsys.readouterr()
    terminal = captured.out + captured.err
    normalized_terminal = " ".join(terminal.split())
    compact_terminal = "".join(terminal.split())
    assert exit_code == 1
    assert "Traceback" not in terminal
    assert "LDSC heritability stopped" in normalized_terminal
    assert "Existing LDSC outputs were found" in normalized_terminal
    assert "stopped before changing any files" in normalized_terminal
    assert "--overwrite" in normalized_terminal
    assert "different --output-directory" in normalized_terminal
    assert str(existing) in compact_terminal
    assert "No new LDSC results were published" in normalized_terminal
    assert existing.read_text(encoding="utf-8") == "previous validated result\n"

    service_log = output / "logs" / "STUDY_ldsc_service.log"
    assert "FAILED" in service_log.read_text(encoding="utf-8")


def test_explicit_cli_overrides_module_yaml(tmp_path):
    run_config = tmp_path / "ldsc.yaml"
    run_config.write_text(
        "minimum_info: 0.8\nn_blocks: 80\nuse_m_5_50: false\n",
        encoding="utf-8",
    )
    configuration = resolve_ldsc_configuration(Namespace(
        run_config=str(run_config),
        ldsc_minimum_info=0.95,
        ldsc_n_blocks=120,
        ldsc_use_m_5_50=True,
    ))
    module = configuration.modules.ldsc
    assert module.minimum_info == 0.95
    assert module.n_blocks == 120
    assert module.use_m_5_50 is True
    assert module.minimum_maf == 0.01


def test_direct_cli_values_override_upstream_defaults(tmp_path):
    sumstats, merge, reference, weights = _write_inputs(tmp_path)
    args = build_parser().parse_args([
        "--ldsc-input", str(sumstats),
        "--merge-alleles", str(merge),
        "--ref-ld-chr", str(reference),
        "--w-ld-chr", str(weights),
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path / "output"),
        "--ldsc-population", "AFR",
        "--ldsc-genome-build", "GRCh37",
        "--info-min", "0.85",
        "--maf-min", "0.02",
        "--n-min", "700",
        "--chunksize", "2000",
        "--keep-maf",
        "--intercept-h2", "1.0",
        "--chisq-max", "90",
        "--n-blocks", "100",
        "--not-M-5-50",
        "--print-cov",
        "--print-delete-vals",
    ])
    module = resolve_ldsc_configuration(args).modules.ldsc
    assert module.population.value == "AFR"
    assert module.genome_build.value == "GRCh37"
    assert module.minimum_info == 0.85
    assert module.minimum_maf == 0.02
    assert module.minimum_n == 700
    assert module.chunksize == 2000
    assert module.keep_maf is True
    assert module.intercept == 1.0
    assert module.chisq_max == 90
    assert module.n_blocks == 100
    assert module.use_m_5_50 is False
    assert module.print_covariance is True
    assert module.print_delete_values is True


def test_command_builders_apply_every_supported_resolved_setting(tmp_path):
    munge = build_munge_sumstats_command(
        "munge_sumstats.py",
        sumstats=tmp_path / "input.tsv",
        output_prefix=tmp_path / "munged",
        merge_alleles=tmp_path / "hm3",
        minimum_info=0.8,
        minimum_maf=0.02,
        minimum_n=900,
        chunksize=1234,
        keep_maf=True,
    )
    assert munge[munge.index("--chunksize") + 1] == "1234"
    assert munge[munge.index("--n-min") + 1] == "900"
    assert "--keep-maf" in munge
    assert munge[munge.index("--info-min") + 1] == "0.8"
    assert munge[munge.index("--maf-min") + 1] == "0.02"

    h2 = build_h2_command(
        "ldsc.py",
        sumstats=tmp_path / "munged.sumstats.gz",
        reference_ld_directory=tmp_path / "reference",
        weights_ld_directory=tmp_path / "weights",
        output_prefix=tmp_path / "result",
        intercept=1.0,
        two_step=None,
        chisq_max=75,
        n_blocks=150,
        use_m_5_50=False,
        print_covariance=True,
        print_delete_values=True,
        sample_prevalence=0.2,
        population_prevalence=0.01,
    )
    assert h2[h2.index("--n-blocks") + 1] == "150"
    assert h2[h2.index("--intercept-h2") + 1] == "1.0"
    assert h2[h2.index("--chisq-max") + 1] == "75"
    assert "--not-M-5-50" in h2
    assert "--print-cov" in h2
    assert "--print-delete-vals" in h2
    assert h2[h2.index("--samp-prev") + 1] == "0.2"
    assert h2[h2.index("--pop-prev") + 1] == "0.01"
    assert h2[h2.index("--ref-ld-chr") + 1] != h2[h2.index("--w-ld-chr") + 1]


def test_default_commands_preserve_upstream_data_dependent_options(tmp_path):
    module = load_configuration().modules.ldsc
    munge = build_munge_sumstats_command(
        "munge_sumstats.py",
        sumstats=tmp_path / "input.tsv",
        output_prefix=tmp_path / "munged",
        merge_alleles=tmp_path / "hm3",
        minimum_info=module.minimum_info,
        minimum_maf=module.minimum_maf,
        minimum_n=module.minimum_n,
        chunksize=module.chunksize,
        keep_maf=module.keep_maf,
    )
    assert "--n-min" not in munge
    assert "--keep-maf" not in munge

    h2 = build_h2_command(
        "ldsc.py",
        sumstats=tmp_path / "munged.sumstats.gz",
        reference_ld_directory=tmp_path / "reference",
        weights_ld_directory=tmp_path / "weights",
        output_prefix=tmp_path / "result",
        intercept=module.intercept,
        two_step=module.two_step,
        chisq_max=module.chisq_max,
        n_blocks=module.n_blocks,
        use_m_5_50=module.use_m_5_50,
        print_covariance=module.print_covariance,
        print_delete_values=module.print_delete_values,
    )
    assert h2[h2.index("--n-blocks") + 1] == "200"
    for option in (
        "--intercept-h2", "--two-step", "--chisq-max", "--not-M-5-50",
        "--print-cov", "--print-delete-vals", "--samp-prev", "--pop-prev",
    ):
        assert option not in h2


def test_service_uses_formatter_returned_sample_prevalence(
    tmp_path, monkeypatch, capsys,
):
    args = _service_args(tmp_path, pop_prev=0.01)
    captured = {}

    def fake_run(**kwargs):
        captured["sample_prevalence"] = kwargs["configuration"].sample_prevalence
        return _fake_ldsc_result(**kwargs)

    monkeypatch.setattr("postgwas.modules.ldsc.service.run_ldsc", fake_run)
    ctx = {"formatter": {"ldsc": {
        "sample_prev": 0.25, "trait_type": "binary",
    }}}
    result = run_ldsc_direct(args, ctx)

    assert captured["sample_prevalence"] == 0.25
    assert result["sample_prevalence"] == 0.25
    assert result["sample_prevalence_source"] == "formatter_return_value"
    assert Path(result["h2_observed"]).is_file()
    assert Path(result["h2_liability"]).is_file()
    assert ctx["heritability"] == result
    resolved = yaml.safe_load(
        Path(result["resolved_configuration"]).read_text(encoding="utf-8")
    )
    assert resolved["modules"]["ldsc"]["sample_prevalence"] == 0.25
    assert resolved["modules"]["ldsc"]["population_prevalence"] == 0.01
    assert resolved["modules"]["ldsc"]["population"] is None
    assert resolved["modules"]["ldsc"]["genome_build"] is None

    screen = capsys.readouterr().out
    for finding in (
        "LDSC heritability summary",
        "Main findings",
        "Observed-scale h²",
        "0.20 (0.03)",
        "Observed intercept",
        "Observed ratio",
        "Liability-scale h²",
        "0.35 (0.05)",
        "Sample prevalence",
        "0.25 (formatter-returned value)",
        "Population prevalence",
        "Full PostGWAS log",
    ):
        assert finding in screen
    assert "Scientific findings" not in screen

    summary_lines = screen.splitlines()
    finding_labels = (
        "Observed-scale h²",
        "Observed intercept",
        "Observed ratio",
        "Liability-scale h²",
        "Liability intercept",
        "Liability ratio",
        "Sample prevalence",
        "Population prevalence",
    )
    finding_lines = {
        label: next(line for line in summary_lines if label in line)
        for label in finding_labels
    }
    assert len({line.index(":") for line in finding_lines.values()}) == 1
    observed_ratio_index = summary_lines.index(finding_lines["Observed ratio"])
    liability_ratio_index = summary_lines.index(finding_lines["Liability ratio"])
    assert summary_lines[observed_ratio_index + 1] == ""
    assert summary_lines[liability_ratio_index + 1] == ""
    service_log = Path(result["service_log"]).read_text(encoding="utf-8")
    assert "RESULT   ldsc_observed_scale" in service_log
    assert "RESULT   ldsc_liability_scale" in service_log


def test_explicit_sample_prevalence_overrides_formatter_return(tmp_path, monkeypatch):
    args = _service_args(tmp_path, samp_prev=0.3, pop_prev=0.01)
    monkeypatch.setattr(
        "postgwas.modules.ldsc.service.run_ldsc", _fake_ldsc_result,
    )
    result = run_ldsc_direct(args, {
        "formatter": {"ldsc": {"sample_prev": 0.25, "trait_type": "binary"}},
    })
    assert result["sample_prevalence"] == 0.3
    assert result["sample_prevalence_source"] == "configuration_or_cli"


def test_pipeline_runner_passes_formatter_result_without_calculating_prevalence(
    tmp_path, monkeypatch,
):
    from postgwas.pipeline.runners import run_heritability_runner

    captured = {}

    def fake_service(args, ctx):
        captured["ctx"] = ctx
        captured["has_sample_override"] = hasattr(args, "samp_prev")
        captured["input"] = args.ldsc_input
        return {"sample_prevalence_source": "formatter_return_value"}

    monkeypatch.setattr(
        "postgwas.modules.ldsc.service.run_ldsc_direct", fake_service,
    )
    args = Namespace(
        output_directory=str(tmp_path), dataset_id="STUDY", _step_num="04",
    )
    ctx = {"formatter": {"ldsc": {
        "ldsc_file": "formatter-returned.tsv", "sample_prev": 0.25,
    }}}
    result = run_heritability_runner(args, ctx)

    assert result == {"sample_prevalence_source": "formatter_return_value"}
    assert captured == {
        "ctx": ctx,
        "has_sample_override": False,
        "input": "formatter-returned.tsv",
    }
    assert args.output_directory == str(tmp_path)


def test_population_prevalence_without_formatter_value_fails(tmp_path):
    args = _service_args(tmp_path, pop_prev=0.01)
    with pytest.raises(LDSCError, match="formatter returned no sample prevalence"):
        run_ldsc_direct(args, {
            "formatter": {"ldsc": {"sample_prev": None, "trait_type": "quantitative"}},
        })


def test_formatter_returned_prevalence_is_schema_validated(tmp_path):
    args = _service_args(tmp_path, pop_prev=0.01)
    with pytest.raises(ValidationError, match="sample_prevalence"):
        run_ldsc_direct(args, {
            "formatter": {"ldsc": {"sample_prev": 1.2, "trait_type": "binary"}},
        })


def test_sample_prevalence_without_population_runs_observed_only(
    tmp_path, monkeypatch, capsys,
):
    args = _service_args(tmp_path, samp_prev=0.2)
    monkeypatch.setattr(
        "postgwas.modules.ldsc.service.run_ldsc", _fake_ldsc_result,
    )

    result = run_ldsc_direct(args)

    assert Path(result["h2_observed"]).is_file()
    assert result["h2_liability"] is None
    assert result["liability_outputs"] == []
    assert result["sample_prevalence"] == 0.2
    assert result["population_prevalence"] is None
    screen = capsys.readouterr().out
    assert "Observed-scale h²" in screen
    assert "not run because population prevalence was not provided" in screen
    assert "Liability intercept" not in screen


def test_constrained_intercept_rejects_two_step_estimator():
    with pytest.raises(ConfigurationError, match="intercept and two_step"):
        resolve_ldsc_configuration(Namespace(
            ldsc_intercept=1.0, ldsc_two_step=30,
        ))


@pytest.mark.parametrize("name,value", [
    ("samp_prev", 0), ("samp_prev", 1), ("pop_prev", -0.1),
    ("pop_prev", 1.1),
])
def test_prevalence_must_be_strictly_between_zero_and_one(name, value):
    with pytest.raises(ConfigurationError, match="prevalence"):
        resolve_ldsc_configuration(Namespace(**{name: value}))


def test_success_without_required_outputs_is_rejected(tmp_path):
    sumstats, merge, reference, weights = _write_inputs(tmp_path)
    no_op = tmp_path / "no-op"
    no_op.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    no_op.chmod(0o755)
    with pytest.raises(LDSCError, match="expected output is missing or empty"):
        run_ldsc(
            sumstats_tsv=sumstats,
            output_prefix=tmp_path / "result",
            merge_alleles=merge,
            reference_ld_directory=reference,
            weights_ld_directory=weights,
            munge_executable=str(no_op),
            ldsc_executable=str(no_op),
            configuration=load_configuration().modules.ldsc,
        )


def test_runner_skips_liability_without_population_prevalence(
    tmp_path, monkeypatch,
):
    sumstats, merge, reference, weights = _write_inputs(tmp_path)
    module = load_configuration().modules.ldsc.model_copy(update={
        "sample_prevalence": 0.2,
    })
    layout = module.output_layout
    commands = []
    events = []

    class RecordingLogger:
        def record(self, *values, **details):
            events.append((values, details))

    def fake_subprocess_run(command, **_kwargs):
        commands.append(command)
        output_prefix = Path(command[command.index("--out") + 1])
        if "--sumstats" in command:
            Path(str(output_prefix) + layout.munged_sumstats_suffix).write_text(
                "sumstats\n", encoding="utf-8",
            )
            Path(str(output_prefix) + layout.upstream_log_suffix).write_text(
                "munge complete\n", encoding="utf-8",
            )
        else:
            Path(str(output_prefix) + layout.upstream_log_suffix).write_text(
                "Total Observed scale h2: 0.20 (0.03)\n"
                "Intercept: 1.01 (0.01)\n",
                encoding="utf-8",
            )
        return Namespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(
        "postgwas.core.processes.subprocess.run", fake_subprocess_run,
    )
    result = run_ldsc(
        sumstats_tsv=sumstats,
        output_prefix=tmp_path / "result",
        merge_alleles=merge,
        reference_ld_directory=reference,
        weights_ld_directory=weights,
        munge_executable="munge_sumstats.py",
        ldsc_executable="ldsc.py",
        configuration=module,
        logger=RecordingLogger(),
    )

    assert len(commands) == 2
    assert result["h2_liability"] is None
    assert result["liability_outputs"] == []
    assert (("SKIP", "ldsc_liability_scale"), {
        "reason": "population_prevalence_not_provided",
    }) in events


@pytest.mark.parametrize("missing_stage", ["observed", "liability"])
def test_zero_exit_missing_heritability_log_is_rejected(
    tmp_path, monkeypatch, missing_stage,
):
    sumstats, merge, reference, weights = _write_inputs(tmp_path)
    module = load_configuration().modules.ldsc
    if missing_stage == "liability":
        module = module.model_copy(update={
            "sample_prevalence": 0.2,
            "population_prevalence": 0.01,
        })
    layout = module.output_layout

    def fake_subprocess_run(command, **_kwargs):
        output_prefix = Path(command[command.index("--out") + 1])
        if "--sumstats" in command:
            Path(str(output_prefix) + layout.munged_sumstats_suffix).write_text(
                "sumstats\n", encoding="utf-8",
            )
            Path(str(output_prefix) + layout.upstream_log_suffix).write_text(
                "munge complete\n", encoding="utf-8",
            )
        elif (
            missing_stage == "liability"
            and not output_prefix.name.endswith(layout.liability_prefix_suffix)
        ):
            Path(str(output_prefix) + layout.upstream_log_suffix).write_text(
                "Total Observed scale h2: 0.20 (0.03)\n"
                "Intercept: 1.01 (0.01)\n",
                encoding="utf-8",
            )
        return Namespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(
        "postgwas.core.processes.subprocess.run", fake_subprocess_run,
    )
    with pytest.raises(
        LDSCError,
        match=(
            "LDSC %s-scale heritability returned success but expected output "
            "is missing or empty" % missing_stage
        ),
    ):
        run_ldsc(
            sumstats_tsv=sumstats,
            output_prefix=tmp_path / "result",
            merge_alleles=merge,
            reference_ld_directory=reference,
            weights_ld_directory=weights,
            munge_executable="munge_sumstats.py",
            ldsc_executable="ldsc.py",
            configuration=module,
        )


def test_metrics_parser_requires_scientific_results(tmp_path):
    log = tmp_path / "result.log"
    log.write_text("LDSC completed\n", encoding="utf-8")
    with pytest.raises(LDSCError, match="missing required result metric"):
        extract_ldsc_metrics(log)


def test_metrics_parser_accepts_ratio_na(tmp_path):
    log = tmp_path / "result.log"
    log.write_text(
        "Total Observed scale h2: 0.20 (0.03)\n"
        "Intercept: 1.01 (0.01)\nRatio: NA\n",
        encoding="utf-8",
    )
    assert extract_ldsc_metrics(log) == {
        "h2": "0.20 (0.03)",
        "intercept": "1.01 (0.01)",
        "ratio": "NA",
    }


def test_metrics_parser_accepts_constrained_intercept(tmp_path):
    log = tmp_path / "result.log"
    log.write_text(
        "Total Observed scale h2: 0.20 (0.03)\n"
        "Intercept: constrained to 1\n",
        encoding="utf-8",
    )
    assert extract_ldsc_metrics(log)["intercept"] == "constrained to 1"


def test_metrics_parser_preserves_negative_ratio_message(tmp_path):
    log = tmp_path / "result.log"
    log.write_text(
        "Total Observed scale h2: 0.20 (0.03)\n"
        "Intercept: 0.99 (0.01)\n"
        "Ratio < 0 (usually indicates GC correction).\n",
        encoding="utf-8",
    )
    assert extract_ldsc_metrics(log)["ratio"].startswith("Ratio < 0")


def test_metrics_parser_rejects_nonfinite_scientific_result(tmp_path):
    log = tmp_path / "result.log"
    log.write_text(
        "Total Observed scale h2: nan (nan)\n"
        "Intercept: 1.01 (0.01)\n",
        encoding="utf-8",
    )
    with pytest.raises(LDSCError, match="heritability"):
        extract_ldsc_metrics(log)


def test_dockerfile_uses_cbiit_ldsc_in_postgwas_environment():
    dockerfile = (Path(__file__).parents[1] / "Dockerfile").read_text(
        encoding="utf-8",
    )
    assert "https://github.com/CBIIT/ldsc.git" in dockerfile
    assert "LDSC_COMMIT=6c673952cee74bd5c57aef1555a03b1c015399a0" in dockerfile
    assert "pip install --no-deps --no-cache-dir /opt/ldsc" in dockerfile
    assert "python=2.7" not in dockerfile
    assert "micromamba create -y -n ldsc" not in dockerfile
    assert "github.com/bulik/ldsc" not in dockerfile
