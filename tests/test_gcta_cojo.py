import argparse
from argparse import Namespace
from pathlib import Path
import shutil
import subprocess

import polars as pl
import pytest

from postgwas.config import load_configuration, load_module_configuration
from postgwas.core.contracts import RunContext
from postgwas.core.errors import ConfigurationError
from postgwas.modules.gcta_cojo.adapters import (
    build_cojo_command,
    require_supported_gcta,
)
from postgwas.modules.gcta_cojo.cli import build_parser
from postgwas.modules.gcta_cojo.errors import GctaCojoError
from postgwas.modules.gcta_cojo.results import normalize_cojo_results
from postgwas.modules.gcta_cojo.service import (
    _resolved_configuration,
    run_gcta_cojo_direct,
)
from postgwas.pipeline.runners import run_formatter_runner, run_gcta_cojo_runner


def _reference(tmp_path, *, variant="rs1", samples=4):
    prefix = tmp_path / "reference"
    Path(str(prefix) + ".bed").write_bytes(b"BED")
    Path(str(prefix) + ".bim").write_text(
        "1\t%s\t0\t100\tG\tA\n" % variant,
        encoding="utf-8",
    )
    Path(str(prefix) + ".fam").write_text(
        "".join("F%d I%d 0 0 0 -9\n" % (index, index) for index in range(samples)),
        encoding="utf-8",
    )
    return prefix


def _summary(tmp_path, *, variant="rs1"):
    path = tmp_path / "input.ma"
    path.write_text(
        "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\n"
        "%s\tG\tA\t0.2\t0.1\t0.05\t1e-9\t1000\n" % variant,
        encoding="utf-8",
    )
    return path


def _fake_gcta(tmp_path, *, version="1.94.1"):
    path = tmp_path / "fake_gcta"
    path.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        " print('GCTA v%s')\n"
        " raise SystemExit(0)\n"
        "prefix = sys.argv[sys.argv.index('--out') + 1]\n"
        "summary = sys.argv[sys.argv.index('--cojo-file') + 1]\n"
        "variant = Path(summary).read_text().splitlines()[1].split()[0]\n"
        "Path(prefix + '.log').write_text('Analysis finished\\n')\n"
        "if '--cojo-cond' in sys.argv:\n"
        " Path(prefix + '.cma.cojo').write_text("
        "'Chr SNP bp freq refA b se p n freq_geno bC bC_se pC\\n'"
        "+ '1 ' + variant + ' 100 0.2 G 0.1 0.05 1e-9 1000 0.2 0.08 0.04 0.04\\n')\n"
        "else:\n"
        " Path(prefix + '.jma.cojo').write_text("
        "'Chr SNP bp freq refA b se p n freq_geno bJ bJ_se pJ LD_r\\n'"
        "+ '1 ' + variant + ' 100 0.2 G 0.1 0.05 1e-9 1000 0.2 0.09 0.04 0.02 0\\n')\n"
        " Path(prefix + '.jma.ldr.cojo').write_text("
        "'SNP ' + variant + '\\n' + variant + ' 1\\n')\n"
        % version,
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _config(tmp_path, reference, *, mode="slct", condition=None):
    path = tmp_path / "gcta_cojo.yaml"
    lines = [
        "enabled: true",
        "mode: %s" % mode,
        "genome_build: GRCh37",
        "reference:",
        "  prefix: '%s'" % reference,
        "  population: EUR",
        "input_validation:",
        "  reference_sample_size_warning_threshold: 1",
    ]
    if condition is not None:
        lines.extend(["inputs:", "  condition_snps: '%s'" % condition])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _args(tmp_path, *, mode="slct", condition=None, version="1.94.1"):
    reference = _reference(tmp_path)
    return Namespace(
        gcta_cojo_input_file=str(_summary(tmp_path)),
        dataset_id="STUDY",
        output_directory=str(tmp_path / "out"),
        run_config=str(_config(tmp_path, reference, mode=mode, condition=condition)),
        gcta=str(_fake_gcta(tmp_path, version=version)),
        threads=2,
        resume=False,
        overwrite=False,
        dry_run=False,
    )


def test_cli_exposes_every_documented_cojo_control_without_argparse_defaults():
    parser = build_parser()
    help_text = parser.format_help()
    for option in (
        "--cojo-file",
        "--cojo-mode", "--cojo-p", "--cojo-top-snps", "--cojo-window-kb",
        "--cojo-collinear", "--cojo-diff-freq", "--cojo-maf", "--cojo-gc",
        "--cojo-gc-lambda", "--condition-snps", "--joint-snps",
        "--cojo-extract", "--cojo-exclude", "--cojo-chromosome",
    ):
        assert option in help_text
    assert "--vcf" not in help_text
    compact_help = " ".join(help_text.split())
    assert "Configured YAML default: slct" in compact_help
    assert "Known SNPs to adjust for" in compact_help
    assert "any explicit mode must be cond" in compact_help
    assert "effects are estimated together in one multiple-SNP model" in compact_help
    assert "any explicit mode must be joint" in compact_help
    assert "Only test or select SNPs listed in this file" in compact_help
    assert "Do not test or select SNPs listed in this file" in compact_help
    configurable_destinations = {
        "cojo_chromosome", "cojo_collinear", "cojo_diff_freq",
        "cojo_exclude", "cojo_extract", "cojo_finding_threshold", "cojo_gc",
        "cojo_gc_lambda", "cojo_maf", "cojo_minimum_reference_overlap",
        "cojo_mode", "cojo_p", "cojo_p_value_digits",
        "cojo_reference_population", "cojo_reference_prefix", "cojo_top_results",
        "cojo_top_snps", "cojo_window_kb", "condition_snps", "dataset_id",
        "dry_run", "gcta", "gcta_cojo_input_file", "genome_build", "joint_snps",
        "memory_gb",
        "output_directory", "overwrite", "resume", "run_config", "seed",
        "threads",
    }
    observed = {action.dest for action in parser._actions if action.dest != "help"}
    assert observed == configurable_destinations
    for action in parser._actions:
        if action.dest in configurable_destinations:
            assert action.default == argparse.SUPPRESS


def test_condition_list_implies_cond_mode(tmp_path):
    condition = tmp_path / "condition.txt"
    condition.write_text("rs1\n", encoding="utf-8")
    configuration = _resolved_configuration(Namespace(
        condition_snps=str(condition),
    ))
    assert configuration.modules.gcta_cojo.mode == "cond"
    assert configuration.modules.gcta_cojo.inputs.condition_snps == str(condition)


@pytest.mark.parametrize(
    ("mode", "mode_flag"),
    (
        ("slct", "--cojo-slct"),
        ("top_snps", "--cojo-top-SNPs"),
        ("joint", "--cojo-joint"),
        ("cond", "--cojo-cond"),
    ),
)
def test_command_builder_covers_all_cojo_modes(tmp_path, mode, mode_flag):
    snps = tmp_path / "snps.txt"
    snps.write_text("rs1\n", encoding="utf-8")
    overrides = {"modules.gcta_cojo.mode": mode}
    if mode == "joint":
        overrides["modules.gcta_cojo.inputs.joint_snps"] = str(snps)
    if mode == "cond":
        overrides["modules.gcta_cojo.inputs.condition_snps"] = str(snps)
    module = load_configuration(cli_overrides=overrides).modules.gcta_cojo
    command = build_cojo_command(
        "gcta64", "study.ma", "reference", "output", 4, module,
    )
    assert mode_flag in command
    for option in (
        "--cojo-file", "--bfile", "--maf", "--diff-freq", "--cojo-wind",
        "--cojo-collinear", "--thread-num", "--out",
    ):
        assert option in command
    assert ("--cojo-p" in command) is (mode == "slct")


def test_version_probe_isolates_and_removes_gcta_log(monkeypatch):
    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        output_prefix = Path(
            command[command.index("--out") + 1]
        )
        captured["output_prefix"] = output_prefix
        output_prefix.with_suffix(".log").write_text(
            "version probe\n", encoding="utf-8",
        )
        return subprocess.CompletedProcess(
            command, 1, stdout="GCTA v1.95.3\n", stderr="invalid option\n",
        )

    class Logger:
        def record(self, *args, **kwargs):
            captured["record"] = (args, kwargs)

    monkeypatch.setattr("postgwas.core.gcta.subprocess.run", fake_run)
    module = load_configuration().modules.gcta_cojo

    version = require_supported_gcta("gcta64", module, Logger(), 10)

    assert version == "1.95.3"
    assert module.version_probe_output_argument in captured["command"]
    assert not captured["output_prefix"].parent.exists()


def test_conditional_na_estimates_are_preserved_and_labelled(tmp_path):
    condition = tmp_path / "condition.txt"
    condition.write_text("rs1\n", encoding="utf-8")
    module = load_configuration(cli_overrides={
        "modules.gcta_cojo.mode": "cond",
        "modules.gcta_cojo.inputs.condition_snps": str(condition),
    }).modules.gcta_cojo
    source = tmp_path / "result.cma.cojo"
    source.write_text(
        "Chr SNP bp freq refA b se p n freq_geno bC bC_se pC\n"
        "1 rs1 100 0.2 G 0.1 0.05 1e-9 1000 0.2 NA NA NA\n",
        encoding="utf-8",
    )
    destination = tmp_path / "normalized.tsv"
    metrics, _ = normalize_cojo_results(source, destination, module)
    normalized = pl.read_csv(destination, separator="\t", null_values="NA")
    assert metrics["not_estimable_rows"] == 1
    assert normalized["estimation_status"].item() == "not_estimable_collinearity"
    assert normalized["pC"].item() is None


@pytest.mark.parametrize(
    ("replacement", "message"),
    (
        ("invalid invalid invalid", "non-numeric"),
        ("0.08 -0.04 0.04", "standard errors"),
        ("0.08 0.04 1.4", "invalid p-values"),
    ),
)
def test_conditional_invalid_adjusted_statistics_fail(tmp_path, replacement, message):
    condition = tmp_path / "condition.txt"
    condition.write_text("rs1\n", encoding="utf-8")
    module = load_configuration(cli_overrides={
        "modules.gcta_cojo.mode": "cond",
        "modules.gcta_cojo.inputs.condition_snps": str(condition),
    }).modules.gcta_cojo
    source = tmp_path / "result.cma.cojo"
    source.write_text(
        "Chr SNP bp freq refA b se p n freq_geno bC bC_se pC\n"
        "1 rs1 100 0.2 G 0.1 0.05 1e-9 1000 0.2 %s\n" % replacement,
        encoding="utf-8",
    )
    with pytest.raises(GctaCojoError, match=message):
        normalize_cojo_results(source, tmp_path / "normalized.tsv", module)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("bp", "100.5", "non-integer"),
        ("se", "0", "non-positive"),
        ("freq", "1.2", "outside"),
    ),
)
def test_marginal_result_constraints_come_from_yaml(tmp_path, field, value, message):
    module = load_configuration().modules.gcta_cojo
    columns = {
        "Chr": "1", "SNP": "rs1", "bp": "100", "freq": "0.2",
        "refA": "G", "b": "0.1", "se": "0.05", "p": "1e-9",
        "n": "1000", "freq_geno": "0.2", "bJ": "0.08",
        "bJ_se": "0.04", "pJ": "0.04", "LD_r": "0",
    }
    columns[field] = value
    source = tmp_path / "result.jma.cojo"
    source.write_text(
        " ".join(columns) + "\n" + " ".join(columns.values()) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(GctaCojoError, match=message):
        normalize_cojo_results(source, tmp_path / "normalized.tsv", module)


def test_slct_run_validates_normalizes_logs_and_summarizes(tmp_path, capsys):
    result = run_gcta_cojo_direct(_args(tmp_path))
    terminal = capsys.readouterr().out
    assert result.metrics["mode"] == "slct"
    assert result.metrics["reference_overlap_fraction"] == 1.0
    assert result.metrics["variant_id_type"] == "rsid"
    assert result.artifacts["normalized_results"].path.is_file()
    assert result.artifacts["ld_matrix"].path.is_file()
    assert result.artifacts["completion_manifest"].path.is_file()
    assert "GCTA-COJO completed" in terminal
    assert "Independent/joint signals" in terminal
    assert "BIM identifier type" in terminal
    log = tmp_path / "out" / "logs" / "STUDY_slct_gcta_cojo.log"
    assert "gcta_cojo status=COMPLETED" in log.read_text(encoding="utf-8")


def test_conditional_run_reports_findings_without_joint_ld_artifact(tmp_path, capsys):
    condition = tmp_path / "condition.txt"
    condition.write_text("rs1\n", encoding="utf-8")

    result = run_gcta_cojo_direct(
        _args(tmp_path, mode="cond", condition=condition)
    )

    terminal = capsys.readouterr().out
    assert result.metrics["mode"] == "cond"
    assert result.metrics["estimated_rows"] == 1
    assert "ld_matrix" not in result.artifacts
    assert result.artifacts["raw_results"].path.name.endswith(".cma.cojo")
    assert "Conditional findings" in terminal


def test_completed_run_resumes_only_from_validated_manifest(tmp_path):
    args = _args(tmp_path)
    first = run_gcta_cojo_direct(args)
    args.resume = True

    resumed = run_gcta_cojo_direct(args)

    assert first.artifacts["normalized_results"].path == (
        resumed.artifacts["normalized_results"].path
    )
    assert resumed.metrics["resumed"] is True


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_pipeline_vcf_mode_formats_then_runs_cojo_with_bim_identifier_type(tmp_path):
    prefix = tmp_path / "unique_reference"
    Path(str(prefix) + ".bed").write_bytes(b"BED")
    Path(str(prefix) + ".fam").write_text(
        "F1 I1 0 0 0 -9\nF2 I2 0 0 0 -9\n",
        encoding="utf-8",
    )
    Path(str(prefix) + ".bim").write_text(
        "1\t1_100_A_G\t0\t100\tG\tA\n"
        "1\t1_200_C_T\t0\t200\tT\tC\n"
        "2\t2_300_G_A\t0\t300\tA\tG\n"
        "2\t2_400_T_C\t0\t400\tC\tT\n",
        encoding="utf-8",
    )
    fixture = Path(__file__).parent / "fixtures" / "formatting" / "harmonised.vcf"
    output = tmp_path / "out"
    args = Namespace(
        modules=["gcta_cojo"],
        vcf=str(fixture),
        cojo_reference_prefix=str(prefix),
        cojo_reference_population="EUR",
        genome_build="GRCh37",
        dataset_id="STUDY",
        output_directory=str(output),
        gcta=str(_fake_gcta(tmp_path)),
        bcftools=shutil.which("bcftools"),
        threads=2,
        resume=False,
        overwrite=False,
        dry_run=False,
        _step_num="01",
    )
    context = RunContext()

    run_formatter_runner(args, context)
    args._step_num = "02"
    result = run_gcta_cojo_runner(args, context)

    assert result.metrics["variant_id_type"] == "unique"
    assert result.metrics["summary_variants"] == 3
    formatted = output / "01_formatter" / "STUDY_gcta.ma"
    identifiers = pl.read_csv(formatted, separator="\t")["SNP"].to_list()
    assert identifiers == ["1_100_A_G", "1_200_C_T", "2_300_G_A"]
    assert args.variant_id_types["gcta_gene"] == "unique"
    assert result.artifacts["normalized_results"].path.parent == (
        output / "02_gcta_cojo" / "results"
    )


def test_direct_mode_requires_existing_cojo_file(tmp_path, capsys):
    reference = _reference(tmp_path)
    args = Namespace(
        dataset_id="STUDY",
        output_directory=str(tmp_path / "out"),
        cojo_reference_prefix=str(reference),
        cojo_reference_population="EUR",
        genome_build="GRCh37",
        gcta=str(_fake_gcta(tmp_path)),
        threads=2,
        resume=False,
        overwrite=False,
        dry_run=False,
    )

    with pytest.raises(GctaCojoError, match="--cojo-file"):
        run_gcta_cojo_direct(args)

    terminal = capsys.readouterr().out
    assert "GCTA-COJO failed" in terminal
    assert "postgwas pipeline --modules gcta_cojo" in terminal


def test_conditioning_snp_must_exist_in_summary_and_bim(tmp_path, capsys):
    condition = tmp_path / "condition.txt"
    condition.write_text("rs_missing\n", encoding="utf-8")
    args = _args(tmp_path, mode="cond", condition=condition)
    with pytest.raises(GctaCojoError, match="incompatible with the GCTA .ma"):
        run_gcta_cojo_direct(args)
    terminal = capsys.readouterr().out
    assert "GCTA-COJO failed" in terminal
    assert "Reason" in terminal
    log = tmp_path / "out" / "logs" / "STUDY_cond_gcta_cojo.log"
    assert "FAILED" in log.read_text(encoding="utf-8")


def test_old_gcta_version_fails_with_terminal_reason_and_file_log(tmp_path, capsys):
    with pytest.raises(GctaCojoError, match="too old"):
        run_gcta_cojo_direct(_args(tmp_path, version="1.93.3"))
    terminal = capsys.readouterr().out
    assert "GCTA-COJO failed" in terminal
    assert "too old" in terminal
    log = tmp_path / "out" / "logs" / "STUDY_slct_gcta_cojo.log"
    assert "FAILED" in log.read_text(encoding="utf-8")


def test_module_only_configuration_is_schema_validated(tmp_path):
    condition = tmp_path / "condition.txt"
    condition.write_text("rs1\n", encoding="utf-8")
    path = tmp_path / "cojo.yaml"
    path.write_text(
        "mode: cond\ninputs:\n  condition_snps: '%s'\n" % condition,
        encoding="utf-8",
    )
    assert load_module_configuration("gcta_cojo", path).mode == "cond"


@pytest.mark.parametrize(
    ("key", "value"),
    (
        ("modules.gcta_cojo.results.delimiter", "unsupported"),
        ("modules.gcta_cojo.results.delimiter_candidates", ["unsupported"]),
        ("modules.gcta_cojo.results.atomic_output_suffix", "../temporary"),
    ),
)
def test_result_io_configuration_fails_during_schema_validation(key, value):
    with pytest.raises(ConfigurationError):
        load_configuration(cli_overrides={key: value})


def test_pipeline_runner_consumes_formatter_artifact_and_restores_root(tmp_path):
    args = _args(tmp_path)
    root = tmp_path / "pipeline"
    args.output_directory = str(root)
    args._step_num = "02"
    summary = Path(args.gcta_cojo_input_file)
    context = RunContext({
        "formatter": {
            "gcta_gene": {"summary_statistics_input_file": str(summary)},
        },
    })

    result = run_gcta_cojo_runner(args, context)

    assert result.module == "gcta_cojo"
    assert context["gcta_cojo"] is result
    assert args.output_directory == str(root)
    assert result.artifacts["normalized_results"].path.parent == (
        root / "02_gcta_cojo" / "results"
    )


def test_pipeline_formatter_uses_cojo_bim_identifier_type(
    tmp_path,
    monkeypatch,
):
    reference = _reference(tmp_path, variant="1_100_A_G")
    captured = {}

    def fake_formatter(args, ctx):
        captured["target_types"] = dict(args.variant_id_types)
        captured["observed"] = dict(args.variant_id_observations["gcta_gene"])
        captured["formats"] = list(args.format)
        return {"gcta_gene": {"summary_statistics_input_file": "study.ma"}}

    monkeypatch.setattr(
        "postgwas.modules.formatting.service.run_formatter_direct",
        fake_formatter,
    )
    fixture = Path(__file__).parent / "fixtures" / "formatting" / "harmonised.vcf"
    output = tmp_path / "pipeline"
    args = Namespace(
        modules=["gcta_cojo"],
        vcf=str(fixture),
        cojo_reference_prefix=str(reference),
        output_directory=str(output),
        dataset_id="STUDY",
        bcftools=shutil.which("bcftools") or shutil.which("python"),
        _step_num="01",
    )

    result = run_formatter_runner(args, {})

    assert result["gcta_gene"]["summary_statistics_input_file"] == "study.ma"
    assert captured["target_types"] == {"gcta_gene": "unique"}
    assert captured["observed"]["unique_ids"] == 1
    assert captured["formats"] == ["gcta_gene"]
    assert args.output_directory == str(output)
