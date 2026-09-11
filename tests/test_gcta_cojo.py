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
from postgwas.core.screen_logging import ScreenSettings, record_screen
from postgwas.modules.gcta_cojo.adapters import (
    build_cojo_command,
    require_supported_gcta,
)
from postgwas.modules.gcta_cojo.cli import build_parser
from postgwas.modules.gcta_cojo.errors import GctaCojoError
from postgwas.modules.gcta_cojo.results import normalize_cojo_results
from postgwas.modules.gcta_cojo.parallel import (
    numeric_chromosome_workloads,
    plan_chromosome_workers,
    should_parallelize_slct,
)
from postgwas.modules.gcta_cojo.service import (
    GctaCojoPipelineResources,
    gcta_cojo_analysis_label,
    gcta_cojo_pipeline_title,
    preflight_gcta_cojo_pipeline,
    resolve_gcta_cojo_configuration,
    run_gcta_cojo_direct,
)
from postgwas.pipeline.registry import REGISTRY, resolve_reference
from postgwas.pipeline.runners import run_formatter_runner, run_gcta_cojo_runner
from preflight_support import pipeline_input_vcf_evidence


def _reference(tmp_path, *, variant="rs1", samples=4):
    prefix = tmp_path / "reference"
    # One SNP, four two-bit genotypes per byte, plus the SNP-major header.
    Path(str(prefix) + ".bed").write_bytes(bytes((0x6C, 0x1B, 0x01)) + bytes((samples + 3) // 4))
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


def _replace_validation_inputs(args, reference, summary_rows, bim_rows):
    Path(args.gcta_cojo_input_file).write_text(
        "SNP A1 A2 freq b se p N\n"
        + "".join(
            "%s %s %s 0.2 0.1 0.05 1e-9 1000\n" % row
            for row in summary_rows
        ),
        encoding="utf-8",
    )
    Path(str(reference) + ".bim").write_text(
        "".join(
            "1\t%s\t0\t%d\t%s\t%s\n"
            % (variant, position, allele1, allele2)
            for position, (variant, allele1, allele2) in enumerate(bim_rows, 100)
        ),
        encoding="utf-8",
    )


def _fake_gcta(
    tmp_path,
    *,
    version="1.94.1",
    no_signals=False,
    replace_modeled_snp=False,
):
    suffix = "_no_signals" if no_signals else ""
    suffix += "_replaced_snp" if replace_modeled_snp else ""
    path = tmp_path / ("fake_gcta" + suffix)
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
        "gc = '--cojo-gc' in sys.argv\n"
        "p = 'p_GC' if gc else 'p'\n"
        "pj = 'pJ_GC' if gc else 'pJ'\n"
        "pc = 'pC_GC' if gc else 'pC'\n"
        "no_signals = %r\n"
        "replace_modeled_snp = %r\n"
        "if no_signals:\n"
        " Path(prefix + '.log').write_text('No SNPs have been selected.\\n')\n"
        " raise SystemExit(0)\n"
        "Path(prefix + '.log').write_text('Analysis finished\\n')\n"
        "if '--cojo-cond' in sys.argv:\n"
        " requested = Path(sys.argv[sys.argv.index('--cojo-cond') + 1]).read_text().split()[0]\n"
        " modeled = 'rs_replaced' if replace_modeled_snp else requested\n"
        " Path(prefix + '.given.cojo').write_text('SNP\\n' + modeled + '\\n')\n"
        " Path(prefix + '.cma.cojo').write_text("
        "'Chr SNP bp freq refA b se ' + p + ' n freq_geno bC bC_se ' + pc + '\\n'"
        "+ '1 ' + variant + ' 100 0.2 G 0.1 0.05 1e-9 1000 0.2 0.08 0.04 0.04\\n')\n"
        "else:\n"
        " modeled = variant\n"
        " if '--cojo-joint' in sys.argv:\n"
        "  requested = Path(sys.argv[sys.argv.index('--extract') + 1]).read_text().split()[0]\n"
        "  modeled = 'rs_replaced' if replace_modeled_snp else requested\n"
        " Path(prefix + '.jma.cojo').write_text("
        "'Chr SNP bp freq refA b se ' + p + ' n freq_geno bJ bJ_se ' + pj + ' LD_r\\n'"
        "+ '1 ' + modeled + ' 100 0.2 G 0.1 0.05 1e-9 1000 0.2 0.09 0.04 0.02 0\\n')\n"
        " Path(prefix + '.ldr.cojo').write_text("
        "'SNP ' + modeled + '\\n' + modeled + ' 1\\n')\n"
        " if '--cojo-joint' not in sys.argv:\n"
        "  Path(prefix + '.cma.cojo').write_text("
        "'Chr SNP bp freq refA b se ' + p + ' n freq_geno bC bC_se ' + pc + '\\n'"
        "+ '1 ' + modeled + ' 100 0.2 G 0.1 0.05 1e-9 1000 0.2 0.08 0.04 0.04\\n')\n"
        % (version, no_signals, replace_modeled_snp),
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _parallel_fake_gcta(
    tmp_path,
    *,
    fail_once_chromosome=None,
    always_fail_chromosome=None,
    no_signal_chromosomes=(),
    emit_selection_cma=True,
):
    path = tmp_path / (
        "parallel_fake_gcta_%s_%s_%s_%s"
        % (
            fail_once_chromosome or "none",
            always_fail_chromosome or "none",
            "_".join(str(value) for value in no_signal_chromosomes) or "signals",
            "selection_cma" if emit_selection_cma else "no_selection_cma",
        )
    )
    counter = tmp_path / (path.name + ".attempts")
    path.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        " print('GCTA v1.94.1')\n"
        " raise SystemExit(0)\n"
        "prefix = sys.argv[sys.argv.index('--out') + 1]\n"
        "bfile = sys.argv[sys.argv.index('--bfile') + 1]\n"
        "rows = [line.split() for line in Path(bfile + '.bim').read_text().splitlines() if line.strip()]\n"
        "chromosome = int(sys.argv[sys.argv.index('--chr') + 1]) if '--chr' in sys.argv else None\n"
        "fail_once = %r\n"
        "always_fail = %r\n"
        "no_signal = set(%r)\n"
        "emit_selection_cma = %r\n"
        "counter = Path(%r)\n"
        "if chromosome is not None:\n"
        " previous = int(counter.read_text()) if counter.exists() else 0\n"
        " if chromosome == fail_once and previous == 0:\n"
        "  counter.write_text('1')\n"
        "  print('transient chromosome failure', file=sys.stderr)\n"
        "  raise SystemExit(7)\n"
        " if chromosome == always_fail:\n"
        "  counter.write_text(str(previous + 1))\n"
        "  print('persistent chromosome failure', file=sys.stderr)\n"
        "  raise SystemExit(8)\n"
        " selected_rows = [row for row in rows if chromosome is None or int(row[0]) == chromosome]\n"
        "Path(prefix + '.log').write_text('Analysis finished\\n')\n"
        "if '--cojo-cond' in sys.argv:\n"
        " selected = [line.split()[0] for line in Path(sys.argv[sys.argv.index('--cojo-cond') + 1]).read_text().splitlines() if line.strip()]\n"
        " Path(prefix + '.given.cojo').write_text('SNP\\n' + ''.join(value + '\\n' for value in selected))\n"
        " remaining = [row for row in rows if row[1] not in set(selected)]\n"
        " header = 'Chr SNP bp freq refA b se p n freq_geno bC bC_se pC\\n'\n"
        " body = ''.join(f'{row[0]} {row[1]} {row[3]} 0.2 G 0.1 0.05 1e-6 1000 0.2 0.08 0.04 0.04\\n' for row in remaining)\n"
        " Path(prefix + '.cma.cojo').write_text(header + body)\n"
        " raise SystemExit(0)\n"
        "if chromosome in no_signal:\n"
        " Path(prefix + '.log').write_text('No SNPs have been selected.\\n')\n"
        " raise SystemExit(0)\n"
        "modeled = selected_rows[0]\n"
        "variant = modeled[1]\n"
        "header = 'Chr SNP bp freq refA b se p n freq_geno bJ bJ_se pJ LD_r\\n'\n"
        "body = f'{modeled[0]} {variant} {modeled[3]} 0.2 G 0.1 0.05 1e-9 1000 0.2 0.09 0.04 0.02 0\\n'\n"
        "Path(prefix + '.jma.cojo').write_text(header + body)\n"
        "Path(prefix + '.ldr.cojo').write_text(f'SNP {variant}\\n{variant} 1\\n')\n"
        "remaining = selected_rows[1:]\n"
        "c_header = 'Chr SNP bp freq refA b se p n freq_geno bC bC_se pC\\n'\n"
        "c_body = ''.join(f'{row[0]} {row[1]} {row[3]} 0.2 G 0.1 0.05 1e-6 1000 0.2 0.08 0.04 0.04\\n' for row in remaining)\n"
        "if emit_selection_cma:\n"
        " Path(prefix + '.cma.cojo').write_text(c_header + c_body)\n"
        % (
            fail_once_chromosome,
            always_fail_chromosome,
            tuple(no_signal_chromosomes),
            emit_selection_cma,
            str(counter),
        ),
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _set_parallel_inputs(args, reference, chromosomes=(1, 2)):
    summary_rows = []
    bim_rows = []
    for chromosome in chromosomes:
        for offset, suffix in enumerate(("lead", "other"), start=1):
            variant = "rs%d%02d" % (chromosome, offset)
            position = chromosome * 1000 + offset
            summary_rows.append(
                "%s G A 0.2 0.1 0.05 %s 1000\n"
                % (variant, "1e-9" if suffix == "lead" else "1e-6")
            )
            bim_rows.append(
                "%d\t%s\t0\t%d\tG\tA\n"
                % (chromosome, variant, position)
            )
    Path(args.gcta_cojo_input_file).write_text(
        "SNP A1 A2 freq b se p N\n" + "".join(summary_rows),
        encoding="utf-8",
    )
    Path(str(reference) + ".bim").write_text(
        "".join(bim_rows), encoding="utf-8",
    )


def _config(tmp_path, reference, *, mode="slct", condition=None, joint=None):
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
    if joint is not None:
        lines.extend(["inputs:", "  joint_snps: '%s'" % joint])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _args(
    tmp_path,
    *,
    mode="slct",
    condition=None,
    joint=None,
    version="1.94.1",
):
    reference = _reference(tmp_path)
    return Namespace(
        gcta_cojo_input_file=str(_summary(tmp_path)),
        dataset_id="STUDY",
        output_directory=str(tmp_path / "out"),
        run_config=str(_config(
            tmp_path, reference, mode=mode, condition=condition, joint=joint,
        )),
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
    assert "Default: slct" in compact_help
    assert "Known SNPs to adjust for" in compact_help
    assert "any explicit mode must be cond" in compact_help
    assert "effects are estimated together in one multiple-SNP model" in compact_help
    assert "any explicit mode must be joint" in compact_help
    assert "Only test or select SNPs listed in this file" in compact_help
    assert "cannot be combined with joint mode" in compact_help
    assert "Do not test or select SNPs listed in this file" in compact_help
    for example in (
        "Select independent signals (default when no SNP list is supplied):",
        "Select up to a requested number of independent signals:",
        "Estimate a specified SNP set jointly:",
        "Condition on known lead variants:",
    ):
        assert example in help_text
    assert "adjusted P-value" not in help_text
    assert "adjusted-p" not in help_text
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
        "show_screen", "threads",
    }
    observed = {action.dest for action in parser._actions if action.dest != "help"}
    assert observed == configurable_destinations
    for action in parser._actions:
        if action.dest in configurable_destinations:
            assert action.default == argparse.SUPPRESS
    normalized_help = " ".join(help_text.split())
    assert "both SNP ID and complete allele pair must match" in normalized_help
    assert "Default: 0.7" in normalized_help


def test_cojo_genomic_control_is_one_way_opt_in_from_yaml_default():
    parser = build_parser()
    help_text = parser.format_help()

    assert "--cojo-gc" in help_text
    assert "--no-cojo-gc" not in help_text
    assert "omit this flag to leave genomic control disabled" in " ".join(
        help_text.split()
    )
    assert load_configuration().modules.gcta_cojo.analysis.genomic_control is False
    assert not hasattr(parser.parse_args([]), "cojo_gc")
    assert parser.parse_args(["--cojo-gc"]).cojo_gc is True


def test_condition_list_implies_cond_mode(tmp_path):
    condition = tmp_path / "condition.txt"
    condition.write_text("rs1\n", encoding="utf-8")
    configuration = resolve_gcta_cojo_configuration(Namespace(
        condition_snps=str(condition),
    ))
    assert configuration.modules.gcta_cojo.mode == "cond"
    assert configuration.modules.gcta_cojo.inputs.condition_snps == str(condition)


def test_pipeline_preflight_validates_external_resources_before_generated_input(
    tmp_path,
):
    args = _args(tmp_path)
    del args.gcta_cojo_input_file

    evidence = preflight_gcta_cojo_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )

    assert evidence.module == "gcta_cojo"
    assert isinstance(evidence.resources, GctaCojoPipelineResources)
    assert evidence.resources.reference.samples == 4
    assert evidence.resources.reference.files["bim"].is_file()
    assert evidence.resources.identifier_observation["variant_id_type"] == "rsid"
    assert evidence.resources.snp_list_metrics == {}
    assert len(evidence.resources.file_identities) == 4
    assert "formatter-created GCTA .ma" in " ".join(evidence.deferred_checks)
    assert not (tmp_path / "out").exists()


def test_pipeline_preflight_rejects_malformed_external_snp_list_before_output(
    tmp_path,
):
    condition = tmp_path / "condition.txt"
    condition.write_text("rs1\nrs1\n", encoding="utf-8")
    args = _args(tmp_path, mode="cond", condition=condition)

    with pytest.raises(GctaCojoError, match="duplicate identifiers"):
        preflight_gcta_cojo_pipeline(
            args,
            preflight_evidence=pipeline_input_vcf_evidence(),
        )

    assert not (tmp_path / "out").exists()


def test_pipeline_execution_reuses_reference_validation_and_rejects_mutation(
    tmp_path,
    monkeypatch,
):
    args = _args(tmp_path)
    evidence = preflight_gcta_cojo_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )
    resources = evidence.resources
    changed_root = tmp_path / "changed"
    changed_root.mkdir()
    changed_args = _args(changed_root)
    changed_evidence = preflight_gcta_cojo_pipeline(
        changed_args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )
    monkeypatch.setattr(
        "postgwas.modules.gcta_cojo.service.validate_gcta_cojo_reference",
        lambda *_args, **_kwargs: pytest.fail(
            "GCTA-COJO reference validation was repeated"
        ),
    )

    result = run_gcta_cojo_direct(args, pipeline_resources=resources)

    assert result.module == "gcta_cojo"
    assert result.metrics["gcta_version"] == "1.94.1"

    with Path(
        str(changed_evidence.resources.reference.prefix) + ".bim"
    ).open("a", encoding="utf-8") as handle:
        handle.write("1\trs2\t0\t200\tT\tC\n")

    with pytest.raises(GctaCojoError, match="changed after pipeline preflight"):
        run_gcta_cojo_direct(
            changed_args,
            pipeline_resources=changed_evidence.resources,
        )


def test_joint_mode_rejects_cojo_extract_during_cli_configuration_resolution(
    tmp_path,
):
    joint = tmp_path / "joint.txt"
    extract = tmp_path / "extract.txt"
    joint.write_text("rs1\n", encoding="utf-8")
    extract.write_text("rs1\n", encoding="utf-8")
    args = build_parser().parse_args([
        "--cojo-mode", "joint",
        "--joint-snps", str(joint),
        "--cojo-extract", str(extract),
    ])

    with pytest.raises(
        ConfigurationError,
        match="joint mode uses inputs.joint_snps as GCTA --extract",
    ):
        resolve_gcta_cojo_configuration(args)


def test_top_snp_count_implies_top_mode_and_rejects_incompatible_explicit_mode():
    configuration = resolve_gcta_cojo_configuration(
        Namespace(cojo_top_snps=5)
    )
    assert configuration.modules.gcta_cojo.mode == "top_snps"
    assert configuration.modules.gcta_cojo.analysis.top_snp_count == 5
    assert gcta_cojo_analysis_label(
        configuration.modules.gcta_cojo
    ) == "GCTA-COJO top-SNP selection (5 SNPs requested)"
    assert gcta_cojo_pipeline_title(
        Namespace(cojo_top_snps=5)
    ) == "Run GCTA-COJO top-SNP selection (5 SNPs requested)."
    title_factory = REGISTRY.get("gcta_cojo").pipeline_title_factory
    assert title_factory is not None
    assert resolve_reference(title_factory)(Namespace(cojo_top_snps=5)) == (
        "Run GCTA-COJO top-SNP selection (5 SNPs requested)."
    )

    with pytest.raises(GctaCojoError, match="valid only.*top_snps"):
        resolve_gcta_cojo_configuration(Namespace(
            cojo_mode="slct",
            cojo_top_snps=5,
        ))


@pytest.mark.parametrize(
    ("key", "value"),
    (
        ("modules.gcta_cojo.analysis.significance_threshold", 0.06),
        ("modules.gcta_cojo.analysis.top_snp_count", 10001),
        ("modules.gcta_cojo.analysis.window_kb", 100001),
        ("modules.gcta_cojo.analysis.collinearity_cutoff", 1.0),
        ("modules.gcta_cojo.analysis.genomic_control_lambda", 10.1),
        ("modules.gcta_cojo.analysis.chromosome", 101),
        ("modules.gcta_cojo.input_validation.minimum_summary_sample_size", 9),
        (
            "modules.gcta_cojo.input_validation.maximum_allele_mismatch_examples",
            101,
        ),
    ),
)
def test_upstream_gcta_parameter_limits_are_schema_validated(key, value):
    with pytest.raises(ConfigurationError):
        load_configuration(cli_overrides={key: value})


def test_zero_frequency_difference_threshold_matches_upstream_gcta_range():
    module = load_configuration(cli_overrides={
        "modules.gcta_cojo.analysis.frequency_difference_max": 0,
    }).modules.gcta_cojo
    assert module.analysis.frequency_difference_max == 0


def test_chromosome_parallelism_is_default_only_for_genome_wide_slct():
    defaults = load_configuration()
    module = defaults.modules.gcta_cojo

    assert (
        defaults.modules.ld_clumping.cojo.parallel_output_contract
        == "selection_only"
    )
    assert module.chromosome_execution.enabled is True
    assert module.chromosome_execution.max_workers == "auto"
    assert module.chromosome_execution.failed_chromosome_retries == 1
    assert should_parallelize_slct(module, {"1": 10, "2": 8}) is True

    top_module = load_configuration(cli_overrides={
        "modules.gcta_cojo.mode": "top_snps",
    }).modules.gcta_cojo
    assert should_parallelize_slct(top_module, {"1": 10, "2": 8}) is False

    chromosome_module = load_configuration(cli_overrides={
        "modules.gcta_cojo.analysis.chromosome": 1,
    }).modules.gcta_cojo
    assert should_parallelize_slct(
        chromosome_module, {"1": 10, "2": 8},
    ) is False
    assert numeric_chromosome_workloads({"1": 10, "X": 8}) is None


def test_ld_clumping_parallel_output_contract_is_schema_validated():
    with pytest.raises(ConfigurationError):
        load_configuration(cli_overrides={
            "modules.ld_clumping.cojo.parallel_output_contract": "invalid",
        })


def test_parallel_worker_plan_is_bounded_by_cpu_memory_and_user_cap():
    plan = plan_chromosome_workers(
        8,
        total_threads=8,
        threads_per_worker=2,
        memory_budget_gb=10,
        parent_memory_gb=1,
        estimated_memory_gb_per_worker=3,
        configured_max_workers="auto",
    )

    assert plan.cpu_limit == 4
    assert plan.memory_limit == 3
    assert plan.workers == 3

    capped = plan_chromosome_workers(
        8,
        total_threads=8,
        threads_per_worker=2,
        memory_budget_gb=10,
        parent_memory_gb=1,
        estimated_memory_gb_per_worker=3,
        configured_max_workers=2,
    )
    assert capped.workers == 2

    with pytest.raises(GctaCojoError, match="cannot safely start another"):
        plan_chromosome_workers(
            2,
            total_threads=4,
            threads_per_worker=1,
            memory_budget_gb=2,
            parent_memory_gb=1,
            estimated_memory_gb_per_worker=2,
            configured_max_workers="auto",
        )


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
    if mode == "joint":
        assert command.count("--extract") == 1
        assert command[command.index("--extract") + 1] == str(snps.resolve())


def test_command_builder_prefers_validated_effective_exclusion(tmp_path):
    configured = tmp_path / "configured.txt"
    effective = tmp_path / "effective.txt"
    configured.write_text("rs1\n", encoding="utf-8")
    effective.write_text("rs1\nrs2\n", encoding="utf-8")
    module = load_configuration(cli_overrides={
        "modules.gcta_cojo.inputs.exclude_snps": str(configured),
    }).modules.gcta_cojo

    command = build_cojo_command(
        "gcta64",
        "study.ma",
        "reference",
        "output",
        4,
        module,
        exclude_snps=effective,
    )

    assert command.count("--exclude") == 1
    assert command[command.index("--exclude") + 1] == str(effective.resolve())


def test_command_builder_rejects_unvalidated_joint_extract_combination(tmp_path):
    joint = tmp_path / "joint.txt"
    extract = tmp_path / "extract.txt"
    joint.write_text("rs1\n", encoding="utf-8")
    extract.write_text("rs1\n", encoding="utf-8")
    module = load_configuration(cli_overrides={
        "modules.gcta_cojo.mode": "joint",
        "modules.gcta_cojo.inputs.joint_snps": str(joint),
    }).modules.gcta_cojo
    invalid_inputs = module.inputs.model_copy(update={
        "extract_snps": str(extract),
    })
    invalid_module = module.model_copy(update={"inputs": invalid_inputs})

    with pytest.raises(
        ValueError,
        match="joint mode uses inputs.joint_snps as GCTA --extract",
    ):
        build_cojo_command(
            "gcta64", "study.ma", "reference", "output", 4, invalid_module,
        )


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


def test_version_probe_accepts_optional_logger(monkeypatch):
    monkeypatch.setattr(
        "postgwas.core.gcta.subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(
            command, 1, stdout="GCTA v1.95.3\n", stderr="invalid option\n",
        ),
    )
    module = load_configuration().modules.gcta_cojo

    version = require_supported_gcta("gcta64", module, None, 10)

    assert version == "1.95.3"


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


def test_conditional_header_only_output_is_a_valid_zero_test_set(tmp_path):
    condition = tmp_path / "condition.txt"
    condition.write_text("rs1\n", encoding="utf-8")
    module = load_configuration(cli_overrides={
        "modules.gcta_cojo.mode": "cond",
        "modules.gcta_cojo.inputs.condition_snps": str(condition),
    }).modules.gcta_cojo
    source = tmp_path / "result.cma.cojo"
    source.write_text(
        "Chr SNP bp refA freq b se p n freq_geno bC bC_se pC\n",
        encoding="utf-8",
    )

    metrics, frame = normalize_cojo_results(
        source, tmp_path / "normalized.tsv", module,
    )

    assert frame.is_empty()
    assert metrics["result_rows"] == 0
    assert "no_signals" not in metrics


def test_genomic_control_native_p_value_headers_are_preserved(tmp_path):
    module = load_configuration(cli_overrides={
        "modules.gcta_cojo.analysis.genomic_control": True,
    }).modules.gcta_cojo
    source = tmp_path / "result.jma.cojo"
    source.write_text(
        "Chr SNP bp freq refA b se p_GC n freq_geno bJ bJ_se pJ_GC LD_r\n"
        "1 rs1 100 0.2 G 0.1 0.05 1e-9 1000 0.2 0.08 0.04 0.02 0\n",
        encoding="utf-8",
    )

    metrics, frame = normalize_cojo_results(
        source, tmp_path / "normalized.tsv", module,
    )

    assert metrics["cojo_p_value_column"] == "pJ_GC"
    assert "p_GC" in frame.columns
    assert "pJ_GC" in frame.columns


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
    assert result.artifacts["ld_matrix"].path.name.endswith(".ldr.cojo")
    assert result.artifacts["conditional_results"].path.name.endswith(
        ".cma.cojo"
    )
    assert result.artifacts["completion_manifest"].path.is_file()
    assert result.metrics["normalized_result"] == str(
        result.artifacts["normalized_results"].path
    )
    assert result.metrics["postgwas_version"]
    assert result.metrics["python_version"]
    assert (
        result.artifacts["normalized_results"].metadata["effect_allele_origin"]
        == "input_ma"
    )
    assert "GCTA-COJO completed" in terminal
    assert "Independent/joint signals" in terminal
    assert "BIM identifier type" in terminal
    log = tmp_path / "out" / "logs" / "STUDY_slct_gcta_cojo.log"
    assert "gcta_cojo status=COMPLETED" in log.read_text(encoding="utf-8")


def test_small_reference_warning_states_reason_and_consequence(tmp_path):
    args = _args(tmp_path)
    configuration_path = Path(args.run_config)
    configuration_path.write_text(
        configuration_path.read_text(encoding="utf-8").replace(
            "reference_sample_size_warning_threshold: 1",
            "reference_sample_size_warning_threshold: 4000",
        ),
        encoding="utf-8",
    )

    result = run_gcta_cojo_direct(args)

    warning = next(
        warning for warning in result.warnings
        if "LD reference sample size" in warning
    )
    assert "4 samples available; more than 4,000 recommended" in warning
    assert "signals may be less stable" in warning


def test_genome_wide_slct_runs_by_chromosome_and_restores_complete_outputs(
    tmp_path, capsys,
):
    args = _args(tmp_path)
    args.threads = 4
    args.gcta = str(_parallel_fake_gcta(tmp_path))
    _set_parallel_inputs(args, tmp_path / "reference")

    result = run_gcta_cojo_direct(args)
    terminal = capsys.readouterr().out

    assert result.metrics["execution_strategy"] == "chromosome_parallel"
    assert result.metrics["parallel_output_contract"] == "complete"
    assert result.metrics["conditional_reconstruction_status"] == "completed"
    assert result.metrics["chromosomes"] == [1, 2]
    assert result.metrics["pilot_chromosome"] == 1
    assert result.metrics["chromosome_attempts"] == {"1": 1, "2": 1}
    assert result.metrics["retried_chromosomes"] == []
    assert result.metrics["threads_per_worker"] == 1
    assert result.metrics["result_rows"] == 2
    assert "GCTA-COJO chromosome progress" in terminal
    assert "GCTA-COJO chromosome execution" in terminal
    assert "Chromosomes scheduled" in terminal
    assert "Chromosomes running simultaneously" in terminal
    assert "up to 1" in terminal
    assert "Parallel chromosome analyses" in terminal
    assert "up to 1 simultaneously · 2 total" in terminal
    assert "GCTA-COJO memory pilot" not in terminal
    assert "Memory pilot" not in terminal
    assert "Running now" not in terminal
    assert "Remaining workload" not in terminal
    assert "Concurrency decision" not in terminal
    assert "Concurrency limits" not in terminal
    assert "Completed 2/2 · Select chromosomes" in terminal
    assert "GCTA-COJO final model" in terminal

    normalized = pl.read_csv(
        result.artifacts["normalized_results"].path,
        separator="\t",
    )
    assert normalized["Chr"].to_list() == [1.0, 2.0]
    assert normalized["SNP"].to_list() == ["rs101", "rs201"]

    ld_lines = result.artifacts["ld_matrix"].path.read_text(
        encoding="utf-8",
    ).splitlines()
    assert ld_lines == [
        "SNP\trs101\trs201",
        "rs101\t1\t0",
        "rs201\t0\t1",
    ]
    conditional_text = result.artifacts["conditional_results"].path.read_text(
        encoding="utf-8",
    )
    assert "rs102" in conditional_text
    assert "rs202" in conditional_text
    assert "rs101" not in conditional_text
    assert "rs201" not in conditional_text

    commands = result.metrics["command"]
    selection_commands = [command for command in commands if "--cojo-slct" in command]
    reconstruction = [command for command in commands if "--cojo-cond" in command]
    assert len(selection_commands) == 2
    assert len(reconstruction) == 1
    assert sorted(
        int(command[command.index("--chr") + 1])
        for command in selection_commands
    ) == [1, 2]
    assert "--chr" not in reconstruction[0]
    assert all(
        command[command.index("--cojo-file") + 1]
        == str(Path(args.gcta_cojo_input_file).resolve())
        for command in commands
    )
    combined_log = result.artifacts["gcta_log"].path.read_text(encoding="utf-8")
    assert "chromosome 1 attempt 1: COMPLETED" in combined_log
    assert "chromosome 2 attempt 1: COMPLETED" in combined_log
    assert "genome-wide conditional reconstruction attempt 1" in combined_log
    assert not (tmp_path / "out" / ".partial" / "STUDY_slct").exists()
    canonical_log = (
        tmp_path / "out" / "logs" / "STUDY_slct_gcta_cojo.log"
    ).read_text(encoding="utf-8")
    assert "gcta_cojo_chromosome_workers" in canonical_log
    assert "threads_per_worker=1" in canonical_log


def test_selection_only_parallel_contract_preserves_joint_results_without_cma(
    tmp_path, capsys,
):
    complete_root = tmp_path / "complete"
    selection_root = tmp_path / "selection"
    complete_root.mkdir()
    selection_root.mkdir()

    complete_args = _args(complete_root)
    complete_args.threads = 4
    complete_args.gcta = str(_parallel_fake_gcta(complete_root))
    _set_parallel_inputs(complete_args, complete_root / "reference")
    complete = run_gcta_cojo_direct(complete_args)
    capsys.readouterr()

    selection_args = _args(selection_root)
    selection_args.threads = 4
    selection_args.gcta = str(_parallel_fake_gcta(
        selection_root,
        emit_selection_cma=False,
    ))
    _set_parallel_inputs(selection_args, selection_root / "reference")
    selection = run_gcta_cojo_direct(
        selection_args,
        parallel_slct_output_contract="selection_only",
    )
    terminal = capsys.readouterr().out

    complete_frame = pl.read_csv(
        complete.artifacts["normalized_results"].path,
        separator="\t",
    )
    selection_frame = pl.read_csv(
        selection.artifacts["normalized_results"].path,
        separator="\t",
    )
    assert selection_frame.equals(complete_frame)
    assert selection.metrics["result_rows"] == complete.metrics["result_rows"]
    assert selection.metrics["parallel_output_contract"] == "selection_only"
    assert selection.metrics["conditional_reconstruction_status"] == (
        "not_requested"
    )
    assert selection.metrics["conditional_reconstruction_attempts"] == 0
    assert "conditional_results" not in selection.artifacts
    manifest_text = selection.artifacts["completion_manifest"].path.read_text(
        encoding="utf-8",
    )
    assert ".cma.cojo" not in manifest_text
    assert "parallel_output_contract: selection_only" in manifest_text
    assert "conditional_reconstruction_status: not_requested" in manifest_text
    commands = selection.metrics["command"]
    assert len(commands) == 2
    assert all("--cojo-slct" in command for command in commands)
    assert all("--cojo-cond" not in command for command in commands)
    assert "GCTA-COJO final model" not in terminal
    selection_args.resume = True
    resumed = run_gcta_cojo_direct(
        selection_args,
        parallel_slct_output_contract="selection_only",
    )
    capsys.readouterr()
    assert resumed.metrics["resumed"] is True
    assert resumed.metrics["conditional_reconstruction_status"] == (
        "not_requested"
    )
    assert "conditional_results" not in resumed.artifacts
    canonical_log = (
        selection_root / "out" / "logs" / "STUDY_slct_gcta_cojo.log"
    ).read_text(encoding="utf-8")
    assert "gcta_cojo_conditional_reconstruction" in canonical_log
    assert "reason=selection_only_output_contract" in canonical_log
    assert "scientific_effect=" in canonical_log


def test_selection_only_contract_rejects_untracked_stale_cma(tmp_path):
    args = _args(tmp_path)
    args.threads = 4
    args.gcta = str(_parallel_fake_gcta(
        tmp_path,
        emit_selection_cma=False,
    ))
    _set_parallel_inputs(args, tmp_path / "reference")
    layout = load_configuration().modules.gcta_cojo.output_layout
    stale = Path(args.output_directory) / layout.conditional_results.format(
        dataset_id="STUDY",
        mode="slct",
    )
    stale.parent.mkdir(parents=True, exist_ok=True)
    stale.write_text("stale conditional output\n", encoding="utf-8")

    with pytest.raises(
        GctaCojoError,
        match="not part of the resolved selection-only output contract",
    ):
        run_gcta_cojo_direct(
            args,
            parallel_slct_output_contract="selection_only",
        )

    assert stale.is_file()
    args.overwrite = True
    result = run_gcta_cojo_direct(
        args,
        parallel_slct_output_contract="selection_only",
    )
    assert not stale.exists()
    assert "conditional_results" not in result.artifacts


def test_parallel_progress_is_saved_to_hidden_screen_log(tmp_path, capsys):
    args = _args(tmp_path)
    args.threads = 4
    args.gcta = str(_parallel_fake_gcta(tmp_path))
    _set_parallel_inputs(args, tmp_path / "reference")
    transcript = tmp_path / "screen.log"

    with capsys.disabled():
        with record_screen(ScreenSettings(False, transcript)):
            run_gcta_cojo_direct(args)

    text = transcript.read_text(encoding="utf-8")
    assert "GCTA-COJO chromosome execution" in text
    assert "Chromosomes scheduled" in text
    assert "Chromosomes running simultaneously" in text
    assert "up to 1" in text
    assert "GCTA-COJO memory pilot" not in text
    assert "Pilot memory" not in text
    assert "Running now" not in text
    assert "Remaining workload" not in text
    assert "Concurrency decision" not in text
    assert "Concurrency limits" not in text
    assert "Completed 2/2 · Select chromosomes" in text
    assert "GCTA-COJO final model" in text
    assert "━" not in text
    assert "\x1b" not in text


def test_failed_chromosome_is_retried_without_rerunning_successes(tmp_path):
    args = _args(tmp_path)
    args.gcta = str(_parallel_fake_gcta(
        tmp_path, fail_once_chromosome=2,
    ))
    _set_parallel_inputs(args, tmp_path / "reference")

    result = run_gcta_cojo_direct(args)

    assert result.metrics["chromosome_attempts"] == {"1": 1, "2": 2}
    assert result.metrics["retried_chromosomes"] == [2]
    combined_log = result.artifacts["gcta_log"].path.read_text(encoding="utf-8")
    assert "chromosome 2 attempt 1: FAILED" in combined_log
    assert "chromosome 2 attempt 2: COMPLETED" in combined_log


def test_persistent_chromosome_failure_stops_without_partial_final_outputs(
    tmp_path,
):
    args = _args(tmp_path)
    args.gcta = str(_parallel_fake_gcta(
        tmp_path, always_fail_chromosome=2,
    ))
    _set_parallel_inputs(args, tmp_path / "reference")

    with pytest.raises(
        GctaCojoError,
        match="after automatically retrying only the failed chromosome",
    ):
        run_gcta_cojo_direct(args)

    assert not (tmp_path / "out" / "raw").exists()
    assert not (tmp_path / "out" / "results").exists()


def test_no_signal_chromosome_remains_in_conditional_output(tmp_path):
    args = _args(tmp_path)
    args.gcta = str(_parallel_fake_gcta(
        tmp_path, no_signal_chromosomes=(2,),
    ))
    _set_parallel_inputs(args, tmp_path / "reference")

    result = run_gcta_cojo_direct(args)

    assert result.metrics["no_signal_chromosomes"] == [2]
    assert result.metrics["result_rows"] == 1
    conditional_text = result.artifacts["conditional_results"].path.read_text(
        encoding="utf-8",
    )
    assert "rs201" in conditional_text
    assert "rs202" in conditional_text


def test_all_no_signal_chromosomes_produce_resumable_empty_result(tmp_path):
    args = _args(tmp_path)
    args.gcta = str(_parallel_fake_gcta(
        tmp_path, no_signal_chromosomes=(1, 2),
    ))
    _set_parallel_inputs(args, tmp_path / "reference")

    first = run_gcta_cojo_direct(args)

    assert first.metrics["no_signals"] is True
    assert first.metrics["no_signal_chromosomes"] == [1, 2]
    assert first.metrics["result_rows"] == 0
    assert "raw_results" not in first.artifacts
    assert "conditional_results" not in first.artifacts

    args.resume = True
    resumed = run_gcta_cojo_direct(args)
    assert resumed.metrics["resumed"] is True
    assert resumed.metrics["no_signals"] is True


def test_exclude_list_validation_does_not_replace_reference_overlap_metric(tmp_path):
    exclude = tmp_path / "exclude.txt"
    exclude.write_text("rs_absent\n", encoding="utf-8")
    args = _args(tmp_path)
    args.cojo_exclude = str(exclude)

    result = run_gcta_cojo_direct(args)

    assert result.metrics["reference_overlap_variants"] == 1
    assert result.metrics["reference_overlap_fraction"] == 1.0
    assert result.metrics["effective_exclude_variants"] == 1
    assert any(
        "Configured exclusion IDs with no effect: 1"
        in warning
        for warning in result.warnings
    )


def test_reversed_bim_allele_order_is_gcta_usable(tmp_path):
    args = _args(tmp_path)
    reference = tmp_path / "reference"
    _replace_validation_inputs(
        args,
        reference,
        [("rs1", "G", "A")],
        [("rs1", "A", "G")],
    )

    result = run_gcta_cojo_direct(args)

    assert result.metrics["reference_overlap_variants"] == 1
    assert result.metrics["reference_allele_mismatch_variants"] == 0
    assert "allele_mismatch_exclusion" not in result.artifacts


def test_stale_optional_mismatch_exclusion_requires_overwrite(tmp_path):
    args = _args(tmp_path)
    layout = load_configuration().modules.gcta_cojo.output_layout
    stale = Path(args.output_directory) / layout.allele_mismatch_exclusion.format(
        dataset_id="STUDY",
        mode="slct",
    )
    stale.parent.mkdir(parents=True, exist_ok=True)
    stale.write_text("rs_stale\n", encoding="utf-8")

    with pytest.raises(GctaCojoError, match="Output already exists"):
        run_gcta_cojo_direct(args)

    args.overwrite = True
    result = run_gcta_cojo_direct(args)
    assert not stale.exists()
    assert "allele_mismatch_exclusion" not in result.artifacts


def test_usable_overlap_at_seventy_percent_warns_and_excludes_allele_mismatches(
    tmp_path,
    capsys,
):
    args = _args(tmp_path)
    reference = tmp_path / "reference"
    summary_rows = [
        ("rs%d" % index, "G", "A") for index in range(1, 11)
    ]
    bim_rows = [
        ("rs%d" % index, "G", "A") for index in range(1, 8)
    ] + [("rs8", "T", "C")]
    _replace_validation_inputs(args, reference, summary_rows, bim_rows)
    user_exclude = tmp_path / "user_exclude.txt"
    user_exclude.write_text("rs2\n", encoding="utf-8")
    args.cojo_exclude = str(user_exclude)

    result = run_gcta_cojo_direct(args)
    terminal = capsys.readouterr().out

    assert result.metrics["reference_id_overlap_variants"] == 8
    assert result.metrics["reference_id_overlap_fraction"] == 0.8
    assert result.metrics["reference_overlap_variants"] == 7
    assert result.metrics["reference_overlap_fraction"] == 0.7
    assert result.metrics["reference_missing_variants"] == 2
    assert result.metrics["reference_allele_mismatch_variants"] == 1
    assert result.metrics["reference_allele_mismatch_examples"] == [{
        "variant_id": "rs8",
        "summary_alleles": ["G", "A"],
        "bim_alleles": ["T", "C"],
    }]
    assert result.metrics["effective_exclude_variants"] == 2
    exclusion = result.artifacts["allele_mismatch_exclusion"].path
    assert exclusion.read_text(encoding="utf-8").splitlines() == ["rs2", "rs8"]
    command = result.metrics["command"]
    assert "--exclude" in command
    assert Path(command[command.index("--exclude") + 1]).name.endswith(
        ".allele_mismatch.exclude.txt"
    )
    assert "GCTA-usable GWAS/reference variants" in terminal
    assert "7/10 variants (70.00%)" in terminal
    assert "Shared IDs with incompatible alleles" in terminal
    assert "rs8 (.ma G/A; BIM T/C)" in terminal
    compact_terminal = " ".join(terminal.split())
    assert (
        "GWAS variants absent from the LD reference: 2 of 10"
        in compact_terminal
    )
    assert "GWAS/reference allele mismatch: 1 of 8" in compact_terminal
    log = tmp_path / "out" / "logs" / "STUDY_slct_gcta_cojo.log"
    log_text = log.read_text(encoding="utf-8")
    assert "reference_allele_mismatch_variants=1" in log_text
    assert "rs8" in log_text

    args.resume = True
    resumed = run_gcta_cojo_direct(args)
    assert resumed.metrics["resumed"] is True
    assert resumed.artifacts["allele_mismatch_exclusion"].path == exclusion


def test_usable_overlap_below_seventy_percent_fails_before_gcta(tmp_path, capsys):
    args = _args(tmp_path)
    reference = tmp_path / "reference"
    summary_rows = [
        ("rs%d" % index, "G", "A") for index in range(1, 11)
    ]
    bim_rows = [
        ("rs%d" % index, "G", "A") for index in range(1, 7)
    ] + [("rs7", "T", "C")]
    _replace_validation_inputs(args, reference, summary_rows, bim_rows)

    with pytest.raises(GctaCojoError, match="Only 60.00%.*minimum is 70.00%"):
        run_gcta_cojo_direct(args)

    terminal = capsys.readouterr().out
    terminal_text = " ".join(terminal.split())
    assert "usable: 6" in terminal
    assert "IDs absent from the BIM: 3" in terminal_text
    assert "shared variant IDs with incompatible alleles: 1" in terminal_text
    assert "rs7 (.ma G/A; BIM T/C)" in terminal_text
    assert not (tmp_path / "out" / "raw").exists()


@pytest.mark.parametrize("mode", ("cond", "joint", "extract"))
def test_requested_snp_with_incompatible_alleles_always_fails(tmp_path, mode):
    requested = tmp_path / (mode + ".txt")
    requested.write_text("rs4\n", encoding="utf-8")
    kwargs = {"mode": mode} if mode in {"cond", "joint"} else {}
    if mode == "cond":
        kwargs["condition"] = requested
    elif mode == "joint":
        kwargs["joint"] = requested
    args = _args(tmp_path, **kwargs)
    if mode == "extract":
        args.cojo_extract = str(requested)
    reference = tmp_path / "reference"
    summary_rows = [
        ("rs%d" % index, "G", "A") for index in range(1, 5)
    ]
    bim_rows = [
        ("rs%d" % index, "G", "A") for index in range(1, 4)
    ] + [("rs4", "T", "C")]
    _replace_validation_inputs(args, reference, summary_rows, bim_rows)

    with pytest.raises(
        GctaCojoError,
        match="1 have incompatible .ma/BIM allele pairs",
    ):
        run_gcta_cojo_direct(args)


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
    assert result.artifacts["condition_snps"].path.name.endswith(".given.cojo")
    assert result.artifacts["raw_results"].path.name.endswith(".cma.cojo")
    assert "Conditional findings" in terminal


def test_genomic_control_run_accepts_gcta_gc_headers(tmp_path):
    args = _args(tmp_path)
    args.cojo_gc = True

    result = run_gcta_cojo_direct(args)

    assert result.metrics["cojo_p_value_column"] == "pJ_GC"
    normalized = pl.read_csv(
        result.artifacts["normalized_results"].path,
        separator="\t",
    )
    assert "p_GC" in normalized.columns
    assert "pJ_GC" in normalized.columns


def test_no_selected_signals_is_successful_and_resumable(tmp_path, capsys):
    args = _args(tmp_path)
    args.gcta = str(_fake_gcta(tmp_path, no_signals=True))

    first = run_gcta_cojo_direct(args)
    terminal = capsys.readouterr().out

    assert first.metrics["no_signals"] is True
    assert first.metrics["result_rows"] == 0
    assert "raw_results" not in first.artifacts
    assert "ld_matrix" not in first.artifacts
    assert "selected no independent signals" in terminal
    assert first.artifacts["normalized_results"].path.is_file()

    args.resume = True
    resumed = run_gcta_cojo_direct(args)
    assert resumed.metrics["resumed"] is True
    assert resumed.metrics["no_signals"] is True


@pytest.mark.parametrize("mode", ("cond", "joint"))
def test_requested_modeled_snps_must_survive_gcta_exactly(tmp_path, mode):
    snps = tmp_path / (mode + "_snps.txt")
    snps.write_text("rs1\n", encoding="utf-8")
    kwargs = {"mode": mode, "condition" if mode == "cond" else "joint": snps}
    args = _args(tmp_path, **kwargs)
    args.gcta = str(_fake_gcta(tmp_path, replace_modeled_snp=True))

    with pytest.raises(GctaCojoError, match="does not exactly match"):
        run_gcta_cojo_direct(args)


def test_completed_run_resumes_only_from_validated_manifest(tmp_path):
    args = _args(tmp_path)
    first = run_gcta_cojo_direct(args)
    args.resume = True

    resumed = run_gcta_cojo_direct(args)

    assert first.artifacts["normalized_results"].path == (
        resumed.artifacts["normalized_results"].path
    )
    assert resumed.metrics["resumed"] is True


def test_completed_run_restarts_when_every_declared_output_is_missing(tmp_path):
    args = _args(tmp_path)
    first = run_gcta_cojo_direct(args)
    completion = first.artifacts["completion_manifest"].path
    for name, artifact in first.artifacts.items():
        if name != "completion_manifest":
            artifact.path.unlink()
    args.resume = True

    restarted = run_gcta_cojo_direct(args)

    assert restarted.metrics["resumed"] is False
    assert restarted.artifacts["normalized_results"].path.is_file()
    assert completion.is_file()
    log = tmp_path / "out" / "logs" / "STUDY_slct_gcta_cojo.log"
    assert "reason=incomplete_outputs" in log.read_text(
        encoding="utf-8"
    )


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_pipeline_vcf_mode_formats_then_runs_cojo_with_bim_identifier_type(
    tmp_path, capsys, monkeypatch,
):
    prefix = tmp_path / "unique_reference"
    Path(str(prefix) + ".bed").write_bytes(bytes((0x6C, 0x1B, 0x01)) + bytes(4))
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
        cojo_mode="top_snps",
        cojo_top_snps=10,
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
    preflight = preflight_gcta_cojo_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )
    context = RunContext(validations={"gcta_cojo": preflight})
    monkeypatch.setattr(
        "postgwas.pipeline.runners._validate_current_pipeline_vcf",
        lambda *_args: None,
    )

    run_formatter_runner(args, context)
    formatter_screen = capsys.readouterr().out
    args._step_num = "02"
    result = run_gcta_cojo_runner(args, context)

    analysis_label = "GCTA-COJO top-SNP selection (10 SNPs requested)"
    assert analysis_label in formatter_screen
    assert "GCTA COJO / fastBAT / mBAT-combo" not in formatter_screen
    assert args.variant_id_observations["gcta_gene"]["consumers"] == [
        analysis_label
    ]
    formatter_log = output / "01_formatter" / "logs" / "STUDY_formatter.log"
    formatter_log_text = formatter_log.read_text(encoding="utf-8")
    assert "downstream_consumer='GCTA-COJO top-SNP selection" in formatter_log_text
    assert "(10 SNPs requested)'" in formatter_log_text
    formatter_report = (
        output / "01_formatter" / "reports" / "STUDY_formatter_report.html"
    )
    assert analysis_label in formatter_report.read_text(encoding="utf-8")
    assert result.metrics["mode"] == "top_snps"
    assert result.metrics["variant_id_type"] == "unique"
    assert result.metrics["summary_variants"] == 3
    formatted = output / "01_formatter" / "STUDY_gcta.ma"
    identifiers = pl.read_csv(formatted, separator="\t")["SNP"].to_list()
    assert identifiers == ["1_100_A_G", "1_200_C_T", "2_300_G_A"]
    assert args.variant_id_types["gcta_gene"] == "unique"
    assert (
        result.artifacts["normalized_results"].metadata["effect_allele_origin"]
        == "VCF_ALT"
    )
    assert result.artifacts["normalized_results"].path.parent == (
        output / "02_gcta_cojo" / "results"
    )


def test_direct_mode_requires_existing_cojo_file(tmp_path, capsys):
    reference = _reference(tmp_path)
    args = Namespace(
        dataset_id="STUDY",
        output_directory=str(tmp_path / "out"),
        cojo_reference_prefix=str(reference),
        cojo_mode="top_snps",
        cojo_top_snps=10,
        cojo_reference_population="EUR",
        genome_build="GRCh37",
        gcta=str(_fake_gcta(tmp_path)),
        threads=2,
        resume=False,
        overwrite=False,
        dry_run=False,
    )

    with pytest.raises(ConfigurationError, match="--cojo-file"):
        run_gcta_cojo_direct(args)

    terminal = capsys.readouterr().out
    assert "GCTA-COJO failed" in terminal
    assert "Required argument not provided: --cojo-file" in terminal
    assert "modules.gcta_cojo.input_file" in terminal


def test_all_missing_base_requirements_are_reported_together(tmp_path):
    args = Namespace(
        dataset_id="STUDY",
        output_directory=str(tmp_path / "out"),
        resume=False,
        overwrite=False,
        dry_run=False,
    )

    with pytest.raises(ConfigurationError) as captured:
        run_gcta_cojo_direct(args)

    message = str(captured.value)
    for option in (
        "--cojo-file",
        "--cojo-reference-prefix",
        "--genome-build",
        "--cojo-reference-population",
    ):
        assert "Required argument not provided: %s" % option in message


@pytest.mark.parametrize(
    ("mode", "option", "configuration_path"),
    (
        ("cond", "--condition-snps", "modules.gcta_cojo.inputs.condition_snps"),
        ("joint", "--joint-snps", "modules.gcta_cojo.inputs.joint_snps"),
    ),
)
def test_mode_specific_missing_list_uses_shared_required_error(
    tmp_path,
    mode,
    option,
    configuration_path,
):
    args = _args(tmp_path, mode=mode)

    with pytest.raises(ConfigurationError) as captured:
        run_gcta_cojo_direct(args)

    assert str(captured.value) == (
        "Required argument not provided: %s. Provide %s VALUE or set %s in the "
        "run configuration." % (option, option, configuration_path)
    )


def test_official_gcta_ma_header_is_accepted_and_values_are_validated(tmp_path):
    args = _args(tmp_path)
    summary = Path(args.gcta_cojo_input_file)
    summary.write_text(
        "SNP A1 A2 freq b se p N\n"
        "rs1 G A 0.2 0.1 0.05 1e-9 1000.5\n",
        encoding="utf-8",
    )
    result = run_gcta_cojo_direct(args)
    assert result.metrics["summary_variants"] == 1

    invalid_directory = tmp_path / "invalid"
    invalid_directory.mkdir()
    invalid_args = _args(invalid_directory)
    invalid_summary = Path(invalid_args.gcta_cojo_input_file)
    invalid_summary.write_text(
        "SNP A1 A2 freq b se p N\n"
        "rs1 G A 0.2 0.1 0.05 1.2 1000\n",
        encoding="utf-8",
    )
    with pytest.raises(GctaCojoError, match="P-value must be within"):
        run_gcta_cojo_direct(invalid_args)

    low_n_directory = tmp_path / "low_n"
    low_n_directory.mkdir()
    low_n_args = _args(low_n_directory)
    Path(low_n_args.gcta_cojo_input_file).write_text(
        "SNP A1 A2 freq b se p N\n"
        "rs1 G A 0.2 0.1 0.05 1e-9 9.5\n",
        encoding="utf-8",
    )
    with pytest.raises(GctaCojoError, match="sample size must be at least"):
        run_gcta_cojo_direct(low_n_args)


@pytest.mark.parametrize("list_option", ("cojo_exclude", "cojo_extract"))
def test_conditioning_snps_cannot_be_excluded_or_omitted_from_extract(
    tmp_path,
    list_option,
):
    condition = tmp_path / "condition.txt"
    condition.write_text("rs1\n", encoding="utf-8")
    restriction = tmp_path / (list_option + ".txt")
    args = _args(tmp_path, mode="cond", condition=condition)
    if list_option == "cojo_extract":
        with Path(args.gcta_cojo_input_file).open("a", encoding="utf-8") as handle:
            handle.write("rs2\tT\tC\t0.3\t0.2\t0.06\t2e-6\t1000\n")
        with Path(str(tmp_path / "reference") + ".bim").open(
            "a", encoding="utf-8",
        ) as handle:
            handle.write("1\trs2\t0\t200\tT\tC\n")
    restriction.write_text(
        "rs1\n" if list_option == "cojo_exclude" else "rs2\n",
        encoding="utf-8",
    )
    setattr(args, list_option, str(restriction))

    with pytest.raises(GctaCojoError, match="conditioning SNPs"):
        run_gcta_cojo_direct(args)


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
    preflight = preflight_gcta_cojo_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )
    context = RunContext({
        "formatter": {
            "gcta_gene": {"summary_statistics_input_file": str(summary)},
        },
    }, validations={"gcta_cojo": preflight})

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
        cojo_mode="top_snps",
        cojo_top_snps=10,
        cojo_reference_population="EUR",
        genome_build="GRCh37",
        output_directory=str(output),
        dataset_id="STUDY",
        gcta=str(_fake_gcta(tmp_path)),
        bcftools=shutil.which("bcftools") or shutil.which("python"),
        _step_num="01",
    )
    preflight = preflight_gcta_cojo_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )
    context = RunContext(validations={"gcta_cojo": preflight})
    monkeypatch.setattr(
        "postgwas.pipeline.runners._validate_current_pipeline_vcf",
        lambda *_args: None,
    )

    result = run_formatter_runner(args, context)

    assert result["gcta_gene"]["summary_statistics_input_file"] == "study.ma"
    assert captured["target_types"] == {"gcta_gene": "unique"}
    assert captured["observed"]["unique_ids"] == 1
    assert captured["observed"]["consumers"] == [
        "GCTA-COJO top-SNP selection (10 SNPs requested)"
    ]
    assert captured["formats"] == ["gcta_gene"]
    assert args.output_directory == str(output)
