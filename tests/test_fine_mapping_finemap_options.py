"""Regression tests for canonical FINEMAP configuration and command controls."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from unittest.mock import patch

import pytest
import polars as pl

from postgwas.config import load_configuration
from postgwas.cli.common import get_plink_binary_parser
from postgwas.core.errors import ConfigurationError
from postgwas.core.processes import (
    SupervisedProcessTimeout,
    run_supervised_process,
)
from postgwas.modules.fine_mapping.arguments import (
    get_common_susie_arguments,
    get_common_finemap_finemap_arguments,
    get_finemap_common_parser,
)
from postgwas.modules.fine_mapping.engines.finemap.adapter import (
    _effective_max_causal_snps,
    _prepare_loci,
    _reconcile_model_selection_failures,
    generate_tasks,
    process_single_locus,
)
from postgwas.modules.fine_mapping.engines.finemap.merge_results import (
    INVALID_MODEL_PROBABILITY_REASON,
    InvalidModelProbabilityError,
    create_flames_file,
    parse_cred_header,
    process_finemap_output,
    save_summary_csv,
    select_and_copy_best_models,
)
from postgwas.modules.fine_mapping.engines.finemap.runner import (
    FinemapExternalToolTimeout,
    run_finemap_binary,
    run_ldstore,
)
from postgwas.modules.fine_mapping.engines.susie.main import (
    resolve_susie_r_runtime,
    run_susie,
)
from postgwas.modules.fine_mapping.service import (
    _resolve_fine_mapping_arguments,
)
from postgwas.modules.fine_mapping.resource_guard import dense_ld_resource_qc
from postgwas.modules.fine_mapping.progress import write_run_configuration


def _finemap_config(**overrides):
    settings = (
        load_configuration()
        .modules.fine_mapping.engines.finemap.model_dump()
    )
    settings.update(overrides)
    settings["prob_cred_set"] = 0.95
    return settings


def _resource_config():
    module = load_configuration().modules.fine_mapping
    return {
        "maximum_variants_per_locus": (
            module.ld_resource_guard.maximum_variants_per_locus
        ),
        "finemap_ld_peak_matrix_multiplier": (
            module.ld_resource_guard.finemap_peak_matrix_multiplier
        ),
        "minimum_memory_per_worker_gb": module.memory_per_worker_gb,
    }


def test_packaged_finemap_defaults_match_the_documented_baseline():
    module = load_configuration().modules.fine_mapping

    assert module.locus_window_kb == 500
    assert module.credible_set_coverage == 0.95
    assert module.memory_per_worker_gb == 14.0
    assert module.sample_size.model_dump() == {
        "policy": "warn",
        "summary_statistic": "median",
        "relative_range_warning_threshold": 0.05,
    }
    assert module.ld_resource_guard.model_dump() == {
        "maximum_variants_per_locus": 30000,
        "finemap_peak_matrix_multiplier": 4.0,
        "susie_peak_matrix_multiplier": 5.0,
    }
    assert module.engines.susie.recovery_audit_filename == "susie_recovery_audit.tsv"
    assert module.engines.finemap.model_dump() == {
        "algorithm": "sss",
        "ldstore_timeout_seconds": 21600,
        "finemap_timeout_seconds": 43200,
        "termination_grace_seconds": 30,
        "plink_memory_mb": 2000,
        "bgen_bits": 8,
        "external_tool_threads": 1,
        "n_causal_snps": 5,
        "n_iter": 100000,
        "n_conv_sss": 100,
        "prob_conv_sss_tol": 0.001,
        "n_configs_top": 50000,
        "corr_config": 0.95,
        "pvalue_snps": 1.0,
        "cond_pvalue": 5e-8,
        "prior_std": 0.05,
        "prior_k": False,
        "force_n_samples": False,
        "std_effects": False,
        "flames_manifest_filename": "finemap_FLAMES_manifest.tsv",
    }


def test_predefined_range_boundaries_can_be_preserved_or_flanked(tmp_path):
    locus_file = tmp_path / "loci.tsv"
    locus_file.write_text(
        "CHR\tSTART\tEND\tLP\n1\t1000000\t2000000\t8.0\n",
        encoding="utf-8",
    )
    args = argparse.Namespace(
        locus_file=locus_file,
        locus_type="range",
        window_kb=0,
        lp_threshold=7.3,
        finemap_skip_mhc=False,
        finemap_include_mhc=True,
    )

    exact = _prepare_loci(args).row(0, named=True)
    assert (exact["START"], exact["END"]) == (1000000, 2000000)

    args.window_kb = 500
    flanked = _prepare_loci(args).row(0, named=True)
    assert (flanked["START"], flanked["END"]) == (500000, 2500000)


def test_causal_snp_limit_is_a_per_locus_upper_bound():
    assert _effective_max_causal_snps(5, 20) == 5
    assert _effective_max_causal_snps(5, 3) == 3
    with pytest.raises(ValueError, match="at least 1"):
        _effective_max_causal_snps(0, 3)


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        ("# Post-Pr = 0\n", 0.0),
        ("# Post-Pr = 1\n", 1.0),
        ("# Post-Pr=.5\n", 0.5),
        ("# Post-Pr = 1e-3\n", 0.001),
    ],
)
def test_finemap_model_probability_accepts_probability_boundaries(
    tmp_path, header, expected
):
    credible_file = tmp_path / "locus.cred1"
    credible_file.write_text(
        f"{header}cred1 prob1\nrs1 1\n", encoding="utf-8"
    )

    assert parse_cred_header(credible_file) == pytest.approx(expected)


@pytest.mark.parametrize(
    "header",
    [
        "# FINEMAP credible set\n",
        "# Post-Pr = invalid\n",
        "# Post-Pr = NaN\n",
        "# Post-Pr = Inf\n",
        "# Post-Pr = -0.01\n",
        "# Post-Pr = 1.01\n",
        "# Post-Pr = 1e309\n",
        "# Post-Pr = 0.4\n# Post-Pr = 0.6\n",
    ],
)
def test_finemap_model_probability_rejects_invalid_headers(tmp_path, header):
    credible_file = tmp_path / "locus.cred1"
    credible_file.write_text(
        f"{header}cred1 prob1\nrs1 1\n", encoding="utf-8"
    )

    with pytest.raises(InvalidModelProbabilityError):
        parse_cred_header(credible_file)


def test_invalid_model_probability_fails_only_its_locus(tmp_path):
    raw = tmp_path / "raw"
    selected = tmp_path / "selected"
    good = raw / "good"
    bad = raw / "bad"
    good.mkdir(parents=True)
    bad.mkdir()
    (good / "good.z").write_text("z data\n", encoding="utf-8")
    (bad / "bad.z").write_text("z data\n", encoding="utf-8")
    (good / "good.cred1").write_text(
        "# Post-Pr = 0.8\ncred1 prob1\nrs1 1\n", encoding="utf-8"
    )
    (good / "good.cred2").write_text(
        "# Post-Pr = 0.2\ncred1 prob1\nrs1 1\n", encoding="utf-8"
    )
    (bad / "bad.cred1").write_text(
        "# Post-Pr = invalid\ncred1 prob1\nrs1 1\n", encoding="utf-8"
    )
    (bad / "bad.cred2").write_text(
        "# Post-Pr = 0.4\ncred1 prob1\nrs1 1\n", encoding="utf-8"
    )
    failures = []

    tracking = select_and_copy_best_models(
        raw,
        selected,
        allowed_loci={"good", "bad"},
        failure_records=failures,
    )
    reconciled = _reconcile_model_selection_failures(
        [
            {"locus_id": "good", "status": "success", "failure_reason": None},
            {"locus_id": "bad", "status": "success", "failure_reason": None},
        ],
        failures,
    )

    assert (selected / "good.cred1").is_file()
    assert (selected / "good.z").is_file()
    assert not list(selected.glob("bad.cred*"))
    assert len(failures) == 1
    assert failures[0]["locus_id"] == "bad"
    assert (
        failures[0]["failure_reason"]
        == INVALID_MODEL_PROBABILITY_REASON
    )
    assert {row["Locus"] for row in tracking if row["Selected"]} == {"good"}
    assert reconciled[0]["status"] == "success"
    assert reconciled[1]["status"] == "failed"
    assert (
        reconciled[1]["failure_reason"]
        == INVALID_MODEL_PROBABILITY_REASON
    )
    assert "bad.cred1" in reconciled[1]["failure_detail"]
    save_summary_csv(tracking, selected / "finemap_all_models_summary.csv")
    model_summary = (selected / "finemap_all_models_summary.csv").read_text(
        encoding="utf-8"
    )
    assert "Validation_Status" in model_summary
    assert INVALID_MODEL_PROBABILITY_REASON in model_summary
    assert "bad.cred1" in model_summary


def test_all_invalid_finemap_models_return_auditable_failures(tmp_path):
    raw = tmp_path / "raw"
    invalid_locus = raw / "invalid_locus"
    invalid_locus.mkdir(parents=True)
    (invalid_locus / "invalid_locus.cred1").write_text(
        "# Post-Pr = malformed\ncred1 prob1\nrs1 1\n",
        encoding="utf-8",
    )
    intermediate = tmp_path / "intermediate"
    final = tmp_path / "final"
    model_summary = tmp_path / "quality_control" / "finemap_models.csv"

    result = process_finemap_output(
        raw,
        intermediate,
        final,
        model_summary,
        final / "annotations",
        "indexfile.txt",
        "finemap_FLAMES_manifest.tsv",
        "annotated_",
        allowed_loci={"invalid_locus"},
    )

    assert result["selected_loci"] == set()
    assert result["failures"][0]["locus_id"] == "invalid_locus"
    assert (
        result["failures"][0]["failure_reason"]
        == INVALID_MODEL_PROBABILITY_REASON
    )
    assert model_summary.is_file()
    assert not intermediate.exists()
    assert not final.exists()
    assert not (final / "indexfile.txt").exists()


def test_finemap_records_per_locus_nef_warning_without_failing(tmp_path):
    dirs = {
        "root": tmp_path,
        "temp": tmp_path / "temp",
        "loci": tmp_path / "loci",
    }
    dirs["temp"].mkdir()
    dirs["loci"].mkdir()
    loci = pl.DataFrame({
        "CHR": ["1"], "START": [1], "END": [3],
        "GenomicLocus": ["chr1:1-3"],
    })
    sumstats = pl.DataFrame({
        "rsid": ["rs1", "rs2", "rs3"],
        "chromosome": ["1", "1", "1"],
        "position": [1, 2, 3],
        "allele1": ["A", "C", "G"],
        "allele2": ["G", "T", "A"],
        "maf": [0.1, 0.2, 0.3],
        "beta": [0.1, 0.2, 0.3],
        "se": [0.1, 0.1, 0.1],
        "NEF": [60.0, 100.0, 102.0],
    })
    args = argparse.Namespace(
        **_finemap_config(),
        plink="plink2",
        bgenix="bgenix",
        ldstore="ldstore",
        finemap_executable="finemap",
        sample_size_policy="warn",
        sample_size_summary_statistic="median",
        sample_size_relative_range_warning_threshold=0.05,
        **_resource_config(),
    )

    tasks, skipped = generate_tasks(
        loci,
        sumstats,
        "reference",
        dirs,
        args,
        tmp_path / "task_generation_skips.tsv",
    )

    assert skipped == []
    assert tasks[0]["n_samples"] == 100
    assert tasks[0]["sample_size_qc"]["nef_min"] == 60
    assert tasks[0]["sample_size_qc"]["nef_max"] == 102
    assert tasks[0]["sample_size_qc"]["nef_relative_range"] == pytest.approx(0.42)
    assert (
        tasks[0]["sample_size_qc"]["warning_reason"]
        == "nef_relative_range_exceeds_threshold"
    )


def test_finemap_cli_has_no_independent_defaults_or_undocumented_control():
    parsers = (
        get_finemap_common_parser(),
        get_common_susie_arguments(),
        get_common_finemap_finemap_arguments(),
    )

    assert all(vars(parser.parse_args([])) == {} for parser in parsers)
    help_text = parsers[-1].format_help()
    common_help_text = parsers[0].format_help()
    assert "--maximum-variants-per-locus" in common_help_text
    assert "--prior-std" in help_text
    assert "--ldstore-timeout-seconds" in help_text
    assert "--finemap-timeout-seconds" in help_text
    assert "--finemap-termination-grace-seconds" in help_text
    assert "--collinear-tol" not in help_text
    assert vars(get_plink_binary_parser().parse_args([])) == {}


def test_plink_executable_is_resolved_from_yaml_or_explicit_cli(tmp_path):
    config_file = tmp_path / "run.yaml"
    config_file.write_text(
        "resources:\n"
        "  executables:\n"
        "    plink: configured-plink\n"
        "modules:\n"
        "  fine_mapping:\n"
        "    engine: susie\n",
        encoding="utf-8",
    )

    yaml_args = _resolve_fine_mapping_arguments(
        argparse.Namespace(finemap_method="susie", run_config=str(config_file))
    )
    assert yaml_args.plink == "configured-plink"

    cli_args = _resolve_fine_mapping_arguments(
        argparse.Namespace(
            finemap_method="susie",
            run_config=str(config_file),
            plink="explicit-plink",
        )
    )
    assert cli_args.plink == "explicit-plink"


def test_rscript_executable_is_resolved_from_yaml(tmp_path):
    config_file = tmp_path / "run.yaml"
    config_file.write_text(
        "resources:\n"
        "  executables:\n"
        "    rscript: configured-Rscript\n"
        "modules:\n"
        "  fine_mapping:\n"
        "    engine: susie\n",
        encoding="utf-8",
    )

    resolved = _resolve_fine_mapping_arguments(
        argparse.Namespace(finemap_method="susie", run_config=str(config_file))
    )

    assert resolved.rscript == "configured-Rscript"


def test_yaml_is_resolved_then_explicit_cli_values_override_it(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "credible_set_coverage: 0.9\n"
        "memory_per_worker_gb: 12\n"
        "ld_resource_guard:\n"
        "  maximum_variants_per_locus: 25000\n"
        "engines:\n"
        "  finemap:\n"
        "    algorithm: cond\n"
        "    ldstore_timeout_seconds: 1000\n"
        "    finemap_timeout_seconds: 2000\n"
        "    termination_grace_seconds: 20\n"
        "    n_causal_snps: 9\n",
        encoding="utf-8",
    )
    yaml_args = argparse.Namespace(
        finemap_method="finemap",
        run_config=str(config_file),
    )
    resolved_yaml = _resolve_fine_mapping_arguments(yaml_args)
    assert resolved_yaml.algorithm == "cond"
    assert resolved_yaml.cond is True
    assert resolved_yaml.sss is False
    assert resolved_yaml.n_causal_snps == 9
    assert resolved_yaml.ldstore_timeout_seconds == 1000
    assert resolved_yaml.finemap_timeout_seconds == 2000
    assert resolved_yaml.finemap_termination_grace_seconds == 20
    assert resolved_yaml.prob_cred_set == 0.9
    assert resolved_yaml.minimum_memory_per_worker_gb == 12.0
    assert resolved_yaml.maximum_variants_per_locus == 25000

    cli_args = argparse.Namespace(
        finemap_method="finemap",
        run_config=str(config_file),
        sss=True,
        n_causal_snps=7,
        minimum_memory_per_worker_gb=16,
        maximum_variants_per_locus=40000,
        ldstore_timeout_seconds=3000,
        finemap_timeout_seconds=4000,
        finemap_termination_grace_seconds=40,
    )
    resolved_cli = _resolve_fine_mapping_arguments(cli_args)
    assert resolved_cli.algorithm == "sss"
    assert resolved_cli.n_causal_snps == 7
    assert resolved_cli.minimum_memory_per_worker_gb == 16.0
    assert resolved_cli.maximum_variants_per_locus == 40000
    assert resolved_cli.ldstore_timeout_seconds == 3000
    assert resolved_cli.finemap_timeout_seconds == 4000
    assert resolved_cli.finemap_termination_grace_seconds == 40


def test_python_runtime_values_are_derived_from_resolved_yaml(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "runtime:\n"
        "  tool_version_timeout_seconds: 45\n"
        "  fallback_memory_gb: 24\n"
        "validation:\n"
        "  ld_symmetry_tolerance: 0.00002\n"
        "summary_statistics_preparation:\n"
        "  schema_inference_length: 25000\n"
        "engines:\n"
        "  finemap:\n"
        "    plink_memory_mb: 4096\n"
        "    bgen_bits: 12\n"
        "    external_tool_threads: 2\n",
        encoding="utf-8",
    )

    resolved = _resolve_fine_mapping_arguments(
        argparse.Namespace(finemap_method="finemap", run_config=str(config_file))
    )

    runtime = resolved.fine_mapping_runtime_defaults
    assert runtime["tool_version_timeout_seconds"] == 45
    assert runtime["fallback_memory_gb"] == 24
    assert runtime["ld_symmetry_tolerance"] == pytest.approx(0.00002)
    assert runtime["schema_inference_length"] == 25000
    assert runtime["plink_memory_mb"] == 4096
    assert runtime["bgen_bits"] == 12
    assert runtime["external_tool_threads"] == 2


def test_susie_yaml_and_cli_controls_are_resolved_end_to_end(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "locus_type: point\n"
        "locus_window_kb: 2000\n"
        "locus_lp_threshold: 8.0\n"
        "skip_mhc: false\n"
        "mhc_chromosome: 6\n"
        "mhc_start: 26000000\n"
        "mhc_end: 34000000\n"
        "memory_per_worker_gb: 12\n"
        "engines:\n"
        "  susie:\n"
        "    max_causal_components: 8\n"
        "    minimum_purity: 0.8\n"
        "    fitting:\n"
        "      main_max_iter: 120\n"
        "    execution:\n"
        "      ld_timeout_seconds: 200\n"
        "      susie_timeout_seconds: 300\n",
        encoding="utf-8",
    )

    resolved = _resolve_fine_mapping_arguments(
        argparse.Namespace(finemap_method="susie", run_config=str(config_file))
    )
    assert resolved.locus_type == "point"
    assert resolved.window_kb == 2000
    assert resolved.lp_threshold == 8.0
    assert resolved.finemap_skip_mhc is False
    assert resolved.finemap_include_mhc is True
    assert resolved.finemap_mhc_start == 26000000
    assert resolved.finemap_mhc_end == 34000000
    assert resolved.minimum_memory_per_worker_gb == 12.0
    assert resolved.L == 8
    assert resolved.minimum_purity == 0.8
    assert resolved.ld_timeout_seconds == 200
    assert resolved.susie_timeout_seconds == 300
    assert (
        resolved.resolved_fine_mapping_configuration["engines"]["susie"]
        ["fitting"]["main_max_iter"]
        == 120
    )
    assert resolved.recovery_audit_filename == "susie_recovery_audit.tsv"
    assert resolved.sample_size_policy == "warn"
    assert resolved.sample_size_summary_statistic == "median"
    assert resolved.sample_size_relative_range_warning_threshold == 0.05
    assert resolved.maximum_variants_per_locus == 30000

    overridden = _resolve_fine_mapping_arguments(
        argparse.Namespace(
            finemap_method="susie",
            run_config=str(config_file),
            window_kb=500,
            minimum_purity=0.7,
            susie_main_max_iter=150,
            finemap_skip_mhc=True,
        )
    )
    assert overridden.window_kb == 500
    assert overridden.minimum_purity == 0.7
    assert (
        overridden.resolved_fine_mapping_configuration["engines"]["susie"]
        ["fitting"]["main_max_iter"]
        == 150
    )
    assert overridden.finemap_skip_mhc is True


@pytest.mark.parametrize(
    ("namespace", "message"),
    [
        (
            argparse.Namespace(finemap_method="susie", n_causal_snps=7),
            "FINEMAP controls do not apply",
        ),
        (
            argparse.Namespace(finemap_method="finemap", minimum_purity=0.7),
            "SuSiE controls do not apply",
        ),
    ],
)
def test_engine_specific_cli_controls_cannot_be_silently_ignored(
    namespace, message
):
    with pytest.raises(ValueError, match=message):
        _resolve_fine_mapping_arguments(namespace)


def test_susie_adapter_receives_every_resolved_control(tmp_path):
    resolved_configuration = tmp_path / "resolved_susie_configuration.json"
    resolved_configuration.write_text("{}\n", encoding="utf-8")
    completed = subprocess.CompletedProcess([], 0, stdout="", stderr="")
    with patch(
        "postgwas.modules.fine_mapping.engines.susie.main.subprocess.run",
        return_value=completed,
    ) as run:
        run_susie(
            locus_file=tmp_path / "loci.tsv",
            sumstat_manifest=tmp_path / "locus_sumstats_manifest.tsv",
            sample_id="STUDY",
            ld_ref=tmp_path / "reference",
            plink="plink2",
            output_folder=tmp_path / "output",
            resolved_configuration_file=resolved_configuration,
            rscript="configured-Rscript",
            r_environment={"R_LIBS_USER": "validated-library-order"},
        )

    command = run.call_args.args[0]
    assert command[0] == "configured-Rscript"
    assert run.call_args.kwargs["env"] == {
        "R_LIBS_USER": "validated-library-order"
    }
    assert command[command.index("--sumstat_manifest") + 1] == str(
        tmp_path / "locus_sumstats_manifest.tsv"
    )
    assert command[command.index("--resolved_configuration_file") + 1] == str(
        resolved_configuration
    )
    for obsolete_option in (
        "--lp_threshold",
        "--L",
        "--minimum_purity",
        "--timeout_ld_seconds",
        "--timeout_susie_seconds",
        "--maximum_variants_per_locus",
    ):
        assert obsolete_option not in command


def test_susie_r_runtime_prioritizes_configured_r_library(monkeypatch, tmp_path):
    rscript = tmp_path / "Rscript"
    rscript.write_text("executable", encoding="utf-8")
    rscript.chmod(0o755)
    configured_library = tmp_path / "configured-library"
    user_library = tmp_path / "user-library"
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        if len(calls) == 1:
            return subprocess.CompletedProcess(
                command,
                0,
                stdout=(
                    f"{configured_library}\n"
                    f"{user_library}{os.pathsep}{configured_library}\n"
                ),
                stderr="",
            )
        assert kwargs["env"]["R_LIBS_USER"] == os.pathsep.join(
            [str(configured_library), str(user_library)]
        )
        assert kwargs["env"]["PATH"].split(os.pathsep)[0] == str(tmp_path)
        return subprocess.CompletedProcess(
            command, 0, stdout="R version 4.3.3\n", stderr=""
        )

    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.engines.susie.main.subprocess.run",
        fake_run,
    )

    runtime = resolve_susie_r_runtime(rscript, timeout_seconds=20)

    assert runtime["rscript"] == str(rscript)
    assert runtime["version"] == "R version 4.3.3"
    assert runtime["library_paths"] == [
        str(configured_library),
        str(user_library),
    ]


def test_susie_recovery_audit_filename_must_be_a_tsv_basename(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "engines:\n  susie:\n    recovery_audit_filename: ../audit.tsv\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="recovery_audit_filename"):
        _resolve_fine_mapping_arguments(
            argparse.Namespace(finemap_method="susie", run_config=str(config_file))
        )


def test_sample_size_warning_threshold_is_schema_validated(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "sample_size:\n  relative_range_warning_threshold: -0.01\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="relative_range_warning_threshold"):
        _resolve_fine_mapping_arguments(
            argparse.Namespace(finemap_method="susie", run_config=str(config_file))
        )


def test_susie_ld_mismatch_threshold_order_is_schema_validated(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "engines:\n"
        "  susie:\n"
        "    ld_validation:\n"
        "      mismatch_warning_threshold: 0.2\n"
        "      mismatch_failure_threshold: 0.1\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="mismatch_warning_threshold"):
        _resolve_fine_mapping_arguments(
            argparse.Namespace(finemap_method="susie", run_config=str(config_file))
        )


def test_maximum_variants_per_locus_is_schema_validated(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "ld_resource_guard:\n  maximum_variants_per_locus: 0\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="maximum_variants_per_locus"):
        _resolve_fine_mapping_arguments(
            argparse.Namespace(finemap_method="finemap", run_config=str(config_file))
        )


def test_dense_ld_resource_guard_accepts_limit_and_rejects_next_variant():
    at_limit = dense_ld_resource_qc(30000, 30000, 4.0, 14.0)
    over_limit = dense_ld_resource_qc(30001, 30000, 4.0, 14.0)

    assert at_limit["resource_failure_reason"] is None
    assert at_limit["dense_ld_matrix_gb"] == pytest.approx(6.705522537)
    assert at_limit["estimated_peak_ld_memory_gb"] == pytest.approx(26.822090149)
    assert over_limit["resource_failure_reason"] == (
        "maximum_variants_per_locus_exceeded"
    )


def test_finemap_oversized_locus_is_stopped_before_ld_inputs(tmp_path):
    dirs = {
        "root": tmp_path,
        "temp": tmp_path / "temp",
        "loci": tmp_path / "loci",
    }
    dirs["temp"].mkdir()
    dirs["loci"].mkdir()
    loci = pl.DataFrame({
        "CHR": ["1"], "START": [1], "END": [3],
        "GenomicLocus": ["chr1:1-3"],
    })
    sumstats = pl.DataFrame({
        "rsid": ["rs1", "rs2", "rs3"],
        "chromosome": ["1", "1", "1"],
        "position": [1, 2, 3],
        "allele1": ["A", "C", "G"],
        "allele2": ["G", "T", "A"],
        "maf": [0.1, 0.2, 0.3],
        "beta": [0.1, 0.2, 0.3],
        "se": [0.1, 0.1, 0.1],
        "NEF": [100.0, 100.0, 100.0],
    })
    args = argparse.Namespace(
        **_finemap_config(),
        **{
            **_resource_config(),
            "maximum_variants_per_locus": 2,
        },
        plink="plink2",
        sample_size_policy="warn",
        sample_size_summary_statistic="median",
        sample_size_relative_range_warning_threshold=0.05,
    )

    tasks, skipped = generate_tasks(
        loci,
        sumstats,
        "reference",
        dirs,
        args,
        tmp_path / "task_generation_skips.tsv",
    )

    assert tasks == []
    assert skipped[0]["status"] == "resource_limit_exceeded"
    assert skipped[0]["failure_reason"] == "maximum_variants_per_locus_exceeded"
    assert skipped[0]["input_variant_count"] == 3
    assert "--maximum-variants-per-locus" in skipped[0]["failure_detail"]
    assert not list(dirs["loci"].glob("*.z"))
    assert not list(dirs["loci"].glob("*.snp"))


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is unavailable")
def test_susie_recovery_routes_invalid_ld_before_refitting(tmp_path):
    susie_file = (
        Path(__file__).parents[1]
        / "src/postgwas/modules/fine_mapping/engines/susie/susie.r"
    )
    resolved_configuration_file = write_run_configuration(
        tmp_path / "resolved_susie_configuration.json",
        load_configuration().modules.fine_mapping.model_dump(mode="json"),
    )
    script = f"""
    library <- function(package, ...) {{
      package_name <- as.character(substitute(package))
      if (package_name %in% c("susieR", "processx")) return(invisible(TRUE))
      base::library(package_name, character.only = TRUE, ...)
    }}
    source({str(susie_file)!r})
    SUSIE_DEFAULTS <- load_susie_defaults({str(resolved_configuration_file)!r})
    quiet <- function(...) invisible(NULL)
    alignment_input <- data.table::data.table(
      SNP = c("rs1", "rs2"), CHR = c("1", "1"), BP = c(1L, 2L),
      REF = c("A", "A"), ALT = c("G", "T"), EZ = c(2, 1)
    )
    alignment_ld <- diag(2)
    rownames(alignment_ld) <- colnames(alignment_ld) <- c("rs1", "rs2")
    attr(alignment_ld, "variant_map") <- data.frame(
      SNP = c("rs1", "rs2"), CHR = c("1", "1"), BP = c(1L, 2L),
      counted_allele = c("G", "T"), other_allele = c("A", "A")
    )
    aligned_test <- align_sumstats_to_ld(alignment_input, alignment_ld, quiet)
    data.table::set(
      aligned_test$selected, j = "variable_index",
      value = seq_len(nrow(aligned_test$selected))
    )
    stopifnot(
      nrow(aligned_test$selected) == 1L,
      identical(aligned_test$selected$SNP, "rs1"),
      identical(aligned_test$selected$variable_index, 1L)
    )
    non_psd <- matrix(c(1, 2, 2, 1), nrow = 2)
    caught <- tryCatch(
      validate_ld_for_susie(non_psd, c(1, 1), 100, quiet),
      susie_ld_validation_error = function(e) e
    )
    stopifnot(
      inherits(caught, "susie_ld_validation_error"),
      identical(caught$reason, "ld_non_psd"),
      isTRUE(caught$repairable)
    )

    locus <- list(chr = "1", start = 1L, end = 2L, genomic_locus = "chr1:1-2")
    selected <- data.table::data.table(
      SNP = c("rs1", "rs2"), CHR = c("1", "1"), BP = c(1L, 2L),
      REF = c("A", "C"), ALT = c("G", "T"), EZ = c(1, 1),
      LP = c(8, 8), variable_index = 1:2
    )
    repair_calls <- 0L
    fit_calls <- 0L
    fit_input_paths <- character()
    fit_input_states <- character()
    LD_fix_fast <- function(...) {{ repair_calls <<- repair_calls + 1L; diag(2) }}
    run_susie_hard_timeout <- function(fit_input, ...) {{
      validate_susie_fit_input(fit_input)
      fit_input$child_processes <- fit_input$child_processes + 1L
      fit_calls <<- fit_calls + 1L
      fit_input_paths <<- c(fit_input_paths, fit_input$path)
      fit_input_states <<- c(fit_input_states, fit_input$matrix_state)
      susie_fit_failure("susie_nonconvergence", "test non-convergence")
    }}
    validate_ld_for_susie <- function(...) {{
      list(ld = diag(2), lambda = 0, min_eigenvalue = 1, warning_reason = NA_character_)
    }}
    defaults <- utils::modifyList(SUSIE_DEFAULTS, list(
      recovery_timeout_multiplier = 1,
      repaired_timeout_multiplier = 1,
      ld_repair_timeout_multiplier = 1,
      sample_size_policy = "warn",
      sample_size_summary_statistic = "median",
      sample_size_relative_range_warning_threshold = 0.05,
      maximum_variants_per_locus = 30000,
      ld_peak_matrix_multiplier = 5,
      reserved_worker_memory_gb = 14
    ))
    sample_size_qc <- summarize_locus_nef(c(60, 100, 102), "warn", "median", 0.05)
    stopifnot(
      identical(sample_size_qc$warning_reason, "nef_relative_range_exceeds_threshold"),
      sample_size_qc$nef_min == 60,
      sample_size_qc$nef_max == 102,
      sample_size_qc$nef_median == 100,
      abs(sample_size_qc$nef_relative_range - 0.42) < 1e-12
    )
    stopifnot(identical(
      combine_warning_reasons(
        sample_size_qc$warning_reason, "ld_z_mismatch_warning"
      ),
      "nef_relative_range_exceeds_threshold;ld_z_mismatch_warning"
    ))

    fit_input <- create_susie_fit_input(
      selected$EZ, diag(2), 100, "valid_original", "fit_retry", quiet
    )
    fit_input_path <- fit_input$path
    serialized_fit_input <- readRDS(fit_input_path)
    stopifnot(
      identical(serialized_fit_input$z, as.numeric(selected$EZ)),
      identical(serialized_fit_input$R, diag(2)),
      identical(serialized_fit_input$n, 100),
      identical(serialized_fit_input$matrix_state, "valid_original")
    )
    fit_job <- list(
      locus = locus, selected = selected, ld.mat = diag(2), n_eff = 100,
      ld_qc_lambda = 0, ld_min_eigenvalue = 1,
      recovery_route = "fit_retry", audit = list(),
      fit_input = fit_input,
      sample_size_qc = sample_size_qc,
      resource_qc = summarize_locus_ld_resources(2, 30000, 5, 14)
    )
    fit_result <- recover_locus(
      fit_job, "study", {str(tmp_path)!r}, 2L, 1, 1, FALSE,
      defaults = defaults
    )
    stopifnot(
      identical(fit_result$reason, "fit_recovery_exhausted"),
      repair_calls == 0L,
      fit_calls == 2L,
      identical(fit_result$fit_input_qc$fit_input_serializations, 1L),
      identical(fit_result$fit_input_qc$fit_child_processes, 2L),
      identical(fit_result$fit_input_qc$fit_input_reuses, 1L),
      length(unique(fit_input_paths[1:2])) == 1L,
      all(fit_input_states[1:2] == "valid_original"),
      !file.exists(fit_input_path)
    )

    repair_job <- fit_job
    repair_job$ld.mat <- non_psd
    repair_job$ld_min_eigenvalue <- -1
    repair_job$recovery_route <- "ld_repair"
    repair_job$fit_input <- NULL
    repair_result <- recover_locus(
      repair_job, "study", {str(tmp_path)!r}, 2L, 1, 1, FALSE,
      defaults = defaults
    )
    stopifnot(
      repair_calls == 1L,
      fit_calls == 4L,
      identical(repair_result$reason, "repaired_ld_fit_recovery_exhausted"),
      identical(repair_result$fit_input_qc$fit_input_serializations, 1L),
      identical(repair_result$fit_input_qc$fit_child_processes, 2L),
      identical(repair_result$fit_input_qc$fit_input_reuses, 1L),
      length(unique(fit_input_paths[3:4])) == 1L,
      all(fit_input_states[3:4] == "valid_repaired"),
      !file.exists(fit_input_paths[[3L]]),
      any(vapply(repair_result$audit, function(x) identical(x$action, "validate_repaired_ld"), logical(1)))
    )

    audit_file <- file.path({str(tmp_path)!r}, "audit.tsv")
    aggregate_and_write_results(
      list(repair_result), "study", {str(tmp_path)!r}, FALSE,
      target_coverage = SUSIE_DEFAULTS$credible_set_coverage,
      recovery_audit_file = audit_file
    )
    qc_file <- file.path({str(tmp_path)!r}, "study_SuSiE_QC_summary.tsv")
    stopifnot(file.exists(audit_file), file.exists(qc_file))
    audit <- data.table::fread(audit_file)
    qc <- data.table::fread(qc_file)
    stopifnot(
      nrow(audit) == length(repair_result$audit),
      identical(qc$note[[1]], "repaired_ld_fit_recovery_exhausted"),
      identical(qc$converged[[1]], FALSE)
    )
    stopifnot(
      identical(qc$fit_input_serializations[[1]], 1L),
      identical(qc$fit_child_processes[[1]], 2L),
      identical(qc$fit_input_reuses[[1]], 1L),
      identical(qc$fit_input_matrix_state[[1]], "valid_repaired")
    )

    mismatch_fit_calls <- 0L
    wait_for_memory <- function(...) invisible(NULL)
    postgwas_ld_matrix <- function(...) diag(2)
    align_sumstats_to_ld <- function(selected, ld, log_msg) {{
      list(selected = selected, ld = ld)
    }}
    validate_ld_for_susie <- function(...) {{
      stop_susie_ld_validation(
        "ld_z_mismatch", "test severe LD/z mismatch",
        repairable = FALSE, lambda = 0.2, min_eigenvalue = 1
      )
    }}
    run_susie_hard_timeout <- function(...) {{
      mismatch_fit_calls <<- mismatch_fit_calls + 1L
      stop("SuSiE fitting must not run after severe LD/z mismatch")
    }}
    mismatch_df <- data.table::copy(selected)
    mismatch_df[, NEF := 100]
    mismatch_result <- process_locus(
      locus, mismatch_df, "study", "reference", "plink",
      {str(tmp_path)!r}, 7.3, 2L, 1, 1, FALSE, "6", 25000000,
      35000000, FALSE, defaults = defaults
    )
    stopifnot(
      identical(mismatch_result$status, "FAILED"),
      identical(mismatch_result$reason, "ld_z_mismatch"),
      mismatch_fit_calls == 0L,
      identical(mismatch_result$audit[[1]]$repairable, FALSE),
      identical(mismatch_result$ld_min_eigenvalue, 1)
    )

    ld_calls <- 0L
    postgwas_ld_matrix <- function(...) {{
      ld_calls <<- ld_calls + 1L
      stop("LD construction must not run for an oversized locus")
    }}
    oversized_defaults <- utils::modifyList(defaults, list(
      maximum_variants_per_locus = 1L
    ))
    oversized_result <- process_locus(
      locus, mismatch_df, "study", "reference", "plink",
      {str(tmp_path)!r}, 7.3, 2L, 1, 1, FALSE, "6", 25000000,
      35000000, FALSE, defaults = oversized_defaults
    )
    stopifnot(
      identical(oversized_result$status, "FAILED"),
      identical(oversized_result$reason, "maximum_variants_per_locus_exceeded"),
      identical(oversized_result$resource_qc$input_variant_count, 2L),
      grepl("--maximum-variants-per-locus", oversized_result$failure_detail),
      ld_calls == 0L,
      identical(oversized_result$audit[[1]]$ld_status, "not_constructed")
    )
    aggregate_and_write_results(
      list(oversized_result), "oversized", {str(tmp_path)!r}, FALSE,
      target_coverage = SUSIE_DEFAULTS$credible_set_coverage,
      recovery_audit_file = file.path({str(tmp_path)!r}, "oversized_audit.tsv")
    )
    oversized_qc <- data.table::fread(file.path(
      {str(tmp_path)!r}, "oversized_SuSiE_QC_summary.tsv"
    ))
    stopifnot(
      identical(oversized_qc$failure_reason[[1]],
                "maximum_variants_per_locus_exceeded"),
      identical(oversized_qc$resource_failure_reason[[1]],
                "maximum_variants_per_locus_exceeded"),
      identical(oversized_qc$input_variant_count[[1]], 2L),
      grepl("--maximum-variants-per-locus", oversized_qc$failure_detail[[1]])
    )

    annotated_cs <- annotate_sample_size_qc(
      data.table::data.table(SNP = "rs1", cs = 1L), sample_size_qc
    )
    annotated_cs <- annotate_resource_qc(
      annotated_cs, summarize_locus_ld_resources(2, 30000, 5, 14)
    )
    annotated_cs[, `:=`(
      analysis_status = "success", outcome_reason = "primary_fit_converged"
    )]
    success_result <- list(
      status = "OK", reason = "primary_fit_converged", locus = locus,
      vars = annotated_cs, cs = annotated_cs, sample_size_qc = sample_size_qc,
      resource_qc = summarize_locus_ld_resources(2, 30000, 5, 14),
      ld_qc_lambda = 0, ld_min_eigenvalue = 1, audit = list()
    )
    aggregate_and_write_results(
      list(success_result), "success", {str(tmp_path)!r}, FALSE,
      target_coverage = SUSIE_DEFAULTS$credible_set_coverage,
      recovery_audit_file = file.path({str(tmp_path)!r}, "success_audit.tsv")
    )
    combined_cs <- data.table::fread(file.path(
      {str(tmp_path)!r}, "success_SUSIE_combined_credibleset.csv"
    ))
    stopifnot(
      identical(combined_cs$warning_reason[[1]], "nef_relative_range_exceeds_threshold"),
      identical(combined_cs$outcome_reason[[1]], "primary_fit_converged"),
      identical(combined_cs$maximum_variants_per_locus[[1]], 30000L),
      abs(combined_cs$nef_relative_range[[1]] - 0.42) < 1e-12
    )
    """
    completed = subprocess.run(
        ["Rscript", "-e", script],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is unavailable")
def test_susie_reused_fit_payload_is_numerically_equivalent(tmp_path):
    availability = subprocess.run(
        [
            "Rscript",
            "-e",
            (
                "quit(status=if (requireNamespace('susieR', quietly=TRUE) && "
                "requireNamespace('processx', quietly=TRUE)) 0L else 42L)"
            ),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if availability.returncode != 0:
        pytest.skip("susieR and processx are required for numerical equivalence")

    susie_dir = (
        Path(__file__).parents[1]
        / "src/postgwas/modules/fine_mapping/engines/susie"
    )
    module = load_configuration().modules.fine_mapping
    susie = module.engines.susie
    script = f"""
    source({str(susie_dir / 'defaults.r')!r})
    source({str(susie_dir / 'utlities.r')!r})
    SUSIE_DEFAULTS <- utils::modifyList(SUSIE_DEFAULTS, list(
      credible_set_coverage = {module.credible_set_coverage!r},
      min_abs_corr = {susie.minimum_purity!r},
      process_poll_seconds = {susie.execution.process_poll_seconds!r},
      process_terminate_grace_seconds = {susie.execution.termination_grace_seconds!r},
      stderr_tail_lines = {susie.execution.stderr_tail_lines!r}
    ))
    quiet <- function(...) invisible(NULL)
    z <- c(4.0, 2.5, -1.5)
    ld <- matrix(c(
      1.0, 0.3, 0.1,
      0.3, 1.0, 0.2,
      0.1, 0.2, 1.0
    ), nrow = 3, byrow = TRUE)
    direct <- susieR::susie_rss(
      z = z, R = ld, n = 1000, L = 2L, max_iter = 1000L,
      coverage = SUSIE_DEFAULTS$credible_set_coverage,
      min_abs_corr = SUSIE_DEFAULTS$min_abs_corr,
      return_correlation = TRUE, verbose = FALSE
    )
    fit_input <- create_susie_fit_input(
      z, ld, 1000, "valid_original", "equivalence", quiet
    )
    input_path <- fit_input$path
    timed_out <- run_susie_hard_timeout(
      fit_input, 2L, 1000L, 1e-9, "equivalence_timeout", quiet,
      susie_verbose = FALSE
    )
    child_one <- run_susie_hard_timeout(
      fit_input, 2L, 1000L, 30, "equivalence_one", quiet,
      susie_verbose = FALSE
    )
    child_two <- run_susie_hard_timeout(
      fit_input, 2L, 1000L, 30, "equivalence_two", quiet,
      susie_verbose = FALSE
    )
    qc <- summarize_susie_fit_input(fit_input)
    stopifnot(
      is_susie_fit_success(child_one),
      is_susie_fit_success(child_two),
      inherits(timed_out, "susie_fit_failure"),
      identical(timed_out$reason, "susie_timeout"),
      isTRUE(all.equal(direct$pip, child_one$pip, tolerance = 1e-12)),
      isTRUE(all.equal(direct$alpha, child_one$alpha, tolerance = 1e-12)),
      isTRUE(all.equal(child_one$pip, child_two$pip, tolerance = 1e-12)),
      isTRUE(all.equal(child_one$alpha, child_two$alpha, tolerance = 1e-12)),
      identical(qc$fit_input_serializations, 1L),
      identical(qc$fit_child_processes, 3L),
      identical(qc$fit_input_reuses, 2L)
    )
    cleanup_susie_fit_input(fit_input, quiet)
    stopifnot(!file.exists(input_path))
    """
    completed = subprocess.run(
        ["Rscript", "-e", script],
        check=False,
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert completed.returncode == 0, completed.stderr


def test_invalid_finemap_probability_is_rejected_during_configuration(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "engines:\n  finemap:\n    corr_config: 1.1\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="corr_config"):
        _resolve_fine_mapping_arguments(
            argparse.Namespace(finemap_method="finemap", run_config=str(config_file))
        )


@pytest.mark.parametrize(
    "setting",
    [
        "ldstore_timeout_seconds",
        "finemap_timeout_seconds",
        "termination_grace_seconds",
    ],
)
def test_finemap_execution_timeouts_must_be_positive(tmp_path, setting):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        f"engines:\n  finemap:\n    {setting}: 0\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match=setting):
        _resolve_fine_mapping_arguments(
            argparse.Namespace(
                finemap_method="finemap", run_config=str(config_file)
            )
        )


def test_prior_k_fails_before_analysis_without_a_k_file(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "engines:\n  finemap:\n    prior_k: true\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="requires a K file"):
        _resolve_fine_mapping_arguments(
            argparse.Namespace(finemap_method="finemap", run_config=str(config_file))
        )


def _timeout_error(command, timeout_seconds, forced=True):
    return SupervisedProcessTimeout(
        command=command,
        timeout_seconds=timeout_seconds,
        elapsed_seconds=timeout_seconds + 0.1,
        stdout="partial stdout\n",
        stderr="partial stderr\n",
        forced_termination=forced,
    )


def test_finemap_timeout_is_typed_and_logged(tmp_path):
    master = tmp_path / "locus.master"
    master.write_text("master", encoding="utf-8")
    timeout = _timeout_error(["finemap"], 2)

    with patch(
        "postgwas.modules.fine_mapping.engines.finemap.runner.run_supervised_process",
        side_effect=timeout,
    ):
        with pytest.raises(FinemapExternalToolTimeout) as caught:
            run_finemap_binary(
                master,
                _finemap_config(),
                finemap_binary="finemap",
                timeout_seconds=2,
                termination_grace_seconds=1,
            )

    assert caught.value.reason == "finemap_timeout"
    assert caught.value.stage == "FINEMAP"
    assert caught.value.forced_termination is True
    log = master.with_suffix(".finemap.log").read_text(encoding="utf-8")
    assert "Status: TIMEOUT" in log
    assert "partial stderr" in log
    assert "Forced Termination: True" in log


def test_finemap_no_causal_configuration_remains_a_locus_outcome(tmp_path):
    master = tmp_path / "locus.master"
    master.write_text("master", encoding="utf-8")
    completed = subprocess.CompletedProcess(
        [],
        1,
        stdout="No causal configuration found\n",
        stderr="",
    )

    with patch(
        "postgwas.modules.fine_mapping.engines.finemap.runner.run_supervised_process",
        return_value=completed,
    ):
        success, reason = run_finemap_binary(
            master,
            _finemap_config(),
            finemap_binary="finemap",
            timeout_seconds=2,
            termination_grace_seconds=1,
        )

    assert success is False
    assert reason == "No causal SNPs found"


def test_ldstore_timeout_is_typed_and_logged(tmp_path):
    master = tmp_path / "locus.ldstore.master"
    master.write_text("master", encoding="utf-8")
    timeout = _timeout_error(["ldstore"], 2, forced=False)

    with patch(
        "postgwas.modules.fine_mapping.engines.finemap.runner.run_supervised_process",
        side_effect=timeout,
    ):
        with pytest.raises(FinemapExternalToolTimeout) as caught:
            run_ldstore(
                master,
                ldstore_binary="ldstore",
                timeout_seconds=2,
                termination_grace_seconds=1,
            )

    assert caught.value.reason == "ldstore_timeout"
    assert caught.value.stage == "LDstore"
    assert caught.value.forced_termination is False
    log = master.with_suffix(".ldstore.log").read_text(encoding="utf-8")
    assert "Status: TIMEOUT" in log


def test_ldstore_timeout_marks_only_the_locus_and_removes_partial_outputs(
    tmp_path,
):
    source_z = tmp_path / "source.z"
    source_snp = tmp_path / "source.snp"
    source_z.write_text("z\n", encoding="utf-8")
    source_snp.write_text("rs1\n", encoding="utf-8")
    locus_id = "chr1_1_2"
    locus_dir = tmp_path / "loci" / locus_id

    def fake_plink_extraction(_reference, _snps, out_prefix, **_kwargs):
        out_prefix.with_suffix(".fam").write_text(
            "F I 0 0 0 -9\n", encoding="utf-8"
        )

    def timeout_ldstore(*_args, **_kwargs):
        locus_dir.joinpath(f"{locus_id}.bcor").write_text(
            "partial", encoding="utf-8"
        )
        locus_dir.joinpath(f"{locus_id}.ld").write_text(
            "partial", encoding="utf-8"
        )
        raise FinemapExternalToolTimeout(
            "ldstore_timeout",
            "LDstore",
            _timeout_error(["ldstore"], 2),
            locus_dir / f"{locus_id}.ldstore.log",
        )

    task = {
        "locus_id": locus_id,
        "z_file": source_z,
        "snp_file": source_snp,
        "out_dir": tmp_path / "loci",
        "ld_ref": "reference",
        "n_samples": 100,
        "finemap_config": _finemap_config(),
        "ldstore_timeout_seconds": 2,
        "finemap_timeout_seconds": 3,
        "termination_grace_seconds": 1,
        "plink": "plink2",
        "bgenix": "bgenix",
        "ldstore": "ldstore",
        "finemap_executable": "finemap",
        "genomic_locus": "chr1:1-2",
        "chromosome": "1",
        "start": 1,
        "end": 2,
        "runtime_defaults": {},
        "file_log_level": "INFO",
        "sample_size_qc": {},
        "resource_qc": {},
    }

    with (
        patch(
            "postgwas.modules.fine_mapping.engines.finemap.adapter."
            "run_plink_extraction",
            side_effect=fake_plink_extraction,
        ),
        patch(
            "postgwas.modules.fine_mapping.engines.finemap.adapter."
            "_reconcile_extracted_variants",
            return_value=1,
        ),
        patch(
            "postgwas.modules.fine_mapping.engines.finemap.adapter."
            "run_plink_to_bgen"
        ),
        patch(
            "postgwas.modules.fine_mapping.engines.finemap.adapter."
            "run_bgen_indexing"
        ),
        patch(
            "postgwas.modules.fine_mapping.engines.finemap.adapter."
            "create_ldstore_master"
        ),
        patch(
            "postgwas.modules.fine_mapping.engines.finemap.adapter.run_ldstore",
            side_effect=timeout_ldstore,
        ),
    ):
        result = process_single_locus(task)

    assert result["status"] == "failed"
    assert result["failure_reason"] == "ldstore_timeout"
    assert result["timeout_stage"] == "LDstore"
    assert result["timeout_seconds"] == 2
    assert result["timeout_forced_termination"] is True
    assert not (locus_dir / f"{locus_id}.bcor").exists()
    assert not (locus_dir / f"{locus_id}.ld").exists()
    assert (locus_dir / f"{locus_id}_QC.tsv").is_file()


@pytest.mark.skipif(os.name != "posix", reason="POSIX process groups required")
def test_supervised_timeout_terminates_descendant_process_tree(tmp_path):
    child_pid_file = tmp_path / "child.pid"
    parent_script = (
        "import pathlib, subprocess, sys, time; "
        "child=subprocess.Popen([sys.executable, '-c', "
        "'import signal, time; signal.signal(signal.SIGTERM, "
        "signal.SIG_IGN); time.sleep(60)']); "
        f"pathlib.Path({str(child_pid_file)!r}).write_text(str(child.pid)); "
        "time.sleep(60)"
    )

    with pytest.raises(SupervisedProcessTimeout) as caught:
        run_supervised_process(
            [sys.executable, "-c", parent_script],
            timeout_seconds=0.5,
            termination_grace_seconds=0.2,
        )

    assert caught.value.forced_termination is True
    child_pid = int(child_pid_file.read_text(encoding="utf-8"))
    deadline = time.monotonic() + 2
    state = ""
    while time.monotonic() < deadline:
        status = subprocess.run(
            ["ps", "-o", "stat=", "-p", str(child_pid)],
            capture_output=True,
            text=True,
            check=False,
        )
        state = status.stdout.strip()
        if not state or state.startswith("Z"):
            break
        time.sleep(0.05)
    assert not state or state.startswith("Z")


def test_sss_command_receives_every_applicable_configured_control(tmp_path):
    master = tmp_path / "locus.master"
    master.write_text("master", encoding="utf-8")
    completed = subprocess.CompletedProcess([], 0, stdout="", stderr="")

    with patch(
        "postgwas.modules.fine_mapping.engines.finemap.runner.run_supervised_process",
        return_value=completed,
    ) as run:
        run_finemap_binary(
            master,
            _finemap_config(),
            threads=3,
            finemap_binary="finemap",
            timeout_seconds=43200,
            termination_grace_seconds=30,
        )

    command = run.call_args.args[0]
    assert run.call_args.kwargs["timeout_seconds"] == 43200
    assert run.call_args.kwargs["termination_grace_seconds"] == 30
    assert command[:2] == ["finemap", "--sss"]
    assert command[command.index("--n-causal-snps") + 1] == "5"
    assert command[command.index("--n-conv-sss") + 1] == "100"
    assert command[command.index("--n-configs-top") + 1] == "50000"
    assert command[command.index("--corr-config") + 1] == "0.95"
    assert "--pvalue-snps" not in command
    assert "native unfiltered endpoint" in master.with_suffix(".finemap.log").read_text()
    assert command[command.index("--prior-std") + 1] == "0.05"
    assert "--cond-pvalue" not in command


def test_conditional_command_excludes_sss_controls(tmp_path):
    master = tmp_path / "locus.master"
    master.write_text("master", encoding="utf-8")
    completed = subprocess.CompletedProcess([], 0, stdout="", stderr="")
    config = _finemap_config(
        algorithm="cond",
        cond_pvalue=1e-6,
        force_n_samples=True,
        std_effects=True,
    )

    with patch(
        "postgwas.modules.fine_mapping.engines.finemap.runner.run_supervised_process",
        return_value=completed,
    ) as run:
        run_finemap_binary(
            master,
            config,
            threads=2,
            finemap_binary="finemap",
            timeout_seconds=43200,
            termination_grace_seconds=30,
        )

    command = run.call_args.args[0]
    assert command[:2] == ["finemap", "--cond"]
    assert command[command.index("--cond-pvalue") + 1] == "1e-06"
    assert "--n-iter" not in command
    assert "--n-configs-top" not in command
    assert "--corr-config" not in command
    assert "--force-n-samples" in command
    assert "--std-effects" in command


def test_explicit_credible_set_coverage_is_validated_and_recorded(tmp_path):
    inputs = tmp_path / "selected"
    outputs = tmp_path / "flames"
    inputs.mkdir()
    (inputs / "locus.z").write_text(
        "rsid chromosome position allele1 allele2\n"
        "rs1 1 100 A G\n",
        encoding="utf-8",
    )
    (inputs / "locus.cred1").write_text(
        "# Post-Pr = 1\n"
        "cred1 prob1\n"
        "rs1 0.9\n",
        encoding="utf-8",
    )

    create_flames_file(
        inputs,
        outputs,
        outputs / "annotations",
        "indexfile.txt",
        "finemap_FLAMES_manifest.tsv",
        "annotated_",
        target_coverage=0.9,
        locus_metadata={
            "locus": {
                "GenomicLocus": "chr1:1-200",
                "analysis_status": "success",
                "outcome_reason": "finemap_completed",
                "warning_reason": "nef_relative_range_exceeds_threshold",
                "failure_reason": None,
                "nef_min": 60,
                "nef_max": 102,
                "nef_median": 100,
                "nef_relative_range": 0.42,
                "selected_nef": 100,
                "nef_policy": "warn",
                "nef_summary_statistic": "median",
                "nef_warning_threshold": 0.05,
            }
        },
    )

    manifest = (outputs / "finemap_FLAMES_manifest.tsv").read_text(
        encoding="utf-8"
    )
    assert "target_coverage" in manifest
    assert "\t0.9\t0.9\t" in manifest
    assert "finemap_0.9_coverage_credible_set" in manifest
    assert "nef_relative_range_exceeds_threshold" in manifest
    assert "nef_relative_range" in manifest
    assert "outcome_reason" in manifest
