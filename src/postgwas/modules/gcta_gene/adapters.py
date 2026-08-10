"""GCTA command construction and version validation."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

from postgwas.core.gcta import require_supported_gcta as _require_supported_gcta
from postgwas.core.processes import run_checked_command
from postgwas.modules.gcta_gene.errors import GctaGeneError


def require_supported_gcta(
    executable: str,
    module_config,
    logger,
    timeout_seconds: float,
) -> str:
    return _require_supported_gcta(
        executable,
        module_config,
        logger,
        timeout_seconds,
        error_type=GctaGeneError,
        configuration_path="modules.gcta_gene",
        analysis_name="fastBAT/mBAT-combo",
    )


def build_gcta_command(
    executable: str,
    input_file: str | Path,
    reference_prefix: str | Path,
    annotation_file: str | Path | None,
    output_prefix: str | Path,
    threads: int,
    module_config,
) -> list[str]:
    """Build one shell-free command following the official GCTA contract."""
    common = [
        executable,
        "--bfile", str(reference_prefix),
        "--maf", str(module_config.reference_maf_min),
        "--thread-num", str(threads),
        "--out", str(output_prefix),
    ]
    if module_config.method.startswith("fastbat_"):
        command = [
            *common,
            "--fastBAT", str(input_file),
            "--fastBAT-ld-cutoff", str(module_config.fastbat_ld_cutoff),
        ]
        if module_config.method == "fastbat_gene":
            command.extend([
                "--fastBAT-gene-list", str(annotation_file),
                "--fastBAT-wind", str(module_config.gene_window_kb),
            ])
        elif module_config.method == "fastbat_segment":
            command.extend(["--fastBAT-seg", str(module_config.segment_size_kb)])
        else:
            command.extend(["--fastBAT-set-list", str(annotation_file)])
        if module_config.write_snpset:
            command.append("--fastBAT-write-snpset")
        return command
    command = [
        *common,
        "--mBAT-combo", str(input_file),
        "--mBAT-gene-list", str(annotation_file),
        "--mBAT-wind", str(module_config.gene_window_kb),
        "--mBAT-svd-gamma", str(module_config.mbat_svd_gamma),
        "--diff-freq", str(module_config.frequency_difference_max),
        "--fastBAT-ld-cutoff", str(module_config.fastbat_ld_cutoff),
    ]
    if module_config.print_component_p_values:
        command.append("--mBAT-print-all-p")
    if module_config.write_snpset:
        command.append("--mBAT-write-snpset")
    return command


def run_gcta_command(
    command: Sequence[str],
    expected_result: Path,
    configuration,
    logger,
    *,
    dry_run: bool,
) -> None:
    run_checked_command(
        command,
        "%s association" % configuration.modules.gcta_gene.method,
        logger=logger,
        error_type=GctaGeneError,
        timeout_seconds=configuration.execution.timeout_seconds,
        expected_outputs=[expected_result],
        dry_run=dry_run,
    )


__all__ = ["build_gcta_command", "require_supported_gcta", "run_gcta_command"]
