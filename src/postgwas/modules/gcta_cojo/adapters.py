"""Shell-free GCTA-COJO command construction and checked execution."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

from postgwas.core.gcta import require_supported_gcta as _require_supported_gcta
from postgwas.core.processes import run_checked_command
from postgwas.modules.gcta_cojo.errors import GctaCojoError


def require_supported_gcta(executable, module, logger, timeout_seconds) -> str:
    return _require_supported_gcta(
        executable,
        module,
        logger,
        timeout_seconds,
        error_type=GctaCojoError,
        configuration_path="modules.gcta_cojo",
        analysis_name="GCTA-COJO",
    )


def build_cojo_command(
    executable: str,
    summary_statistics: str | Path,
    reference_prefix: str | Path,
    output_prefix: str | Path,
    threads: int,
    module,
) -> list[str]:
    """Build one GCTA-COJO command from the resolved canonical configuration."""
    analysis = module.analysis
    inputs = module.inputs
    command = [
        executable,
        "--bfile", str(reference_prefix),
        "--cojo-file", str(summary_statistics),
        "--maf", str(analysis.reference_maf_min),
        "--diff-freq", str(analysis.frequency_difference_max),
        "--cojo-wind", str(analysis.window_kb),
        "--cojo-collinear", str(analysis.collinearity_cutoff),
        "--thread-num", str(threads),
        "--out", str(output_prefix),
    ]
    if analysis.chromosome is not None:
        command.extend(["--chr", analysis.chromosome])
    if inputs.exclude_snps is not None:
        command.extend([
            "--exclude", str(Path(inputs.exclude_snps).expanduser().resolve()),
        ])
    if inputs.extract_snps is not None:
        command.extend([
            "--extract", str(Path(inputs.extract_snps).expanduser().resolve()),
        ])
    if analysis.genomic_control:
        command.append("--cojo-gc")
        if analysis.genomic_control_lambda is not None:
            command.append(str(analysis.genomic_control_lambda))

    if module.mode == "slct":
        command.extend(["--cojo-slct", "--cojo-p", str(analysis.significance_threshold)])
    elif module.mode == "top_snps":
        command.extend(["--cojo-top-SNPs", str(analysis.top_snp_count)])
    elif module.mode == "joint":
        command.extend([
            "--extract", str(Path(inputs.joint_snps).expanduser().resolve()),
            "--cojo-joint",
        ])
    else:
        command.extend([
            "--cojo-cond", str(Path(inputs.condition_snps).expanduser().resolve()),
        ])
    return command


def run_cojo_command(
    command: Sequence[str],
    expected_result: Path,
    configuration,
    logger,
    *,
    dry_run: bool,
) -> None:
    run_checked_command(
        command,
        "GCTA-COJO %s analysis" % configuration.modules.gcta_cojo.mode,
        logger=logger,
        error_type=GctaCojoError,
        timeout_seconds=configuration.execution.timeout_seconds,
        expected_outputs=[expected_result],
        dry_run=dry_run,
    )


__all__ = ["build_cojo_command", "require_supported_gcta", "run_cojo_command"]
