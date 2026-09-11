"""Validated command construction and execution for CBIIT LDSC."""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any

from postgwas.core.paths import require_nonempty_file
from postgwas.core.processes import run_checked_command
from postgwas.core.reference_resources import require_file_inventory
from postgwas.core.validation_reporting import register_file_availability_bundle


class LDSCError(RuntimeError):
    """Raised when an LDSC input, command, or required result is invalid."""


@dataclass(frozen=True)
class LDSCReferenceValidation:
    """Exact external files consumed by one resolved LDSC analysis."""

    merge_alleles: Path | None
    reference_directory: Path
    weights_directory: Path
    required_files: tuple[Path, ...]


_ESTIMATE_WITH_SE = re.compile(
    r"^(?P<estimate>[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
    r"\s+\((?P<standard_error>[-+]?(?:\d+(?:\.\d*)?|\.\d+)"
    r"(?:[eE][-+]?\d+)?)\)$"
)


def _validate_estimate(value: str, label: str, *, constrained: bool = False) -> None:
    if constrained and value.startswith("constrained to "):
        try:
            estimate = float(value.removeprefix("constrained to "))
        except ValueError as exc:
            raise LDSCError("LDSC %s is not numeric: %s" % (label, value)) from exc
        if not math.isfinite(estimate):
            raise LDSCError("LDSC %s is not finite: %s" % (label, value))
        return
    match = _ESTIMATE_WITH_SE.fullmatch(value)
    if match is None:
        raise LDSCError(
            "LDSC %s does not contain an estimate and standard error: %s"
            % (label, value)
        )
    estimate = float(match.group("estimate"))
    standard_error = float(match.group("standard_error"))
    if not math.isfinite(estimate) or not math.isfinite(standard_error):
        raise LDSCError("LDSC %s is not finite: %s" % (label, value))
    if standard_error < 0:
        raise LDSCError("LDSC %s has a negative standard error: %s" % (label, value))


def build_munge_sumstats_command(
    executable: str,
    *,
    sumstats: str | Path,
    output_prefix: str | Path,
    merge_alleles: str | Path,
    minimum_info: float,
    minimum_maf: float,
    minimum_n: float | None,
    chunksize: int,
    keep_maf: bool,
) -> list[str]:
    """Build the shared HapMap3 munging command from resolved configuration."""
    command = [
        str(executable),
        "--sumstats", str(sumstats),
        "--out", str(output_prefix),
        "--merge-alleles", str(merge_alleles),
        "--info-min", str(minimum_info),
        "--maf-min", str(minimum_maf),
        "--chunksize", str(chunksize),
    ]
    if minimum_n is not None:
        command.extend(["--n-min", str(minimum_n)])
    if keep_maf:
        command.append("--keep-maf")
    return command


def build_h2_command(
    executable: str,
    *,
    sumstats: str | Path,
    reference_ld_directory: str | Path,
    weights_ld_directory: str | Path,
    output_prefix: str | Path,
    intercept: float | None,
    two_step: float | None,
    chisq_max: float | None,
    n_blocks: int,
    use_m_5_50: bool,
    print_covariance: bool,
    print_delete_values: bool,
    sample_prevalence: float | None = None,
    population_prevalence: float | None = None,
) -> list[str]:
    """Build one single-trait ``ldsc.py --h2`` command."""
    if (sample_prevalence is None) != (population_prevalence is None):
        raise LDSCError(
            "An LDSC h2 command requires both sample and population prevalence"
        )
    if intercept is not None and two_step is not None:
        raise LDSCError(
            "An LDSC h2 command cannot combine a constrained intercept and two-step"
        )
    command = [
        str(executable),
        "--h2", str(sumstats),
        "--ref-ld-chr", _chromosome_prefix(reference_ld_directory),
        "--w-ld-chr", _chromosome_prefix(weights_ld_directory),
        "--out", str(output_prefix),
        "--n-blocks", str(n_blocks),
    ]
    if intercept is not None:
        command.extend(["--intercept-h2", str(intercept)])
    if two_step is not None:
        command.extend(["--two-step", str(two_step)])
    if chisq_max is not None:
        command.extend(["--chisq-max", str(chisq_max)])
    if not use_m_5_50:
        command.append("--not-M-5-50")
    if print_covariance:
        command.append("--print-cov")
    if print_delete_values:
        command.append("--print-delete-vals")
    if sample_prevalence is not None and population_prevalence is not None:
        command.extend([
            "--samp-prev", str(sample_prevalence),
            "--pop-prev", str(population_prevalence),
        ])
    return command


def build_h2_cts_command(
    executable: str,
    *,
    sumstats: str | Path,
    baseline_ld_prefixes: list[str] | tuple[str, ...],
    weights_ld_prefix: str,
    ldcts_file: str | Path,
    output_prefix: str | Path,
    prefix_separator: str,
) -> list[str]:
    """Build the official LDSC ``--h2-cts`` regression command."""
    return [
        str(executable),
        "--h2-cts", str(sumstats),
        "--ref-ld-chr", prefix_separator.join(baseline_ld_prefixes),
        "--w-ld-chr", str(weights_ld_prefix),
        "--ref-ld-chr-cts", str(ldcts_file),
        "--out", str(output_prefix),
    ]


def _chromosome_prefix(directory: str | Path) -> str:
    """Render the directory prefix expected by CBIIT's 1..22 reader."""
    return str(Path(directory)) + os.sep


def _required_directory(value: str | Path, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise LDSCError("%s directory does not exist: %s" % (label, path))
    return path


def validate_reference_files(
    reference_directory: str | Path,
    weights_directory: str | Path,
    configuration,
) -> tuple[Path, Path]:
    """Require every chromosome-split file consumed by the pinned LDSC reader."""
    validated = validate_ldsc_reference_resources(
        None,
        reference_directory,
        weights_directory,
        configuration,
    )
    return validated.reference_directory, validated.weights_directory


def validate_ldsc_reference_resources(
    merge_alleles: str | Path | None,
    reference_directory: str | Path,
    weights_directory: str | Path,
    configuration,
) -> LDSCReferenceValidation:
    """Validate and enumerate the resolved external LDSC resource contract."""
    reference = _required_directory(reference_directory, "Reference LD-score")
    weights = _required_directory(weights_directory, "Regression-weight LD-score")
    layout = configuration.reference_layout
    m_suffix = layout.m_5_50_suffix if configuration.use_m_5_50 else layout.m_suffix
    required_files: list[Path] = []
    reference_ld_files: list[Path] = []
    reference_snp_count_files: list[Path] = []
    weight_ld_files: list[Path] = []
    for chromosome in range(1, configuration.chromosomes + 1):
        reference_ld_files.append(
            reference / (str(chromosome) + layout.ld_score_suffix)
        )
        reference_snp_count_files.append(
            reference / (str(chromosome) + m_suffix)
        )
        weight_ld_files.append(
            weights / (str(chromosome) + layout.ld_score_suffix)
        )
        required_files.extend((
            reference_ld_files[-1],
            reference_snp_count_files[-1],
            weight_ld_files[-1],
        ))
    require_file_inventory(
        required_files, "LDSC chromosome-split reference",
        missing_message="LDSC chromosome-split reference files are missing or empty",
        error_type=LDSCError,
    )
    directory_fields = (
        (("info", "ldsc_shared_reference_directory", reference, True),)
        if reference == weights
        else (
            ("info", "ldsc_reference_directory", reference, True),
            ("info", "ldsc_weight_directory", weights, True),
        )
    )
    register_file_availability_bundle(
        required_files,
        "LDSC chromosome reference bundle",
        (
            ("count", "chromosomes", [
                str(chromosome)
                for chromosome in range(1, configuration.chromosomes + 1)
            ]),
            ("success", "ldsc_reference_ld_files", "%d / %d" % (
                len(reference_ld_files), configuration.chromosomes,
            )),
            ("success", "ldsc_reference_snp_count_files", "%d / %d" % (
                len(reference_snp_count_files), configuration.chromosomes,
            )),
            ("success", "ldsc_weight_ld_files", "%d / %d" % (
                len(weight_ld_files), configuration.chromosomes,
            )),
            ("info", "ldsc_snp_count_suffix", m_suffix),
            ("count", "ldsc_unique_bundle_files", len(set(required_files))),
            *directory_fields,
        ),
    )
    merge_file = (
        None
        if merge_alleles is None
        else require_nonempty_file(
            merge_alleles, "LDSC merge-alleles file", error_type=LDSCError,
        )
    )
    return LDSCReferenceValidation(
        merge_alleles=merge_file,
        reference_directory=reference,
        weights_directory=weights,
        required_files=tuple(required_files),
    )


def extract_ldsc_metrics(log_file: str | Path) -> dict[str, str]:
    """Parse required h²/intercept and optional attenuation ratio from an LDSC log."""
    path = require_nonempty_file(
        log_file, "LDSC heritability log", error_type=LDSCError,
    )
    metrics: dict[str, str | None] = {
        "h2": None,
        "intercept": None,
        "ratio": None,
    }
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if line.startswith((
                    "Total Liability scale h2:", "Total Observed scale h2:",
                )):
                    metrics["h2"] = line.split(":", 1)[1].strip()
                elif line.startswith("Intercept:"):
                    metrics["intercept"] = line.split(":", 1)[1].strip()
                elif line.startswith("Ratio:"):
                    metrics["ratio"] = line.split(":", 1)[1].strip()
                elif line.startswith("Ratio < 0"):
                    metrics["ratio"] = line
    except OSError as exc:
        raise LDSCError(
            "Cannot read LDSC heritability log %s: %s" % (path, exc)
        ) from exc
    missing = [name for name in ("h2", "intercept") if metrics[name] is None]
    if missing:
        raise LDSCError(
            "LDSC log is missing required result metric(s) %s: %s"
            % (", ".join(missing), path)
        )
    _validate_estimate(str(metrics["h2"]), "heritability")
    _validate_estimate(str(metrics["intercept"]), "intercept", constrained=True)
    if (
        metrics["ratio"] is not None
        and not str(metrics["ratio"]).startswith(("NA", "Ratio < 0"))
    ):
        _validate_estimate(str(metrics["ratio"]), "attenuation ratio")
    return {
        "h2": str(metrics["h2"]),
        "intercept": str(metrics["intercept"]),
        "ratio": "N/A" if metrics["ratio"] is None else str(metrics["ratio"]),
    }


def _validate_liability_conversion_invariant(
    observed_metrics: dict[str, str],
    liability_metrics: dict[str, str],
    *,
    observed_log: Path,
    liability_log: Path,
    logger=None,
) -> None:
    """Require prevalence to leave the fitted LDSC regression unchanged.

    The pinned CBIIT LDSC implementation passes prevalence only to the summary
    conversion after fitting the regression. It may therefore rescale h² and
    its standard error, but the intercept and attenuation ratio must be shared
    by the observed- and liability-scale reports.
    """
    shared_metrics = ("intercept", "ratio")
    mismatches = [
        name for name in shared_metrics
        if observed_metrics[name] != liability_metrics[name]
    ]
    differences = "; ".join(
        "%s (observed=%r, liability=%r)"
        % (name, observed_metrics[name], liability_metrics[name])
        for name in mismatches
    )
    unverifiable_ratio_scales = [
        scale
        for scale, metrics in (
            ("observed", observed_metrics),
            ("liability", liability_metrics),
        )
        if (
            metrics["ratio"] == "N/A"
            and not metrics["intercept"].startswith("constrained to ")
        )
    ]
    problems = []
    if mismatches:
        problems.append(
            "the prevalence run changed shared regression metric(s): %s"
            % differences
        )
    if unverifiable_ratio_scales:
        problems.append(
            "the attenuation ratio is missing despite an unconstrained "
            "intercept in the %s log(s)"
            % ", ".join(unverifiable_ratio_scales)
        )
    if problems:
        validation_problem = "; ".join(problems)
        if logger is not None and hasattr(logger, "record"):
            logger.record(
                "FAILED", "ldsc_liability_conversion_invariant",
                mismatched_metrics=mismatches,
                unverifiable_ratio_scales=unverifiable_ratio_scales,
                observed_log=str(observed_log),
                liability_log=str(liability_log),
                problem=validation_problem,
            )
        raise LDSCError(
            "LDSC liability-scale validation failed because %s. In the "
            "pinned CBIIT LDSC implementation, --samp-prev and --pop-prev "
            "may rescale only h² and its standard error. PostGWAS compared "
            "observed log %s with liability log %s and stopped before "
            "publication."
            % (validation_problem, observed_log, liability_log)
        )
    if logger is not None and hasattr(logger, "record"):
        logger.record(
            "PASS", "ldsc_liability_conversion_invariant",
            result="intercept_and_ratio_unchanged",
            shared_intercept=observed_metrics["intercept"],
            shared_ratio=observed_metrics["ratio"],
            observed_log=str(observed_log),
            liability_log=str(liability_log),
        )


def _analysis_outputs(
    prefix: Path, configuration, *, include_all_optional: bool = False,
) -> list[Path]:
    layout = configuration.output_layout
    outputs = [Path(str(prefix) + layout.upstream_log_suffix)]
    if include_all_optional or configuration.print_covariance:
        outputs.append(Path(str(prefix) + layout.covariance_suffix))
    if include_all_optional or configuration.print_delete_values:
        outputs.extend([
            Path(str(prefix) + layout.delete_values_suffix),
            Path(str(prefix) + layout.partitioned_delete_values_suffix),
        ])
    return outputs


def _run_h2_analysis(
    command: list[str],
    purpose: str,
    outputs: list[Path],
    *,
    logger=None,
) -> dict[str, str]:
    """Run one h² command and validate its durable scientific report.

    The pinned CBIIT logger sends each report line to both stdout and ``.log``.
    On affected executions its unclosed buffered file can remain zero bytes
    despite a complete zero-exit stdout report. In that one condition, preserve
    the exact captured stdout as the staged log and apply the normal strict
    scientific parser before accepting the run.
    """
    log_path = outputs[0]
    stdout = run_checked_command(
        command,
        purpose,
        logger=logger,
        error_type=LDSCError,
        expected_outputs=outputs[1:],
    )
    if log_path.is_file() and log_path.stat().st_size > 0:
        return extract_ldsc_metrics(log_path)
    if not stdout:
        raise LDSCError(
            "%s returned success but expected output is missing or empty: %s. "
            "Captured stdout was also empty, so no scientific result can be "
            "validated."
            % (purpose, log_path)
        )
    try:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(stdout, encoding="utf-8")
        metrics = extract_ldsc_metrics(log_path)
    except (OSError, LDSCError) as exc:
        log_path.unlink(missing_ok=True)
        raise LDSCError(
            "%s returned success but expected output is missing or empty: %s. "
            "Captured stdout was not a complete, parseable LDSC heritability "
            "report: %s"
            % (purpose, log_path, exc)
        ) from exc
    if logger is not None and hasattr(logger, "record"):
        logger.record(
            "WARNING", "ldsc_log_recovered_from_stdout",
            purpose=purpose,
            path=str(log_path),
            captured_stdout_bytes=len(stdout.encode("utf-8")),
            reason="upstream_log_missing_or_empty_after_zero_exit",
            validation="required_heritability_metrics_parsed",
        )
    return metrics


def ldsc_owned_output_paths(
    output_prefix: str | Path, configuration,
) -> list[Path]:
    """Return every configured scientific path that one LDSC run may own."""
    prefix = Path(output_prefix)
    layout = configuration.output_layout
    paths = [
        Path(str(prefix) + layout.munged_sumstats_suffix),
        Path(str(prefix) + layout.upstream_log_suffix),
    ]
    for prefix_suffix in (
        layout.observed_prefix_suffix, layout.liability_prefix_suffix,
    ):
        paths.extend(_analysis_outputs(
            Path(str(prefix) + prefix_suffix),
            configuration,
            include_all_optional=True,
        ))
    return paths


def run_ldsc(
    *,
    sumstats_tsv: str | Path,
    output_prefix: str | Path,
    merge_alleles: str | Path,
    reference_ld_directory: str | Path,
    weights_ld_directory: str | Path,
    munge_executable: str,
    ldsc_executable: str,
    configuration,
    logger=None,
    reference_validation: LDSCReferenceValidation | None = None,
) -> dict[str, Any]:
    """Munge formatter output and run validated observed/liability h² analyses."""
    sumstats = require_nonempty_file(
        sumstats_tsv, "LDSC formatter input", error_type=LDSCError,
    )
    validated_reference = reference_validation or validate_ldsc_reference_resources(
        merge_alleles,
        reference_ld_directory,
        weights_ld_directory,
        configuration,
    )
    merge_file = validated_reference.merge_alleles
    if (
        merge_file is None
        or not merge_file.is_file()
        or merge_file.stat().st_size <= 0
    ):
        raise LDSCError(
            "Prevalidated LDSC merge-alleles file is missing or empty: %s"
            % merge_file
        )
    reference = validated_reference.reference_directory
    weights = validated_reference.weights_directory
    if logger is not None and hasattr(logger, "record"):
        logger.record(
            "INPUT", "ldsc_inputs",
            sumstats=str(sumstats),
            merge_alleles=str(merge_file),
            reference_ld_directory=str(reference),
            weights_ld_directory=str(weights),
        )
        logger.record(
            "PASS", "ldsc_reference_validation",
            chromosomes=configuration.chromosomes,
            m_file_suffix=(
                configuration.reference_layout.m_5_50_suffix
                if configuration.use_m_5_50
                else configuration.reference_layout.m_suffix
            ),
        )
    prefix = Path(output_prefix).expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    layout = configuration.output_layout

    munged_sumstats = Path(str(prefix) + layout.munged_sumstats_suffix)
    munge_log = Path(str(prefix) + layout.upstream_log_suffix)
    run_checked_command(
        build_munge_sumstats_command(
            munge_executable,
            sumstats=sumstats,
            output_prefix=prefix,
            merge_alleles=merge_file,
            minimum_info=configuration.minimum_info,
            minimum_maf=configuration.minimum_maf,
            minimum_n=configuration.minimum_n,
            chunksize=configuration.chunksize,
            keep_maf=configuration.keep_maf,
        ),
        "LDSC munge_sumstats",
        logger=logger,
        error_type=LDSCError,
        expected_outputs=[munged_sumstats, munge_log],
    )

    observed_prefix = Path(str(prefix) + layout.observed_prefix_suffix)
    observed_outputs = _analysis_outputs(observed_prefix, configuration)
    observed_metrics = _run_h2_analysis(
        build_h2_command(
            ldsc_executable,
            sumstats=munged_sumstats,
            reference_ld_directory=reference,
            weights_ld_directory=weights,
            output_prefix=observed_prefix,
            intercept=configuration.intercept,
            two_step=configuration.two_step,
            chisq_max=configuration.chisq_max,
            n_blocks=configuration.n_blocks,
            use_m_5_50=configuration.use_m_5_50,
            print_covariance=configuration.print_covariance,
            print_delete_values=configuration.print_delete_values,
        ),
        "LDSC observed-scale heritability",
        logger=logger,
        outputs=observed_outputs,
    )
    observed_log = observed_outputs[0]

    result: dict[str, Any] = {
        "munged_sumstats": str(munged_sumstats),
        "munge_log": str(munge_log),
        "h2_observed": str(observed_log),
        "observed_outputs": [str(path) for path in observed_outputs],
        "observed_metrics": observed_metrics,
        "h2_liability": None,
        "liability_outputs": [],
        "liability_metrics": None,
    }
    if configuration.population_prevalence is not None:
        liability_prefix = Path(str(prefix) + layout.liability_prefix_suffix)
        liability_outputs = _analysis_outputs(liability_prefix, configuration)
        liability_metrics = _run_h2_analysis(
            build_h2_command(
                ldsc_executable,
                sumstats=munged_sumstats,
                reference_ld_directory=reference,
                weights_ld_directory=weights,
                output_prefix=liability_prefix,
                intercept=configuration.intercept,
                two_step=configuration.two_step,
                chisq_max=configuration.chisq_max,
                n_blocks=configuration.n_blocks,
                use_m_5_50=configuration.use_m_5_50,
                print_covariance=configuration.print_covariance,
                print_delete_values=configuration.print_delete_values,
                sample_prevalence=configuration.sample_prevalence,
                population_prevalence=configuration.population_prevalence,
            ),
            "LDSC liability-scale heritability",
            logger=logger,
            outputs=liability_outputs,
        )
        liability_log = liability_outputs[0]
        _validate_liability_conversion_invariant(
            observed_metrics,
            liability_metrics,
            observed_log=observed_log,
            liability_log=liability_log,
            logger=logger,
        )
        result.update({
            "h2_liability": str(liability_log),
            "liability_outputs": [str(path) for path in liability_outputs],
            "liability_metrics": liability_metrics,
        })
    elif logger is not None and hasattr(logger, "record"):
        logger.record(
            "SKIP", "ldsc_liability_scale",
            reason="population_prevalence_not_provided",
        )
    return result


__all__ = [
    "LDSCError", "LDSCReferenceValidation", "build_h2_command", "build_h2_cts_command",
    "build_munge_sumstats_command", "extract_ldsc_metrics", "run_ldsc",
    "ldsc_owned_output_paths", "validate_ldsc_reference_resources",
    "validate_reference_files",
]
