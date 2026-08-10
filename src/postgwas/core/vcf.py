"""Reusable, checked access to single-sample VCF data."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Mapping, Type

from postgwas.core.processes import run_checked_command


VCF_TAG = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


class VcfQueryError(RuntimeError):
    """A VCF could not be queried without ambiguity."""


def select_vcf_sample(
    vcf_path: str | Path,
    dataset_id: str,
    bcftools: str,
    *,
    logger=None,
    error_type: Type[Exception] = VcfQueryError,
) -> str:
    """Select ``dataset_id`` or the only sample; never guess among samples."""
    samples = [
        value
        for value in run_checked_command(
            [bcftools, "query", "-l", str(vcf_path)],
            "Reading VCF sample names",
            logger=logger,
            error_type=error_type,
        ).splitlines()
        if value
    ]
    if dataset_id in samples:
        return dataset_id
    if len(samples) == 1:
        sample = samples[0]
        if logger is not None and sample != dataset_id:
            logger.warning(
                "VCF sample %r differs from dataset_id %r; the only sample was used."
                % (sample, dataset_id)
            )
        return sample
    raise error_type(
        "VCF must contain dataset sample %r or exactly one sample; found %s."
        % (dataset_id, ", ".join(samples) or "none")
    )


def extract_vcf_table(
    vcf_path: str | Path,
    table_path: str | Path,
    dataset_id: str,
    columns: Mapping[str, str],
    bcftools: str,
    *,
    delimiter: str,
    io_buffer_bytes: int,
    include_expression: str | None = None,
    allow_undefined_tags: bool = False,
    logger=None,
    error_type: Type[Exception] = VcfQueryError,
    purpose: str = "Extracting VCF fields",
) -> str:
    """Extract one tabular VCF projection with an explicit header.

    ``columns`` maps output column names to bcftools query expressions. FORMAT
    expressions must include their brackets, for example ``[%ES]``.
    """
    vcf = Path(vcf_path).expanduser().resolve()
    destination = Path(table_path).expanduser().resolve()
    if not vcf.is_file() or vcf.stat().st_size <= 0:
        raise error_type("VCF does not exist or is empty: %s" % vcf)
    if not columns:
        raise error_type("At least one VCF field is required for extraction.")

    sample = select_vcf_sample(
        vcf, dataset_id, bcftools, logger=logger, error_type=error_type,
    )
    if len(delimiter) != 1:
        raise error_type("VCF table delimiter must be exactly one character.")
    if io_buffer_bytes <= 0:
        raise error_type("VCF table I/O buffer must be greater than zero.")
    query = delimiter.join(columns.values()) + "\n"
    arguments = [bcftools, "query"]
    if allow_undefined_tags:
        arguments.append("--allow-undef-tags")
    arguments.extend(["--samples", sample])
    if include_expression:
        arguments.extend(["--include", include_expression])
    arguments.extend(["--format", query, str(vcf)])

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".body")
    try:
        run_checked_command(
            arguments,
            purpose,
            logger=logger,
            error_type=error_type,
            stdout_path=temporary,
        )
        with destination.open("w", encoding="utf-8", newline="") as output:
            output.write(delimiter.join(columns) + "\n")
            with temporary.open("r", encoding="utf-8") as body:
                for block in iter(lambda: body.read(io_buffer_bytes), ""):
                    output.write(block)
    except OSError as exc:
        raise error_type("Cannot write extracted VCF table %s: %s" % (destination, exc)) from exc
    finally:
        temporary.unlink(missing_ok=True)

    if logger is not None and hasattr(logger, "record"):
        logger.record("OUTPUT", "vcf_table", path=str(destination), sample=sample)
    return str(destination)


__all__ = [
    "VCF_TAG",
    "VcfQueryError",
    "extract_vcf_table",
    "select_vcf_sample",
]
