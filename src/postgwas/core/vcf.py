"""Reusable, checked access to single-sample VCF data."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Iterable, Mapping, Sequence, Type

from postgwas.core.processes import run_checked_command


VCF_TAG = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


class VcfQueryError(RuntimeError):
    """A VCF could not be queried without ambiguity."""


def read_vcf_header(
    vcf_path: str | Path,
    bcftools: str,
    *,
    logger=None,
    error_type: Type[Exception] = VcfQueryError,
) -> str:
    """Read one complete VCF header with the configured bcftools executable."""
    header = run_checked_command(
        [bcftools, "view", "--header-only", str(vcf_path)],
        "Reading VCF header",
        logger=logger,
        error_type=error_type,
    )
    if not header.strip():
        raise error_type("bcftools returned an empty VCF header for %s" % vcf_path)
    return header


def declared_vcf_tags(header: str, category: str) -> set[str]:
    """Return INFO or FORMAT IDs declared by one VCF header."""
    normalized = str(category).upper()
    if normalized not in {"INFO", "FORMAT"}:
        raise ValueError("VCF tag category must be INFO or FORMAT")
    return set(
        re.findall(r"^##%s=<ID=([^,>]+)" % normalized, header, re.MULTILINE)
    )


def required_vcf_query_tags(expressions: Iterable[str]) -> list[str]:
    """Derive declared INFO/FORMAT tags required by bcftools query expressions."""
    required: list[str] = []
    for expression in expressions:
        text = str(expression)
        required.extend(
            "%s/%s" % match
            for match in re.findall(
                r"%(INFO|FORMAT)/([A-Za-z][A-Za-z0-9_.-]*)", text
            )
        )
        required.extend(
            "FORMAT/%s" % tag
            for tag in re.findall(
                r"\[%([A-Za-z][A-Za-z0-9_.-]*)\]", text
            )
        )
    return list(dict.fromkeys(required))


def validate_vcf_header_contract(
    *,
    header: str,
    genome_build_header: str,
    supported_genome_builds: Sequence[str],
    required_fields: Sequence[str],
) -> tuple[str, list[str]]:
    """Infer one exact declared build and validate active INFO/FORMAT fields."""
    template = str(genome_build_header).strip()
    if template.count("{build}") != 1 or "=" not in template:
        raise ValueError(
            "The configured VCF genome-build header must contain {build} exactly once"
        )
    builds = [str(build).strip() for build in supported_genome_builds]
    if not builds or any(not build for build in builds):
        raise ValueError("At least one non-empty supported genome build is required")
    rendered_headers = {build: template.format(build=build) for build in builds}
    if len(set(rendered_headers.values())) != len(rendered_headers):
        raise ValueError("The configured VCF genome-build header is ambiguous")
    metadata_prefix = template.split("=", 1)[0] + "="
    build_declarations = [
        line.strip()
        for line in header.splitlines()
        if line.strip().startswith(metadata_prefix)
    ]
    if len(build_declarations) != 1:
        raise ValueError(
            "Input VCF header must contain exactly one %s<build> declaration; "
            "found %d" % (metadata_prefix, len(build_declarations))
        )
    matches = [
        build
        for build, rendered in rendered_headers.items()
        if build_declarations[0] == rendered
    ]
    if len(matches) != 1:
        raise ValueError(
            "Input VCF declares unsupported genome build metadata %r. Expected "
            "exactly one of: %s"
            % (build_declarations[0], ", ".join(rendered_headers.values()))
        )
    declared = {
        "FORMAT": declared_vcf_tags(header, "FORMAT"),
        "INFO": declared_vcf_tags(header, "INFO"),
    }
    missing = []
    for field in dict.fromkeys(required_fields):
        try:
            category, tag = field.split("/", 1)
        except ValueError as exc:
            raise ValueError(
                "Required VCF fields must use INFO/TAG or FORMAT/TAG: %s" % field
            ) from exc
        if category not in declared:
            raise ValueError("Unsupported required VCF field category: %s" % category)
        if tag not in declared[category]:
            missing.append(field)
    if missing:
        raise ValueError(
            "Input VCF header is missing fields required by active operations: %s"
            % ", ".join(missing)
        )
    contigs = re.findall(r"^##contig=<ID=([^,>]+)", header, re.MULTILINE)
    return matches[0], contigs


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
    "declared_vcf_tags",
    "extract_vcf_table",
    "read_vcf_header",
    "required_vcf_query_tags",
    "select_vcf_sample",
    "validate_vcf_header_contract",
]
