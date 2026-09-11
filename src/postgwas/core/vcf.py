"""Reusable, checked access to single-sample VCF data."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable, Mapping, Sequence, Type

from postgwas.core.errors import FormattingError

from postgwas.core.input_validation import (
    current_validation_session,
    record_file_validation,
    validate_once,
)
from postgwas.core.paths import require_nonempty_file
from postgwas.core.processes import run_checked_command


VCF_TAG = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


class VcfQueryError(RuntimeError):
    """A VCF could not be queried without ambiguity."""


@dataclass(frozen=True)
class IndexedVcfValidation:
    """Reusable evidence for one unchanged, indexed, single-sample VCF."""

    vcf: Path
    bcftools: str
    header: str
    genome_build: str
    contigs: tuple[str, ...]
    dataset_id: str
    sample: str
    variant_count: int
    required_fields: tuple[str, ...]
    size_bytes: int
    modified_ns: int
    index_identities: tuple[tuple[Path, int, int], ...]
    reused: bool = False


def _vcf_file_identity(
    vcf_path: str | Path,
    *,
    error_type: Type[Exception],
) -> tuple[Path, int, int]:
    """Return the minimum identity needed to reject stale in-run evidence."""
    vcf = require_nonempty_file(
        vcf_path,
        "input VCF",
        error_type=error_type,
    )
    try:
        stat = vcf.stat()
    except OSError as exc:
        raise error_type("Cannot inspect input VCF %s: %s" % (vcf, exc)) from exc
    return vcf, int(stat.st_size), int(stat.st_mtime_ns)


def _vcf_index_identities(
    vcf: Path,
    *,
    error_type: Type[Exception],
) -> tuple[tuple[Path, int, int], ...]:
    """Describe any standard tabix/CSI indexes beside a compressed VCF."""
    identities = []
    for suffix in (".tbi", ".csi"):
        path = Path(str(vcf) + suffix)
        if not path.exists():
            continue
        try:
            stat = path.stat()
        except OSError as exc:
            raise error_type("Cannot inspect VCF index %s: %s" % (path, exc)) from exc
        if not path.is_file() or stat.st_size <= 0:
            message = "VCF index does not exist or is empty: %s" % path
            record_file_validation(path, "VCF index", status="failed", message=message)
            raise error_type(message)
        identities.append((path.resolve(), int(stat.st_size), int(stat.st_mtime_ns)))
    return tuple(identities)


def validate_indexed_vcf(
    vcf_path: str | Path,
    dataset_id: str,
    bcftools: str,
    *,
    genome_build_header: str,
    supported_genome_builds: Sequence[str],
    required_fields: Sequence[str] = (),
    provenance_headers: Mapping[str, str] | None = None,
    expected_genome_build: str | None = None,
    cached: IndexedVcfValidation | None = None,
    logger=None,
    error_type: Type[Exception] = VcfQueryError,
) -> IndexedVcfValidation:
    """Validate one indexed GWAS-VCF, reusing unchanged in-run evidence.

    File size and nanosecond modification time guard reuse within a single
    invocation. Durable checkpoint reuse continues to use the stronger shared
    checkpoint fingerprints; this transient evidence is never checkpointed.
    Study consumers supply ``provenance_headers``; reference-panel validators
    leave it unset. The origin gate reuses the header, not a record scan.
    """
    vcf, size_bytes, modified_ns = _vcf_file_identity(
        vcf_path,
        error_type=error_type,
    )
    index_identities = _vcf_index_identities(vcf, error_type=error_type)
    fields = tuple(dict.fromkeys(str(field) for field in required_fields))
    requested_dataset = str(dataset_id).strip()
    if not requested_dataset:
        raise error_type("A non-empty dataset identifier is required for VCF validation")
    resolved_bcftools = str(bcftools).strip()
    if not resolved_bcftools:
        raise error_type("A resolved bcftools executable is required for VCF validation")

    can_reuse = (
        isinstance(cached, IndexedVcfValidation)
        and cached.vcf == vcf
        and cached.size_bytes == size_bytes
        and cached.modified_ns == modified_ns
        and bool(index_identities)
        and cached.index_identities == index_identities
        and cached.dataset_id == requested_dataset
        and cached.bcftools == resolved_bcftools
    )
    # During a pipeline, all reuse goes through the session's stronger identity
    # checks. Explicit evidence must not bypass replacement/change detection.
    reuse_explicit = can_reuse and current_validation_session() is None
    if reuse_explicit:
        header = cached.header
    else:
        header = read_vcf_header(
            vcf,
            resolved_bcftools,
            logger=logger,
            error_type=error_type,
        )

    if provenance_headers is not None:
        validate_postgwas_vcf_provenance(
            header, provenance_headers, vcf_path=vcf,
            logger=logger, error_type=error_type,
        )

    try:
        genome_build, contigs = validate_vcf_header_contract(
            header=header,
            genome_build_header=genome_build_header,
            supported_genome_builds=supported_genome_builds,
            required_fields=fields,
        )
    except ValueError as exc:
        record_file_validation(vcf, "VCF header contract", status="failed", message=str(exc))
        raise error_type(str(exc)) from exc
    if (
        expected_genome_build is not None
        and genome_build != str(expected_genome_build)
    ):
        message = "Input VCF genome build %s does not match the required genome build %s." % (
            genome_build, expected_genome_build,
        )
        record_file_validation(vcf, "VCF genome build", status="failed", message=message)
        raise error_type(message)
    if reuse_explicit:
        sample = cached.sample
        variant_count = cached.variant_count
    else:
        variant_count = count_indexed_vcf_records(
            vcf,
            resolved_bcftools,
            logger=logger,
            error_type=error_type,
        )
        if variant_count == 0:
            message = "Input VCF contains no variant records: %s" % vcf
            record_file_validation(vcf, "Indexed VCF", status="failed", message=message)
            raise error_type(message)
        sample = select_vcf_sample(
            vcf,
            requested_dataset,
            resolved_bcftools,
            logger=logger,
            error_type=error_type,
        )
    evidence = IndexedVcfValidation(
        vcf=vcf,
        bcftools=resolved_bcftools,
        header=header,
        genome_build=genome_build,
        contigs=tuple(contigs),
        dataset_id=requested_dataset,
        sample=sample,
        variant_count=variant_count,
        required_fields=fields,
        size_bytes=size_bytes,
        modified_ns=modified_ns,
        index_identities=index_identities,
        reused=can_reuse,
    )
    record_file_validation(
        vcf,
        "Indexed VCF",
        checks=(
            "declared genome build",
            "required tag declarations",
            "exactly one sample column",
        ),
        metrics={
            "declared_genome_build": genome_build,
            "contigs": list(contigs),
            "variants_from_index": variant_count,
            "sample_count": 1,
            "sample": sample,
        },
    )
    if logger is not None and hasattr(logger, "record"):
        logger.record(
            "VALIDATE",
            "indexed_vcf_contract",
            status="PASSED",
            input_vcf=str(vcf),
            variants=variant_count,
            genome_build=genome_build,
            sample_count=1,
            sample=sample,
            contigs=evidence.contigs,
            required_fields=fields,
            validation_reused=can_reuse,
        )
    return evidence


def validate_postgwas_vcf_provenance(
    header: str,
    provenance_headers: Mapping[str, str],
    *,
    vcf_path: str | Path | None = None,
    logger=None,
    error_type: Type[Exception] = VcfQueryError,
) -> dict[str, str]:
    """Require declared PostGWAS origin, without certifying variant-level QC.

    Header names come from the resolved input contract. These declarations are
    not a digital signature, and the producing version need not equal the
    installed version. Keep this separate from consumer-specific field checks
    so a study-origin requirement never becomes a reference-panel requirement.
    """
    try:
        metadata = declared_vcf_metadata_values(header, provenance_headers.values())
        missing = [name for name in provenance_headers.values() if name not in metadata]
        empty = [name for name, value in metadata.items() if not value.strip()]
        invalid = [
            name for name, value in metadata.items()
            if value.strip().startswith("<") or any(c in value for c in "\r\n\t\0")
        ]
        issues = []
        for label, names in (("missing", missing), ("empty", empty), ("invalid text in", invalid)):
            if names:
                issues.append("%s PostGWAS provenance header(s): %s" % (label, ", ".join(names)))
        if issues:
            raise ValueError("; ".join(issues))
    except ValueError as exc:
        location = " %s" % vcf_path if vcf_path is not None else ""
        message = (
            "Input VCF%s is not a valid PostGWAS-harmonised GWAS-VCF: %s. "
            "This command accepts PostGWAS-harmonised study VCFs only. "
            "Re-run PostGWAS harmonisation on the original summary statistics; "
            "do not add provenance headers manually."
            % (location, exc)
        )
        if vcf_path is not None:
            record_file_validation(vcf_path, "GWAS-VCF provenance", status="failed", message=message)
        if logger is not None and hasattr(logger, "record"):
            logger.record("VALIDATE", "postgwas_vcf_provenance", status="FAILED", message=message)
        raise error_type(message) from exc
    evidence = {
        "postgwas_dataset_id": metadata[provenance_headers["dataset_id"]],
        "postgwas_version": metadata[provenance_headers["version"]],
        "postgwas_status": metadata[provenance_headers["status"]],
    }
    if vcf_path is not None:
        record_file_validation(
            vcf_path, "GWAS-VCF provenance", checks=("PostGWAS provenance",), metrics=evidence,
        )
    if logger is not None and hasattr(logger, "record"):
        logger.record(
            "VALIDATE", "postgwas_vcf_provenance", status="PASSED",
            input_vcf=str(vcf_path) if vcf_path is not None else None, **evidence,
        )
    return evidence


def validate_harmonised_vcf_header(
    header: str,
    config,
    *,
    vcf_path: str | Path | None = None,
    logger=None,
) -> dict[str, str | list[str]]:
    """Validate the formatter contract and declared PostGWAS origin."""
    contract = config.input_contract
    provenance = validate_postgwas_vcf_provenance(
        header, contract.provenance_headers.model_dump(), vcf_path=vcf_path,
        logger=logger, error_type=FormattingError,
    )
    required_fields = required_vcf_query_tags(config.vcf_fields.root.values())
    try:
        genome_build, contigs = validate_vcf_header_contract(
            header=header,
            genome_build_header="##" + contract.genome_build_metadata,
            supported_genome_builds=contract.supported_genome_builds,
            required_fields=required_fields,
        )
    except ValueError as exc:
        raise FormattingError(
            "Input is not a valid PostGWAS-harmonised GWAS-VCF: %s" % exc
        ) from exc

    return {
        "genome_build": genome_build,
        "contigs": contigs,
        "required_fields": required_fields,
        **provenance,
    }


def validate_harmonised_indexed_vcf(
    vcf_path: str | Path,
    dataset_id: str,
    bcftools: str,
    config,
    *,
    cached: IndexedVcfValidation | None = None,
) -> dict[str, object]:
    """Establish the reusable indexed and PostGWAS header contracts together."""
    indexed = validate_indexed_vcf(
        vcf_path,
        dataset_id,
        bcftools,
        genome_build_header="##" + config.input_contract.genome_build_metadata,
        supported_genome_builds=config.input_contract.supported_genome_builds,
        required_fields=required_vcf_query_tags(config.vcf_fields.root.values()),
        cached=cached,
        error_type=FormattingError,
    )
    harmonised = validate_harmonised_vcf_header(indexed.header, config, vcf_path=indexed.vcf)
    return {"indexed": indexed, "harmonised": harmonised}


def count_indexed_vcf_records(
    vcf_path: str | Path,
    bcftools: str,
    *,
    logger=None,
    error_type: Type[Exception] = VcfQueryError,
) -> int:
    """Return the validated record count stored in an indexed VCF."""
    # Index contents, not just the VCF, determine this answer. Never reuse a
    # count after either standard companion has been added, replaced or removed.
    indexes = _vcf_index_identities(Path(vcf_path).expanduser().resolve(), error_type=error_type)
    value = _read_vcf_validation_command(
        vcf_path, bcftools, ("index", "-n"), "Counting indexed VCF records",
        dependencies=tuple(item[0] for item in indexes),
        logger=logger, error_type=error_type,
    ).strip()
    try:
        count = int(value)
    except ValueError as exc:
        message = "bcftools returned an invalid VCF record count for %s: %r" % (vcf_path, value)
        record_file_validation(vcf_path, "VCF record count", status="failed", message=message)
        raise error_type(message) from exc
    if count < 0:
        message = "bcftools returned a negative VCF record count for %s: %s" % (vcf_path, count)
        record_file_validation(vcf_path, "VCF record count", status="failed", message=message)
        raise error_type(message)
    record_file_validation(
        vcf_path, "VCF", checks=("indexed record count",), metrics={"variants_from_index": count},
    )
    for index_path, _, _ in indexes:
        record_file_validation(index_path, "VCF index", checks=("index available for VCF record counting",))
    return count


def read_vcf_header(
    vcf_path: str | Path,
    bcftools: str,
    *,
    logger=None,
    error_type: Type[Exception] = VcfQueryError,
) -> str:
    """Read one complete VCF header with the configured bcftools executable."""
    header = _read_vcf_validation_command(
        vcf_path, bcftools, ("view", "--header-only"), "Reading VCF header",
        logger=logger, error_type=error_type,
    )
    if not header.strip():
        message = "bcftools returned an empty VCF header for %s" % vcf_path
        record_file_validation(vcf_path, "VCF header", status="failed", message=message)
        raise error_type(message)
    return header


def _read_vcf_validation_command(
    vcf_path, bcftools, arguments, purpose, *, dependencies=(), logger=None,
    error_type=VcfQueryError,
) -> str:
    """Share read-only VCF evidence, never arbitrary analysis subprocesses.

    Header, sample-list and index-count queries are format-level observations.
    Field requirements and scientific policies remain explicit in the callers.
    The cache is confined to the active run and tracks executable identity when
    an executable path is available, as well as each data/index dependency.
    """
    vcf = Path(vcf_path).expanduser().resolve()
    executable = Path(bcftools).expanduser()
    paths = (vcf, *dependencies)
    if executable.is_file():
        paths = (*paths, executable.resolve())

    def inspect():
        try:
            result = run_checked_command(
                [bcftools, *arguments, str(vcf_path)], purpose,
                logger=logger, error_type=error_type,
            )
        except Exception as exc:
            record_file_validation(vcf, purpose, checks=(purpose,), status="failed", message=str(exc))
            raise
        record_file_validation(
            vcf, purpose, checks=(purpose,),
            message="Read-only format evidence; record-level values are not certified by this query.",
        )
        return result

    return validate_once(
        paths,
        {"validator": "vcf_query", "executable": str(bcftools), "arguments": arguments},
        inspect, error_type=error_type,
    )


def declared_vcf_tags(header: str, category: str) -> set[str]:
    """Return INFO, FORMAT, or FILTER IDs declared by one VCF header."""
    normalized = str(category).upper()
    if normalized not in {"INFO", "FORMAT", "FILTER"}:
        raise ValueError("VCF tag category must be INFO, FORMAT, or FILTER")
    return set(
        re.findall(r"^##%s=<ID=([^,>]+)" % normalized, header, re.MULTILINE)
    )


def declared_vcf_tag_definitions(
    header: str,
    category: str,
) -> dict[str, dict[str, str]]:
    """Return Number and Type metadata for declared INFO or FORMAT fields.

    Duplicate declarations are invalid because callers cannot prove which
    metadata an external VCF tool will apply to the field.
    """
    normalized = str(category).upper()
    if normalized not in {"INFO", "FORMAT"}:
        raise ValueError("VCF tag category must be INFO or FORMAT")
    definitions: dict[str, dict[str, str]] = {}
    for payload in re.findall(
        r"^##%s=<(.*)>$" % normalized,
        header,
        re.MULTILINE,
    ):
        attributes: dict[str, str] = {}
        for name in ("ID", "Number", "Type"):
            match = re.search(r"(?:^|,)%s=([^,>]+)" % name, payload)
            if match:
                attributes[name.lower()] = match.group(1)
        tag = attributes.pop("id", "")
        if not tag:
            continue
        if tag in definitions:
            raise ValueError(
                "VCF header declares %s/%s more than once" % (normalized, tag)
            )
        definitions[tag] = attributes
    return definitions


def declared_vcf_metadata_values(
    header: str,
    names: Iterable[str],
) -> dict[str, str]:
    """Return uniquely declared values for selected simple metadata headers.

    PostGWAS provenance is written as JSON-compatible quoted text. Unquoted
    values remain supported for ordinary VCF metadata such as genome-build
    declarations, while structured ``<...>`` declarations are outside this
    helper's contract.
    """
    requested = list(dict.fromkeys(str(name) for name in names))
    invalid = [name for name in requested if not VCF_TAG.fullmatch(name)]
    if invalid:
        raise ValueError(
            "VCF metadata names are invalid: %s" % ", ".join(invalid)
        )
    requested_set = set(requested)
    values: dict[str, str] = {}
    for line in header.splitlines():
        match = re.match(r"^##([A-Za-z][A-Za-z0-9_.-]*)=(.*)$", line.strip())
        if match is None or match.group(1) not in requested_set:
            continue
        name, payload = match.groups()
        if name in values:
            raise ValueError(
                "VCF header declares metadata %s more than once" % name
            )
        if payload.startswith('"'):
            try:
                parsed = json.loads(payload)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    "VCF metadata %s has an invalid quoted value" % name
                ) from exc
            if not isinstance(parsed, str):
                raise ValueError(
                    "VCF metadata %s must contain a text value" % name
                )
            values[name] = parsed
        else:
            values[name] = payload
    return values


def annotation_vcf_fields(
    columns: Iterable[str],
    category: str,
) -> list[str]:
    """Extract qualified source fields from bcftools annotation columns."""
    normalized = str(category).upper()
    if normalized not in {"INFO", "FORMAT"}:
        raise ValueError("VCF tag category must be INFO or FORMAT")
    fields: list[str] = []
    for configured in columns:
        source = str(configured).split(":=", 1)[-1].lstrip("+.-")
        prefix = normalized + "/"
        if not source.startswith(prefix):
            continue
        tag = source[len(prefix):]
        if VCF_TAG.fullmatch(tag):
            fields.append(prefix + tag)
    return list(dict.fromkeys(fields))


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


def vcf_query_field_label(expression: str) -> str:
    """Render one configured bcftools query expression for user-facing provenance."""
    query = str(expression).strip()
    match = re.fullmatch(r"\[%([A-Za-z][A-Za-z0-9_.-]*)\]", query)
    if match:
        return "FORMAT/%s" % match.group(1)
    structural = re.fullmatch(r"%([A-Za-z][A-Za-z0-9_.-]*)", query)
    return structural.group(1) if structural else query


def required_vcf_field_presence(
    header: str,
    required_fields: Sequence[str],
) -> dict[str, bool]:
    """Return ordered declaration status for required INFO/FORMAT fields.

    This is the shared source for both contract enforcement and user-facing
    validation evidence. A declared field may still be missing on individual
    records; callers must report those record-level values separately.
    """
    declared = {
        "FORMAT": declared_vcf_tags(header, "FORMAT"),
        "INFO": declared_vcf_tags(header, "INFO"),
    }
    presence: dict[str, bool] = {}
    for field in dict.fromkeys(required_fields):
        try:
            category, tag = field.split("/", 1)
        except ValueError as exc:
            raise ValueError(
                "Required VCF fields must use INFO/TAG or FORMAT/TAG: %s" % field
            ) from exc
        if category not in declared:
            raise ValueError("Unsupported required VCF field category: %s" % category)
        presence[field] = tag in declared[category]
    return presence


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
    field_presence = required_vcf_field_presence(header, required_fields)
    missing = [field for field, present in field_presence.items() if not present]
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
    """Require one VCF sample and report a lone dataset-ID mismatch."""
    samples = [
        value
        for value in _read_vcf_validation_command(
            vcf_path, bcftools, ("query", "-l"), "Reading VCF sample names",
            logger=logger, error_type=error_type,
        ).splitlines()
        if value
    ]
    sample_count = len(samples)
    if sample_count != 1:
        message = (
            "Input GWAS-VCF must contain exactly one sample column; "
            "bcftools reported %d. Supply a single-sample PostGWAS GWAS-VCF."
            % sample_count
        )
        record_file_validation(
            vcf_path,
            "VCF sample",
            checks=("exactly one sample column",),
            metrics={"sample_count": sample_count},
            status="failed",
            message=message,
        )
        raise error_type(message)

    sample = samples[0]
    if sample != dataset_id:
        message = (
            "VCF sample %r differs from dataset ID %r; the only VCF sample "
            "will be used. Verify that this is the intended GWAS study."
            % (sample, dataset_id)
        )
        if logger is not None:
            logger.warning(message)
        record_file_validation(
            vcf_path,
            "VCF sample",
            checks=("exactly one sample column", "dataset identifier comparison"),
            metrics={
                "sample_count": sample_count,
                "sample": sample,
                "dataset_id": dataset_id,
            },
            status="warning",
            message=message,
        )
        return sample

    record_file_validation(
        vcf_path,
        "VCF sample",
        checks=("exactly one sample column", "dataset identifier match"),
        metrics={
            "sample_count": sample_count,
            "sample": sample,
            "dataset_id": dataset_id,
        },
    )
    return sample


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
    validated_sample: str | None = None,
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

    if validated_sample is None:
        sample = select_vcf_sample(
            vcf, dataset_id, bcftools, logger=logger, error_type=error_type,
        )
    else:
        sample = str(validated_sample).strip()
        if not sample:
            raise error_type("Validated VCF sample must not be empty")
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
    "IndexedVcfValidation",
    "VCF_TAG",
    "VcfQueryError",
    "annotation_vcf_fields",
    "count_indexed_vcf_records",
    "declared_vcf_metadata_values",
    "declared_vcf_tag_definitions",
    "declared_vcf_tags",
    "extract_vcf_table",
    "read_vcf_header",
    "required_vcf_field_presence",
    "required_vcf_query_tags",
    "select_vcf_sample",
    "validate_harmonised_vcf_header",
    "validate_harmonised_indexed_vcf",
    "validate_indexed_vcf",
    "validate_postgwas_vcf_provenance",
    "validate_vcf_header_contract",
    "vcf_query_field_label",
]
