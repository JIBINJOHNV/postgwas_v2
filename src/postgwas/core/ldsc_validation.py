"""Read-only LDSC manifests, configured reference inventories and input headers.

The two-field manifest contract follows the official LDSC cell-type tutorial:
https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses
Reference inventory passes establish availability, not numeric contents.
Summary-statistic header passes do not certify the remaining records.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import io
from pathlib import Path
from typing import Any

from postgwas.core.input_validation import record_file_validation, validate_once
from postgwas.core.io.delimiters import open_binary
from postgwas.core.paths import require_nonempty_file
from postgwas.core.validation_reporting import register_file_availability_bundle


@dataclass(frozen=True)
class LdscCelltypeManifestEntry:
    """One tested annotation and its required control annotation prefixes."""

    label: str
    prefixes: tuple[str, ...]


def resolve_ldscore_prefix(value: str, *, relative_to: Path) -> str:
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = relative_to / candidate
    return str(candidate.resolve())


def _inspect_ldcts(
    path: Path, *, prefix_separator: str, minimum_prefixes_per_cell_type: int,
    error_type,
) -> tuple[LdscCelltypeManifestEntry, ...]:
    """Parse the exact two-field syntax consumed by LDSC ``--ref-ld-chr-cts``."""
    entries: list[LdscCelltypeManifestEntry] = []
    labels: set[str] = set()
    tested_prefixes: set[str] = set()
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise error_type(
            "Cannot read LDSC .ldcts file %s: %s" % (path, exc)
        ) from exc
    for line_number, line in enumerate(lines, start=1):
        fields = line.split()
        if len(fields) != 2:
            raise error_type(
                "LDSC .ldcts line %d must contain exactly a cell-type name and "
                "a comma-delimited prefix list" % line_number
            )
        label, raw_prefixes = fields
        if label in labels:
            raise error_type(
                "Duplicate LDSC .ldcts cell-type name: %s" % label
            )
        labels.add(label)
        values = raw_prefixes.split(prefix_separator)
        if any(not value.strip() for value in values):
            raise error_type(
                "LDSC .ldcts line %d contains an empty LD-score prefix" % line_number
            )
        if len(values) < minimum_prefixes_per_cell_type:
            raise error_type(
                "LDSC .ldcts line %d requires at least %d prefixes: the tested "
                "annotation followed by its all-genes control"
                % (line_number, minimum_prefixes_per_cell_type)
            )
        prefixes = tuple(
            resolve_ldscore_prefix(value.strip(), relative_to=path.parent)
            for value in values
        )
        if len(prefixes) != len(set(prefixes)):
            raise error_type(
                "LDSC .ldcts line %d contains duplicate LD-score prefixes"
                % line_number
            )
        if prefixes[0] in tested_prefixes:
            raise error_type(
                "LDSC .ldcts tested annotation prefix is repeated: %s"
                % prefixes[0]
            )
        tested_prefixes.add(prefixes[0])
        entries.append(LdscCelltypeManifestEntry(label, prefixes))
    if not entries:
        raise error_type("LDSC .ldcts file must contain at least one cell type")
    return tuple(entries)


def read_ldcts_manifest(
    path: str | Path, *, prefix_separator: str,
    minimum_prefixes_per_cell_type: int, error_type=ValueError,
) -> tuple[LdscCelltypeManifestEntry, ...]:
    """Read every manifest row once for its exact native parsing contract."""
    source = require_nonempty_file(path, "LDSC .ldcts manifest", error_type=error_type)
    options = dict(
        prefix_separator=prefix_separator,
        minimum_prefixes_per_cell_type=minimum_prefixes_per_cell_type,
    )

    def inspect():
        entries = _inspect_ldcts(source, **options, error_type=error_type)
        record_file_validation(
            source, "LDSC cell-type manifest",
            checks=(
                "complete two-field manifest", "unique labels and tested prefixes",
                "configured tested/control prefix structure",
            ),
            metrics={
                "cell_types": len(entries),
                "prefix_entries": sum(len(row.prefixes) for row in entries),
            },
        )
        return entries

    return validate_once(
        (source,), {"validator": "ldsc_celltype_manifest", **options},
        inspect, error_type=error_type,
    )


def _chromosome_prefix(
    prefix: str, chromosome: int, placeholder: str, *, error_type,
) -> str:
    occurrences = prefix.count(placeholder)
    if occurrences > 1:
        raise error_type(
            "LDSC reference prefix contains the chromosome placeholder more "
            "than once: %s" % prefix
        )
    if occurrences:
        return prefix.replace(placeholder, str(chromosome))
    return "%s%d" % (prefix, chromosome)


def validate_ldscore_reference_inventory(
    baseline_ld_prefixes: tuple[str, ...], weights_ld_prefix: str,
    entries: tuple[LdscCelltypeManifestEntry, ...], *, chromosomes: int,
    chromosome_placeholder: str, ldscore_suffix: str, m_suffix: str,
    error_type=ValueError,
) -> tuple[dict[str, Any], ...]:
    """Check only the configured companion inventory; never certify table values."""
    requirements = [(prefix, "baseline", True) for prefix in baseline_ld_prefixes]
    requirements.append((weights_ld_prefix, "regression_weights", False))
    seen = set()
    for entry in entries:
        for index, prefix in enumerate(entry.prefixes):
            if prefix not in seen:
                seen.add(prefix)
                role = "cell_type_annotation" if index == 0 else "all_genes_control"
                requirements.append((prefix, role, True))
    expected = []
    for prefix, role, require_m_file in requirements:
        for chromosome in range(1, chromosomes + 1):
            expanded = _chromosome_prefix(
                prefix, chromosome, chromosome_placeholder, error_type=error_type,
            )
            companions = [(ldscore_suffix, "ld_score")]
            if require_m_file:
                companions.append((m_suffix, "regression_snp_count"))
            for suffix, resource in companions:
                expected.append((
                    Path(expanded + suffix),
                    dict(role=role, prefix=prefix, chromosome=chromosome, resource=resource),
                ))

    def inspect():
        records = []
        for path, metadata in expected:
            source = require_nonempty_file(path, "LDSC reference file", error_type=error_type)
            stat = source.stat()
            records.append({
                **metadata, "path": str(source), "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            })
            record_file_validation(
                source, "LDSC reference companion",
                checks=("regular file", "non-empty file"), metrics=metadata,
                message=(
                    "Availability only; LD-score/SNP-count contents are not "
                    "validated by this check."
                ),
            )
        return tuple(records)

    contract = {
        "validator": "ldscore_reference_inventory", "requirements": requirements,
        "chromosomes": chromosomes, "chromosome_placeholder": chromosome_placeholder,
        "ldscore_suffix": ldscore_suffix, "m_suffix": m_suffix,
    }
    inventory = deepcopy(validate_once(
        (path for path, _ in expected), contract, inspect, error_type=error_type,
    ))
    prefixes_by_role: dict[str, set[str]] = {}
    for item in inventory:
        prefixes_by_role.setdefault(item["role"], set()).add(item["prefix"])
    paths = tuple(dict.fromkeys(item["path"] for item in inventory))
    register_file_availability_bundle(
        paths,
        "LDSC cell-type reference bundle",
        (
            (
                "count",
                "chromosomes",
                sorted({str(item["chromosome"]) for item in inventory}, key=int),
            ),
            (
                "count",
                "ldsc_celltype_prefixes",
                len({item["prefix"] for item in inventory}),
            ),
            (
                "success",
                "ldsc_celltype_ld_files",
                sum(item["resource"] == "ld_score" for item in inventory),
            ),
            (
                "success",
                "ldsc_celltype_snp_count_files",
                sum(
                    item["resource"] == "regression_snp_count"
                    for item in inventory
                ),
            ),
            (
                "count",
                "ldsc_celltype_roles",
                {
                    role: len(prefixes)
                    for role, prefixes in prefixes_by_role.items()
                },
            ),
            ("count", "ldsc_unique_bundle_files", len(paths)),
            (
                "info",
                "validation_scope",
                "presence and non-empty files; numeric contents are not "
                "validated",
            ),
        ),
        covered_metric_keys=("role", "prefix", "chromosome", "resource"),
    )
    return inventory


def _inspect_munged_header(
    path: Path, *, required_columns: tuple[str, ...], error_type,
) -> dict[str, Any]:
    """Validate columns required by the LDSC regression before execution."""
    sumstats = path
    try:
        with open_binary(sumstats) as raw, io.TextIOWrapper(raw, encoding="utf-8") as handle:
            header_line = handle.readline()
            first_data_line = next((line for line in handle if line.strip()), "")
    except (OSError, UnicodeError) as exc:
        raise error_type(
            "Cannot read munged LDSC summary statistics %s: %s" % (sumstats, exc)
        ) from exc
    columns = header_line.split()
    required = required_columns
    missing = [column for column in required if column not in columns]
    if missing:
        raise error_type(
            "Munged LDSC summary statistics are missing required columns: %s"
            % ", ".join(missing)
        )
    if len(columns) != len(set(columns)):
        raise error_type("Munged LDSC summary-statistic columns must be unique")
    if not first_data_line:
        raise error_type("Munged LDSC summary statistics contain no variants")
    return {"columns": columns, "required_columns": list(required)}


def validate_ldsc_sumstats_header(
    path: str | Path, *, required_columns, error_type=ValueError,
) -> dict[str, Any]:
    """Validate the existing header/first-record gate, not every statistic value."""
    source = require_nonempty_file(
        path, "munged LDSC summary statistics", error_type=error_type,
    )
    required = tuple(required_columns)

    def inspect():
        result = _inspect_munged_header(
            source, required_columns=required, error_type=error_type,
        )
        record_file_validation(
            source, "Munged LDSC summary-statistic header",
            checks=("unique required header columns", "at least one nonempty data record"),
            metrics=result,
            message="Header/first-record check only; remaining records are not scanned.",
        )
        return result

    return deepcopy(validate_once(
        (source,), ("ldsc_munged_header", required), inspect, error_type=error_type,
    ))
