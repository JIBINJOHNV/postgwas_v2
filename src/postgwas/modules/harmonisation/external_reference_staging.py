"""One-pass staging for user reference tables shared by chromosome workers.

Per-chromosome path templates already give each worker a bounded input and are
therefore left unchanged.  A single user file without ``{chromosome}``, on the
other hand, would otherwise be parsed independently by every spawned worker.
This module projects only the configured variant/value fields and writes one
atomic Parquet partition per observed chromosome before the process pool starts.
Parquet is an internal algorithm invariant here: it preserves the projected
schema without a second delimiter decision and supports configured compression
while each worker reads only its chromosome-sized partition.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import polars as pl
import pyarrow.parquet as pq

from postgwas.core.io.delimiters import open_text, resolve_delimiter
from postgwas.core.paths import configured_output_path
from postgwas.core.values import optional_text

from .resource_paths import external_resource_template_fields
from .shared.variant_columns import canonical_chromosome_expression


class ExternalReferenceStagingError(RuntimeError):
    """A shared user reference could not be staged safely."""


@dataclass(frozen=True)
class _ReferenceRole:
    resource_key: str
    label: str
    specification: str | None
    value_column: str | None
    mapping: Mapping[str, Any]
    output_layout_key: str


@dataclass(frozen=True)
class _StagingGroup:
    source: Path
    roles: tuple[_ReferenceRole, ...]
    columns: tuple[str, ...]
    chromosome_column: str
    configured_delimiter: str
    output_layout_key: str


def _selected_columns(
    role: _ReferenceRole,
) -> tuple[str, ...]:
    mapping = role.mapping
    missing = [
        key for key in ("chr", "pos", "a1", "a2", "delimiter")
        if optional_text(mapping.get(key)) is None
    ]
    if missing:
        raise ExternalReferenceStagingError(
            "%s mapping is missing: %s."
            % (role.label, ", ".join(missing))
        )
    if optional_text(role.value_column) is None:
        raise ExternalReferenceStagingError(
            "%s has no configured value column." % role.label
        )
    return tuple(dict.fromkeys([
        str(mapping[key]) for key in ("chr", "pos", "a1", "a2")
    ] + [str(role.value_column)]))


def _eligible_roles(
    resource_maps: Mapping[str, Mapping[str, Any]],
    roles: Sequence[_ReferenceRole],
    chromosomes: Sequence[str],
) -> list[tuple[_ReferenceRole, Path, tuple[str, ...]]]:
    if len(chromosomes) <= 1:
        return []
    eligible = []
    for role in roles:
        specification = optional_text(role.specification)
        if specification is None:
            continue
        if "chromosome" in external_resource_template_fields(specification):
            continue
        resolved = {
            optional_text(resource_maps[chromosome].get(role.resource_key))
            for chromosome in chromosomes
        }
        if None in resolved or len(resolved) != 1:
            raise ExternalReferenceStagingError(
                "%s does not resolve to one shared file for all observed "
                "chromosomes: %s."
                % (role.label, ", ".join(sorted(str(value) for value in resolved)))
            )
        source = Path(next(iter(resolved))).expanduser().resolve()
        eligible.append((role, source, _selected_columns(role)))
    return eligible


def _group_roles(
    eligible: Sequence[tuple[_ReferenceRole, Path, tuple[str, ...]]],
) -> tuple[_StagingGroup, ...]:
    grouped: dict[tuple[Any, ...], list[tuple[_ReferenceRole, tuple[str, ...]]]] = {}
    sources: dict[tuple[Any, ...], Path] = {}
    for role, source, columns in eligible:
        mapping = role.mapping
        signature = (
            str(source),
            str(mapping["chr"]),
            str(mapping["pos"]),
            str(mapping["a1"]),
            str(mapping["a2"]),
        )
        grouped.setdefault(signature, []).append((role, columns))
        sources[signature] = source

    results = []
    for signature, members in grouped.items():
        roles = tuple(member[0] for member in members)
        columns = tuple(dict.fromkeys(
            column for _role, selected in members for column in selected
        ))
        results.append(_StagingGroup(
            source=sources[signature],
            roles=roles,
            columns=columns,
            chromosome_column=signature[1],
            configured_delimiter=str(roles[0].mapping["delimiter"]),
            output_layout_key=roles[0].output_layout_key,
        ))
    return tuple(results)


def _reference_batches(
    source: Path,
    *,
    columns: Sequence[str],
    separator: str,
    null_values: Sequence[str],
    infer_schema_length: int,
    batch_rows: int,
    compressed_suffixes: Sequence[str],
) -> Iterable[pl.DataFrame]:
    """Yield a bounded projected scan using the established table semantics."""
    leading_rows = 0
    with open_text(source) as handle:
        header_line = None
        for line in handle:
            if not line.strip() or line.startswith("##"):
                leading_rows += 1
                continue
            header_line = line.rstrip("\r\n")
            break
    if header_line is None:
        raise ExternalReferenceStagingError(
            "Shared external reference has no readable header: %s" % source
        )
    header = (
        header_line.split()
        if separator == " "
        else next(csv.reader([header_line], delimiter=separator))
    )
    absent = [column for column in columns if column not in header]
    if absent:
        raise ExternalReferenceStagingError(
            "Shared external reference '%s' is missing configured column(s): %s. "
            "Available columns: %s."
            % (source, ", ".join(absent), ", ".join(header))
        )

    compressed = source.name.lower().endswith(tuple(compressed_suffixes))
    if separator == " " or compressed:
        try:
            import pandas as pd
        except ImportError as exc:
            raise ExternalReferenceStagingError(
                "Compressed or whitespace-separated external-reference staging "
                "requires pandas."
            ) from exc
        with open_text(source) as handle:
            reader = pd.read_csv(
                handle,
                sep=r"\s+" if separator == " " else separator,
                usecols=list(columns),
                dtype=str,
                na_values=list(null_values),
                keep_default_na=False,
                chunksize=int(batch_rows),
                skiprows=leading_rows,
            )
            for batch in reader:
                yield pl.from_pandas(batch).select(list(columns))
        return

    reader = pl.read_csv_batched(
        source,
        separator=separator,
        columns=list(columns),
        null_values=list(null_values),
        infer_schema_length=int(infer_schema_length),
        ignore_errors=False,
        low_memory=True,
        batch_size=int(batch_rows),
        skip_rows=leading_rows,
    )
    while True:
        batches = reader.next_batches(1)
        if not batches:
            return
        yield from batches


def projected_reference_batches(
    source: str | Path,
    *,
    columns: Sequence[str],
    separator: str | None,
    null_values: Sequence[str],
    infer_schema_length: int,
    batch_rows: int,
    compressed_suffixes: Sequence[str],
) -> Iterable[pl.DataFrame]:
    """Yield projected batches from a delimited source or staged Parquet file.

    Dataset-level consumers use the same bounded reader as reference staging,
    while staged partitions take the native Parquet batch path.  This keeps
    source parsing and memory bounds identical instead of adding a second eager
    whole-file implementation for dataset-wide validation.
    """
    path = Path(source).expanduser().resolve()
    if path.suffix.lower() == ".parquet":
        try:
            parquet = pq.ParquetFile(path)
            absent = [
                column for column in columns
                if column not in parquet.schema_arrow.names
            ]
            if absent:
                raise ExternalReferenceStagingError(
                    "Projected external-reference partition '%s' is missing "
                    "configured column(s): %s."
                    % (path, ", ".join(absent))
                )
            for batch in parquet.iter_batches(
                batch_size=int(batch_rows), columns=list(columns),
            ):
                yield pl.from_arrow(batch)
        except ExternalReferenceStagingError:
            raise
        except Exception as exc:
            raise ExternalReferenceStagingError(
                "Could not read projected external-reference partition '%s': %s"
                % (path, exc)
            ) from exc
        return

    if separator is None:
        raise ExternalReferenceStagingError(
            "A resolved delimiter is required for external reference %s." % path
        )
    yield from _reference_batches(
        path,
        columns=columns,
        separator=separator,
        null_values=null_values,
        infer_schema_length=infer_schema_length,
        batch_rows=batch_rows,
        compressed_suffixes=compressed_suffixes,
    )


def _partition_group(
    group: _StagingGroup,
    *,
    chromosomes: Sequence[str],
    dataset_id: str,
    output_directory: Path,
    output_layout: Mapping[str, str],
    settings: Mapping[str, Any],
    policies,
) -> tuple[dict[str, str], dict[str, Any]]:
    detected = resolve_delimiter(
        group.source,
        group.configured_delimiter,
        candidates=list(policies.get("input.delimiter_candidates")),
        minimum_columns=int(policies.get("input.delimiter_min_columns")),
        maximum_columns=int(policies.get("input.delimiter_max_columns")),
        sample_lines=int(policies.get("input.delimiter_sample_rows")),
    )
    if detected.value is None:
        raise ExternalReferenceStagingError(
            "Could not resolve the delimiter for shared external reference %s."
            % group.source
        )

    batch_rows = int(settings["batch_rows"])
    compression = str(settings["compression"])
    atomic_suffix = str(settings["atomic_output_suffix"])
    compressed_suffixes = tuple(settings["compressed_suffixes"])
    parquet_compression = None if compression == "none" else compression
    partition_column = "__postgwas_external_reference_chromosome"
    while partition_column in group.columns:
        partition_column += "_"

    destinations = {
        chromosome: configured_output_path(
            output_directory,
            output_layout[group.output_layout_key],
            error_type=ExternalReferenceStagingError,
            dataset_id=dataset_id,
            chromosome=chromosome,
        )
        for chromosome in chromosomes
    }
    if len(set(destinations.values())) != len(destinations):
        raise ExternalReferenceStagingError(
            "Configured external-reference partition paths are not unique by chromosome."
        )
    temporary = {
        chromosome: path.with_name(path.name + atomic_suffix)
        for chromosome, path in destinations.items()
    }
    for chromosome, path in destinations.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.is_symlink() or temporary[chromosome].is_symlink():
            raise ExternalReferenceStagingError(
                "Refusing to replace a symlinked external-reference partition: %s"
                % path
            )
        temporary[chromosome].unlink(missing_ok=True)

    writers: dict[str, pq.ParquetWriter] = {}
    rows_read = 0
    rows_by_chromosome = {chromosome: 0 for chromosome in chromosomes}
    read_started = False
    failure: Exception | None = None
    try:
        for batch in projected_reference_batches(
            group.source,
            columns=group.columns,
            separator=detected.value,
            null_values=list(policies.get("input.null_values")),
            infer_schema_length=int(policies.get("input.schema_inference_rows")),
            batch_rows=batch_rows,
            compressed_suffixes=compressed_suffixes,
        ):
            read_started = True
            selected = batch.select(list(group.columns))
            if not writers:
                arrow_schema = selected.head(0).to_arrow().schema
                for chromosome in chromosomes:
                    writers[chromosome] = pq.ParquetWriter(
                        temporary[chromosome],
                        arrow_schema,
                        compression=parquet_compression,
                    )
            rows_read += selected.height
            routed = (
                selected.with_columns(
                    canonical_chromosome_expression(
                        pl.col(group.chromosome_column), policies,
                    ).alias(partition_column)
                )
                .filter(pl.col(partition_column).is_in(list(chromosomes)))
            )
            for key, subset in routed.partition_by(
                [partition_column], as_dict=True, maintain_order=True,
            ).items():
                chromosome = str(key[0] if isinstance(key, tuple) else key)
                output = subset.select(list(group.columns))
                writers[chromosome].write_table(output.to_arrow())
                rows_by_chromosome[chromosome] += output.height
        if not read_started:
            raise ExternalReferenceStagingError(
                "Shared external reference contains no data rows: %s"
                % group.source
            )
    except Exception as exc:  # normalized below after every writer is closed
        failure = exc
    finally:
        for writer in writers.values():
            try:
                writer.close()
            except Exception as exc:
                if failure is None:
                    failure = exc

    if failure is not None:
        for path in temporary.values():
            path.unlink(missing_ok=True)
        if isinstance(failure, ExternalReferenceStagingError):
            raise failure
        raise ExternalReferenceStagingError(
            "Could not stage shared external reference '%s': %s"
            % (group.source, failure)
        ) from failure

    try:
        for chromosome in chromosomes:
            temporary[chromosome].replace(destinations[chromosome])
    except OSError as exc:
        for path in temporary.values():
            path.unlink(missing_ok=True)
        raise ExternalReferenceStagingError(
            "Could not publish external-reference partitions for '%s': %s"
            % (group.source, exc)
        ) from exc

    return (
        {chromosome: str(path) for chromosome, path in destinations.items()},
        {
            "source": str(group.source),
            "roles": [role.resource_key for role in group.roles],
            "selected_columns": list(group.columns),
            "rows_read": rows_read,
            "rows_by_chromosome": rows_by_chromosome,
            "rows_not_in_observed_chromosomes": (
                rows_read - sum(rows_by_chromosome.values())
            ),
            "delimiter_method": detected.method,
            "batch_rows": batch_rows,
            "compression": compression,
            "partitions": {
                chromosome: str(path)
                for chromosome, path in destinations.items()
            },
        },
    )


def stage_shared_external_reference_files(
    resource_maps: Mapping[str, Mapping[str, Any]],
    *,
    chromosomes: Sequence[str],
    dataset_id: str,
    output_directory: str | Path,
    output_layout: Mapping[str, str],
    user_eaf_specification: str | None,
    user_eaf_column: str | None,
    external_eaf_mapping: Mapping[str, Any],
    user_info_specification: str | None,
    user_info_column: str | None,
    external_info_mapping: Mapping[str, Any],
    settings: Mapping[str, Any],
    policies,
    logger=None,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    """Stage only whole-file user EAF/INFO inputs and return worker maps.

    A user path containing ``{chromosome}`` and every single-chromosome run are
    deliberately no-ops.  When EAF and INFO use the same file and identical
    structural mapping, their selected value fields are projected in one read
    and both worker resource keys point at the same chromosome partition.
    """
    ordered_chromosomes = tuple(str(chromosome) for chromosome in chromosomes)
    maps = {
        chromosome: dict(resource_maps[chromosome])
        for chromosome in ordered_chromosomes
    }
    roles = (
        _ReferenceRole(
            "user_eaf_file",
            "external allele-frequency table",
            user_eaf_specification,
            user_eaf_column,
            external_eaf_mapping,
            "external_eaf_partition",
        ),
        _ReferenceRole(
            "user_info_file",
            "external INFO table",
            user_info_specification,
            user_info_column,
            external_info_mapping,
            "external_info_partition",
        ),
    )
    groups = _group_roles(_eligible_roles(maps, roles, ordered_chromosomes))
    if not groups:
        return maps, {
            "status": "not_needed",
            "source_files_read": 0,
            "reason": (
                "no user EAF/INFO whole-genome file was shared by multiple "
                "observed chromosomes"
            ),
            "resources": [],
        }

    results = []
    for group in groups:
        partitions, result = _partition_group(
            group,
            chromosomes=ordered_chromosomes,
            dataset_id=dataset_id,
            output_directory=Path(output_directory).expanduser().resolve(),
            output_layout=output_layout,
            settings=settings,
            policies=policies,
        )
        for role in group.roles:
            for chromosome in ordered_chromosomes:
                maps[chromosome][role.resource_key] = partitions[chromosome]
        results.append(result)
        if logger is not None:
            logger.info(
                "Shared external reference staged once before chromosome fan-out: "
                "%s; roles=%s; rows=%s; selected_columns=%s; partitions=%s."
                % (
                    group.source,
                    ",".join(role.resource_key for role in group.roles),
                    result["rows_read"],
                    len(group.columns),
                    len(ordered_chromosomes),
                )
            )

    return maps, {
        "status": "staged",
        "source_files_read": len(groups),
        "resources": results,
    }


__all__ = [
    "ExternalReferenceStagingError",
    "projected_reference_batches",
    "stage_shared_external_reference_files",
]
