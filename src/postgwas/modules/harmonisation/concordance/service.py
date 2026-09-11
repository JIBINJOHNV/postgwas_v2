"""I/O, logging and report orchestration for concordance validation."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import re
import shutil
import tempfile
import traceback
from typing import Any, Mapping

import polars as pl
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

from postgwas.core.io.delimiters import open_binary, open_text, resolve_delimiter
from postgwas.core.paths import configured_output_matches, configured_output_path
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.core.processes import run_checked_command
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.vcf import (
    extract_vcf_table,
    read_vcf_header,
    validate_postgwas_vcf_provenance,
)
from postgwas.modules.harmonisation.resource_paths import (
    external_resource_template_fields,
    resolve_resource_file,
)
from postgwas.modules.harmonisation.study_properties import (
    finalise_eaf_decision_from_chromosomes,
)
from postgwas.modules.harmonisation.strand import STRAND_ACTION_COLUMN

from .analysis import (
    INPUT_SOURCE_ROW_COLUMN,
    PARTITION_COLUMN,
    VCF_SOURCE_ROW_COLUMN,
    ConcordanceAnalysis,
    combine_concordance_partitions,
    compare_input_to_vcf,
    input_chromosome_expression,
    reference_chromosome_expression,
    vcf_chromosome_expression,
)
from .errors import ConcordanceValidationError


def _validation_logger(
    path: Path,
    dataset_id: str,
    policies,
    file_level: str,
    screen_level: str,
) -> PipelineLogger:
    """Create the canonical structured logger at the validation report path."""
    return PipelineLogger(
        sample_id=dataset_id,
        scope="validation",
        log_dir=str(path.parent),
        policies=policies,
        level=file_level,
        screen_level=screen_level,
        log_path=str(path),
    )


def _dataset_harmonisation_directory(
    output_root: str | Path,
    dataset_id: str,
    output_layout,
) -> Path:
    """Resolve the canonical result root shared with harmonisation."""
    return configured_output_path(
        output_root,
        output_layout["dataset_directory"],
        dataset_id=dataset_id,
    )


def write_validation_preflight_failure(
    output_root: str | Path,
    dataset_id: str,
    message: Any,
    output_layout,
    *,
    file_level: str,
    screen_level: str,
) -> str:
    """Write an actionable log even when failure happens before data loading."""
    base = _dataset_harmonisation_directory(
        output_root, dataset_id, output_layout,
    )
    path = configured_output_path(
        base, output_layout["concordance_log"], dataset_id=dataset_id,
    )
    log = _validation_logger(
        path, dataset_id, None, file_level, screen_level,
    )
    try:
        log.error(message)
        log.info("Concordance validation finished during preflight; log: %s" % path)
    finally:
        log.close()
    return str(path)


def _read_run_manifest(path: str | Path | None) -> tuple[Path | None, dict[str, Any]]:
    if path is None:
        return None, {}
    candidate = Path(path).expanduser().resolve()
    if not candidate.is_file():
        return candidate, {}
    try:
        value = json.loads(candidate.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConcordanceValidationError(
            "Cannot read the harmonisation run manifest %s: %s" % (candidate, exc)
        ) from exc
    if not isinstance(value, dict):
        raise ConcordanceValidationError(
            "Harmonisation run manifest must contain a JSON object: %s" % candidate
        )
    return candidate, value


def _manifest_context(manifest: dict[str, Any]) -> tuple[str | None, dict[str, Any]]:
    dataset = manifest.get("dataset") if isinstance(manifest.get("dataset"), dict) else {}
    build_info = dataset.get("genome_build")
    if isinstance(build_info, dict):
        input_build = build_info.get("inferred_build")
    else:
        input_build = None
    decisions = dataset.get("study_decisions")
    decisions = decisions if isinstance(decisions, dict) else {}
    return input_build, finalise_eaf_decision_from_chromosomes(
        decisions,
        dataset.get("chromosome_summaries"),
        dataset.get("completed"),
    )


def _vcf_build(header: str, supported_builds: list[str]) -> str | None:
    # Command-history headers may mention both builds after liftover. Contig
    # assembly metadata describes the records themselves and therefore wins.
    configured = {build.lower(): build for build in supported_builds}
    aliases = {
        alias.lower(): build
        for build in supported_builds
        for alias in (build, build.removeprefix("GRCh"))
    }
    assemblies = {
        aliases[value.strip().lower()]
        for value in re.findall(
            r"(?im)^##contig=<[^>\n]*assembly=([^,>]+)", header,
        )
        if value.strip().lower() in aliases
    }
    if len(assemblies) == 1:
        return assemblies.pop()
    if len(assemblies) > 1:
        return None
    explicit = {
        canonical
        for lowered, canonical in configured.items()
        if re.search(r"\b%s\b" % re.escape(lowered), header, flags=re.I)
    }
    return explicit.pop() if len(explicit) == 1 else None


def _separator(path: str | Path, configured: str, policies) -> str:
    try:
        return resolve_delimiter(
            path,
            configured,
            candidates=list(policies.get("input.delimiter_candidates")),
            minimum_columns=int(policies.get("input.delimiter_min_columns")),
            maximum_columns=int(policies.get("input.delimiter_max_columns")),
            sample_lines=int(policies.get("input.delimiter_sample_rows")),
        ).value
    except ValueError as exc:
        raise ConcordanceValidationError(str(exc)) from exc


def _input_columns(row) -> list[str]:
    columns = [
        row.chromosome_column,
        row.position_column,
        row.chromosome_position_column,
        row.effect_allele_column,
        row.other_allele_column,
        row.effect_allele_frequency_column,
        row.effect_column,
        row.standard_error_column,
        row.z_score_column,
        row.p_value_column,
    ]
    return list(dict.fromkeys(column for column in columns if column))


def _temporary_path(directory: Path) -> Path:
    handle = tempfile.NamedTemporaryFile(dir=directory, delete=False)
    path = Path(handle.name)
    handle.close()
    path.unlink()
    return path


def _stage_csv(
    source: str | Path,
    destination: Path,
    *,
    separator: str,
    columns: list[str],
    partition_expression: pl.Expr,
    null_values: list[str],
    schema_inference_rows: int,
    row_index_name: str | None,
    comment_prefix: str | None,
    description: str,
    log: PipelineLogger,
    batch_rows: int | None = None,
    compression: str | None = None,
) -> int:
    """Stream selected columns to a partition-addressable Parquet table."""
    try:
        if row_index_name or batch_rows is not None or compression is not None:
            if batch_rows is None or compression is None:
                raise ConcordanceValidationError(
                    "Batched concordance staging requires the resolved "
                    "staging_batch_rows and staging_compression settings."
                )
            rows = _stage_csv_batches(
                source,
                destination,
                separator=separator,
                columns=columns,
                partition_expression=partition_expression,
                null_values=null_values,
                schema_inference_rows=schema_inference_rows,
                batch_rows=batch_rows,
                compression=compression,
                row_index_name=row_index_name,
                comment_prefix=comment_prefix,
            )
        else:
            lazy = _partitioned_csv_scan(
                source,
                separator=separator,
                columns=columns,
                partition_expression=partition_expression,
                null_values=null_values,
                schema_inference_rows=schema_inference_rows,
                row_index_name=None,
                comment_prefix=comment_prefix,
            )
            sink_options: dict[str, Any] = {"maintain_order": True}
            if compression is not None:
                sink_options["compression"] = (
                    "uncompressed" if compression == "none" else compression
                )
            lazy.sink_parquet(destination, **sink_options)
            rows = int(
                pl.scan_parquet(destination)
                .select(pl.len().alias("rows"))
                .collect()
                .item()
                or 0
            )
    except (OSError, pl.exceptions.PolarsError, pa.ArrowException) as exc:
        destination.unlink(missing_ok=True)
        raise ConcordanceValidationError(
            "Cannot stage %s %s: %s" % (description, source, exc)
        ) from exc
    log.info(
        "%s staged: %s rows and %s selected columns; chromosome partitions "
        "will be read one at a time."
        % (description.capitalize(), f"{rows:,}", len(columns))
    )
    return rows


def _stage_csv_batches(
    source: str | Path,
    destination: Path,
    *,
    separator: str,
    columns: list[str],
    partition_expression: pl.Expr,
    null_values: list[str],
    schema_inference_rows: int,
    batch_rows: int,
    compression: str,
    row_index_name: str | None,
    comment_prefix: str | None,
) -> int:
    """Stream plain or compressed source text with optional stable row IDs.

    Polars 0.20's batched CSV reader opens compressed paths as raw bytes.  The
    shared input stream performs the configured format-independent
    decompression, while Arrow keeps every projected source token as text.
    Downstream concordance then applies the same explicit numerical and
    coordinate conversions as harmonisation, including exact extreme-P text.
    """
    del schema_inference_rows  # Projected source tokens intentionally remain text.

    leading_rows = 0
    if comment_prefix:
        with open_text(source) as handle:
            for line in handle:
                if not line.strip() or line.startswith(comment_prefix):
                    leading_rows += 1
                    continue
                break

    def invalid_row_handler(invalid_row) -> str:
        if comment_prefix and invalid_row.text.startswith(comment_prefix):
            return "skip"
        return "error"

    writer = None
    rows = 0
    try:
        with open_binary(source) as source_handle:
            with pa.PythonFile(source_handle, mode="r") as stream:
                reader = pacsv.open_csv(
                    stream,
                    read_options=pacsv.ReadOptions(skip_rows=leading_rows),
                    parse_options=pacsv.ParseOptions(
                        delimiter=separator,
                        invalid_row_handler=invalid_row_handler,
                    ),
                    convert_options=pacsv.ConvertOptions(
                        include_columns=columns,
                        column_types={column: pa.string() for column in columns},
                        null_values=null_values,
                        strings_can_be_null=True,
                    ),
                )
                for record_batch in reader:
                    batch = pl.from_arrow(record_batch)
                    for offset in range(0, batch.height, int(batch_rows)):
                        piece = batch.slice(offset, int(batch_rows))
                        indexed = (
                            piece.with_row_index(row_index_name, offset=rows + 1)
                            if row_index_name else piece
                        )
                        selected = (
                            [row_index_name] if row_index_name else []
                        ) + columns + [PARTITION_COLUMN]
                        staged = (
                            indexed.with_columns(
                                partition_expression.alias(PARTITION_COLUMN)
                            )
                            .select(selected)
                        )
                        table = staged.to_arrow()
                        if writer is None:
                            writer = pq.ParquetWriter(
                                destination,
                                table.schema,
                                compression=(
                                    None if compression == "none" else compression
                                ),
                            )
                        writer.write_table(table)
                        rows += piece.height
        if writer is None:
            empty = pl.DataFrame(
                schema={column: pl.String for column in columns}
            )
            indexed = (
                empty.with_row_index(row_index_name, offset=1)
                if row_index_name else empty
            )
            selected = (
                [row_index_name] if row_index_name else []
            ) + columns + [PARTITION_COLUMN]
            empty = (
                indexed.with_columns(partition_expression.alias(PARTITION_COLUMN))
                .select(selected)
            )
            empty.write_parquet(
                destination,
                compression="uncompressed" if compression == "none" else compression,
            )
    finally:
        if writer is not None:
            writer.close()
    return rows


def _partitioned_csv_scan(
    source: str | Path,
    *,
    separator: str,
    columns: list[str],
    partition_expression: pl.Expr,
    null_values: list[str],
    schema_inference_rows: int,
    row_index_name: str | None,
    comment_prefix: str | None,
) -> pl.LazyFrame:
    """Scan selected CSV fields with the canonical chromosome partition."""
    selected = ([row_index_name] if row_index_name else []) + columns
    return (
        pl.scan_csv(
            source,
            separator=separator,
            comment_prefix=comment_prefix,
            null_values=null_values,
            infer_schema_length=int(schema_inference_rows),
            ignore_errors=False,
            low_memory=True,
            row_index_name=row_index_name,
            row_index_offset=1 if row_index_name else 0,
        )
        .select(selected)
        .with_columns(partition_expression.alias(PARTITION_COLUMN))
    )


def _stage_input(
    row,
    destination: Path,
    policies,
    log: PipelineLogger,
    *,
    batch_rows: int,
    compression: str,
) -> int:
    separator = _separator(row.input_file, row.delimiter, policies)
    columns = _input_columns(row)
    log.info(
        "Streaming input fields [%s] from %s with separator %r."
        % (", ".join(columns), row.input_file, separator)
    )
    return _stage_csv(
        row.input_file,
        destination,
        separator=separator,
        columns=columns,
        partition_expression=input_chromosome_expression(row, policies),
        null_values=list(policies.get("input.null_values")),
        schema_inference_rows=int(policies.get("input.schema_inference_rows")),
        batch_rows=batch_rows,
        compression=compression,
        row_index_name=INPUT_SOURCE_ROW_COLUMN,
        comment_prefix="##",
        description="input summary statistics",
        log=log,
    )


def _stage_strand_actions(
    base: Path,
    destination: Path,
    dataset_id: str,
    output_layout,
    vcf_config,
    policies,
    log: PipelineLogger,
    *,
    batch_rows: int,
    compression: str,
) -> Path | None:
    """Stage archived reference-oriented rows for per-variant concordance."""
    archive_directory = configured_output_path(
        base,
        output_layout["adapter_input_archive"],
        dataset_id=dataset_id,
        chromosome="*",
    )
    archive_pattern = Path(output_layout["adapter_input"]).name + "*"
    sources = sorted(
        path
        for path in configured_output_matches(
            archive_directory,
            archive_pattern,
            dataset_id=dataset_id,
            chromosome="*",
        )
        if path.is_file() and path.stat().st_size > 0
    )
    mapping_path = configured_output_path(
        base,
        output_layout["adapter_merged_mapping"],
        dataset_id=dataset_id,
        chromosome="*",
    )
    if not sources or not mapping_path.is_file():
        log.warning(
            "Archived per-chromosome GWAS-to-VCF inputs or their merged mapping "
            "are unavailable; palindromic concordance will fall back to the "
            "study-wide strand decision."
        )
        return None
    try:
        mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConcordanceValidationError(
            "Cannot read archived GWAS-to-VCF column mapping %s: %s"
            % (mapping_path, exc)
        ) from exc
    required_keys = ("chr_col", "pos_col", "ea_col", "oa_col")
    missing_keys = [key for key in required_keys if key not in mapping]
    if missing_keys:
        raise ConcordanceValidationError(
            "Archived GWAS-to-VCF column mapping is missing: %s"
            % ", ".join(missing_keys)
        )
    separator = vcf_config["table_delimiter"]
    lazy_tables = []
    staged_sources: list[Path] = []
    rows = 0
    try:
        for source in sources:
            try:
                with open_text(source) as handle:
                    header = next(csv.reader(handle, delimiter=separator))
            except (OSError, StopIteration, csv.Error) as exc:
                raise ConcordanceValidationError(
                    "Cannot read archived GWAS-to-VCF header %s: %s" % (source, exc)
                ) from exc
            try:
                names = {
                    key: header[int(mapping[key])]
                    for key in required_keys
                }
            except (IndexError, TypeError, ValueError) as exc:
                raise ConcordanceValidationError(
                    "Archived GWAS-to-VCF mapping indexes do not match %s: %s"
                    % (source, exc)
                ) from exc
            if STRAND_ACTION_COLUMN not in header:
                raise ConcordanceValidationError(
                    "Archived GWAS-to-VCF input lacks required audit column %s: %s"
                    % (STRAND_ACTION_COLUMN, source)
                )
            source_columns = list(dict.fromkeys([
                names["chr_col"],
                names["pos_col"],
                names["oa_col"],
                names["ea_col"],
                STRAND_ACTION_COLUMN,
            ]))
            staged_source = _temporary_path(destination.parent)
            staged_sources.append(staged_source)
            rows += _stage_csv(
                source,
                staged_source,
                separator=separator,
                columns=source_columns,
                partition_expression=reference_chromosome_expression(
                    names["chr_col"], policies,
                ),
                null_values=list(policies.get("input.null_values")),
                schema_inference_rows=int(
                    policies.get("input.schema_inference_rows")
                ),
                row_index_name=None,
                comment_prefix=None,
                description="archived GWAS-to-VCF strand-action table",
                log=log,
                batch_rows=batch_rows,
                compression=compression,
            )
            lazy_tables.append(
                pl.scan_parquet(staged_source).select(
                    pl.col(names["chr_col"]).alias("CHROM"),
                    pl.col(names["pos_col"]).alias("POS"),
                    pl.col(names["oa_col"]).alias("REF"),
                    pl.col(names["ea_col"]).alias("ALT"),
                    pl.col(STRAND_ACTION_COLUMN),
                    pl.col(PARTITION_COLUMN),
                )
            )
        pl.concat(lazy_tables, how="vertical_relaxed").sink_parquet(
            destination,
            compression="uncompressed" if compression == "none" else compression,
            maintain_order=True,
        )
        staged_rows = int(
            pl.scan_parquet(destination).select(pl.len()).collect().item() or 0
        )
        if staged_rows != rows:
            raise ConcordanceValidationError(
                "Archived GWAS-to-VCF strand-action staging changed the row "
                "count (%s input rows; %s staged rows)."
                % (f"{rows:,}", f"{staged_rows:,}")
            )
    finally:
        for staged_source in staged_sources:
            staged_source.unlink(missing_ok=True)
    log.info(
        "Archived per-row strand actions staged: %s rows from %s chromosome "
        "input file(s)." % (f"{rows:,}", f"{len(sources):,}")
    )
    return destination


def _stage_duplicate_report(
    base: Path,
    destination: Path,
    row,
    output_layout,
    delimiter: str,
    policies,
    log: PipelineLogger,
    *,
    batch_rows: int | None = None,
    compression: str | None = None,
) -> Path | None:
    path = configured_output_path(
        base, output_layout["duplicates"], dataset_id=row.dataset_id,
    )
    if not path.is_file():
        return None
    columns = _input_columns(row) + ["duplicate_input_row", "duplicate_action"]
    rows = _stage_csv(
        path,
        destination,
        separator=delimiter,
        columns=columns,
        partition_expression=input_chromosome_expression(row, policies),
        null_values=list(policies.get("input.null_values")),
        schema_inference_rows=int(policies.get("input.schema_inference_rows")),
        batch_rows=batch_rows,
        compression=compression,
        row_index_name=None,
        comment_prefix=None,
        description="harmonisation duplicate report",
        log=log,
    )
    log.info("Applying %s recorded duplicate decisions from %s." % (f"{rows:,}", path))
    return destination


def _stage_external_eaf(
    row,
    mapping,
    destination: Path,
    policies,
    log: PipelineLogger,
    reference_build: str,
    *,
    batch_rows: int | None = None,
    compression: str | None = None,
) -> Path | None:
    if row.effect_allele_frequency_column:
        return None
    if not row.external_eaf_file or not row.external_eaf_column or mapping is None:
        raise ConcordanceValidationError(
            "The sample sheet does not provide a complete internal or external EAF source."
        )
    columns = _external_eaf_columns(row, mapping)
    try:
        template_fields = external_resource_template_fields(row.external_eaf_file)
    except ValueError as exc:
        raise ConcordanceValidationError(
            "Invalid external EAF path template %r: %s"
            % (row.external_eaf_file, exc)
        ) from exc
    if "chromosome" in template_fields:
        log.info(
            "External effect-allele frequencies will be resolved and read once "
            "per concordance chromosome from template %s."
            % row.external_eaf_file
        )
        return None
    try:
        source = resolve_resource_file(
            row.external_eaf_file,
            reference_build,
            "all",
        )
    except (FileNotFoundError, ValueError) as exc:
        raise ConcordanceValidationError(
            "Cannot resolve the external EAF file for concordance: %s" % exc
        ) from exc
    separator = _separator(source, mapping.delimiter, policies)
    log.info(
        "Streaming external effect-allele frequencies [%s] from %s with separator %r."
        % (", ".join(columns), source, separator)
    )
    _stage_csv(
        source,
        destination,
        separator=separator,
        columns=columns,
        partition_expression=reference_chromosome_expression(
            mapping.chromosome, policies,
        ),
        null_values=list(policies.get("input.null_values")),
        schema_inference_rows=int(policies.get("input.schema_inference_rows")),
        batch_rows=batch_rows,
        compression=compression,
        row_index_name=None,
        comment_prefix="##",
        description="external effect-allele-frequency table",
        log=log,
    )
    return destination


def _external_eaf_columns(row, mapping) -> list[str]:
    return list(dict.fromkeys([
        mapping.chromosome,
        mapping.position,
        mapping.effect_allele,
        mapping.other_allele,
        row.external_eaf_column,
    ]))


def _read_external_eaf_partition(
    row,
    mapping,
    partition: str,
    reference_build: str,
    policies,
    log: PipelineLogger,
    *,
    workspace: Path,
    batch_rows: int,
    compression: str,
) -> pl.DataFrame:
    """Read one resolved external-EAF chromosome without staging the full panel."""
    try:
        source = resolve_resource_file(
            row.external_eaf_file,
            reference_build,
            partition,
        )
    except (FileNotFoundError, ValueError) as exc:
        raise ConcordanceValidationError(
            "Cannot resolve the external EAF file for concordance chromosome %s "
            "from template %s: %s"
            % (partition, row.external_eaf_file, exc)
        ) from exc
    columns = _external_eaf_columns(row, mapping)
    separator = _separator(source, mapping.delimiter, policies)
    log.info(
        "Concordance chromosome %s: reading external effect-allele frequencies "
        "[%s] from %s with separator %r."
        % (partition, ", ".join(columns), source, separator)
    )
    staged_path = _temporary_path(workspace)
    try:
        _stage_csv(
            source,
            staged_path,
            separator=separator,
            columns=columns,
            partition_expression=reference_chromosome_expression(
                mapping.chromosome, policies,
            ),
            null_values=list(policies.get("input.null_values")),
            schema_inference_rows=int(
                policies.get("input.schema_inference_rows")
            ),
            row_index_name=None,
            comment_prefix="##",
            description=(
                "external effect-allele-frequency chromosome %s" % partition
            ),
            log=log,
            batch_rows=batch_rows,
            compression=compression,
        )
        frame = _read_partition(staged_path, partition)
    except (OSError, pl.exceptions.PolarsError) as exc:
        raise ConcordanceValidationError(
            "Cannot read the external EAF file for concordance chromosome %s "
            "(%s): %s" % (partition, source, exc)
        ) from exc
    finally:
        staged_path.unlink(missing_ok=True)
    log.info(
        "Concordance chromosome %s: external EAF rows read=%s."
        % (partition, f"{frame.height:,}")
    )
    return frame


def _extract_vcf(
    vcf_path: Path,
    dataset_id: str,
    temporary_directory: Path,
    destination: Path,
    bcftools: str,
    vcf_config,
    schema_inference_rows: int,
    policies,
    log: PipelineLogger,
    *,
    batch_rows: int,
    compression: str,
    header: str | None = None,
) -> tuple[int, str]:
    if header is None:
        header = run_checked_command(
            [bcftools, "view", "-h", str(vcf_path)],
            "Reading VCF header",
            logger=log,
            error_type=ConcordanceValidationError,
        )
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=vcf_config["temporary_table_suffix"],
        dir=temporary_directory, delete=False, encoding="utf-8",
    )
    temporary = Path(handle.name)
    handle.close()
    log.info(
        "Extracting configured VCF concordance fields once: %s."
        % ", ".join(vcf_config["concordance_fields"])
    )
    try:
        extract_vcf_table(
            vcf_path,
            temporary,
            dataset_id,
            vcf_config["concordance_fields"],
            bcftools,
            delimiter=vcf_config["table_delimiter"],
            io_buffer_bytes=vcf_config["io_buffer_bytes"],
            logger=log,
            error_type=ConcordanceValidationError,
            purpose="Extracting VCF fields for concordance validation",
        )
        rows = _stage_csv(
            temporary,
            separator=vcf_config["table_delimiter"],
            destination=destination,
            columns=list(vcf_config["concordance_fields"]),
            partition_expression=vcf_chromosome_expression(policies),
            null_values=vcf_config["table_null_values"],
            schema_inference_rows=int(schema_inference_rows),
            batch_rows=batch_rows,
            compression=compression,
            row_index_name=VCF_SOURCE_ROW_COLUMN,
            comment_prefix=None,
            description="extracted VCF concordance table",
            log=log,
        )
        return rows, header
    except (OSError, pl.exceptions.PolarsError) as exc:
        raise ConcordanceValidationError("Cannot extract VCF fields: %s" % exc) from exc
    finally:
        if not handle.closed:
            handle.close()
        try:
            temporary.unlink()
        except OSError:
            pass


def _partition_values(input_path: Path, vcf_path: Path) -> list[str | None]:
    partitions = (
        pl.concat([
            pl.scan_parquet(input_path).select(PARTITION_COLUMN),
            pl.scan_parquet(vcf_path).select(PARTITION_COLUMN),
        ])
        .unique()
        .sort(PARTITION_COLUMN, nulls_last=True)
        .collect()
    )
    values = partitions.get_column(PARTITION_COLUMN).to_list()
    return values or [None]


def _read_partition(path: Path, partition: str | None) -> pl.DataFrame:
    predicate = (
        pl.col(PARTITION_COLUMN).is_null()
        if partition is None
        else pl.col(PARTITION_COLUMN) == partition
    )
    return (
        pl.scan_parquet(path)
        .filter(predicate)
        .drop(PARTITION_COLUMN)
        .collect()
    )


def _spool_frame(
    frame: pl.DataFrame,
    directory: Path,
    paths: list[Path],
) -> pl.DataFrame:
    empty = frame.head(0)
    if frame.height:
        path = _temporary_path(directory)
        frame.write_parquet(path)
        paths.append(path)
    return empty


def _spooled_lazy(paths: list[Path], empty: pl.DataFrame) -> pl.LazyFrame:
    return pl.scan_parquet(paths) if paths else empty.lazy()


def _compare_staged_partitions(
    *,
    input_path: Path,
    vcf_path: Path,
    duplicate_path: Path | None,
    external_eaf_path: Path | None,
    strand_action_path: Path | None = None,
    workspace: Path,
    row,
    effect_type: str,
    se_scale: str | None = None,
    p_value_type: str,
    eaf_is_maf: bool | None,
    strand_consensus: str | None = None,
    settings,
    policies,
    external_eaf_mapping,
    external_eaf_build: str,
    log: PipelineLogger,
) -> ConcordanceAnalysis:
    partitions = _partition_values(input_path, vcf_path)
    log.info(
        "Concordance comparison will process %s canonical chromosome partition(s) "
        "sequentially: %s."
        % (
            f"{len(partitions):,}",
            ", ".join("invalid/null" if value is None else str(value) for value in partitions),
        )
    )
    partition_summaries: list[dict[str, Any]] = []
    report_paths: dict[str, list[Path]] = {
        "matched": [],
        "position_matches": [],
        "input_only": [],
        "vcf_only": [],
        "vcf_duplicates": [],
    }
    empty_reports: dict[str, pl.DataFrame] = {}
    for partition in partitions:
        input_frame = _read_partition(input_path, partition)
        vcf_frame = _read_partition(vcf_path, partition)
        duplicate_frame = (
            _read_partition(duplicate_path, partition)
            if duplicate_path is not None else None
        )
        strand_action_frame = (
            _read_partition(strand_action_path, partition)
            if strand_action_path is not None else None
        )
        if external_eaf_path is not None:
            external_eaf_frame = _read_partition(external_eaf_path, partition)
        elif (
            not row.effect_allele_frequency_column
            and row.external_eaf_file
            and input_frame.height
            and partition is not None
        ):
            external_eaf_frame = _read_external_eaf_partition(
                row,
                external_eaf_mapping,
                partition,
                external_eaf_build,
                policies,
                log,
                workspace=workspace,
                batch_rows=int(settings.staging_batch_rows),
                compression=str(settings.staging_compression),
            )
        else:
            external_eaf_frame = None
        analysis = compare_input_to_vcf(
            input_frame,
            vcf_frame,
            row,
            effect_type=effect_type,
            se_scale=se_scale,
            p_value_type=p_value_type,
            eaf_is_maf=eaf_is_maf,
            strand_consensus=strand_consensus,
            settings=settings,
            policies=policies,
            external_eaf_frame=external_eaf_frame,
            external_eaf_mapping=external_eaf_mapping,
            duplicate_report_frame=duplicate_frame,
            strand_action_frame=strand_action_frame,
        )
        partition_summaries.append(analysis.summary)
        for name in report_paths:
            frame = getattr(analysis, name)
            empty_reports[name] = _spool_frame(
                frame, workspace, report_paths[name],
            )
        log.info(
            "Concordance partition %s: input=%s, VCF=%s, matched=%s."
            % (
                "invalid/null" if partition is None else partition,
                f"{analysis.summary['input_rows']:,}",
                f"{analysis.summary['vcf_records']:,}",
                f"{analysis.summary['matched_variants']:,}",
            )
        )
        del (
            analysis,
            input_frame,
            vcf_frame,
            duplicate_frame,
            external_eaf_frame,
            strand_action_frame,
        )

    return combine_concordance_partitions(
        partition_summaries,
        matched=_spooled_lazy(report_paths["matched"], empty_reports["matched"]),
        position_matches=_spooled_lazy(
            report_paths["position_matches"], empty_reports["position_matches"],
        ),
        input_only=_spooled_lazy(
            report_paths["input_only"], empty_reports["input_only"],
        ),
        vcf_only=_spooled_lazy(report_paths["vcf_only"], empty_reports["vcf_only"]),
        vcf_duplicates=_spooled_lazy(
            report_paths["vcf_duplicates"], empty_reports["vcf_duplicates"],
        ),
        eaf_is_maf=eaf_is_maf,
        settings=settings,
    )


def _write_frame(
    frame: pl.DataFrame | pl.LazyFrame,
    path: Path,
    delimiter: str,
    null_value: str,
) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    compressed = path.suffix.lower() == ".gz"
    temporary = path.with_name(
        ".%s.tmp%s" % (path.stem, path.suffix)
        if compressed else ".%s.tmp" % path.name
    )
    uncompressed = path.with_name(".%s.uncompressed.tmp" % path.name)
    table_path = uncompressed if compressed else temporary
    try:
        if isinstance(frame, pl.LazyFrame):
            frame.sink_csv(
                table_path,
                separator=delimiter,
                null_value=null_value,
                maintain_order=True,
            )
        else:
            frame.write_csv(
                table_path,
                separator=delimiter,
                null_value=null_value,
            )
        if compressed:
            with uncompressed.open("rb") as source, gzip.open(temporary, "wb") as target:
                shutil.copyfileobj(source, target)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
        uncompressed.unlink(missing_ok=True)
    return str(path)


def _summary_rows(analysis: ConcordanceAnalysis, context: dict[str, Any]) -> list[dict[str, Any]]:
    rows = [
        {"section": "run", "metric": key, "value": value}
        for key, value in context.items()
    ]
    for key, value in analysis.summary.items():
        if key == "match_types":
            for match_type, count in sorted(value.items()):
                rows.append({"section": "matching", "metric": match_type, "value": count})
        elif key == "variant_types":
            for variant_type, values in value.items():
                for metric, metric_value in values.items():
                    rows.append({
                        "section": "matching_%s" % variant_type,
                        "metric": metric,
                        "value": metric_value,
                    })
        elif key == "position_diagnostics":
            for metric, metric_value in value.items():
                rows.append({
                    "section": "position_diagnostics",
                    "metric": metric,
                    "value": metric_value,
                })
        else:
            rows.append({"section": "retention", "metric": key, "value": value})
    for metric, values in analysis.metric_summary.items():
        for key, value in values.items():
            rows.append({"section": metric, "metric": key, "value": value})
    for variant_type, metrics in analysis.metric_summary_by_variant_type.items():
        for metric, values in metrics.items():
            for key, value in values.items():
                rows.append({
                    "section": "values_%s_%s" % (variant_type, metric),
                    "metric": key,
                    "value": value,
                })
    for variant_type, metrics in (
        analysis.position_metric_summary_by_variant_type.items()
    ):
        for metric, values in metrics.items():
            for key, value in values.items():
                rows.append({
                    "section": "position_%s_%s" % (variant_type, metric),
                    "metric": key,
                    "value": value,
                })
    return rows


def _write_summary(
    rows: list[dict[str, Any]],
    path: Path,
    delimiter: str,
    null_value: str,
) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(".%s.tmp" % path.name)
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["section", "metric", "value"],
            delimiter=delimiter,
        )
        writer.writeheader()
        for row in rows:
            value = row.get("value")
            if isinstance(value, (list, tuple)):
                value = "; ".join(str(item) for item in value)
            writer.writerow({
                "section": row.get("section"),
                "metric": row.get("metric"),
                "value": null_value if value is None else value,
            })
    temporary.replace(path)
    return str(path)


def _failure_summary(
    path: Path,
    dataset_id: str,
    error: Exception,
    delimiter: str,
    null_value: str,
) -> None:
    _write_summary([
        {"section": "run", "metric": "dataset_id", "value": dataset_id},
        {"section": "run", "metric": "status", "value": "FAIL"},
        {"section": "run", "metric": "error_type", "value": type(error).__name__},
        {"section": "run", "metric": "error", "value": str(error)},
    ], path, delimiter, null_value)


def _update_run_manifest(path: Path | None, result: dict[str, Any], log: PipelineLogger) -> None:
    if path is None or not path.is_file():
        return
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
        document["concordance_validation"] = result
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(path)
    except (OSError, json.JSONDecodeError, TypeError) as exc:
        log.warning("Could not add validation results to run manifest %s: %s" % (path, exc))


def _screen_summary(dataset_id: str, analysis: ConcordanceAnalysis, reports: dict[str, str]) -> None:
    summary = analysis.summary
    snps = summary["variant_types"]["snps"]
    indels = summary["variant_types"]["indels"]
    other_variants = summary["variant_types"]["other_variants"]
    position = summary["position_diagnostics"]
    width = 36

    def fraction_text(numerator: int, denominator: int, unit: str = "") -> str:
        unit_text = " %s" % unit if unit else ""
        if denominator <= 0:
            return "%s / %s%s (not applicable)" % (
                f"{numerator:,}", f"{denominator:,}", unit_text,
            )
        return "%s / %s%s (%0.2f%%)" % (
            f"{numerator:,}", f"{denominator:,}", unit_text,
            100.0 * numerator / denominator,
        )

    def metric_text(metric):
        if metric.get("frequency_type") == "unresolved":
            return "not compared; final frequency type unresolved"
        if not metric["checked"]:
            return "0 / 0 (not applicable)"
        text = fraction_text(metric["concordant"], metric["checked"])
        if metric["unavailable_in_input"]:
            text += "; %s input values unavailable" % f"{metric['unavailable_in_input']:,}"
        return text

    variant_groups = (
        ("SNPs", "snps", snps),
        ("Indels", "indels", indels),
        ("Other variants", "other_variants", other_variants),
    )
    position_labels = {
        "snps": "SNP",
        "indels": "Indel",
        "other_variants": "Other",
    }
    common_labels = {
        "snps": "SNPs",
        "indels": "indels",
        "other_variants": "other variants",
    }
    statistic_labels = {
        "snps": "SNP",
        "indels": "Indel",
        "other_variants": "Other-variant",
    }
    metric_groups = (
        ("Effect", "effect", "analysis"),
        ("Standard error", "standard_error", "analysis"),
        ("Allele frequency", "allele_frequency", "genetic"),
        ("Z score", "z_score", "analysis"),
        ("P-value (-log10)", "p_value", "analysis"),
    )
    lines = [
        screen_line("analysis", "Harmonisation concordance validation", indent=4),
        screen_field("info", "Dataset", dataset_id, indent=6, label_width=24),
        screen_field(
            "success" if analysis.status == "PASS" else "warning" if analysis.status == "WARNING" else "error",
            "Status", analysis.status, indent=6, label_width=24,
        ),
        screen_field(
            "info", "Reason", "; ".join(summary["status_reasons"]),
            indent=6, label_width=24,
        ),
        "",
        screen_line("genetic", "1. Input composition", indent=6),
        screen_field(
            "count", "Total unique variants",
            f"{summary['input_unique_variants']:,}", indent=8, label_width=width,
        ),
    ]
    for label, _key, values in variant_groups:
        lines.append(screen_field(
            "count", label,
            fraction_text(
                values["input_unique_variants"], summary["input_unique_variants"],
            ), indent=10, label_width=width,
        ))
    lines.extend([
        screen_field(
            "warning" if summary["input_invalid_rows"] else "info",
            "Invalid variants",
            fraction_text(summary["input_invalid_rows"], summary["input_rows"]),
            indent=8, label_width=width,
        ),
        screen_field(
            "warning" if summary["input_duplicate_rows"] else "info",
            "Duplicate variants",
            fraction_text(summary["input_duplicate_rows"], summary["input_rows"]),
            indent=8, label_width=width,
        ),
        "",
        screen_line("genetic", "2. Final VCF composition", indent=6),
        screen_field(
            "count", "Total unique variants",
            f"{summary['vcf_unique_variants']:,}", indent=8, label_width=width,
        ),
    ])
    for label, _key, values in variant_groups:
        lines.append(screen_field(
            "count", label,
            fraction_text(
                values["vcf_unique_variants"], summary["vcf_unique_variants"],
            ), indent=10, label_width=width,
        ))
    lines.extend([
        screen_field(
            "error" if summary["vcf_invalid_rows"] else "info",
            "Invalid VCF records",
            fraction_text(summary["vcf_invalid_rows"], summary["vcf_records"]),
            indent=8, label_width=width,
        ),
        screen_field(
            "error" if summary["vcf_duplicate_records"] else "info",
            "Duplicate VCF records",
            fraction_text(summary["vcf_duplicate_records"], summary["vcf_records"]),
            indent=8, label_width=width,
        ),
        "",
        screen_line("genetic", "3. Allele-aware comparison", indent=6),
        screen_field(
            "info", "Match key",
            "chromosome, minimally represented position and unordered allele pair",
            indent=8, label_width=width,
        ),
        screen_line("genetic", "3.1 Common allele-aware variants", indent=8),
        screen_field(
            "success", "Total common variants",
            fraction_text(
                summary["matched_variants"], summary["variant_union"],
                "union variants",
            ), indent=10, label_width=width,
        ),
    ])
    for label, _key, values in variant_groups:
        lines.append(screen_field(
            "success" if values["exact_matched_variants"] else "info",
            "Common %s" % common_labels[_key],
            fraction_text(
                values["exact_matched_variants"], values["variant_union"],
                "union %s" % common_labels[_key],
            ), indent=12, label_width=width,
        ))

    match_types = summary["match_types"]
    palindromic = sum(
        int(count)
        for match_type, count in match_types.items()
        if str(match_type).startswith("palindromic_")
    )
    palindromic_compared = int(summary["palindromic_variants_compared"])
    palindromic_excluded = int(
        summary["palindromic_variants_excluded_from_orientation"]
    )
    if not palindromic:
        palindromic_message = fraction_text(0, summary["matched_variants"])
    elif palindromic_excluded:
        palindromic_message = (
            "%s; %s compared and %s excluded because a safe strand orientation "
            "was unavailable or disabled by policy; SE and p-value remain comparable"
            % (
                fraction_text(palindromic, summary["matched_variants"]),
                f"{palindromic_compared:,}",
                f"{palindromic_excluded:,}",
            )
        )
    elif summary["palindromic_comparison_basis"] == "listed_allele_order":
        palindromic_message = (
            fraction_text(palindromic, summary["matched_variants"])
            + "; included from listed allele order by explicit policy"
        )
    else:
        palindromic_message = (
            fraction_text(palindromic, summary["matched_variants"])
            + "; included in effect, EAF and Z using the study-wide strand consensus"
        )
    lines.extend([
        screen_line("genetic", "3.2 Orientation of common variants", indent=8),
        screen_field("info", "Direct", fraction_text(
            int(match_types.get("direct", 0)), summary["matched_variants"],
        ), indent=10, label_width=width),
        screen_field("info", "Allele-swapped", fraction_text(
            int(match_types.get("allele_swapped", 0)), summary["matched_variants"],
        ), indent=10, label_width=width),
        screen_field("info", "Strand-complement", fraction_text(
            int(match_types.get("strand_complement", 0)), summary["matched_variants"],
        ), indent=10, label_width=width),
        screen_field("info", "Complement-swapped", fraction_text(
            int(match_types.get("strand_complement_swapped", 0)),
            summary["matched_variants"],
        ), indent=10, label_width=width),
        screen_field(
            "warning" if palindromic_excluded else "info", "Palindromic",
            palindromic_message,
            indent=10, label_width=width,
        ),
        screen_line("analysis", "3.3 Statistical concordance", indent=8),
    ])
    for label, key, kind in metric_groups:
        metric = analysis.metric_summary[key]
        lines.append(screen_field(
            "warning" if metric["mismatches"] else kind,
            label, metric_text(metric), indent=10, label_width=width,
        ))
    for _variant_label, variant_key, _values in variant_groups:
        lines.append(screen_line(
            "analysis", "%s statistics" % statistic_labels[variant_key],
            indent=10,
        ))
        for label, key, kind in metric_groups:
            metric = analysis.metric_summary_by_variant_type[variant_key][key]
            lines.append(screen_field(
                "warning" if metric["mismatches"] else kind,
                label, metric_text(metric), indent=12, label_width=width,
            ))

    lines.append(screen_line(
        "loss", "3.4 Variants failing allele-aware matching", indent=8,
    ))
    for source, total_key, count_key in (
        ("Input unmatched", "input_unique_variants", "input_only_variants"),
        ("VCF unmatched", "vcf_unique_variants", "vcf_only_variants"),
    ):
        lines.append(screen_line("loss", source, indent=10))
        lines.append(screen_field(
            "loss" if summary[count_key] else "success", "Total",
            fraction_text(summary[count_key], summary[total_key]),
            indent=12, label_width=width,
        ))
        type_count_key = (
            "input_only_variants" if source.startswith("Input")
            else "vcf_only_variants"
        )
        type_total_key = (
            "input_unique_variants" if source.startswith("Input")
            else "vcf_unique_variants"
        )
        for label, _key, values in variant_groups:
            lines.append(screen_field(
                "loss" if values[type_count_key] else "info", label,
                fraction_text(values[type_count_key], values[type_total_key]),
                indent=14, label_width=width,
            ))

    lines.extend([
        "",
        screen_line(
            "analysis",
            "4. Position diagnostics for allele-unmatched variants only",
            indent=6,
        ),
        screen_field(
            "info", "Method",
            "runs only after allele-aware matching fails; position matches do not establish variant identity",
            indent=8, label_width=width,
        ),
        screen_field(
            "info", "Variant-type strata",
            "position counts can overlap when a mixed or multiallelic position contains more than one variant type",
            indent=8, label_width=width,
        ),
        screen_field(
            "count", "Unmatched input positions",
            f"{position['input_unmatched_positions']:,}",
            indent=8, label_width=width,
        ),
    ])
    for label, _key, values in variant_groups:
        lines.append(screen_field(
            "count", "%s positions" % position_labels[_key],
            fraction_text(
                values["input_unmatched_positions"],
                position["input_unmatched_positions"],
            ), indent=10, label_width=width,
        ))
    lines.append(screen_field(
        "count", "Unmatched VCF positions",
        f"{position['vcf_unmatched_positions']:,}",
        indent=8, label_width=width,
    ))
    for label, _key, values in variant_groups:
        lines.append(screen_field(
            "count", "%s positions" % position_labels[_key],
            fraction_text(
                values["vcf_unmatched_positions"],
                position["vcf_unmatched_positions"],
            ), indent=10, label_width=width,
        ))
    lines.append(screen_field(
        "genetic", "Shared unmatched positions",
        fraction_text(
            position["shared_unmatched_positions"],
            position["unmatched_position_union"],
            "union positions",
        ), indent=8, label_width=width,
    ))
    for label, _key, values in variant_groups:
        lines.append(screen_field(
            "genetic", "%s positions" % position_labels[_key],
            fraction_text(
                values["shared_unmatched_positions"],
                values["unmatched_position_union"],
                "union %s positions" % position_labels[_key],
            ), indent=10, label_width=width,
        ))
    lines.append(screen_field(
        "loss", "Input-only unmatched positions",
        fraction_text(
            position["input_specific_positions"],
            position["input_unmatched_positions"],
        ), indent=8, label_width=width,
    ))
    for label, _key, values in variant_groups:
        lines.append(screen_field(
            "loss" if values["input_specific_positions"] else "info",
            "%s positions" % position_labels[_key],
            fraction_text(
                values["input_specific_positions"],
                values["input_unmatched_positions"],
            ), indent=10, label_width=width,
        ))
    lines.append(screen_field(
        "loss", "VCF-only unmatched positions",
        fraction_text(
            position["vcf_specific_positions"],
            position["vcf_unmatched_positions"],
        ), indent=8, label_width=width,
    ))
    for label, _key, values in variant_groups:
        lines.append(screen_field(
            "loss" if values["vcf_specific_positions"] else "info",
            "%s positions" % position_labels[_key],
            fraction_text(
                values["vcf_specific_positions"],
                values["vcf_unmatched_positions"],
            ), indent=10, label_width=width,
        ))
    lines.extend([
        screen_field(
            "analysis", "One-to-one diagnostic pairs",
            fraction_text(
                position["unambiguous_position_pairs"],
                position["shared_unmatched_positions"],
                "shared positions",
            ), indent=8, label_width=width,
        ),
        screen_field(
            "warning" if position["ambiguous_shared_positions"] else "info",
            "Ambiguous shared positions",
            fraction_text(
                position["ambiguous_shared_positions"],
                position["shared_unmatched_positions"],
            ), indent=8, label_width=width,
        ),
        screen_field(
            "warning" if position["variant_type_mismatch_positions"] else "info",
            "Variant-type mismatches",
            fraction_text(
                position["variant_type_mismatch_positions"],
                position["shared_unmatched_positions"],
            ), indent=10, label_width=width,
        ),
        screen_field(
            "warning" if position["multiallelic_ambiguous_positions"] else "info",
            "Multiallelic positions",
            fraction_text(
                position["multiallelic_ambiguous_positions"],
                position["shared_unmatched_positions"],
            ), indent=10, label_width=width,
        ),
    ])
    for _variant_label, variant_key, values in variant_groups:
        lines.append(screen_line(
            "analysis", "%s diagnostic statistics" % statistic_labels[variant_key],
            indent=8,
        ))
        lines.append(screen_field(
            "analysis", "Eligible pairs",
            fraction_text(
                values["unambiguous_position_pairs"],
                values["shared_unmatched_positions"],
            ), indent=10, label_width=width,
        ))
        for label, key, kind in metric_groups:
            diagnostic_label = {
                "effect": "Effect magnitude",
                "allele_frequency": "Folded allele frequency",
                "z_score": "Absolute Z magnitude",
            }.get(key, label)
            metric = analysis.position_metric_summary_by_variant_type[variant_key][key]
            lines.append(screen_field(
                "warning" if metric["mismatches"] else kind,
                diagnostic_label, metric_text(metric),
                indent=10, label_width=width,
            ))

    lines.extend(["", screen_line("info", "5. Concordance reports", indent=6)])
    for key, label in (
        ("summary", "Detailed summary"),
        ("input_only", "Input-unmatched variants"),
        ("vcf_only", "VCF-unmatched variants"),
        ("mismatches", "Statistical mismatches"),
        ("same_position_matches", "Same-position diagnostics"),
        ("vcf_duplicates", "Duplicate VCF records"),
        ("all_matches", "All allele-aware matches"),
    ):
        if key in reports:
            lines.append(screen_field(
                "info", label, Path(reports[key]).name,
                indent=8, label_width=width,
            ))
    print("\n" + "\n".join(lines) + "\n")


def run_concordance_validation(
    *,
    row,
    vcf_path: str | Path,
    output_root: str | Path,
    settings,
    policies,
    threads: int,
    bcftools: str,
    output_layout,
    vcf_config,
    file_log_level: str,
    screen_log_level: str,
    run_manifest: str | Path | None = None,
    external_eaf_mapping=None,
    show_screen: bool = True,
    provenance_headers: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Run validation with a complete log and status summary.

    Standalone callers supply the resolved PostGWAS provenance contract for
    their primary study VCF. Internal harmonisation retains its existing
    producer-owned input path; external reference tables are not subject to
    this study-origin check.
    """
    vcf = Path(vcf_path).expanduser().resolve()
    base = _dataset_harmonisation_directory(
        output_root, row.dataset_id, output_layout,
    )
    log_path = configured_output_path(
        base, output_layout["concordance_log"], dataset_id=row.dataset_id,
    )
    summary_path = configured_output_path(
        base, output_layout["concordance_summary"], dataset_id=row.dataset_id,
    )
    log = _validation_logger(
        log_path, row.dataset_id, policies, file_log_level, screen_log_level,
    )
    manifest_path = None
    workspace_context = None
    try:
        log.info("Concordance validation started for dataset %s." % row.dataset_id)
        log.info("Input summary statistics: %s" % row.input_file)
        log.info("VCF: %s" % vcf)
        log.info(
            "Concordance execution: %s configured thread(s); chromosome "
            "partitions are processed sequentially." % int(threads)
        )
        log.info("Resolved validation settings: %s" % json.dumps(settings.model_dump(mode="json"), sort_keys=True))
        log.info(
            "Shared input and p-value policies: %s"
            % json.dumps({
                "input.null_values": policies.get("input.null_values"),
                "input.schema_inference_rows": policies.get("input.schema_inference_rows"),
                "pvalue.clip_low": policies.get("pvalue.clip_low"),
                "pvalue.clip_high": policies.get("pvalue.clip_high"),
                "pvalue.se_tail": policies.get("pvalue.se_tail"),
            }, sort_keys=True)
        )
        if not vcf.is_file() or vcf.stat().st_size <= 0:
            raise ConcordanceValidationError("VCF does not exist or is empty: %s" % vcf)

        header = None
        if provenance_headers is not None:
            header = read_vcf_header(
                vcf, bcftools, logger=log, error_type=ConcordanceValidationError,
            )
            validate_postgwas_vcf_provenance(
                header, provenance_headers, vcf_path=vcf, logger=log,
                error_type=ConcordanceValidationError,
            )

        if run_manifest is None:
            candidate = configured_output_path(
                base, output_layout["run_manifest"], dataset_id=row.dataset_id,
            )
            run_manifest = candidate if candidate.is_file() else None
        manifest_path, manifest = _read_run_manifest(run_manifest)
        input_build, decisions = _manifest_context(manifest)
        effect_type = decisions.get("effect_type") or (
            row.effect_type if row.effect_type != "auto" else None
        )
        p_value_type = decisions.get("pvalue_type") or (
            row.p_value_type if row.p_value_type != "auto" else None
        )
        if effect_type not in ("beta", "odds_ratio"):
            raise ConcordanceValidationError(
                "Effect type is unresolved. Keep the harmonisation run manifest beside the VCF "
                "or set effect_type to beta/odds_ratio in the sample sheet."
            )
        if p_value_type not in ("raw", "neglog10"):
            raise ConcordanceValidationError(
                "P-value type is unresolved. Keep the harmonisation run manifest beside the VCF "
                "or set p_value_type in the sample sheet."
            )
        se_scale = decisions.get("se_scale")
        if effect_type == "odds_ratio" and row.standard_error_column:
            configured_se_scale = policies.get("effect.se_scale")
            if se_scale not in ("log_odds", "as_given"):
                se_scale = (
                    configured_se_scale
                    if configured_se_scale in ("log_odds", "as_given")
                    else None
                )
            if se_scale is None:
                raise ConcordanceValidationError(
                    "Odds-ratio concordance cannot compare the supplied standard "
                    "error because its study-level scale is unresolved. Keep the "
                    "harmonisation run manifest beside the VCF, or explicitly set "
                    "effect.se_scale to 'log_odds' or 'as_given'."
                )
        eaf_is_maf = decisions.get("eaf_is_maf")
        strand_consensus = decisions.get("strand")
        if eaf_is_maf is None:
            log.warning(
                "The final frequency type is unresolved after chromosome reference "
                "validation. Allele-frequency concordance will be skipped; variant, "
                "effect, and Z-score concordance will still run."
            )

        temporary_prefix = configured_output_path(
            base,
            output_layout["concordance_temporary"],
            dataset_id=row.dataset_id,
        )
        temporary_prefix.parent.mkdir(parents=True, exist_ok=True)
        workspace_context = tempfile.TemporaryDirectory(
            prefix=temporary_prefix.name,
            dir=temporary_prefix.parent,
        )
        workspace = Path(workspace_context.name)
        staging_batch_rows = int(settings.staging_batch_rows)
        staging_compression = str(settings.staging_compression)
        input_path = _temporary_path(workspace)
        duplicate_path = _temporary_path(workspace)
        external_eaf_path = _temporary_path(workspace)
        strand_action_path = _temporary_path(workspace)
        vcf_path_staged = _temporary_path(workspace)
        _stage_input(
            row,
            input_path,
            policies,
            log,
            batch_rows=staging_batch_rows,
            compression=staging_compression,
        )
        duplicate_path = _stage_duplicate_report(
            base,
            duplicate_path,
            row,
            output_layout,
            vcf_config["table_delimiter"],
            policies,
            log,
            batch_rows=staging_batch_rows,
            compression=staging_compression,
        )
        strand_action_path = _stage_strand_actions(
            base,
            strand_action_path,
            row.dataset_id,
            output_layout,
            vcf_config,
            policies,
            log,
            batch_rows=staging_batch_rows,
            compression=staging_compression,
        )
        _, header = _extract_vcf(
            vcf,
            row.dataset_id,
            workspace,
            vcf_path_staged,
            bcftools,
            vcf_config,
            policies.get("input.schema_inference_rows"),
            policies,
            log,
            batch_rows=staging_batch_rows,
            compression=staging_compression,
            header=header,
        )
        observed_build = _vcf_build(
            header, list(vcf_config["target_builds"]),
        )
        if observed_build is None:
            raise ConcordanceValidationError(
                "The VCF header does not identify one configured genome build (%s): %s"
                % (", ".join(vcf_config["target_builds"]), vcf)
            )
        if input_build and input_build != observed_build:
            raise ConcordanceValidationError(
                "Genome-build mismatch: input was harmonised as %s but the supplied VCF is %s."
                % (input_build, observed_build)
            )
        if input_build is None:
            log.warning(
                "No adjacent run manifest supplied an independently inferred input build; "
                "the VCF header build %s is recorded but cannot be cross-checked." % observed_build
            )

        external_eaf_path = _stage_external_eaf(
            row,
            external_eaf_mapping,
            external_eaf_path,
            policies,
            log,
            observed_build,
            batch_rows=staging_batch_rows,
            compression=staging_compression,
        )

        log.info(
            "Study decisions: build=%s, effect_type=%s, se_scale=%s, "
            "p_value_type=%s, frequency_type=%s, strand_consensus=%s."
            % (
                observed_build,
                effect_type,
                se_scale or "not_applicable",
                p_value_type,
                "MAF" if eaf_is_maf is True else "EAF" if eaf_is_maf is False else "unresolved",
                strand_consensus or "unavailable",
            )
        )
        analysis = _compare_staged_partitions(
            input_path=input_path,
            vcf_path=vcf_path_staged,
            duplicate_path=duplicate_path,
            external_eaf_path=external_eaf_path,
            strand_action_path=strand_action_path,
            workspace=workspace,
            row=row,
            effect_type=effect_type, se_scale=se_scale,
            p_value_type=p_value_type,
            eaf_is_maf=eaf_is_maf,
            strand_consensus=strand_consensus,
            settings=settings,
            policies=policies,
            external_eaf_mapping=external_eaf_mapping,
            external_eaf_build=observed_build,
            log=log,
        )
        reports = {
            "summary": _write_summary(_summary_rows(analysis, {
                "dataset_id": row.dataset_id,
                "input_file": str(row.input_file),
                "vcf": str(vcf),
                "genome_build": observed_build,
                "effect_type": effect_type,
                "se_scale": se_scale,
                "p_value_type": p_value_type,
                "frequency_type": (
                    "minor_allele_frequency"
                    if eaf_is_maf is True
                    else "effect_allele_frequency"
                    if eaf_is_maf is False
                    else "unresolved"
                ),
                "status": analysis.status,
            }), summary_path, vcf_config["table_delimiter"], vcf_config["table_null_output"]),
            "mismatches": _write_frame(
                analysis.mismatches,
                configured_output_path(
                    base, output_layout["concordance_mismatches"],
                    dataset_id=row.dataset_id,
                ), vcf_config["table_delimiter"], vcf_config["table_null_output"],
            ),
            "input_only": _write_frame(
                analysis.input_only,
                configured_output_path(
                    base, output_layout["concordance_input_only"],
                    dataset_id=row.dataset_id,
                ), vcf_config["table_delimiter"], vcf_config["table_null_output"],
            ),
            "vcf_only": _write_frame(
                analysis.vcf_only,
                configured_output_path(
                    base, output_layout["concordance_vcf_only"],
                    dataset_id=row.dataset_id,
                ), vcf_config["table_delimiter"], vcf_config["table_null_output"],
            ),
            "vcf_duplicates": _write_frame(
                analysis.vcf_duplicates,
                configured_output_path(
                    base, output_layout["concordance_vcf_duplicates"],
                    dataset_id=row.dataset_id,
                ), vcf_config["table_delimiter"], vcf_config["table_null_output"],
            ),
            "same_position_matches": _write_frame(
                analysis.position_matches,
                configured_output_path(
                    base, output_layout["concordance_position_matches"],
                    dataset_id=row.dataset_id,
                ), vcf_config["table_delimiter"], vcf_config["table_null_output"],
            ),
            "log": str(log_path),
        }
        if settings.write_all_matches:
            reports["all_matches"] = _write_frame(
                analysis.matched,
                configured_output_path(
                    base, output_layout["concordance_all_matches"],
                    dataset_id=row.dataset_id,
                ), vcf_config["table_delimiter"], vcf_config["table_null_output"],
            )
        result = {
            "status": analysis.status,
            "summary": analysis.summary,
            "metrics": analysis.metric_summary,
            "metrics_by_variant_type": analysis.metric_summary_by_variant_type,
            "position_metrics_by_variant_type": (
                analysis.position_metric_summary_by_variant_type
            ),
            "reports": reports,
        }
        log.info("Validation result: %s" % analysis.status)
        log.info("Validation summary: %s" % json.dumps(analysis.summary, sort_keys=True))
        for name, path in reports.items():
            log.info("Output %s: %s" % (name, path))
        _update_run_manifest(manifest_path, result, log)
        if show_screen:
            _screen_summary(row.dataset_id, analysis, reports)
        return result
    except Exception as exc:
        error = exc if isinstance(exc, ConcordanceValidationError) else ConcordanceValidationError(str(exc))
        log.error("Validation failed: %s: %s" % (type(exc).__name__, exc))
        log.error(traceback.format_exc())
        try:
            _failure_summary(
                summary_path,
                row.dataset_id,
                error,
                vcf_config["table_delimiter"],
                vcf_config["table_null_output"],
            )
        except Exception as report_error:
            log.error("Could not write failure summary: %s" % report_error)
        if error is exc:
            raise
        raise error from exc
    finally:
        if workspace_context is not None:
            workspace_context.cleanup()
        log.info("Concordance validation finished; log: %s" % log_path)
        log.close()
