"""I/O, logging and report orchestration for concordance validation."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import re
import tempfile
import traceback
from typing import Any

import polars as pl

from postgwas.core.io.delimiters import resolve_delimiter
from postgwas.core.paths import configured_output_path
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.core.processes import run_checked_command
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.vcf import extract_vcf_table
from postgwas.modules.harmonisation.study_properties import (
    finalise_eaf_decision_from_chromosomes,
)

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
    dataset_root = configured_output_path(
        output_root, output_layout["dataset_directory"], dataset_id=dataset_id,
    )
    base = configured_output_path(
        dataset_root, output_layout["analysis_directory"], dataset_id=dataset_id,
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
) -> int:
    """Stream selected columns to a partition-addressable Parquet table."""
    selected = ([row_index_name] if row_index_name else []) + columns
    try:
        lazy = pl.scan_csv(
            source,
            separator=separator,
            comment_prefix=comment_prefix,
            null_values=null_values,
            infer_schema_length=int(schema_inference_rows),
            ignore_errors=False,
            low_memory=True,
            row_index_name=row_index_name,
            row_index_offset=1 if row_index_name else 0,
        ).select(selected).with_columns(
            partition_expression.alias(PARTITION_COLUMN)
        )
        lazy.sink_parquet(destination, maintain_order=True)
        rows = int(
            pl.scan_parquet(destination)
            .select(pl.len().alias("rows"))
            .collect()
            .item()
            or 0
        )
    except (OSError, pl.exceptions.PolarsError) as exc:
        raise ConcordanceValidationError(
            "Cannot stage %s %s: %s" % (description, source, exc)
        ) from exc
    log.info(
        "%s staged: %s rows and %s selected columns; chromosome partitions "
        "will be read one at a time."
        % (description.capitalize(), f"{rows:,}", len(columns))
    )
    return rows


def _stage_input(
    row,
    destination: Path,
    policies,
    log: PipelineLogger,
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
        row_index_name=INPUT_SOURCE_ROW_COLUMN,
        comment_prefix="##",
        description="input summary statistics",
        log=log,
    )


def _stage_duplicate_report(
    base: Path,
    destination: Path,
    row,
    output_layout,
    delimiter: str,
    policies,
    log: PipelineLogger,
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
) -> Path | None:
    if row.effect_allele_frequency_column:
        return None
    if not row.external_eaf_file or not row.external_eaf_column or mapping is None:
        raise ConcordanceValidationError(
            "The sample sheet does not provide a complete internal or external EAF source."
        )
    columns = list(dict.fromkeys([
        mapping.chromosome,
        mapping.position,
        mapping.effect_allele,
        mapping.other_allele,
        row.external_eaf_column,
    ]))
    separator = _separator(row.external_eaf_file, "auto", policies)
    log.info(
        "Streaming external effect-allele frequencies [%s] from %s with separator %r."
        % (", ".join(columns), row.external_eaf_file, separator)
    )
    _stage_csv(
        row.external_eaf_file,
        destination,
        separator=separator,
        columns=columns,
        partition_expression=reference_chromosome_expression(mapping.chromosome),
        null_values=list(policies.get("input.null_values")),
        schema_inference_rows=int(policies.get("input.schema_inference_rows")),
        row_index_name=None,
        comment_prefix="##",
        description="external effect-allele-frequency table",
        log=log,
    )
    return destination


def _extract_vcf(
    vcf_path: Path,
    dataset_id: str,
    temporary_directory: Path,
    destination: Path,
    bcftools: str,
    vcf_config,
    schema_inference_rows: int,
    log: PipelineLogger,
) -> tuple[int, str]:
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
            partition_expression=vcf_chromosome_expression(),
            null_values=vcf_config["table_null_values"],
            schema_inference_rows=int(schema_inference_rows),
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
    workspace: Path,
    row,
    effect_type: str,
    p_value_type: str,
    eaf_is_maf: bool | None,
    settings,
    policies,
    external_eaf_mapping,
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
        external_eaf_frame = (
            _read_partition(external_eaf_path, partition)
            if external_eaf_path is not None else None
        )
        analysis = compare_input_to_vcf(
            input_frame,
            vcf_frame,
            row,
            effect_type=effect_type,
            p_value_type=p_value_type,
            eaf_is_maf=eaf_is_maf,
            settings=settings,
            policies=policies,
            external_eaf_frame=external_eaf_frame,
            external_eaf_mapping=external_eaf_mapping,
            duplicate_report_frame=duplicate_frame,
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
        del analysis, input_frame, vcf_frame, duplicate_frame, external_eaf_frame

    return combine_concordance_partitions(
        partition_summaries,
        matched=_spooled_lazy(report_paths["matched"], empty_reports["matched"]),
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
    compression = "gzip" if path.suffix == ".gz" else "uncompressed"
    temporary = path.with_name(
        ".%s.tmp%s" % (path.stem, path.suffix)
        if path.suffix == ".gz" else ".%s.tmp" % path.name
    )
    if isinstance(frame, pl.LazyFrame):
        frame.sink_csv(
            temporary,
            separator=delimiter,
            null_value=null_value,
            compression=compression,
            maintain_order=True,
        )
    else:
        frame.write_csv(
            temporary,
            separator=delimiter,
            null_value=null_value,
            compression=compression,
        )
    temporary.replace(path)
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
        else:
            rows.append({"section": "retention", "metric": key, "value": value})
    for metric, values in analysis.metric_summary.items():
        for key, value in values.items():
            rows.append({"section": metric, "metric": key, "value": value})
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
    effect = analysis.metric_summary["effect"]
    frequency = analysis.metric_summary["allele_frequency"]
    z_score = analysis.metric_summary["z_score"]
    snps = summary["variant_types"]["snps"]
    indels = summary["variant_types"]["indels"]
    other_variants = summary["variant_types"]["other_variants"]

    def concordance_text(metric):
        if metric.get("frequency_type") == "unresolved":
            return "not compared; final frequency type unresolved"
        if not metric["checked"]:
            return "not comparable; input value unavailable"
        text = "%s/%s concordant" % (
            f"{metric['concordant']:,}", f"{metric['checked']:,}",
        )
        if metric["unavailable_in_input"]:
            text += "; %s input values unavailable" % f"{metric['unavailable_in_input']:,}"
        return text

    def matching_text(values):
        input_count = values["input_unique_variants"]
        fraction = values["exact_match_fraction"]
        exact = "%s/%s exact" % (
            f"{values['exact_matched_variants']:,}", f"{input_count:,}",
        )
        if fraction is not None:
            exact += " (%0.2f%%)" % (100.0 * fraction)
        return "%s; input-only %s; VCF-only %s" % (
            exact,
            f"{values['input_only_variants']:,}",
            f"{values['vcf_only_variants']:,}",
        )

    print("\n" + "\n".join([
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
        screen_line("genetic", "Variant matching", indent=6),
        screen_field(
            "count", "Input unique variants",
            f"{summary['input_unique_variants']:,}", indent=8, label_width=28,
        ),
        screen_field(
            "warning"
            if summary["input_invalid_rows"] or summary["input_duplicate_rows"]
            else "info",
            "Invalid / duplicate input",
            "%s / %s" % (
                f"{summary['input_invalid_rows']:,}",
                f"{summary['input_duplicate_rows']:,}",
            ),
            indent=8, label_width=28,
        ),
        screen_field(
            "success", "Matched variants",
            f"{summary['matched_variants']:,} ({summary['retained_fraction']:.2%})",
            indent=8, label_width=28,
        ),
        screen_field("loss" if summary["input_only_variants"] else "info", "Unmatched input variants", f"{summary['input_only_variants']:,}", indent=8, label_width=28),
        screen_field("loss" if summary["vcf_only_variants"] else "info", "Unmatched VCF variants", f"{summary['vcf_only_variants']:,}", indent=8, label_width=28),
        screen_field("warning" if summary["vcf_duplicate_records"] else "info", "Duplicate VCF records", f"{summary['vcf_duplicate_records']:,}", indent=8, label_width=28),
        screen_field(
            "success" if snps["status"] == "PASS" else "warning" if snps["status"] in ("WARNING", "NOT_APPLICABLE") else "error",
            "SNP concordance",
            "%s — %s" % (snps["status"], matching_text(snps)),
            indent=8,
            label_width=28,
        ),
        screen_field(
            "success" if indels["status"] == "PASS" else "warning" if indels["status"] in ("WARNING", "NOT_APPLICABLE") else "error",
            "Indel concordance",
            "%s — %s" % (indels["status"], matching_text(indels)),
            indent=8,
            label_width=28,
        ),
        screen_field(
            "warning" if indels["maximum_possible_representation_pairs"] else "info",
            "Possible normalized indels",
            "%s maximum; not individually proven"
            % f"{indels['maximum_possible_representation_pairs']:,}",
            indent=8,
            label_width=28,
        ),
        screen_field(
            "success"
            if other_variants["status"] == "PASS"
            else "warning"
            if other_variants["status"] in ("WARNING", "NOT_APPLICABLE")
            else "error",
            "Other-variant concordance",
            "%s — %s" % (other_variants["status"], matching_text(other_variants)),
            indent=8,
            label_width=28,
        ),
        screen_field(
            "warning" if summary["palindromic_variants_excluded_from_orientation"] else "info",
            "Palindromic effect/Z",
            "%s excluded by policy"
            % f"{summary['palindromic_variants_excluded_from_orientation']:,}",
            indent=8, label_width=28,
        ),
        "",
        screen_line("analysis", "Value concordance", indent=6),
        screen_field("warning" if effect["unavailable_in_input"] else "analysis", "Effect", concordance_text(effect), indent=8, label_width=28),
        screen_field("warning" if frequency["unavailable_in_input"] else "genetic", "Allele frequency", concordance_text(frequency), indent=8, label_width=28),
        screen_field("warning" if z_score["unavailable_in_input"] else "analysis", "Z score", concordance_text(z_score), indent=8, label_width=28),
        screen_field("info", "Detailed summary", Path(reports["summary"]).name, indent=8, label_width=28),
    ]) + "\n")


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
) -> dict[str, Any]:
    """Run one validation and always leave a complete log and status summary."""
    vcf = Path(vcf_path).expanduser().resolve()
    dataset_root = configured_output_path(
        output_root,
        output_layout["dataset_directory"],
        dataset_id=row.dataset_id,
    )
    base = configured_output_path(
        dataset_root,
        output_layout["analysis_directory"],
        dataset_id=row.dataset_id,
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
        if p_value_type not in ("raw", "neglog10", "negln"):
            raise ConcordanceValidationError(
                "P-value type is unresolved. Keep the harmonisation run manifest beside the VCF "
                "or set p_value_type in the sample sheet."
            )
        eaf_is_maf = decisions.get("eaf_is_maf")
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
        input_path = _temporary_path(workspace)
        duplicate_path = _temporary_path(workspace)
        external_eaf_path = _temporary_path(workspace)
        vcf_path_staged = _temporary_path(workspace)
        _stage_input(row, input_path, policies, log)
        duplicate_path = _stage_duplicate_report(
            base,
            duplicate_path,
            row,
            output_layout,
            vcf_config["table_delimiter"],
            policies,
            log,
        )
        external_eaf_path = _stage_external_eaf(
            row,
            external_eaf_mapping,
            external_eaf_path,
            policies,
            log,
        )
        _, header = _extract_vcf(
            vcf,
            row.dataset_id,
            workspace,
            vcf_path_staged,
            bcftools,
            vcf_config,
            policies.get("input.schema_inference_rows"),
            log,
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

        log.info(
            "Study decisions: build=%s, effect_type=%s, p_value_type=%s, frequency_type=%s."
            % (
                observed_build,
                effect_type,
                p_value_type,
                "MAF" if eaf_is_maf is True else "EAF" if eaf_is_maf is False else "unresolved",
            )
        )
        analysis = _compare_staged_partitions(
            input_path=input_path,
            vcf_path=vcf_path_staged,
            duplicate_path=duplicate_path,
            external_eaf_path=external_eaf_path,
            workspace=workspace,
            row=row,
            effect_type=effect_type, p_value_type=p_value_type,
            eaf_is_maf=eaf_is_maf,
            settings=settings,
            policies=policies,
            external_eaf_mapping=external_eaf_mapping,
            log=log,
        )
        reports = {
            "summary": _write_summary(_summary_rows(analysis, {
                "dataset_id": row.dataset_id,
                "input_file": str(row.input_file),
                "vcf": str(vcf),
                "genome_build": observed_build,
                "effect_type": effect_type,
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
