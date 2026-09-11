"""GCTA-COJO selection handoff and paper-style physical locus grouping."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any, Mapping

import polars as pl

from postgwas.core.io.reports import (
    write_delimited_report,
    write_text_lines_report,
)
from postgwas.modules.ld_clumping.common import (
    LDClumpingError,
    chromosome_sort_key,
    normalise_chromosome,
)


COJO_LOCUS_BASE_COLUMNS = (
    "COJO_locus",
    "Chr",
    "start",
    "end",
    "span_bp",
    "index_SNP",
    "index_bp",
    "selected_signals",
    "selected_signal_ids",
)

COJO_SIGNAL_GROUPING_COLUMNS = (
    "COJO_locus",
    "is_locus_index",
    "locus_start",
    "locus_end",
    "distance_to_previous_signal_bp",
)


def _one_column_identifiers(path: str | Path, delimiter_pattern: str) -> set[str]:
    pattern = re.compile(delimiter_pattern)
    identifiers: set[str] = set()
    with Path(path).expanduser().resolve().open(
        "r", encoding="utf-8", errors="replace",
    ) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = pattern.split(raw.strip())
            if len(fields) != 1 or not fields[0]:
                raise LDClumpingError(
                    "COJO exclusion list line %d must contain exactly one SNP "
                    "identifier: %s" % (line_number, path)
                )
            identifiers.add(fields[0])
    return identifiers


def prepare_cojo_exclusion_list(
    *,
    reference_prefix: str | Path,
    destination: str | Path,
    ld_module,
    cojo_module,
) -> dict[str, Any]:
    """Create the exact GCTA exclusion list, including configured MHC removal."""
    configured = set()
    if cojo_module.inputs.exclude_snps is not None:
        configured = _one_column_identifiers(
            cojo_module.inputs.exclude_snps,
            cojo_module.input_validation.table_delimiter_pattern,
        )

    mhc_identifiers: set[str] = set()
    if ld_module.remove_mhc:
        bim_suffixes = [
            suffix for suffix in cojo_module.reference.required_extensions
            if suffix.lower() == ".bim"
        ]
        if len(bim_suffixes) != 1:
            raise LDClumpingError(
                "modules.gcta_cojo.reference.required_extensions must contain "
                "exactly one .bim suffix."
            )
        bim_path = Path(str(reference_prefix) + bim_suffixes[0])
        roles = {
            name: index
            for index, name in enumerate(cojo_module.reference.bim_columns)
        }
        required_roles = {"chromosome", "variant_id", "position"}
        if not required_roles.issubset(roles):
            raise LDClumpingError(
                "The configured COJO BIM schema must define chromosome, "
                "variant_id, and position roles."
            )
        region = ld_module.mhc_regions[ld_module.genome_build]
        pattern = re.compile(cojo_module.reference.table_delimiter_pattern)
        with bim_path.open("r", encoding="utf-8", errors="replace") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip():
                    continue
                fields = pattern.split(raw.strip())
                if len(fields) != len(roles):
                    raise LDClumpingError(
                        "PLINK BIM line %d has %d fields; expected %d: %s"
                        % (line_number, len(fields), len(roles), bim_path)
                    )
                try:
                    position = int(fields[roles["position"]])
                except ValueError as exc:
                    raise LDClumpingError(
                        "PLINK BIM position is not an integer at line %d: %s"
                        % (line_number, bim_path)
                    ) from exc
                if (
                    normalise_chromosome(fields[roles["chromosome"]])
                    == normalise_chromosome(region.chromosome)
                    and region.start <= position <= region.end
                ):
                    mhc_identifiers.add(fields[roles["variant_id"]])
    combined = configured | mhc_identifiers
    output_path = None
    if combined:
        output_path = write_text_lines_report(sorted(combined), destination)
    return {
        "path": None if output_path is None else str(output_path),
        "configured_exclusions": len(configured),
        "mhc_exclusions": len(mhc_identifiers),
        "combined_exclusions": len(combined),
    }


def _resolved_pvalue_columns(cojo_module) -> tuple[str, str]:
    schema = cojo_module.results.schemas["slct"]
    suffix = (
        cojo_module.results.genomic_control_p_value_suffix
        if cojo_module.analysis.genomic_control
        else ""
    )
    return (
        schema.marginal_p_value_column + suffix,
        schema.model_p_value_column + suffix,
    )


def _locus_columns(marginal_p: str, joint_p: str) -> tuple[str, ...]:
    """Preserve GCTA's exact p-value labels, including genomic control."""
    return (
        *COJO_LOCUS_BASE_COLUMNS[:7],
        "index_%s" % marginal_p,
        "index_%s" % joint_p,
        *COJO_LOCUS_BASE_COLUMNS[7:],
    )


def _index_sort_key(row: Mapping[str, Any], cojo_module):
    schema = cojo_module.results.schemas["slct"]
    marginal_p, joint_p = _resolved_pvalue_columns(cojo_module)
    if row["_index_policy"] == "joint":
        p_column = joint_p
        effect_column = schema.model_effect_column
        se_column = schema.model_standard_error_column
    else:
        p_column = marginal_p
        effect_column = schema.marginal_effect_column
        se_column = schema.marginal_standard_error_column
    z_score = abs(float(row[effect_column]) / float(row[se_column]))
    return (
        float(row[p_column]),
        -z_score,
        int(row[schema.position_column]),
        str(row[schema.identifier_column]),
    )


def _report_table(
    records: list[dict[str, Any]],
    columns: list[str],
    *,
    output_file: str | Path,
    aliases: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    aliases = dict(aliases or {})
    display_columns = [aliases.get(column, column) for column in columns]
    rows = [
        tuple(record.get(aliases.get(column, column)) for column in columns)
        for record in records
    ]
    return {
        "columns": display_columns,
        "rows": rows,
        "total_rows": len(records),
        "output_file": str(output_file),
    }


def group_cojo_selected_signals(
    *,
    normalized_result: str | Path,
    selected_signals_path: str | Path,
    loci_path: str | Path,
    ld_module,
    cojo_module,
) -> dict[str, Any]:
    """Group all GCTA-selected signals by consecutive physical separation."""
    if cojo_module.mode != "slct":
        raise LDClumpingError(
            "LD-clumping method cojo-slct requires modules.gcta_cojo.mode=slct."
        )
    try:
        frame = pl.read_csv(
            normalized_result,
            separator=cojo_module.results.normalized_delimiter,
            null_values=cojo_module.results.null_values,
            infer_schema_length=cojo_module.results.infer_schema_length,
        )
    except (OSError, pl.exceptions.PolarsError) as exc:
        raise LDClumpingError(
            "Cannot read normalized GCTA-COJO result %s: %s"
            % (normalized_result, exc)
        ) from exc

    schema = cojo_module.results.schemas["slct"]
    required = list(schema.required_columns) + [
        cojo_module.results.estimation_status_column,
    ]
    suffix = (
        cojo_module.results.genomic_control_p_value_suffix
        if cojo_module.analysis.genomic_control
        else ""
    )
    required = [
        column + suffix
        if suffix and column in {
            schema.marginal_p_value_column,
            schema.model_p_value_column,
        }
        else column
        for column in required
    ]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise LDClumpingError(
            "Normalized GCTA-COJO result is missing required columns: %s"
            % ", ".join(missing)
        )

    chromosome_column = schema.chromosome_column
    position_column = schema.position_column
    identifier_column = schema.identifier_column
    records = frame.select(required).to_dicts()
    for record in records:
        record[chromosome_column] = normalise_chromosome(
            record[chromosome_column]
        )
        record[position_column] = int(record[position_column])
    records.sort(key=lambda row: (
        chromosome_sort_key(row[chromosome_column]),
        row[position_column],
        str(row[identifier_column]),
    ))

    if ld_module.remove_mhc:
        region = ld_module.mhc_regions[ld_module.genome_build]
        retained = [
            row for row in records
            if not (
                row[chromosome_column] == normalise_chromosome(region.chromosome)
                and region.start <= row[position_column] <= region.end
            )
        ]
        if len(retained) != len(records):
            raise LDClumpingError(
                "GCTA selected %d signal(s) inside the configured MHC despite "
                "the validated exclusion list; COJO locus reporting stopped."
                % (len(records) - len(retained))
            )

    grouped: list[list[dict[str, Any]]] = []
    previous = None
    for record in records:
        if (
            previous is None
            or record[chromosome_column] != previous[chromosome_column]
            or record[position_column] - previous[position_column]
            > ld_module.cojo.merge_distance_bp
        ):
            grouped.append([])
        grouped[-1].append(record)
        previous = record

    marginal_p, joint_p = _resolved_pvalue_columns(cojo_module)
    marginal_locus_p = "index_%s" % marginal_p
    joint_locus_p = "index_%s" % joint_p
    locus_columns = _locus_columns(marginal_p, joint_p)
    signal_records: list[dict[str, Any]] = []
    locus_records: list[dict[str, Any]] = []
    for locus_number, members in enumerate(grouped, 1):
        locus_id = locus_number
        start = min(row[position_column] for row in members)
        end = max(row[position_column] for row in members)
        ranked = []
        for row in members:
            ranked_row = dict(row)
            ranked_row["_index_policy"] = ld_module.cojo.index_pvalue
            ranked.append(ranked_row)
        index = min(ranked, key=lambda row: _index_sort_key(row, cojo_module))
        index_id = index[identifier_column]
        previous_position = None
        member_ids = []
        for member in members:
            member_ids.append(str(member[identifier_column]))
            signal_records.append({
                "COJO_locus": locus_id,
                "is_locus_index": member[identifier_column] == index_id,
                "locus_start": start,
                "locus_end": end,
                "distance_to_previous_signal_bp": (
                    None
                    if previous_position is None
                    else member[position_column] - previous_position
                ),
                **member,
            })
            previous_position = member[position_column]
        locus_records.append({
            "COJO_locus": locus_id,
            "Chr": members[0][chromosome_column],
            "start": start,
            "end": end,
            "span_bp": end - start,
            "index_SNP": index_id,
            "index_bp": index[position_column],
            marginal_locus_p: index[marginal_p],
            joint_locus_p: index[joint_p],
            "selected_signals": len(members),
            "selected_signal_ids": ",".join(member_ids),
        })

    signal_columns = list(COJO_SIGNAL_GROUPING_COLUMNS) + required
    write_delimited_report(
        signal_records,
        selected_signals_path,
        fieldnames=signal_columns,
        delimiter=ld_module.table.delimiter,
        null_value=ld_module.table.null_values[0],
    )
    write_delimited_report(
        locus_records,
        loci_path,
        fieldnames=locus_columns,
        delimiter=ld_module.table.delimiter,
        null_value=ld_module.table.null_values[0],
    )

    report_columns = ld_module.reporting.result_table_columns
    aliases = {
        schema.marginal_p_value_column: marginal_p,
        schema.model_p_value_column: joint_p,
    }
    return {
        "status": "completed",
        "selected_signals": len(signal_records),
        "genomic_loci": len(locus_records),
        "merge_distance_bp": ld_module.cojo.merge_distance_bp,
        "index_pvalue": ld_module.cojo.index_pvalue,
        "selected_signals_file": str(selected_signals_path),
        "loci_file": str(loci_path),
        "output_files": {
            "selected_signals": str(selected_signals_path),
            "loci": str(loci_path),
        },
        "_report_tables": {
            "cojo_loci": _report_table(
                locus_records,
                report_columns["cojo_loci"],
                output_file=loci_path,
                aliases={
                    "index_p": marginal_locus_p,
                    "index_pJ": joint_locus_p,
                },
            ),
            "cojo_selected_signals": _report_table(
                signal_records,
                report_columns["cojo_selected_signals"],
                output_file=selected_signals_path,
                aliases=aliases,
            ),
        },
    }


__all__ = [
    "COJO_LOCUS_BASE_COLUMNS",
    "COJO_SIGNAL_GROUPING_COLUMNS",
    "group_cojo_selected_signals",
    "prepare_cojo_exclusion_list",
]
