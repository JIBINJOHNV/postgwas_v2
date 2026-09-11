"""Export the validated chromosome table consumed by the gwas2vcf adapter."""

import json
from pathlib import Path
from typing import Dict

import polars as pl
from postgwas.core.paths import configured_output_path
from postgwas.core.values import optional_text

from .p_values import P_VALUE_EXACT_RAW_KEY
from .strand import EXPORTABLE_STRAND_ACTIONS, STRAND_ACTION_COLUMN


# This is the internal audit contract shared by the chromosome exporter and
# the dataset-level atomic merge. These names are protocol invariants rather
# than user-selectable study columns.
GWAS2VCF_SUMMARY_COLUMNS = (
    "chromosome",
    "key",
    "column_name",
    "dtype",
    "n_missing",
    "min",
    "max",
    "mean",
    "median",
    "std",
    "num_rows",
    "num_cols",
    "status",
)
GWAS2VCF_SUMMARY_CHROMOSOME_COLUMN = "chromosome"
GWAS2VCF_SUMMARY_KEY_COLUMN = "key"
GWAS2VCF_SUMMARY_MISSING_COLUMN = "n_missing"
GWAS2VCF_SUMMARY_ROW_COUNT_COLUMN = "num_rows"
GWAS2VCF_SUMMARY_COLUMN_COUNT_COLUMN = "num_cols"
GWAS2VCF_SUMMARY_STATUS_COLUMN = "status"
GWAS2VCF_SUMMARY_SUCCESS_STATUS = "success"


def summarise_gwas2vcf_columns(
    df: pl.DataFrame, exported_column_mappings: Dict, chromosome: str
) -> pl.DataFrame:
    """Summarize only the scientific mappings exported to GWAS-to-VCF."""
    numeric_types = {
        pl.Float32, pl.Float64, pl.Int8, pl.Int16, pl.Int32, pl.Int64,
        pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64,
    }
    mappings = []
    for key, value in exported_column_mappings.items():
        column = optional_text(value)
        if column is not None and column in df.columns:
            mappings.append((key, column))

    unique_columns = list(dict.fromkeys(column for _key, column in mappings))
    expressions = []
    for index, column in enumerate(unique_columns):
        prefix = "__column_%d" % index
        expressions.append(
            pl.col(column).null_count().alias(prefix + "_missing")
        )
        if df.schema[column] in numeric_types:
            expressions.extend([
                pl.col(column).min().alias(prefix + "_min"),
                pl.col(column).max().alias(prefix + "_max"),
                pl.col(column).mean().alias(prefix + "_mean"),
                pl.col(column).median().alias(prefix + "_median"),
                pl.col(column).std().alias(prefix + "_std"),
            ])
    values = df.select(expressions).to_dicts()[0] if expressions else {}
    column_index = {column: index for index, column in enumerate(unique_columns)}

    records = []
    for key, column in mappings:
        prefix = "__column_%d" % column_index[column]
        record = {
            "chromosome": chromosome,
            "key": key,
            "column_name": column,
            "dtype": str(df.schema[column]),
            "n_missing": values[prefix + "_missing"],
            "min": None, "max": None, "mean": None, "median": None, "std": None,
        }
        if df.schema[column] in numeric_types:
            record.update({
                statistic: values[prefix + "_" + statistic]
                for statistic in ("min", "max", "mean", "median", "std")
            })
        records.append(record)
    return pl.DataFrame(records)


def export_gwas2vcf_input(
    df: pl.DataFrame,
    sample_column_dict: Dict,
    output_dir: str,
    gwas_outputname: str,
    chromosome: str,
    genome_build: str,
    layout: Dict[str, str],
    input_config: Dict,
    logger=None,
) -> pl.DataFrame:
    """Write one adapter input TSV, its index-based mapping, and QC summary."""
    required_keys = tuple(input_config["required_column_keys"])
    optional_keys = tuple(input_config["optional_column_keys"])
    audit_columns = tuple(input_config["audit_columns"])
    missing = [
        key for key in required_keys
        if optional_text(sample_column_dict.get(key)) is None
    ]
    if missing:
        raise ValueError(
            "Cannot export chromosome %s; required mappings are missing: %s"
            % (chromosome, ", ".join(missing))
        )

    pairs = {
        key: sample_column_dict[key]
        for key in required_keys + optional_keys
        if optional_text(sample_column_dict.get(key)) is not None
    }
    absent = [column for column in pairs.values() if column not in df.columns]
    if absent:
        raise ValueError(
            "Cannot export chromosome %s; mapped columns are absent: %s"
            % (chromosome, ", ".join(absent))
        )
    absent_audit = [column for column in audit_columns if column not in df.columns]
    if absent_audit:
        raise ValueError(
            "Cannot export chromosome %s; configured audit columns are absent: %s"
            % (chromosome, ", ".join(absent_audit))
        )
    audit_null_counts = df.select([
        pl.col(column).null_count().alias(column) for column in audit_columns
    ]).row(0, named=True)
    incomplete_audit = {
        column: int(count)
        for column, count in audit_null_counts.items()
        if count
    }
    if incomplete_audit:
        raise ValueError(
            "Cannot export chromosome %s; configured audit columns contain missing "
            "values: %s"
            % (
                chromosome,
                ", ".join(
                    "%s=%s" % (column, count)
                    for column, count in incomplete_audit.items()
                ),
            )
        )
    if STRAND_ACTION_COLUMN in audit_columns:
        invalid_actions = (
            df.filter(
                ~pl.col(STRAND_ACTION_COLUMN)
                .is_in(list(EXPORTABLE_STRAND_ACTIONS))
                .fill_null(False)
            )
            .get_column(STRAND_ACTION_COLUMN)
            .drop_nulls()
            .unique()
            .sort()
            .to_list()
        )
        if invalid_actions:
            raise ValueError(
                "Cannot export chromosome %s; strand_action contains values that "
                "are neither reference-resolved nor explicitly retained for "
                "mandatory genome-FASTA validation: %s"
                % (chromosome, ", ".join(str(value) for value in invalid_actions))
            )
    mapped_audit = sorted(set(pairs.values()) & set(audit_columns))
    if mapped_audit:
        raise ValueError(
            "Cannot export chromosome %s; audit columns must not also be mapped "
            "GWAS-to-VCF columns: %s"
            % (chromosome, ", ".join(mapped_audit))
        )

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    export_columns = list(pairs.values()) + list(audit_columns)
    pvalue_column = optional_text(pairs.get("pval_col"))
    exact_raw_column = optional_text(
        sample_column_dict.get(P_VALUE_EXACT_RAW_KEY)
    )
    use_exact_raw = bool(
        exact_raw_column is not None
        and exact_raw_column in df.columns
        and df.get_column(exact_raw_column).null_count() < df.height
    )
    expressions = []
    for column in export_columns:
        if (
            column == pvalue_column
            and use_exact_raw
        ):
            # The vendored adapter accepts raw P as Decimal text and creates
            # FORMAT/LP itself. Prefer exact underflow text at this boundary;
            # changing the adapter contract to LP would be incorrect.
            expressions.append(
                pl.coalesce(
                    pl.col(exact_raw_column).cast(pl.String, strict=False),
                    pl.col(column).cast(pl.String, strict=False),
                ).alias(column)
            )
        else:
            expressions.append(pl.col(column))
    export_df = df.select(expressions)
    if pvalue_column is not None:
        missing_export_p = export_df.get_column(pvalue_column).null_count()
        if missing_export_p:
            raise ValueError(
                "Cannot export chromosome %s; %s p-values have neither a raw "
                "Float64 representation nor preserved exact raw text."
                % (chromosome, missing_export_p)
            )
    values = {"dataset_id": gwas_outputname, "chromosome": chromosome}
    tsv_path = configured_output_path(
        output, layout["adapter_input"], **values,
    )
    dict_path = configured_output_path(
        output, layout["adapter_mapping"], **values,
    )
    summary_path = configured_output_path(
        output, layout["adapter_summary"], **values,
    )

    export_df.write_csv(tsv_path, separator=input_config["delimiter"])
    renamed_keys = input_config["renamed_keys"]
    column_positions = {
        renamed_keys.get(key, key): export_df.columns.index(column)
        for key, column in pairs.items()
    }
    column_positions.update(
        delimiter=input_config["delimiter"],
        header=input_config["header"],
        build=genome_build,
    )
    with dict_path.open("w", encoding="utf-8") as handle:
        json.dump(column_positions, handle, indent=2)
        handle.write("\n")

    summary = summarise_gwas2vcf_columns(export_df, pairs, chromosome)
    summary = summary.with_columns(
        pl.lit(export_df.height).alias(GWAS2VCF_SUMMARY_ROW_COUNT_COLUMN),
        pl.lit(export_df.width).alias(GWAS2VCF_SUMMARY_COLUMN_COUNT_COLUMN),
        pl.lit(GWAS2VCF_SUMMARY_SUCCESS_STATUS).alias(
            GWAS2VCF_SUMMARY_STATUS_COLUMN
        ),
    ).select(list(GWAS2VCF_SUMMARY_COLUMNS))
    summary.write_csv(summary_path, separator=input_config["delimiter"])
    if logger is not None:
        logger.record(
            "OUTPUT", "adapter_input", rows=export_df.height,
            table=str(tsv_path), mapping=str(dict_path), summary=str(summary_path),
            audit_columns=",".join(audit_columns),
        )
    return summary
