"""Export the validated chromosome table consumed by the gwas2vcf adapter."""

import json
from pathlib import Path
from typing import Dict

import polars as pl
from postgwas.core.paths import configured_output_path
from postgwas.core.values import optional_text


def summarise_gwas2vcf_columns(
    df: pl.DataFrame, sample_column_dict: Dict, chromosome: str
) -> pl.DataFrame:
    """Summarize all exported columns with one Polars aggregation."""
    numeric_types = {
        pl.Float32, pl.Float64, pl.Int8, pl.Int16, pl.Int32, pl.Int64,
        pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64,
    }
    mappings = []
    for key, value in sample_column_dict.items():
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

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    export_df = df.select(list(pairs.values()))
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

    summary = summarise_gwas2vcf_columns(export_df, sample_column_dict, chromosome)
    summary = summary.with_columns(
        pl.lit(str(tsv_path)).alias("tsv_path"),
        pl.lit(str(dict_path)).alias("dict_path"),
        pl.lit(export_df.height).alias("num_rows"),
        pl.lit(export_df.width).alias("num_cols"),
        pl.lit("success").alias("status"),
    )
    summary.write_csv(summary_path, separator=input_config["delimiter"])
    if logger is not None:
        logger.record(
            "OUTPUT", "adapter_input", rows=export_df.height,
            table=str(tsv_path), mapping=str(dict_path), summary=str(summary_path),
        )
    return summary
