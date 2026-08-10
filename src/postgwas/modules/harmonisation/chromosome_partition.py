from pathlib import Path
from typing import Optional

import polars as pl

from postgwas.core.paths import configured_output_path

from .rejects import SOURCE_INPUT_ROW_COLUMN
from .shared.variant_columns import has_canonical_variant_columns


def write_chromosome_partitions(
    df: pl.DataFrame,
    sample_gwas_dict: dict,
    output_layout: dict[str, str],
    delimiter: str,
    source_snapshot: Optional[pl.DataFrame] = None,
) -> dict[str, str]:
    """
    Split a GWAS table into per-chromosome files in one partitioning pass.

    Parameters
    ----------
    df : pl.DataFrame
        Full GWAS summary statistics dataframe.
    sample_gwas_dict : dict
        Dictionary with 'chr_col' and 'gwas_outputname'.
    output_layout : dict[str, str]
        Canonical output-path patterns from the resolved YAML configuration.
    delimiter : str
        Configured output-table delimiter.
    source_snapshot : pl.DataFrame, optional
        Immutable parsed study rows keyed by the stable source-row ID. When
        supplied, a compressed Parquet snapshot is written for each chromosome
        so rejected rows can recover pre-harmonisation values efficiently.

    Returns
    -------
    dict[str, str]
        Chromosome-to-file mapping for downstream parallel processing.
    """
    output_folder = Path(sample_gwas_dict["output_folder"])
    output_folder.mkdir(parents=True, exist_ok=True)
    chr_col = sample_gwas_dict.get("chr_col")
    gwas_name = sample_gwas_dict["gwas_outputname"]

    if chr_col not in df.columns:
        raise KeyError(f"❌ Chromosome column '{chr_col}' not found in dataframe.")
    if not has_canonical_variant_columns(df, sample_gwas_dict):
        # Direct library callers may supply a frame outside the validated
        # dataset pipeline; retain their existing compatibility normalization.
        df = df.with_columns(pl.col(chr_col).cast(pl.Utf8).str.strip_chars())
    df = df.filter(pl.col(chr_col).is_not_null())
    source_counts = None
    if source_snapshot is not None:
        if SOURCE_INPUT_ROW_COLUMN not in df.columns:
            raise KeyError(
                "Working data has no stable source-row ID column %r."
                % SOURCE_INPUT_ROW_COLUMN
            )
        if SOURCE_INPUT_ROW_COLUMN not in source_snapshot.columns:
            raise KeyError(
                "Source snapshot has no stable source-row ID column %r."
                % SOURCE_INPUT_ROW_COLUMN
            )
        partition_col = "__postgwas_source_partition_chromosome"
        while partition_col in source_snapshot.columns or partition_col in df.columns:
            partition_col += "_"
        assignments = df.select([
            SOURCE_INPUT_ROW_COLUMN,
            pl.col(chr_col).alias(partition_col),
        ])
        if (
            assignments.get_column(SOURCE_INPUT_ROW_COLUMN).null_count()
            or assignments.get_column(SOURCE_INPUT_ROW_COLUMN).n_unique()
            != assignments.height
        ):
            raise ValueError(
                "Stable source-row IDs must be non-null and unique before chromosome partitioning."
            )
        source_for_partition = source_snapshot.join(
            assignments,
            on=SOURCE_INPUT_ROW_COLUMN,
            how="inner",
        )
        if source_for_partition.height != df.height:
            raise ValueError(
                "Immutable source snapshot matched %d of %d retained rows before "
                "chromosome partitioning; every row must match exactly once."
                % (source_for_partition.height, df.height)
            )
        source_counts = {}
        for key, subset in source_for_partition.partition_by(
            partition_col, as_dict=True, maintain_order=True,
        ).items():
            chrom = str(key[0])
            source_subset = subset.drop(partition_col)
            source_path = configured_output_path(
                output_folder,
                output_layout["chromosome_source_snapshot"],
                dataset_id=gwas_name,
                chromosome=chrom,
            )
            source_path.parent.mkdir(parents=True, exist_ok=True)
            source_subset.write_parquet(source_path)
            source_counts[chrom] = source_subset.height

    output_files = {}
    partitions = df.partition_by(chr_col, as_dict=True, maintain_order=True)
    for key, subset in partitions.items():
        chrom = str(key[0])
        out_path = configured_output_path(
            output_folder,
            output_layout["chromosome_table"],
            dataset_id=gwas_name,
            chromosome=chrom,
        )
        subset.write_csv(out_path, separator=delimiter)
        if source_counts is not None:
            source_count = source_counts.get(chrom)
            if source_count != subset.height:
                raise ValueError(
                    "Chromosome %s has %d working rows but %s immutable source rows."
                    % (
                        chrom,
                        subset.height,
                        0 if source_count is None else source_count,
                    )
                )
        output_files[str(chrom)] = str(out_path)
    return output_files
