import os
from dataclasses import dataclass
from pathlib import Path
import tempfile
from typing import Optional

import polars as pl

from postgwas.core.paths import configured_output_path

from .rejects import SOURCE_INPUT_ROW_COLUMN
from .shared.variant_columns import has_canonical_variant_columns


@dataclass(frozen=True)
class ChromosomePartitions:
    """Validated chromosome paths and their exact retained-row counts."""

    paths: dict[str, str]
    rows: dict[str, int]


def _validate_source_row_ids(frame: pl.DataFrame, label: str) -> None:
    """Require the one-to-one provenance key used by every chromosome."""
    if SOURCE_INPUT_ROW_COLUMN not in frame.columns:
        raise KeyError(
            "%s has no stable source-row ID column %r."
            % (label, SOURCE_INPUT_ROW_COLUMN)
        )
    identifiers = frame.get_column(SOURCE_INPUT_ROW_COLUMN)
    if identifiers.null_count() or identifiers.n_unique() != frame.height:
        raise ValueError(
            "%s stable source-row IDs must be non-null and unique before "
            "chromosome partitioning." % label
        )


def _write_parquet_atomically(
    frame: pl.DataFrame,
    destination: Path,
    *,
    compression: str,
) -> None:
    """Publish one validated Parquet shard without exposing a partial file."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=".%s." % destination.name,
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
        frame.write_parquet(temporary_path, compression=compression)

        scan = pl.scan_parquet(temporary_path)
        observed_schema = list(dict(scan.schema).items())
        expected_schema = list(frame.schema.items())
        if observed_schema != expected_schema:
            raise ValueError(
                "Parquet schema validation failed for chromosome work file %s: "
                "observed %r, expected %r."
                % (destination, observed_schema, expected_schema)
            )
        observed_rows = int(
            scan.select(pl.len().alias("rows")).collect().item()
        )
        if observed_rows != frame.height:
            raise ValueError(
                "Parquet row-count validation failed for chromosome work file "
                "%s: observed %d, expected %d."
                % (destination, observed_rows, frame.height)
            )
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def write_chromosome_partitions(
    df: pl.DataFrame,
    sample_gwas_dict: dict,
    output_layout: dict[str, str],
    compression: str,
    source_snapshot: Optional[pl.DataFrame] = None,
) -> ChromosomePartitions:
    """
    Write one bounded-memory, typed Parquet shard per observed chromosome.

    Parameters
    ----------
    df : pl.DataFrame
        Full GWAS summary statistics dataframe.
    sample_gwas_dict : dict
        Dictionary with 'chr_col' and 'gwas_outputname'.
    output_layout : dict[str, str]
        Canonical output-path patterns from the resolved YAML configuration.
    compression : str
        Schema-validated Parquet compression from the canonical policy set.
    source_snapshot : pl.DataFrame, optional
        Immutable parsed study rows keyed by the stable source-row ID. When
        supplied, a compressed Parquet snapshot is written for each chromosome
        so rejected rows can recover pre-harmonisation values efficiently.

    Returns
    -------
    ChromosomePartitions
        Chromosome-to-file mapping and the exact number of retained study rows
        written to each validated partition. The counts are returned from the
        same write loop so dataset reconciliation does not rescan the files.
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
    null_chromosomes = df.get_column(chr_col).null_count()
    if null_chromosomes:
        raise ValueError(
            "Working data contains %d row(s) with a null chromosome before "
            "partitioning; PostGWAS will not discard them silently."
            % null_chromosomes
        )

    observed_chromosomes = [
        str(value)
        for value in df.get_column(chr_col).unique(maintain_order=True).to_list()
    ]

    if source_snapshot is not None:
        _validate_source_row_ids(df, "Working data")
        _validate_source_row_ids(source_snapshot, "Source snapshot")
        if (
            df.schema[SOURCE_INPUT_ROW_COLUMN]
            != source_snapshot.schema[SOURCE_INPUT_ROW_COLUMN]
        ):
            raise ValueError(
                "Working data and source snapshot use different stable source-row "
                "ID types (%s versus %s)."
                % (
                    df.schema[SOURCE_INPUT_ROW_COLUMN],
                    source_snapshot.schema[SOURCE_INPUT_ROW_COLUMN],
                )
            )

    output_files = {}
    rows_by_chromosome = {}
    working_rows_written = 0
    source_rows_written = 0
    for chrom in observed_chromosomes:
        working_subset = df.filter(pl.col(chr_col) == chrom)
        if working_subset.is_empty():
            raise ValueError(
                "Observed chromosome %s produced an empty working partition."
                % chrom
            )

        if source_snapshot is not None:
            identifiers = working_subset.select(SOURCE_INPUT_ROW_COLUMN)
            source_subset = (
                source_snapshot.join(
                    identifiers,
                    on=SOURCE_INPUT_ROW_COLUMN,
                    how="semi",
                )
                .sort(SOURCE_INPUT_ROW_COLUMN)
            )
            if source_subset.height != working_subset.height:
                raise ValueError(
                    "Chromosome %s has %d working rows but %d matching immutable "
                    "source rows; every retained row must match exactly once."
                    % (chrom, working_subset.height, source_subset.height)
                )
            source_path = configured_output_path(
                output_folder,
                output_layout["chromosome_source_snapshot"],
                dataset_id=gwas_name,
                chromosome=chrom,
            )
            _write_parquet_atomically(
                source_subset,
                source_path,
                compression=compression,
            )
            source_rows_written += source_subset.height
            del source_subset

        out_path = configured_output_path(
            output_folder,
            output_layout["chromosome_table"],
            dataset_id=gwas_name,
            chromosome=chrom,
        )
        _write_parquet_atomically(
            working_subset,
            out_path,
            compression=compression,
        )
        working_rows_written += working_subset.height
        output_files[str(chrom)] = str(out_path)
        rows_by_chromosome[str(chrom)] = int(working_subset.height)
        del working_subset

    if working_rows_written != df.height:
        raise ValueError(
            "Chromosome working partitions contain %d of %d retained rows."
            % (working_rows_written, df.height)
        )
    if source_snapshot is not None and source_rows_written != df.height:
        raise ValueError(
            "Chromosome source snapshots contain %d of %d retained rows."
            % (source_rows_written, df.height)
        )
    return ChromosomePartitions(
        paths=output_files,
        rows=rows_by_chromosome,
    )
