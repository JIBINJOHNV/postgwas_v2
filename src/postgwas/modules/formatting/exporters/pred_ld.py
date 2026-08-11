"""Chromosome-partitioned PRED-LD input export."""

import polars as pl
from postgwas.core.paths import configured_output_path
from postgwas.core.io.tables import write_dataframe_table

from postgwas.modules.formatting.table import (
    complete_rows,
    mapped_table,
    validation_expression,
    FormattingError,
)


def export_pred_ld(frame, output_directory, dataset_id, config, *, overwrite):
    """Write PRED-LD inputs, with NC consistently meaning control count."""
    schema = config.exports["pred_ld"]
    partition_column = schema.partition_by
    if not partition_column or partition_column not in schema.columns:
        raise FormattingError(
            "PRED-LD partition_by must name a canonical column in the "
            "configured PRED-LD column mapping."
        )
    if partition_column not in frame.columns:
        raise FormattingError(
            "Required PRED-LD chromosome field is missing: %s"
            % partition_column
        )

    configured_chromosomes = list(config.chromosomes)
    configured_mask = (
        pl.col(partition_column)
        .cast(pl.String)
        .is_in(configured_chromosomes)
        .fill_null(False)
    )
    usable_chromosomes = frame.filter(configured_mask)
    unconfigured_counts = (
        frame
        .select(
            pl.col(partition_column)
            .cast(pl.String)
            .fill_null(config.runtime.output_null_value),
            configured_mask.alias("_configured"),
        )
        .filter(~pl.col("_configured"))
        .group_by(partition_column)
        .len(name="rows")
        .sort(partition_column)
    )
    excluded_chromosomes = {
        str(row[partition_column]): int(row["rows"])
        for row in unconfigured_counts.iter_rows(named=True)
    }
    excluded_chromosome_rows = sum(excluded_chromosomes.values())
    if usable_chromosomes.is_empty():
        configured = ", ".join(configured_chromosomes) or "<none>"
        excluded_labels = ", ".join(
            "%s=%s" % (label, f"{count:,}")
            for label, count in excluded_chromosomes.items()
        ) or "<none>"
        raise FormattingError(
            "No PRED-LD variants remain on the configured chromosomes: %s. "
            "Excluded unconfigured chromosome counts: %s. Check "
            "modules.formatting.chromosomes and chromosome-label "
            "normalization." % (configured, excluded_labels)
        )

    usable, excluded = complete_rows(
        usable_chromosomes,
        schema.validation.required_columns,
        validation_expression(schema.validation),
    )
    output = mapped_table(usable, schema, config.minimum_p_value)
    folder = configured_output_path(
        output_directory, schema.output_directory, dataset_id=dataset_id,
    )
    partition_output_column = schema.columns[partition_column]
    files = []
    for chromosome in config.chromosomes:
        path = write_dataframe_table(
            output.filter(pl.col(partition_output_column) == str(chromosome)),
            configured_output_path(
                folder,
                schema.partition_file,
                dataset_id=dataset_id,
                partition=chromosome,
            ),
            overwrite=overwrite,
            runtime=config.runtime,
            error_type=FormattingError,
        )
        files.append(path)
    return {
        "pred_ld_folder": str(folder),
        "files": files,
        "rows_in": frame.height,
        "rows_out": output.height,
        "rows_excluded": excluded,
        "rows_excluded_unconfigured_chromosome": excluded_chromosome_rows,
        "excluded_chromosomes": excluded_chromosomes,
        "columns": output.columns,
    }


__all__ = ["export_pred_ld"]
