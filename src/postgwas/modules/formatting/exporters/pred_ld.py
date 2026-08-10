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
    usable, excluded = complete_rows(
        frame,
        schema.validation.required_columns,
        validation_expression(schema.validation),
    )
    output = mapped_table(usable, schema, config.minimum_p_value)
    folder = configured_output_path(
        output_directory, schema.output_directory, dataset_id=dataset_id,
    )
    partition_output_column = schema.columns[schema.partition_by]
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
        "columns": output.columns,
    }


__all__ = ["export_pred_ld"]
