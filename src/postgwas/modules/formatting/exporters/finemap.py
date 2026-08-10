"""FINEMAP summary-statistics export."""

from postgwas.modules.formatting.table import (
    FormattingError,
    complete_rows,
    mapped_table,
    validation_expression,
)
from postgwas.core.io.tables import write_dataframe_table
from postgwas.core.paths import configured_output_path


def export_finemap(frame, output_directory, dataset_id, config, *, overwrite):
    """Write FINEMAP fields with ALT as effect allele and true MAF."""
    schema = config.exports["finemap"]
    usable, excluded = complete_rows(
        frame,
        schema.validation.required_columns,
        validation_expression(schema.validation),
    )
    output = mapped_table(usable, schema, config.minimum_p_value)
    path = write_dataframe_table(
        output,
        configured_output_path(
            output_directory, schema.output_file, dataset_id=dataset_id,
        ),
        overwrite=overwrite,
        runtime=config.runtime,
        error_type=FormattingError,
    )
    return {
        "finemap_input": path,
        "rows_in": frame.height,
        "rows_out": output.height,
        "rows_excluded": excluded,
        "columns": output.columns,
    }


__all__ = ["export_finemap"]
