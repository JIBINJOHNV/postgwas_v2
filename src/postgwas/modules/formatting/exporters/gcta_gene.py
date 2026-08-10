"""GCTA fastBAT and mBAT-combo input export."""

from postgwas.core.paths import configured_output_path
from postgwas.core.io.tables import write_dataframe_table
from postgwas.modules.formatting.contracts import (
    GCTA_SAMPLE_SIZE_MODE,
    NAMED_OUTPUT_RESULT_KEYS,
)
from postgwas.modules.formatting.table import (
    FormattingError,
    complete_rows,
    count_bounded_negative_log10_values,
    mapped_table,
    transformation_source,
    validation_expression,
)


def _export_table(
    frame, output_directory, dataset_id, config, schema, *, overwrite,
):
    usable, excluded = complete_rows(
        frame,
        schema.validation.required_columns,
        validation_expression(schema.validation),
    )
    p_source = transformation_source(schema, "negative_log10_to_raw_p")
    bounded = count_bounded_negative_log10_values(
        usable, p_source, config.minimum_p_value,
    )
    table = mapped_table(usable, schema, config.minimum_p_value)
    path = write_dataframe_table(
        table,
        configured_output_path(
            output_directory, schema.output_file, dataset_id=dataset_id,
        ),
        overwrite=overwrite,
        runtime=config.runtime,
        error_type=FormattingError,
    )
    return {
        "path": path,
        "rows_out": usable.height,
        "rows_excluded": excluded,
        "p_values_bounded": bounded,
        "columns": table.columns,
    }


def export_gcta_gene(frame, output_directory, dataset_id, config, *, overwrite):
    """Write the official GCTA .ma input shared by fastBAT and mBAT-combo."""
    schema = config.exports["gcta_gene"]
    result_keys = NAMED_OUTPUT_RESULT_KEYS["gcta_gene"]
    summary_statistics = _export_table(
        frame,
        output_directory,
        dataset_id,
        config,
        schema.outputs["summary_statistics"],
        overwrite=overwrite,
    )
    return {
        result_keys["summary_statistics"]: summary_statistics["path"],
        "rows_in": frame.height,
        "rows_out": summary_statistics["rows_out"],
        "rows_excluded": summary_statistics["rows_excluded"],
        "p_values_bounded": summary_statistics["p_values_bounded"],
        "outputs": {"summary_statistics": summary_statistics},
        "columns": {"summary_statistics": summary_statistics["columns"]},
        "sample_size_mode": GCTA_SAMPLE_SIZE_MODE,
    }


__all__ = ["export_gcta_gene"]
