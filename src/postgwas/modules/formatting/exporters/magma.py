"""MAGMA p-value and SNP-location export."""

from postgwas.core.io.tables import write_dataframe_table
from postgwas.core.paths import configured_output_path
from postgwas.modules.formatting.contracts import NAMED_OUTPUT_RESULT_KEYS
from postgwas.modules.formatting.table import (
    complete_rows,
    count_bounded_negative_log10_values,
    mapped_table,
    transformation_source,
    validation_expression,
    FormattingError,
)


def export_magma(frame, output_directory, dataset_id, config, *, overwrite):
    """Write the two inputs consumed by PostGWAS' MAGMA runner."""
    schema = config.exports["magma"]
    result_keys = NAMED_OUTPUT_RESULT_KEYS["magma"]
    usable, excluded = complete_rows(
        frame,
        schema.validation.required_columns,
        validation_expression(schema.validation),
    )
    location_schema = schema.outputs["snp_location"]
    p_value_schema = schema.outputs["p_values"]
    p_source = transformation_source(
        p_value_schema, "negative_log10_to_raw_p",
    )
    p_bounded = count_bounded_negative_log10_values(
        usable, p_source, config.minimum_p_value,
    )
    location_table = mapped_table(
        usable, location_schema, config.minimum_p_value,
    )
    p_value_table = mapped_table(
        usable, p_value_schema, config.minimum_p_value,
    )
    location = write_dataframe_table(
        location_table,
        configured_output_path(
            output_directory, location_schema.output_file, dataset_id=dataset_id,
        ),
        overwrite=overwrite,
        runtime=config.runtime,
        error_type=FormattingError,
    )
    p_values = write_dataframe_table(
        p_value_table,
        configured_output_path(
            output_directory, p_value_schema.output_file, dataset_id=dataset_id,
        ),
        overwrite=overwrite,
        runtime=config.runtime,
        error_type=FormattingError,
    )
    return {
        result_keys["snp_location"]: location,
        result_keys["p_values"]: p_values,
        "rows_in": frame.height,
        "rows_out": usable.height,
        "rows_excluded": excluded,
        "p_values_bounded": p_bounded,
        "columns": {
            "snp_location": location_table.columns,
            "p_values": p_value_table.columns,
        },
    }


__all__ = ["export_magma"]
