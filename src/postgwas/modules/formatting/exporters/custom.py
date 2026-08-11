"""Additional user-named summary-statistics table export."""

from postgwas.config.models.modules.formatting import FormattingValidationConfig
from postgwas.core.io.tables import write_dataframe_table
from postgwas.core.paths import configured_output_path
from postgwas.modules.formatting.table import (
    FormattingColumnSpec,
    FormattingError,
    complete_rows,
    count_bounded_negative_log10_values,
    mapped_table,
    validation_expression,
)


def export_custom(frame, output_directory, dataset_id, config, *, overwrite):
    """Write configured CLI fields in their command-line order."""
    custom = config.custom_output
    if not custom.active:
        raise FormattingError("Custom output was requested without an output filename.")

    contracts = [
        (role, custom.field_contracts[role], destination)
        for role, destination in custom.columns.items()
    ]
    required = list(dict.fromkeys(contract.source for _, contract, _ in contracts))
    validation_values = {
        "required_columns": required,
        "positive_columns": [],
        "nonnegative_columns": [],
        "open_unit_interval_columns": [],
        "closed_unit_interval_columns": [],
    }
    for _, contract, _ in contracts:
        if contract.constraint == "none":
            continue
        destination = "%s_columns" % contract.constraint
        if contract.source not in validation_values[destination]:
            validation_values[destination].append(contract.source)
    validation = FormattingValidationConfig.model_validate(validation_values)
    usable, excluded = complete_rows(
        frame,
        validation.required_columns,
        validation_expression(validation),
    )
    specs = [
        FormattingColumnSpec(
            source=contract.source,
            destination=destination,
            transformation=contract.transformation,
        )
        for _, contract, destination in contracts
    ]
    output = mapped_table(
        usable,
        custom,
        config.minimum_p_value,
        column_specs=specs,
    )
    bounded_sources = {
        contract.source
        for _, contract, _ in contracts
        if contract.transformation == "negative_log10_to_raw_p"
    }
    bounded = sum(
        count_bounded_negative_log10_values(
            usable, source, config.minimum_p_value,
        )
        for source in bounded_sources
    )
    path = write_dataframe_table(
        output,
        configured_output_path(
            output_directory,
            custom.output_file,
            error_type=FormattingError,
            dataset_id=dataset_id,
        ),
        overwrite=overwrite,
        runtime=config.runtime,
        error_type=FormattingError,
    )
    return {
        "custom_output_file": path,
        "rows_in": frame.height,
        "rows_out": output.height,
        "rows_excluded": excluded,
        "p_values_bounded": bounded,
        "columns": output.columns,
        "field_roles": dict(custom.columns),
    }


__all__ = ["export_custom"]
