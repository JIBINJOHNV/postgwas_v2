"""LDSC munge-sumstats input export."""

import polars as pl
from postgwas.core.io.tables import write_dataframe_table
from postgwas.core.paths import configured_output_path

from postgwas.modules.formatting.table import (
    complete_rows,
    count_bounded_negative_log10_values,
    infer_study_design,
    mapped_table,
    transformation_source,
    validation_expression,
    FormattingError,
)


def export_ldsc(frame, output_directory, dataset_id, config, *, overwrite):
    """Write the sample-size schema understood by CBIIT ``munge_sumstats``."""
    schema = config.exports["ldsc"]
    design = infer_study_design(
        frame,
        config.study_design.case_count_column,
        config.study_design.control_count_column,
    )
    case_control = design.trait_type == "binary"
    sample_columns = schema.trait_required_columns[design.trait_type]
    required = [*schema.validation.required_columns, *sample_columns]
    usable, excluded = complete_rows(
        frame,
        required,
        validation_expression(schema.validation, positive_columns=sample_columns),
    )
    p_source = transformation_source(
        schema, "negative_log10_to_raw_p",
    )
    p_bounded = count_bounded_negative_log10_values(
        usable, p_source, config.minimum_p_value,
    )

    sample_prevalence = None
    if case_control:
        case_column = config.study_design.case_count_column
        control_column = config.study_design.control_count_column
        sizes = usable.select(
            pl.col(case_column).filter(pl.col(case_column) > 0).median().alias("cases"),
            pl.col(control_column)
            .filter(pl.col(control_column) > 0)
            .median()
            .alias("controls"),
        ).to_dicts()[0]
        cases, controls = sizes["cases"], sizes["controls"]
        if cases and controls:
            sample_prevalence = float(cases / (cases + controls))

    mapping = dict(schema.columns)
    mapping.update(schema.trait_columns[design.trait_type])
    mapping.update(schema.trailing_columns)
    if not case_control:
        sample_size_mode = "quantitative_n_from_nco"
    else:
        sample_size_mode = "case_control_counts"
    output = mapped_table(
        usable, schema, config.minimum_p_value, columns=mapping,
    )
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
        "ldsc_file": path,
        "sample_prev": sample_prevalence,
        "trait_type": design.trait_type,
        "sample_size_mode": sample_size_mode,
        "columns": output.columns,
        "rows_in": frame.height,
        "rows_out": output.height,
        "rows_excluded": excluded,
        "p_values_bounded": p_bounded,
    }


__all__ = ["export_ldsc"]
