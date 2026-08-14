"""Typed configuration for GCTA-COJO conditional and joint analysis."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Literal, get_args

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import DelimitedTableReadConfig, StrictModel
from postgwas.config.models.gcta import GctaBackedModuleConfig, GctaReferenceConfig


GctaCojoMode = Literal["slct", "top_snps", "joint", "cond"]


class GctaCojoInputValidationConfig(StrictModel):
    minimum_reference_overlap_fraction: float = Field(
        gt=0, le=1, allow_inf_nan=False,
    )
    reference_sample_size_warning_threshold: int = Field(ge=1)
    sqlite_batch_size: int = Field(ge=1)
    table_delimiter_pattern: str
    temporary_database_prefix: str
    temporary_database_suffix: str

    @field_validator("table_delimiter_pattern")
    @classmethod
    def valid_delimiter_pattern(cls, value: str) -> str:
        try:
            re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        return value

    @field_validator("temporary_database_prefix", "temporary_database_suffix")
    @classmethod
    def nonempty_filename_fragment(cls, value: str) -> str:
        if not value.strip() or "/" in value or "\\" in value:
            raise ValueError("must be a non-empty filename fragment")
        return value


class GctaCojoAnalysisConfig(StrictModel):
    significance_threshold: float = Field(gt=0, le=1, allow_inf_nan=False)
    top_snp_count: int = Field(ge=1)
    window_kb: int = Field(ge=1)
    collinearity_cutoff: float = Field(gt=0, le=1, allow_inf_nan=False)
    frequency_difference_max: float = Field(gt=0, le=1, allow_inf_nan=False)
    reference_maf_min: float = Field(ge=0, le=0.5, allow_inf_nan=False)
    genomic_control: bool
    genomic_control_lambda: float | None = Field(
        default=None, gt=0, allow_inf_nan=False,
    )
    chromosome: str | None = None

    @field_validator("chromosome")
    @classmethod
    def optional_nonempty_chromosome(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value

    @model_validator(mode="after")
    def genomic_control_lambda_requires_control(self):
        if self.genomic_control_lambda is not None and not self.genomic_control:
            raise ValueError(
                "genomic_control_lambda requires genomic_control: true"
            )
        return self


class GctaCojoInputsConfig(StrictModel):
    condition_snps: str | None = None
    joint_snps: str | None = None
    extract_snps: str | None = None
    exclude_snps: str | None = None

    @field_validator("*")
    @classmethod
    def optional_nonempty_path(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value


class GctaCojoResultSchema(StrictModel):
    required_columns: list[str]
    numeric_columns: list[str]
    required_numeric_columns: list[str]
    positive_columns: list[str]
    integer_columns: list[str]
    closed_unit_interval_columns: list[str]
    identifier_column: str
    chromosome_column: str
    position_column: str
    marginal_p_value_column: str
    p_value_column: str
    adjusted_effect_column: str
    adjusted_standard_error_column: str

    @field_validator(
        "required_columns",
        "numeric_columns",
        "required_numeric_columns",
        "positive_columns",
        "integer_columns",
        "closed_unit_interval_columns",
    )
    @classmethod
    def unique_required_columns(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique columns")
        return values

    @model_validator(mode="after")
    def result_columns_are_required(self):
        for field in (
            self.identifier_column,
            self.chromosome_column,
            self.position_column,
            self.marginal_p_value_column,
            self.p_value_column,
            self.adjusted_effect_column,
            self.adjusted_standard_error_column,
        ):
            if field not in self.required_columns:
                raise ValueError("result statistic columns must be required")
        if not set(self.required_numeric_columns).issubset(self.numeric_columns):
            raise ValueError("required_numeric_columns must be numeric_columns")
        if not set(self.numeric_columns).issubset(self.required_columns):
            raise ValueError("numeric_columns must be required_columns")
        constrained = (
            set(self.positive_columns)
            | set(self.integer_columns)
            | set(self.closed_unit_interval_columns)
        )
        if not constrained.issubset(self.required_numeric_columns):
            raise ValueError(
                "constrained numeric columns must be required_numeric_columns"
            )
        return self


class GctaCojoResultsConfig(DelimitedTableReadConfig):
    normalized_delimiter: str
    normalized_null_value: str
    atomic_output_suffix: str
    estimation_status_column: str
    estimated_status: str
    not_estimable_status: str
    schemas: dict[GctaCojoMode, GctaCojoResultSchema]

    @field_validator("normalized_delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1 or value in "\r\n":
            raise ValueError("must contain one non-newline character")
        return value

    @field_validator(
        "estimation_status_column", "estimated_status", "not_estimable_status",
    )
    @classmethod
    def nonempty_text(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @field_validator("atomic_output_suffix")
    @classmethod
    def safe_atomic_output_suffix(cls, value: str) -> str:
        if not value.strip() or Path(value).name != value:
            raise ValueError("must be a non-empty filename suffix")
        return value

    @model_validator(mode="after")
    def every_mode_has_a_schema(self):
        expected = set(get_args(GctaCojoMode))
        if set(self.schemas) != expected:
            raise ValueError(
                "schemas must contain every COJO mode exactly once: %s"
                % ", ".join(sorted(expected))
            )
        return self


class GctaCojoOutputLayout(StrictModel):
    log_file: str
    resolved_config_file: str
    completion_manifest: str
    output_prefix: str
    staging_directory: str
    normalized_result: str
    primary_results: dict[GctaCojoMode, str]
    ld_matrix_result: str
    gcta_log_result: str

    @field_validator("*")
    @classmethod
    def safe_relative_paths(cls, value):
        values = value.values() if isinstance(value, dict) else (value,)
        for item in values:
            path = Path(str(item).strip())
            if not str(item).strip() or path.is_absolute() or ".." in path.parts:
                raise ValueError(
                    "output paths must be relative and remain below the output directory"
                )
        return value

    @model_validator(mode="after")
    def required_output_tokens(self):
        token_fields = {
            "log_file": ("{dataset_id}", "{mode}"),
            "completion_manifest": ("{dataset_id}", "{mode}"),
            "output_prefix": ("{dataset_id}", "{mode}"),
            "staging_directory": ("{dataset_id}", "{mode}"),
            "normalized_result": ("{dataset_id}", "{mode}"),
            "ld_matrix_result": ("{dataset_id}", "{mode}"),
            "gcta_log_result": ("{dataset_id}", "{mode}"),
        }
        for field, tokens in token_fields.items():
            missing = [token for token in tokens if token not in getattr(self, field)]
            if missing:
                raise ValueError(
                    "%s must contain %s"
                    % (field, ", ".join(repr(token) for token in missing))
                )
        expected = set(get_args(GctaCojoMode))
        if set(self.primary_results) != expected:
            raise ValueError(
                "primary_results must contain every COJO mode exactly once"
            )
        for mode, pattern in self.primary_results.items():
            if "{dataset_id}" not in pattern or "{mode}" not in pattern:
                raise ValueError(
                    "primary_results.%s must contain {dataset_id} and {mode}" % mode
                )
        return self


class GctaCojoReportingConfig(StrictModel):
    top_result_count: int = Field(ge=1)
    p_value_significant_digits: int = Field(ge=1, le=10)
    finding_threshold: float = Field(gt=0, le=1, allow_inf_nan=False)


class GctaCojoLogParsingConfig(StrictModel):
    warning_pattern: str
    maximum_warning_messages: int = Field(ge=1)

    @field_validator("warning_pattern")
    @classmethod
    def valid_warning_pattern(cls, value: str) -> str:
        try:
            re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        return value


class GctaCojoConfig(GctaBackedModuleConfig):
    mode: GctaCojoMode
    input_file: str | None = None
    genome_build: str | None = None
    reference: GctaReferenceConfig
    inputs: GctaCojoInputsConfig
    input_validation: GctaCojoInputValidationConfig
    analysis: GctaCojoAnalysisConfig
    reporting: GctaCojoReportingConfig
    log_parsing: GctaCojoLogParsingConfig
    output_layout: GctaCojoOutputLayout
    results: GctaCojoResultsConfig

    @field_validator("input_file", "genome_build")
    @classmethod
    def optional_nonempty_text(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value

    @model_validator(mode="after")
    def mode_specific_inputs(self):
        if self.mode == "cond" and self.inputs.condition_snps is None:
            raise ValueError("cond mode requires inputs.condition_snps")
        if self.mode != "cond" and self.inputs.condition_snps is not None:
            raise ValueError("inputs.condition_snps is valid only in cond mode")
        if self.mode == "joint" and self.inputs.joint_snps is None:
            raise ValueError("joint mode requires inputs.joint_snps")
        if self.mode != "joint" and self.inputs.joint_snps is not None:
            raise ValueError("inputs.joint_snps is valid only in joint mode")
        if self.mode == "joint" and self.inputs.extract_snps is not None:
            raise ValueError(
                "joint mode uses inputs.joint_snps as GCTA --extract; do not also "
                "set inputs.extract_snps"
            )
        return self


__all__ = ["GctaCojoConfig", "GctaCojoMode"]
