import re
from pathlib import Path
from typing import Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import (
    ModuleConfig,
    MultipleTestingSelectionConfig,
    PValueCorrectionMethod,
    StrictModel,
)


MagmaCovarDirection = Literal["two-sided", "greater", "smaller"]
MagmaCovarMissingValues = Literal["drop", "median", "mean"]
MagmaCovarMissingGenes = Literal["drop", "fill"]


class MagmaCovarInputConfig(StrictModel):
    gene_results_file: Path | None = None
    covariates_file: Path | None = None
    missing_values: MagmaCovarMissingValues
    maximum_missing_fraction: float = Field(ge=0, le=0.2)
    missing_genes: MagmaCovarMissingGenes


class MagmaCovarOutputLayout(StrictModel):
    output_prefix: str
    results_file: str
    corrected_results_file: str
    native_log_file: str
    service_log_file: str
    resolved_config_file: str
    completion_manifest: str
    staging_directory: str


class MagmaCovarMultipleTestingConfig(MultipleTestingSelectionConfig):
    reporting_method_labels: dict[PValueCorrectionMethod, str]

    @field_validator("reporting_method_labels")
    @classmethod
    def nonempty_reporting_labels(cls, values: dict[str, str]) -> dict[str, str]:
        if not values or any(not label.strip() for label in values.values()):
            raise ValueError("must provide non-empty reporting labels")
        return values

    @model_validator(mode="after")
    def labels_cover_selected_methods(self):
        missing = sorted(set(self.methods) - set(self.reporting_method_labels))
        if missing:
            raise ValueError(
                "reporting_method_labels is missing configured methods: %s"
                % ", ".join(missing)
            )
        return self


class MagmaCovarResultSchema(StrictModel):
    variable_column: str
    type_column: str
    gene_count_column: str
    beta_column: str
    standardized_beta_column: str
    standard_error_column: str
    p_value_column: str
    adjusted_p_value_column_pattern: str
    primary_correction_method_column: str
    primary_adjusted_p_value_column: str
    primary_significant_column: str
    delimiter: Literal["\t", ","]
    null_value: str

    @field_validator(
        "variable_column", "type_column", "gene_count_column", "beta_column",
        "standardized_beta_column", "standard_error_column", "p_value_column",
        "primary_correction_method_column", "primary_adjusted_p_value_column",
        "primary_significant_column", "null_value",
    )
    @classmethod
    def nonempty_value(cls, value: str) -> str:
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("must not be empty")
        return cleaned

    @field_validator("adjusted_p_value_column_pattern")
    @classmethod
    def valid_adjusted_p_value_pattern(cls, value: str) -> str:
        cleaned = value.strip()
        if "{method}" not in cleaned:
            raise ValueError("must contain '{method}'")
        try:
            cleaned.format(method="bonferroni")
        except (KeyError, ValueError) as exc:
            raise ValueError("must be a valid format pattern") from exc
        return cleaned


class MagmaCovarReportingConfig(StrictModel):
    highlight_method: PValueCorrectionMethod
    top_property_count: int = Field(ge=1)
    p_value_significant_digits: int = Field(ge=1, le=10)
    effect_significant_digits: int = Field(ge=1, le=10)


class MagmaCovarModelModifierConfig(StrictModel):
    syntax: list[str] = Field(min_length=1)
    requires_value: bool
    description: str = Field(min_length=1)
    published_example: str = Field(min_length=1)
    example: str = Field(min_length=1)

    @field_validator("syntax")
    @classmethod
    def validate_syntax(cls, value: list[str]) -> list[str]:
        cleaned = [item.strip() for item in value]
        if any(not item for item in cleaned):
            raise ValueError("model modifier syntax entries must not be empty")
        if len(set(cleaned)) != len(cleaned):
            raise ValueError("model modifier syntax entries must be unique")
        return cleaned

    @field_validator("description", "published_example", "example")
    @classmethod
    def validate_text(cls, value: str) -> str:
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("model modifier help text must not be empty")
        return cleaned


class MagmaCovarModelUseCaseConfig(StrictModel):
    label: str = Field(min_length=1)
    question: str = Field(min_length=1)
    model: list[str]
    direction: MagmaCovarDirection

    @field_validator("label", "question")
    @classmethod
    def validate_text(cls, value: str) -> str:
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("model use-case text must not be empty")
        return cleaned

    @field_validator("model")
    @classmethod
    def validate_model_entries(cls, value: list[str]) -> list[str]:
        return _clean_model_entries(value)


def _clean_model_entries(value: list[str]) -> list[str]:
    cleaned = [item.strip() for item in value]
    if any(not item or item.startswith("-") for item in cleaned):
        raise ValueError(
            "entries must be MAGMA --model modifiers without a leading dash"
        )
    return cleaned


class MagmaCovarConfig(ModuleConfig):
    model: list[str]
    model_use_cases: dict[str, MagmaCovarModelUseCaseConfig]
    model_placeholders: dict[str, str]
    model_modifiers: dict[str, MagmaCovarModelModifierConfig]
    direction: MagmaCovarDirection
    minimum_genes: int = Field(ge=1)
    input: MagmaCovarInputConfig
    multiple_testing: MagmaCovarMultipleTestingConfig
    result_schema: MagmaCovarResultSchema
    reporting: MagmaCovarReportingConfig
    output_layout: MagmaCovarOutputLayout

    @field_validator("model")
    @classmethod
    def validate_model_entries(cls, value: list[str]) -> list[str]:
        return _clean_model_entries(value)

    @field_validator("model_placeholders")
    @classmethod
    def validate_model_placeholders(cls, value: dict[str, str]) -> dict[str, str]:
        if not value:
            raise ValueError("model_placeholders must declare at least one value")
        cleaned: dict[str, str] = {}
        for name, description in value.items():
            if not re.fullmatch(r"[A-Z][A-Z0-9_]*", name):
                raise ValueError("model placeholder names must be uppercase")
            text = description.strip()
            if not text:
                raise ValueError("model placeholder descriptions must not be empty")
            cleaned[name] = text
        return cleaned

    @model_validator(mode="after")
    def validate_model_modifier_catalog(self):
        if not self.model_modifiers:
            raise ValueError("model_modifiers must declare at least one modifier")
        for name, specification in self.model_modifiers.items():
            if not re.fullmatch(r"[a-z][a-z0-9-]*", name):
                raise ValueError(
                    "model_modifiers keys must be lowercase MAGMA modifier names"
                )
            if any(
                syntax != name and not syntax.startswith(name + "=")
                for syntax in specification.syntax
            ):
                raise ValueError(
                    "model_modifiers.%s syntax must start with %s" % (name, name)
                )

        used_placeholders = {
            placeholder
            for specification in self.model_modifiers.values()
            for syntax in specification.syntax
            for placeholder in re.findall(r"<([A-Z][A-Z0-9_]*)>", syntax)
        }
        missing_placeholders = used_placeholders - self.model_placeholders.keys()
        if missing_placeholders:
            raise ValueError(
                "model_placeholders must define: %s"
                % ", ".join(sorted(missing_placeholders))
            )

        if not self.model_use_cases:
            raise ValueError("model_use_cases must declare at least one use case")
        selections = [("model", self.model)]
        for use_case, specification in self.model_use_cases.items():
            if not re.fullmatch(r"[a-z][a-z0-9_]*", use_case):
                raise ValueError(
                    "model_use_cases keys must be lowercase names with underscores"
                )
            selections.append(
                ("model_use_cases.%s.model" % use_case, specification.model)
            )
        for location, options in selections:
            self._validate_model_selection(options, location)
        return self

    def _validate_model_selection(self, options: list[str], location: str) -> None:
        available = ", ".join(self.model_modifiers)
        for option in options:
            name, separator, raw_value = option.partition("=")
            specification = self.model_modifiers.get(name)
            if specification is None:
                raise ValueError(
                    "%s entry %r is not available for gene-property analysis; "
                    "available modifiers are: %s" % (location, option, available)
                )
            if specification.requires_value and (not separator or not raw_value):
                raise ValueError("%s modifier %r requires a value" % (location, name))
            if not specification.requires_value and separator:
                raise ValueError(
                    "%s modifier %r does not accept a value" % (location, name)
                )

    @model_validator(mode="after")
    def output_columns_do_not_collide(self):
        schema = self.result_schema
        fixed = [
            schema.variable_column,
            schema.type_column,
            schema.gene_count_column,
            schema.beta_column,
            schema.standardized_beta_column,
            schema.standard_error_column,
            schema.p_value_column,
            schema.primary_correction_method_column,
            schema.primary_adjusted_p_value_column,
            schema.primary_significant_column,
        ]
        generated = [
            schema.adjusted_p_value_column_pattern.format(method=method)
            for method in self.multiple_testing.methods
        ]
        columns = fixed + generated
        if len(columns) != len(set(columns)):
            raise ValueError(
                "configured MAGMAcovar source, correction, and primary result "
                "columns must be unique"
            )
        if self.reporting.highlight_method not in self.multiple_testing.methods:
            raise ValueError(
                "reporting.highlight_method must occur in "
                "multiple_testing.methods"
            )
        return self
