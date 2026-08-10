"""Typed configuration for the upstream FLAMES integration."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import GenomeBuild, ModuleConfig, StrictModel


class FlamesInputSchema(StrictModel):
    table_delimiter: str
    whitespace_delimiter_pattern: str
    index_file: str
    index_filename_column: str
    index_locus_column: str
    index_annotation_column: str
    credible_set_index_column: Literal["index"]
    credible_set_variant_column: Literal["cred1"]
    credible_set_probability_column: Literal["prob1"]
    variant_identifier_pattern: str
    chromosome_minimum: int = Field(ge=1)
    chromosome_maximum: int = Field(ge=1)
    position_maximum: int = Field(ge=1)
    magma_gene_column: str
    magma_z_column: str
    magma_covariate_column: str
    magma_covariate_p_column: str
    pops_gene_column: str
    pops_score_column: str
    ensembl_gene_pattern: str

    @field_validator("table_delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @field_validator("index_file")
    @classmethod
    def plain_index_filename(cls, value: str) -> str:
        if not value or Path(value).name != value:
            raise ValueError("must be a plain filename")
        return value

    @field_validator(
        "variant_identifier_pattern",
        "ensembl_gene_pattern",
        "whitespace_delimiter_pattern",
    )
    @classmethod
    def valid_pattern(cls, value: str) -> str:
        try:
            re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        return value

    @field_validator("variant_identifier_pattern")
    @classmethod
    def variant_pattern_groups(cls, value: str) -> str:
        if re.compile(value).groups != 4:
            raise ValueError(
                "must capture chromosome, position, allele 1, and allele 2 in order"
            )
        return value

    @field_validator(
        "index_filename_column", "index_locus_column", "index_annotation_column",
        "magma_gene_column", "magma_z_column", "magma_covariate_column",
        "magma_covariate_p_column", "pops_gene_column", "pops_score_column",
    )
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @model_validator(mode="after")
    def valid_chromosome_range(self):
        if self.chromosome_maximum < self.chromosome_minimum:
            raise ValueError("chromosome_maximum must be at least chromosome_minimum")
        return self


class FlamesUpstreamResources(StrictModel):
    script_relative_path: str
    model_directory_relative_path: str
    model_file: str
    feature_file: str
    cadd_index_suffix: str
    required_annotation_directories: list[str]
    required_annotation_file_patterns: list[str]
    runtime_imports: list[str]
    annotation_bundle_id: str

    @field_validator("script_relative_path", "model_directory_relative_path")
    @classmethod
    def safe_relative_path(cls, value: str) -> str:
        path = Path(value)
        if path.is_absolute() or ".." in path.parts or not value:
            raise ValueError("must be a non-empty relative path")
        return value

    @field_validator("model_file", "feature_file")
    @classmethod
    def plain_filename(cls, value: str) -> str:
        if not value or Path(value).name != value:
            raise ValueError("must be a plain filename")
        return value

    @field_validator("cadd_index_suffix")
    @classmethod
    def nonempty_suffix(cls, value: str) -> str:
        if not value.startswith(".") or Path(value).name != value:
            raise ValueError("must be a filename suffix beginning with a period")
        return value

    @field_validator(
        "required_annotation_directories",
        "required_annotation_file_patterns",
        "runtime_imports",
    )
    @classmethod
    def unique_nonempty_values(cls, values: list[str]) -> list[str]:
        if (
            not values
            or len(values) != len(set(values))
            or any(not value for value in values)
        ):
            raise ValueError("must contain unique non-empty values")
        return values

    @field_validator("required_annotation_file_patterns")
    @classmethod
    def valid_build_patterns(cls, values: list[str]) -> list[str]:
        for value in values:
            path = Path(value)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError(
                    "annotation file patterns must remain below the resource directory"
                )
            if "{" in value and value.count("{genome_build}") != 1:
                raise ValueError(
                    "only {genome_build} is supported in annotation file patterns"
                )
        return values

    @field_validator("annotation_bundle_id")
    @classmethod
    def nonempty_bundle_id(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value


class FlamesResultSchema(StrictModel):
    filename_column: str
    locus_column: str
    symbol_column: str
    gene_column: str
    xgboost_score_column: str
    pops_score_column: str
    raw_score_column: str
    scaled_score_column: str
    precision_column: str
    highest_column: str
    causal_column: str

    @field_validator("*")
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @model_validator(mode="after")
    def unique_columns(self):
        values = list(self.model_dump().values())
        if len(values) != len(set(values)):
            raise ValueError("result columns must be unique")
        return self


class FlamesOutputLayout(StrictModel):
    score_basename: str
    raw_score_suffix: str
    prediction_suffix: str
    annotation_file: str
    index_file: str
    staging_directory: str
    resolved_config_file: str
    completion_manifest: str
    service_log_file: str

    @field_validator("*")
    @classmethod
    def safe_relative_pattern(cls, value: str) -> str:
        path = Path(value)
        if path.is_absolute() or ".." in path.parts or not value:
            raise ValueError("must be a non-empty path below the output directory")
        return value

    @model_validator(mode="after")
    def required_tokens(self):
        dataset_fields = (
            "score_basename", "annotation_file", "index_file",
            "staging_directory", "completion_manifest", "service_log_file",
        )
        for field in dataset_fields:
            if "{dataset_id}" not in getattr(self, field):
                raise ValueError("%s must contain {dataset_id}" % field)
        if "{row_number}" not in self.annotation_file:
            raise ValueError("annotation_file must contain {row_number}")
        for field in ("raw_score_suffix", "prediction_suffix"):
            value = getattr(self, field)
            if not value.startswith(".") or Path(value).name != value:
                raise ValueError("%s must be a filename suffix" % field)
        return self


class FlamesConfig(ModuleConfig):
    genome_build: GenomeBuild
    # The bundled classifier and published calibration use a 750 kb window.
    locus_window_bp: Literal[750000]
    minimum_credible_set_coverage: Literal[0.95]
    probability_tolerance: float = Field(gt=0, le=0.01)
    credible_sets_directory: str | None = None
    magma_gene_results_file: str | None = None
    magma_covariate_results_file: str | None = None
    pops_scores_file: str | None = None
    annotation_resource_directory: str | None = None
    model_directory: str | None = None
    vep_mode: Literal["api", "local"]
    vep_command: str | None = None
    vep_cache: str | None = None
    cadd_mode: Literal["api", "local"]
    cadd_file: str | None = None
    input_schema: FlamesInputSchema
    upstream: FlamesUpstreamResources
    result_schema: FlamesResultSchema
    output_layout: FlamesOutputLayout

    @field_validator(
        "credible_sets_directory", "magma_gene_results_file",
        "magma_covariate_results_file", "pops_scores_file",
        "annotation_resource_directory", "model_directory", "vep_command",
        "vep_cache", "cadd_file",
    )
    @classmethod
    def nonempty_optional_path(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value

    @model_validator(mode="after")
    def local_annotation_requirements(self):
        if self.vep_mode == "local" and (
            self.vep_command is None or self.vep_cache is None
        ):
            raise ValueError("local VEP mode requires vep_command and vep_cache")
        if self.cadd_mode == "local" and self.cadd_file is None:
            raise ValueError("local CADD mode requires cadd_file")
        return self


__all__ = ["FlamesConfig", "FlamesMagmaCovarContract"]
