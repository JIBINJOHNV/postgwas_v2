"""Typed configuration for GCTA fastBAT and mBAT-combo gene tests."""

from __future__ import annotations

from pathlib import Path
from typing import get_args, Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import (
    ChromosomeAnalysisScopeConfig,
    DelimitedTableReadConfig,
    MHCAnalysisScopeConfig,
    StrictModel,
)
from postgwas.config.models.gcta import GctaBackedModuleConfig, GctaReferenceConfig


GctaGeneMethod = Literal[
    "fastbat_gene", "fastbat_segment", "fastbat_set", "mbat_combo",
]
GctaChromosomeLabelPolicy = Literal["exact", "strip_chr_prefix"]
GctaDuplicateGenePolicy = Literal["error", "deduplicate"]
GctaUnmappedGenePolicy = Literal["error", "report"]
GctaEmptyPathwayPolicy = Literal["error", "omit"]
GctaOversizedSetPolicy = Literal["error", "omit"]
GctaPathwayAuditLevel = Literal["summary", "normalized", "expanded"]
GctaPathwayAuditFormat = Literal["parquet"]
GctaPathwayAuditCompression = Literal["zstd", "snappy", "gzip", "none"]


def _validate_output_filenames(model: StrictModel, label: str) -> None:
    values = [getattr(model, name) for name in type(model).model_fields]
    if len(values) != len(set(values)):
        raise ValueError("%s output names must be unique" % label)
    for value in values:
        path = Path(value.strip())
        if not value.strip() or path.is_absolute() or len(path.parts) != 1:
            raise ValueError("%s output names must be simple filenames" % label)


class GctaVariantHarmonisationConfig(StrictModel):
    chromosome_label_policy: GctaChromosomeLabelPolicy
    minimum_overlap_fraction: float = Field(
        gt=0, le=1, allow_inf_nan=False,
    )
    allow_strand_complement: bool


class GctaGeneResourceOutputNames(StrictModel):
    manifest: str
    readme: str
    checksums: str

    @model_validator(mode="after")
    def safe_unique_filenames(self):
        _validate_output_filenames(self, "gene resource")
        return self


class GctaGeneAnnotationConfig(StrictModel):
    file: str | None = None
    columns: list[str]
    resource_output_names: GctaGeneResourceOutputNames

    @field_validator("file")
    @classmethod
    def optional_nonempty_text(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value

    @field_validator("columns")
    @classmethod
    def four_unique_columns(cls, values: list[str]) -> list[str]:
        if len(values) != 4 or len(values) != len(set(values)):
            raise ValueError("must contain four unique column names")
        if any(not str(value).strip() for value in values):
            raise ValueError("must not contain empty values")
        required = {"chromosome", "start", "end", "gene"}
        if set(values) != required:
            raise ValueError(
                "must map the gene-list roles: %s" % ", ".join(sorted(required))
            )
        return values


class GctaPathwayOutputNames(StrictModel):
    set_list: str
    pathway_mapping: str
    pathway_gene_mapping: str
    gene_variant_mapping: str
    expanded_mapping: str
    unmapped_genes: str
    manifest: str
    readme: str
    checksums: str

    @model_validator(mode="after")
    def safe_unique_filenames(self):
        _validate_output_filenames(self, "pathway conversion")
        return self


class GctaPathwayParallelismConfig(StrictModel):
    worker_memory_multiplier: float = Field(
        gt=1, allow_inf_nan=False,
    )


class GctaPathwayAuditConfig(StrictModel):
    level: GctaPathwayAuditLevel
    format: GctaPathwayAuditFormat
    compression: GctaPathwayAuditCompression
    batch_rows: int = Field(gt=0)


class GctaPathwayDiskConfig(StrictModel):
    minimum_free_gb: float = Field(ge=0, allow_inf_nan=False)
    estimation_safety_factor: float = Field(ge=1, allow_inf_nan=False)


class GctaPathwayConversionConfig(StrictModel):
    allowed_chromosomes: list[str]
    chromosome_label_policy: GctaChromosomeLabelPolicy
    duplicate_gene_policy: GctaDuplicateGenePolicy
    unmapped_gene_policy: GctaUnmappedGenePolicy
    minimum_gene_id_overlap_fraction: float = Field(
        gt=0, le=1, allow_inf_nan=False,
    )
    empty_pathway_policy: GctaEmptyPathwayPolicy
    parallelism: GctaPathwayParallelismConfig
    audit: GctaPathwayAuditConfig
    disk: GctaPathwayDiskConfig
    output_names: GctaPathwayOutputNames

    @field_validator("allowed_chromosomes")
    @classmethod
    def unique_nonempty_chromosomes(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique chromosome labels")
        if any(not str(value).strip() for value in values):
            raise ValueError("must not contain empty chromosome labels")
        return values


class GctaSetAnnotationConfig(StrictModel):
    file: str | None = None
    gmt_file: str | None = None
    maximum_set_variants: int = Field(gt=0)
    oversized_set_policy: GctaOversizedSetPolicy
    conversion: GctaPathwayConversionConfig

    @field_validator("file", "gmt_file")
    @classmethod
    def optional_nonempty_text(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value

    @model_validator(mode="after")
    def one_set_source(self):
        if self.file is not None and self.gmt_file is not None:
            raise ValueError("file and gmt_file are mutually exclusive")
        return self


class GctaGeneOutputLayout(StrictModel):
    log_file: str
    resolved_config_file: str
    completion_manifest: str
    output_prefix: str
    staging_directory: str
    prepared_set_directory: str
    analysis_set_list: str
    variant_exclusion_list: str
    scoped_gene_list: str
    excluded_gene_list: str
    primary_results: dict[GctaGeneMethod, str]
    normalized_result: str
    html_report: str
    frequency_qc_result: str
    mbat_snpset_result: str

    @field_validator(
        "log_file", "resolved_config_file", "completion_manifest", "output_prefix",
        "staging_directory", "prepared_set_directory", "analysis_set_list",
        "variant_exclusion_list", "scoped_gene_list", "excluded_gene_list",
        "normalized_result", "html_report", "frequency_qc_result",
        "mbat_snpset_result",
    )
    @classmethod
    def safe_relative_pattern(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("must not be empty")
        path = Path(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("must be a relative path below the output directory")
        return value

    @model_validator(mode="after")
    def required_tokens(self):
        token_requirements = {
            "log_file": ("{dataset_id}", "{method}"),
            "completion_manifest": ("{dataset_id}", "{method}"),
            "output_prefix": ("{dataset_id}",),
            "staging_directory": ("{dataset_id}", "{method}"),
            "prepared_set_directory": ("{dataset_id}",),
            "analysis_set_list": ("{dataset_id}",),
            "variant_exclusion_list": ("{dataset_id}", "{method}"),
            "scoped_gene_list": ("{dataset_id}", "{method}"),
            "excluded_gene_list": ("{dataset_id}", "{method}"),
            "normalized_result": ("{dataset_id}", "{method}"),
            "html_report": ("{dataset_id}", "{method}"),
            "frequency_qc_result": ("{dataset_id}",),
            "mbat_snpset_result": ("{dataset_id}",),
        }
        for field, tokens in token_requirements.items():
            missing = [token for token in tokens if token not in getattr(self, field)]
            if missing:
                raise ValueError(
                    "%s must contain %s"
                    % (field, ", ".join(repr(token) for token in missing))
                )
        expected = set(get_args(GctaGeneMethod))
        if set(self.primary_results) != expected:
            raise ValueError(
                "primary_results must contain every method exactly once: %s"
                % ", ".join(sorted(expected))
            )
        for method, pattern in self.primary_results.items():
            self.safe_relative_pattern(pattern)
            if "{dataset_id}" not in pattern:
                raise ValueError(
                    "primary_results.%s must contain {dataset_id}" % method
                )
        if not self.html_report.lower().endswith(".html"):
            raise ValueError("html_report must end with .html")
        return self


class GctaResultSchema(StrictModel):
    required_columns: list[str]
    reportable_columns: list[str]
    p_value_column: str
    identifier_columns: list[str]
    component_p_value_columns: list[str]
    unit_label: str

    @field_validator(
        "required_columns", "reportable_columns", "identifier_columns",
        "component_p_value_columns",
    )
    @classmethod
    def unique_columns(cls, values: list[str]) -> list[str]:
        if len(values) != len(set(values)):
            raise ValueError("must not contain duplicate columns")
        return values

    @model_validator(mode="after")
    def referenced_columns_are_required(self):
        required = set(self.required_columns)
        reportable = set(self.reportable_columns)
        if not required.issubset(reportable):
            raise ValueError("required_columns must all be reportable")
        if self.p_value_column not in required:
            raise ValueError("p_value_column must be required")
        if not self.identifier_columns or not set(self.identifier_columns).issubset(required):
            raise ValueError("identifier_columns must be nonempty required columns")
        if required & set(self.component_p_value_columns):
            raise ValueError("component p-value columns must not duplicate required columns")
        if not set(self.component_p_value_columns).issubset(reportable):
            raise ValueError("component p-value columns must all be reportable")
        if not self.unit_label.strip():
            raise ValueError("unit_label must not be empty")
        return self


class GctaGeneResultConfig(DelimitedTableReadConfig):
    normalized_schema_version: str
    normalized_delimiter: str
    normalized_null_value: str
    schemas: dict[GctaGeneMethod, GctaResultSchema]

    @field_validator("normalized_schema_version")
    @classmethod
    def nonempty_schema_version(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @field_validator("normalized_delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @model_validator(mode="after")
    def every_method_has_a_schema(self):
        expected = set(get_args(GctaGeneMethod))
        if set(self.schemas) != expected:
            raise ValueError(
                "schemas must contain every method exactly once: %s"
                % ", ".join(sorted(expected))
            )
        for method, schema in self.schemas.items():
            if method == "mbat_combo" and not schema.component_p_value_columns:
                raise ValueError(
                    "schemas.mbat_combo must define component p-value columns"
                )
            if method != "mbat_combo" and schema.component_p_value_columns:
                raise ValueError(
                    "only schemas.mbat_combo may define component p-value columns"
                )
        return self


class GctaCorrectionColumns(StrictModel):
    nominal_significant: str
    bonferroni_adjusted_p: str
    bonferroni_significant: str
    fdr_bh_adjusted_p: str
    fdr_bh_significant: str

    @model_validator(mode="after")
    def unique_nonempty_column_names(self):
        values = [getattr(self, name) for name in type(self).model_fields]
        if any(not value.strip() for value in values):
            raise ValueError("correction column names must not be empty")
        if len(values) != len(set(values)):
            raise ValueError("correction column names must be unique")
        return self


class GctaReportingConfig(StrictModel):
    top_result_count: int = Field(ge=1)
    p_value_significant_digits: int = Field(ge=1, le=10)
    nominal_alpha: float = Field(gt=0, le=1, allow_inf_nan=False)
    familywise_alpha: float = Field(gt=0, le=1, allow_inf_nan=False)
    fdr_alpha: float = Field(gt=0, le=1, allow_inf_nan=False)
    correction_columns: GctaCorrectionColumns
    chromosome_columns: dict[GctaGeneMethod, str | None]

    @model_validator(mode="after")
    def every_method_has_a_chromosome_mapping(self):
        expected = set(get_args(GctaGeneMethod))
        if set(self.chromosome_columns) != expected:
            raise ValueError(
                "chromosome_columns must contain every method exactly once: %s"
                % ", ".join(sorted(expected))
            )
        for method, column in self.chromosome_columns.items():
            if column is not None and not column.strip():
                raise ValueError(
                    "chromosome_columns.%s must be null or nonempty" % method
                )
        return self


class GctaHtmlReportConfig(StrictModel):
    """Presentation settings for complete method-specific GCTA result tables."""

    page_size: int = Field(gt=0)
    columns: dict[GctaGeneMethod, list[str]]

    @model_validator(mode="after")
    def every_method_has_unique_columns(self):
        expected = set(get_args(GctaGeneMethod))
        if set(self.columns) != expected:
            raise ValueError(
                "columns must contain every method exactly once: %s"
                % ", ".join(sorted(expected))
            )
        for method, columns in self.columns.items():
            if not columns or len(columns) != len(set(columns)):
                raise ValueError(
                    "columns.%s must contain one or more unique column names"
                    % method
                )
            if any(not str(column).strip() for column in columns):
                raise ValueError(
                    "columns.%s must not contain empty column names" % method
                )
        return self


class GctaGeneConfig(GctaBackedModuleConfig):
    method: GctaGeneMethod
    input_file: str | None = None
    genome_build: str | None = None
    reference: GctaReferenceConfig
    variant_harmonisation: GctaVariantHarmonisationConfig
    mhc: MHCAnalysisScopeConfig
    chromosomes: ChromosomeAnalysisScopeConfig
    gene_annotation: GctaGeneAnnotationConfig
    set_annotation: GctaSetAnnotationConfig
    gene_window_kb: int = Field(ge=0)
    segment_size_kb: int = Field(gt=0)
    reference_maf_min: float = Field(ge=0, le=0.5, allow_inf_nan=False)
    fastbat_ld_cutoff: float = Field(gt=0, le=1, allow_inf_nan=False)
    mbat_svd_gamma: float = Field(gt=0, le=1, allow_inf_nan=False)
    frequency_difference_max: float = Field(gt=0, le=1, allow_inf_nan=False)
    print_component_p_values: bool
    write_snpset: bool
    reporting: GctaReportingConfig
    html_report: GctaHtmlReportConfig
    output_layout: GctaGeneOutputLayout
    results: GctaGeneResultConfig

    @field_validator("input_file", "genome_build")
    @classmethod
    def optional_nonempty_text(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value

    @model_validator(mode="after")
    def reporting_columns_reference_result_schemas(self):
        correction_columns = set(
            self.reporting.correction_columns.model_dump().values()
        )
        for method, schema in self.results.schemas.items():
            conflicts = correction_columns & (
                set(schema.required_columns)
                | set(schema.component_p_value_columns)
            )
            if conflicts:
                raise ValueError(
                    "reporting correction columns conflict with %s result "
                    "columns: %s" % (method, ", ".join(sorted(conflicts)))
                )
        for method, column in self.reporting.chromosome_columns.items():
            if method == "fastbat_set" and column is not None:
                raise ValueError(
                    "reporting.chromosome_columns.fastbat_set must be null"
                )
            if method != "fastbat_set" and column is None:
                raise ValueError(
                    "gene and segment reporting requires a chromosome column"
                )
            if column is not None and column not in self.results.schemas[method].required_columns:
                raise ValueError(
                    "reporting chromosome column for %s must be required by its "
                    "result schema" % method
                )
        correction_order = list(
            self.reporting.correction_columns.model_dump().values()
        )
        for method, columns in self.html_report.columns.items():
            schema = self.results.schemas[method]
            available = set(
                schema.reportable_columns
                + correction_order
            )
            unknown = [column for column in columns if column not in available]
            if unknown:
                raise ValueError(
                    "html_report.columns.%s contains columns not guaranteed by "
                    "the normalized result schema: %s"
                    % (method, ", ".join(unknown))
                )
            required = [
                *schema.identifier_columns,
                schema.p_value_column,
                *correction_order,
            ]
            missing = [column for column in required if column not in columns]
            if missing:
                raise ValueError(
                    "html_report.columns.%s must contain identifier, primary "
                    "p-value, and correction columns: %s"
                    % (method, ", ".join(missing))
                )
        return self


__all__ = [
    "GctaGeneConfig",
    "GctaGeneMethod",
    "GctaOversizedSetPolicy",
]
