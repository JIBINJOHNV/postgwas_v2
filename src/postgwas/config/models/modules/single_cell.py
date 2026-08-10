"""Typed configuration for GWAS and single-cell integration analyses."""

from __future__ import annotations

from pathlib import Path
import re
from typing import get_args, Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import ModuleConfig, StrictModel


SingleCellTool = Literal["magma_celltype", "scdrs"]
SINGLE_CELL_TOOLS = get_args(SingleCellTool)
SingleCellCorrectionMethod = Literal["bonferroni", "sidak", "holm", "fdr_bh"]
ScdrsSpecies = Literal["human", "hsapiens", "mouse", "mmusculus"]
ScdrsMatrixState = Literal["raw_counts", "normalized_log1p"]


def _safe_relative_pattern(value: str) -> str:
    cleaned = value.strip()
    if not cleaned:
        raise ValueError("must not be empty")
    path = Path(cleaned)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError("must be a relative path below the output directory")
    return cleaned


class MagmaCelltypeInputConfig(StrictModel):
    gene_results_file: Path | None = None
    covariates_file: Path | None = None


class MagmaCelltypeConfig(StrictModel):
    workflow: Literal["base"]
    magmacovar_use_case: str
    average_property: str
    input: MagmaCelltypeInputConfig

    @field_validator("magmacovar_use_case", "average_property")
    @classmethod
    def nonempty_name(cls, value: str) -> str:
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("must not be empty")
        if any(character.isspace() for character in cleaned):
            raise ValueError("must not contain whitespace")
        return cleaned


class ScdrsInputConfig(StrictModel):
    h5ad_file: Path | None = None
    gene_set_file: Path | None = None
    magma_gene_results_file: Path | None = None
    gene_identifier_map_file: Path | None = None
    covariate_file: Path | None = None


class ScdrsMagmaGeneSetConfig(StrictModel):
    source: Literal["file", "magma"]
    gene_id_column: str
    z_score_column: str
    statistics_gene_column: str
    statistics_null_value: str
    table_delimiter_pattern: str
    source_identifier_type: str
    target_identifier_type: str
    mapping_mode: Literal["exact", "crosswalk"]
    mapping_source_column: str
    mapping_target_column: str
    mapping_delimiter: str
    unmapped_policy: Literal["error", "exclude"]
    ambiguous_policy: Literal["error"]
    trait_pattern: str
    weight: Literal["zscore", "uniform"]
    minimum_genes: int = Field(ge=1)
    maximum_genes: int = Field(ge=1)

    @field_validator(
        "gene_id_column", "z_score_column", "statistics_gene_column",
        "statistics_null_value",
        "source_identifier_type",
        "target_identifier_type", "mapping_source_column",
        "mapping_target_column",
    )
    @classmethod
    def nonempty_name(cls, value: str) -> str:
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("must not be empty")
        return cleaned

    @field_validator("table_delimiter_pattern")
    @classmethod
    def valid_delimiter_pattern(cls, value: str) -> str:
        try:
            re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        return value

    @field_validator("mapping_delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @field_validator("trait_pattern")
    @classmethod
    def dataset_trait_pattern(cls, value: str) -> str:
        cleaned = value.strip()
        if "{dataset_id}" not in cleaned:
            raise ValueError("must contain the {dataset_id} placeholder")
        try:
            cleaned.format(dataset_id="STUDY")
        except (KeyError, ValueError) as exc:
            raise ValueError("must only use the {dataset_id} placeholder") from exc
        return cleaned

    @model_validator(mode="after")
    def valid_gene_range(self):
        if self.minimum_genes > self.maximum_genes:
            raise ValueError("minimum_genes must not exceed maximum_genes")
        if self.mapping_source_column == self.mapping_target_column:
            raise ValueError("mapping source and target columns must be distinct")
        if (
            self.mapping_mode == "exact"
            and self.source_identifier_type != self.target_identifier_type
        ):
            raise ValueError(
                "mapping_mode exact requires identical source and target "
                "identifier types"
            )
        return self


class ScdrsGeneSetFormatConfig(StrictModel):
    trait_column: str
    gene_set_column: str
    delimiter: Literal["\t"]
    gene_separator: Literal[","]
    weight_separator: Literal[":"]

    @field_validator("trait_column", "gene_set_column")
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("must not be empty")
        return cleaned

    @model_validator(mode="after")
    def distinct_values(self):
        if self.trait_column == self.gene_set_column:
            raise ValueError("trait and gene-set columns must be distinct")
        separators = [
            self.delimiter, self.gene_separator, self.weight_separator,
        ]
        if len(separators) != len(set(separators)):
            raise ValueError("gene-set separators must be distinct")
        return self


class ScdrsCovariateFormatConfig(StrictModel):
    delimiter: Literal["\t"]
    require_exact_cell_ids: bool
    minimum_cell_overlap_fraction: float = Field(
        gt=0, lt=1, allow_inf_nan=False,
    )
    constant_column: str | None = None
    constant_tolerance: float = Field(ge=0, allow_inf_nan=False)

    @field_validator("constant_column")
    @classmethod
    def optional_nonempty_name(cls, value: str | None) -> str | None:
        if value is None:
            return None
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("must not be empty")
        return cleaned


class ScdrsValidationConfig(StrictModel):
    matrix_chunk_rows: int = Field(ge=1)
    raw_count_integer_tolerance: float = Field(
        ge=0, allow_inf_nan=False,
    )
    minimum_effective_genes: int = Field(ge=1)
    maximum_effective_gene_fraction: float = Field(
        gt=0, le=1, allow_inf_nan=False,
    )
    allow_missing_annotation_values: bool


class ScdrsDownstreamConfig(StrictModel):
    group_analysis: list[str]
    correlation_analysis: list[str]
    gene_analysis: bool
    knn_neighbors: int = Field(ge=1)
    knn_principal_components: int = Field(ge=1)

    @field_validator("group_analysis", "correlation_analysis")
    @classmethod
    def unique_annotations(cls, values: list[str]) -> list[str]:
        cleaned = [value.strip() for value in values]
        if any(not value for value in cleaned):
            raise ValueError("annotation names must not be empty")
        if len(cleaned) != len(set(cleaned)):
            raise ValueError("annotation names must be unique")
        return cleaned


class ScdrsConfig(StrictModel):
    input: ScdrsInputConfig
    magma_gene_set: ScdrsMagmaGeneSetConfig
    h5ad_species: ScdrsSpecies
    gene_set_species: ScdrsSpecies
    matrix_state: ScdrsMatrixState
    weight_option: Literal["vs", "uniform"]
    adjust_proportion_column: str | None = None
    filter_data: bool
    control_gene_sets: int = Field(ge=1)
    minimum_genes_per_cell: int = Field(ge=1)
    minimum_cells_per_gene: int = Field(ge=1)
    return_control_raw_score: bool
    return_control_normalized_score: bool
    version_probe_arguments: list[str] = Field(min_length=1)
    supported_versions: list[str] = Field(min_length=1)
    downstream: ScdrsDownstreamConfig
    gene_set_format: ScdrsGeneSetFormatConfig
    covariate_format: ScdrsCovariateFormatConfig
    validation: ScdrsValidationConfig

    @field_validator("adjust_proportion_column")
    @classmethod
    def optional_nonempty_name(cls, value: str | None) -> str | None:
        if value is None:
            return None
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("must not be empty")
        return cleaned

    @field_validator("version_probe_arguments", "supported_versions")
    @classmethod
    def unique_nonempty_values(cls, values: list[str]) -> list[str]:
        cleaned = [value.strip() for value in values]
        if any(not value for value in cleaned):
            raise ValueError("values must not be empty")
        if len(cleaned) != len(set(cleaned)):
            raise ValueError("values must be unique")
        return cleaned

    @model_validator(mode="after")
    def supported_native_contract(self):
        species_family = {
            "human": "human",
            "hsapiens": "human",
            "mouse": "mouse",
            "mmusculus": "mouse",
        }
        if (
            species_family[self.h5ad_species]
            != species_family[self.gene_set_species]
        ):
            raise ValueError(
                "cross-species scDRS gene-set conversion is not supported by "
                "the current PostGWAS preflight; use matching species"
            )
        if not (
            self.return_control_raw_score
            or self.return_control_normalized_score
        ):
            raise ValueError(
                "scDRS requires at least one control-score output so the "
                "expected full_score file is produced"
            )
        downstream_requested = bool(
            self.downstream.group_analysis
            or self.downstream.correlation_analysis
            or self.downstream.gene_analysis
        )
        if downstream_requested and not self.return_control_normalized_score:
            raise ValueError(
                "scDRS downstream analyses require normalized control scores"
            )
        return self


class SingleCellMultipleTestingConfig(StrictModel):
    methods: list[SingleCellCorrectionMethod] = Field(min_length=1)
    primary_method: SingleCellCorrectionMethod
    significance_threshold: float = Field(gt=0, le=1, allow_inf_nan=False)

    @field_validator("methods")
    @classmethod
    def unique_methods(
        cls, values: list[SingleCellCorrectionMethod],
    ) -> list[SingleCellCorrectionMethod]:
        if len(values) != len(set(values)):
            raise ValueError("must contain unique correction methods")
        return values

    @model_validator(mode="after")
    def primary_is_selected(self):
        if self.primary_method not in self.methods:
            raise ValueError("primary_method must occur in methods")
        return self


class SingleCellResultSchema(StrictModel):
    normalized_schema_version: int = Field(ge=1)
    dataset_column: str
    method_column: str
    workflow_column: str
    cell_type_column: str
    gene_count_column: str
    beta_column: str
    standardized_beta_column: str
    standard_error_column: str
    p_value_column: str
    adjusted_p_value_column_pattern: str
    significance_column_pattern: str
    delimiter: str
    null_value: str

    @field_validator(
        "dataset_column", "method_column", "workflow_column", "cell_type_column",
        "gene_count_column", "beta_column", "standardized_beta_column",
        "standard_error_column", "p_value_column", "null_value",
    )
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @field_validator(
        "adjusted_p_value_column_pattern", "significance_column_pattern",
    )
    @classmethod
    def correction_pattern(cls, value: str) -> str:
        try:
            rendered = value.format(method="bonferroni")
        except (KeyError, ValueError) as exc:
            raise ValueError("must support the {method} placeholder") from exc
        if "{method}" not in value or not rendered.strip():
            raise ValueError("must contain the {method} placeholder")
        return value

    @field_validator("delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @model_validator(mode="after")
    def distinct_output_columns(self):
        names = [
            self.dataset_column,
            self.method_column,
            self.workflow_column,
            self.cell_type_column,
            self.gene_count_column,
            self.beta_column,
            self.standardized_beta_column,
            self.standard_error_column,
            self.p_value_column,
        ]
        if len(names) != len(set(names)):
            raise ValueError("fixed result column names must be unique")
        return self


class SingleCellOutputLayout(StrictModel):
    engine_directory: str
    results_file: str
    log_file: str
    resolved_config_file: str
    completion_manifest: str
    staging_directory: str
    scdrs_engine_directory: str
    scdrs_qc_report: str
    scdrs_completion_manifest: str
    scdrs_staging_directory: str
    scdrs_score_file_pattern: str
    scdrs_full_score_file_pattern: str
    scdrs_group_file_pattern: str
    scdrs_correlation_file_pattern: str
    scdrs_gene_file_pattern: str
    scdrs_generated_gene_statistics: str
    scdrs_generated_gene_set: str
    scdrs_gene_mapping_report: str

    @field_validator("*")
    @classmethod
    def safe_relative_pattern(cls, value: str) -> str:
        return _safe_relative_pattern(value)

    @model_validator(mode="after")
    def required_dataset_tokens(self):
        trait_patterns = {
            "scdrs_score_file_pattern",
            "scdrs_full_score_file_pattern",
            "scdrs_group_file_pattern",
            "scdrs_correlation_file_pattern",
            "scdrs_gene_file_pattern",
        }
        annotation_patterns = {"scdrs_group_file_pattern"}
        for name, value in self.model_dump().items():
            if name in trait_patterns:
                if "{trait}" not in value:
                    raise ValueError("%s must contain '{trait}'" % name)
                if name in annotation_patterns and "{annotation}" not in value:
                    raise ValueError("%s must contain '{annotation}'" % name)
                continue
            if "{dataset_id}" not in value:
                raise ValueError("%s must contain '{dataset_id}'" % name)
        return self


class SingleCellConfig(ModuleConfig):
    tools: list[SingleCellTool] = Field(min_length=1)
    magma_celltype: MagmaCelltypeConfig
    scdrs: ScdrsConfig
    multiple_testing: SingleCellMultipleTestingConfig
    result_schema: SingleCellResultSchema
    output_layout: SingleCellOutputLayout

    @field_validator("tools")
    @classmethod
    def unique_tools(cls, values: list[SingleCellTool]) -> list[SingleCellTool]:
        if len(values) != len(set(values)):
            raise ValueError("must contain unique tool names")
        return values

    @model_validator(mode="after")
    def output_columns_do_not_collide(self):
        fixed = set(
            value
            for name, value in self.result_schema.model_dump().items()
            if name.endswith("_column")
        )
        generated = []
        for method in self.multiple_testing.methods:
            generated.extend(
                (
                    self.result_schema.adjusted_p_value_column_pattern.format(
                        method=method,
                    ),
                    self.result_schema.significance_column_pattern.format(
                        method=method,
                    ),
                )
            )
        if len(generated) != len(set(generated)) or fixed.intersection(generated):
            raise ValueError(
                "configured fixed and multiple-testing result columns must be unique"
            )
        return self


__all__ = [
    "SINGLE_CELL_TOOLS", "SingleCellConfig", "SingleCellCorrectionMethod",
    "SingleCellTool",
]
