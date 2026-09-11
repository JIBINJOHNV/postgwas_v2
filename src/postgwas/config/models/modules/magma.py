"""Typed configuration for MAGMA gene and competitive gene-set analysis."""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from typing import Literal, get_args

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import (
    ChromosomeAnalysisScopeConfig,
    GenomeBuild,
    MHCAnalysisScopeConfig,
    MHCPolicy,
    MHC_POLICIES,
    ModuleConfig,
    PValueCorrectionMethod,
    Population,
    StrictModel,
)


MagmaDuplicatePolicy = Literal["err", "lowest_p", "remove"]
MAGMA_DUPLICATE_POLICIES = get_args(MagmaDuplicatePolicy)
MagmaMHCPolicy = MHCPolicy
MAGMA_MHC_POLICIES = MHC_POLICIES
MagmaCorrectionMethod = PValueCorrectionMethod
MagmaMappingMethod = Literal[
    "positional", "emagma", "h_magma", "n_magma", "chrom_magma",
]
MagmaGeneIdentifierType = Literal["entrez", "ensembl", "symbol", "mixed"]
MAGMA_GENE_IDENTIFIER_TYPES = get_args(MagmaGeneIdentifierType)
MagmaGeneSetInputFormat = Literal["auto", "gmt", "magma", "membership"]
MagmaGeneSetMismatchAction = Literal["skip", "error"]
MAGMA_GENE_SET_MISMATCH_ACTIONS = get_args(MagmaGeneSetMismatchAction)
MagmaAlternateGeneIdDuplicatePolicy = Literal["longest_interval", "error"]
MAGMA_ALTERNATE_GENE_ID_DUPLICATE_POLICIES = get_args(
    MagmaAlternateGeneIdDuplicatePolicy
)
MagmaResultStatisticType = Literal[
    "calibrated_gene_p_value", "minimum_regulatory_element_p_value",
]


def _safe_relative_pattern(value: str) -> str:
    value = value.strip()
    if not value:
        raise ValueError("must not be empty")
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError("must be a relative path below the output directory")
    return value


def _validated_delimiter_pattern(value: str) -> str:
    """Require a regex separator that consumes at least one character."""
    if not value:
        raise ValueError("must not be empty")
    try:
        compiled = re.compile(value)
    except re.error as exc:
        raise ValueError("must be a valid regular expression") from exc
    if compiled.match("") is not None:
        raise ValueError("must not match empty text")
    return value


class MagmaInputConfig(StrictModel):
    snp_location_file: str | None = None
    p_value_file: str | None = None
    ld_reference_prefix: str | None = None
    gene_location_file: str | None = None
    gene_set_file: str | None = None
    sample_size_column: str
    table_delimiter_pattern: str
    bim_delimiter: Literal["\t", " "]
    output_table_delimiter: Literal["\t", " "]
    chromosome_prefix_pattern: str
    chromosome_aliases: dict[str, str]
    invalid_chromosome_labels: list[str]
    invalid_allele_labels: list[str]
    variant_id_column: str
    chromosome_column: str
    position_column: str
    reference_allele_column: str
    alternate_allele_column: str
    p_value_column: str
    required_reference_extensions: list[str]
    bim_extension: str
    bim_columns: list[str]
    gene_location_has_header: bool
    gene_location_columns: list[str]
    alternate_gene_id_type: MagmaGeneIdentifierType

    @field_validator(
        "snp_location_file", "p_value_file", "ld_reference_prefix",
        "gene_location_file", "gene_set_file",
    )
    @classmethod
    def optional_nonempty_path(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value

    @field_validator(
        "sample_size_column", "variant_id_column", "chromosome_column",
        "position_column", "reference_allele_column", "alternate_allele_column",
        "p_value_column",
    )
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @field_validator("table_delimiter_pattern")
    @classmethod
    def valid_delimiter_pattern(cls, value: str) -> str:
        return _validated_delimiter_pattern(value)

    @field_validator("chromosome_prefix_pattern")
    @classmethod
    def valid_chromosome_prefix_pattern(cls, value: str) -> str:
        try:
            re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        return value

    @field_validator("chromosome_aliases")
    @classmethod
    def nonempty_chromosome_aliases(
        cls, values: dict[str, str],
    ) -> dict[str, str]:
        if any(
            not source.strip() or not target.strip()
            for source, target in values.items()
        ):
            raise ValueError("must not contain empty chromosome labels")
        normalized = {
            source.strip().upper(): target.strip().upper()
            for source, target in values.items()
        }
        if len(normalized) != len(values):
            raise ValueError("must remain unique after chromosome-label normalization")
        return normalized

    @field_validator("required_reference_extensions")
    @classmethod
    def unique_reference_extensions(cls, values: list[str]) -> list[str]:
        if (
            not values
            or len(values) != len(set(values))
            or any(not value.startswith(".") or Path(value).name != value for value in values)
        ):
            raise ValueError("must contain unique filename extensions beginning with '.'")
        return values

    @field_validator("invalid_chromosome_labels", "invalid_allele_labels")
    @classmethod
    def unique_invalid_labels(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique labels")
        return values

    @field_validator("bim_extension")
    @classmethod
    def valid_bim_extension(cls, value: str) -> str:
        if not value.startswith(".") or Path(value).name != value:
            raise ValueError("must be a filename extension beginning with '.'")
        return value

    @field_validator("bim_columns")
    @classmethod
    def complete_bim_roles(cls, values: list[str]) -> list[str]:
        required = {
            "chromosome", "variant_id", "genetic_distance", "position",
            "allele1", "allele2",
        }
        if len(values) != 6 or set(values) != required:
            raise ValueError(
                "must map the six PLINK BIM roles: %s"
                % ", ".join(sorted(required))
            )
        return values

    @field_validator("gene_location_columns")
    @classmethod
    def complete_gene_location_roles(cls, values: list[str]) -> list[str]:
        required = ["gene_id", "chromosome", "start", "end", "strand"]
        allowed = required + ["alternate_gene_id"]
        if values not in (required, allowed):
            raise ValueError(
                "must map, in order, gene_id, chromosome, start, end, strand, "
                "with optional alternate_gene_id as column six"
            )
        return values


class MagmaSnpHarmonisationConfig(StrictModel):
    resolve_variants_to_reference: bool
    minimum_overlap_fraction: float = Field(gt=0, le=1, allow_inf_nan=False)
    duplicate_policy: MagmaDuplicatePolicy


class MagmaMHCConfig(MHCAnalysisScopeConfig):
    """MAGMA specialization of the shared MHC analysis-scope schema."""


class MagmaChromosomeConfig(ChromosomeAnalysisScopeConfig):
    """MAGMA specialization of the shared chromosome-scope schema."""


class MagmaExclusionReportingConfig(StrictModel):
    chromosome_reason: str
    mhc_reason: str
    annotated_units_scope: str
    reason_column: str
    scope_column: str
    excluded_count_column: str
    input_count_column: str
    retained_count_column: str
    start_column: str
    end_column: str

    @field_validator("*")
    @classmethod
    def nonempty_values(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @model_validator(mode="after")
    def distinct_reasons(self):
        if self.chromosome_reason == self.mhc_reason:
            raise ValueError("chromosome_reason and mhc_reason must differ")
        return self


class MagmaGeneSetConfig(StrictModel):
    input_format: MagmaGeneSetInputFormat
    minimum_gene_id_overlap_fraction: float = Field(
        gt=0, le=1, allow_inf_nan=False,
    )
    identifier_mismatch_action: MagmaGeneSetMismatchAction
    alternate_id_duplicate_policy: MagmaAlternateGeneIdDuplicatePolicy
    membership_delimiter_pattern: str
    membership_has_header: bool
    membership_set_column: int = Field(ge=0)
    membership_gene_column: int = Field(ge=0)

    @field_validator("membership_delimiter_pattern")
    @classmethod
    def valid_membership_delimiter(cls, value: str) -> str:
        return _validated_delimiter_pattern(value)

    @model_validator(mode="after")
    def distinct_membership_columns(self):
        if self.membership_set_column == self.membership_gene_column:
            raise ValueError(
                "membership_set_column and membership_gene_column must differ"
            )
        return self


class MagmaAnnotationValidationConfig(StrictModel):
    minimum_bim_variant_overlap_fraction: float = Field(
        gt=0, le=1, allow_inf_nan=False,
    )
    coordinate_pattern: str
    comment_prefix: str
    invalid_variant_identifiers: list[str]

    @field_validator("coordinate_pattern")
    @classmethod
    def valid_coordinate_pattern(cls, value: str) -> str:
        try:
            compiled = re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        required = {"chromosome", "start", "end"}
        if set(compiled.groupindex) != required:
            raise ValueError(
                "must define chromosome, start, and end named groups"
            )
        return value

    @field_validator("comment_prefix")
    @classmethod
    def nonempty_comment_prefix(cls, value: str) -> str:
        if not value:
            raise ValueError("must not be empty")
        return value

    @field_validator("invalid_variant_identifiers")
    @classmethod
    def unique_invalid_variant_identifiers(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)) or any(
            not value for value in values
        ):
            raise ValueError("must contain unique non-empty identifiers")
        return values


class MagmaMappingDefinition(StrictModel):
    method: MagmaMappingMethod
    display_name: str
    genome_build: GenomeBuild
    population: Population
    gene_id_type: MagmaGeneIdentifierType
    context: str | None = None
    source_name: str
    source_version: str
    source_url: str
    result_statistic_type: MagmaResultStatisticType
    result_statistic_interpretation: str
    gene_location_file: str | None = None
    gene_annotation_file: str | None = None
    gene_annotation_files: list[str] = Field(default_factory=list)
    gene_set_file: str | None = None
    gene_set_format: MagmaGeneSetInputFormat = "auto"
    minimum_gene_id_overlap_fraction: float | None = Field(
        default=None, gt=0, le=1, allow_inf_nan=False,
    )
    gene_settings: list[str] = Field(default_factory=list)
    annotation_window_upstream_kb: int | None = Field(default=None, ge=0)
    annotation_window_downstream_kb: int | None = Field(default=None, ge=0)
    regulatory_element_location_file: str | None = None
    element_to_gene_file: str | None = None

    @field_validator(
        "display_name", "source_name", "source_version", "source_url",
        "result_statistic_interpretation",
    )
    @classmethod
    def nonempty_metadata(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value.strip()

    @field_validator("source_url")
    @classmethod
    def https_source(cls, value: str) -> str:
        if not value.startswith("https://"):
            raise ValueError("must use an https URL")
        return value

    @field_validator("gene_annotation_files")
    @classmethod
    def valid_annotation_files(cls, values: list[str]) -> list[str]:
        if len(values) != len(set(values)) or any(not value.strip() for value in values):
            raise ValueError("must contain unique non-empty paths")
        return values

    @field_validator("context")
    @classmethod
    def optional_context(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value.strip() if value is not None else None

    @field_validator(
        "gene_location_file", "gene_annotation_file", "gene_set_file",
        "regulatory_element_location_file", "element_to_gene_file",
    )
    @classmethod
    def optional_path(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value

    @field_validator("gene_settings")
    @classmethod
    def valid_gene_settings(cls, values: list[str]) -> list[str]:
        if len(values) != len(set(values)) or any(
            not value.strip() or value.lstrip().startswith("--") for value in values
        ):
            raise ValueError(
                "must contain unique non-empty MAGMA gene-setting values without flags"
            )
        return values

    @model_validator(mode="after")
    def method_specific_inputs(self):
        external = {"emagma", "h_magma"}
        if self.method in external and not self.gene_annotation_file:
            raise ValueError(
                "%s requires gene_annotation_file" % self.method
            )
        if self.method == "positional" and self.gene_annotation_file is not None:
            raise ValueError("positional mapping cannot use gene_annotation_file")
        if self.method == "n_magma" and (
            not self.gene_location_file or not self.gene_annotation_files
        ):
            raise ValueError(
                "n_magma requires gene_location_file and gene_annotation_files"
            )
        if self.method == "n_magma" and (
            self.annotation_window_upstream_kb is None
            or self.annotation_window_downstream_kb is None
        ):
            raise ValueError(
                "n_magma requires explicit annotation windows so its positional "
                "component cannot inherit conventional MAGMA windows"
            )
        if self.method != "n_magma" and self.gene_annotation_files:
            raise ValueError("gene_annotation_files is valid only for n_magma")
        if self.method == "chrom_magma":
            missing = [
                name
                for name in (
                    "regulatory_element_location_file", "element_to_gene_file",
                )
                if getattr(self, name) is None
            ]
            if missing:
                raise ValueError(
                    "chrom_magma requires %s" % ", ".join(missing)
                )
            if self.gene_set_file is not None:
                raise ValueError(
                    "chrom_magma does not support MAGMA competitive gene-set "
                    "testing because its tested units are regulatory elements"
                )
            if self.result_statistic_type != "minimum_regulatory_element_p_value":
                raise ValueError(
                    "chrom_magma must label its gene-level statistic as the "
                    "minimum regulatory-element p-value"
                )
            if (
                self.annotation_window_upstream_kb is None
                or self.annotation_window_downstream_kb is None
            ):
                raise ValueError(
                    "chrom_magma requires explicit annotation windows"
                )
        elif (
            self.regulatory_element_location_file is not None
            or self.element_to_gene_file is not None
        ):
            raise ValueError(
                "regulatory-element inputs are valid only for chrom_magma"
            )
        elif self.result_statistic_type != "calibrated_gene_p_value":
            raise ValueError(
                "%s must use calibrated_gene_p_value" % self.method
            )
        return self


class MagmaMappingConfig(StrictModel):
    selected: list[str]
    primary: str
    definitions: dict[str, MagmaMappingDefinition]

    @field_validator("selected")
    @classmethod
    def unique_selected(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique analysis names")
        return values

    @field_validator("definitions")
    @classmethod
    def safe_definition_names(
        cls, values: dict[str, MagmaMappingDefinition],
    ) -> dict[str, MagmaMappingDefinition]:
        if not values:
            raise ValueError("must contain at least one mapping definition")
        invalid = [
            name for name in values
            if not re.fullmatch(r"[a-z][a-z0-9_]*", name)
        ]
        if invalid:
            raise ValueError(
                "definition names must use lowercase snake_case: %s"
                % ", ".join(invalid)
            )
        return values

    @model_validator(mode="after")
    def selected_definitions_exist(self):
        unknown = sorted(set(self.selected) - set(self.definitions))
        if unknown:
            raise ValueError(
                "selected mapping definitions do not exist: %s"
                % ", ".join(unknown)
            )
        if self.primary not in self.selected:
            raise ValueError("primary must occur in selected")
        return self


class ChromMagmaMappingSchema(StrictModel):
    mapping_delimiter_pattern: str
    mapping_has_header: bool
    location_delimiter_pattern: str
    location_has_header: bool
    element_id_column: int = Field(ge=0)
    chromosome_column: int = Field(ge=0)
    start_column: int = Field(ge=0)
    end_column: int = Field(ge=0)
    gene_id_column: int = Field(ge=0)
    score_column: int | None = Field(default=None, ge=0)
    minimum_element_mapping_fraction: float = Field(
        gt=0, le=1, allow_inf_nan=False,
    )
    gene_assignment_policy: Literal["lowest_p"]
    result_element_column: str
    result_p_value_column: str
    report_element_column: str
    report_element_count_column: str
    report_source_score_column: str
    report_original_column_prefix: str
    tie_breaker: Literal["element_id"]

    @field_validator("mapping_delimiter_pattern", "location_delimiter_pattern")
    @classmethod
    def valid_delimiter(cls, value: str) -> str:
        return _validated_delimiter_pattern(value)

    @field_validator(
        "result_element_column", "result_p_value_column", "report_element_column",
        "report_element_count_column", "report_source_score_column",
        "report_original_column_prefix",
    )
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @model_validator(mode="after")
    def distinct_input_columns(self):
        required = [
            self.element_id_column, self.chromosome_column, self.start_column,
            self.end_column, self.gene_id_column,
        ]
        if len(required) != len(set(required)):
            raise ValueError("required element-to-gene columns must be distinct")
        return self


class MagmaBatchingConfig(StrictModel):
    enabled: bool
    minimum_genes_per_batch: int = Field(ge=1)
    memory_per_process_gb: float = Field(gt=0, allow_inf_nan=False)


class MagmaCorrectionFamily(StrictModel):
    pattern: str
    methods: list[MagmaCorrectionMethod]

    @field_validator("pattern")
    @classmethod
    def valid_pattern(cls, value: str) -> str:
        try:
            compiled = re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        if compiled.groups:
            raise ValueError(
                "must not contain capturing groups; use (?:...) for grouping"
            )
        return value

    @field_validator("methods")
    @classmethod
    def unique_methods(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique correction methods")
        return values


class MagmaMultipleTestingConfig(StrictModel):
    reporting_significance_threshold: float = Field(
        gt=0, le=1, allow_inf_nan=False,
    )
    primary_method: MagmaCorrectionMethod
    reporting_method_labels: dict[MagmaCorrectionMethod, str]
    gene_methods: list[MagmaCorrectionMethod]
    global_methods: list[MagmaCorrectionMethod]
    families: dict[str, MagmaCorrectionFamily]

    @field_validator("gene_methods", "global_methods")
    @classmethod
    def unique_methods(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique correction methods")
        return values

    @field_validator("families")
    @classmethod
    def safe_family_names(
        cls, values: dict[str, MagmaCorrectionFamily],
    ) -> dict[str, MagmaCorrectionFamily]:
        if any(not re.fullmatch(r"[a-z][a-z0-9_]*", name) for name in values):
            raise ValueError("family names must use lowercase snake_case")
        return values

    @field_validator("reporting_method_labels")
    @classmethod
    def nonempty_reporting_labels(cls, values: dict[str, str]) -> dict[str, str]:
        if not values or any(not label.strip() for label in values.values()):
            raise ValueError("must provide non-empty reporting labels")
        return values

    @model_validator(mode="after")
    def reporting_labels_cover_configured_methods(self):
        required = set(self.gene_methods) | set(self.global_methods)
        missing = sorted(required - set(self.reporting_method_labels))
        if missing:
            raise ValueError(
                "reporting_method_labels is missing configured methods: %s"
                % ", ".join(missing)
            )
        if self.primary_method not in self.global_methods:
            raise ValueError(
                "primary_method must occur in global_methods because the primary "
                "gene-set decision is calculated across all tested gene sets"
            )
        return self


class MagmaResultSchema(StrictModel):
    table_delimiter_pattern: str
    gene_set_name_column: str
    gene_set_p_value_column: str
    gene_set_full_name_column: str
    gene_id_column: str
    gene_p_value_column: str
    gene_reference_chromosome_column: str
    gene_reference_start_column: str
    gene_reference_end_column: str
    gene_reference_strand_column: str
    gene_reference_alternate_id_column: str
    global_correction_column_pattern: str
    family_correction_column_pattern: str
    primary_correction_method_column: str
    primary_adjusted_p_value_column: str
    primary_significant_column: str
    report_dataset_column: str
    report_gene_set_description_column: str
    report_source_input_genes_column: str
    report_input_genes_column: str
    report_common_genes_column: str
    report_common_gene_p_values_column: str
    report_total_genes_column: str
    report_common_gene_count_column: str
    report_untested_genes_column: str
    report_untested_gene_count_column: str
    report_mapping_name_column: str
    report_mapping_method_column: str
    report_mapping_context_column: str
    report_gene_id_type_column: str
    report_source_name_column: str
    report_source_version_column: str
    report_statistic_type_column: str
    report_statistic_interpretation_column: str
    report_mhc_policy_column: str
    report_mhc_region_column: str
    report_excluded_chromosomes_column: str
    report_excluded_units_column: str
    report_gene_set_status_column: str
    report_gene_set_reason_column: str
    report_delimiter: str
    report_null_value: str

    @field_validator("table_delimiter_pattern")
    @classmethod
    def valid_table_delimiter_pattern(cls, value: str) -> str:
        return _validated_delimiter_pattern(value)

    @field_validator(
        "gene_set_name_column", "gene_set_p_value_column",
        "gene_set_full_name_column", "gene_id_column", "gene_p_value_column",
        "gene_reference_chromosome_column", "gene_reference_start_column",
        "gene_reference_end_column", "gene_reference_strand_column",
        "gene_reference_alternate_id_column",
        "primary_correction_method_column", "primary_adjusted_p_value_column",
        "primary_significant_column",
        "report_dataset_column", "report_gene_set_description_column",
        "report_source_input_genes_column",
        "report_input_genes_column", "report_common_genes_column",
        "report_common_gene_p_values_column", "report_total_genes_column",
        "report_common_gene_count_column", "report_null_value",
        "report_untested_genes_column", "report_untested_gene_count_column",
        "report_mapping_name_column", "report_mapping_method_column",
        "report_mapping_context_column", "report_gene_id_type_column",
        "report_source_name_column", "report_source_version_column",
        "report_statistic_type_column", "report_statistic_interpretation_column",
        "report_mhc_policy_column", "report_mhc_region_column",
        "report_excluded_chromosomes_column", "report_excluded_units_column",
        "report_gene_set_status_column", "report_gene_set_reason_column",
    )
    @classmethod
    def nonempty_name(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @model_validator(mode="after")
    def distinct_primary_result_columns(self):
        columns = [
            self.primary_correction_method_column,
            self.primary_adjusted_p_value_column,
            self.primary_significant_column,
        ]
        if len(columns) != len(set(columns)):
            raise ValueError("primary result columns must be distinct")
        non_column_fields = {
            "table_delimiter_pattern",
            "global_correction_column_pattern",
            "family_correction_column_pattern",
            "primary_correction_method_column",
            "primary_adjusted_p_value_column",
            "primary_significant_column",
            "report_delimiter",
            "report_null_value",
        }
        other_columns = {
            value
            for name, value in self.model_dump().items()
            if name not in non_column_fields
        }
        collisions = sorted(set(columns) & other_columns)
        if collisions:
            raise ValueError(
                "primary result columns conflict with other configured result "
                "columns: %s" % ", ".join(collisions)
            )
        return self

    @field_validator("global_correction_column_pattern")
    @classmethod
    def valid_global_pattern(cls, value: str) -> str:
        try:
            rendered = value.format(method="bonferroni")
        except (KeyError, ValueError) as exc:
            raise ValueError("must support the {method} placeholder") from exc
        if "{method}" not in value or not rendered.strip():
            raise ValueError("must contain the {method} placeholder")
        return value

    @field_validator("family_correction_column_pattern")
    @classmethod
    def valid_family_pattern(cls, value: str) -> str:
        try:
            rendered = value.format(family="family", method="bonferroni")
        except (KeyError, ValueError) as exc:
            raise ValueError("must support {family} and {method} placeholders") from exc
        if "{family}" not in value or "{method}" not in value or not rendered.strip():
            raise ValueError("must contain {family} and {method} placeholders")
        return value

    @field_validator("report_delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

class MagmaOutputLayout(StrictModel):
    log_file: str
    resolved_config_file: str
    completion_manifest: str
    staging_directory: str
    harmonised_p_values: str
    harmonised_snp_locations: str
    excluded_variants: str
    annotation_prefix: str
    scoped_annotation: str
    component_annotation_prefix: str
    gene_result_prefix: str
    gene_batch_prefix: str
    gene_set_prefix: str
    corrected_genes: str
    corrected_gene_sets: str
    annotated_gene_sets: str
    pathway_compatible_gene_locations: str
    chrom_magma_genes: str
    excluded_genes: str
    exclusion_summary: str
    mapping_comparison: str
    pipeline_summary_csv: str
    pipeline_summary_html: str

    @field_validator("*")
    @classmethod
    def safe_relative_pattern(cls, value: str) -> str:
        return _safe_relative_pattern(value)

    @model_validator(mode="after")
    def required_dataset_tokens(self):
        mapping_fields = {
            "annotation_prefix",
            "scoped_annotation",
            "component_annotation_prefix",
            "gene_result_prefix",
            "gene_batch_prefix",
            "gene_set_prefix",
            "corrected_genes",
            "corrected_gene_sets",
            "annotated_gene_sets",
            "pathway_compatible_gene_locations",
            "chrom_magma_genes",
            "excluded_genes",
            "exclusion_summary",
        }
        for field, value in self.model_dump().items():
            if field == "resolved_config_file":
                continue
            if "{dataset_id}" not in value:
                raise ValueError("%s must contain '{dataset_id}'" % field)
            if field in mapping_fields and "{analysis}" not in value:
                raise ValueError("%s must contain '{analysis}'" % field)
        return self


class MagmaPipelineSummarySchema(StrictModel):
    """Configured machine-readable columns for the ordered pipeline summary."""

    version: int = Field(gt=0)
    step_column: str
    stage_column: str
    status_column: str
    summary_column: str
    details_column: str
    output_column: str
    csv_delimiter: str
    null_value: str

    @field_validator(
        "step_column", "stage_column", "status_column", "summary_column",
        "details_column", "output_column",
    )
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @field_validator("csv_delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @model_validator(mode="after")
    def unique_columns(self):
        columns = [
            self.step_column,
            self.stage_column,
            self.status_column,
            self.summary_column,
            self.details_column,
            self.output_column,
        ]
        duplicates = sorted(
            name for name, count in Counter(columns).items() if count > 1
        )
        if duplicates:
            raise ValueError(
                "pipeline summary columns must be unique: %s"
                % ", ".join(duplicates)
            )
        return self


class MagmaHtmlReportConfig(StrictModel):
    """Pagination and columns for the complete scientific HTML result tables."""

    page_size: int = Field(gt=0)
    gene_columns: list[str]
    pathway_columns: list[str]

    @field_validator("gene_columns", "pathway_columns")
    @classmethod
    def unique_nonempty_columns(cls, values: list[str]) -> list[str]:
        if not values or any(not value.strip() for value in values):
            raise ValueError("must contain one or more non-empty column names")
        if len(values) != len(set(values)):
            raise ValueError("must contain unique column names")
        return values


class MagmaScreenSummaryConfig(StrictModel):
    """Presentation settings for ranked associations shown after a MAGMA run."""

    top_gene_rows: int = Field(gt=0)
    top_pathway_rows: int = Field(gt=0)
    p_value_significant_digits: int = Field(ge=2, le=10)


class MagmaConfig(ModuleConfig):
    genome_build: GenomeBuild
    population: Population
    gene_window_upstream_kb: int = Field(ge=0)
    gene_window_downstream_kb: int = Field(ge=0)
    gene_model: str
    input: MagmaInputConfig
    snp_harmonisation: MagmaSnpHarmonisationConfig
    mhc: MagmaMHCConfig
    chromosomes: MagmaChromosomeConfig
    exclusion_reporting: MagmaExclusionReportingConfig
    gene_sets: MagmaGeneSetConfig
    annotation_validation: MagmaAnnotationValidationConfig
    mapping: MagmaMappingConfig
    chrom_magma_mapping: ChromMagmaMappingSchema
    batching: MagmaBatchingConfig
    multiple_testing: MagmaMultipleTestingConfig
    result_schema: MagmaResultSchema
    pipeline_summary_schema: MagmaPipelineSummarySchema
    html_report: MagmaHtmlReportConfig
    screen_summary: MagmaScreenSummaryConfig
    minimum_magma_version: str
    version_arguments: list[str]
    version_pattern: str
    output_layout: MagmaOutputLayout

    @field_validator("minimum_magma_version")
    @classmethod
    def semantic_version(cls, value: str) -> str:
        if not re.fullmatch(r"\d+(?:\.\d+){1,2}", value):
            raise ValueError("must be a dotted numeric version")
        return value

    @field_validator("gene_model")
    @classmethod
    def p_value_compatible_gene_model(cls, value: str) -> str:
        fixed = {
            "snp-wise", "snp-wise=mean", "snp-wise=top",
            "multi=snp-wise", "snp-wise=multi", "multi",
        }
        if value in fixed:
            return value
        match = re.fullmatch(
            r"snp-wise=top,((?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)",
            value,
        )
        if match is not None and float(match.group(1)) > 0:
            return value
        raise ValueError(
            "must be a MAGMA model supported with SNP p-value input: "
            "snp-wise[=mean|=top[,P]], multi=snp-wise, snp-wise=multi, or multi"
        )

    @field_validator("version_arguments")
    @classmethod
    def safe_version_arguments(cls, values: list[str]) -> list[str]:
        if any(not value.strip() for value in values):
            raise ValueError("must not contain empty arguments")
        return values

    @field_validator("version_pattern")
    @classmethod
    def compilable_version_pattern(cls, value: str) -> str:
        try:
            compiled = re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        if compiled.groups < 1:
            raise ValueError("must capture the numeric version in group 1")
        return value

    @model_validator(mode="after")
    def cross_field_contracts(self):
        if self.input.bim_extension not in self.input.required_reference_extensions:
            raise ValueError(
                "input.bim_extension must occur in input.required_reference_extensions"
            )
        for name in self.mapping.selected:
            definition = self.mapping.definitions[name]
            if definition.genome_build != self.genome_build:
                raise ValueError(
                    "mapping definition %s uses %s but modules.magma.genome_build "
                    "is %s"
                    % (name, definition.genome_build.value, self.genome_build.value)
                )
            if definition.population != self.population:
                raise ValueError(
                    "mapping definition %s uses %s but modules.magma.population is %s"
                    % (name, definition.population.value, self.population.value)
                )
        if (
            self.mhc.region_override is not None
            and self.mhc.region_override.start < 1
        ):
            raise ValueError("mhc.region_override.start must be one-based")
        correction_columns = [
            self.result_schema.global_correction_column_pattern.format(method=method)
            for method in self.multiple_testing.global_methods
        ]
        correction_columns.extend(
            self.result_schema.family_correction_column_pattern.format(
                family=family_name, method=method,
            )
            for family_name, family in self.multiple_testing.families.items()
            for method in family.methods
        )
        correction_columns.extend((
            self.result_schema.primary_correction_method_column,
            self.result_schema.primary_adjusted_p_value_column,
            self.result_schema.primary_significant_column,
        ))
        counts = Counter(correction_columns)
        duplicates = sorted(
            column for column, count in counts.items() if count > 1
        )
        if duplicates:
            raise ValueError(
                "configured correction result columns must be unique: %s"
                % ", ".join(duplicates)
            )
        source_columns = {
            self.result_schema.gene_set_name_column,
            self.result_schema.gene_set_p_value_column,
            self.result_schema.gene_set_full_name_column,
            self.result_schema.gene_id_column,
            self.result_schema.gene_p_value_column,
        }
        collisions = sorted(source_columns & set(correction_columns))
        if collisions:
            raise ValueError(
                "configured correction result columns conflict with MAGMA source "
                "columns: %s" % ", ".join(collisions)
            )
        return self


__all__ = [
    "MAGMA_DUPLICATE_POLICIES", "MAGMA_GENE_SET_MISMATCH_ACTIONS",
    "MAGMA_MHC_POLICIES", "MagmaConfig",
    "MagmaCorrectionMethod", "MagmaDuplicatePolicy", "MagmaMappingMethod",
]
