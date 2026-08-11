"""Typed values controlling harmonisation after sample-sheet normalization."""

from __future__ import annotations

import re
from typing import Any, Literal

from pydantic import Field, RootModel, field_validator, model_validator

from postgwas.config.models.common import ModuleConfig, StrictModel


class FrequencyReferenceConfig(StrictModel):
    source: str
    column: str
    available_sources: list[str]
    resource_examples: dict[str, str]

    @field_validator("source", "column")
    @classmethod
    def nonempty(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("must not be empty")
        return value

    @model_validator(mode="after")
    def supported_source(self):
        if (
            not self.available_sources
            or len(self.available_sources) != len(set(self.available_sources))
        ):
            raise ValueError("available_sources must contain unique values")
        if self.source not in self.available_sources:
            raise ValueError("source must be one of available_sources")
        if set(self.resource_examples) != set(self.available_sources):
            raise ValueError("resource_examples must define every available source")
        if any(
            not key.strip() or not value.strip()
            for key, value in self.resource_examples.items()
        ):
            raise ValueError("available sources and resource examples must not be empty")
        return self


class ComparisonAFConfig(FrequencyReferenceConfig):
    """Frequency reference stored as an indexed annotation VCF."""


class HarmonisationReferenceConfig(StrictModel):
    dbsnp_source: str


class VariantReferenceMappingConfig(StrictModel):
    """Structural columns and delimiter for a variant reference table."""

    chromosome: str
    position: str
    effect_allele: str
    other_allele: str
    delimiter: Literal["auto", "tab", "comma", "semicolon", "space", "whitespace", "pipe"]

    @field_validator("chromosome", "position", "effect_allele", "other_allele")
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("reference column names must not be empty")
        return value

    @model_validator(mode="after")
    def unique_columns(self):
        columns = [
            self.chromosome,
            self.position,
            self.effect_allele,
            self.other_allele,
        ]
        if len(columns) != len(set(columns)):
            raise ValueError("reference structural columns must be unique")
        return self


class FixedInfoConfig(StrictModel):
    """CLI-only constant INFO fallback and its collision-safe working column."""

    value: float | None = Field(default=None, ge=0, le=1)
    column: str

    @field_validator("column")
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("fixed INFO working column must not be empty")
        return value


class SampleSheetColumnAliases(StrictModel):
    """Configured input-header aliases used only to draft sample sheets."""

    chromosome_column: list[str]
    position_column: list[str]
    chromosome_position_column: list[str]
    variant_id_column: list[str]
    effect_allele_column: list[str]
    other_allele_column: list[str]
    effect_allele_frequency_column: list[str]
    maf_frequency_column: list[str]
    ambiguous_frequency_column: list[str]
    reference_frequency_column: list[str]
    beta_effect_column: list[str]
    odds_ratio_effect_column: list[str]
    ambiguous_effect_column: list[str]
    standard_error_column: list[str]
    z_score_column: list[str]
    p_value_column: list[str]
    control_count_column: list[str]
    case_count_column: list[str]
    total_sample_size_column: list[str]
    effective_sample_size_column: list[str]
    invalid_sample_size_column: list[str]
    imputation_info_column: list[str]

    @field_validator("*")
    @classmethod
    def nonempty_unique_aliases(cls, values: list[str]) -> list[str]:
        cleaned = [str(value).strip() for value in values]
        normalized = [re.sub(r"[^A-Z0-9]+", "_", value.upper()).strip("_") for value in cleaned]
        if not cleaned or any(not value for value in cleaned):
            raise ValueError("column-alias lists must not be empty")
        if len(normalized) != len(set(normalized)):
            raise ValueError("column aliases must be unique after normalization")
        return cleaned

    @model_validator(mode="after")
    def aliases_do_not_cross_scientific_fields(self):
        owners: dict[str, str] = {}
        for field, values in self.model_dump().items():
            for value in values:
                normalized = re.sub(r"[^A-Z0-9]+", "_", value.upper()).strip("_")
                previous = owners.get(normalized)
                if previous is not None:
                    raise ValueError(
                        "column alias %r is assigned to both %s and %s"
                        % (value, previous, field)
                    )
                owners[normalized] = field
        return self


class SampleSheetAlternativeFrequencyAliases(StrictModel):
    """ALT-labelled allele and frequency headers that must remain aligned."""

    allele_columns: list[str]
    frequency_columns: list[str]

    @field_validator("*")
    @classmethod
    def nonempty_unique_aliases(cls, values: list[str]) -> list[str]:
        cleaned = [str(value).strip() for value in values]
        normalized = [
            re.sub(r"[^A-Z0-9]+", "_", value.upper()).strip("_")
            for value in cleaned
        ]
        if not cleaned or any(not value for value in cleaned):
            raise ValueError("alternative-frequency alias lists must not be empty")
        if len(normalized) != len(set(normalized)):
            raise ValueError(
                "alternative-frequency aliases must be unique after normalization"
            )
        return cleaned


class SampleSheetGeneratorConfig(StrictModel):
    """File discovery and header semantics for v2 sample-sheet generation."""

    supported_suffixes: list[str]
    trait_type: Literal["auto"]
    effect_type: Literal["auto"]
    p_value_type: Literal["auto"]
    on_missing_sample_size: Literal["write_draft", "fail"]
    effect_frequency_prefixes: list[str]
    alternative_frequency: SampleSheetAlternativeFrequencyAliases
    column_aliases: SampleSheetColumnAliases

    @field_validator("supported_suffixes")
    @classmethod
    def valid_suffixes(cls, values: list[str]) -> list[str]:
        cleaned = [str(value).strip().lower() for value in values]
        if not cleaned or any(not value.startswith(".") for value in cleaned):
            raise ValueError("supported_suffixes must contain non-empty dot-prefixed suffixes")
        if len(cleaned) != len(set(cleaned)):
            raise ValueError("supported_suffixes must be unique")
        return cleaned

    @field_validator("effect_frequency_prefixes")
    @classmethod
    def valid_frequency_prefixes(cls, values: list[str]) -> list[str]:
        cleaned = [
            re.sub(r"[^A-Z0-9]+", "_", str(value).upper()).strip("_")
            for value in values
        ]
        if not cleaned or any(not value for value in cleaned):
            raise ValueError("effect_frequency_prefixes must not be empty")
        if len(cleaned) != len(set(cleaned)):
            raise ValueError("effect_frequency_prefixes must be unique")
        return cleaned

    @model_validator(mode="after")
    def alternative_frequency_aliases_are_consistent(self):
        def normalize(value: str) -> str:
            return re.sub(
                r"[^A-Z0-9]+", "_", str(value).upper()
            ).strip("_")

        effect_alleles = {
            normalize(value) for value in self.column_aliases.effect_allele_column
        }
        alternative_alleles = {
            normalize(value) for value in self.alternative_frequency.allele_columns
        }
        unknown_alleles = sorted(alternative_alleles - effect_alleles)
        if unknown_alleles:
            raise ValueError(
                "alternative-frequency allele aliases must also be configured as "
                "effect-allele aliases: %s" % ", ".join(unknown_alleles)
            )

        ordinary_aliases = {
            normalize(value)
            for values in self.column_aliases.model_dump().values()
            for value in values
        }
        alternative_frequencies = {
            normalize(value) for value in self.alternative_frequency.frequency_columns
        }
        overlap = sorted(alternative_frequencies & ordinary_aliases)
        if overlap:
            raise ValueError(
                "alternative-frequency aliases must not also appear in ordinary "
                "column aliases: %s" % ", ".join(overlap)
            )
        return self


class HarmonisationResourceLayout(StrictModel):
    """Relative resource paths, formatted below the configured resource root."""

    default_eaf: str
    comparison_af: str
    dbsnp: str
    fasta: str
    annotation: str
    chain: str
    build_check: str

    @field_validator(
        "default_eaf", "comparison_af", "dbsnp", "fasta", "annotation", "chain", "build_check",
    )
    @classmethod
    def valid_relative_template(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("must not be empty")
        if value.startswith(("/", "~")) or ".." in value.split("/"):
            raise ValueError("must be a relative path below the resource directory")
        return value


class HarmonisationOutputLayout(RootModel[dict[str, str]]):
    """Canonical dataset root and relative harmonisation result paths."""

    @model_validator(mode="after")
    def required_output_paths(self):
        required = {
            "dataset_directory", "logs_directory", "rejected_directory",
            "frequency_qc_directory", "qc_directory", "concordance_directory",
            "adapter_log_directory", "chromosome_table", "chromosome_source_snapshot",
            "chromosome_log",
            "chromosome_reject", "input_reject", "duplicates", "adapter_input",
            "adapter_mapping", "adapter_summary", "adapter_transcript",
            "adapter_output_vcf",
            "bcftools_transcript", "chromosome_raw_vcf", "chromosome_original_vcf",
            "chromosome_normalized_vcf", "chromosome_id_vcf",
            "chromosome_frequency_vcf", "chromosome_annotated_vcf",
            "chromosome_lifted_vcf", "chromosome_not_lifted_vcf",
            "sort_temp_directory", "merged_build_vcf", "merged_raw_vcf",
            "merged_not_lifted_vcf", "dataset_reject", "reject_summary",
            "run_manifest", "qc_summary", "gwas2vcf_summary", "dataset_log",
            "qc_assessment_summary", "qc_filter_rules", "qc_assessment_json",
            "qc_assessment_temporary",
            "population_frequency_qc", "population_frequency_temporary",
            "screen_report",
            "combined_log", "concordance_log", "concordance_summary",
            "concordance_mismatches", "concordance_input_only",
            "concordance_vcf_only", "concordance_vcf_duplicates",
            "concordance_position_matches", "concordance_all_matches",
            "missing_eaf", "out_of_range_eaf",
            "concordance_temporary",
            "adapter_merged_mapping", "adapter_input_archive",
        }
        missing = sorted(required - set(self.root))
        if missing:
            raise ValueError(
                "output_layout is missing paths: %s" % ", ".join(missing)
            )
        if any(not key.strip() or not value.strip() for key, value in self.root.items()):
            raise ValueError("output_layout keys and patterns must not be empty")
        return self


class Gwas2VcfInputConfig(StrictModel):
    required_column_keys: list[str]
    optional_column_keys: list[str]
    audit_columns: list[str]
    renamed_keys: dict[str, str]
    delimiter: str
    header: bool

    @field_validator(
        "required_column_keys", "optional_column_keys", "audit_columns",
    )
    @classmethod
    def unique_nonempty_entries(cls, values: list[str]) -> list[str]:
        cleaned = [str(value).strip() for value in values]
        if (
            not cleaned
            or any(not value for value in cleaned)
            or len(cleaned) != len(set(cleaned))
        ):
            raise ValueError("must contain unique non-empty values")
        return cleaned

    @model_validator(mode="after")
    def consistent_mapping(self):
        overlap = set(self.required_column_keys) & set(self.optional_column_keys)
        if overlap:
            raise ValueError(
                "required and optional column keys overlap: %s"
                % ", ".join(sorted(overlap))
            )
        known = set(self.required_column_keys) | set(self.optional_column_keys)
        unknown = sorted(set(self.renamed_keys) - known)
        if unknown:
            raise ValueError(
                "renamed_keys contains unknown column keys: %s" % ", ".join(unknown)
            )
        if not self.delimiter:
            raise ValueError("delimiter must not be empty")
        return self


class HarmonisationQCVcfFields(StrictModel):
    chromosome: str
    position: str
    reference_allele: str
    alternate_allele: str
    study_info_af: str
    external_info_af: str
    study_format_af: str
    imputation_format: str
    log_pvalue_format: str
    effective_sample_size_format: str

    @field_validator("*")
    @classmethod
    def nonempty_query_field(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("VCF query fields must not be empty")
        return value

    @model_validator(mode="after")
    def external_field_is_parameterised(self):
        if "{external_af}" not in self.external_info_af:
            raise ValueError("external_info_af must contain {external_af}")
        return self


class HarmonisationVcfConfig(StrictModel):
    """Tool-facing VCF fields and build transitions used by bcftools."""

    external_frequency_columns: list[str]
    missing_id_format: str
    target_builds: dict[str, str]
    required_merge_groups: list[str]
    liftover_plugin: str
    qc_fields: HarmonisationQCVcfFields
    concordance_fields: dict[str, str]
    table_delimiter: str
    table_null_values: list[str]
    table_null_output: str
    temporary_table_suffix: str
    io_buffer_bytes: int = Field(gt=0)

    @field_validator("external_frequency_columns")
    @classmethod
    def unique_frequency_columns(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain unique VCF columns and must not be empty")
        return values

    @field_validator("required_merge_groups")
    @classmethod
    def complete_required_merge_groups(cls, values: list[str]) -> list[str]:
        required = {"input_build", "target_build", "raw_gwas2vcf"}
        if len(values) != len(set(values)) or set(values) != required:
            raise ValueError(
                "required_merge_groups must contain exactly: %s"
                % ", ".join(sorted(required))
            )
        return values

    @field_validator("concordance_fields")
    @classmethod
    def complete_concordance_fields(cls, values: dict[str, str]) -> dict[str, str]:
        required = {"CHROM", "POS", "ID", "REF", "ALT", "ES", "SE", "EZ", "AF"}
        if set(values) != required or any(not str(value).strip() for value in values.values()):
            raise ValueError(
                "concordance_fields must define non-empty queries for %s"
                % ", ".join(sorted(required))
            )
        return values

    @field_validator("table_delimiter")
    @classmethod
    def single_table_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("table_delimiter must be exactly one character")
        return value

    @field_validator("table_null_values")
    @classmethod
    def unique_null_values(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("table_null_values must contain unique values")
        return values

    @field_validator("temporary_table_suffix")
    @classmethod
    def nonempty_temporary_suffix(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("temporary_table_suffix must not be empty")
        return value

    @field_validator("liftover_plugin")
    @classmethod
    def nonempty_plugin(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("must not be empty")
        return value

    @field_validator("missing_id_format")
    @classmethod
    def valid_missing_id_format(cls, value: str) -> str:
        value = value.strip()
        required = ("%CHROM", "%POS", "%REF", "%ALT")
        if not value.startswith("+") or any(field not in value for field in required):
            raise ValueError(
                "must start with '+' and contain %CHROM, %POS, %REF and %ALT"
            )
        return value

    @model_validator(mode="after")
    def complete_build_map(self):
        if len(self.target_builds) < 2:
            raise ValueError("target_builds must define at least two genome builds")
        if any(not source.strip() or not target.strip() for source, target in self.target_builds.items()):
            raise ValueError("target_builds names must not be empty")
        if any(source == target for source, target in self.target_builds.items()):
            raise ValueError("target_builds must map each build to a different build")
        if any(
            target not in self.target_builds
            or self.target_builds[target] != source
            for source, target in self.target_builds.items()
        ):
            raise ValueError("target_builds must contain reversible build pairs")
        return self


class PopulationFrequencyQCConfig(StrictModel):
    """Dataset-level population-frequency similarity on the raw merged VCF."""

    enabled: bool
    study_field: str
    population_fields: dict[str, str]
    minimum_comparable_variants: int = Field(ge=2)
    minimum_correlation: float = Field(ge=-1, le=1)
    minimum_correlation_gap: float = Field(ge=0, le=2)
    require_mae_agreement: bool
    warn_on_filename_mismatch: bool
    on_error: Literal["warn", "fail"]

    @field_validator("study_field")
    @classmethod
    def valid_study_info_field(cls, value: str) -> str:
        value = value.strip()
        if not re.fullmatch(r"%INFO/[A-Za-z][A-Za-z0-9_.-]*", value):
            raise ValueError("must be a bcftools INFO query such as %INFO/AF")
        return value

    @field_validator("population_fields")
    @classmethod
    def valid_population_fields(cls, values: dict[str, str]) -> dict[str, str]:
        if len(values) < 2:
            raise ValueError("must define at least two reference populations")
        if len(values.values()) != len(set(values.values())):
            raise ValueError("population VCF query fields must be unique")
        for population, query in values.items():
            if not re.fullmatch(r"[A-Z][A-Z0-9]*", population):
                raise ValueError(
                    "population labels must be uppercase alphanumeric tokens"
                )
            if not re.fullmatch(r"%INFO/[A-Za-z][A-Za-z0-9_.-]*", query.strip()):
                raise ValueError(
                    "population fields must be bcftools INFO queries such as %INFO/EUR"
                )
        return {population: query.strip() for population, query in values.items()}

    @model_validator(mode="after")
    def study_field_is_distinct(self):
        if self.study_field in self.population_fields.values():
            raise ValueError("study_field must differ from every population field")
        return self


class HarmonisationRuntimeConfig(StrictModel):
    display_screen: bool
    metadata_directory: str
    top_metadata_directory: str
    resolved_config_file: str
    run_summary_file: str
    sample_sheet_row_file: str
    command_file: str
    supplied_config_file: str

    @field_validator(
        "metadata_directory", "top_metadata_directory", "resolved_config_file",
        "run_summary_file",
        "sample_sheet_row_file", "command_file", "supplied_config_file",
    )
    @classmethod
    def relative_nonempty_path(cls, value: str) -> str:
        value = value.strip()
        if not value or value.startswith(("/", "~")) or ".." in value.split("/"):
            raise ValueError("must be a non-empty relative path")
        return value


class ConcordanceToleranceConfig(StrictModel):
    absolute: float = Field(ge=0)
    relative: float = Field(ge=0)


class ConcordanceFailureConfig(StrictModel):
    maximum_value_mismatch_fraction: float = Field(ge=0, le=1)
    maximum_vcf_duplicate_records: int = Field(ge=0)
    maximum_invalid_vcf_records: int = Field(ge=0)


class ConcordanceValidationConfig(StrictModel):
    enabled: bool
    allow_strand_complement: bool
    palindromic_action: Literal[
        "exclude", "compare_resolved", "compare_as_listed",
    ]
    write_all_matches: bool
    effect: ConcordanceToleranceConfig
    standard_error: ConcordanceToleranceConfig
    allele_frequency: ConcordanceToleranceConfig
    z_score: ConcordanceToleranceConfig
    failure: ConcordanceFailureConfig


class HarmonisationConfig(ModuleConfig):
    comparison_af: ComparisonAFConfig
    default_eaf: FrequencyReferenceConfig
    default_eaf_mapping: VariantReferenceMappingConfig
    reference: HarmonisationReferenceConfig
    external_eaf_mapping: VariantReferenceMappingConfig
    external_info_mapping: VariantReferenceMappingConfig
    fixed_info: FixedInfoConfig
    sample_sheet_generator: SampleSheetGeneratorConfig
    build_check_mapping: VariantReferenceMappingConfig
    resource_layout: HarmonisationResourceLayout
    output_layout: HarmonisationOutputLayout
    gwas2vcf_input: Gwas2VcfInputConfig
    vcf_processing: HarmonisationVcfConfig
    population_frequency_qc: PopulationFrequencyQCConfig
    runtime: HarmonisationRuntimeConfig
    concordance_validation: ConcordanceValidationConfig
    policies: dict[str, Any]
