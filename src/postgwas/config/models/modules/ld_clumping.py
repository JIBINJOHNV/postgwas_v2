"""Validated configuration for region and LD-reference clumping."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import (
    GenomeBuild,
    GenomicRegion,
    ModuleConfig,
    Population,
    StrictModel,
)
from postgwas.core.paths import validate_filename_component


LDClumpingMethod = Literal["region", "standard"]


class LDClumpingInputsConfig(StrictModel):
    vcf: Path | None = None
    dataset_id: str | None = None

    @field_validator("dataset_id")
    @classmethod
    def safe_dataset_id(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return validate_filename_component(
            value, "modules.ld_clumping.inputs.dataset_id"
        )


class LDClumpingReferenceConfig(StrictModel):
    directory: Path | None = None
    manifest_filename: str
    file_pattern: str
    index_suffix: str
    format_version: int = Field(ge=1)
    orientation: Literal["symmetric_first_endpoint"]
    columns: list[str]

    @field_validator("manifest_filename")
    @classmethod
    def safe_manifest_filename(cls, value: str) -> str:
        return validate_filename_component(
            value, "modules.ld_clumping.reference.manifest_filename"
        )

    @field_validator("file_pattern")
    @classmethod
    def valid_file_pattern(cls, value: str) -> str:
        value = value.strip()
        required = ("{population}", "{chromosome}")
        if any(value.count(token) != 1 for token in required):
            raise ValueError(
                "must contain {population} and {chromosome} exactly once"
            )
        try:
            rendered = Path(value.format(population="EUR", chromosome="1"))
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if rendered.is_absolute() or ".." in rendered.parts:
            raise ValueError("must stay inside the configured reference directory")
        if not str(rendered).endswith(".ld.gz"):
            raise ValueError("must render a .ld.gz file")
        return value

    @field_validator("index_suffix")
    @classmethod
    def valid_index_suffix(cls, value: str) -> str:
        if not value or Path(value).name != value:
            raise ValueError("must be a non-empty filename suffix")
        return value

    @field_validator("columns")
    @classmethod
    def valid_columns(cls, values: list[str]) -> list[str]:
        required = [
            "chromosome_a",
            "position_a",
            "variant_a",
            "chromosome_b",
            "position_b",
            "variant_b",
            "r2",
        ]
        if values != required:
            raise ValueError(
                "must describe the seven-column PLINK --r2 contract in order: %s"
                % ", ".join(required)
            )
        return values


class LDClumpingVcfFieldsConfig(StrictModel):
    chromosome: str
    position: str
    reference_allele: str
    alternate_allele: str
    variant_id: str
    effect: str
    standard_error: str
    allele_frequency: str
    log_pvalue: str
    ld_block: str

    @field_validator("*")
    @classmethod
    def nonempty_expression(cls, value: str) -> str:
        value = value.strip()
        if not value or "\t" in value or "\n" in value:
            raise ValueError("must be one non-empty bcftools query expression")
        return value

    @field_validator("ld_block")
    @classmethod
    def valid_ld_block_template(cls, value: str) -> str:
        if value.count("{population}") != 1:
            raise ValueError("must contain {population} exactly once")
        try:
            rendered = value.format(population="EUR")
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if not re.fullmatch(r"%INFO/[A-Za-z][A-Za-z0-9_.-]*", rendered):
            raise ValueError("must render one bcftools INFO-field expression")
        return value


class LDClumpingTableConfig(StrictModel):
    delimiter: str
    null_values: list[str]
    infer_schema_length: int = Field(ge=1)
    io_buffer_bytes: int = Field(ge=1)
    biallelic_include_expression: str

    @field_validator("delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @field_validator("null_values")
    @classmethod
    def unique_null_values(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique values")
        return values

    @field_validator("biallelic_include_expression")
    @classmethod
    def nonempty_include_expression(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("must not be empty")
        return value


class LDClumpingComputeConfig(StrictModel):
    minimum_worker_memory_gb: float = Field(gt=0, allow_inf_nan=False)
    input_memory_multiplier: float = Field(ge=1, allow_inf_nan=False)


class LDClumpingOutputLayoutConfig(StrictModel):
    region_working_table: str
    region_raw_table: str
    region_pruned_table: str
    region_significant_table: str
    region_log: str
    standard_working_table: str
    standard_formatted_table: str
    standard_log: str
    standard_summary: str
    standard_hierarchy: str
    standard_independent_clusters: str
    standard_lead_clusters: str
    canonical_log: str
    resolved_configuration: str

    @field_validator("*")
    @classmethod
    def safe_output_pattern(cls, value: str) -> str:
        value = value.strip()
        if "{dataset_id}" not in value:
            raise ValueError("must contain {dataset_id}")
        try:
            rendered = Path(value.format(dataset_id="dataset", population="EUR"))
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if rendered.is_absolute() or ".." in rendered.parts:
            raise ValueError("must stay inside the configured output directory")
        return value

    @model_validator(mode="after")
    def unique_destinations(self):
        values = list(self.model_dump().values())
        if len(values) != len(set(values)):
            raise ValueError("LD-clumping output path patterns must be unique")
        population_specific = (
            self.region_pruned_table,
            self.region_significant_table,
        )
        if any("{population}" not in value for value in population_specific):
            raise ValueError("region output tables must contain {population}")
        return self


class LDClumpingConfig(ModuleConfig):
    inputs: LDClumpingInputsConfig
    output_directory: Path | None = None
    methods: list[LDClumpingMethod]
    genome_build: GenomeBuild
    population: Population
    lead_pvalue: float = Field(gt=0, le=1, allow_inf_nan=False)
    candidate_pvalue: float = Field(gt=0, le=1, allow_inf_nan=False)
    clump_r2: float = Field(ge=0, le=1, allow_inf_nan=False)
    lead_r2: float = Field(ge=0, le=1, allow_inf_nan=False)
    window_kb: int = Field(gt=0)
    merge_distance_bp: int = Field(ge=0)
    missing_index_action: Literal["error", "self_only"]
    remove_mhc: bool
    mhc_regions: dict[GenomeBuild, GenomicRegion]
    summary_pvalue_thresholds: dict[str, float]
    reference: LDClumpingReferenceConfig
    vcf_fields: LDClumpingVcfFieldsConfig
    table: LDClumpingTableConfig
    compute: LDClumpingComputeConfig
    output_layout: LDClumpingOutputLayoutConfig

    @field_validator("methods")
    @classmethod
    def unique_methods(cls, values: list[LDClumpingMethod]) -> list[LDClumpingMethod]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique clumping methods")
        return values

    @field_validator("mhc_regions")
    @classmethod
    def configured_mhc_regions(
        cls, values: dict[GenomeBuild, GenomicRegion]
    ) -> dict[GenomeBuild, GenomicRegion]:
        if not values:
            raise ValueError("must define at least one build-specific MHC region")
        return values

    @field_validator("summary_pvalue_thresholds")
    @classmethod
    def valid_summary_thresholds(cls, values: dict[str, float]) -> dict[str, float]:
        if not values:
            raise ValueError("must contain at least one named threshold")
        for name, threshold in values.items():
            validate_filename_component(
                name, "modules.ld_clumping.summary_pvalue_thresholds key"
            )
            if not 0 < threshold <= 1:
                raise ValueError("all summary P-value thresholds must be in (0, 1]")
        return values

    @model_validator(mode="after")
    def scientifically_consistent_thresholds(self):
        if self.candidate_pvalue < self.lead_pvalue:
            raise ValueError(
                "candidate_pvalue must be greater than or equal to lead_pvalue"
            )
        if self.lead_r2 > self.clump_r2:
            raise ValueError("lead_r2 must be less than or equal to clump_r2")
        if self.remove_mhc and self.genome_build not in self.mhc_regions:
            raise ValueError(
                "mhc_regions must define the selected genome_build when MHC removal "
                "is enabled"
            )
        return self
