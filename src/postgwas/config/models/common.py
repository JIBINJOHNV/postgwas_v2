"""Strict configuration primitives shared by scientific modules."""

from enum import Enum
from pathlib import Path
from typing import get_args, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from postgwas.core.io.delimiters import NAMED_DELIMITERS


PValueCorrectionMethod = Literal["bonferroni", "sidak", "holm", "fdr_bh"]
MHCPolicy = Literal[
    "include", "exclude_snps", "exclude_genes", "exclude_both",
]
MHC_POLICIES = get_args(MHCPolicy)


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)


class MultipleTestingSelectionConfig(StrictModel):
    """Shared selection and primary-decision contract for adjusted p-values."""

    methods: list[PValueCorrectionMethod] = Field(min_length=1)
    primary_method: PValueCorrectionMethod
    significance_threshold: float = Field(gt=0, le=1, allow_inf_nan=False)

    @field_validator("methods")
    @classmethod
    def unique_methods(
        cls, values: list[PValueCorrectionMethod],
    ) -> list[PValueCorrectionMethod]:
        if len(values) != len(set(values)):
            raise ValueError("must contain unique correction methods")
        return values

    @model_validator(mode="after")
    def primary_is_selected(self):
        if self.primary_method not in self.methods:
            raise ValueError("primary_method must occur in methods")
        return self


class GenomeBuild(str, Enum):
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"


class Population(str, Enum):
    EUR = "EUR"
    AFR = "AFR"
    EAS = "EAS"
    SAS = "SAS"
    AMR = "AMR"


class ModuleConfig(StrictModel):
    enabled: bool = False


class DelimitedTableReadConfig(StrictModel):
    """Shared schema for configuration-driven headed-table readers."""

    delimiter: str
    delimiter_candidates: list[str]
    sample_lines: int = Field(ge=1)
    infer_schema_length: int = Field(ge=1)
    maximum_columns: int = Field(ge=1)
    null_values: list[str]

    @field_validator("delimiter")
    @classmethod
    def known_delimiter(cls, value: str) -> str:
        if value != "auto" and value not in NAMED_DELIMITERS:
            raise ValueError("is not a supported delimiter name")
        return value

    @field_validator("delimiter_candidates", "null_values")
    @classmethod
    def unique_values(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique values")
        return values

    @field_validator("delimiter_candidates")
    @classmethod
    def known_delimiter_candidates(cls, values: list[str]) -> list[str]:
        unknown = [value for value in values if value not in NAMED_DELIMITERS]
        if unknown:
            raise ValueError(
                "contains unsupported delimiter names: %s" % ", ".join(unknown)
            )
        return values


class GeneRankingReportingConfig(StrictModel):
    top_gene_count: int = Field(ge=1)
    score_decimal_places: int = Field(ge=1, le=10)


class GenomicRegion(StrictModel):
    chromosome: str
    start: int = Field(ge=0)
    end: int = Field(gt=0)

    @model_validator(mode="after")
    def end_must_follow_start(self):
        if self.end <= self.start:
            raise ValueError("end must be greater than start")
        return self


class MHCAnalysisScopeConfig(StrictModel):
    """Shared policy for excluding MHC variants and coordinate-defined units."""

    policy: MHCPolicy
    region_override: GenomicRegion | None = None

    @property
    def excludes_snps(self) -> bool:
        return self.policy in {"exclude_snps", "exclude_both"}

    @property
    def excludes_genes(self) -> bool:
        return self.policy in {"exclude_genes", "exclude_both"}


class ChromosomeAnalysisScopeConfig(StrictModel):
    """Shared normalized chromosome-exclusion list."""

    exclude: list[str]

    @field_validator("exclude")
    @classmethod
    def normalized_unique_labels(cls, values: list[str]) -> list[str]:
        normalized = [str(value).strip().upper() for value in values]
        if any(not value for value in normalized):
            raise ValueError("must not contain empty chromosome labels")
        if len(normalized) != len(set(normalized)):
            raise ValueError("must contain unique chromosome labels")
        return normalized


class InputOutputConfig(StrictModel):
    input: Path | None = None
    output_directory: Path | None = None
    dataset_id: str | None = None
