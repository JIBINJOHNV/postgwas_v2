"""Strict configuration primitives shared by scientific modules."""

from enum import Enum
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from postgwas.core.io.delimiters import NAMED_DELIMITERS


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)


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


class InputOutputConfig(StrictModel):
    input: Path | None = None
    output_directory: Path | None = None
    dataset_id: str | None = None
