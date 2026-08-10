"""Strict configuration primitives shared by scientific modules."""

from enum import Enum
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, model_validator


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
