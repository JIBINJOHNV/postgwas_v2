from typing import Literal

from pydantic import Field

from postgwas.config.models.common import ModuleConfig, Population, StrictModel


class PredLDConfig(StrictModel):
    minimum_r2: float = Field(ge=0, le=1)
    minimum_maf: float = Field(ge=0, le=0.5)
    mode: str
    correlation_method: Literal["pearson", "spearman"]


class RaissConfig(StrictModel):
    minimum_r2: float = Field(ge=0, le=1)


class SSImpConfig(StrictModel):
    minimum_r2: float = Field(ge=0, le=1)


class SummaryGWASImputationConfig(StrictModel):
    minimum_r2: float = Field(ge=0, le=1)


class ImputationEnginesConfig(StrictModel):
    pred_ld: PredLDConfig
    raiss: RaissConfig
    ssimp: SSImpConfig
    summary_gwas_imputation: SummaryGWASImputationConfig


class ImputationConfig(ModuleConfig):
    engine: Literal["pred_ld", "raiss", "ssimp", "summary_gwas_imputation"]
    population: Population
    engines: ImputationEnginesConfig
