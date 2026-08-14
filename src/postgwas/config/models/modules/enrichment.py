from pydantic import Field

from postgwas.config.models.common import ModuleConfig


class EnrichmentConfig(ModuleConfig):
    providers: list[str]
    organism: str
    significance_threshold: float = Field(gt=0, le=1)
    correction_method: str
    minimum_genes: int = Field(ge=1)
    reference_set: str = Field(min_length=1)
