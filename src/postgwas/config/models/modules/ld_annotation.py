from pydantic import field_validator

from postgwas.config.models.common import GenomeBuild, ModuleConfig, Population


class LDAnnotationConfig(ModuleConfig):
    genome_build: GenomeBuild
    populations: list[Population]
    include_unassigned: bool

    @field_validator("populations")
    @classmethod
    def unique_populations(cls, values: list[Population]) -> list[Population]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique populations")
        return values
