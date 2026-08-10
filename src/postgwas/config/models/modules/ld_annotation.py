from postgwas.config.models.common import GenomeBuild, ModuleConfig, Population


class LDAnnotationConfig(ModuleConfig):
    genome_build: GenomeBuild
    population: Population
    include_unassigned: bool
