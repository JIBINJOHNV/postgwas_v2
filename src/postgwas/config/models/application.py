from pydantic import Field, model_validator

from postgwas.config.models.common import StrictModel
from postgwas.config.models.execution import ExecutionConfig
from postgwas.config.models.logging import LoggingConfig
from postgwas.config.models.pipeline import PipelineConfig
from postgwas.config.models.resources import ResourcesConfig
from postgwas.config.models.run import RunConfig
from postgwas.config.models.modules import (
    AlleleOrientationConfig,
    CalderaConfig,
    EnrichmentConfig,
    FilteringConfig,
    FineMappingConfig,
    FlamesConfig,
    FormattingConfig,
    GctaCojoConfig,
    GctaGeneConfig,
    HarmonisationConfig,
    ImputationConfig,
    KPopsConfig,
    LDAnnotationConfig,
    LDClumpingConfig,
    LDSCConfig,
    MagmaConfig,
    MagmaCovarConfig,
    ManhattanConfig,
    MixerConfig,
    PopsConfig,
    QCSummaryConfig,
    SingleCellConfig,
)


class ModulesConfig(StrictModel):
    allele_orientation: AlleleOrientationConfig
    caldera: CalderaConfig
    enrichment: EnrichmentConfig
    filtering: FilteringConfig
    fine_mapping: FineMappingConfig
    flames: FlamesConfig
    formatting: FormattingConfig
    gcta_cojo: GctaCojoConfig
    gcta_gene: GctaGeneConfig
    harmonisation: HarmonisationConfig
    imputation: ImputationConfig
    kpops: KPopsConfig
    ld_annotation: LDAnnotationConfig
    ld_clumping: LDClumpingConfig
    ldsc: LDSCConfig
    magma: MagmaConfig
    magmacovar: MagmaCovarConfig
    manhattan: ManhattanConfig
    mixer: MixerConfig
    pops: PopsConfig
    qc_summary: QCSummaryConfig
    single_cell: SingleCellConfig


class PostGWASConfig(StrictModel):
    config_version: int = Field(ge=1)
    run: RunConfig
    execution: ExecutionConfig
    logging: LoggingConfig
    resources: ResourcesConfig
    pipeline: PipelineConfig
    modules: ModulesConfig

    @model_validator(mode="after")
    def validate_pipeline_modules(self):
        known = set(ModulesConfig.model_fields)
        unknown = sorted(set(self.pipeline.modules) - known)
        if unknown:
            raise ValueError("Unknown pipeline modules: " + ", ".join(unknown))
        if self.pipeline.stop_after and self.pipeline.stop_after not in self.pipeline.modules:
            raise ValueError("pipeline.stop_after must name an enabled pipeline module")
        configured_builds = set(self.resources.genomes)
        if self.modules.filtering.genome_build.value not in configured_builds:
            raise ValueError(
                "modules.filtering.genome_build must be defined in resources.genomes"
            )
        if self.modules.mixer.genome_build not in configured_builds:
            raise ValueError(
                "modules.mixer.genome_build must be defined in resources.genomes"
            )
        gcta_gene = self.modules.gcta_gene
        if (
            gcta_gene.genome_build is not None
            and gcta_gene.genome_build not in configured_builds
        ):
            raise ValueError(
                "modules.gcta_gene.genome_build must be defined in resources.genomes"
            )
        if (
            gcta_gene.reference.population is not None
            and gcta_gene.reference.population not in self.resources.populations
        ):
            raise ValueError(
                "modules.gcta_gene.reference.population must be defined in "
                "resources.populations"
            )
        gcta_cojo = self.modules.gcta_cojo
        if (
            gcta_cojo.genome_build is not None
            and gcta_cojo.genome_build not in configured_builds
        ):
            raise ValueError(
                "modules.gcta_cojo.genome_build must be defined in resources.genomes"
            )
        if (
            gcta_cojo.reference.population is not None
            and gcta_cojo.reference.population not in self.resources.populations
        ):
            raise ValueError(
                "modules.gcta_cojo.reference.population must be defined in "
                "resources.populations"
            )
        harmonisation_builds = set(
            self.modules.harmonisation.vcf_processing.target_builds
        )
        if not harmonisation_builds.issubset(configured_builds):
            raise ValueError(
                "harmonisation target_builds must be defined in resources.genomes"
            )
        return self
