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
        if {"ld_annotation", "ld_clumping"}.issubset(self.pipeline.modules):
            annotation_field = (
                "%INFO/" + self.modules.ld_annotation.info_field_template
            )
            clumping_field = self.modules.ld_clumping.vcf_fields.ld_block
            if annotation_field != clumping_field:
                raise ValueError(
                    "modules.ld_annotation.info_field_template and "
                    "modules.ld_clumping.vcf_fields.ld_block must describe the "
                    "same population INFO field when both modules are in the pipeline"
                )
        configured_builds = set(self.resources.genomes)
        filtering_mhc_builds = {
            build.value for build in self.modules.filtering.mhc_regions
        }
        missing_filtering_mhc_builds = configured_builds - filtering_mhc_builds
        if self.modules.filtering.remove_mhc and missing_filtering_mhc_builds:
            raise ValueError(
                "modules.filtering.mhc_regions must define every configured genome "
                "build when MHC removal is enabled; missing: %s"
                % ", ".join(sorted(missing_filtering_mhc_builds))
            )
        ld_clumping = self.modules.ld_clumping
        if ld_clumping.genome_build.value not in configured_builds:
            raise ValueError(
                "modules.ld_clumping.genome_build must be defined in "
                "resources.genomes"
            )
        if ld_clumping.population.value not in self.resources.populations:
            raise ValueError(
                "modules.ld_clumping.population must be defined in "
                "resources.populations"
            )
        ld_clumping_mhc_builds = {
            build.value for build in ld_clumping.mhc_regions
        }
        missing_ld_clumping_mhc_builds = (
            configured_builds - ld_clumping_mhc_builds
        )
        if ld_clumping.remove_mhc and missing_ld_clumping_mhc_builds:
            raise ValueError(
                "modules.ld_clumping.mhc_regions must define every configured "
                "genome build when MHC removal is enabled; missing: %s"
                % ", ".join(sorted(missing_ld_clumping_mhc_builds))
            )
        qc_summary = self.modules.qc_summary
        qc_mhc_builds = {
            build.value for build in qc_summary.rules.mhc_regions
        }
        missing_qc_mhc_builds = configured_builds - qc_mhc_builds
        if qc_summary.rules.remove_mhc and missing_qc_mhc_builds:
            raise ValueError(
                "modules.qc_summary.rules.mhc_regions must define every configured "
                "genome build when MHC assessment is enabled; missing: %s"
                % ", ".join(sorted(missing_qc_mhc_builds))
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
