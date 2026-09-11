from pathlib import Path

from pydantic import Field

from postgwas.config.models.common import GenomicRegion, StrictModel


class ExecutableResources(StrictModel):
    bash: str
    bcftools: str
    tabix: str
    plink: str
    plink2: str
    magma: str
    scdrs: str
    ldsc: str
    munge_sumstats: str
    bgenix: str
    ldstore: str
    finemap: str
    gcta: str
    rscript: str
    python: str
    mixer: str
    mixer_figures: str
    pigz: str | None = None


class MixerContainerResource(StrictModel):
    runtime: str
    image: str
    platform: str | None = None
    python: str
    mixer: str
    mixer_figures: str
    use_host_user: bool


class ContainerResources(StrictModel):
    mixer: MixerContainerResource


class GenomeResources(StrictModel):
    fasta: Path | None = None
    dbsnp: Path | None = None
    gene_locations: Path | None = None
    regions: dict[str, GenomicRegion] = Field(default_factory=dict)


class PopulationResources(StrictModel):
    plink_prefix: Path | None = None
    allele_frequencies: Path | None = None
    ld_scores: Path | None = None


class ResourcesConfig(StrictModel):
    root: Path | None = None
    executables: ExecutableResources
    containers: ContainerResources
    genomes: dict[str, GenomeResources]
    populations: dict[str, PopulationResources]
