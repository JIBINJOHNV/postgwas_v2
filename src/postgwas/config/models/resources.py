from pathlib import Path

from postgwas.config.models.common import StrictModel


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


class PopulationResources(StrictModel):
    plink_prefix: Path | None = None
    allele_frequencies: Path | None = None
    ld_scores: Path | None = None
    ld_blocks: Path | None = None


class ResourcesConfig(StrictModel):
    root: Path | None = None
    executables: ExecutableResources
    containers: ContainerResources
    genomes: dict[str, GenomeResources]
    populations: dict[str, PopulationResources]
