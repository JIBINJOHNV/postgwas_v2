from postgwas.config.models.common import GenomeBuild, ModuleConfig


class QCSummaryConfig(ModuleConfig):
    target_build: GenomeBuild
    include_plots: bool
    include_rejection_reasons: bool
