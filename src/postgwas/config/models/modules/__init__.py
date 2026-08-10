"""Typed configuration models owned by individual scientific modules."""

from postgwas.config.models.modules.allele_orientation import AlleleOrientationConfig
from postgwas.config.models.modules.caldera import CalderaConfig
from postgwas.config.models.modules.enrichment import EnrichmentConfig
from postgwas.config.models.modules.filtering import FilteringConfig
from postgwas.config.models.modules.fine_mapping import FineMappingConfig
from postgwas.config.models.modules.flames import FlamesConfig
from postgwas.config.models.modules.formatting import FormattingConfig
from postgwas.config.models.modules.gcta_cojo import GctaCojoConfig
from postgwas.config.models.modules.gcta_gene import GctaGeneConfig
from postgwas.config.models.modules.harmonisation import HarmonisationConfig
from postgwas.config.models.modules.imputation import ImputationConfig
from postgwas.config.models.modules.kpops import KPopsConfig
from postgwas.config.models.modules.ld_annotation import LDAnnotationConfig
from postgwas.config.models.modules.ld_clumping import LDClumpingConfig
from postgwas.config.models.modules.ldsc import LDSCConfig
from postgwas.config.models.modules.magma import MagmaConfig
from postgwas.config.models.modules.magmacovar import MagmaCovarConfig
from postgwas.config.models.modules.manhattan import ManhattanConfig
from postgwas.config.models.modules.mixer import MixerConfig
from postgwas.config.models.modules.pops import PopsConfig
from postgwas.config.models.modules.qc_summary import QCSummaryConfig
from postgwas.config.models.modules.single_cell import SingleCellConfig

__all__ = [name for name in globals() if name.endswith("Config")]
