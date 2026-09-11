"""Validated configuration for region and LD-reference clumping."""

from __future__ import annotations

from pathlib import Path
import re
from typing import ClassVar, Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import (
    GenomeBuild,
    GenomicRegion,
    ModuleConfig,
    Population,
    StrictModel,
)
from postgwas.config.models.modules.gcta_cojo import (
    GctaCojoParallelOutputContract,
)
from postgwas.core.paths import validate_filename_component


LDClumpingMethod = Literal["region", "standard", "cojo-slct"]


class LDClumpingCojoConfig(StrictModel):
    """Post-COJO physical-locus definition; GCTA model settings stay canonical."""

    merge_distance_bp: int = Field(ge=0)
    index_pvalue: Literal["marginal", "joint"]
    parallel_output_contract: GctaCojoParallelOutputContract
    subordinate_dataset_id: str
    formatter_resolved_config_file: str
    formatter_completion_manifest_file: str
    gcta_resolved_config_file: str

    @field_validator("subordinate_dataset_id")
    @classmethod
    def valid_subordinate_dataset_id(cls, value: str) -> str:
        value = value.strip()
        required = ("{dataset_id}", "{population}")
        if any(value.count(token) != 1 for token in required):
            raise ValueError(
                "must contain {dataset_id} and {population} exactly once"
            )
        try:
            rendered = value.format(dataset_id="dataset", population="EUR")
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        validate_filename_component(
            rendered, "modules.ld_clumping.cojo.subordinate_dataset_id"
        )
        return value

    @field_validator(
        "formatter_resolved_config_file",
        "formatter_completion_manifest_file",
        "gcta_resolved_config_file",
    )
    @classmethod
    def valid_subordinate_metadata_path(cls, value: str, info) -> str:
        value = value.strip()
        required = ["{dataset_id}"]
        if info.field_name == "gcta_resolved_config_file":
            required.append("{mode}")
        if any(value.count(token) != 1 for token in required):
            raise ValueError(
                "must contain %s exactly once"
                % " and ".join(required)
            )
        try:
            rendered = Path(value.format(dataset_id="dataset_EUR", mode="slct"))
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if rendered.is_absolute() or ".." in rendered.parts:
            raise ValueError("must stay inside the subordinate output directory")
        if rendered.suffix != ".yaml":
            raise ValueError("must render a .yaml file")
        return value


class LDClumpingInputsConfig(StrictModel):
    vcf: Path | None = None
    dataset_id: str | None = None

    @field_validator("dataset_id")
    @classmethod
    def safe_dataset_id(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return validate_filename_component(
            value, "modules.ld_clumping.inputs.dataset_id"
        )


class LDClumpingReferenceConfig(StrictModel):
    directory: Path | None = None
    manifest_filename: str
    file_pattern: str
    reverse_file_pattern: str
    variant_inventory_pattern: str
    index_suffix: str
    format_version: int = Field(ge=1)
    orientation: Literal[
        "symmetric_first_endpoint", "upper_triangle_dual_index"
    ]
    columns: list[str]
    variant_inventory_columns: list[str]

    @field_validator("manifest_filename")
    @classmethod
    def safe_manifest_filename(cls, value: str) -> str:
        return validate_filename_component(
            value, "modules.ld_clumping.reference.manifest_filename"
        )

    @field_validator(
        "file_pattern", "reverse_file_pattern", "variant_inventory_pattern"
    )
    @classmethod
    def valid_file_pattern(cls, value: str, info) -> str:
        value = value.strip()
        required = ("{population}", "{chromosome}")
        if any(value.count(token) != 1 for token in required):
            raise ValueError(
                "must contain {population} and {chromosome} exactly once"
            )
        try:
            rendered = Path(value.format(population="EUR", chromosome="1"))
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if rendered.is_absolute() or ".." in rendered.parts:
            raise ValueError("must stay inside the configured reference directory")
        suffix = (
            ".variants.tsv.gz"
            if info.field_name == "variant_inventory_pattern"
            else ".ld.gz"
        )
        if not str(rendered).endswith(suffix):
            raise ValueError("must render a %s file" % suffix)
        return value

    @field_validator("index_suffix")
    @classmethod
    def valid_index_suffix(cls, value: str) -> str:
        if not value or Path(value).name != value:
            raise ValueError("must be a non-empty filename suffix")
        return value

    @field_validator("columns")
    @classmethod
    def valid_columns(cls, values: list[str]) -> list[str]:
        required = [
            "chromosome_a",
            "position_a",
            "variant_a",
            "chromosome_b",
            "position_b",
            "variant_b",
            "r2",
        ]
        if values != required:
            raise ValueError(
                "must describe the seven-column PLINK --r2 contract in order: %s"
                % ", ".join(required)
            )
        return values

    @field_validator("variant_inventory_columns")
    @classmethod
    def valid_variant_inventory_columns(cls, values: list[str]) -> list[str]:
        required = [
            "chromosome",
            "position",
            "reference_id",
            "canonical_id",
            "allele_1",
            "allele_2",
            "minor_allele_frequency",
        ]
        if values != required:
            raise ValueError(
                "must describe the seven-column variant inventory contract in "
                "order: %s" % ", ".join(required)
            )
        return values

    @model_validator(mode="after")
    def orientation_files(self):
        if self.orientation == "upper_triangle_dual_index":
            if self.reverse_file_pattern == self.file_pattern:
                raise ValueError(
                    "reverse_file_pattern must differ from file_pattern"
                )
        return self


class LDClumpingVcfFieldsConfig(StrictModel):
    chromosome: str
    position: str
    reference_allele: str
    alternate_allele: str
    variant_id: str
    effect: str
    standard_error: str
    allele_frequency: str
    log_pvalue: str
    ld_block: str

    @field_validator("*")
    @classmethod
    def nonempty_expression(cls, value: str) -> str:
        value = value.strip()
        if not value or "\t" in value or "\n" in value:
            raise ValueError("must be one non-empty bcftools query expression")
        return value

    @field_validator("ld_block")
    @classmethod
    def valid_ld_block_template(cls, value: str) -> str:
        if value.count("{population}") != 1:
            raise ValueError("must contain {population} exactly once")
        try:
            rendered = value.format(population="EUR")
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if not re.fullmatch(r"%INFO/[A-Za-z][A-Za-z0-9_.-]*", rendered):
            raise ValueError("must render one bcftools INFO-field expression")
        return value


class LDClumpingTableConfig(StrictModel):
    delimiter: str
    null_values: list[str]
    infer_schema_length: int = Field(ge=1)
    io_buffer_bytes: int = Field(ge=1)
    biallelic_include_expression: str
    compressor: str | None = None

    @field_validator("compressor")
    @classmethod
    def executable_name(cls, value: str | None) -> str | None:
        """Name a parallel gzip-compatible compressor, or null for Python gzip."""
        if value is None or not str(value).strip():
            return None
        return validate_filename_component(
            value, "modules.ld_clumping.table.compressor"
        )

    @field_validator("delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @field_validator("null_values")
    @classmethod
    def unique_null_values(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique values")
        return values

    @field_validator("biallelic_include_expression")
    @classmethod
    def nonempty_include_expression(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("must not be empty")
        return value


class LDClumpingComputeConfig(StrictModel):
    minimum_worker_memory_gb: float = Field(gt=0, allow_inf_nan=False)
    input_memory_multiplier: float = Field(ge=1, allow_inf_nan=False)


LDClumpingResultTable = Literal[
    "genomic_risk_loci",
    "lead_snps",
    "independent_significant_snps",
    "cojo_loci",
    "cojo_selected_signals",
]


# Stable columns produced by the standard and COJO result frames. Keeping this
# contract in the schema makes an invalid HTML-column override fail during
# configuration validation, before VCF extraction or chromosome processing.
_LD_CLUMPING_RESULT_COLUMNS = {
    "genomic_risk_loci": {
        "Genomic_locus", "Locus_Length_kb", "Merged_by_Distance",
        "Lead_uniqID", "Lead_rsID", "CHR", "POS", "START", "END", "ea",
        "nea", "eaf", "LP", "P_value", "Beta", "SE", "n_refsnps",
        "n_members", "n_5e_8", "n_5e_5", "n_0_05", "nIndSigSNPs",
        "nLeadSNPs", "IndSig_LD_Groups", "Lead_LD_Groups",
    },
    "lead_snps": {
        "lead_SNP_id", "chr", "start", "end", "p", "lp", "l_pos",
        "l_beta", "l_se", "l_rsid", "l_ea", "l_nea", "l_eaf", "is_list",
        "is_groups", "candidate_list", "gwas_candidate_list", "l_ld_list",
        "sum_refsnps", "sum_members", "Lead_Group", "Genomic_locus",
    },
    "independent_significant_snps": {
        "ind_sig_SNP_id", "chr", "start", "end", "pos", "candidate_list",
        "gwas_candidate_list", "p", "lp", "beta", "se", "rsID", "ea",
        "nea", "eaf", "ld_list", "n_refsnps", "n_members", "IndSig_Group",
        "Genomic_locus", "lead_SNP_id", "r2_with_Lead",
    },
    "cojo_loci": {
        "COJO_locus", "Chr", "start", "end", "span_bp", "index_SNP",
        "index_bp", "index_p", "index_pJ", "selected_signals",
        "selected_signal_ids",
    },
    "cojo_selected_signals": {
        "COJO_locus", "is_locus_index", "locus_start", "locus_end",
        "distance_to_previous_signal_bp", "Chr", "SNP", "bp", "freq",
        "refA", "b", "se", "p", "n", "freq_geno", "bJ", "bJ_se",
        "pJ", "LD_r", "estimation_status",
    },
}


class LDClumpingReportingConfig(StrictModel):
    result_table_columns: dict[LDClumpingResultTable, list[str]]

    @field_validator("result_table_columns")
    @classmethod
    def complete_unique_result_columns(
        cls,
        values: dict[LDClumpingResultTable, list[str]],
    ) -> dict[LDClumpingResultTable, list[str]]:
        required = {
            "genomic_risk_loci",
            "lead_snps",
            "independent_significant_snps",
            "cojo_loci",
            "cojo_selected_signals",
        }
        if set(values) != required:
            raise ValueError(
                "must define exactly these result tables: %s"
                % ", ".join(sorted(required))
            )
        for table_name, columns in values.items():
            if not columns or len(columns) != len(set(columns)):
                raise ValueError(
                    "%s must contain one or more unique column names"
                    % table_name
                )
            if any(
                not isinstance(column, str)
                or not column.strip()
                or "\t" in column
                or "\n" in column
                for column in columns
            ):
                raise ValueError(
                    "%s contains an invalid result column name" % table_name
                )
            unknown = sorted(
                set(columns) - _LD_CLUMPING_RESULT_COLUMNS[table_name]
            )
            if unknown:
                raise ValueError(
                    "%s contains columns not produced by that result "
                    "table: %s" % (table_name, ", ".join(unknown))
                )
        return values


class LDClumpingOutputLayoutConfig(StrictModel):
    shared_cojo_directory_fields: ClassVar[frozenset[str]] = frozenset({
        "cojo_formatter_directory",
        "cojo_root_directory",
    })

    region_working_table: str
    region_raw_table: str
    region_pruned_table: str
    region_significant_table: str
    region_significant_outside_ld_regions: str
    region_log: str
    standard_working_table: str
    standard_formatted_table: str
    standard_log: str
    standard_summary: str
    standard_hierarchy: str
    standard_independent_clusters: str
    standard_lead_clusters: str
    standard_reference_exclusions: str
    cojo_formatter_directory: str
    cojo_root_directory: str
    cojo_selected_signals: str
    cojo_loci: str
    cojo_excluded_variants: str
    summary_csv: str
    html_report: str
    canonical_log: str
    resolved_configuration: str

    @field_validator("*")
    @classmethod
    def safe_output_pattern(cls, value: str, info) -> str:
        value = value.strip()
        if (
            info.field_name not in cls.shared_cojo_directory_fields
            and "{dataset_id}" not in value
        ):
            raise ValueError("must contain {dataset_id}")
        try:
            rendered = Path(value.format(dataset_id="dataset", population="EUR"))
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if rendered.is_absolute() or ".." in rendered.parts:
            raise ValueError("must stay inside the configured output directory")
        return value

    @model_validator(mode="after")
    def unique_destinations(self):
        values = self.model_dump()
        if len(values) != len(set(values.values())):
            raise ValueError("LD-clumping output path patterns must be unique")
        # Shared COJO containers are safe because the validated subordinate
        # dataset ID is applied to every formatter and GCTA artifact filename.
        # Every direct scientific result and report remains population-specific.
        missing = sorted(
            name
            for name, pattern in values.items()
            if name not in self.shared_cojo_directory_fields
            and "{population}" not in pattern
        )
        if missing:
            raise ValueError(
                "scientific output patterns must contain {population}: %s"
                % ", ".join(missing)
            )
        rendered_root = Path(self.cojo_root_directory.format(
            dataset_id="dataset", population="EUR",
        ))
        rendered_formatter = Path(self.cojo_formatter_directory.format(
            dataset_id="dataset", population="EUR",
        ))
        if rendered_root not in rendered_formatter.parents:
            raise ValueError(
                "cojo_formatter_directory must be below cojo_root_directory"
            )
        if not self.summary_csv.endswith(".csv"):
            raise ValueError("summary_csv must end with .csv")
        if not self.html_report.endswith(".html"):
            raise ValueError("html_report must end with .html")
        return self


class LDClumpingConfig(ModuleConfig):
    inputs: LDClumpingInputsConfig
    output_directory: Path | None = None
    methods: list[LDClumpingMethod]
    genome_build: GenomeBuild
    population: Population
    lead_pvalue: float = Field(gt=0, le=1, allow_inf_nan=False)
    candidate_pvalue: float = Field(gt=0, le=1, allow_inf_nan=False)
    clump_r2: float = Field(ge=0, le=1, allow_inf_nan=False)
    lead_r2: float = Field(ge=0, le=1, allow_inf_nan=False)
    window_kb: int = Field(gt=0)
    merge_distance_bp: int = Field(ge=0)
    missing_index_action: Literal["error", "warning_skip"]
    minimum_reference_maf: float = Field(ge=0, le=0.5, allow_inf_nan=False)
    missing_chromosome_action: Literal["error", "warning"]
    remove_mhc: bool
    mhc_regions: dict[GenomeBuild, GenomicRegion]
    summary_pvalue_thresholds: dict[str, float]
    reference: LDClumpingReferenceConfig
    vcf_fields: LDClumpingVcfFieldsConfig
    table: LDClumpingTableConfig
    compute: LDClumpingComputeConfig
    cojo: LDClumpingCojoConfig
    reporting: LDClumpingReportingConfig
    output_layout: LDClumpingOutputLayoutConfig

    @field_validator("methods")
    @classmethod
    def unique_methods(cls, values: list[LDClumpingMethod]) -> list[LDClumpingMethod]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique clumping methods")
        return values

    @field_validator("mhc_regions")
    @classmethod
    def configured_mhc_regions(
        cls, values: dict[GenomeBuild, GenomicRegion]
    ) -> dict[GenomeBuild, GenomicRegion]:
        if not values:
            raise ValueError("must define at least one build-specific MHC region")
        return values

    @field_validator("summary_pvalue_thresholds")
    @classmethod
    def valid_summary_thresholds(cls, values: dict[str, float]) -> dict[str, float]:
        if not values:
            raise ValueError("must contain at least one named threshold")
        for name, threshold in values.items():
            validate_filename_component(
                name, "modules.ld_clumping.summary_pvalue_thresholds key"
            )
            if not 0 < threshold <= 1:
                raise ValueError("all summary P-value thresholds must be in (0, 1]")
        return values

    @model_validator(mode="after")
    def scientifically_consistent_thresholds(self):
        if self.candidate_pvalue < self.lead_pvalue:
            raise ValueError(
                "candidate_pvalue must be greater than or equal to lead_pvalue"
            )
        if self.lead_r2 > self.clump_r2:
            raise ValueError("lead_r2 must be less than or equal to clump_r2")
        if self.remove_mhc and self.genome_build not in self.mhc_regions:
            raise ValueError(
                "mhc_regions must define the selected genome_build when MHC removal "
                "is enabled"
            )
        return self
