from pathlib import Path
from typing import Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import GenomeBuild, ModuleConfig, StrictModel


FineMappingEngine = Literal["susie", "finemap"]


class SusieFittingConfig(StrictModel):
    main_max_iter: int = Field(ge=1)


class SusieRecoveryConfig(StrictModel):
    max_iter: int = Field(ge=1)
    repaired_max_iter: int = Field(ge=1)
    reduced_l_max: int = Field(ge=1)
    ld_retry_multiplier: float = Field(ge=1)
    ld_retry_add_seconds: int = Field(ge=0)
    timeout_multiplier: float = Field(ge=1)
    ld_repair_timeout_multiplier: float = Field(ge=1)
    repaired_timeout_multiplier: float = Field(ge=1)


class SusieLdValidationConfig(StrictModel):
    eigenvalue_tolerance: float = Field(gt=0)
    eigenvalue_floor: float = Field(gt=0)
    mismatch_warning_threshold: float = Field(ge=0, le=1)
    mismatch_failure_threshold: float = Field(gt=0, le=1)
    repair_maximum_change: float = Field(gt=0, le=1)

    @model_validator(mode="after")
    def validate_threshold_order(self):
        if self.eigenvalue_floor > self.eigenvalue_tolerance:
            raise ValueError(
                "eigenvalue_floor must not exceed eigenvalue_tolerance"
            )
        if self.mismatch_warning_threshold >= self.mismatch_failure_threshold:
            raise ValueError(
                "mismatch_warning_threshold must be less than "
                "mismatch_failure_threshold"
            )
        return self


class SusieExecutionConfig(StrictModel):
    ld_timeout_seconds: int = Field(gt=0)
    susie_timeout_seconds: int = Field(gt=0)
    plink_timeout_seconds: int = Field(gt=0)
    process_poll_seconds: float = Field(gt=0)
    termination_grace_seconds: float = Field(gt=0)
    stderr_tail_lines: int = Field(ge=1)
    verbose: bool


class SusieMemoryConfig(StrictModel):
    used_threshold_percent: float = Field(gt=0, le=100)
    check_interval_seconds: float = Field(gt=0)
    maximum_wait_seconds: int = Field(gt=0)
    meminfo_read_lines: int = Field(ge=1)


class SusiePlottingConfig(StrictModel):
    minimum_pip_to_label: float = Field(ge=0, le=1)
    background_point_size: float = Field(gt=0)
    background_alpha: float = Field(ge=0, le=1)
    credible_set_point_size: float = Field(gt=0)
    credible_set_alpha: float = Field(ge=0, le=1)
    maximum_ld_dimension: int = Field(ge=1)
    width: float = Field(gt=0)
    height: float = Field(gt=0)
    dpi: int = Field(ge=1)
    ld_palette_colors: int = Field(ge=2)
    legend_gradient_points: int = Field(ge=2)
    legend_title_size: float = Field(gt=0)
    legend_text_size: float = Field(gt=0)
    legend_key_size_cm: float = Field(gt=0)
    legend_scale: float = Field(gt=0)
    credible_set_legend_size: float = Field(gt=0)
    labels_per_credible_set: int = Field(ge=0)
    large_locus_threshold: int = Field(ge=1)
    medium_locus_threshold: int = Field(ge=1)
    large_locus_scale: float = Field(gt=0)
    medium_locus_scale: float = Field(gt=0)
    label_lift: float = Field(ge=0)
    label_size: float = Field(gt=0)
    label_segment_alpha: float = Field(ge=0, le=1)
    label_segment_width: float = Field(gt=0)
    panel_heights: list[float]
    panel_widths: list[float]

    @model_validator(mode="after")
    def validate_plotting_layout(self):
        if self.medium_locus_threshold >= self.large_locus_threshold:
            raise ValueError(
                "medium_locus_threshold must be less than large_locus_threshold"
            )
        if len(self.panel_heights) != 3 or any(x <= 0 for x in self.panel_heights):
            raise ValueError("panel_heights must contain three positive values")
        if len(self.panel_widths) != 2 or any(x <= 0 for x in self.panel_widths):
            raise ValueError("panel_widths must contain two positive values")
        return self


class SusieEngineConfig(StrictModel):
    max_causal_components: int = Field(ge=1)
    minimum_purity: float = Field(ge=0, le=1)
    credible_set_tolerance: float = Field(gt=0, lt=1)
    recovery_audit_filename: str = Field(pattern=r"^[^/\\]+\.tsv$")
    fitting: SusieFittingConfig
    recovery: SusieRecoveryConfig
    ld_validation: SusieLdValidationConfig
    execution: SusieExecutionConfig
    memory: SusieMemoryConfig
    plotting: SusiePlottingConfig

    @model_validator(mode="after")
    def validate_recovery_policy(self):
        if self.recovery.max_iter < self.fitting.main_max_iter:
            raise ValueError("recovery.max_iter must be at least fitting.main_max_iter")
        if self.recovery.repaired_max_iter < self.fitting.main_max_iter:
            raise ValueError(
                "recovery.repaired_max_iter must be at least fitting.main_max_iter"
            )
        if self.recovery.reduced_l_max > self.max_causal_components:
            raise ValueError(
                "recovery.reduced_l_max must not exceed max_causal_components"
            )
        return self


class FinemapEngineConfig(StrictModel):
    algorithm: Literal["sss", "cond"]
    ldstore_timeout_seconds: int = Field(gt=0)
    finemap_timeout_seconds: int = Field(gt=0)
    termination_grace_seconds: int = Field(gt=0)
    plink_memory_mb: int = Field(gt=0)
    bgen_bits: int = Field(ge=1, le=32)
    external_tool_threads: int = Field(ge=1)
    n_causal_snps: int = Field(ge=1)
    n_iter: int = Field(ge=1)
    n_conv_sss: int = Field(ge=1)
    prob_conv_sss_tol: float = Field(gt=0, le=1)
    n_configs_top: int = Field(ge=1)
    corr_config: float = Field(gt=0, le=1)
    pvalue_snps: float = Field(gt=0, le=1)
    cond_pvalue: float = Field(gt=0, le=1)
    prior_std: float = Field(gt=0)
    prior_k: bool
    force_n_samples: bool
    std_effects: bool
    flames_manifest_filename: str = Field(pattern=r"^[^/\\]+\.tsv$")


class FineMappingEnginesConfig(StrictModel):
    susie: SusieEngineConfig
    finemap: FinemapEngineConfig


class FineMappingSampleSizeConfig(StrictModel):
    policy: Literal["warn"]
    summary_statistic: Literal["median"]
    relative_range_warning_threshold: float = Field(ge=0)


class FineMappingLdResourceGuardConfig(StrictModel):
    maximum_variants_per_locus: int = Field(ge=1)
    finemap_peak_matrix_multiplier: float = Field(ge=1)
    susie_peak_matrix_multiplier: float = Field(ge=1)


class FineMappingValidationConfig(StrictModel):
    credible_set_coverage_tolerance: float = Field(gt=0, lt=1)
    ld_correlation_tolerance: float = Field(gt=0)
    ld_symmetry_tolerance: float = Field(gt=0)
    ld_diagonal_tolerance: float = Field(gt=0)
    ld_eigenvalue_tolerance: float = Field(gt=0)


class FineMappingInputConfig(StrictModel):
    """Direct inputs; pipeline execution replaces generated entries in memory."""

    locus_file: Path | None = None
    susie_summary_statistics_file: Path | None = None
    finemap_summary_statistics_file: Path | None = None
    ld_reference_prefix: Path | None = None


class FineMappingRuntimeConfig(StrictModel):
    tool_version_timeout_seconds: int = Field(gt=0)
    software_version_timeout_seconds: int = Field(gt=0)
    fallback_memory_gb: float = Field(gt=0)


class FineMappingCredibleSetFileColumnsConfig(StrictModel):
    """Column contract for the authoritative FLAMES-compatible set files."""

    rank: str
    variant_id: str
    posterior_probability: str

    @field_validator("*")
    @classmethod
    def nonempty_column_name(cls, value: str) -> str:
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("must be a non-empty column name")
        return cleaned

    @model_validator(mode="after")
    def distinct_columns(self):
        values = (self.rank, self.variant_id, self.posterior_probability)
        if len(set(values)) != len(values):
            raise ValueError("credible-set file column names must be distinct")
        return self


class FineMappingHtmlReportConfig(StrictModel):
    """Presentation-only controls for the self-contained scientific report."""

    page_size: int = Field(ge=1, le=500)
    probability_significant_digits: int = Field(ge=2, le=12)
    credible_set_file_columns: FineMappingCredibleSetFileColumnsConfig
    credible_set_columns: list[str] = Field(min_length=1)
    locus_status_columns: list[str] = Field(min_length=1)
    overlap_columns: list[str] = Field(min_length=1)
    preflight_columns: list[str] = Field(min_length=1)
    susie_recovery_columns: list[str] = Field(min_length=1)
    finemap_model_columns: list[str] = Field(min_length=1)

    @field_validator(
        "credible_set_columns",
        "locus_status_columns",
        "overlap_columns",
        "preflight_columns",
        "susie_recovery_columns",
        "finemap_model_columns",
    )
    @classmethod
    def unique_nonempty_columns(cls, values: list[str]) -> list[str]:
        cleaned = [value.strip() for value in values]
        if any(not value for value in cleaned):
            raise ValueError("must contain only non-empty column names")
        if len(cleaned) != len(set(cleaned)):
            raise ValueError("must not contain duplicate column names")
        return cleaned


class FineMappingOutputLayoutConfig(StrictModel):
    """Validated, single-source layout for published and intermediate outputs."""

    results_directory: str
    combined_results_directory: str
    primary_credible_sets_directory: str
    diagnostic_plots_directory: str
    fitted_models_directory: str
    quality_control_directory: str
    locus_logs_directory: str
    ld_diagnostics_directory: str
    run_metadata_directory: str
    downstream_inputs_directory: str
    downstream_flames_directory: str
    downstream_flames_annotations_directory: str
    intermediate_directory: str
    workers_directory: str
    locus_chunks_directory: str
    prepared_loci_directory: str
    summary_statistics_directory: str
    chromosome_cache_directory: str
    locus_summary_statistics_directory: str
    primary_flames_directory: str
    primary_flames_annotations_directory: str
    finemap_inputs_directory: str
    finemap_locus_work_directory: str
    finemap_selected_models_directory: str
    joint_rerun_directory: str
    pipeline_log_file: str
    pipeline_progress_file: str
    run_configuration_file: str
    resolved_susie_configuration_file: str
    susie_software_versions_file: str
    susie_r_configuration_file: str
    susie_r_software_versions_file: str
    susie_combined_log_file: str
    susie_combined_results_file: str
    susie_combined_credible_sets_file: str
    susie_qc_file: str
    susie_failed_loci_file: str
    susie_locus_progress_file: str
    preflight_validation_file: str
    html_report_file: str
    finemap_input_qc_file: str
    finemap_task_skips_file: str
    finemap_locus_status_file: str
    finemap_model_summary_file: str
    finemap_software_versions_file: str
    point_locus_windows_file: str
    range_locus_windows_file: str
    cleanup_successful_workers: bool

    @field_validator("*", mode="before")
    @classmethod
    def safe_relative_paths(cls, value, info):
        if info.field_name == "cleanup_successful_workers":
            return value
        path = Path(str(value))
        if path.is_absolute() or ".." in path.parts or not str(value).strip():
            raise ValueError("must be a non-empty path below the output directory")
        return str(value)

    @model_validator(mode="after")
    def validate_layout_relationships(self):
        top_level = [
            self.results_directory,
            self.quality_control_directory,
            self.run_metadata_directory,
            self.downstream_inputs_directory,
            self.intermediate_directory,
        ]
        if len(top_level) != len(set(top_level)):
            raise ValueError("top-level fine-mapping output directories must differ")

        relationships = {
            self.results_directory: [
                self.combined_results_directory,
                self.primary_credible_sets_directory,
                self.diagnostic_plots_directory,
            ],
            self.quality_control_directory: [
                self.locus_logs_directory,
                self.ld_diagnostics_directory,
            ],
            self.downstream_inputs_directory: [self.downstream_flames_directory],
            self.downstream_flames_directory: [
                self.downstream_flames_annotations_directory,
            ],
            self.intermediate_directory: [
                self.workers_directory,
                self.locus_chunks_directory,
                self.prepared_loci_directory,
                self.summary_statistics_directory,
                self.primary_flames_directory,
                self.finemap_inputs_directory,
                self.finemap_locus_work_directory,
                self.finemap_selected_models_directory,
                self.fitted_models_directory,
                self.joint_rerun_directory,
            ],
            self.primary_flames_directory: [
                self.primary_flames_annotations_directory,
            ],
            self.summary_statistics_directory: [
                self.chromosome_cache_directory,
                self.locus_summary_statistics_directory,
            ],
        }
        for parent_value, child_values in relationships.items():
            parent = Path(parent_value)
            for child_value in child_values:
                child = Path(child_value)
                if parent not in child.parents:
                    raise ValueError(
                        f"{child_value!r} must be below {parent_value!r}"
                    )

        directory_values = [
            value
            for name, value in self.model_dump().items()
            if name.endswith("_directory")
        ]
        if len(directory_values) != len(set(directory_values)):
            raise ValueError("fine-mapping output directory paths must be unique")

        file_relationships = {
            self.run_metadata_directory: [
                self.pipeline_log_file,
                self.pipeline_progress_file,
                self.run_configuration_file,
                self.resolved_susie_configuration_file,
                self.susie_software_versions_file,
                self.susie_r_configuration_file,
                self.susie_r_software_versions_file,
                self.susie_combined_log_file,
                self.finemap_software_versions_file,
            ],
            self.combined_results_directory: [
                self.susie_combined_results_file,
                self.susie_combined_credible_sets_file,
            ],
            self.results_directory: [self.html_report_file],
            self.quality_control_directory: [
                self.preflight_validation_file,
                self.susie_qc_file,
                self.susie_failed_loci_file,
                self.susie_locus_progress_file,
                self.finemap_input_qc_file,
                self.finemap_task_skips_file,
                self.finemap_locus_status_file,
                self.finemap_model_summary_file,
            ],
            self.prepared_loci_directory: [
                self.point_locus_windows_file,
                self.range_locus_windows_file,
            ],
        }
        for parent_value, file_values in file_relationships.items():
            parent = Path(parent_value)
            for file_value in file_values:
                if parent not in Path(file_value).parents:
                    raise ValueError(
                        f"{file_value!r} must be below {parent_value!r}"
                    )

        dataset_patterns = [
            self.susie_software_versions_file,
            self.susie_r_configuration_file,
            self.susie_r_software_versions_file,
            self.susie_combined_log_file,
            self.susie_combined_results_file,
            self.susie_combined_credible_sets_file,
            self.susie_qc_file,
            self.susie_failed_loci_file,
            self.susie_locus_progress_file,
            self.html_report_file,
        ]
        if any("{dataset_id}" not in value for value in dataset_patterns):
            raise ValueError(
                "Dataset-specific fine-mapping output file patterns must contain "
                "{dataset_id}"
            )
        if not self.html_report_file.lower().endswith(".html"):
            raise ValueError("html_report_file must end with .html")
        file_values = [
            value
            for name, value in self.model_dump().items()
            if name.endswith("_file")
        ]
        if len(file_values) != len(set(file_values)):
            raise ValueError("fine-mapping output file paths must be unique")
        for value in file_values:
            try:
                value.format(dataset_id="dataset")
            except (KeyError, ValueError) as exc:
                raise ValueError(
                    "fine-mapping output files support only {dataset_id}"
                ) from exc
        return self


class FineMappingSummaryStatisticsPreparationConfig(StrictModel):
    policy: Literal["exact_file_per_locus"]
    source_read_policy: Literal["single_streaming_pass"]
    validation_scope: Literal["analysed_chromosomes"]
    chunk_rows: int = Field(gt=0)
    schema_inference_length: int = Field(gt=0)
    separator: Literal["\t"]
    compression: Literal["infer", "gzip", "bz2", "zip", "xz", "zstd", "none"]
    chromosome_filename_template: str
    locus_filename_template: str
    chromosome_manifest_filename: str = Field(pattern=r"^[^/\\]+\.tsv$")
    locus_manifest_filename: str = Field(pattern=r"^[^/\\]+\.tsv$")
    preparation_summary_filename: str = Field(pattern=r"^[^/\\]+\.tsv$")

    @model_validator(mode="after")
    def validate_owned_paths_and_templates(self):
        try:
            chromosome_name = self.chromosome_filename_template.format(chromosome="1")
            locus_name = self.locus_filename_template.format(index=1)
        except (KeyError, ValueError) as exc:
            raise ValueError(
                "summary-statistics filename templates must accept {chromosome} "
                "and {index}, respectively"
            ) from exc
        if "{chromosome" not in self.chromosome_filename_template:
            raise ValueError("chromosome_filename_template must contain {chromosome}")
        if "{index" not in self.locus_filename_template:
            raise ValueError("locus_filename_template must contain {index}")
        for value in (chromosome_name, locus_name):
            if value in {"", ".", ".."} or "/" in value or "\\" in value:
                raise ValueError(
                    "summary-statistics filename templates must produce basenames"
                )
        return self


class FineMappingOverlapResolutionConfig(StrictModel):
    policy: Literal["rerun_connected_groups"]
    maximum_joint_region_kb: int | None = Field(default=None, gt=0)
    joint_failure_policy: Literal["exclude_and_report"]
    plan_filename: str = Field(pattern=r"^[^/\\]+\.tsv$")
    summary_filename: str = Field(pattern=r"^[^/\\]+\.tsv$")
    joint_locus_filename: str = Field(pattern=r"^[^/\\]+\.tsv$")
    final_combined_filename: str = Field(pattern=r"^[^/\\]+\.tsv$")
    configuration_filename: str = Field(pattern=r"^[^/\\]+\.json$")
    index_filename: str = Field(pattern=r"^[^/\\]+\.txt$")
    genomic_loci_filename: str = Field(pattern=r"^[^/\\]+\.tsv$")
    annotation_prefix: str = Field(pattern=r"^[A-Za-z0-9_.-]+$")
    overlap_group_prefix: str = Field(pattern=r"^[A-Za-z0-9_.-]+$")
    provenance_delimiter: Literal[";"]
    primary_file_prefix: str = Field(pattern=r"^[A-Za-z0-9_.-]+$")
    joint_file_prefix: str = Field(pattern=r"^[A-Za-z0-9_.-]+$")

    @model_validator(mode="after")
    def validate_owned_paths(self):
        if self.primary_file_prefix == self.joint_file_prefix:
            raise ValueError("primary and joint file prefixes must differ")
        return self


class FineMappingConfig(ModuleConfig):
    engine: FineMappingEngine
    genome_build: GenomeBuild
    locus_type: Literal["range", "point"]
    locus_window_kb: int = Field(ge=0)
    locus_lp_threshold: float = Field(ge=0)
    skip_mhc: bool
    mhc_chromosome: int = Field(ge=1, le=26)
    mhc_start: int = Field(ge=0)
    mhc_end: int = Field(gt=0)
    memory_per_worker_gb: float = Field(gt=0)
    credible_set_coverage: float = Field(gt=0, le=1)
    input: FineMappingInputConfig
    sample_size: FineMappingSampleSizeConfig
    ld_resource_guard: FineMappingLdResourceGuardConfig
    validation: FineMappingValidationConfig
    runtime: FineMappingRuntimeConfig
    html_report: FineMappingHtmlReportConfig
    output_layout: FineMappingOutputLayoutConfig
    summary_statistics_preparation: FineMappingSummaryStatisticsPreparationConfig
    overlap_resolution: FineMappingOverlapResolutionConfig
    engines: FineMappingEnginesConfig

    @model_validator(mode="after")
    def validate_mhc_interval(self):
        if self.mhc_end <= self.mhc_start:
            raise ValueError("mhc_end must be greater than mhc_start")
        return self
