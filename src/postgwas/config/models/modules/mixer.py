from datetime import datetime, timezone
import math
from pathlib import Path
from typing import Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import ModuleConfig, StrictModel


MixerAnalysis = Literal["univariate", "gsa", "all"]
MixerExecutionBackend = Literal["auto", "native", "docker"]
MixerQualityAction = Literal["ignore", "warning", "error"]
MixerSummaryStatistic = Literal[
    "point_estimate", "mean", "median", "std", "min", "max",
]
MixerFigureExtension = Literal["png", "svg"]
MixerLoglikeMethod = Literal["fast", "full"]


class MixerWorkflowConfig(StrictModel):
    fit_command: str
    test_command: str
    combine_command: str
    split_sumstats_command: str
    gsa_command: str
    figure_command: str
    chromosome_placeholder: str
    replicate_placeholder: str
    run_id_format: str

    @field_validator("*")
    @classmethod
    def nonempty_workflow_value(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @field_validator("run_id_format")
    @classmethod
    def safe_run_id_format(cls, value: str) -> str:
        rendered = datetime(2000, 1, 1, tzinfo=timezone.utc).strftime(value)
        if not rendered or rendered in {".", ".."} or Path(rendered).name != rendered:
            raise ValueError("must produce one safe directory name")
        return value


class MixerOutputLayout(StrictModel):
    log_file: str
    resolved_config_file: str
    fit_prefix: str
    test_prefix: str
    fit_replicate_prefix: str
    test_replicate_prefix: str
    summary_yaml: str
    summary_tsv: str
    figure_prefix: str
    gsa_split_pattern: str
    gsa_baseline_prefix: str
    gsa_full_prefix: str
    gsa_summary_yaml: str
    gsa_top_results_tsv: str

    @field_validator("*")
    @classmethod
    def nonempty_output_pattern(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("must not be empty")
        path = Path(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("must be a relative path below the output directory")
        return value


class MixerQualityConfig(StrictModel):
    missing_chromosome_action: MixerQualityAction
    unexpected_chromosome_action: MixerQualityAction
    seed_mismatch_action: MixerQualityAction
    nonfinite_qq_action: MixerQualityAction
    missing_log_metrics_action: MixerQualityAction
    missing_uncertainty_action: MixerQualityAction
    insufficient_model_power_action: MixerQualityAction


class MixerReportingConfig(StrictModel):
    enabled: bool
    generate_figures: bool
    figure_extensions: list[MixerFigureExtension]
    figure_statistics: list[MixerSummaryStatistic]
    table_delimiter: str
    null_value: str
    quality: MixerQualityConfig

    @field_validator("table_delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @field_validator("null_value")
    @classmethod
    def nonempty_null_value(cls, value: str) -> str:
        if not value:
            raise ValueError("must not be empty")
        return value

    @field_validator("figure_extensions", "figure_statistics")
    @classmethod
    def nonempty_unique_values(cls, value: list[str]) -> list[str]:
        if not value or len(value) != len(set(value)):
            raise ValueError("must contain one or more unique values")
        return value

    @model_validator(mode="after")
    def figures_require_reporting(self):
        if self.generate_figures and not self.enabled:
            raise ValueError("generate_figures requires reporting.enabled: true")
        return self


class MixerGsaResultColumns(StrictModel):
    identifier: str
    gene_count: str
    heritability: str
    heritability_standard_error: str
    baseline_heritability: str
    enrichment: str
    enrichment_standard_error: str
    snp_count: str
    log_likelihood_difference: str
    degrees_of_freedom: str
    model_aic: str

    @field_validator("*")
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @model_validator(mode="after")
    def unique_columns(self):
        values = list(self.model_dump().values())
        if len(values) != len(set(values)):
            raise ValueError("result column names must be unique")
        return self


class MixerUnivariateConfig(StrictModel):
    fit_extract_file_pattern: str | None = None
    replicate_indices: list[int]

    @field_validator("fit_extract_file_pattern")
    @classmethod
    def nonempty_optional_pattern(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or a non-empty path pattern")
        return value

    @field_validator("replicate_indices")
    @classmethod
    def positive_unique_replicates(cls, values: list[int]) -> list[int]:
        if (
            not values
            or len(values) != len(set(values))
            or any(value < 1 for value in values)
        ):
            raise ValueError(
                "must contain one or more unique positive replicate indices"
            )
        return values


class MixerGsaConfig(StrictModel):
    annotation_file_pattern: str | None = None
    loadlib_file_pattern: str | None = None
    baseline_go_file: str | None = None
    model_go_file: str | None = None
    test_go_file: str | None = None
    use_complete_tag_indices: bool
    exclude_ranges: list[str]
    hardprune_maf: float = Field(ge=0, le=0.5)
    hardprune_r2: float = Field(ge=0, le=1)
    z_max: float | None = Field(default=None, gt=0, allow_inf_nan=False)
    go_extend_bp: int = Field(ge=0)
    go_all_genes_label: str
    loglike_difference_method: MixerLoglikeMethod
    standard_error_samples: int = Field(ge=1)
    adam_epoch: list[int] | None = None
    adam_step: list[float] | None = None
    go_file_delimiter: str
    go_file_required_columns: list[str]
    go_identifier_column: str
    result_delimiter: str
    result_columns: MixerGsaResultColumns
    evidence_aic_threshold: float = Field(allow_inf_nan=False)
    top_results: int = Field(ge=1)

    @field_validator(
        "annotation_file_pattern", "loadlib_file_pattern", "baseline_go_file",
        "model_go_file", "test_go_file",
    )
    @classmethod
    def nonempty_optional_path(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or a non-empty path")
        return value

    @field_validator("go_all_genes_label", "go_identifier_column")
    @classmethod
    def nonempty_label(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @field_validator("go_file_delimiter", "result_delimiter")
    @classmethod
    def one_character_table_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @field_validator("exclude_ranges", "go_file_required_columns")
    @classmethod
    def nonempty_unique_text(cls, values: list[str]) -> list[str]:
        if (
            not values
            or len(values) != len(set(values))
            or any(not str(value).strip() for value in values)
        ):
            raise ValueError("must contain one or more unique non-empty values")
        return values

    @field_validator("adam_epoch")
    @classmethod
    def positive_epochs(cls, values: list[int] | None) -> list[int] | None:
        if values is not None and (not values or any(value < 1 for value in values)):
            raise ValueError("must be null or contain positive integers")
        return values

    @field_validator("adam_step")
    @classmethod
    def positive_steps(cls, values: list[float] | None) -> list[float] | None:
        if values is not None and (
            not values
            or any(not math.isfinite(value) or value <= 0 for value in values)
        ):
            raise ValueError("must be null or contain positive values")
        return values

    @model_validator(mode="after")
    def paired_adam_schedule(self):
        if (self.adam_epoch is None) != (self.adam_step is None):
            raise ValueError("adam_epoch and adam_step must both be null or both supplied")
        if self.adam_epoch is not None and len(self.adam_epoch) != len(self.adam_step or []):
            raise ValueError("adam_epoch and adam_step must have the same length")
        return self


class MixerConfig(ModuleConfig):
    analysis: MixerAnalysis
    execution_backend: MixerExecutionBackend
    genome_build: str
    bim_file_pattern: str | None = None
    ld_file_pattern: str | None = None
    chromosomes: list[str]
    fit_arguments: list[str]
    test_arguments: list[str]
    univariate: MixerUnivariateConfig
    gsa: MixerGsaConfig
    workflow: MixerWorkflowConfig
    output_layout: MixerOutputLayout
    reporting: MixerReportingConfig

    @model_validator(mode="after")
    def validate_patterns(self):
        if not self.genome_build.strip():
            raise ValueError("genome_build must not be empty")
        placeholder = self.workflow.chromosome_placeholder
        for name, value in (
            ("bim_file_pattern", self.bim_file_pattern),
            ("ld_file_pattern", self.ld_file_pattern),
            ("gsa.annotation_file_pattern", self.gsa.annotation_file_pattern),
            ("gsa.loadlib_file_pattern", self.gsa.loadlib_file_pattern),
        ):
            if value is not None and placeholder not in value:
                raise ValueError(
                    "%s must contain the chromosome placeholder %r"
                    % (name, placeholder)
                )
        replicate_placeholder = self.workflow.replicate_placeholder
        extract_pattern = self.univariate.fit_extract_file_pattern
        if extract_pattern is not None and replicate_placeholder not in extract_pattern:
            raise ValueError(
                "univariate.fit_extract_file_pattern must contain the replicate "
                "placeholder %r" % replicate_placeholder
            )
        required_output_tokens = {
            "log_file": ("{dataset_id}", "{run_id}"),
            "fit_prefix": ("{dataset_id}",),
            "test_prefix": ("{dataset_id}",),
            "fit_replicate_prefix": ("{dataset_id}", "{replicate}"),
            "test_replicate_prefix": ("{dataset_id}", "{replicate}"),
            "summary_yaml": ("{dataset_id}",),
            "summary_tsv": ("{dataset_id}",),
            "figure_prefix": ("{dataset_id}",),
            "gsa_split_pattern": ("{dataset_id}", placeholder),
            "gsa_baseline_prefix": ("{dataset_id}",),
            "gsa_full_prefix": ("{dataset_id}",),
            "gsa_summary_yaml": ("{dataset_id}",),
            "gsa_top_results_tsv": ("{dataset_id}",),
        }
        for name, tokens in required_output_tokens.items():
            pattern = getattr(self.output_layout, name)
            missing = [token for token in tokens if token not in pattern]
            if missing:
                raise ValueError(
                    "output_layout.%s must contain %s"
                    % (name, ", ".join(repr(token) for token in missing))
                )
        if len(self.chromosomes) != len(set(self.chromosomes)):
            raise ValueError("chromosomes must not contain duplicate values")
        if (
            self.analysis in {"univariate", "all"}
            and self.reporting.figure_statistics[0] not in {"mean", "median"}
        ):
            raise ValueError(
                "reporting.figure_statistics must begin with an aggregate "
                "location statistic for replicated univariate MiXeR results"
            )
        try:
            chromosomes = [int(value) for value in self.chromosomes]
        except (TypeError, ValueError) as exc:
            raise ValueError("chromosomes must contain positive integer labels") from exc
        if (
            not chromosomes
            or any(value < 1 for value in chromosomes)
            or chromosomes != list(range(chromosomes[0], chromosomes[-1] + 1))
        ):
            raise ValueError(
                "chromosomes must be one ordered, consecutive range of positive "
                "integer labels"
            )

        protected = {
            "--bim-file", "--chr2use", "--ld-file", "--load-params", "--out",
            "--seed", "--threads", "--trait1-file", "--extract", "--json",
            "--rep2use",
        }
        for field in ("fit_arguments", "test_arguments"):
            arguments = getattr(self, field)
            invalid = [
                token for token in arguments
                if not str(token).strip()
                or token in protected
                or token.startswith("--trait2")
                or token in {"fit2", "test2"}
            ]
            if invalid:
                raise ValueError(
                    "%s contains protected or multi-trait arguments: %s"
                    % (field, ", ".join(repr(value) for value in invalid))
                )
        return self


__all__ = [
    "MixerAnalysis",
    "MixerExecutionBackend",
    "MixerFigureExtension",
    "MixerQualityAction",
    "MixerSummaryStatistic",
]
