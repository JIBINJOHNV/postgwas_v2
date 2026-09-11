"""Standalone command-line boundary for GWAS harmonisation."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import shutil
import sys
import traceback
from pathlib import Path
from typing import Any

from postgwas.config import load_configuration, write_resolved_configuration
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.config.merger import deep_merge
from postgwas.cli.compute import get_compute_parser
from postgwas.core.errors import ConfigurationError
from postgwas.core.paths import configured_output_path
from postgwas.core.screen_logging import screen_recording_active
from postgwas.modules.harmonisation.sample_sheet import (
    HarmonisationSampleSheetRow,
    load_harmonisation_sample_sheet,
    to_harmonisation_input,
)
from postgwas.modules.harmonisation.policies import POLICIES, load_policies
from postgwas.modules.harmonisation.qc_reporting import (
    effect_scale_description,
    summarise_strand_orientation,
)
from postgwas.modules.harmonisation.html_report import (
    write_dataset_report as write_dataset_html_report,
    write_run_report as write_run_html_report,
)
from postgwas.core.pipeline_logging import write_log_record as _append_run_log
from postgwas.modules.harmonisation.service import (
    ConfigError,
    PipelineError,
    run_harmonisation_pipeline,
)
from postgwas.modules.harmonisation.concordance.errors import (
    ConcordanceValidationError,
)
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_default,
    format_cli_examples,
    help_with_default,
    screen_field,
    screen_line,
)
from postgwas.modules.harmonisation.vcf_processing import require_binaries


_RESOLVED_RESOURCE_PATHS = (
    "root",
    "executables.bash",
    "executables.bcftools",
    "executables.python",
    "executables.tabix",
    "executables.pigz",
)

_PREPARATION_TITLE = "Preparing PostGWAS harmonisation"
_PREPARATION_STAGE = "Validating configuration and required resources"

_RUN_SUMMARY_FIELDS = (
    # Dataset identity and terminal state.
    "dataset_id",
    "input_file",
    "status",
    "failure_reason",
    # Dataset step 3: initial input validation.
    "input_variants",
    "input_ready_variants",
    "input_ready_snps",
    "input_ready_indels_or_other_variants",
    "parsed_variants",
    # Dataset step 5: genome-build inference.
    "genome_build",
    "genome_build_detected",
    "genome_build_source",
    # Dataset steps 6-7: study-wide statistical and strand decisions.
    "effect_type",
    "effect_type_detected",
    "effect_type_source",
    "effect_type_automatically_detected",
    "effect_scale",
    "p_value_type",
    "p_value_type_detected",
    "p_value_type_source",
    "p_value_type_automatically_detected",
    "frequency_type",
    "frequency_type_source",
    "strand_consensus",
    # Chromosome processing: combined orientation evidence.
    "strand_status",
    "strand_metadata_complete",
    "strand_reference_panel",
    "strand_reference_population",
    "strand_reference_file_count",
    "strand_reference_files",
    "strand_chromosomes_expected",
    "strand_chromosomes_completed",
    "strand_chromosomes_summarized",
    "strand_variants_evaluated",
    "strand_matched",
    "strand_forward",
    "strand_forward_swapped",
    "strand_reverse_complement",
    "strand_reverse_complement_swapped",
    "strand_reference_unmatched_detected",
    "strand_reference_unmatched_retained",
    "strand_reference_unmatched",
    "strand_palindromic_frequency_conflict",
    "strand_palindromic_orientation_unavailable",
    "strand_palindromic_frequency_discordant",
    "strand_palindromic_ambiguous",
    "strand_reference_ambiguous",
    "strand_removed_total",
    "strand_accounting_balanced",
    # Chromosome processing and GWAS-to-VCF export.
    "pre_vcf_invalid_effect_statistic_removals",
    "rejected_variants",
    "unprocessed_failed_chromosome_variants",
    "row_accounting_balanced",
    "row_accounting_complete",
    "harmonised_variants",
    # Post-merge step 2: population-frequency similarity.
    "closest_population",
    "comparison_af_panel",
    "comparison_af_population",
    # Post-merge step 4: raw merged VCF and virtual QC.
    "final_vcf_variants",
    "final_vcf_snps",
    "final_vcf_indels_or_other_variants",
    "qc_active_rule_count",
    "qc_passed_variants",
    "qc_passed_snps",
    "qc_failed_snps",
    "qc_passed_snp_percent",
    "qc_failed_snp_percent",
    "af_comparable_variants",
    "af_concordant_variants",
    "af_concordance_percent",
    "af_mismatched_variants",
    "af_mismatch_percent",
    "study_af_missing",
    "reference_af_missing",
    "af_difference_cutoff",
    "qc_passed_missing_or_invalid_neff",
    "qc_passed_missing_imputation_score",
    "neff_reference_quantile",
    "neff_reference_value",
    "neff_minimum_fraction_of_reference",
    "neff_minimum_threshold",
    "raw_low_neff_variants",
    "raw_low_neff_percent",
    "qc_passed_low_neff_variants",
    "qc_passed_low_neff_percent",
    "qc_passed_neff_upper_outliers",
    "manifest",
    "screen_report",
    "html_report",
)


CLI_OVERRIDE_PATHS = {
    "resource_directory": "resources.root",
    "output_directory": "run.output_directory",
    "comparison_af_source": "modules.harmonisation.comparison_af.source",
    "comparison_af_column": "modules.harmonisation.comparison_af.column",
    "threads": "execution.threads",
    "memory_gb": "execution.memory_gb",
    "seed": "execution.random_seed",
    "validate": "modules.harmonisation.concordance_validation.enabled",
    "fixed_info": "modules.harmonisation.fixed_info.value",
    "zero_p_se_action": (
        "modules.harmonisation.policies.pvalue.zero_missing_se"
    ),
    "keep_gwas2vcf_intermediate": (
        "modules.harmonisation.policies.vcf.keep_gwas2vcf_intermediate"
    ),
}


class RunSummaryError(PipelineError):
    """The required multi-dataset harmonisation summary could not be written."""


class HtmlReportError(RunSummaryError):
    """A required harmonisation HTML report could not be written."""


class ScreenReportError(PipelineError):
    """A required per-dataset screen transcript could not be written."""


_ANSI_ESCAPE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")


def _dataset_output_paths(
    config,
    rows: list[HarmonisationSampleSheetRow],
    *,
    output_key: str,
    error_type: type[PipelineError],
) -> dict[str, Path]:
    """Resolve one configured dataset-root output for every selected row."""
    output_root = Path(config.run.output_directory).expanduser().resolve()
    layout = config.modules.harmonisation.output_layout.root
    paths = {}
    for row in rows:
        dataset_root = configured_output_path(
            output_root,
            layout["dataset_directory"],
            error_type=error_type,
            dataset_id=row.dataset_id,
        )
        paths[row.dataset_id] = configured_output_path(
            dataset_root,
            layout[output_key],
            error_type=error_type,
            dataset_id=row.dataset_id,
        )
    return paths


def _screen_report_paths(
    config,
    rows: list[HarmonisationSampleSheetRow],
) -> dict[str, Path]:
    """Resolve the configured screen-report path for every selected dataset."""
    return _dataset_output_paths(
        config, rows, output_key="screen_report", error_type=ScreenReportError,
    )


def _html_report_paths(
    config,
    rows: list[HarmonisationSampleSheetRow],
) -> tuple[dict[str, Path], Path]:
    """Resolve the per-dataset and combined HTML-report paths."""
    dataset_paths = _dataset_output_paths(
        config, rows, output_key="html_report", error_type=HtmlReportError,
    )
    output_root = Path(config.run.output_directory).expanduser().resolve()
    metadata = _runtime_path(config, output_root, "top_metadata_directory")
    run_path = _runtime_path(config, metadata, "run_html_report")
    return dataset_paths, run_path


def _html_report_context(
    config,
    row: HarmonisationSampleSheetRow,
    *,
    sample_sheet: str | Path | None = None,
) -> dict[str, Any]:
    """Reuse resolved inputs for presentation, including before analysis starts.

    This context belongs to the HTML view only. It does not modify the engine
    manifest, resolved policies, input mappings or the run-summary CSV schema.
    """
    harmonisation = config.modules.harmonisation
    return {
        "input_mapping": row.model_dump(mode="json"),
        "sample_sheet_warnings": list(row.normalisation_warnings),
        "sample_sheet": str(sample_sheet) if sample_sheet is not None else None,
        "concordance_requested": harmonisation.concordance_validation.enabled,
        "configuration": {
            "harmonisation": harmonisation.model_dump(mode="json"),
            "qc_summary": config.modules.qc_summary.model_dump(mode="json"),
        },
        "resource_directory": (
            str(config.resources.root) if config.resources.root is not None else None
        ),
        "output_directory": str(config.run.output_directory),
    }


class _DatasetScreenRouter:
    """Copy stdout to one required report per dataset and optionally the terminal."""

    def __init__(self, paths: dict[str, Path], *, display: bool):
        self.paths = dict(paths)
        self.display = bool(display)
        self._terminal = None
        self._active_dataset: str | None = None
        self._active_handle = None

    def __enter__(self):
        self._terminal = sys.stdout
        try:
            for path in self.paths.values():
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("", encoding="utf-8")
        except OSError as exc:
            raise ScreenReportError(
                "Cannot initialize the required harmonisation screen report %s: %s"
                % (path, exc)
            ) from exc
        sys.stdout = self
        return self

    def __exit__(self, exc_type, exc, traceback_value):
        sys.stdout = self._terminal
        self._close_active()
        return False

    @property
    def encoding(self):
        return getattr(self._terminal, "encoding", "utf-8")

    def isatty(self) -> bool:
        return bool(
            self.display
            and self._terminal is not None
            and getattr(self._terminal, "isatty", lambda: False)()
        )

    def fileno(self):
        if self._terminal is None:
            raise OSError("screen report router is not active")
        return self._terminal.fileno()

    def _close_active(self) -> None:
        if self._active_handle is not None:
            try:
                self._active_handle.close()
            finally:
                self._active_handle = None

    def select(self, dataset_id: str | None) -> None:
        """Route subsequent dataset output to one report, or shared output to all."""
        if dataset_id == self._active_dataset:
            return
        self._close_active()
        self._active_dataset = dataset_id
        if dataset_id is None:
            return
        try:
            path = self.paths[dataset_id]
        except KeyError as exc:
            raise ScreenReportError(
                "No screen-report path was resolved for dataset %s." % dataset_id
            ) from exc
        try:
            self._active_handle = path.open("a", encoding="utf-8")
        except OSError as exc:
            raise ScreenReportError(
                "Cannot append to the required harmonisation screen report %s: %s"
                % (path, exc)
            ) from exc

    def report_path(self, dataset_id: str) -> Path:
        try:
            return self.paths[dataset_id]
        except KeyError as exc:
            raise ScreenReportError(
                "No screen-report path was resolved for dataset %s." % dataset_id
            ) from exc

    @staticmethod
    def _report_text(text: str) -> str:
        """Keep terminal colour controls out of the persistent plain-text report."""
        return _ANSI_ESCAPE.sub("", text)

    def _append_shared(self, text: str) -> None:
        for path in self.paths.values():
            try:
                with path.open("a", encoding="utf-8") as handle:
                    handle.write(text)
                    handle.flush()
            except OSError as exc:
                raise ScreenReportError(
                    "Cannot append to the required harmonisation screen report "
                    "%s: %s" % (path, exc)
                ) from exc

    def write(self, value) -> int:
        text = str(value)
        report_text = self._report_text(text)
        try:
            if self._active_dataset is None:
                self._append_shared(report_text)
            elif self._active_handle is None:
                raise ScreenReportError(
                    "The screen report for dataset %s is not open."
                    % self._active_dataset
                )
            else:
                self._active_handle.write(report_text)
                self._active_handle.flush()
        except OSError as exc:
            path = self.report_path(self._active_dataset)
            raise ScreenReportError(
                "Cannot write the required harmonisation screen report %s: %s"
                % (path, exc)
            ) from exc
        if (
            self._terminal is not None
            and (self.display or screen_recording_active())
        ):
            self._terminal.write(text)
        return len(text)

    def flush(self) -> None:
        try:
            if self._active_handle is not None:
                self._active_handle.flush()
        except OSError as exc:
            path = self.report_path(self._active_dataset)
            raise ScreenReportError(
                "Cannot flush the required harmonisation screen report %s: %s"
                % (path, exc)
            ) from exc
        if (
            self._terminal is not None
            and (self.display or screen_recording_active())
        ):
            self._terminal.flush()


def _unit_interval_float(value: str) -> float:
    """Argparse type for a finite score in the closed unit interval."""
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("must be a number between 0 and 1") from exc
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise argparse.ArgumentTypeError("must be between 0 and 1 inclusive")
    return number


def get_harmonisation_parser(add_help: bool = False) -> argparse.ArgumentParser:
    """Build a parser whose defaults are displayed but never owned by argparse."""
    defaults = load_configuration()
    comparison = defaults.modules.harmonisation.comparison_af
    policy_defaults = load_policies(defaults.modules.harmonisation.policies)
    zero_p_se_policy = POLICIES["pvalue.zero_missing_se"]
    parser = argparse.ArgumentParser(
        add_help=add_help,
        formatter_class=AlignedRichHelpFormatter,
        parents=[get_compute_parser()],
    )
    files = parser.add_argument_group("Files and folders")
    files.add_argument(
        "--sample-sheet",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=(
            "REQUIRED. CSV or TSV file listing each GWAS dataset and the names "
            "of its columns. See examples/configs/harmonisation for templates."
        ),
    )
    files.add_argument(
        "--dataset-id",
        default=argparse.SUPPRESS,
        metavar="ID",
        help=(
            "Process only this dataset_id from the sample sheet. Without this option, "
            "all rows are processed."
        ),
    )
    files.add_argument(
        "--run-config",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=(
            "Optional YAML file with settings for this run. Command-line values "
            "override values in this file."
        ),
    )
    files.add_argument(
        "--resource-directory",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=(
            "Folder containing the PostGWAS reference data. REQUIRED unless "
            "resources.root is set in --run-config."
        ),
    )
    files.add_argument(
        "--output-directory",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=help_with_default(
            "Folder where results, logs, and the resolved configuration will be saved",
            defaults.run.output_directory,
        ),
    )
    frequency = parser.add_argument_group("Reference frequency check")
    frequency.add_argument(
        "--comparison-af-source",
        default=argparse.SUPPRESS,
        metavar="SOURCE",
        choices=tuple(comparison.available_sources),
        help=(
            "Reference VCF used for population-frequency annotation and QC; it "
            "is not used for strand orientation or to fill missing study "
            "frequencies.\n"
            "Options: %s.\n"
            "The resource directory must contain one VCF per chromosome:\n"
            "%s\n"
            "%s."
            % (
                ", ".join(comparison.available_sources),
                "\n".join(
                    "  %s: %s" % (source, comparison.resource_examples[source])
                    for source in comparison.available_sources
                ),
                format_cli_default(comparison.source),
            )
        ),
    )
    validation = parser.add_argument_group("Optional accuracy validation")
    validation.add_argument(
        "--validate",
        action="store_true",
        default=argparse.SUPPRESS,
        help=(
            "After harmonisation, compare each input dataset with its automatically "
            "selected same-build merged VCF. Disabled by default."
        ),
    )
    files.add_argument(
        "--keep_gwas2vcf_intermediate",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Keep the merged raw GWAS-to-VCF adapter VCF and its index after "
            "successful final validation. Without this option, that intermediate "
            "is deleted; the final GRCh37 and GRCh38 VCFs are always retained",
            policy_defaults.get("vcf.keep_gwas2vcf_intermediate"),
        ),
    )
    frequency.add_argument(
        "--comparison-af-column",
        default=argparse.SUPPRESS,
        metavar="COLUMN",
        help=help_with_default(
            "Population column to use in that reference, such as EUR, AFR, or EAS",
            defaults.modules.harmonisation.comparison_af.column,
        ),
    )
    quality = parser.add_argument_group("Imputation quality")
    quality.add_argument(
        "--fixed-info",
        type=_unit_interval_float,
        default=argparse.SUPPRESS,
        metavar="VALUE",
        help=(
            "Explicitly assign one INFO value to every variant only for datasets "
            "that provide neither an internal INFO column nor an external INFO "
            "file. Must be between 0 and 1. Internal INFO has first priority and "
            "external INFO has second priority. The assigned value is logged as "
            "user-provided, not measured per variant."
        ),
    )
    statistics = parser.add_argument_group("Effect reconstruction")
    statistics.add_argument(
        "--zero-p-se-action",
        dest="zero_p_se_action",
        choices=zero_p_se_policy.validator.members,
        default=argparse.SUPPRESS,
        metavar="ACTION",
        help=help_with_default(
            "What to do when a raw p-value is exactly zero and SE cannot be "
            "obtained from supplied SE or a usable Z score: 'fail' stops without "
            "discarding the variants, 'approximate' derives SE using the configured "
            "p-value floor and records the approximation, and 'reject' writes the "
            "variants to the rejected-variant output",
            policy_defaults.get("pvalue.zero_missing_se"),
        ),
    )
    return parser


def _module_policy_block(config) -> dict[str, Any]:
    harmonisation = config.modules.harmonisation
    configured = {
        "execution": {
            "total_cpu_budget": config.execution.threads,
            "memory_budget_gb": config.execution.memory_gb,
        },
        "logging": {"level": config.logging.file_level},
    }
    return deep_merge(configured, harmonisation.policies)


def harmonisation_executable_requirements(config) -> dict[str, str]:
    """Return the configured binaries required by harmonisation execution."""
    executables = config.resources.executables
    return {
        "bash": executables.bash,
        "bcftools": executables.bcftools,
        "python": executables.python,
        "tabix": executables.tabix,
    }


def _engine_defaults(config, resolved_executables=None) -> dict[str, Any]:
    """Build the engine configuration from the resolved run config."""
    harmonisation = config.modules.harmonisation
    mapping = harmonisation.external_eaf_mapping

    def reference_mapping(value) -> dict[str, str]:
        return {
            "chr": value.chromosome,
            "pos": value.position,
            "a1": value.effect_allele,
            "a2": value.other_allele,
            "delimiter": value.delimiter,
        }

    policy_block = _module_policy_block(config)
    return {
        "executables": dict(
            resolved_executables
            if resolved_executables is not None
            else harmonisation_executable_requirements(config)
        ),
        "executables_prevalidated": resolved_executables is not None,
        "compression_executable": config.resources.executables.pigz,
        "default_eaf_reference_source": harmonisation.default_eaf.source,
        "default_eaf_reference_column": harmonisation.default_eaf.column,
        "default_comparison_af_file": harmonisation.comparison_af.source,
        "default_comparison_af_column": harmonisation.comparison_af.column,
        "default_dbsnp": harmonisation.reference.dbsnp_source,
        "default_eaf_colmap": reference_mapping(
            harmonisation.default_eaf_mapping
        ),
        "external_eaf_colmap": reference_mapping(mapping),
        "external_info_colmap": reference_mapping(
            harmonisation.external_info_mapping
        ),
        "build_check_colmap": reference_mapping(
            harmonisation.build_check_mapping
        ),
        "resource_layout": harmonisation.resource_layout.model_dump(),
        "output_layout": dict(harmonisation.output_layout.root),
        "gwas2vcf_input": harmonisation.gwas2vcf_input.model_dump(),
        "vcf_processing": harmonisation.vcf_processing.model_dump(),
        "qc_summary": config.modules.qc_summary.model_dump(mode="json"),
        "population_frequency_qc": harmonisation.population_frequency_qc.model_dump(),
        "external_reference_staging": (
            harmonisation.external_reference_staging.model_dump()
        ),
        # Resolved global presentation settings are passed through to the
        # engine; this is not a second user-configurable default.
        "terminal_progress": {
            "enabled": config.logging.show_progress,
            "outcome_label_width": config.logging.terminal_label_width,
        },
        "policies": policy_block,
    }


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def _initial_run_summary_record(
    row: HarmonisationSampleSheetRow,
) -> dict[str, Any]:
    """Create one row before analysis so every selected dataset is represented."""
    record = {field: None for field in _RUN_SUMMARY_FIELDS}
    record.update({
        "dataset_id": row.dataset_id,
        "input_file": str(row.input_file),
        "status": "NOT_RUN",
    })
    return record


def _read_summary_manifest(path: Path, *, required: bool) -> dict[str, Any]:
    """Read one dataset manifest for run-level reporting."""
    if not path.is_file():
        if required:
            raise RunSummaryError(
                "Dataset completed without its required run manifest: %s" % path
            )
        return {}
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        if required:
            raise RunSummaryError(
                "Cannot read the completed dataset manifest %s: %s" % (path, exc)
            ) from exc
        return {}
    if not isinstance(document, dict):
        if required:
            raise RunSummaryError(
                "Dataset run manifest must contain a JSON object: %s" % path
            )
        return {}
    return document


def _single_line(value: Any) -> str | None:
    """Keep failure details readable as one CSV cell and one terminal line."""
    if value is None:
        return None
    text = " ".join(str(value).split())
    return text or None


def _percentage_value(numerator: Any, denominator: Any) -> float | None:
    """Return a numeric percentage while preserving unavailable denominators."""
    if numerator is None or denominator is None:
        return None
    try:
        denominator_value = float(denominator)
        if denominator_value <= 0:
            return None
        return 100.0 * float(numerator) / denominator_value
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def _automatically_detected(source: Any) -> bool | None:
    """Expose detector provenance explicitly without altering its raw label."""
    if source is None:
        return None
    return str(source).strip().lower() == "detector"


def _read_qc_assessment_report(
    manifest_path: Path,
    manifest: dict[str, Any],
    *,
    required: bool,
) -> dict[str, Any]:
    """Load the persisted QC metrics without rereading the merged VCF."""
    summary = manifest.get("qc_assessment") or {}
    report_value = (summary.get("reports") or {}).get("json")
    if not report_value:
        if required:
            raise RunSummaryError(
                "Completed dataset manifest does not identify its required QC "
                "assessment JSON: %s" % manifest_path
            )
        return {}
    report_path = Path(report_value).expanduser()
    if not report_path.is_absolute():
        report_path = manifest_path.parent / report_path
    report_path = report_path.resolve()
    if not report_path.is_file():
        if required:
            raise RunSummaryError(
                "Completed dataset is missing its required QC assessment JSON: %s"
                % report_path
            )
        return {}
    try:
        document = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RunSummaryError(
            "Cannot read the QC assessment JSON %s: %s" % (report_path, exc)
        ) from exc
    if not isinstance(document, dict):
        raise RunSummaryError(
            "QC assessment JSON must contain an object: %s" % report_path
        )
    return document


def _qc_run_summary_fields(
    pre_vcf: dict[str, Any],
    assessment: dict[str, Any],
) -> dict[str, Any]:
    """Flatten already-calculated final-QC evidence into CSV-ready values."""
    raw = assessment.get("raw") or {}
    passed = assessment.get("qc_passed") or {}
    total_snps = raw.get("num_snps")
    passed_snps = passed.get("num_snps")
    comparable = raw.get("af_comparable")
    mismatched = raw.get("af_difference_above_cutoff")

    def nonnegative_difference(total: Any, subset: Any) -> int | None:
        if total is None or subset is None:
            return None
        try:
            difference = int(total) - int(subset)
        except (TypeError, ValueError):
            return None
        return difference if difference >= 0 else None

    failed_snps = nonnegative_difference(total_snps, passed_snps)
    concordant = nonnegative_difference(comparable, mismatched)
    return {
        "final_vcf_snps": total_snps,
        "final_vcf_indels_or_other_variants": raw.get("num_non_snps"),
        "qc_active_rule_count": assessment.get("active_rule_count"),
        "qc_passed_snps": passed_snps,
        "qc_failed_snps": failed_snps,
        "qc_passed_snp_percent": _percentage_value(passed_snps, total_snps),
        "qc_failed_snp_percent": _percentage_value(failed_snps, total_snps),
        "af_comparable_variants": comparable,
        "af_concordant_variants": concordant,
        "af_concordance_percent": _percentage_value(concordant, comparable),
        "af_mismatched_variants": mismatched,
        "af_mismatch_percent": _percentage_value(mismatched, comparable),
        "study_af_missing": raw.get("study_af_missing"),
        "reference_af_missing": raw.get("external_af_missing"),
        "af_difference_cutoff": assessment.get("af_difference_cutoff"),
        "pre_vcf_invalid_effect_statistic_removals": pre_vcf.get(
            "total_variant_with_invalid_beta_se"
        ),
        "qc_passed_missing_or_invalid_neff": passed.get(
            "effective_sample_size_missing_or_invalid"
        ),
        "qc_passed_missing_imputation_score": passed.get("format_si_missing"),
        "neff_reference_quantile": assessment.get(
            "sample_size_reference_quantile"
        ),
        "neff_reference_value": raw.get(
            "effective_sample_size_reference_quantile_value"
        ),
        "neff_minimum_fraction_of_reference": assessment.get(
            "sample_size_minimum_fraction_of_reference"
        ),
        "neff_minimum_threshold": raw.get(
            "effective_sample_size_minimum_threshold"
        ),
        "raw_low_neff_variants": raw.get(
            "effective_sample_size_below_minimum_threshold"
        ),
        "raw_low_neff_percent": _percentage_value(
            raw.get("effective_sample_size_below_minimum_threshold"),
            raw.get("effective_sample_size_available"),
        ),
        "qc_passed_low_neff_variants": passed.get(
            "effective_sample_size_below_minimum_threshold"
        ),
        "qc_passed_low_neff_percent": _percentage_value(
            passed.get("effective_sample_size_below_minimum_threshold"),
            passed.get("effective_sample_size_available"),
        ),
        "qc_passed_neff_upper_outliers": passed.get(
            "effective_sample_size_above_outlier_threshold"
        ),
    }


def _run_summary_record(
    row: HarmonisationSampleSheetRow,
    *,
    status: str,
    manifest_path: Path,
    manifest: dict[str, Any],
    strand_reference_panel: str,
    strand_reference_population: str,
    comparison_af_panel: str,
    comparison_af_population: str,
    failure: Any = None,
    report_context: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Flatten already-calculated dataset and chromosome evidence into one row."""
    record = _initial_run_summary_record(row)
    dataset = manifest.get("dataset") or {}
    build = dataset.get("genome_build") or {}
    decisions = dataset.get("study_decisions") or {}
    strand = summarise_strand_orientation(dataset)
    population = manifest.get("population_frequency_qc") or {}
    pre_vcf = manifest.get("pre_vcf_summary") or {}
    assessment_summary = manifest.get("qc_assessment") or {}
    assessment = _read_qc_assessment_report(
        manifest_path,
        manifest,
        required=str(status).upper() in {"OK", "PARTIAL", "CONCORDANCE_FAILED"},
    )
    if report_context is not None:
        # Preserve the already-loaded rule and diagnostic evidence for the HTML
        # view; do not reopen the QC JSON or reread any scientific data.
        report_context["qc_assessment"] = assessment
    reconciliation = dataset.get("reconciliation") or {}
    inferred_build = build.get("inferred_build")
    forced_build = build.get("forced") is True
    manifest_failure = manifest.get("failure") or {}
    failure_reason = failure
    if failure_reason is None and isinstance(manifest_failure, dict):
        failure_reason = manifest_failure.get("message")

    record.update({
        "status": str(status).upper(),
        "failure_reason": _single_line(failure_reason),
        "genome_build": inferred_build,
        "genome_build_detected": None if forced_build else inferred_build,
        "genome_build_source": (
            "policy" if forced_build else "detector" if inferred_build else None
        ),
        "effect_type": decisions.get("effect_type"),
        "effect_type_detected": decisions.get("effect_type_detected"),
        "effect_type_source": decisions.get("effect_type_source"),
        "effect_type_automatically_detected": _automatically_detected(
            decisions.get("effect_type_source")
        ),
        "effect_scale": effect_scale_description(decisions.get("effect_type")),
        "p_value_type": decisions.get("pvalue_type"),
        "p_value_type_detected": decisions.get("pvalue_type_detected"),
        "p_value_type_source": decisions.get("pvalue_type_source"),
        "p_value_type_automatically_detected": _automatically_detected(
            decisions.get("pvalue_type_source")
        ),
        "frequency_type": decisions.get("frequency_type"),
        "frequency_type_source": decisions.get("eaf_is_maf_source"),
        "strand_consensus": decisions.get("strand") or strand.get("study_consensus"),
        "strand_status": strand.get("status"),
        "strand_metadata_complete": strand.get("metadata_complete"),
        "strand_reference_panel": strand_reference_panel,
        "strand_reference_population": (
            strand.get("reference_population") or strand_reference_population
        ),
        "strand_reference_file_count": strand.get("reference_file_count"),
        "strand_reference_files": json.dumps(
            strand.get("reference_files") or [], separators=(",", ":"),
        ),
        "strand_chromosomes_expected": strand.get("chromosomes_expected"),
        "strand_chromosomes_completed": strand.get("chromosomes_completed"),
        "strand_chromosomes_summarized": strand.get("chromosomes_summarized"),
        "strand_variants_evaluated": strand.get("variants_evaluated"),
        "strand_matched": strand.get("variants_matched"),
        "strand_forward": strand.get("forward"),
        "strand_forward_swapped": strand.get("forward_swapped"),
        "strand_reverse_complement": strand.get("reverse_complement"),
        "strand_reverse_complement_swapped": strand.get(
            "reverse_complement_swapped"
        ),
        "strand_reference_unmatched_detected": strand.get(
            "reference_unmatched_detected"
        ),
        "strand_reference_unmatched_retained": strand.get(
            "reference_unmatched_retained"
        ),
        "strand_reference_unmatched": strand.get("reference_unmatched"),
        "strand_palindromic_frequency_conflict": strand.get(
            "palindromic_frequency_conflict"
        ),
        "strand_palindromic_orientation_unavailable": strand.get(
            "palindromic_orientation_unavailable"
        ),
        "strand_palindromic_frequency_discordant": strand.get(
            "palindromic_frequency_discordant"
        ),
        "strand_palindromic_ambiguous": strand.get("palindromic_ambiguous"),
        "strand_reference_ambiguous": strand.get("reference_ambiguous"),
        "strand_removed_total": strand.get("removed_total"),
        "strand_accounting_balanced": strand.get("accounting_balanced"),
        "closest_population": population.get("closest_population"),
        "comparison_af_panel": comparison_af_panel,
        "comparison_af_population": comparison_af_population,
        "input_variants": (
            pre_vcf.get("total_variant_infile")
            if pre_vcf.get("total_variant_infile") is not None
            else build.get("input_variants")
        ),
        "input_ready_variants": (
            pre_vcf.get("total_variant_remaining_for_harmonisation")
            if pre_vcf.get("total_variant_remaining_for_harmonisation")
            is not None
            else build.get("input_variants")
        ),
        "input_ready_snps": pre_vcf.get("total_variant_ready_snps"),
        "input_ready_indels_or_other_variants": pre_vcf.get(
            "total_variant_ready_indels_or_other"
        ),
        "parsed_variants": reconciliation.get("rows_read"),
        "rejected_variants": reconciliation.get("rejected"),
        "unprocessed_failed_chromosome_variants": reconciliation.get(
            "unprocessed_failed_chromosome_rows"
        ),
        "row_accounting_balanced": reconciliation.get("balanced"),
        "row_accounting_complete": reconciliation.get("complete"),
        "harmonised_variants": (
            pre_vcf.get("total_variant_passed_chromosome_harmonisation")
            if pre_vcf.get("total_variant_passed_chromosome_harmonisation")
            is not None
            else reconciliation.get("rows_exported")
        ),
        "final_vcf_variants": assessment_summary.get("raw_variants"),
        "qc_passed_variants": assessment_summary.get("qc_passed_variants"),
        "manifest": str(manifest_path) if manifest_path.is_file() else None,
    })
    record.update(_qc_run_summary_fields(pre_vcf, assessment))
    return record


def _write_run_summary(
    path: Path,
    records: list[dict[str, Any]],
) -> None:
    """Atomically write the required one-row-per-dataset run summary."""
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=_RUN_SUMMARY_FIELDS)
    writer.writeheader()
    for record in records:
        writer.writerow({
            field: (
                "true" if record.get(field) is True
                else "false" if record.get(field) is False
                else "" if record.get(field) is None
                else record.get(field)
            )
            for field in _RUN_SUMMARY_FIELDS
        })
    try:
        _write_text(path, stream.getvalue())
    except OSError as exc:
        raise RunSummaryError(
            "Cannot write the required harmonisation run summary %s: %s"
            % (path, exc)
        ) from exc


def _write_html_reports(
    *,
    run_path: Path,
    records: dict[str, dict[str, Any]],
    manifests: dict[str, dict[str, Any]],
    dataset_paths: dict[str, Path],
    changed_dataset: str | None = None,
) -> None:
    """Write HTML views of existing records without rereading scientific data."""
    dataset_ids = (
        [changed_dataset]
        if changed_dataset is not None
        else list(records)
    )
    for dataset_id in dataset_ids:
        try:
            write_dataset_html_report(
                dataset_paths[dataset_id],
                records[dataset_id],
                manifests.get(dataset_id) or {},
            )
        except (KeyError, OSError, TypeError, ValueError) as exc:
            path = dataset_paths.get(dataset_id, dataset_id)
            raise HtmlReportError(
                "Cannot write the required dataset harmonisation HTML report "
                "%s: %s" % (path, exc)
            ) from exc
    try:
        write_run_html_report(
            run_path,
            list(records.values()),
            manifests,
        )
    except (OSError, TypeError, ValueError) as exc:
        raise HtmlReportError(
            "Cannot write the required combined harmonisation HTML report %s: %s"
            % (run_path, exc)
        ) from exc


def _run_summary_block(
    path: Path,
    html_path: Path,
    records: list[dict[str, Any]],
) -> str:
    """Show final report locations and dataset-status composition."""
    successful = sum(
        str(record.get("status") or "").upper() == "OK" for record in records
    )
    partial = sum(
        str(record.get("status") or "").upper() == "PARTIAL"
        for record in records
    )
    unfinished = sum(
        str(record.get("status") or "").upper() in {"NOT_RUN", "RUNNING"}
        for record in records
    )
    failed = sum(
        str(record.get("status") or "").upper()
        not in {"OK", "PARTIAL", "NOT_RUN", "RUNNING"}
        for record in records
    )
    return "\n".join([
        screen_line("decision", "Harmonisation run summary", indent=4),
        screen_field(
            "count", "Datasets",
            "%d total · %d OK · %d partial · %d failed · %d not finished"
            % (len(records), successful, partial, failed, unfinished),
            indent=8, label_width=18,
        ),
        screen_field(
            "info", "CSV", path,
            indent=8, label_width=18,
        ),
        screen_field(
            "info", "HTML", html_path,
            indent=8, label_width=18,
        ),
    ])


def _write_preflight_failure_run_summary(
    config,
    rows: list[HarmonisationSampleSheetRow],
    failure: Any,
    *,
    sample_sheet: str | Path | None = None,
) -> Path:
    """Represent every selected dataset when command preflight stops the run."""
    output_root = Path(config.run.output_directory).expanduser().resolve()
    report_paths = _screen_report_paths(config, rows)
    html_paths, run_html_path = _html_report_paths(config, rows)
    top_metadata = _runtime_path(
        config, output_root, "top_metadata_directory",
    )
    path = _runtime_path(config, top_metadata, "run_summary_file")
    records = []
    for row in rows:
        record = _initial_run_summary_record(row)
        record.update({
            "status": "PREFLIGHT_FAILED",
            "failure_reason": _single_line(failure),
            "screen_report": str(report_paths[row.dataset_id]),
            "html_report": str(html_paths[row.dataset_id]),
        })
        records.append(record)
        report = "\n".join([
            screen_line("run", _PREPARATION_TITLE, indent=4),
            screen_field(
                "info", "Dataset", row.dataset_id,
                indent=8, label_width=20,
            ),
            screen_field(
                "error", "Preflight failed", failure,
                indent=8, label_width=20,
            ),
        ]) + "\n"
        try:
            _write_text(report_paths[row.dataset_id], report)
        except OSError as exc:
            raise ScreenReportError(
                "Cannot write the required harmonisation screen report %s: %s"
                % (report_paths[row.dataset_id], exc)
            ) from exc
    _write_run_summary(path, records)
    record_map = {record["dataset_id"]: record for record in records}
    _write_html_reports(
        run_path=run_html_path,
        records=record_map,
        manifests={
            row.dataset_id: {
                "report_context": _html_report_context(
                    config, row, sample_sheet=sample_sheet,
                ),
            }
            for row in rows
        },
        dataset_paths=html_paths,
    )
    return path


def _runtime_path(config, root, setting: str, **values) -> Path:
    pattern = getattr(config.modules.harmonisation.runtime, setting)
    return configured_output_path(root, pattern, **values)


def _top_run_log(config, root) -> Path:
    metadata = _runtime_path(config, root, "top_metadata_directory")
    return configured_output_path(metadata, config.logging.filename)


def _fixed_info_cli_value(args, config) -> float | None:
    """Return the explicit fallback while rejecting run-config activation."""
    if hasattr(args, "fixed_info"):
        return float(args.fixed_info)
    if config.modules.harmonisation.fixed_info.value is not None:
        raise ConfigurationError(
            "fixed_info.value cannot be set in a run configuration. Supply "
            "the explicit command-line option --fixed-info VALUE instead."
        )
    return None


def _validate_info_fallback(
    rows: list[HarmonisationSampleSheetRow], fixed_info: float | None,
) -> None:
    """Fail preflight when any selected dataset has no usable INFO source."""
    missing = [
        row.dataset_id
        for row in rows
        if row.imputation_info_column is None
        and not (row.external_info_file and row.external_info_column)
    ]
    if missing and fixed_info is None:
        raise ConfigurationError(
            "No imputation-quality source is available for dataset(s): %s. "
            "Provide imputation_info_column, provide external_info_file together "
            "with external_info_column, or explicitly supply --fixed-info VALUE."
            % ", ".join(missing)
        )


def _sample_sheet_display(sample_sheet: str | Path | None) -> str | Path:
    """Resolve a CLI sheet path while labelling validated direct API input."""
    if sample_sheet is None:
        return "Direct API input — no sample-sheet path"
    return Path(sample_sheet).expanduser().resolve()


def _preparation_block(
    *,
    sample_sheet: str | Path | None,
    dataset_count: int,
    resource_directory: str | Path,
    output_directory: str | Path,
) -> str:
    """Render the run-level context before per-dataset output begins."""
    sample_sheet_display = _sample_sheet_display(sample_sheet)
    return "\n".join([
        screen_line("run", _PREPARATION_TITLE, indent=4),
        screen_field(
            "info", "Sample sheet", sample_sheet_display,
            indent=8, label_width=20,
        ),
        screen_field(
            "count", "Datasets selected", dataset_count,
            indent=8, label_width=20,
        ),
        screen_field(
            "info", "Resource directory", resource_directory,
            indent=8, label_width=20,
        ),
        screen_field(
            "info", "Output directory", output_directory,
            indent=8, label_width=20,
        ),
        screen_field(
            "info", "Current stage", _PREPARATION_STAGE,
            indent=8, label_width=20,
        ),
    ])


def _dataset_start_block(
    *,
    dataset_index: int,
    dataset_count: int,
    dataset_id: str,
    input_file: str | Path,
) -> str:
    """Render the identity and progress position of the next dataset."""
    return "\n".join([
        screen_line(
            "run",
            "Starting dataset %d of %d — %s"
            % (dataset_index, dataset_count, dataset_id),
            indent=4,
        ),
        screen_field(
            "info", "Summary statistics", input_file,
            indent=8, label_width=20,
        ),
    ])


def _write_run_metadata(config, args, row: HarmonisationSampleSheetRow) -> Path:
    root = Path(config.run.output_directory).expanduser().resolve()
    metadata = _runtime_path(
        config, root, "metadata_directory", dataset_id=row.dataset_id,
    )
    metadata.mkdir(parents=True, exist_ok=True)
    write_resolved_configuration(
        config,
        _runtime_path(config, metadata, "resolved_config_file"),
        modules=("harmonisation", "qc_summary"),
        resource_paths=_RESOLVED_RESOURCE_PATHS,
    )
    _write_text(
        _runtime_path(config, metadata, "sample_sheet_row_file"),
        json.dumps(row.model_dump(mode="json"), indent=2, sort_keys=True) + "\n",
    )
    _write_text(
        _runtime_path(config, metadata, "command_file"),
        " ".join(sys.argv) + "\n",
    )
    if hasattr(args, "run_config"):
        shutil.copyfile(
            args.run_config,
            _runtime_path(config, metadata, "supplied_config_file"),
        )
    return metadata


def _run_validated_rows(
    args,
    config,
    rows: list[HarmonisationSampleSheetRow],
    resolved_executables=None,
    screen_router: _DatasetScreenRouter | None = None,
):
    if screen_router is None:
        report_paths = _screen_report_paths(config, rows)
        with _DatasetScreenRouter(
            report_paths,
            display=config.logging.show_screen,
        ) as managed_router:
            return _run_validated_rows(
                args,
                config,
                rows,
                resolved_executables=resolved_executables,
                screen_router=managed_router,
            )

    fixed_info_from_cli = _fixed_info_cli_value(args, config)
    _validate_info_fallback(rows, fixed_info_from_cli)
    resource_root = config.resources.root
    if resource_root is None:
        raise ConfigurationError(
            "resources.root is required; supply --resource-directory or resources.root in --run-config"
        )
    resource_root = Path(resource_root).expanduser().resolve()
    if not resource_root.is_dir():
        raise ConfigurationError("Resource directory does not exist or is not a directory: %s" % resource_root)
    output_root = Path(config.run.output_directory).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    top_log = _top_run_log(config, output_root)
    sample_sheet = getattr(args, "sample_sheet", None)
    sample_sheet_display = _sample_sheet_display(sample_sheet)
    preparation = _preparation_block(
        sample_sheet=sample_sheet,
        dataset_count=len(rows),
        resource_directory=resource_root,
        output_directory=output_root,
    )
    print("\n" + preparation)
    _append_run_log(
        top_log,
        "INFO",
        "%s: sample_sheet=%s; datasets_selected=%d; "
        "resource_directory=%s; output_directory=%s; current_stage=%s."
        % (
            _PREPARATION_TITLE,
            sample_sheet_display,
            len(rows),
            resource_root,
            output_root,
            _PREPARATION_STAGE.lower(),
        ),
    )
    top_metadata = _runtime_path(
        config, output_root, "top_metadata_directory",
    )
    write_resolved_configuration(
        config,
        _runtime_path(config, top_metadata, "resolved_config_file"),
        modules=("harmonisation", "qc_summary"),
        resource_paths=_RESOLVED_RESOURCE_PATHS,
    )
    run_summary_path = _runtime_path(
        config, top_metadata, "run_summary_file",
    )
    html_report_paths, run_html_report_path = _html_report_paths(config, rows)
    run_summary_records = {}
    run_summary_manifests: dict[str, dict[str, Any]] = {}
    for row in rows:
        record = _initial_run_summary_record(row)
        record["screen_report"] = str(screen_router.report_path(row.dataset_id))
        record["html_report"] = str(html_report_paths[row.dataset_id])
        run_summary_records[row.dataset_id] = record
        run_summary_manifests[row.dataset_id] = {
            "report_context": _html_report_context(
                config, row, sample_sheet=sample_sheet,
            ),
        }
    _write_run_summary(run_summary_path, list(run_summary_records.values()))
    _write_html_reports(
        run_path=run_html_report_path,
        records=run_summary_records,
        manifests=run_summary_manifests,
        dataset_paths=html_report_paths,
    )

    defaults = _engine_defaults(config, resolved_executables=resolved_executables)
    results: dict[str, Any] = {}
    failures = []
    dataset_count = len(rows)
    strand_reference = config.modules.harmonisation.default_eaf
    comparison_reference = config.modules.harmonisation.comparison_af

    def record_dataset(
        row: HarmonisationSampleSheetRow,
        status: str,
        manifest_path: Path,
        *,
        manifest: dict[str, Any] | None = None,
        failure: Any = None,
    ) -> None:
        """Persist one dataset state while retaining all other selected rows."""
        manifest_document = (
            manifest
            if manifest is not None
            else _read_summary_manifest(manifest_path, required=False)
        )
        report_context = run_summary_manifests[row.dataset_id]["report_context"]
        record = _run_summary_record(
            row,
            status=status,
            manifest_path=manifest_path,
            manifest=manifest_document,
            strand_reference_panel=strand_reference.source,
            strand_reference_population=strand_reference.column,
            comparison_af_panel=comparison_reference.source,
            comparison_af_population=comparison_reference.column,
            failure=failure,
            report_context=report_context,
        )
        record["screen_report"] = str(screen_router.report_path(row.dataset_id))
        record["html_report"] = str(html_report_paths[row.dataset_id])
        run_summary_records[row.dataset_id] = record
        run_summary_manifests[row.dataset_id] = {
            **manifest_document,
            "report_context": report_context,
        }
        _write_run_summary(run_summary_path, list(run_summary_records.values()))
        _write_html_reports(
            run_path=run_html_report_path,
            records=run_summary_records,
            manifests=run_summary_manifests,
            dataset_paths=html_report_paths,
            changed_dataset=row.dataset_id,
        )

    for dataset_index, row in enumerate(rows, start=1):
        screen_router.select(row.dataset_id)
        run_summary_records[row.dataset_id]["status"] = "RUNNING"
        _write_run_summary(run_summary_path, list(run_summary_records.values()))
        _write_html_reports(
            run_path=run_html_report_path,
            records=run_summary_records,
            manifests=run_summary_manifests,
            dataset_paths=html_report_paths,
            changed_dataset=row.dataset_id,
        )
        dataset_start = _dataset_start_block(
            dataset_index=dataset_index,
            dataset_count=dataset_count,
            dataset_id=row.dataset_id,
            input_file=row.input_file,
        )
        print("\n" + dataset_start)
        start_log_message = (
            "Starting dataset %d of %d: dataset_id=%s; input_file=%s."
            % (dataset_index, dataset_count, row.dataset_id, row.input_file)
        )
        _append_run_log(top_log, "INFO", start_log_message)
        metadata = _runtime_path(
            config, output_root, "metadata_directory", dataset_id=row.dataset_id,
        )
        dataset_log = metadata / config.logging.filename
        dataset_output = configured_output_path(
            output_root,
            defaults["output_layout"]["dataset_directory"],
            error_type=RunSummaryError,
            dataset_id=row.dataset_id,
        )
        manifest_path = configured_output_path(
            dataset_output,
            defaults["output_layout"]["run_manifest"],
            error_type=RunSummaryError,
            dataset_id=row.dataset_id,
        )
        active_manifest_path = manifest_path
        try:
            _write_run_metadata(config, args, row)
            for warning in row.normalisation_warnings:
                _append_run_log(dataset_log, "WARNING", warning)
                print(screen_field(
                    "warning", "Sample-sheet value", warning,
                    indent=4, label_width=20,
                ))
            engine_input = to_harmonisation_input(
                row,
                resource_directory=resource_root,
                output_directory=output_root,
                output_layout=defaults["output_layout"],
                fixed_info=fixed_info_from_cli,
                fixed_info_column=config.modules.harmonisation.fixed_info.column,
            )
            if (
                fixed_info_from_cli is not None
                and engine_input["info_source"] != "fixed_cli"
            ):
                _append_run_log(
                    dataset_log,
                    "INFO",
                    "--fixed-info %g was not used for dataset %s because %s has "
                    "priority."
                    % (
                        fixed_info_from_cli,
                        row.dataset_id,
                        engine_input["info_source_detail"],
                    ),
                )
            elif engine_input["info_source"] == "fixed_cli":
                message = (
                    "Dataset %s assigns INFO=%g to every variant from the explicit "
                    "--fixed-info option. This is a user-assigned constant, not a "
                    "variant-level measured imputation-quality score."
                    % (row.dataset_id, fixed_info_from_cli)
                )
                _append_run_log(dataset_log, "WARNING", message)
                print(screen_field(
                    "warning", "Fixed INFO", message,
                    indent=4, label_width=20,
                ))
        except KeyboardInterrupt:
            record_dataset(
                row, "INTERRUPTED", manifest_path,
                failure="Run interrupted by the user during dataset setup.",
            )
            _append_run_log(
                top_log,
                "WARNING",
                "Run interrupted by the user while preparing dataset %s."
                % row.dataset_id,
            )
            raise
        except (RunSummaryError, ScreenReportError):
            raise
        except Exception as exc:
            record_dataset(row, "FAILED", manifest_path, failure=exc)
            _append_run_log(
                dataset_log,
                "ERROR",
                "Dataset setup failed: %s: %s"
                % (type(exc).__name__, exc),
            )
            _append_run_log(
                top_log,
                "ERROR",
                "Dataset %s failed before analysis; see %s."
                % (row.dataset_id, dataset_log),
            )
            failures.append((row.dataset_id, exc))
            continue
        _append_run_log(dataset_log, "INFO", start_log_message)
        attempts = config.execution.retries + 1
        for attempt in range(1, attempts + 1):
            try:
                result = run_harmonisation_pipeline(
                    sample_column_dict=engine_input,
                    default_cfg=defaults,
                    threads=config.execution.threads,
                )
                validation_manifest = result.get("manifest") or str(manifest_path)
                reported_manifest = Path(validation_manifest).expanduser().resolve()
                active_manifest_path = reported_manifest
                manifest_document = _read_summary_manifest(
                    reported_manifest, required=True,
                )
                if config.modules.harmonisation.concordance_validation.enabled:
                    from postgwas.modules.harmonisation.concordance.service import run_concordance_validation

                    dataset_context = manifest_document.get("dataset") or {}
                    build_context = dataset_context.get("genome_build") or {}
                    input_build = build_context.get("inferred_build")
                    supported_builds = set(
                        defaults["vcf_processing"]["target_builds"]
                    )
                    if input_build not in supported_builds:
                        raise ConcordanceValidationError(
                            "Harmonisation result did not record a resolved input genome build."
                        )
                    validation = run_concordance_validation(
                        row=row,
                        vcf_path=result[input_build],
                        output_root=output_root,
                        settings=config.modules.harmonisation.concordance_validation,
                        policies=load_policies(config.modules.harmonisation.policies),
                        threads=config.execution.threads,
                        bcftools=config.resources.executables.bcftools,
                        output_layout=defaults["output_layout"],
                        vcf_config=defaults["vcf_processing"],
                        file_log_level=config.logging.file_level,
                        screen_log_level=config.logging.console_level,
                        run_manifest=validation_manifest,
                        external_eaf_mapping=config.modules.harmonisation.external_eaf_mapping,
                    )
                    result["concordance_validation"] = validation
                    if validation["status"] == "FAIL":
                        raise ConcordanceValidationError(
                            "Concordance validation failed for dataset %s; see %s"
                            % (row.dataset_id, validation["reports"]["summary"])
                        )
                    manifest_document = _read_summary_manifest(
                        reported_manifest, required=True,
                    )
                dataset_status = str(
                    result.get("status")
                    or manifest_document.get("status")
                    or "OK"
                ).upper()
                if dataset_status not in {"OK", "PARTIAL"}:
                    raise PipelineError(
                        "Dataset %s returned unsupported terminal status %r; "
                        "expected OK or PARTIAL."
                        % (row.dataset_id, dataset_status)
                    )
                record_dataset(
                    row,
                    dataset_status,
                    reported_manifest,
                    manifest=manifest_document,
                )
                result["screen_report"] = str(
                    screen_router.report_path(row.dataset_id)
                )
                result["html_report"] = str(html_report_paths[row.dataset_id])
                results[row.dataset_id] = result
                if dataset_status == "OK":
                    _append_run_log(
                        dataset_log, "INFO", "Dataset completed successfully.",
                    )
                    _append_run_log(
                        top_log,
                        "INFO",
                        "Dataset %s completed successfully." % row.dataset_id,
                    )
                    print(screen_line(
                        "success", "Dataset %s completed" % row.dataset_id,
                        indent=4,
                    ))
                else:
                    _append_run_log(
                        dataset_log,
                        "WARNING",
                        "Dataset completed with status %s." % dataset_status,
                    )
                    _append_run_log(
                        top_log,
                        "WARNING",
                        "Dataset %s completed with status %s."
                        % (row.dataset_id, dataset_status),
                    )
                    print(screen_line(
                        "warning",
                        "Dataset %s completed with status %s"
                        % (row.dataset_id, dataset_status),
                        indent=4,
                    ))
                break
            except ScreenReportError:
                raise
            except RunSummaryError as exc:
                try:
                    record_dataset(
                        row, "FAILED", active_manifest_path, failure=exc,
                    )
                except RunSummaryError as update_exc:
                    _append_run_log(
                        top_log,
                        "ERROR",
                        "The dataset failure also could not be saved to the required "
                        "run reports: %s" % update_exc,
                    )
                _append_run_log(
                    dataset_log, "ERROR", "Required reporting failure: %s" % exc,
                )
                _append_run_log(
                    top_log, "ERROR", "Required reporting failure: %s" % exc,
                )
                raise
            except KeyboardInterrupt:
                try:
                    record_dataset(
                        row, "INTERRUPTED", active_manifest_path,
                        failure="Run interrupted by the user.",
                    )
                except RunSummaryError as summary_exc:
                    _append_run_log(
                        top_log, "ERROR",
                        "Run interruption was recorded, but the run summary could "
                        "not be updated: %s" % summary_exc,
                    )
                _append_run_log(dataset_log, "WARNING", "Dataset interrupted by the user.")
                _append_run_log(
                    top_log,
                    "WARNING",
                    "Run interrupted by the user while processing dataset %s."
                    % row.dataset_id,
                )
                raise
            except ConfigError as exc:
                record_dataset(row, "FAILED", active_manifest_path, failure=exc)
                _append_run_log(dataset_log, "ERROR", "Configuration failure: %s" % exc)
                _append_run_log(
                    top_log,
                    "ERROR",
                    "Dataset %s failed configuration validation; see %s."
                    % (row.dataset_id, dataset_log),
                )
                failures.append((row.dataset_id, exc))
                break
            except ConcordanceValidationError as exc:
                record_dataset(
                    row, "CONCORDANCE_FAILED", active_manifest_path, failure=exc,
                )
                _append_run_log(
                    dataset_log,
                    "ERROR",
                    "Concordance validation failed after harmonisation: %s" % exc,
                )
                _append_run_log(
                    top_log,
                    "ERROR",
                    "Dataset %s harmonised but failed concordance validation; see %s."
                    % (row.dataset_id, dataset_log),
                )
                failures.append((row.dataset_id, exc))
                break
            except Exception as exc:
                _append_run_log(
                    dataset_log,
                    "ERROR",
                    "Attempt %d/%d failed: %s: %s"
                    % (attempt, attempts, type(exc).__name__, exc),
                )
                if attempt == attempts:
                    record_dataset(
                        row, "FAILED", active_manifest_path, failure=exc,
                    )
                    _append_run_log(
                        dataset_log,
                        "ERROR",
                        "Dataset failed permanently.\n%s" % traceback.format_exc(),
                    )
                    _append_run_log(
                        top_log,
                        "ERROR",
                        "Dataset %s failed after %d attempt(s); see %s."
                        % (row.dataset_id, attempts, dataset_log),
                    )
                    failures.append((row.dataset_id, exc))
    screen_router.select(None)
    summary_block = _run_summary_block(
        run_summary_path,
        run_html_report_path,
        list(run_summary_records.values()),
    )
    print("\n" + summary_block)
    _append_run_log(
        top_log,
        "INFO",
        "Harmonisation run summary: %s" % run_summary_path,
    )
    _append_run_log(
        top_log,
        "INFO",
        "Harmonisation HTML run report: %s" % run_html_report_path,
    )
    ok_count = sum(
        str(record.get("status") or "").upper() == "OK"
        for record in run_summary_records.values()
    )
    partial_count = sum(
        str(record.get("status") or "").upper() == "PARTIAL"
        for record in run_summary_records.values()
    )
    if failures and not results:
        _append_run_log(
            top_log,
            "ERROR",
            "Run failed: all %d dataset(s) failed." % len(failures),
        )
        raise PipelineError("All harmonisation datasets failed; see %s" % top_log)
    if failures:
        _append_run_log(
            top_log,
            "WARNING",
            "Run completed with %d OK, %d partial, and %d failed dataset(s)."
            % (ok_count, partial_count, len(failures)),
        )
        print(screen_field(
            "warning", "Datasets failed", "%d; see %s" % (len(failures), top_log),
            indent=4, label_width=18,
        ))
        if partial_count:
            print(screen_field(
                "warning", "Datasets partial", partial_count,
                indent=4, label_width=18,
            ))
    elif partial_count:
        _append_run_log(
            top_log,
            "WARNING",
            "Run completed with %d OK and %d partial dataset(s)."
            % (ok_count, partial_count),
        )
        print(screen_field(
            "warning", "Datasets partial", partial_count,
            indent=4, label_width=18,
        ))
    else:
        _append_run_log(
            top_log,
            "INFO",
            "Run completed successfully: %d dataset(s)." % ok_count,
        )
    print(screen_line(
        "success", "Harmonisation complete — %d dataset(s)" % len(results),
        indent=4,
    ))
    return results


def run_harmonisation(
    args,
    *,
    configuration=None,
    resolved_executables=None,
):
    if not hasattr(args, "sample_sheet"):
        raise ConfigurationError("--sample-sheet is required")

    overrides = explicit_overrides(args, CLI_OVERRIDE_PATHS)
    config = configuration or load_configuration(
        getattr(args, "run_config", None), cli_overrides=overrides,
    )
    preflight_log = _top_run_log(config, config.run.output_directory)
    rows: list[HarmonisationSampleSheetRow] = []
    try:
        rows = load_harmonisation_sample_sheet(args.sample_sheet)
        if hasattr(args, "dataset_id"):
            rows = [row for row in rows if row.dataset_id == args.dataset_id]
            if not rows:
                raise ConfigurationError(
                    "dataset_id %r is not present in sample sheet %s"
                    % (args.dataset_id, args.sample_sheet)
                )
        _validate_info_fallback(rows, _fixed_info_cli_value(args, config))
        if resolved_executables is None:
            resolved_executables = require_binaries(
                harmonisation_executable_requirements(config),
                plugins=(
                    config.modules.harmonisation.vcf_processing.liftover_plugin,
                ),
            )
    except Exception as exc:
        _append_run_log(preflight_log, "ERROR", "Harmonisation preflight failed: %s" % exc)
        if rows:
            try:
                summary_path = _write_preflight_failure_run_summary(
                    config, rows, exc, sample_sheet=args.sample_sheet,
                )
                _append_run_log(
                    preflight_log,
                    "INFO",
                    "Preflight-failure dataset summary: %s" % summary_path,
                )
            except (RunSummaryError, ScreenReportError) as summary_exc:
                _append_run_log(
                    preflight_log,
                    "ERROR",
                    "The preflight failure could not be written to the required "
                    "run reports: %s" % summary_exc,
                )
        raise
    return _run_validated_rows(
        args, config, rows, resolved_executables=resolved_executables,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="postgwas harmonisation",
        usage="postgwas harmonisation --sample-sheet PATH [options]",
        description=(
            "Prepare one or more GWAS summary-statistics files for PostGWAS. "
            "The command checks the sample sheet, harmonises alleles and genome "
            "coordinates, and writes analysis-ready GWAS-VCF files."
        ),
        epilog=format_cli_examples(
            (
                "Harmonise every dataset in a sample sheet:",
                "postgwas harmonisation",
                (
                    "--sample-sheet studies.csv",
                    "--run-config harmonisation.yaml",
                    "--resource-directory resources",
                    "--output-directory results",
                ),
            ),
            (
                "Harmonise one selected dataset and validate its merged VCF:",
                "postgwas harmonisation",
                (
                    "--sample-sheet studies.csv",
                    "--dataset-id STUDY",
                    "--resource-directory resources",
                    "--output-directory results",
                    "--validate",
                ),
            ),
            (
                "Export a reusable harmonisation configuration:",
                "postgwas config export",
                (
                    "--module harmonisation",
                    "--style full",
                    "--output harmonisation.yaml",
                ),
            ),
            notes=(
                "Sample-sheet templates are available in examples/configs/harmonisation/.",
            ),
        ),
        parents=[get_harmonisation_parser()],
        formatter_class=AlignedRichHelpFormatter,
    )
    if len(sys.argv) == 1:
        parser.print_help()
        return 0
    args = parser.parse_args()
    try:
        run_harmonisation(args)
    except KeyboardInterrupt:
        try:
            fallback = load_configuration()
            output = Path(getattr(args, "output_directory", fallback.run.output_directory))
            interruption_log = _top_run_log(fallback, output)
            _append_run_log(
                interruption_log,
                "WARNING",
                "Harmonisation interrupted by the user; the run did not complete.",
            )
        except Exception:
            interruption_log = None
        suffix = "; see %s" % interruption_log if interruption_log else ""
        print(screen_field(
            "warning", "Interrupted", "harmonisation stopped%s" % suffix,
            indent=4, label_width=18,
        ), file=sys.stderr)
        return 130
    except (ConfigurationError, ConfigError) as exc:
        try:
            fallback = load_configuration()
            output = Path(getattr(args, "output_directory", fallback.run.output_directory))
            _append_run_log(
                _top_run_log(fallback, output),
                "ERROR",
                "Configuration error: %s" % exc,
            )
        except Exception:
            pass
        print(screen_field(
            "error", "Configuration", exc,
            indent=4, label_width=18,
        ), file=sys.stderr)
        return 2
    except Exception as exc:
        try:
            fallback = load_configuration()
            output = Path(getattr(args, "output_directory", fallback.run.output_directory))
            failure_log = _top_run_log(fallback, output)
            _append_run_log(
                failure_log,
                "ERROR",
                "Harmonisation failed: %s: %s\n%s"
                % (type(exc).__name__, exc, traceback.format_exc()),
            )
        except Exception:
            failure_log = None
        suffix = (
            "; see %s" % failure_log
            if failure_log and str(failure_log) not in str(exc)
            else ""
        )
        print(screen_field(
            "error", "Harmonisation", "%s%s" % (exc, suffix),
            indent=4, label_width=18,
        ), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
