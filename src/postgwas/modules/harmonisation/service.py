"""PostGWAS — harmonisation orchestration.

This module owns the two loops the pipeline is built from:

* ``run_harmonisation_pipeline``  — one dataset, eight dataset-level steps and
  then everything after the per-chromosome merge.
* ``harmonise_chromosomes``       — the per-chromosome fan-out, with retry
  rounds, and ``process_one_chromosome`` — the sixteen steps a single
  chromosome goes through inside a worker process.

Three rules shape the code, all from plan v3:

1. **Nothing in library code calls ``sys.exit``.**  ``SystemExit`` derives from
   ``BaseException``, so ``except Exception`` in ``cli.py`` cannot see it and a
   single bad dataset used to kill an entire multi-dataset batch.  Everything
   that used to exit now raises :class:`ConfigError` (a bad configuration) or
   :class:`PipelineError` (a run that could not be completed), and ``cli.py``
   decides what the process exit status should be.

2. **Workers never write to stdout.**  Each chromosome worker buffers its
   screen text in its own :class:`~postgwas.core.pipeline_logging.PipelineLogger`
   and returns it; the parent prints one complete block per chromosome from its
   single-threaded ``as_completed`` loop, so interleaving is structurally
   impossible.  This is why ``print_lock`` and ``safe_print`` are gone: under
   ``ProcessPoolExecutor(spawn)`` every process gets its **own** copy of a
   ``threading.Lock``, so the old "process-safe print" serialised nothing at all.

3. **Every number is a policy.**  The resolved policy set is loaded from the
   canonical YAML and written to the dataset log and run manifest before
   anything runs.
"""

import os
import sys
import json
import time
import argparse
import platform
import traceback
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from concurrent.futures import Future, as_completed
from concurrent.futures.process import BrokenProcessPool
from multiprocessing import get_context

import polars as pl
import pandas as pd


# ----------------------------------------------------------------------
# PostGWAS imports
# ----------------------------------------------------------------------
from postgwas.core.execution.runtime import safe_thread_count, validate_path
from postgwas.core.execution.scheduling import ResourceQueue
from postgwas.modules.harmonisation.scheduling import plan_chromosome_schedule
from postgwas.core.polars_runtime import (
    POLARS_THREAD_ENVIRONMENT_VARIABLE as _POLARS_THREAD_ENVIRONMENT_VARIABLE,
    bounded_polars_executor as _bounded_chromosome_executor,
    initialise_polars_worker as _initialise_chromosome_worker,
)

from postgwas.modules.harmonisation.policies import (
    load_policies,
    default_policies,
    PolicyError,
)
from postgwas.core.pipeline_logging import (
    PipelineLogger,
    print_chromosome_summary,
)
from postgwas.core.ui import MeasuredProgress, StageProgress
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.paths import configured_output_matches, configured_output_path
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    RejectOutputError,
    concat_reject_files,
    reconcile,
    resolved_reject_output_path,
    REASONS,
    REASON_STEPS,
)
from postgwas.modules.harmonisation.summary_statistics_io import (
    finalise_field_completeness_report,
    read_summary_statistics,
    inspect_summary_statistics_file,
    write_field_completeness_report,
)
from postgwas.modules.harmonisation.resource_paths import resolve_resource_file
from postgwas.modules.harmonisation.external_reference_staging import (
    ExternalReferenceStagingError,
    stage_shared_external_reference_files,
)
from postgwas.modules.harmonisation.input_validation import (
    validate_config,
    validate_header,
    validate_content,
    format_problems,
    ConfigError,
)
from postgwas.modules.harmonisation.chromosome_partition import write_chromosome_partitions
from postgwas.modules.harmonisation.genome_build import infer_genome_build
from postgwas.modules.harmonisation.strand import resolve_strand_consensus
from postgwas.modules.harmonisation.study_properties import (
    finalise_eaf_decision_from_chromosomes,
    resolve_study_properties,
)
from postgwas.modules.harmonisation.allele_frequency import harmonise_allele_frequency
from postgwas.modules.harmonisation.sample_size import (
    harmonise_sample_sizes,
    prepare_missing_sample_sizes,
)
from postgwas.modules.harmonisation.effect_type import harmonise_effect_estimates
from postgwas.modules.harmonisation.effect_from_z import (
    XChromosomeZOnlyReconstructionError,
    derive_effect_and_standard_error_from_z,
    enforce_x_chromosome_z_only_policy,
)
from postgwas.modules.harmonisation.z_score import derive_z_score_from_effect_and_standard_error
from postgwas.modules.harmonisation.p_values import (
    harmonise_p_values,
)
from postgwas.modules.harmonisation.standard_error import (
    derive_standard_error_from_effect_and_p_value,
)
from postgwas.modules.harmonisation.effect_validation import (
    validate_effect_statistics,
    final_completeness_check,
)
from postgwas.modules.harmonisation.imputation_quality import (
    harmonise_imputation_quality,
    resolve_dataset_info_score_type,
)
from postgwas.modules.harmonisation.variant_identifiers import harmonise_variant_identifiers
from postgwas.modules.harmonisation.vcf_processing import (
    annotate_and_liftover_vcf,
    concat_vcfs_by_build,
    require_binaries,
)
from postgwas.modules.harmonisation.vcf_provenance import (
    build_harmonisation_vcf_provenance,
)
from postgwas.modules.harmonisation.population_frequency import (
    PopulationFrequencyQCError,
    run_population_frequency_qc,
)
from postgwas.modules.harmonisation.gwas2vcf_export import (
    export_gwas2vcf_input,
)
from postgwas.modules.harmonisation.gwas2vcf_runner import run_gwas2vcf
from postgwas.modules.harmonisation.cleanup import (
    finalise_harmonisation_outputs,
    remove_merged_gwas2vcf_intermediate,
    remove_partial_chromosome_outputs,
)
from postgwas.modules.harmonisation.resource_preflight import (
    ResourcePreflightError,
    recheck_preflighted_resource_map,
    validate_harmonisation_resource_maps,
)
from postgwas.modules.harmonisation.qc_results import qc_results_to_dataframe
from postgwas.core.values import optional_text

from postgwas.config.models.modules.qc_summary import QCSummaryConfig
from postgwas.modules.qc_summary.service import run_qc_assessment
from postgwas.modules.harmonisation.qc_reporting import (
    harmonisation_qc_summary_lines,
    harmonisation_qc_takeaway_lines,
)


# ===============================================================
# Typed exceptions  (plan Part 7 prerequisite)
# ===============================================================

class PipelineError(RuntimeError):
    """A run could not be completed.

    Raised instead of ``sys.exit()`` from anywhere inside the pipeline, so that
    ``except Exception`` in ``cli.py`` sees it and the remaining datasets in a
    batch still run.  Derived from ``RuntimeError`` so that any caller which
    already catches ``RuntimeError`` keeps working.

    ``results``     whatever QC was collected before the failure, so the caller
                    can still write the QC summary for a partial run.
    ``chromosomes`` the chromosomes that did not complete.
    ``status``      ``FAILED`` (nothing usable) or ``PARTIAL`` (some
                    chromosomes completed).
    """

    def __init__(self, message, results=None, chromosomes=None, status="FAILED"):
        # type: (str, Optional[Dict[str, Any]], Optional[List[str]], str) -> None
        super(PipelineError, self).__init__(message)
        self.results = results
        self.chromosomes = list(chromosomes or [])
        self.status = status


#: Re-exported so ``cli.py`` can import both typed exceptions from one place.
__all__ = [
    "ConfigError",
    "PipelineError",
    "run_harmonisation_pipeline",
    "harmonise_chromosomes",
    "process_one_chromosome",
    "build_resource_map",
    "preflight_harmonisation_resources",
]


# ===============================================================
# Small helpers
# ===============================================================

#: The sixteen per-chromosome steps. Coordinates are already harmonised before
#: the validated dataset is split, so workers do not repeat that dataset step.
CHR_STEP_TOTAL = 16

#: The eight dataset-level steps of plan Part 1.
DATASET_STEP_TOTAL = 8

#: The five post-merge steps: merge, population AF QC, cleanup, virtual QC,
#: and final reports.
POST_MERGE_STEP_TOTAL = 5


class _HarmonisationProgress:
    """Route validated harmonisation boundaries through shared progress UI.

    Dataset and post-merge steps have independent numbering, while chromosome
    work is parallel and retryable. Keeping those as three shared progress
    objects prevents a reset from 8/8 to 1/5 and lets only the parent process
    update the measured chromosome count.
    """

    def __init__(self, *, enabled, outcome_label_width, console=None):
        self.enabled = bool(enabled)
        self.outcome_label_width = int(outcome_label_width)
        self.console = console
        self.dataset_stages = StageProgress(
            "Harmonisation dataset stages",
            enabled=self.enabled,
            outcome_label_width=self.outcome_label_width,
            console=self.console,
        )
        self._chromosome_progress = None
        self._chromosomes = ()
        self._completed_chromosomes = set()
        self._chromosome_phase_finished = False
        self._post_merge_stages = None

    def start_chromosomes(self, chromosomes):
        ordered = tuple(str(value) for value in chromosomes)
        if not ordered:
            raise ValueError("Chromosome progress requires at least one chromosome")
        if len(ordered) != len(set(ordered)):
            raise ValueError("Chromosome progress received duplicate chromosomes")
        if self._chromosome_progress is not None:
            raise RuntimeError("Chromosome progress has already started")
        self._chromosomes = ordered
        self._chromosome_progress = MeasuredProgress(
            "Harmonisation chromosome progress",
            enabled=self.enabled,
            console=self.console,
        )
        self._chromosome_progress.start(
            "Complete chromosome harmonisation",
            total=len(ordered),
        )

    def record_chromosome_result(self, chromosome, status, attempt):
        if self._chromosome_progress is None:
            return
        chromosome = str(chromosome)
        if chromosome not in self._chromosomes:
            raise ValueError(
                "Chromosome progress received an unexpected chromosome: %s"
                % chromosome
            )
        if str(status).lower() == "ok":
            self._completed_chromosomes.add(chromosome)
            title = "Chromosome %s completed on attempt %d" % (
                chromosome,
                int(attempt),
            )
        else:
            title = "Chromosome %s failed on attempt %d" % (
                chromosome,
                int(attempt),
            )
            self._chromosome_progress.set_phase(title)
            return
        self._chromosome_progress.update(
            len(self._completed_chromosomes),
            total=len(self._chromosomes),
            title=title,
        )

    def finish_chromosomes(self, failed=()):
        if (
            self._chromosome_progress is None
            or self._chromosome_phase_finished
        ):
            return
        failed = tuple(str(value) for value in failed)
        if failed:
            self._chromosome_progress.fail(
                title="%d chromosome%s incomplete: %s"
                % (
                    len(failed),
                    "" if len(failed) == 1 else "s",
                    ", ".join(failed),
                )
            )
        else:
            total = len(self._chromosomes)
            if len(self._completed_chromosomes) != total:
                raise RuntimeError(
                    "Chromosome progress cannot complete: %d of %d validated "
                    "chromosomes returned success"
                    % (len(self._completed_chromosomes), total)
                )
            self._chromosome_progress.complete(
                total,
                total=total,
                title="All chromosome outputs and required counts validated",
            )
        self._chromosome_phase_finished = True

    def fail_active_chromosomes(self):
        if (
            self._chromosome_progress is None
            or self._chromosome_phase_finished
        ):
            return
        self._chromosome_progress.fail(
            title="Chromosome harmonisation stopped before validation completed"
        )
        self._chromosome_phase_finished = True

    def start_post_merge(self):
        if self._post_merge_stages is not None:
            raise RuntimeError("Post-merge progress has already started")
        self._post_merge_stages = StageProgress(
            "Harmonisation post-merge stages",
            enabled=self.enabled,
            outcome_label_width=self.outcome_label_width,
            console=self.console,
        )
        return self._post_merge_stages

    def close(self):
        self.dataset_stages.close()
        if self._chromosome_progress is not None:
            self._chromosome_progress.close()
        if self._post_merge_stages is not None:
            self._post_merge_stages.close()


@dataclass(frozen=True)
class _DatasetPreparation:
    """Validated, partitioned dataset state needed by chromosome fan-out.

    The full Polars frame deliberately is not part of this result.  Once the
    helper returns, the parent process can release that genome-wide frame
    before chromosome workers are spawned.
    """

    sample_columns: Dict[str, Any]
    chromosome_files: Dict[str, str]
    resource_maps: Dict[str, Dict[str, Any]]
    genome_build: str
    genome_build_info: Dict[str, Any]
    study_decisions: Dict[str, Any]
    dataset_study_decisions: Dict[str, Any]
    missing_sample_size_plan: Dict[str, Any]
    resource_preflight: Dict[str, Any]
    field_completeness: Dict[str, Any]
    field_completeness_path: Path
    input_variants: int
    rows_read: int
    removed_invalid_coordinates: int
    removed_unsupported_chromosomes: int
    removed_non_standard_alleles: int
    removed_duplicates: int
    removed_missing_values: int
    ready_variants: int
    ready_snps: int
    ready_indels_or_other: int
    partition_rows: Dict[str, int]


@dataclass(frozen=True)
class _PostMergeOutcome:
    """Validated outputs from the five post-merge stages."""

    primary_outputs: Dict[str, str]
    assessment: Dict[str, Any]
    population_frequency: Dict[str, Any]
    concatenated_vcfs: Dict[str, Any]
    pre_vcf_summary: Dict[str, Any]
    status: str
    combined_log: Optional[str]


def _announce(logger, message, marker="", indent=0):
    # type: (Any, Any, str, int) -> None
    """A dataset-level message: recorded in the dataset log AND shown on screen.

    Only the parent process calls this.  Chromosome workers buffer their text
    and the parent prints it as one block, so nothing can interleave.
    """
    if logger is not None:
        try:
            logger.log(marker, message, indent=indent)
        except Exception:
            pass
    sys.stdout.write("%s\n" % (message,))
    sys.stdout.flush()


def _study_decisions_block(study_decisions):
    """Render resolved study properties in a short, aligned screen summary."""
    frequency_type = (
        "MAF-like; reference confirmation required"
        if study_decisions.get("eaf_is_maf") is True
        else "effect allele frequency"
        if study_decisions.get("eaf_is_maf") is False
        else "not decided"
    )
    se_scale = (
        study_decisions.get("se_scale") or "not decided"
        if study_decisions.get("effect_type") == "odds_ratio"
        else "not applicable to beta"
    )
    return "\n".join([
        screen_line("decision", "Study-wide decisions", indent=4),
        screen_field(
            "analysis", "Effect type",
            study_decisions.get("effect_type") or "not decided",
            indent=8, label_width=16,
        ),
        screen_field(
            "analysis", "SE scale", se_scale,
            indent=8, label_width=16,
        ),
        screen_field(
            "analysis", "P-value type",
            study_decisions.get("pvalue_type") or "not decided",
            indent=8, label_width=16,
        ),
        screen_field(
            "genetic", "Frequency type", frequency_type,
            indent=8, label_width=16,
        ),
        screen_field(
            "genetic", "Strand consensus",
            study_decisions.get("strand") or "not decided",
            indent=8, label_width=16,
        ),
    ])


def _population_frequency_block(result):
    """Render the raw merged-VCF frequency comparison without hiding missingness."""
    lines = [screen_line("genetic", "Population-frequency similarity", indent=4)]
    status = str(result.get("status") or "unavailable")
    if result.get("raw_merged_vcf"):
        lines.append(screen_field(
            "info", "Raw merged VCF", result["raw_merged_vcf"],
            indent=8, label_width=24,
        ))
    if result.get("total_records") is not None:
        lines.append(screen_field(
            "count", "VCF records", _count(result["total_records"]),
            indent=8, label_width=24,
        ))
    if result.get("comparable_variants") is not None:
        lines.append(screen_field(
            "count", "Complete AF rows",
            _count(result["comparable_variants"]),
            indent=8, label_width=24,
        ))
    vcf_fields = result.get("vcf_fields") or {}
    for field, metrics in (result.get("fields") or {}).items():
        missing = metrics.get("missing")
        fraction = metrics.get("missing_fraction")
        value = (
            "%s (%.2f%%)" % (_count(missing), 100.0 * float(fraction))
            if fraction is not None else _count(missing)
        )
        lines.append(screen_field(
            "warning" if missing else "success",
            "Missing %s" % str(vcf_fields.get(field, field)).lstrip("%"),
            value,
            indent=8,
            label_width=24,
        ))
    for population, metrics in (result.get("populations") or {}).items():
        correlation = metrics.get("pearson_correlation")
        difference = metrics.get("mean_absolute_difference")
        inverted_difference = metrics.get("inverted_mean_absolute_difference")
        value = (
            "r=%s; mean |AF difference|=%s; after 1-AF=%s"
            % (
                "unavailable" if correlation is None else "%.6f" % correlation,
                "unavailable" if difference is None else "%.6f" % difference,
                (
                    "unavailable"
                    if inverted_difference is None
                    else "%.6f" % inverted_difference
                ),
            )
        )
        lines.append(screen_field(
            "analysis", population, value, indent=8, label_width=24,
        ))
    closest = result.get("closest_population")
    lines.append(screen_field(
        "decision" if closest else "warning",
        "Closest population",
        closest or "inconclusive — %s" % result.get("decision_reason", status),
        indent=8,
        label_width=24,
    ))
    selected = result.get("selected_population_check") or {}
    selected_status = selected.get("status")
    selected_column = selected.get("selected_column")
    selected_display = selected.get("normalized_population") or selected_column
    if selected_column:
        if selected_status == "match":
            selected_value = "%s matches closest population %s" % (
                selected_display, closest,
            )
        elif selected_status == "mismatch":
            selected_value = "%s does not match closest population %s" % (
                selected_display, closest,
            )
        elif selected_status == "not_compared":
            selected_value = "%s not compared because the result is inconclusive" % (
                selected_display,
            )
        else:
            selected_value = (
                "%s is not one of the configured populations; not compared"
                % selected_display
            )
        lines.append(screen_field(
            "success" if selected_status == "match" else
            "warning" if selected_status == "mismatch" else "info",
            "Selected AF column", selected_value,
            indent=8, label_width=24,
        ))
    for check in result.get("external_file_checks") or []:
        check_status = check.get("status")
        filename = check.get("filename")
        indicated = check.get("indicated_population")
        source = str(check.get("source") or "external")
        source_label = source[:1].upper() + source[1:]
        if check_status == "match":
            check_value = "%s indicates %s and matches closest population %s" % (
                filename, indicated, closest,
            )
        elif check_status == "mismatch":
            check_value = "%s indicates %s but closest population is %s" % (
                filename, indicated, closest,
            )
        elif check_status == "ambiguous":
            check_value = "%s has multiple population tokens after normalization" % (
                filename,
            )
        elif check_status == "not_compared":
            check_value = "%s not compared because the result is inconclusive" % (
                filename,
            )
        else:
            check_value = "%s has no configured population token after normalization" % (
                filename,
            )
        lines.append(screen_field(
            "success" if check_status == "match" else
            "warning" if check_status in {"mismatch", "ambiguous"} else "info",
            "%s filename" % source_label, check_value,
            indent=8, label_width=24,
        ))
    if result.get("report"):
        lines.append(screen_field(
            "info", "QC report", result["report"], indent=8, label_width=24,
        ))
    return "\n".join(lines)


def _validated_column_mapping_block(sample_columns, study_decisions, width=104):
    """Show every validated sample-sheet field using its public field name."""
    def display_value(key):
        value = sample_columns.get(key)
        if optional_text(value) is None:
            return "Not provided — will not be used"
        if str(value).strip().lower() == "auto":
            return "auto — not explicitly provided"
        return str(value)

    groups = (
        ("Dataset", "info", (
            ("dataset_id", "gwas_outputname"),
            ("input_file", "sumstat_file"),
            ("trait_type", "trait_type"),
            ("delimiter", "delimiter"),
        )),
        ("Variant columns", "genetic", (
            ("chromosome_column", "chr_col"),
            ("position_column", "pos_col"),
            ("chromosome_position_column", "chr_pos_col"),
            ("variant_id_column", "snp_id_col"),
            ("effect_allele_column", "ea_col"),
            ("other_allele_column", "oa_col"),
            ("effect_allele_frequency_column", "eaf_col"),
        )),
        ("Association columns", "analysis", (
            ("effect_column", "beta_or_col"),
            ("effect_type", "declared_effect_type"),
            ("standard_error_column", "se_col"),
            ("z_score_column", "imp_z_col"),
            ("p_value_column", "pval_col"),
            ("p_value_type", "declared_pvalue_type"),
        )),
        ("Sample size", "count", (
            ("control_count_column", "ncontrol_col"),
            ("case_count_column", "ncase_col"),
            ("control_count", "ncontrol"),
            ("case_count", "ncase"),
        )),
        ("Quality and external sources", "analysis", (
            ("resolved_info_source", "info_source_detail"),
            ("imputation_info_column", "imp_info_col"),
            ("external_info_file", "infofile"),
            ("external_info_column", "infocolumn"),
            ("fixed_info", "fixed_info"),
            ("external_eaf_file", "eaffile"),
            ("external_eaf_column", "eafcolumn"),
        )),
    )

    lines = [screen_line("analysis", "Validated sample-sheet values", indent=4)]
    label_width = max(len(label) for _, _, entries in groups for label, _ in entries)
    for heading, kind, entries in groups:
        lines.extend(["", screen_line(kind, heading, indent=8)])
        for label, key in entries:
            value = display_value(key)
            lines.extend(screen_field(
                "info", label, value,
                width=max(60, int(width)), indent=10, label_width=label_width,
            ).splitlines())
    return "\n".join(lines)


def _count(value):
    # type: (Any) -> str
    try:
        return "{:,}".format(int(value))
    except (TypeError, ValueError):
        return str(value)


def _field_completeness_block(report):
    """Keep original missingness separate from later numeric conversion counts."""
    labels = {
        "chr": "Chromosome", "pos": "Position", "ea": "Effect allele",
        "oa": "Other allele", "pval": "P-value", "beta": "Effect estimate",
        "zscore": "Z score", "se": "Standard error", "eaf": "Allele frequency",
        "info": "Imputation score", "snp": "Variant ID", "n": "Sample count",
        "ncase": "Case count", "ncontrol": "Control count",
    }

    def counted(count, total):
        if not total:
            return "%s / %s (percentage unavailable)" % (_count(count), _count(total))
        percentage = 100.0 * count / total
        percent_text = "<0.01" if 0 < percentage < 0.01 else "%.2f" % percentage
        return "%s / %s (%s%%)" % (_count(count), _count(total), percent_text)

    def detail(label, value):
        return screen_field("info", label, value, indent=12, label_width=22).splitlines()

    lines = [
        screen_line("analysis", "Input completeness and numeric conversion", indent=4),
        screen_field("count", "Variants read", _count(report.get("input_rows")), indent=8, label_width=22),
        screen_field(
            "info", "Count timing",
            "Missing: all input rows after missing-token parsing. Numeric conversion: "
            "rows remaining after initial missing-field and coordinate filtering.",
            indent=8, label_width=22,
        ),
    ]
    for field in report.get("fields") or []:
        name = field.get("field")
        missing = field.get("post_parse_missing")
        column = field.get("input_column")
        conversion = (report.get("numeric_conversion_by_column") or {}).get(column) or {}
        source = column or "Not supplied"
        if field.get("source_kind") == "configured_column_absent":
            source = field["source_label"]
        elif not column and name not in {"beta", "zscore", "se"}:
            values = field.get("source_values") or {}
            if values:
                source = "; ".join(str(value) for value in values.values())
        lines.append("")
        lines.extend(screen_field(
            "warning" if missing or conversion.get("invalid") or field.get("source_kind") == "configured_column_absent"
            else "success" if column else "analysis",
            labels.get(name, str(name)), source, indent=8, label_width=22,
        ).splitlines())
        if missing is not None:
            lines.extend(detail("Missing in input", counted(missing, field.get("input_rows"))))
        if conversion.get("status") == "assessed":
            text = counted(conversion["invalid"], conversion["non_missing"]) + " of non-missing values assessed"
            if conversion["invalid"]:
                text += "; could not become numbers; set to missing (not automatically removed here)"
            lines.extend(detail("Invalid numeric", text))
        elif conversion.get("status") == "compound_values_deferred":
            lines.extend(detail("Numeric conversion", "Contains multi-value cells; handled by the INFO-specific parser"))
        elif name == "pos":
            lines.extend(detail("Numeric conversion", "Handled during coordinate validation; not counted in scalar conversion"))
        action = None
        if name == "snp" and (missing or not column):
            action = "Generate chromosome_position_effect_other IDs for affected variants retained at the ID-processing step"
        elif not column and name == "zscore":
            action = "Calculate from harmonised BETA and SE when valid inputs are available"
        elif not column and name in {"beta", "se"}:
            action = field.get("lifecycle_note")
        elif name == "info" and missing:
            mode = report.get("info_missing_action")
            action = {"keep": "Keep variants with missing INFO", "reject": "Remove variants with missing INFO", "fail": "Stop if INFO remains missing"}.get(mode, "Missing INFO action: %s" % mode)
        if action:
            lines.extend(detail("Planned action", action))
    for group in report.get("read_alternative_groups") or []:
        missing = group.get("post_parse_missing_all")
        value = (
            "unavailable because no member maps to a study column"
            if missing is None else
            "%s variants missing all of: %s"
            % (
                counted(missing, group.get("input_rows")),
                ", ".join(labels.get(key, key) for key in group.get("fields") or []),
            )
        )
        lines.extend(screen_field(
            "warning" if missing else "analysis",
            "Input alternatives",
            value,
            indent=8, label_width=24,
        ).splitlines())
    lines.extend(screen_field(
        "loss" if report.get("rows_removed_at_read_mandatory_gate") else "success",
        "Initial removals",
        "%s variants removed for missing required inputs" % _count(report.get("rows_removed_at_read_mandatory_gate")),
        indent=8, label_width=24,
    ).splitlines())
    lines.extend(screen_field(
        "info", "Further checks",
        "Present values and successful numeric conversion do not confirm valid ranges "
        "or statistical consistency; these are checked separately.",
        indent=8, label_width=22,
    ).splitlines())
    return "\n".join(lines)


def _final_field_completeness_block(report):
    """Render post-recovery completeness with its configured count semantics."""
    final = report.get("final_check") or {}
    lines = [
        screen_line("analysis", "Final field completeness", indent=4),
        screen_field(
            "info", "Assessment",
            "after recovery; %s" % final.get("measurement"),
            indent=8, label_width=24,
        ),
        screen_field(
            "info", "Chromosome coverage",
            "%s assessed / %s expected; unassessed %s"
            % (
                _count(len(final.get("completed_chromosomes") or [])),
                _count(len(final.get("expected_chromosomes") or [])),
                ",".join(final.get("unassessed_chromosomes") or []) or "none",
            ),
            indent=8, label_width=24,
        ),
        screen_field(
            "count", "Assessed rows entering",
            (
                "unavailable" if final.get("rows_entering_gate") is None
                else _count(final.get("rows_entering_gate"))
            ),
            indent=8, label_width=24,
        ),
        screen_field(
            "success", "Assessed rows retained",
            (
                "unavailable" if final.get("rows_retained_after_gate") is None
                else _count(final.get("rows_retained_after_gate"))
            ),
            indent=8, label_width=24,
        ),
    ]
    final_fields = sorted(
        (
            field for field in report.get("fields") or []
            if field.get("final_required")
        ),
        key=lambda field: int(field["final_order"]),
    )
    for field in final_fields:
        missing = field.get("final_missing")
        lines.extend(screen_field(
            "warning" if missing is None else "loss" if missing else "success",
            str(field.get("field")),
            "%s attributed missing; action %s; status %s"
            % (
                "unavailable" if missing is None else _count(missing),
                field.get("final_action"),
                field.get("final_status"),
            ),
            indent=8, label_width=24,
        ).splitlines())
    if report.get("report_path"):
        lines.extend(screen_field(
            "info", "QC table", report["report_path"],
            indent=8, label_width=24,
        ).splitlines())
    return "\n".join(lines)


def _ready_variant_type_counts(df, sample_column_dict):
    # type: (pl.DataFrame, Dict[str, Any]) -> Dict[str, int]
    """Count SNPs versus indels/other among input-QC-retained study rows.

    A SNP requires two single-base DNA alleles. Every retained row outside that
    definition is reported as indel/other so the two categories always reconcile
    exactly to ``Ready for harmonisation``.
    """
    ea_col = optional_text(sample_column_dict.get("ea_col"))
    oa_col = optional_text(sample_column_dict.get("oa_col"))
    missing = [
        column
        for column in (ea_col, oa_col)
        if column is None or column not in df.columns
    ]
    if missing:
        raise PipelineError(
            "Cannot classify input variant types because the validated effect- "
            "and other-allele columns are unavailable: %s."
            % ", ".join(str(column) for column in missing)
        )
    bases = ["A", "C", "G", "T"]
    effect = pl.col(ea_col).cast(pl.String, strict=False).str.to_uppercase()
    other = pl.col(oa_col).cast(pl.String, strict=False).str.to_uppercase()
    is_snp = (
        (effect.str.len_chars() == 1)
        & (other.str.len_chars() == 1)
        & effect.is_in(bases)
        & other.is_in(bases)
    ).fill_null(False)
    snps = int(df.select(is_snp.sum().alias("snps")).item() or 0)
    return {
        "snps": snps,
        "indels_or_other_variants": int(df.height) - snps,
    }


def _input_validation_summary_block(
    *,
    input_variants,
    variants_read,
    missing_required,
    invalid_coordinates,
    unsupported_chromosomes,
    non_standard_alleles,
    duplicate_variants,
    ready_variants,
    ready_snps,
    ready_indels_or_other,
):
    """Render input-QC accounting and its ready-row variant-type breakdown."""
    return "\n".join([
        screen_line("count", "Input validation summary", indent=4),
        screen_field(
            "count", "Input variants", _count(input_variants),
            indent=8, label_width=30,
        ),
        screen_field(
            "count", "Variants read", _count(variants_read),
            indent=8, label_width=30,
        ),
        screen_field(
            "loss", "Missing read-stage mandatory values", _count(missing_required),
            indent=8, label_width=30,
        ),
        screen_field(
            "loss", "Invalid coordinates", _count(invalid_coordinates),
            indent=8, label_width=30,
        ),
        screen_field(
            "loss", "Unsupported chromosomes", _count(unsupported_chromosomes),
            indent=8, label_width=30,
        ),
        screen_field(
            "loss", "Non-standard alleles", _count(non_standard_alleles),
            indent=8, label_width=30,
        ),
        screen_field(
            "loss", "Duplicate variants", _count(duplicate_variants),
            indent=8, label_width=30,
        ),
        screen_field(
            "success", "Ready for harmonisation", _count(ready_variants),
            indent=8, label_width=30,
        ),
        screen_field(
            "genetic", "Ready SNPs", _count(ready_snps),
            indent=8, label_width=30,
        ),
        screen_field(
            "genetic", "Ready indels / other variants",
            _count(ready_indels_or_other),
            indent=8, label_width=30,
        ),
    ])


def _chromosome_sort_key(chromosome):
    # type: (Any) -> Tuple[int, str]
    text = str(chromosome).strip().upper()
    if text.isdigit():
        return (int(text), "")
    return ({"X": 100, "Y": 101, "XY": 102, "MT": 103}.get(text, 200), text)


def _safe_glob_name(value, what):
    # type: (Any, str) -> str
    """A config-supplied name that is about to be interpolated into a glob.

    The five ``os.system("rm -f ...")`` calls this replaces interpolated
    ``gwas_outputname`` and ``outdir`` straight from the user's config CSV,
    unquoted, into a shell command.
    """
    text = str(value or "").strip()
    if not text:
        raise PipelineError(
            "%s is empty, so the cleanup pattern would match every file in the "
            "output directory. Refusing to continue." % what
        )
    if os.sep in text or (os.altsep and os.altsep in text) or text in (".", ".."):
        raise PipelineError(
            "%s contains a path separator (%r); it is used as a filename prefix, "
            "not a path." % (what, text)
        )
    return text


def _deep_merge_policy_block(base, override):
    # type: (Dict[str, Any], Dict[str, Any]) -> Dict[str, Any]
    """Merge one policy block over another, group by group."""
    merged = dict(base)
    for group, entries in (override or {}).items():
        if isinstance(entries, dict) and isinstance(merged.get(group), dict):
            inner = dict(merged[group])
            inner.update(entries)
            merged[group] = inner
        else:
            merged[group] = entries
    return merged


def _resolve_policies(default_cfg_obj, sample_column_dict):
    # type: (Dict[str, Any], Dict[str, Any]) -> Any
    """Build the one :class:`Policies` object for this dataset.

    The defaults YAML may carry a ``policies:`` block; a single config row may
    override individual keys through its own ``policies`` cell (a dict, or a
    JSON object as text). Explicit scientific declarations in the sample sheet
    have final precedence. ``load_policies`` reports **every** problem at once,
    which is turned into a :class:`ConfigError`.
    """
    block = {}  # type: Dict[str, Any]
    required_override = {}  # type: Dict[str, Any]

    if isinstance(default_cfg_obj, dict):
        yaml_required = default_cfg_obj.get("required_inputs")
        if isinstance(yaml_required, dict):
            required_override.update(yaml_required)
        candidate = default_cfg_obj.get("policies")
        if isinstance(candidate, dict):
            block = _deep_merge_policy_block(block, candidate)
        elif optional_text(candidate) is not None:
            raise ConfigError(
                [],
                message=(
                    "The 'policies:' entry in the defaults YAML must be a block of "
                    "settings, not %r." % (candidate,)
                ),
            )

    try:
        policies = load_policies(block)
    except PolicyError as exc:
        raise ConfigError(
            [],
            message=(
                "The policy block in the defaults YAML has %d problem(s):\n\n%s"
                % (len(exc.problems), exc)
            ),
        )

    declared_trait = optional_text((sample_column_dict or {}).get("trait_type"))
    trait_policy = {
        "case_control": "binary",
        "binary": "binary",
        "quantitative": "quantitative",
    }.get(str(declared_trait).lower() if declared_trait is not None else "")
    if trait_policy is None and declared_trait is not None and str(declared_trait).lower() != "auto":
        raise ConfigError(
            [],
            message=(
                "The dataset trait_type must be 'auto', 'case_control' (or 'binary') "
                "or 'quantitative'; received %r." % declared_trait
            ),
        )

    override = (sample_column_dict or {}).get("policies")
    if isinstance(override, str):
        text = override.strip()
        if optional_text(text) is not None:
            try:
                override = json.loads(text)
            except ValueError as exc:
                raise ConfigError(
                    [],
                    message=(
                        "The 'policies' column for this dataset is not valid JSON: %s\n"
                        "Value was: %s" % (exc, text)
                    ),
                )
        else:
            override = None
    if isinstance(override, dict) and override:
        try:
            policies = policies.with_overrides(override)
        except PolicyError as exc:
            raise ConfigError(
                [],
                message=(
                    "The per-dataset policy overrides have %d problem(s):\n\n%s"
                    % (len(exc.problems), exc)
                ),
            )
    declarations = {}  # type: Dict[str, Any]
    for sample_key, policy_key in (
        ("declared_effect_type", "effect.type"),
        ("declared_pvalue_type", "pvalue.type"),
        ("delimiter", "input.delimiter"),
    ):
        declared = optional_text((sample_column_dict or {}).get(sample_key))
        if declared is not None and str(declared).strip().lower() != "auto":
            declarations[policy_key] = str(declared).strip().lower()
    if trait_policy is not None:
        declarations["sample_size.trait_type"] = trait_policy
    if declarations:
        try:
            policies = policies.with_overrides(declarations)
        except PolicyError as exc:
            raise ConfigError(
                [],
                message=(
                    "The declared sample-sheet values have %d problem(s):\n\n%s"
                    % (len(exc.problems), exc)
                ),
            )
    if required_override:
        policies = policies.with_required_inputs(required_override)
    return policies


#: Keys the defaults YAML supplies that the validation stages look for in the
#: dataset config. What a dataset must supply is declared in the registry's
#: ``required_inputs:`` block. These defaults describe the comparison-frequency
#: and dbSNP resources only; study EAF and INFO must still come from the dataset.
_DEFAULTS_KEYS_FOR_VALIDATION = (
    "default_comparison_af_file",
    "default_comparison_af_column",
    "default_dbsnp",
)


def _validation_view(sample_column_dict, default_cfg_obj):
    # type: (Dict[str, Any], Dict[str, Any]) -> Dict[str, Any]
    """The dataset config as the validation stages need to see it.

    A copy, not the real dict: the merged keys are only for the three
    validation stages, and must not leak into the dictionary every harmonisation
    step receives.  A per-dataset value always wins over the defaults file.
    """
    view = dict(sample_column_dict or {})
    for key in _DEFAULTS_KEYS_FOR_VALIDATION:
        if key in view and optional_text(view.get(key)) is not None:
            continue
        if isinstance(default_cfg_obj, dict) and key in default_cfg_obj:
            view[key] = default_cfg_obj[key]
    return view


def _log_resolved_policies(logger, policies):
    # type: (Any, Any) -> None
    """Summarise resolution; the manifest stores the complete policy set."""
    if logger is None or policies is None:
        return
    changed = policies.changed_from_default()
    logger.record(
        "OBSERVED", "configuration_resolved",
        parameters=len(policies.keys()), overrides=len(changed),
    )
    for key in sorted(changed):
        was, now = changed[key]
        logger.record("PARAM", key, value=now, canonical_default=was)


def _write_run_manifest(path, manifest):
    # type: (Any, Dict[str, Any]) -> str
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, default=str)
    return str(path)


def _finalise_failed_run_manifest(
    manifest_path, manifest, error, started_at, logger, status="FAILED",
):
    # type: (Any, Optional[Dict[str, Any]], BaseException, float, Any, str) -> bool
    """Persist a terminal manifest state without hiding the original failure.

    ``PARTIAL`` and ``FAILED`` states already assigned by chromosome processing
    remain meaningful. Any later failure after an apparent ``OK`` result is a
    dataset failure, while an explicit interruption always becomes
    ``INTERRUPTED``. Returning ``False`` tells tests and callers that provenance
    could not be persisted; the original exception is deliberately not raised
    or replaced here.
    """
    if manifest is None:
        return False

    requested_status = str(status).upper()
    current_status = str(manifest.get("status") or "").upper()
    if requested_status == "FAILED" and current_status in ("FAILED", "PARTIAL"):
        terminal_status = current_status
    else:
        terminal_status = requested_status

    message = str(error).strip()
    if not message:
        message = (
            "Run interrupted by the user."
            if terminal_status == "INTERRUPTED"
            else "The exception did not provide an error message."
        )
    manifest["status"] = terminal_status
    manifest["failure"] = {
        "exception_type": type(error).__name__,
        "message": message,
    }
    manifest["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    manifest["elapsed_seconds"] = round(time.time() - started_at, 2)

    try:
        _write_run_manifest(manifest_path, manifest)
    except Exception as manifest_error:
        try:
            logger.error(
                "The dataset failed with %s: %s. PostGWAS also could not update "
                "the run manifest %s: %s: %s. The original exception will be "
                "preserved."
                % (
                    type(error).__name__, message, manifest_path,
                    type(manifest_error).__name__, manifest_error,
                )
            )
        except Exception:
            pass
        return False

    try:
        logger.error(
            "Dataset terminated with %s: %s. Run manifest status=%s."
            % (type(error).__name__, message, terminal_status)
        )
    except Exception:
        pass
    return True


def _resolved_harmonisation_outputs(
    outdir, sample_id, input_build, output_layout, vcf_config,
):
    """Return the three required merged VCFs after validating their presence."""
    output_dir = Path(outdir)
    paths = {
        build: configured_output_path(
            output_dir, output_layout["merged_build_vcf"],
            dataset_id=sample_id, build=build,
        )
        for build in vcf_config["target_builds"]
    }
    paths["gwas2vcf"] = configured_output_path(
        output_dir, output_layout["merged_raw_vcf"],
        dataset_id=sample_id, build=input_build,
    )
    missing = [
        "%s=%s" % (name, path)
        for name, path in paths.items()
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing:
        raise PipelineError(
            "Harmonisation finished its analysis but cannot validate all required "
            "merged VCF outputs. Missing or empty: %s" % "; ".join(missing)
        )
    return {name: str(path) for name, path in paths.items()}


def _combined_dataset_status(chromosome_status, merge_status):
    """Combine chromosome and required-merge completeness into one status."""
    chromosome_status = str(chromosome_status).upper()
    merge_status = str(merge_status).upper()
    if chromosome_status not in ("OK", "PARTIAL", "FAILED"):
        raise PipelineError(
            "Unknown chromosome-processing dataset status: %r" % chromosome_status
        )
    if merge_status not in ("OK", "PARTIAL"):
        raise PipelineError("Unknown required-merge status: %r" % merge_status)
    if chromosome_status == "FAILED":
        return "FAILED"
    if chromosome_status == "PARTIAL" or merge_status == "PARTIAL":
        return "PARTIAL"
    return "OK"


def _write_reject_reason_table(path, counts_by_source, delimiter, logger=None):
    # type: (Any, Dict[str, Dict[str, int]], Any, Any) -> Tuple[str, Dict[str, int]]
    """``qc_summary/{sample}_reject_reasons.tsv`` — reason x chromosome counts.

    Every registered reason is written even when it never fired, so a zero is
    positive evidence that the check ran; a missing row would be ambiguous.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    sources = sorted(counts_by_source, key=lambda name: (name != "input", _chromosome_sort_key(name)))
    reasons = sorted(REASONS, key=lambda code: (REASON_STEPS.get(code, "99"), code))

    lines = [delimiter.join(["step", "reject_reason"] + sources + ["total", "description"])]
    column_totals = dict((name, 0) for name in sources)
    reason_totals = {}  # type: Dict[str, int]
    grand_total = 0
    for reason in reasons:
        row_total = 0
        cells = []
        for name in sources:
            value = int(counts_by_source.get(name, {}).get(reason, 0))
            cells.append(str(value))
            row_total += value
            column_totals[name] += value
        grand_total += row_total
        reason_totals[reason] = row_total
        lines.append(
            delimiter.join(
                [REASON_STEPS.get(reason, ""), reason]
                + cells
                + [str(row_total), REASONS[reason]]
            )
        )
    lines.append(
        delimiter.join(
            ["", "TOTAL"]
            + [str(column_totals[name]) for name in sources]
            + [str(grand_total), "every variant removed by a recorded rule"]
        )
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    if logger is not None:
        logger.info(
            "Rejected-variant reasons written to %s (%s variants across %d source%s)."
            % (path, _count(grand_total), len(sources), "" if len(sources) == 1 else "s")
        )
    return str(path), reason_totals


def _required_row_count(value, label):
    # type: (Any, str) -> int
    """Require an internal row count; never coerce missing evidence to zero."""
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise PipelineError(
            "Dataset row reconciliation requires %s to be a non-negative "
            "integer; found %r." % (label, value)
        )
    return value


def _reconcile_dataset_rows(
    *,
    parsed_rows,
    ready_rows,
    partition_rows,
    completed,
    failed,
    chromosome_summaries,
    reject_rows_by_source,
    combined_reject_rows,
    logger=None,
):
    """Hard-assert every parsed row reaches one scientifically honest bucket.

    A row may be exported, rejected by a registered scientific/QC rule, or
    left unprocessed because its chromosome failed. Failed-chromosome
    survivors are deliberately not relabelled as rejects: no variant-level
    rule rejected them.
    """
    parsed_rows = _required_row_count(parsed_rows, "parsed_rows")
    ready_rows = _required_row_count(ready_rows, "ready_rows")
    combined_reject_rows = _required_row_count(
        combined_reject_rows, "combined_reject_rows",
    )
    partitions = {
        str(chromosome): _required_row_count(
            count, "partition_rows[%s]" % chromosome,
        )
        for chromosome, count in dict(partition_rows).items()
    }
    completed = [str(chromosome) for chromosome in completed]
    failed = [str(chromosome) for chromosome in failed]
    if len(completed) != len(set(completed)) or len(failed) != len(set(failed)):
        raise PipelineError(
            "Dataset row reconciliation received duplicate completed or failed "
            "chromosome identifiers."
        )
    if set(completed) & set(failed):
        raise PipelineError(
            "Dataset row reconciliation received chromosomes marked both "
            "completed and failed: %s."
            % ", ".join(sorted(set(completed) & set(failed)))
        )
    status_chromosomes = set(completed) | set(failed)
    if status_chromosomes != set(partitions):
        raise PipelineError(
            "Dataset row reconciliation requires completed and failed "
            "chromosomes to cover exactly the validated partitions. "
            "Partitions=%s; completed=%s; failed=%s."
            % (
                ", ".join(sorted(partitions, key=_chromosome_sort_key)),
                ", ".join(completed),
                ", ".join(failed),
            )
        )

    rejected_by_source = {
        str(source): _required_row_count(
            count, "reject_rows_by_source[%s]" % source,
        )
        for source, count in dict(reject_rows_by_source).items()
    }
    expected_reject_sources = {"input"} | set(partitions)
    if set(rejected_by_source) != expected_reject_sources:
        raise PipelineError(
            "Dataset row reconciliation requires reject counts for input and "
            "every chromosome. Expected=%s; observed=%s."
            % (
                ", ".join(sorted(expected_reject_sources)),
                ", ".join(sorted(rejected_by_source)),
            )
        )
    if sum(rejected_by_source.values()) != combined_reject_rows:
        raise PipelineError(
            "Combined rejected-variant output contains %d row(s), but its "
            "validated source shards contain %d."
            % (combined_reject_rows, sum(rejected_by_source.values()))
        )

    dataset_stage_rejected = rejected_by_source["input"]
    if dataset_stage_rejected + ready_rows != parsed_rows:
        raise PipelineError(
            "End-to-end row reconciliation failed before chromosome work: "
            "%d parsed row(s) != %d input-stage reject(s) + %d ready row(s)."
            % (parsed_rows, dataset_stage_rejected, ready_rows)
        )
    partition_total = sum(partitions.values())
    if partition_total != ready_rows:
        raise PipelineError(
            "End-to-end row reconciliation failed at chromosome partitioning: "
            "%d ready row(s) != %d validated partition row(s)."
            % (ready_rows, partition_total)
        )

    exported_by_chromosome = {}
    for chromosome in completed:
        summary = chromosome_summaries.get(chromosome)
        if not isinstance(summary, dict):
            raise PipelineError(
                "Completed chromosome %s has no row-accounting summary."
                % chromosome
            )
        rows_in = _required_row_count(
            summary.get("rows_in"),
            "chromosome_summaries[%s].rows_in" % chromosome,
        )
        rows_out = _required_row_count(
            summary.get("rows_out"),
            "chromosome_summaries[%s].rows_out" % chromosome,
        )
        summary_rejected = _required_row_count(
            summary.get("rejected"),
            "chromosome_summaries[%s].rejected" % chromosome,
        )
        source_rejected = rejected_by_source[chromosome]
        if rows_in != partitions[chromosome]:
            raise PipelineError(
                "Completed chromosome %s reports %d input row(s), but its "
                "validated partition contains %d."
                % (chromosome, rows_in, partitions[chromosome])
            )
        if summary_rejected != source_rejected:
            raise PipelineError(
                "Completed chromosome %s reports %d rejected row(s), but its "
                "validated rejection shard contains %d."
                % (chromosome, summary_rejected, source_rejected)
            )
        if rows_out + source_rejected != rows_in:
            raise PipelineError(
                "Completed chromosome %s does not reconcile: %d input row(s) "
                "!= %d rejected + %d exported."
                % (chromosome, rows_in, source_rejected, rows_out)
            )
        exported_by_chromosome[chromosome] = rows_out

    unprocessed_by_chromosome = {}
    for chromosome in failed:
        source_rejected = rejected_by_source[chromosome]
        partition_count = partitions[chromosome]
        if source_rejected > partition_count:
            raise PipelineError(
                "Failed chromosome %s has %d recorded reject(s) for only %d "
                "partition row(s)."
                % (chromosome, source_rejected, partition_count)
            )
        summary = chromosome_summaries.get(chromosome)
        if isinstance(summary, dict) and summary.get("rejected") is not None:
            summary_rejected = _required_row_count(
                summary.get("rejected"),
                "chromosome_summaries[%s].rejected" % chromosome,
            )
            if summary_rejected != source_rejected:
                raise PipelineError(
                    "Failed chromosome %s reports %d rejected row(s), but its "
                    "validated current-attempt rejection shard contains %d."
                    % (chromosome, summary_rejected, source_rejected)
                )
        unprocessed_by_chromosome[chromosome] = (
            partition_count - source_rejected
        )

    rows_exported = sum(exported_by_chromosome.values())
    chromosome_stage_rejected = sum(
        rejected_by_source[chromosome] for chromosome in partitions
    )
    rejected_total = dataset_stage_rejected + chromosome_stage_rejected
    unprocessed_total = sum(unprocessed_by_chromosome.values())
    accounted_rows = rows_exported + rejected_total + unprocessed_total
    if accounted_rows != parsed_rows:
        raise PipelineError(
            "End-to-end dataset row reconciliation failed: %d parsed row(s) "
            "!= %d rejected + %d exported + %d unprocessed because a "
            "chromosome failed."
            % (
                parsed_rows,
                rejected_total,
                rows_exported,
                unprocessed_total,
            )
        )
    if not failed and unprocessed_total:
        raise PipelineError(
            "Dataset row reconciliation found %d unprocessed row(s) although "
            "every chromosome completed." % unprocessed_total
        )

    result = {
        "scope": "end_to_end_parsed_rows",
        "chromosomes": sorted(partitions, key=_chromosome_sort_key),
        "completed_chromosomes": completed,
        "failed_chromosomes": failed,
        "rows_read": parsed_rows,
        "dataset_stage_rejected": dataset_stage_rejected,
        "rows_ready_for_chromosomes": ready_rows,
        "partition_rows_by_chromosome": dict(partitions),
        "chromosome_stage_rejected": chromosome_stage_rejected,
        "rejected": rejected_total,
        "rows_exported_by_chromosome": exported_by_chromosome,
        "rows_exported": rows_exported,
        "unprocessed_rows_by_chromosome": unprocessed_by_chromosome,
        "unprocessed_failed_chromosome_rows": unprocessed_total,
        "accounted_rows": accounted_rows,
        "balanced": True,
        "complete": not failed,
    }
    if logger is not None:
        logger.info(
            "End-to-end row accounting balanced: %s parsed = %s rejected + "
            "%s exported + %s unprocessed because a chromosome failed."
            % (
                _count(parsed_rows),
                _count(rejected_total),
                _count(rows_exported),
                _count(unprocessed_total),
            )
        )
    return result


def _log_tail(path, max_lines=400):
    # type: (Any, int) -> str
    """The end of a chromosome log, for the parent to print when a worker died.

    The log file is the source of truth; a worker that was killed outright
    returns nothing at all, but its log is complete up to the moment it died
    because every record is flushed as it is written.
    """
    try:
        with open(str(path), "r", encoding="utf-8", errors="replace") as handle:
            lines = handle.read().splitlines()
    except (IOError, OSError) as exc:
        return "(the log file %s could not be read: %s)" % (path, exc)
    if not lines:
        return "(the log file %s is empty: the worker died before its first record)" % (path,)
    if len(lines) > max_lines:
        omitted = len(lines) - max_lines
        lines = ["(... %d earlier lines are in %s ...)" % (omitted, path)] + lines[-max_lines:]
    return "\n".join(lines)


def _record_chromosome_exception(
    logger, chromosome, attempt, exc, *, non_retryable=False,
):
    # type: (Any, str, int, BaseException, bool) -> Dict[str, str]
    """Preserve an unhandled worker traceback without duplicating step logs."""
    error = "%s: %s" % (type(exc).__name__, exc)
    trace = traceback.format_exc().rstrip()
    suffix = (
        " This memory failure is non-retryable; PostGWAS will stop the dataset "
        "instead of repeating the chromosome with unchanged memory settings."
        if non_retryable else ""
    )
    logger.error(
        "Chromosome %s did not complete on attempt %d. %s%s"
        % (chromosome, attempt, error, suffix)
    )

    # logger.step() already writes the full traceback for failures raised from
    # a numbered step. Orchestration between steps has no such wrapper, so its
    # traceback must be written here.
    steps = logger.summary().get("steps") or []
    if not steps or steps[-1].get("status") != "failed":
        logger.log(
            "", "For developers - traceback:",
            indent=1, level="ERROR", screen=False,
        )
        for line in trace.split("\n"):
            logger.log(
                "", line, indent=2, level="ERROR", wrap=False, screen=False,
            )
    return {"error": error, "traceback": trace}


def _derive_parallelism(policies, requested_workers, n_chromosomes, logger=None):
    # type: (Any, int, int, Any) -> Tuple[int, int]
    """How many workers, and how many threads each may use.

    ``execution.threads_per_chromosome`` is the per-process ceiling shared by
    Polars and bcftools. When ``execution.total_cpu_budget`` is set, that
    per-worker allocation is preserved and worker count becomes
    ``budget // threads``. Only a budget smaller than one worker's requested
    allocation reduces the per-worker value. The independently resolved memory
    budget can reduce worker count further.
    """
    requested_threads = max(
        1, int(policies.get("execution.threads_per_chromosome"))
    )
    workers_requested = max(
        1,
        min(
            int(requested_workers or 1),
            max(1, int(n_chromosomes or 1)),
        ),
    )
    workers = workers_requested
    threads = requested_threads
    cpu_budget_value = policies.get("execution.total_cpu_budget")
    memory_budget_value = policies.get("execution.memory_budget_gb")
    memory_per_worker = float(
        policies.get("execution.memory_gb_per_chromosome")
    )

    cpu_budget = None
    if cpu_budget_value not in (None, ""):
        cpu_budget = max(1, int(cpu_budget_value))
        threads = min(requested_threads, cpu_budget)
        workers = min(workers, max(1, cpu_budget // threads))

    memory_budget = None
    memory_worker_limit = None
    if memory_budget_value not in (None, ""):
        memory_budget = float(memory_budget_value)
        memory_worker_limit = safe_thread_count(
            workers,
            memory_per_worker,
            available_ram_gb=memory_budget,
            reporter=None,
        )
        workers = min(workers, memory_worker_limit)

    if logger is not None:
        logger.decide(
            "How much chromosome parallelism to use",
            {
                "workers requested": workers_requested,
                "threads per chromosome requested": requested_threads,
                "total CPU budget": cpu_budget,
                "memory budget GB": memory_budget,
                "memory reservation GB per chromosome": memory_per_worker,
                "memory-limited workers": memory_worker_limit,
                "Polars worker variable": _POLARS_THREAD_ENVIRONMENT_VARIABLE,
            },
            "%d worker%s x %d thread%s = at most %d CPU thread%s; "
            "each worker reserves %.1f GB"
            % (
                workers,
                "" if workers == 1 else "s",
                threads,
                "" if threads == 1 else "s",
                workers * threads,
                "" if workers * threads == 1 else "s",
                memory_per_worker,
            ),
        )
    return workers, threads


# ===============================================================
# 1️⃣ Resource Map Builder
# ===============================================================

def build_resource_map(
    chromosome: str,
    grch_version: str,
    target_build: str,
    resource_folder: str,
    user_eaf_file: str,
    default_eaf_reference_source: str,
    default_comparison_af_file: str,
    resource_layout: Dict[str, str],
    user_info_file: Optional[str],
    dbsnp: str,
    user_eaf_column: Optional[str],
    default_eaf_reference_column: Optional[str],
    default_comparison_af_column: Optional[str],
    user_info_column: Optional[str],
    require_default_eaf: bool = False,
    validate_files: bool = True,
) -> Dict[str, Any]:
    """Construct one chromosome's resource map.

    Direct callers retain the historical existence checks. Dataset preflight
    disables that narrow check because its aggregated validator performs the
    same file checks plus indexes, schemas, and cross-file contig validation.
    """
    fields = {
        "build": grch_version,
        "target_build": target_build,
        "source_build": grch_version,
        "chromosome": chromosome,
    }

    def configured_path(name: str, **overrides: str) -> str:
        template = resource_layout.get(name)
        if not template:
            raise PipelineError(
                "The harmonisation resource_layout.%s setting is missing." % name
            )
        values = dict(fields)
        values.update(overrides)
        try:
            relative = Path(template.format(**values))
        except (KeyError, ValueError) as exc:
            raise PipelineError(
                "Invalid resource_layout.%s template %r: %s"
                % (name, template, exc)
            )
        if relative.is_absolute() or ".." in relative.parts:
            raise PipelineError(
                "resource_layout.%s must stay below the resource directory: %s"
                % (name, relative)
            )
        return str(Path(resource_folder) / relative)

    user_eaf_path = resolve_resource_file(
        input_file=user_eaf_file,
        grch_version=grch_version,
        chromosome=chromosome,
        must_exist=validate_files,
    )
    default_eaf_path = configured_path(
        "default_eaf", source=default_eaf_reference_source,
    )
    default_comparison_af_path = configured_path(
        "comparison_af", source=default_comparison_af_file,
    )
    user_info_path = resolve_resource_file(
        input_file=user_info_file,
        grch_version=grch_version,
        chromosome=chromosome,
        must_exist=validate_files,
    )
    # --- Core genome references ---
    genome_fasta_path = configured_path("fasta")
    dbsnp_path = configured_path("dbsnp", source=dbsnp)
    annot_path = configured_path("annotation")
    target_fasta = configured_path("fasta", build=target_build)
    chain_file = configured_path("chain")
    # --- Validate core default files exist & are non-empty ---
    required_resource_paths = [
        default_comparison_af_path,
        dbsnp_path,
        genome_fasta_path,
        target_fasta,
        annot_path,
        chain_file,
    ]
    if require_default_eaf:
        required_resource_paths.insert(0, default_eaf_path)
    if validate_files:
        file_validator = validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
        )
        errors = []
        for path in required_resource_paths:
            try:
                file_validator(path)
            except argparse.ArgumentTypeError as e:
                errors.append(f"❌ File invalid: {path}\n   → {e}")
        if errors:
            raise PipelineError(
                "\n❌ One or more default files are invalid:\n" + "\n".join(errors)
            )
    return {
        # EAF
        "user_eaf_file": user_eaf_path,
        "default_eaf_file": default_eaf_path,
        "default_comparison_af_file": default_comparison_af_path,
        "user_eaf_column": user_eaf_column,
        "default_eaf_reference_column": default_eaf_reference_column,
        "default_comparison_af_column": default_comparison_af_column,
        "user_info_file": user_info_path,
        "user_info_column": user_info_column,
        # Genome, dbSNP, chain
        "genome_fasta_file": genome_fasta_path,
        "dbsnp_file": dbsnp_path,
        "annot_path": annot_path,
        "target_fasta": target_fasta,
        "chain_file": chain_file,
    }


def preflight_harmonisation_resources(
    chromosomes,
    *,
    grch_version,
    target_build,
    resource_folder,
    user_eaf_file,
    default_eaf_reference_source,
    default_comparison_af_file,
    resource_layout,
    user_info_file,
    dbsnp,
    user_eaf_column,
    default_eaf_reference_column,
    default_comparison_af_column,
    user_info_column,
    default_eaf_colmap,
    external_eaf_colmap,
    external_info_colmap,
    executables,
    vcf_config,
    policies,
    require_default_eaf,
    logger=None,
):
    """Resolve and validate every resource before chromosome workers start."""
    resource_maps = {}
    issues = []
    for chromosome in chromosomes:
        try:
            resource_maps[str(chromosome)] = build_resource_map(
                chromosome=str(chromosome),
                grch_version=grch_version,
                target_build=target_build,
                resource_folder=resource_folder,
                user_eaf_file=user_eaf_file,
                default_eaf_reference_source=default_eaf_reference_source,
                default_comparison_af_file=default_comparison_af_file,
                resource_layout=resource_layout,
                user_info_file=user_info_file,
                user_eaf_column=user_eaf_column,
                default_eaf_reference_column=default_eaf_reference_column,
                default_comparison_af_column=default_comparison_af_column,
                user_info_column=user_info_column,
                dbsnp=dbsnp,
                require_default_eaf=require_default_eaf,
                validate_files=False,
            )
        except (OSError, ValueError, PipelineError) as exc:
            issues.append({
                "resource": "resource path resolution",
                "path": str(resource_folder),
                "problem": "%s: %s" % (type(exc).__name__, exc),
                "chromosomes": [str(chromosome)],
            })

    if resource_maps:
        try:
            summary = validate_harmonisation_resource_maps(
                resource_maps,
                bcftools=executables["bcftools"],
                vcf_config=vcf_config,
                default_eaf_colmap=default_eaf_colmap,
                external_eaf_colmap=external_eaf_colmap,
                external_info_colmap=external_info_colmap,
                policies=policies,
                require_default_eaf=require_default_eaf,
                logger=logger,
            )
        except ResourcePreflightError as exc:
            issues.extend(exc.issues)
            summary = None
    else:
        summary = None

    if issues:
        raise ResourcePreflightError(issues)
    return resource_maps, summary


# ===============================================================
# 2️⃣ Per-Chromosome Harmonization Pipeline  (16 steps)
# ===============================================================

def _finalize_chromosome_rejects(rejects, status, qc_dict):
    # type: (RejectCollector, str, Dict[str, Any]) -> str
    """Flush required provenance without replacing an earlier analysis error."""
    try:
        rejects.flush()
    except Exception as exc:
        message = "%s: %s" % (type(exc).__name__, exc)
        qc_dict["reject_output_error"] = message
        if status == "ok":
            qc_dict["error"] = message
        status = "failed"
    finally:
        qc_dict.setdefault("rejects", rejects.summary())
    return status


def process_one_chromosome(
    chromosome: str,
    chr_file: str,
    sample_column_dict: Dict[str, Any],
    resource_folder: str,
    grch_version: str,
    user_eaf_file: str,
    default_eaf_reference_source: str,
    default_comparison_af_file: str,
    resource_layout: Dict[str, str],
    output_layout: Dict[str, str],
    gwas2vcf_input: Dict[str, Any],
    executables: Dict[str, str],
    vcf_config: Dict[str, Any],
    gwas2vcf_main_script_path: str,
    user_info_file: Optional[str],
    user_eaf_column: Optional[str],
    default_eaf_reference_column: Optional[str],
    default_comparison_af_column: Optional[str],
    user_info_column: Optional[str],
    dbsnp: str,
    output_dir: str,
    threads: int,
    external_eaf_colmap,
    default_eaf_colmap,
    external_info_colmap,
    policies,
    study_decisions,
    missing_sample_size_plan,
    attempt: int,
    prevalidated_resource_map=None,
) -> Tuple[str, Dict[str, Any], str, str]:
    """Run the sixteen harmonisation steps for one chromosome, in a worker.

    ``chr_file`` must come from this pipeline's dataset-level validated split.
    Raw, independently supplied chromosome files are outside this function's
    contract and are never routed here by the CLI or pipeline.

    Returns ``(chromosome, qc_dict, screen_text, status)``.

    ``screen_text`` is captured in a ``finally`` and returned on both paths, so
    a chromosome that fails at step 9 still hands the parent everything it
    logged.  ``status`` is ``'ok'`` or ``'failed'``. Ordinary failures are
    returned to the retry controller; ``MemoryError`` is deliberately raised
    so that repeating the chromosome under unchanged memory pressure is
    impossible. The worker never writes to stdout: the parent prints blocks.
    """
    sample_id = sample_column_dict["gwas_outputname"]
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    log_dir = configured_output_path(
        output_dir_path, output_layout["logs_directory"], error_type=PipelineError,
    )
    log_dir.mkdir(parents=True, exist_ok=True)

    pol = policies if policies is not None else default_policies()

    # The log is opened in append mode so every attempt for this chromosome
    # ends up in one file.  Only an explicit policy turns that off.
    if attempt > 1 and not bool(pol.get("logging.append_across_retries")):
        try:
            configured_output_path(
                output_dir_path,
                output_layout["chromosome_log"],
                error_type=PipelineError,
                dataset_id=sample_id,
                chromosome=chromosome,
            ).unlink()
        except (OSError, FileNotFoundError):
            pass

    logger = PipelineLogger(
        sample_id=sample_id,
        scope=str(chromosome),
        log_dir=str(log_dir),
        policies=pol,
        level=None,
        screen_level=None,
        log_path=str(configured_output_path(
            output_dir_path,
            output_layout["chromosome_log"],
            error_type=PipelineError,
            dataset_id=sample_id,
            chromosome=chromosome,
        )),
    )

    qc_dict = {}  # type: Dict[str, Any]
    rejects = None
    reconciliation = None  # type: Optional[Dict[str, Any]]
    rows_read = 0
    status = "failed"
    screen_text = ""
    gwas2vcf_exit_code = None

    try:
        logger.blank()
        logger.info(
            "Chromosome %s of %s, attempt %d. Working directory %s."
            % (chromosome, sample_id, attempt, output_dir_path)
        )

        # --------------------------------------------------------------
        # 01  load_chromosome
        # --------------------------------------------------------------
        with logger.step(
            1, CHR_STEP_TOTAL, "Load the chromosome file", "polars.read_parquet"
        ) as ctx:
            df = pl.read_parquet(chr_file)
            rows_read = df.height
            ctx.info(
                "Read %s variants and %d columns from %s."
                % (_count(df.height), df.width, chr_file)
            )
            ctx.set_rows(df.height, removed=0)

        if bool(pol.get("rejects.enabled")):
            reject_dir = configured_output_path(
                output_dir_path,
                output_layout["rejected_directory"],
                error_type=PipelineError,
            )
            reject_dir.mkdir(parents=True, exist_ok=True)
            rejects = RejectCollector(
                source_snapshot=str(configured_output_path(
                    output_dir_path,
                    output_layout["chromosome_source_snapshot"],
                    error_type=PipelineError,
                    dataset_id=sample_id,
                    chromosome=chromosome,
                )),
                logger=logger,
                out_path=str(configured_output_path(
                    output_dir_path,
                    output_layout["chromosome_reject"],
                    error_type=PipelineError,
                    dataset_id=sample_id,
                    chromosome=chromosome,
                )),
                delimiter=gwas2vcf_input["delimiter"],
                compress=bool(pol.get("rejects.compress")),
            )
        else:
            logger.warn(
                "rejects.enabled is false, so removed variants are counted in the log "
                "but not written to a file, and the reconciliation check cannot run."
            )

        # --------------------------------------------------------------
        # 02  build_resource_map
        # --------------------------------------------------------------
        with logger.step(
            2, CHR_STEP_TOTAL, "Reference resources", "main.build_resource_map",
            rows_in=df.height,
        ) as ctx:
            require_default_eaf = True
            if prevalidated_resource_map is None:
                res = build_resource_map(
                    chromosome=chromosome,
                    grch_version=grch_version,
                    target_build=vcf_config["target_builds"][grch_version],
                    resource_folder=resource_folder,
                    user_eaf_file=user_eaf_file,
                    default_eaf_reference_source=default_eaf_reference_source,
                    default_comparison_af_file=default_comparison_af_file,
                    resource_layout=resource_layout,
                    user_info_file=user_info_file,
                    user_eaf_column=user_eaf_column,
                    default_eaf_reference_column=default_eaf_reference_column,
                    default_comparison_af_column=default_comparison_af_column,
                    user_info_column=user_info_column,
                    dbsnp=dbsnp,
                    require_default_eaf=require_default_eaf,
                )
                resource_source = "resolved in chromosome worker"
            else:
                res = dict(prevalidated_resource_map)
                recheck_preflighted_resource_map(
                    chromosome,
                    res,
                    require_default_eaf=require_default_eaf,
                )
                resource_source = "dataset preflight; availability rechecked"
            ctx.info(
                "Resources from %s. Default EAF %s, comparison AF %s, "
                "dbSNP %s, FASTA %s."
                % (resource_source, res["default_eaf_file"],
                   res["default_comparison_af_file"], res["dbsnp_file"],
                   res["genome_fasta_file"])
            )
            ctx.extra["resource_source"] = resource_source
            ctx.info(
                "Target FASTA %s, annotation %s, chain %s."
                % (res["target_fasta"], res["annot_path"], res["chain_file"])
            )
            ctx.extra["resources"] = dict(res)
            ctx.set_rows(df.height, removed=0)

        # --------------------------------------------------------------
        # 03  effect_type
        #     Normalize OR and any raw OR-scale SE before allele swapping.
        #     On log-odds units a swap negates beta and leaves SE unchanged.
        # --------------------------------------------------------------
        effect_col = optional_text(sample_column_dict.get("beta_or_col"))
        if effect_col is not None and effect_col in df.columns:
            df, or_qc, sample_column_dict = harmonise_effect_estimates(
                chromosome=chromosome,
                df=df,
                sample_column_dict=sample_column_dict,
                logger=logger,
                policies=pol,
                rejects=rejects,
                decision=study_decisions,
                step_number=3,
                step_total=CHR_STEP_TOTAL,
            )
        else:
            with logger.step(
                3, CHR_STEP_TOTAL, "Effect type (beta or odds ratio)",
                "effect_type.harmonise_effect_estimates",
                rows_in=df.height,
            ) as ctx:
                ctx.skip(
                    "No input effect column is present. The effect-from-Z step "
                    "will create a beta directly, so no OR conversion is needed."
                )
                ctx.set_rows(df.height, removed=0)
            or_qc = {
                "initial_variants": df.height,
                "status": "skipped_no_input_effect_column",
                "effect_type": "beta",
                "effect_decision_source": "derived_from_z",
                "final_total": df.height,
            }

        orientation_decisions = dict(study_decisions or {})
        if sample_column_dict.get("harmonised_effect_type") == "beta":
            orientation_decisions["input_effect_type"] = orientation_decisions.get(
                "effect_type"
            )
            orientation_decisions["effect_type"] = "beta"

        # --------------------------------------------------------------
        # 04  eaf_harmonisation and reference allele orientation
        # --------------------------------------------------------------
        df, eaf_qc, sample_column_dict = harmonise_allele_frequency(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
            output_layout=output_layout,
            output_delimiter=gwas2vcf_input["delimiter"],
            eaffile=res["user_eaf_file"],
            external_eaf_colmap=external_eaf_colmap,
            default_eaf_file=res["default_eaf_file"],
            default_eaf_column=res["default_eaf_reference_column"],
            default_eaf_colmap=default_eaf_colmap,
            policies=pol,
            logger=logger,
            rejects=rejects,
            step_number=4,
            step_total=CHR_STEP_TOTAL,
            study_decision=orientation_decisions,
        )

        # --------------------------------------------------------------
        # 05  sample_size
        # --------------------------------------------------------------
        df, size_qc, sample_column_dict = harmonise_sample_sizes(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
            logger=logger,
            policies=pol,
            rejects=rejects,
            missing_plan=missing_sample_size_plan,
            step_number=5,
            step_total=CHR_STEP_TOTAL,
        )

        # --------------------------------------------------------------
        # 06  beta_se_from_z
        #     The effect mapping now always points to log-odds beta.
        # --------------------------------------------------------------
        with logger.step(
            6, CHR_STEP_TOTAL, "Effect size and standard error from Z",
            "effect_from_z.derive_effect_and_standard_error_from_z",
            rows_in=df.height,
        ) as ctx:
            df, beta_qc, sample_column_dict = derive_effect_and_standard_error_from_z(
                chromosome=chromosome,
                df=df,
                sample_column_dict=sample_column_dict,
                logger=logger,
                policies=pol,
                rejects=rejects,
            )
            ctx.set_rows(df.height)

        # --------------------------------------------------------------
        # 07  pvalue_detection
        # --------------------------------------------------------------
        df, pval_qc, sample_column_dict = harmonise_p_values(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
            logger=logger,
            policies=pol,
            rejects=rejects,
            decision=study_decisions,
            step_number=7,
            step_total=CHR_STEP_TOTAL,
        )

        # --------------------------------------------------------------
        # 08  se_from_beta_pval
        # --------------------------------------------------------------
        with logger.step(
            8, CHR_STEP_TOTAL, "Standard error from effect size and p-value",
            "standard_error.derive_standard_error_from_effect_and_p_value",
            rows_in=df.height,
        ) as ctx:
            df, se_qc, sample_column_dict = derive_standard_error_from_effect_and_p_value(
                chromosome=chromosome,
                df=df,
                sample_column_dict=sample_column_dict,
                logger=logger,
                policies=pol,
                rejects=rejects,
            )
            ctx.set_rows(df.height)

        # --------------------------------------------------------------
        # 09  z_from_beta_se
        # --------------------------------------------------------------
        with logger.step(
            9, CHR_STEP_TOTAL, "Z score from effect size and standard error",
            "z_score.derive_z_score_from_effect_and_standard_error",
            rows_in=df.height,
        ) as ctx:
            df, z_qc, sample_column_dict = derive_z_score_from_effect_and_standard_error(
                chromosome=chromosome,
                df=df,
                sample_column_dict=sample_column_dict,
                logger=logger,
                policies=pol,
            )
            ctx.set_rows(df.height)

        # --------------------------------------------------------------
        # 10  effect_statistics_validation
        # --------------------------------------------------------------
        df, effect_validation_qc, sample_column_dict = validate_effect_statistics(
            chromosome,
            df,
            sample_column_dict,
            policies=pol,
            logger=logger,
            rejects=rejects,
            step_number=10,
            step_total=CHR_STEP_TOTAL,
        )

        # --------------------------------------------------------------
        # 11  info_harmonisation
        # --------------------------------------------------------------
        df, info_qc, sample_column_dict = harmonise_imputation_quality(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
            info_file=res["user_info_file"],
            info_column=res["user_info_column"],
            external_info_colmap=external_info_colmap,
            logger=logger,
            policies=pol,
            rejects=rejects,
            step_number=11,
            step_total=CHR_STEP_TOTAL,
        )

        # --------------------------------------------------------------
        # 12  snp_column
        # --------------------------------------------------------------
        df, sample_column_dict = harmonise_variant_identifiers(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
            logger=logger,
            policies=pol,
            rejects=rejects,
            step_number=12,
            step_total=CHR_STEP_TOTAL,
        )

        # --------------------------------------------------------------
        # 13  final_completeness_check
        # --------------------------------------------------------------
        df, final_check_qc, sample_column_dict = final_completeness_check(
            chromosome,
            df,
            sample_column_dict,
            policies=pol,
            logger=logger,
            rejects=rejects,
            step_number=13,
            step_total=CHR_STEP_TOTAL,
        )

        # --------------------------------------------------------------
        # Reconciliation — a hard assertion (plan Part 4).
        # If rows read minus rows rejected does not equal the rows about to be
        # exported, a variant leaked or was counted twice.  Stop.
        # --------------------------------------------------------------
        if rejects is not None:
            reconciliation = reconcile(
                rows_read, df.height, rejects,
                logger=logger, label="chromosome %s" % chromosome,
            )

        # --------------------------------------------------------------
        # 14  export_sumstat
        # --------------------------------------------------------------
        with logger.step(
            14, CHR_STEP_TOTAL, "Export the harmonised summary statistics",
            "gwas2vcf_export.export_gwas2vcf_input",
            rows_in=df.height,
        ) as ctx:
            export_gwas2vcf_input(
                df=df,
                sample_column_dict=sample_column_dict,
                output_dir=str(output_dir_path),
                gwas_outputname=sample_column_dict["gwas_outputname"],
                chromosome=chromosome,
                genome_build=grch_version,
                layout=output_layout,
                input_config=gwas2vcf_input,
                logger=logger,
            )
            ctx.set_rows(df.height, removed=0)

        # --------------------------------------------------------------
        # 15  gwas_to_vcf
        # --------------------------------------------------------------
        with logger.step(
            15, CHR_STEP_TOTAL, "Convert to VCF", "gwas2vcf_runner.run_gwas2vcf",
            rows_in=df.height,
        ) as ctx:
            _gwas2vcf_command, gwas2vcf_exit_code = run_gwas2vcf(
                gwas_outputname=sample_column_dict["gwas_outputname"],
                chromosome=chromosome,
                grch_version=grch_version,
                output_folder=str(output_dir_path),
                fasta=res["genome_fasta_file"],
                main_script_path=gwas2vcf_main_script_path,
                output_layout=output_layout,
                python_bin=executables["python"],
                bcftools_bin=executables["bcftools"],
                expected_variants=df.height,
                logger=logger,
            )
            ctx.extra["gwas2vcf_exit_code"] = gwas2vcf_exit_code
            ctx.set_rows(df.height, removed=0)

        # --------------------------------------------------------------
        # 16  bcftools_annotate_liftover
        # --------------------------------------------------------------
        with logger.step(
            16, CHR_STEP_TOTAL, "Annotate and lift over",
            "vcf_processing.annotate_and_liftover_vcf",
            rows_in=df.height,
            policy_keys=[
                "vcf.liftover_swap",
                "vcf.liftover_warn_fraction",
                "vcf.liftover_critical_fraction",
                "vcf.liftover_fail_fraction",
                "vcf.variant_drop_warn_fraction",
                "vcf.sort_memory_mb_per_thread",
            ],
        ) as ctx:
            target_vcf = annotate_and_liftover_vcf(
                output_dir=str(output_dir_path),
                gwas_outputname=sample_column_dict["gwas_outputname"],
                chromosome=chromosome,
                external_eaf_file=res["default_comparison_af_file"],
                default_dbsnp_file=res["dbsnp_file"],
                genome_fasta_file=res["genome_fasta_file"],
                target_genome_fasta_file=res["target_fasta"],
                gff_file=res["annot_path"],
                chain_file=res["chain_file"],
                grch_version=grch_version,
                threads=threads,
                policies=pol,
                output_layout=output_layout,
                executables=executables,
                vcf_config=vcf_config,
                logger=logger,
                qc_info=ctx.extra.setdefault("liftover_accounting", {}),
            )
            ctx.extra["target_vcf"] = target_vcf
            ctx.set_rows(df.height, removed=0)

        qc_dict = {
            "eaf_qc": eaf_qc,
            "sample_size_qc": size_qc,
            "effect_from_z_qc": beta_qc,
            "beta_or_oddsratio_qc": or_qc,
            "pval_detection_and_conversion_qc": pval_qc,
            "se_from_beta_pvalue_qc": se_qc,
            "z_from_beta_se_qc": z_qc,
            "effect_statistics_validation_qc": effect_validation_qc,
            "info_qc": info_qc,
            "final_completeness_qc": final_check_qc,
            "gwas2vcf_exit_code": gwas2vcf_exit_code,
        }
        status = "ok"
        logger.info("Chromosome %s completed all %d steps." % (chromosome, CHR_STEP_TOTAL))

    except KeyboardInterrupt:
        raise
    except MemoryError as exc:
        qc_dict = _record_chromosome_exception(
            logger, chromosome, attempt, exc, non_retryable=True,
        )
        raise
    except Exception as exc:
        status = "failed"
        qc_dict = _record_chromosome_exception(
            logger, chromosome, attempt, exc,
        )
    finally:
        # Flushed from a finally: a chromosome that failed at step 9 still has
        # steps 1-8 of its rejected variants on disk.
        if rejects is not None:
            status = _finalize_chromosome_rejects(rejects, status, qc_dict)

        summary = logger.summary()
        summary["stage_qc"] = dict(qc_dict)
        if reconciliation:
            summary.update(reconciliation)
        else:
            summary["rows_in"] = rows_read
            if rejects is not None:
                summary["rejected"] = rejects.total()
                summary["reject_counts"] = rejects.counts()
        summary["attempt"] = attempt
        summary["status"] = status
        qc_dict["__block_summary__"] = summary

        # Captured on both paths, so the parent can print the block even when
        # the chromosome failed.
        screen_text = logger.screen_text()
        logger.close()

    return chromosome, qc_dict, screen_text, status


# ===============================================================
#  MULTIPROCESSING VERSION — Python 3.8 / Docker SAFE
# ===============================================================

def _run_one_round(
    pending,
    chr_file_by_chrom,
    resource_maps,
    round_no,
    attempts,
    worker_kwargs,
    max_workers,
    sample_id,
    output_dir,
    output_layout,
    screen_order,
    logger,
    progress=None,
    schedule=None,
):
    """Run one pass over `pending`, printing exactly one block per chromosome.

    The ``as_completed`` loop is single-threaded and is the only place that
    writes to stdout, which is what makes interleaving impossible.
    """
    results = {}   # type: Dict[str, Dict[str, Any]]
    blocks = []    # type: List[Tuple[Any, Dict[str, Any]]]
    mp_ctx = get_context("spawn")
    workers = max(1, min(int(max_workers), len(pending)))

    def emit_block(chromosome, payload):
        args = dict(
            chromosome=chromosome,
            sample_id=sample_id,
            attempt=payload["attempt"],
            status=payload["status_label"],
            elapsed=payload.get("elapsed"),
            screen_text=payload.get("screen_text", ""),
            summary=payload.get("summary"),
        )
        if screen_order == "sorted":
            blocks.append((chromosome, args))
        else:
            print_chromosome_summary(**args)

    started = {}   # type: Dict[str, float]
    try:
        with _bounded_chromosome_executor(
            workers,
            worker_kwargs["threads"],
            mp_ctx,
        ) as executor:
            future_to_chr = {}
            queue = ResourceQueue(
                {chrom: schedule["estimated_memory_gb_by_chromosome"][chrom] for chrom in pending}
                if schedule else {chrom: 1.0 for chrom in pending},
                schedule["worker_memory_budget_gb"] if schedule else float(workers),
                workers,
                mix_sizes=schedule is not None,
            )

            def submit_chromosome(chrom):
                chromosome_args = dict(worker_kwargs)
                chromosome_args.update(
                    chromosome=chrom,
                    chr_file=chr_file_by_chrom[chrom],
                    sample_column_dict=dict(worker_kwargs["sample_column_dict"]),
                    attempt=attempts.get(chrom, 0) + 1,
                    prevalidated_resource_map=dict(resource_maps[chrom]),
                )
                try:
                    future = executor.submit(
                        process_one_chromosome,
                        **chromosome_args,
                    )
                except BrokenProcessPool as exc:
                    # Route a broken-pool submission through the same failure
                    # accounting as a process that dies after submission.
                    future = Future()
                    future.set_exception(exc)
                future_to_chr[future] = chrom
                started[chrom] = time.time()

            while queue.pending or future_to_chr:
                while (chrom := queue.admit()) is not None:
                    submit_chromosome(chrom)
                    if schedule and logger is not None:
                        logger.info(
                            "Chromosome scheduling: started %s; active %s; "
                            "%d/%d workers; %d CPU threads; estimated worker "
                            "memory %.2f/%.2f GB."
                            % (chrom, ", ".join(queue.active), len(queue.active), workers,
                               len(queue.active) * worker_kwargs["threads"],
                               queue.reserved, queue.budget)
                        )
                future = next(iter(as_completed(future_to_chr)))
                chrom = future_to_chr.pop(future)
                queue.complete(chrom)
                attempt = attempts.get(chrom, 0) + 1
                elapsed = time.time() - started.get(chrom, time.time())
                try:
                    chrom_res, qc, screen_text, status = future.result()
                except MemoryError as exc:
                    for other_future in future_to_chr:
                        if other_future is not future:
                            other_future.cancel()
                    message = (
                        "Chromosome %s ran out of memory on attempt %d. PostGWAS "
                        "stopped this dataset without retrying that chromosome, "
                        "because the same memory settings would likely repeat the "
                        "failure. Reduce concurrent chromosome workers, increase "
                        "execution.memory_safety_factor (adaptive mode) or "
                        "execution.memory_gb_per_chromosome (fixed mode), or provide more memory."
                        % (chrom, attempt)
                    )
                    if logger is not None:
                        logger.error(message)
                    emit_block(chrom, {
                        "attempt": attempt,
                        "status_label": "OUT OF MEMORY",
                        "elapsed": elapsed,
                        "screen_text": _log_tail(
                            configured_output_path(
                                output_dir,
                                output_layout["chromosome_log"],
                                dataset_id=sample_id,
                                chromosome=chrom,
                            )
                        ),
                        "summary": None,
                    })
                    if progress is not None:
                        progress.record_chromosome_result(
                            chrom, "failed", attempt,
                        )
                    raise MemoryError(message) from exc
                except BrokenProcessPool as exc:
                    # The worker died so hard it returned nothing.  Its log on
                    # disk is the source of truth, not the return value.
                    results[chrom] = {
                        "status": "failed",
                        "qc": {"error": "worker process died: %s" % exc},
                        "summary": None,
                    }
                    emit_block(chrom, {
                        "attempt": attempt,
                        "status_label": "WORKER DIED",
                        "elapsed": elapsed,
                        "screen_text": _log_tail(
                            configured_output_path(
                                output_dir,
                                output_layout["chromosome_log"],
                                dataset_id=sample_id,
                                chromosome=chrom,
                            )
                        ),
                        "summary": None,
                    })
                    if progress is not None:
                        progress.record_chromosome_result(
                            chrom, "failed", attempt,
                        )
                    # No further jobs can be submitted to a broken executor.
                    # The existing missing-result accounting covers queued jobs.
                    break
                except Exception as exc:
                    results[chrom] = {
                        "status": "failed",
                        "qc": {"error": "%s: %s" % (type(exc).__name__, exc)},
                        "summary": None,
                    }
                    emit_block(chrom, {
                        "attempt": attempt,
                        "status_label": "WORKER DIED",
                        "elapsed": elapsed,
                        "screen_text": _log_tail(
                            configured_output_path(
                                output_dir,
                                output_layout["chromosome_log"],
                                dataset_id=sample_id,
                                chromosome=chrom,
                            )
                        ),
                        "summary": None,
                    })
                    if progress is not None:
                        progress.record_chromosome_result(
                            chrom, "failed", attempt,
                        )
                    continue

                summary = qc.pop("__block_summary__", None)
                results[chrom_res] = {
                    "status": status,
                    "qc": qc,
                    "summary": summary,
                }
                emit_block(chrom_res, {
                    "attempt": attempt,
                    "status_label": "OK" if status == "ok" else "FAILED",
                    "elapsed": elapsed,
                    "screen_text": screen_text,
                    "summary": summary,
                })
                if progress is not None:
                    progress.record_chromosome_result(
                        chrom_res, status, attempt,
                    )
    finally:
        if screen_order == "sorted" and blocks:
            for _chrom, args in sorted(blocks, key=lambda item: _chromosome_sort_key(item[0])):
                print_chromosome_summary(**args)

    # A chromosome whose future never came back at all (a pool torn down
    # mid-round) must not silently look successful.
    for chrom in pending:
        if chrom not in results:
            results[chrom] = {
                "status": "failed",
                "qc": {"error": "the worker never returned a result"},
                "summary": None,
            }
            if logger is not None:
                logger.error(
                    "Chromosome %s produced no result at all in round %d."
                    % (chrom, round_no)
                )
            if progress is not None:
                progress.record_chromosome_result(
                    chrom, "failed", attempts.get(chrom, 0) + 1,
                )
    return results


def _prepare_dataset_for_chromosome_processing(
    *,
    sumstat_file,
    sample_column_dict,
    output_dir_path,
    resource_folder,
    default_eaf_reference_source,
    default_comparison_af_file,
    dbsnp,
    resource_layout,
    output_layout,
    gwas2vcf_input,
    executables,
    vcf_config,
    user_eaf_file,
    user_eaf_column,
    user_info_file,
    user_info_column,
    default_eaf_reference_column,
    default_comparison_af_column,
    build_reference_files,
    build_check_colmap,
    default_eaf_colmap,
    external_eaf_colmap,
    external_info_colmap,
    external_reference_staging,
    policies,
    logger,
    input_line_count,
):
    # type: (...) -> _DatasetPreparation
    """Run dataset steps 03-08 and return only chromosome-worker inputs.

    Scientific order is invariant: read, content validation, chromosome-X
    reconstruction safety, build inference, strand consensus, study-property
    resolution, exact resource preflight, and chromosome partitioning. The
    genome-wide frame is intentionally kept local and is released when this
    function returns.
    """
    sample_id = sample_column_dict["gwas_outputname"]

    # ------------------------------------------------------------------
    # 03/08  read_summary_statistics — the expensive read
    # ------------------------------------------------------------------
    with logger.step(
        3, DATASET_STEP_TOTAL, "Read the summary statistics",
        "summary_statistics_io.read_summary_statistics",
        policy_keys=[
            "input.delimiter", "input.null_values", "input.strip_double_hash_lines",
            "columns.mandatory",
            "chromosome.allowed", "allele.pattern",
            "info.multi_value_delimiter",
            "info.multi_value_aggregation",
            "info.multi_value_output_column",
            "info.multi_value_invalid_token_action",
            "duplicates.key", "duplicates.consistency_fields",
            "duplicates.quality_fields", "duplicates.selection_order",
            "duplicates.conflicting_action",
        ],
    ) as ctx:
        (
            df,
            source_snapshot,
            file_cvariant_count,
            polars_rows,
            sample_column_dict,
            removed_coords,
            removed_unsupported_chromosomes,
            non_standard_allele_count,
            removed_duplicates_count,
            removed_missing_count,
            field_completeness,
        ) = read_summary_statistics(
            sumstat_file=sumstat_file,
            output_dir=str(output_dir_path),
            sample_column_dict=sample_column_dict,
            output_layout=output_layout,
            table_delimiter=gwas2vcf_input["delimiter"],
            policies=policies,
            logger=logger,
            input_line_count=input_line_count,
        )
        ctx.set_rows(df.height)
        ctx.rows_in = polars_rows

    field_completeness_path = configured_output_path(
        output_dir_path,
        output_layout["field_completeness"],
        error_type=PipelineError,
        dataset_id=sample_id,
    )
    write_field_completeness_report(
        field_completeness_path,
        field_completeness,
        delimiter=gwas2vcf_input["delimiter"],
    )
    logger.info(
        "Initial field-completeness QC written to %s."
        % field_completeness_path
    )
    _announce(logger, "\n" + _field_completeness_block(field_completeness))

    ready_variant_types = _ready_variant_type_counts(df, sample_column_dict)
    _announce(
        logger,
        "\n" + _input_validation_summary_block(
            input_variants=file_cvariant_count,
            variants_read=polars_rows,
            missing_required=removed_missing_count,
            invalid_coordinates=removed_coords,
            unsupported_chromosomes=removed_unsupported_chromosomes,
            non_standard_alleles=non_standard_allele_count,
            duplicate_variants=removed_duplicates_count,
            ready_variants=df.height,
            ready_snps=ready_variant_types["snps"],
            ready_indels_or_other=ready_variant_types[
                "indels_or_other_variants"
            ],
        ),
    )

    # ------------------------------------------------------------------
    # 04/08  validate_content
    # ------------------------------------------------------------------
    with logger.step(
        4, DATASET_STEP_TOTAL, "Check the data behind the configuration",
        "validator.validate_content", rows_in=df.height,
        policy_keys=[
            "sample_size.missing_action",
            "sample_size.max_missing_fraction",
        ],
    ) as ctx:
        missing_sample_size_plan = prepare_missing_sample_sizes(
            df,
            sample_column_dict,
            policies=policies,
            ctx=ctx,
            scope="dataset '%s'" % sample_id,
        )
        ok, problems = validate_content(
            sample_column_dict, df, policies=policies, logger=logger,
        )
        if not ok:
            report = format_problems(
                problems,
                title="Content check for %s" % sample_id,
                columns_present=list(df.columns),
                source_name=Path(str(sumstat_file)).name,
                footer="Fix the config and re-run.",
            )
            ctx.failure_hint = (
                "The configuration names columns that the parsed data cannot support."
            )
            raise ConfigError(problems, message=report)
        ctx.info("Every configured column holds usable data.")
        ctx.extra["missing_sample_size"] = dict(missing_sample_size_plan)
        ctx.set_rows(df.height, removed=0)

    observed_chromosomes = sorted(
        (
            str(value)
            for value in df.get_column(sample_column_dict["chr_col"])
            .drop_nulls()
            .unique()
            .to_list()
        ),
        key=_chromosome_sort_key,
    )
    if not observed_chromosomes:
        raise PipelineError(
            "No chromosomes remain after dataset-level validation; resource "
            "preflight and chromosome processing cannot start."
        )

    # Z-only effect reconstruction is a study-schema decision. Enforce the
    # chromosome-X policy immediately after parsed-content validation, before
    # build/strand analysis, resource preflight, partition I/O, or worker
    # launch. The same shared guard remains inside chromosome step 06 as a
    # defence for direct and alternative callers.
    logger.settings(["effect_from_z.x_chromosome_z_only_action"])
    x_z_only_policy = None
    try:
        if "X" in observed_chromosomes:
            x_z_only_policy = enforce_x_chromosome_z_only_policy(
                "X",
                df.columns,
                sample_column_dict,
                policies=policies,
                logger=logger,
                scope="Dataset '{}'".format(sample_id),
            )
    except XChromosomeZOnlyReconstructionError as exc:
        message = "%s No chromosome analysis was started." % exc
        logger.error(message)
        _announce(
            logger,
            "\n" + "\n".join([
                screen_line(
                    "error", "Chromosome-X effect reconstruction blocked",
                    indent=4,
                ),
                screen_field(
                    "error", "Reason", str(exc), indent=8, label_width=22,
                ),
                screen_field(
                    "error", "Chromosome analysis",
                    "not started; no chromosome-X BETA or SE was calculated",
                    indent=8, label_width=22,
                ),
            ]),
            marker="ERROR",
        )
        raise PipelineError(message) from exc

    if x_z_only_policy is not None and x_z_only_policy["applies"]:
        _announce(
            logger,
            "\n" + "\n".join([
                screen_line(
                    "warn", "Chromosome-X reconstruction override", indent=4,
                ),
                screen_field(
                    "warn", "Configured action",
                    x_z_only_policy["action"], indent=8, label_width=22,
                ),
                screen_field(
                    "warn", "Scientific meaning",
                    "autosomal 2*EAF*(1-EAF) variance is being assumed for X; "
                    "the reconstructed BETA and SE are assumption-dependent",
                    indent=8, label_width=22,
                ),
            ]),
            marker="WARNING",
        )

    # ------------------------------------------------------------------
    # 05/08  infer_genome_build
    # ------------------------------------------------------------------
    genome_build_info = infer_genome_build(
        df,
        build_reference_files,
        build_check_colmap,
        sample_column_dict,
        logger=logger,
        policies=policies,
        step_number=5,
        step_total=DATASET_STEP_TOTAL,
    )
    grch_version = genome_build_info["inferred_build"]
    build_label_width = 26
    build_match_lines = [
        screen_field(
            "genetic",
            "%s matches" % build,
            "%s (%s%% input; %s%% testable)" % (
                _count(genome_build_info["matches"][build]),
                genome_build_info["input_percentages"][build],
                genome_build_info["percentages"][build],
            ),
            indent=8,
            label_width=build_label_width,
        )
        for build in build_reference_files
    ]
    build_reference_lines = [
        screen_field(
            "genetic",
            "%s reference coverage" % build,
            "%s / %s markers (%s%%)" % (
                _count(
                    genome_build_info["matched_reference_marker_counts"][build]
                ),
                _count(genome_build_info["reference_marker_counts"][build]),
                genome_build_info["reference_percentages"][build],
            ),
            indent=8,
            label_width=build_label_width,
        )
        for build in build_reference_files
    ]
    _announce(
        logger,
        "\n" + "\n".join([
            screen_line("genetic", "Genome-build inference", indent=4),
            screen_field(
                "success", "Selected build", grch_version,
                indent=8, label_width=build_label_width,
            ),
            screen_field(
                "count", "Input variants",
                _count(genome_build_info["input_variants"]),
                indent=8, label_width=build_label_width,
            ),
            screen_field(
                "count", "Coordinate-testable",
                _count(genome_build_info["testable_variants"]),
                indent=8, label_width=build_label_width,
            ),
            screen_field(
                "info", "Not build-testable",
                _count(genome_build_info["untestable_variants"]),
                indent=8, label_width=build_label_width,
            ),
            screen_field(
                "info", "Reference file",
                build_reference_files.get(grch_version, "not selected"),
                indent=8, label_width=build_label_width,
            ),
        ] + build_match_lines + build_reference_lines),
    )
    if grch_version == "Ambiguous":
        raise PipelineError(
            "❌ Genome build inference failed (Ambiguous).\n"
            "   Reason: %s\n"
            "   Set 'build.mode' to one of the configured genome builds in YAML.\n"
            "   Aborting pipeline to avoid invalid gwas2vcf resource usage."
            % genome_build_info["ambiguous_reason"]
        )

    # ------------------------------------------------------------------
    # 06/08  strand_detection
    # ------------------------------------------------------------------
    with logger.step(
        6, DATASET_STEP_TOTAL, "Strand detection",
        "strand.resolve_strand_consensus", rows_in=df.height,
        policy_keys=[
            "strand.mode", "strand.consensus_threshold",
            "strand.min_informative_variants",
        ],
    ) as ctx:
        strand_decision = resolve_strand_consensus(
            genome_build_info, grch_version, policies=policies,
        )
        ctx.decide(
            "Study strand consensus",
            {
                "forward non-palindromic matches": strand_decision["forward"],
                "reverse non-palindromic matches": strand_decision["reverse"],
                "ambiguous matches excluded": strand_decision["ambiguous"],
                "dominant fraction": strand_decision["dominant_fraction"],
            },
            "%s - %s" % (strand_decision["strand"], strand_decision["reason"]),
        )
        ctx.extra["strand_consensus"] = dict(strand_decision)
        ctx.set_rows(df.height, removed=0)

    # ------------------------------------------------------------------
    # 07/08  resolve_study_properties
    # ------------------------------------------------------------------
    study_decisions = resolve_study_properties(
        df,
        sample_column_dict,
        policies=policies,
        logger=logger,
        step_number=7,
        step_total=DATASET_STEP_TOTAL,
    )
    study_decisions.update(strand_decision)
    dataset_study_decisions = {
        "effect_type": study_decisions.get("effect_type"),
        "effect_type_source": study_decisions.get("effect_type_source"),
        "effect_type_detected": study_decisions.get("effect_type_detected"),
        "effect_type_matches_declaration": study_decisions.get(
            "effect_type_matches_declaration"
        ),
        "se_scale": study_decisions.get("se_scale"),
        "se_scale_source": study_decisions.get("se_scale_source"),
        "se_scale_detected": study_decisions.get("se_scale_detected"),
        "se_scale_matches_declaration": study_decisions.get(
            "se_scale_matches_declaration"
        ),
        "se_scale_evidence": study_decisions.get("se_scale_evidence") or {},
        "pvalue_type": study_decisions.get("pvalue_type"),
        "pvalue_type_source": study_decisions.get("pvalue_type_source"),
        "pvalue_type_detected": study_decisions.get("pvalue_type_detected"),
        "pvalue_type_matches_declaration": study_decisions.get(
            "pvalue_type_matches_declaration"
        ),
        "declaration_mismatches": list(
            study_decisions.get("declaration_mismatches") or []
        ),
        "eaf_is_maf": study_decisions.get("eaf_is_maf"),
        "eaf_is_maf_source": study_decisions.get("eaf_is_maf_source"),
        "strand_mode": study_decisions.get("strand_mode"),
        "strand": study_decisions.get("strand"),
        "strand_consensus": {
            key: study_decisions.get(key)
            for key in (
                "forward", "reverse", "ambiguous", "informative",
                "dominant_fraction", "consensus_threshold", "reason",
            )
        },
    }
    _announce(logger, "\n" + _study_decisions_block(study_decisions))
    _announce(
        logger,
        "\n" + _validated_column_mapping_block(
            sample_column_dict, study_decisions,
        ),
    )

    if x_z_only_policy is not None:
        dataset_study_decisions["x_chromosome_z_only_policy"] = dict(
            x_z_only_policy
        )

    # The exact resources cannot be known at command preflight because both the
    # genome build and chromosomes present are data-derived. Validate all of
    # them before partition I/O or worker launch.
    try:
        resource_maps, resource_preflight = preflight_harmonisation_resources(
            observed_chromosomes,
            grch_version=grch_version,
            target_build=vcf_config["target_builds"][grch_version],
            resource_folder=resource_folder,
            user_eaf_file=user_eaf_file,
            default_eaf_reference_source=default_eaf_reference_source,
            default_comparison_af_file=default_comparison_af_file,
            resource_layout=resource_layout,
            user_info_file=user_info_file,
            dbsnp=dbsnp,
            user_eaf_column=user_eaf_column,
            default_eaf_reference_column=default_eaf_reference_column,
            default_comparison_af_column=default_comparison_af_column,
            user_info_column=user_info_column,
            default_eaf_colmap=default_eaf_colmap,
            external_eaf_colmap=external_eaf_colmap,
            external_info_colmap=external_info_colmap,
            executables=executables,
            vcf_config=vcf_config,
            policies=policies,
            require_default_eaf=True,
            logger=logger,
        )
    except ResourcePreflightError as exc:
        message = (
            "%s\nNo chromosome analysis was started. Fix the listed resource "
            "files, indexes, schemas, or path templates and rerun."
            % exc
        )
        logger.error(message)
        failure_lines = [
            screen_line("error", "Resource preflight failed", indent=4),
            screen_field(
                "error", "Problems", _count(len(exc.issues)),
                indent=8, label_width=26,
            ),
        ]
        for issue in exc.issues:
            failure_lines.append(screen_field(
                "error",
                issue["resource"],
                "chromosome%s %s; %s; %s"
                % (
                    "" if len(issue["chromosomes"]) == 1 else "s",
                    ", ".join(issue["chromosomes"]),
                    issue["path"],
                    issue["problem"],
                ),
                indent=8,
                label_width=26,
            ))
        failure_lines.append(screen_field(
            "error", "Chromosome analysis",
            "not started; complete details are in %s" % logger.log_path,
            indent=8, label_width=26,
        ))
        _announce(logger, "\n" + "\n".join(failure_lines), marker="ERROR")
        raise ConfigError([], message=message) from exc

    logger.info(
        "Chromosome resource preflight passed: %s"
        % json.dumps(resource_preflight, sort_keys=True)
    )
    _announce(
        logger,
        "\n" + "\n".join([
            screen_line("genetic", "Resource preflight", indent=4),
            screen_field(
                "success", "Chromosomes covered",
                "%d (%s)" % (
                    len(observed_chromosomes),
                    ", ".join(observed_chromosomes),
                ),
                indent=8, label_width=24,
            ),
            screen_field(
                "success", "Resource files",
                "%d present, non-empty and readable"
                % resource_preflight["unique_resource_files_checked"],
                indent=8, label_width=24,
            ),
            screen_field(
                "success", "VCF resources",
                "%d readable indexes; %d frequency headers checked"
                % (
                    resource_preflight["indexed_vcfs_checked"],
                    resource_preflight["frequency_vcf_headers_checked"],
                ),
                indent=8, label_width=24,
            ),
            screen_field(
                "success", "FASTA indexes",
                "%d readable and compatible with study/chain chromosome labels"
                % resource_preflight["fasta_indexes_checked"],
                indent=8, label_width=24,
            ),
            screen_field(
                "success", "Reference schemas",
                "%d table header%s and %d GFF/chain file%s checked"
                % (
                    resource_preflight["table_headers_checked"],
                    "" if resource_preflight["table_headers_checked"] == 1 else "s",
                    resource_preflight["structured_text_files_checked"],
                    "" if resource_preflight["structured_text_files_checked"] == 1 else "s",
                ),
                indent=8, label_width=24,
            ),
        ]),
    )

    # ------------------------------------------------------------------
    # 08/08  split by chromosome and stage shared user references
    # ------------------------------------------------------------------
    with logger.step(
        8, DATASET_STEP_TOTAL, "Split by chromosome",
        "chromosome_partition.write_chromosome_partitions + "
        "external_reference_staging.stage_shared_external_reference_files + "
        "imputation_quality.resolve_dataset_info_score_type",
        rows_in=df.height,
        policy_keys=[
            "input.chromosome_partition_compression",
            "info.score_type",
            "info.clip_min",
            "info.clip_max",
            "info.clip_tolerance",
            "info.mach_rsq_max",
            "info.auto_mach_rsq_fraction",
            "info.maximum_invalid_fraction",
        ],
    ) as ctx:
        info_score_decision = None
        if sample_column_dict.get("info_source") in ("internal", "fixed_cli"):
            info_score_decision = resolve_dataset_info_score_type(
                df=df,
                sample_column_dict=sample_column_dict,
                resource_maps=resource_maps,
                chromosomes=observed_chromosomes,
                external_info_column=user_info_column,
                external_info_colmap=external_info_colmap,
                external_reference_staging=external_reference_staging,
                policies=policies,
            )
        chromosome_partitions = write_chromosome_partitions(
            df,
            sample_gwas_dict=sample_column_dict,
            output_layout=output_layout,
            compression=str(
                policies.get("input.chromosome_partition_compression")
            ),
            source_snapshot=(
                source_snapshot if bool(policies.get("rejects.enabled")) else None
            ),
        )
        chr_file_by_chrom = chromosome_partitions.paths
        partition_rows = chromosome_partitions.rows
        ready_variants = df.height
        del source_snapshot
        del df
        ctx.info(
            "Wrote %d typed Parquet chromosome work file%s: %s"
            % (
                len(chr_file_by_chrom),
                "" if len(chr_file_by_chrom) == 1 else "s",
                ", ".join(sorted(chr_file_by_chrom, key=_chromosome_sort_key)),
            )
        )
        try:
            resource_maps, staging_summary = (
                stage_shared_external_reference_files(
                    resource_maps,
                    chromosomes=observed_chromosomes,
                    dataset_id=sample_id,
                    output_directory=output_dir_path,
                    output_layout=output_layout,
                    user_eaf_specification=user_eaf_file,
                    user_eaf_column=user_eaf_column,
                    external_eaf_mapping=external_eaf_colmap,
                    user_info_specification=user_info_file,
                    user_info_column=user_info_column,
                    external_info_mapping=external_info_colmap,
                    settings=external_reference_staging,
                    policies=policies,
                    logger=logger,
                )
            )
        except ExternalReferenceStagingError as exc:
            raise PipelineError(
                "%s No chromosome worker was started; fix the user EAF/INFO "
                "file or external-reference staging configuration and rerun."
                % exc
            ) from exc
        resource_preflight["external_reference_staging"] = staging_summary
        if staging_summary["status"] == "staged":
            ctx.info(
                "Read %d shared whole-genome external reference file%s once "
                "and wrote projected partitions for %d chromosomes."
                % (
                    staging_summary["source_files_read"],
                    "" if staging_summary["source_files_read"] == 1 else "s",
                    len(observed_chromosomes),
                )
            )
        else:
            ctx.info(
                "External-reference staging was not needed; chromosome path "
                "templates and single-chromosome inputs remain unchanged."
            )
        if info_score_decision is None:
            info_score_decision = resolve_dataset_info_score_type(
                df=None,
                sample_column_dict=sample_column_dict,
                resource_maps=resource_maps,
                chromosomes=observed_chromosomes,
                external_info_column=user_info_column,
                external_info_colmap=external_info_colmap,
                external_reference_staging=external_reference_staging,
                policies=policies,
            )
        resolved_info_score_type = info_score_decision["resolved_type"]
        info_score_type_source = info_score_decision["decision_source"]
        sample_column_dict["resolved_info_score_type"] = resolved_info_score_type
        sample_column_dict["info_score_type_source"] = info_score_type_source
        study_decisions["info_score_type"] = resolved_info_score_type
        study_decisions["info_score_type_source"] = info_score_type_source
        study_decisions["info_score_type_evidence"] = dict(info_score_decision)
        dataset_study_decisions["info_score_type"] = resolved_info_score_type
        dataset_study_decisions[
            "info_score_type_source"
        ] = info_score_type_source
        dataset_study_decisions[
            "info_score_type_evidence"
        ] = dict(info_score_decision)
        ctx.decide(
            "Dataset-wide INFO score type",
            {
                "selected source": info_score_decision["source"],
                "finite numeric values": info_score_decision["finite_values"],
                "standard INFO range": info_score_decision[
                    "standard_range_values"
                ],
                "possible MaCH Rsq range": info_score_decision[
                    "mach_range_values"
                ],
                "MaCH-range fraction": info_score_decision[
                    "mach_range_fraction"
                ],
                "above MaCH maximum": info_score_decision[
                    "above_mach_maximum"
                ],
                "configured type": info_score_decision["configured_type"],
            },
            "%s (%s)" % (resolved_info_score_type, info_score_type_source),
        )
        if info_score_decision["mach_range_values"]:
            ctx.warn(
                "{:,} finite INFO-source values ({:.6%}) are above the "
                "standard-INFO tolerance {:g} and no higher than the configured "
                "MaCH Rsq maximum {:g}; the single dataset-wide decision is {}."
                .format(
                    info_score_decision["mach_range_values"],
                    info_score_decision["mach_range_fraction"],
                    info_score_decision["standard_info_maximum"],
                    info_score_decision["mach_rsq_maximum"],
                    resolved_info_score_type,
                )
            )
        if info_score_decision["unusable_values"]:
            ctx.warn(
                "{:,} rows in the selected INFO source were missing, non-numeric, "
                "NaN, or infinite and therefore did not enter the finite-value "
                "score-type denominator. Chromosome step 11 converts any matched "
                "unusable score to missing and applies info.on_missing."
                .format(info_score_decision["unusable_values"])
            )
        if info_score_decision["above_mach_maximum"]:
            ctx.warn(
                "{:,} finite INFO-source values ({:.6%}) exceed the configured "
                "MaCH Rsq maximum {:g}. Their fraction did not exceed "
                "info.maximum_invalid_fraction, so chromosome step 11 will "
                "apply info.out_of_range to those rows."
                .format(
                    info_score_decision["above_mach_maximum"],
                    info_score_decision["invalid_fraction"],
                    info_score_decision["mach_rsq_maximum"],
                )
            )
        ctx.extra["external_reference_staging"] = staging_summary
        ctx.extra["info_score_type"] = dict(info_score_decision)
        ctx.set_rows(ready_variants, removed=0)

    all_chromosomes = sorted(chr_file_by_chrom, key=_chromosome_sort_key)
    if not all_chromosomes:
        raise PipelineError(
            "No per-chromosome files found after splitting. "
            "Check that 'CHR' and 'POS' columns are correctly mapped."
        )
    if all_chromosomes != observed_chromosomes:
        raise PipelineError(
            "Chromosome partitioning produced a different chromosome set from "
            "resource preflight. Preflight=%s; partitions=%s. No chromosome "
            "worker was started."
            % (", ".join(observed_chromosomes), ", ".join(all_chromosomes))
        )

    return _DatasetPreparation(
        sample_columns=sample_column_dict,
        chromosome_files=chr_file_by_chrom,
        resource_maps=resource_maps,
        genome_build=grch_version,
        genome_build_info=genome_build_info,
        study_decisions=study_decisions,
        dataset_study_decisions=dataset_study_decisions,
        missing_sample_size_plan=dict(missing_sample_size_plan),
        resource_preflight=resource_preflight,
        field_completeness=field_completeness,
        field_completeness_path=Path(field_completeness_path),
        input_variants=int(file_cvariant_count),
        rows_read=int(polars_rows),
        removed_invalid_coordinates=int(removed_coords),
        removed_unsupported_chromosomes=int(
            removed_unsupported_chromosomes
        ),
        removed_non_standard_alleles=int(non_standard_allele_count),
        removed_duplicates=int(removed_duplicates_count),
        removed_missing_values=int(removed_missing_count),
        ready_variants=int(ready_variants),
        ready_snps=int(ready_variant_types["snps"]),
        ready_indels_or_other=int(
            ready_variant_types["indels_or_other_variants"]
        ),
        partition_rows=dict(partition_rows),
    )


def harmonise_chromosomes(
    sumstat_file,
    sample_column_dict,
    output_dir,
    resource_folder,
    default_eaf_reference_source,
    default_comparison_af_file,
    dbsnp,
    resource_layout,
    output_layout,
    gwas2vcf_input,
    executables,
    vcf_config,
    max_workers,
    user_eaf_file,
    user_eaf_column,
    user_info_file,
    user_info_column,
    default_eaf_reference_column,
    default_comparison_af_column,
    gwas2vcf_main_script_path,
    build_reference_files,
    build_check_colmap,
    default_eaf_colmap,
    external_eaf_colmap,
    external_info_colmap,
    external_reference_staging,
    policies,
    logger,
    input_line_count=None,
    progress=None,
):
    """Dataset steps 03-08, then the per-chromosome fan-out with retry rounds.

    Returns ``(per_chr_qc, grch_version)`` exactly as before.  The dataset-level
    extras (status, reject files, study decisions, parallelism) are attached
    under the ``"__dataset__"`` key, which ``run_harmonisation_pipeline`` pops
    before the QC JSON is written.

    Raises :class:`PipelineError` when a chromosome cannot be completed and
    ``execution.fail_dataset_on_chr_error`` is true (the default).  Previously
    the errors were collected, printed, and the function returned normally, so
    the pipeline reported success on an incomplete genome.
    """
    pol = policies if policies is not None else default_policies()
    sample_id = sample_column_dict["gwas_outputname"]
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    log_dir = configured_output_path(
        output_dir_path, output_layout["logs_directory"], error_type=PipelineError,
    )
    log_dir.mkdir(parents=True, exist_ok=True)

    # This function is normally called with the dataset logger already built by
    # run_harmonisation_pipeline; building one here keeps it usable on its own
    # rather than crashing on `None.step(...)`.
    if logger is None:
        logger = PipelineLogger(
            sample_id=sample_id,
            scope="dataset",
            log_dir=str(log_dir),
            policies=pol,
            level=None,
            screen_level=None,
        )

    _announce(
        logger,
        "\n" + "\n".join([
            screen_line("run", "Starting harmonisation and GWAS-to-VCF", indent=4),
            screen_field(
                "info", "Input file", sumstat_file,
                indent=8, label_width=18,
            ),
        ]),
    )

    per_chr_qc = {}  # type: Dict[str, Any]
    dataset_extra = {
        "status": "FAILED",
        "chromosomes": {},
        "reject_files": [],
    }  # type: Dict[str, Any]
    prepared = _prepare_dataset_for_chromosome_processing(
        sumstat_file=sumstat_file,
        sample_column_dict=sample_column_dict,
        output_dir_path=output_dir_path,
        resource_folder=resource_folder,
        default_eaf_reference_source=default_eaf_reference_source,
        default_comparison_af_file=default_comparison_af_file,
        dbsnp=dbsnp,
        resource_layout=resource_layout,
        output_layout=output_layout,
        gwas2vcf_input=gwas2vcf_input,
        executables=executables,
        vcf_config=vcf_config,
        user_eaf_file=user_eaf_file,
        user_eaf_column=user_eaf_column,
        user_info_file=user_info_file,
        user_info_column=user_info_column,
        default_eaf_reference_column=default_eaf_reference_column,
        default_comparison_af_column=default_comparison_af_column,
        build_reference_files=build_reference_files,
        build_check_colmap=build_check_colmap,
        default_eaf_colmap=default_eaf_colmap,
        external_eaf_colmap=external_eaf_colmap,
        external_info_colmap=external_info_colmap,
        external_reference_staging=external_reference_staging,
        policies=pol,
        logger=logger,
        input_line_count=input_line_count,
    )
    sample_column_dict = prepared.sample_columns
    chr_file_by_chrom = prepared.chromosome_files
    resource_maps = prepared.resource_maps
    grch_version = prepared.genome_build
    genome_build_info = prepared.genome_build_info
    study_decisions = prepared.study_decisions
    missing_sample_size_plan = prepared.missing_sample_size_plan
    field_completeness = prepared.field_completeness
    field_completeness_path = prepared.field_completeness_path
    partition_rows = prepared.partition_rows
    all_chromosomes = sorted(chr_file_by_chrom, key=_chromosome_sort_key)
    n_chr = len(all_chromosomes)
    if progress is not None:
        progress.start_chromosomes(all_chromosomes)

    dataset_extra["field_completeness"] = field_completeness
    dataset_extra["missing_sample_size"] = dict(missing_sample_size_plan)
    dataset_extra["study_decisions"] = dict(
        prepared.dataset_study_decisions
    )
    dataset_extra["resource_preflight"] = prepared.resource_preflight
    schedule = None
    if pol.get("execution.scheduling_mode") == "adaptive":
        try:
            schedule = plan_chromosome_schedule(
                pol, chr_file_by_chrom, partition_rows, resource_maps,
                requested_workers=max_workers, output_dir=output_dir_path,
                output_layout=output_layout, sample_id=sample_id,
            )
        except (ValueError, OSError) as exc:
            raise PipelineError("Cannot schedule chromosome processing: %s" % exc) from exc
        workers = schedule["workers"]
        threads_per_chromosome = schedule["threads_per_chromosome"]
        logger.decide(
            "Adaptive chromosome scheduling", schedule,
            "large/small admission within CPU and estimated memory budgets",
        )
    else:
        workers, threads_per_chromosome = _derive_parallelism(
            pol, max_workers, n_chr, logger=logger
        )
    dataset_extra["parallelism"] = {
        "mode": pol.get("execution.scheduling_mode"),
        "workers": workers,
        "threads_per_chromosome": threads_per_chromosome,
        "polars_threads_per_worker": threads_per_chromosome,
        "requested_workers": max_workers,
        "total_cpu_budget": pol.get("execution.total_cpu_budget"),
        "memory_budget_gb": pol.get("execution.memory_budget_gb"),
        "memory_gb_per_chromosome": pol.get(
            "execution.memory_gb_per_chromosome"
        ),
        **(schedule or {}),
    }
    _announce(
        logger,
        "\n" + "\n".join([
            screen_line(
                "genetic", "Chromosome processing", indent=4,
            ),
            screen_field(
                "count", "Chromosomes", "%d (%s)" % (
                    n_chr, ", ".join(all_chromosomes),
                ), indent=8, label_width=20,
            ),
            screen_field(
                "count", "Parallel work",
                "up to %d at a time; %d Polars/bcftools thread%s each" % (
                    workers, threads_per_chromosome,
                    "" if threads_per_chromosome == 1 else "s",
                ), indent=8, label_width=20,
            ),
            screen_field(
                "count", "Memory scheduling",
                ("adaptive; %.2f GB for workers; %.2f GB parent/headroom; "
                 "large and small chromosomes mixed" % (
                     schedule["worker_memory_budget_gb"], schedule["memory_headroom_gb"],
                 )) if schedule else "fixed per-chromosome reservations",
                indent=8, label_width=20,
            ),
        ]),
    )

    # ------------------------------------------------------------------
    # Retry rounds (plan Part 7)
    # ------------------------------------------------------------------
    worker_kwargs = {
        "sample_column_dict": sample_column_dict,
        "resource_folder": resource_folder,
        "grch_version": grch_version,
        "user_eaf_file": user_eaf_file,
        "default_eaf_reference_source": default_eaf_reference_source,
        "default_comparison_af_file": default_comparison_af_file,
        "user_info_file": user_info_file,
        "user_eaf_column": user_eaf_column,
        "default_eaf_reference_column": default_eaf_reference_column,
        "default_comparison_af_column": default_comparison_af_column,
        "user_info_column": user_info_column,
        "dbsnp": dbsnp,
        "resource_layout": resource_layout,
        "output_layout": output_layout,
        "gwas2vcf_input": gwas2vcf_input,
        "executables": executables,
        "vcf_config": vcf_config,
        "output_dir": str(output_dir_path),
        "threads": threads_per_chromosome,
        "gwas2vcf_main_script_path": gwas2vcf_main_script_path,
        "external_eaf_colmap": external_eaf_colmap,
        "default_eaf_colmap": default_eaf_colmap,
        "external_info_colmap": external_info_colmap,
        "policies": pol,
        "study_decisions": study_decisions,
        "missing_sample_size_plan": missing_sample_size_plan,
    }

    max_retry_rounds = int(pol.get("execution.max_retry_rounds"))
    clean_before_retry = bool(pol.get("execution.clean_before_retry"))
    retry_requires_progress = bool(pol.get("execution.retry_requires_progress"))
    screen_order = str(pol.get("logging.screen_order"))

    pending = list(all_chromosomes)
    attempts = dict((chrom, 0) for chrom in all_chromosomes)
    statuses = dict((chrom, "not attempted") for chrom in all_chromosomes)
    summaries = {}  # type: Dict[str, Any]
    round_no = 0

    while pending:
        round_no += 1
        if round_no == 1:
            _announce(
                logger,
                "\n" + screen_line(
                    "run", "Round 1 — all %d chromosomes" % len(pending), indent=4,
                ),
            )
        else:
            _announce(
                logger,
                "\n" + screen_line(
                    "retry", "Retry round %d — %d chromosome%s (%s)" % (
                        round_no, len(pending),
                        "" if len(pending) == 1 else "s", ", ".join(pending),
                    ), indent=4,
                ),
            )
            if clean_before_retry:
                for chrom in pending:
                    removed = remove_partial_chromosome_outputs(
                        output_dir_path,
                        sample_id,
                        chrom,
                        output_layout,
                        logger=logger,
                    )
                    _announce(
                        logger,
                        screen_field(
                            "loss", "Retry cleanup",
                            "%d partial file%s removed for chromosome %s" % (
                                removed, "" if removed == 1 else "s", chrom,
                            ), indent=6, label_width=20,
                        ),
                    )
            else:
                _announce(
                    logger,
                    screen_field(
                        "warning", "Retry cleanup",
                        "disabled; partial VCFs from the previous round remain on disk",
                        indent=6, label_width=20,
                    ),
                    marker="WARNING",
                )

        round_results = _run_one_round(
            pending=pending,
            chr_file_by_chrom=chr_file_by_chrom,
            resource_maps=resource_maps,
            round_no=round_no,
            attempts=attempts,
            worker_kwargs=worker_kwargs,
            max_workers=workers,
            sample_id=sample_id,
            output_dir=output_dir_path,
            output_layout=output_layout,
            screen_order=screen_order,
            logger=logger,
            progress=progress,
            schedule=schedule,
        )

        for chrom in pending:
            attempts[chrom] += 1
            outcome = round_results.get(chrom, {})
            per_chr_qc[chrom] = outcome.get("qc", {"error": "no result"})
            if outcome.get("summary") is not None:
                summaries[chrom] = outcome["summary"]
            statuses[chrom] = "ok" if outcome.get("status") == "ok" else "failed"

        succeeded = [c for c in pending if statuses[c] == "ok"]
        failed = [c for c in pending if statuses[c] != "ok"]

        for chrom in succeeded:
            if attempts[chrom] > 1:
                _announce(
                    logger,
                    screen_line(
                        "success", "Chromosome %s recovered on attempt %d" % (
                            chrom, attempts[chrom],
                        ), indent=6,
                    ),
                )

        if not failed:
            break
        if round_no > max_retry_rounds:
            _announce(
                logger,
                screen_field(
                    "error", "Stopping",
                    "maximum retry rounds %d reached after round %d; %d chromosome%s failed: %s"
                    % (max_retry_rounds, round_no, len(failed),
                       "" if len(failed) == 1 else "s", ", ".join(failed)),
                    indent=4, label_width=16,
                ),
                marker="WARNING",
            )
            break
        if round_no > 1 and retry_requires_progress and not succeeded:
            _announce(
                logger,
                screen_field(
                    "error", "Stopping",
                    "retry round %d recovered nothing; another round would repeat the same failure"
                    % round_no,
                    indent=4, label_width=16,
                ),
                marker="WARNING",
            )
            break
        pending = failed

    completed = [c for c in all_chromosomes if statuses[c] == "ok"]
    failed = [c for c in all_chromosomes if statuses[c] != "ok"]

    # ------------------------------------------------------------------
    # Rejected variants: one current, verified file per contributing stage ->
    # one dataset file. Row and reason counts are collected by the same bounded
    # stream that writes the combined provenance output.
    # ------------------------------------------------------------------
    if bool(pol.get("rejects.enabled")):
        provenance_failures = [
            chrom for chrom in all_chromosomes
            if (per_chr_qc.get(chrom) or {}).get("reject_output_error")
        ]
        if provenance_failures:
            raise PipelineError(
                "Required rejected-variant provenance could not be written for "
                "chromosome%s %s. The dataset cannot be finalised as successful "
                "or partial because removed variants would not be auditable."
                % (
                    "" if len(provenance_failures) == 1 else "s",
                    ", ".join(provenance_failures),
                ),
                results=per_chr_qc,
                chromosomes=provenance_failures,
                status="FAILED",
            )

        reject_paths = []  # type: List[str]
        reject_source_labels = {}  # type: Dict[str, str]
        input_reject = resolved_reject_output_path(
            configured_output_path(
                output_dir_path,
                output_layout["input_reject"],
                dataset_id=sample_id,
            ),
            pol.get("rejects.compress"),
        )
        if not Path(input_reject).is_file():
            raise PipelineError(
                "The required current input-stage rejected-variants file is "
                "missing at %s. "
                "The dataset cannot be finalised without rejection provenance."
                % input_reject,
                results=per_chr_qc,
                status="FAILED",
            )
        reject_paths.append(input_reject)
        reject_source_labels[input_reject] = "input"

        for chrom in all_chromosomes:
            reject_metadata = (
                (per_chr_qc.get(chrom) or {}).get("rejects") or {}
            )
            current_reject_flushed = reject_metadata.get("flushed") is True
            expected_reject = resolved_reject_output_path(
                configured_output_path(
                    output_dir_path,
                    output_layout["chromosome_reject"],
                    dataset_id=sample_id,
                    chromosome=chrom,
                ),
                pol.get("rejects.compress"),
            )
            if not current_reject_flushed:
                if chrom in completed:
                    raise PipelineError(
                        "Chromosome %s completed without confirmed current-attempt "
                        "rejection provenance. The dataset cannot be finalised."
                        % chrom,
                        results=per_chr_qc,
                        chromosomes=[chrom],
                        status="FAILED",
                    )
                if Path(expected_reject).exists():
                    logger.warn(
                        "Ignoring unconfirmed reject shard %s for failed "
                        "chromosome %s; it may belong to an earlier attempt."
                        % (expected_reject, chrom)
                    )
                continue
            reported_reject = reject_metadata.get("path")
            if (
                reported_reject is None
                or os.path.abspath(str(reported_reject))
                != os.path.abspath(expected_reject)
            ):
                raise PipelineError(
                    "Chromosome %s reported current rejection provenance at %r; "
                    "the configured current-attempt path is %s."
                    % (chrom, reported_reject, expected_reject),
                    results=per_chr_qc,
                    chromosomes=[chrom],
                    status="FAILED",
                )
            if not Path(expected_reject).is_file():
                raise PipelineError(
                    "Chromosome %s confirmed a flushed rejection shard, but "
                    "the file is missing at %s."
                    % (chrom, expected_reject),
                    results=per_chr_qc,
                    chromosomes=[chrom],
                    status="FAILED",
                )
            reject_paths.append(expected_reject)
            reject_source_labels[expected_reject] = chrom

        try:
            combined = concat_reject_files(
                reject_paths,
                str(configured_output_path(
                    output_dir_path,
                    output_layout["dataset_reject"],
                    dataset_id=sample_id,
                )),
                delimiter=gwas2vcf_input["delimiter"],
                batch_rows=pol.get("rejects.concat_batch_rows"),
                logger=logger,
                compress=True,
                remove_sources=True,
            )
        except RejectOutputError as exc:
            raise PipelineError(
                "%s The dataset cannot be finalised without complete rejection "
                "provenance." % exc,
                results=per_chr_qc,
                status="FAILED",
            ) from exc
        observed_source_paths = set(combined.get("rows_by_source") or {})
        if observed_source_paths != set(reject_source_labels):
            raise PipelineError(
                "Combined rejection provenance returned source counts for %s; "
                "expected %s."
                % (
                    ", ".join(sorted(observed_source_paths)),
                    ", ".join(sorted(reject_source_labels)),
                ),
                results=per_chr_qc,
                status="FAILED",
            )
        reject_rows_by_source = dict(
            (source, 0) for source in ["input"] + all_chromosomes
        )
        counts_by_source = dict(
            (source, {}) for source in ["input"] + all_chromosomes
        )  # type: Dict[str, Dict[str, int]]
        combined_rows_by_path = combined.get("rows_by_source") or {}
        combined_reasons_by_path = (
            combined.get("reject_counts_by_source") or {}
        )
        for path, source in reject_source_labels.items():
            reject_rows_by_source[source] = _required_row_count(
                combined_rows_by_path.get(path),
                "combined reject rows for %s" % source,
            )
            source_reasons = combined_reasons_by_path.get(path)
            if not isinstance(source_reasons, dict):
                raise PipelineError(
                    "Combined rejection provenance did not return reason "
                    "counts for %s (%s)." % (source, path),
                    results=per_chr_qc,
                    status="FAILED",
                )
            validated_reasons = {}
            for reason, count in source_reasons.items():
                if reason not in REASONS:
                    raise PipelineError(
                        "Combined rejection provenance returned the "
                        "unregistered reason %r for %s (%s)."
                        % (reason, source, path),
                        results=per_chr_qc,
                        status="FAILED",
                    )
                validated_reasons[reason] = _required_row_count(
                    count,
                    "combined reject reason %s for %s" % (reason, source),
                )
            if sum(validated_reasons.values()) != reject_rows_by_source[source]:
                raise PipelineError(
                    "Combined rejection provenance reports %d row(s) for %s, "
                    "but its reason counts total %d."
                    % (
                        reject_rows_by_source[source],
                        source,
                        sum(validated_reasons.values()),
                    ),
                    results=per_chr_qc,
                    status="FAILED",
                )
            counts_by_source[source] = validated_reasons

        try:
            dataset_extra["reconciliation"] = _reconcile_dataset_rows(
                parsed_rows=prepared.rows_read,
                ready_rows=prepared.ready_variants,
                partition_rows=partition_rows,
                completed=completed,
                failed=failed,
                chromosome_summaries=summaries,
                reject_rows_by_source=reject_rows_by_source,
                combined_reject_rows=combined.get("rows"),
                logger=logger,
            )
        except PipelineError as exc:
            raise PipelineError(
                str(exc),
                results=per_chr_qc,
                chromosomes=failed,
                status="FAILED",
            ) from exc

        reject_reason_table, reject_counts = _write_reject_reason_table(
            configured_output_path(
                output_dir_path,
                output_layout["reject_summary"],
                dataset_id=sample_id,
            ),
            counts_by_source,
            gwas2vcf_input["delimiter"],
            logger=logger,
        )
        dataset_extra["reject_reason_table"] = reject_reason_table
        dataset_extra["reject_counts"] = reject_counts
        dataset_extra["reject_files"] = [combined.get("path")]
        dataset_extra["rejected_variants_file"] = combined.get("path")
        dataset_extra["rejected_variants_rows"] = combined.get("rows")
        dataset_extra["reject_source_files_removed"] = combined.get(
            "source_files_removed", []
        )
        dataset_extra["reject_source_files_retained"] = combined.get(
            "source_files_retained", []
        )
    else:
        logger.warn(
            "rejects.enabled is false, so no rejected-variants file and no "
            "reason-by-chromosome table were written for this dataset."
        )
        partition_total = sum(
            _required_row_count(
                count, "partition_rows[%s]" % chromosome,
            )
            for chromosome, count in partition_rows.items()
        )
        if partition_total != prepared.ready_variants:
            raise PipelineError(
                "Chromosome partitions contain %d row(s), but %d row(s) were "
                "ready for chromosome processing."
                % (partition_total, prepared.ready_variants),
                results=per_chr_qc,
                chromosomes=failed,
                status="FAILED",
            )
        dataset_extra["reconciliation"] = {
            "scope": "end_to_end_parsed_rows",
            "chromosomes": all_chromosomes,
            "completed_chromosomes": completed,
            "failed_chromosomes": failed,
            "rows_read": prepared.rows_read,
            "rows_ready_for_chromosomes": prepared.ready_variants,
            "partition_rows_by_chromosome": dict(partition_rows),
            "rejected": None,
            "rows_exported": None,
            "unprocessed_failed_chromosome_rows": None,
            "balanced": None,
            "complete": None,
            "note": (
                "rejects.enabled is false, so partition coverage was verified "
                "but the end-to-end terminal buckets cannot be proven"
            ),
        }

    per_chr_qc["total_variant_infile"] = prepared.input_variants
    per_chr_qc["total_variant_read"] = prepared.rows_read
    per_chr_qc["total_variant_removed_null_coords"] = (
        prepared.removed_invalid_coordinates
    )
    per_chr_qc["total_variant_removed_unsupported_chromosomes"] = (
        prepared.removed_unsupported_chromosomes
    )
    per_chr_qc["total_variant_removed_non_standard_alleles"] = (
        prepared.removed_non_standard_alleles
    )
    per_chr_qc["total_variant_remaining_for_harmonisation"] = (
        prepared.ready_variants
    )
    per_chr_qc["total_variant_ready_snps"] = prepared.ready_snps
    per_chr_qc["total_variant_ready_indels_or_other"] = (
        prepared.ready_indels_or_other
    )
    per_chr_qc["total_variant_removed_missing_values"] = (
        prepared.removed_missing_values
    )
    per_chr_qc["total_variant_removed_duplicates"] = prepared.removed_duplicates

    dataset_extra["chromosomes"] = dict(
        (chrom, {"status": statuses[chrom], "attempts": attempts[chrom]})
        for chrom in all_chromosomes
    )
    dataset_extra["completed"] = completed
    dataset_extra["failed"] = failed
    dataset_extra["rounds"] = round_no
    dataset_extra["genome_build"] = genome_build_info
    dataset_extra["chromosome_summaries"] = summaries
    field_completeness = finalise_field_completeness_report(
        field_completeness,
        summaries,
        completed,
        expected_chromosomes=all_chromosomes,
    )
    write_field_completeness_report(
        field_completeness_path,
        field_completeness,
        delimiter=gwas2vcf_input["delimiter"],
    )
    logger.info(
        "Final field-completeness QC written to %s."
        % field_completeness_path
    )
    dataset_extra["field_completeness"] = field_completeness
    _announce(logger, "\n" + _final_field_completeness_block(field_completeness))
    dataset_extra["study_decisions"] = finalise_eaf_decision_from_chromosomes(
        dataset_extra["study_decisions"], summaries, completed,
    )
    final_frequency = dataset_extra["study_decisions"]
    logger.info(
        "Final frequency type: %s (source=%s; chromosome reference decisions=%s)."
        % (
            final_frequency["frequency_type"],
            final_frequency["eaf_is_maf_source"],
            json.dumps(final_frequency["eaf_reference_decisions"], sort_keys=True),
        )
    )

    if not failed:
        dataset_extra["status"] = "OK"
    elif completed:
        dataset_extra["status"] = "PARTIAL"
    else:
        dataset_extra["status"] = "FAILED"
    per_chr_qc["__dataset__"] = dataset_extra

    if failed:
        if progress is not None:
            progress.finish_chromosomes(failed)
        chromosome_log_pattern = configured_output_path(
            output_dir_path,
            output_layout["chromosome_log"],
            dataset_id=sample_id,
            chromosome="{N}",
        )
        message = (
            "%d of %d chromosome%s could not be completed after %d round%s: %s\n"
            "Per-chromosome logs: %s"
            % (len(failed), n_chr, "" if n_chr == 1 else "s", round_no,
               "" if round_no == 1 else "s", ", ".join(failed),
               chromosome_log_pattern)
        )
        _announce(
            logger,
            "\n" + screen_field(
                "warning", "Incomplete dataset", message,
                indent=4, label_width=20,
            ),
            marker="WARNING",
        )
        if bool(pol.get("execution.fail_dataset_on_chr_error")):
            raise PipelineError(
                message
                + "\nexecution.fail_dataset_on_chr_error is true, so this dataset is "
                  "not reported as a success on an incomplete genome. Set it to false "
                  "to merge what did complete.",
                results=per_chr_qc,
                chromosomes=failed,
                status=dataset_extra["status"],
            )
        _announce(
            logger,
            screen_field(
                "warning", "Partial merge",
                "continuing with %d of %d chromosomes because fail-on-error is disabled"
                % (len(completed), n_chr),
                indent=4, label_width=20,
            ),
            marker="WARNING",
        )
    else:
        if progress is not None:
            progress.finish_chromosomes()
        _announce(
            logger,
            "\n" + screen_line(
                "success", "All %d chromosomes completed" % n_chr, indent=4,
            ),
        )

    return per_chr_qc, grch_version


# ===============================================================
#  High-level pipeline wrapper
# ===============================================================

def _save_qc_results(
    qc_results: Dict[str, Any],
    out_file: Path,
    *,
    allowed_chromosomes,
    delimiter: str,
    null_output: str,
):
    """Write chromosome-wise harmonisation metrics as a readable TSV."""
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    table = qc_results_to_dataframe(
        data=qc_results,
        allowed_chromosomes=allowed_chromosomes,
    )
    if table is None or table.empty:
        table = _fallback_qc_frame(qc_results)
    table.to_csv(
        out_file,
        sep=delimiter,
        index=False,
        na_rep=null_output,
        # Seventeen significant digits preserve round-trippable Float64 values
        # while rendering integer-valued counts as 10 rather than 10.0.
        float_format="%.17g",
    )
    return table


def _metric_total(frame, metric):
    # type: (Any, str) -> Optional[int]
    """Sum one QC metric across every chromosome column.

    Returns ``None`` when the metric is genuinely unavailable, so a caller can
    report "unavailable" instead of printing a fabricated zero or failing while
    rendering the final summary.
    """
    if frame is None:
        return None
    try:
        if getattr(frame, "empty", False):
            return None
        if "metric" not in frame.columns:
            return None
        rows = frame.loc[frame["metric"] == metric]
        if rows.empty:
            return None
        value_columns = [c for c in rows.columns if c not in ("section", "metric")]
        if not value_columns:
            return None
        total = (
            rows[value_columns]
            .apply(pd.to_numeric, errors="coerce")
            .sum(axis=1)
            .iloc[0]
        )
    except Exception:
        return None
    if pd.isna(total):
        return None
    return int(total)


def _fallback_qc_frame(qc_results):
    # type: (Dict[str, Any]) -> Any
    """The same schema as ``qc_results_to_dataframe`` produces.

    The old fallback built a one-row wide frame with no ``metric`` column, and
    the block that consumes it selects on ``metric`` — so that path was
    guaranteed to raise.  Producing the success path's schema means the
    consumer needs no special case.
    """
    metrics = [
        "total_variant_infile",
        "total_variant_read",
        "total_variant_removed_null_coords",
        "total_variant_removed_unsupported_chromosomes",
        "total_variant_removed_non_standard_alleles",
        "total_variant_remaining_for_harmonisation",
        "total_variant_ready_snps",
        "total_variant_ready_indels_or_other",
        "total_variant_removed_missing_values",
        "total_variant_removed_duplicates",
    ]
    return pd.DataFrame(
        [
            {
                "section": "GLOBAL",
                "metric": name,
                "GLOBAL_VALUE": qc_results.get(name, pd.NA),
            }
            for name in metrics
        ]
    )


def _completed_adapter_row_expectations(dataset_extra):
    # type: (Dict[str, Any]) -> Dict[str, Optional[int]]
    """Return independently accounted export rows for completed chromosomes."""
    completed = [str(value) for value in dataset_extra.get("completed") or []]
    if not completed:
        raise PipelineError(
            "Cannot finalise GWAS-to-VCF summary audits: no completed "
            "chromosomes were recorded."
        )
    summaries = dataset_extra.get("chromosome_summaries") or {}
    expected = {}
    for chromosome in completed:
        summary = summaries.get(chromosome) or {}
        value = summary.get("rows_out")
        if value is None:
            expected[chromosome] = None
            continue
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise PipelineError(
                "Cannot finalise GWAS-to-VCF summary audits: chromosome %s has "
                "invalid independently accounted rows_out=%r."
                % (chromosome, value)
            )
        expected[chromosome] = value

    reconciled = (
        dataset_extra.get("reconciliation") or {}
    ).get("rows_exported")
    if reconciled is not None:
        if (
            isinstance(reconciled, bool)
            or not isinstance(reconciled, int)
            or reconciled < 0
        ):
            raise PipelineError(
                "Cannot finalise GWAS-to-VCF summary audits: dataset "
                "reconciliation has invalid rows_exported=%r." % reconciled
            )
        unavailable = [
            chromosome for chromosome, value in expected.items()
            if value is None
        ]
        if unavailable:
            raise PipelineError(
                "Cannot finalise GWAS-to-VCF summary audits: dataset "
                "reconciliation reports %d exported rows, but chromosome-level "
                "rows_out is unavailable for: %s."
                % (reconciled, ", ".join(unavailable))
            )
        chromosome_total = sum(expected.values())
        if chromosome_total != reconciled:
            raise PipelineError(
                "Cannot finalise GWAS-to-VCF summary audits: chromosome-level "
                "rows_out totals %d, but dataset reconciliation reports %d."
                % (chromosome_total, reconciled)
            )
    return expected


def _run_post_merge_stages(
    *,
    sample_column_dict,
    sample_id,
    grch_version,
    dataset_extra,
    qc_results,
    chromosome_qc_frame,
    default_comparison_af,
    comparison_cols,
    output_layout,
    gwas2vcf_input,
    vcf_config,
    policies,
    threads,
    executables,
    compression_executable,
    qc_configuration,
    population_frequency_config,
    bcftools_bin,
    manifest,
    manifest_path,
    logger,
):
    # type: (...) -> _PostMergeOutcome
    """Run the five post-merge stages without changing their order.

    ``manifest`` is intentionally updated after the same individual stages as
    before. If a later stage fails, the outer pipeline finalizer therefore
    retains exactly the same partial provenance and terminal status.
    """
    postgwas_qc_df2 = chromosome_qc_frame
    vcf_provenance = build_harmonisation_vcf_provenance(
        sample=sample_column_dict,
        study_decisions=dataset_extra.get("study_decisions") or {},
        resource_preflight=dataset_extra.get("resource_preflight") or {},
        completed_chromosomes=dataset_extra.get("completed") or [],
        policies=policies,
        resource_directory=sample_column_dict["resource_folder"],
        output_directory=sample_column_dict["output_folder"],
        run_manifest=str(manifest_path),
        vcf_config=vcf_config,
    )
    manifest["vcf_provenance"] = vcf_provenance

    # =====================================================
    # POST MERGE 01/05 — concatenate the per-chromosome VCFs
    # =====================================================
    with logger.step(
        1, POST_MERGE_STEP_TOTAL, "Concatenate the per-chromosome VCFs",
        "vcf_processing.concat_vcfs_by_build",
        policy_keys=["chromosome.allowed_after_split", "vcf.concat_max_workers",
                     "vcf.concat_max_attempts", "vcf.on_merge_failure",
                     "vcf.min_valid_size_bytes"],
    ) as ctx:
        concated_vcf_files = concat_vcfs_by_build(
            output_dir=sample_column_dict["output_folder"],
            gwas_outputname=sample_id,
            grch_version=grch_version,
            expected_chromosomes=dataset_extra.get("completed", []),
            threads=threads,
            policies=policies,
            output_layout=output_layout,
            executables=executables,
            vcf_config=vcf_config,
            provenance=vcf_provenance,
            logger=logger,
        )
        manifest["merged_vcfs"] = concated_vcf_files
        ctx.extra["merged_vcfs"] = concated_vcf_files

    # =====================================================
    # POST MERGE 02/05 — population similarity on the raw input-build VCF
    # =====================================================
    with logger.step(
        2, POST_MERGE_STEP_TOTAL,
        "Compare study and reference-population frequencies",
        "population_frequency.run_population_frequency_qc",
    ) as ctx:
        input_build_vcf = concated_vcf_files.get(grch_version.lower())
        incomplete_input_build = grch_version in set(
            concated_vcf_files.get("required_merge_failures") or []
        )
        if not bool(population_frequency_config["enabled"]):
            population_frequency_result = {
                "status": "disabled",
                "closest_population": None,
                "warnings": [],
                "external_file_checks": [],
            }
            ctx.skip("Population-frequency similarity is disabled in configuration.")
        elif input_build_vcf is None or incomplete_input_build:
            population_frequency_result = {
                "status": "skipped_incomplete_input_build_vcf",
                "closest_population": None,
                "warnings": [],
                "external_file_checks": [],
                "decision_reason": (
                    "the required unfiltered input-build merged VCF is missing "
                    "or incomplete"
                ),
            }
            ctx.warn(population_frequency_result["decision_reason"])
        else:
            try:
                population_frequency_result = run_population_frequency_qc(
                    vcf_path=input_build_vcf,
                    output_directory=sample_column_dict["output_folder"],
                    dataset_id=sample_id,
                    genome_build=grch_version,
                    settings=population_frequency_config,
                    output_layout=output_layout,
                    table_delimiter=vcf_config["table_delimiter"],
                    table_null_values=vcf_config["table_null_values"],
                    temporary_table_suffix=vcf_config["temporary_table_suffix"],
                    io_buffer_bytes=vcf_config["io_buffer_bytes"],
                    bcftools_bin=bcftools_bin,
                    selected_population=comparison_cols,
                    external_files={
                        "external EAF": sample_column_dict.get(
                            "provided_external_eaf_file"
                        ),
                        "external INFO": sample_column_dict.get(
                            "provided_external_info_file"
                        ),
                    },
                    logger=logger,
                )
                if (
                    population_frequency_result.get("status")
                    == "frequency_inversion_suspected"
                ):
                    ctx.extra["population_frequency_qc"] = (
                        population_frequency_result
                    )
                    ctx.decide(
                        "Study frequency orientation",
                        population_frequency_result.get("decision_reason"),
                        "failed; the merged AF is not safe for downstream use",
                    )
                    raise PipelineError(
                        "Population-frequency QC detected a likely non-effect-allele "
                        "frequency column: %s. PostGWAS will not treat this run as "
                        "successful. Supply the correct effect-allele-frequency "
                        "column and rerun. Evidence was written to %s."
                        % (
                            population_frequency_result.get("decision_reason"),
                            population_frequency_result.get(
                                "report", "the population-frequency QC report"
                            ),
                        )
                    )
            except PopulationFrequencyQCError as exc:
                if population_frequency_config["on_error"] == "fail":
                    raise PipelineError(
                        "Population-frequency QC failed: %s" % exc
                    ) from exc
                population_frequency_result = {
                    "status": "failed_warning",
                    "closest_population": None,
                    "warnings": [],
                    "external_file_checks": [],
                    "decision_reason": str(exc),
                }
                ctx.warn(
                    "Population-frequency QC could not be completed: %s. "
                    "Processing will continue because on_error='warn'." % exc
                )
        ctx.extra["population_frequency_qc"] = population_frequency_result
        population_records = population_frequency_result.get("total_records")
        if population_records is not None:
            ctx.rows_in = population_records
            ctx.set_rows(population_records, removed=0)
    manifest["population_frequency_qc"] = population_frequency_result
    _announce(
        logger,
        "\n" + _population_frequency_block(population_frequency_result),
    )
    for warning in population_frequency_result.get("warnings") or []:
        logger.warn(warning)

    # =====================================================
    # POST MERGE 03/05 — clean the per-chromosome intermediates
    # =====================================================
    expected_adapter_rows = _completed_adapter_row_expectations(dataset_extra)
    with logger.step(
        3, POST_MERGE_STEP_TOTAL, "Merge the side files and clean up",
        "cleanup.finalise_harmonisation_outputs",
    ) as ctx:
        cleaned = finalise_harmonisation_outputs(
            output_dir=sample_column_dict["output_folder"],
            gwas_outputname=sample_id,
            output_layout=output_layout,
            summary_delimiter=gwas2vcf_input["delimiter"],
            expected_chromosomes=list(expected_adapter_rows),
            expected_rows_by_chromosome=expected_adapter_rows,
            threads=threads,
            compression_executable=compression_executable,
            policies=policies,
            logger=logger,
        )
        ctx.extra["cleaned"] = {
            "moved": len(cleaned.get("moved", [])),
            "removed": len(cleaned.get("removed", [])),
            "adapter_input_rows": cleaned.get("total_adapter_input_rows"),
        }

    observed_adapter_rows = cleaned.get("adapter_input_rows_by_chromosome")
    total_variants_in_gwas2vcf_input = cleaned.get("total_adapter_input_rows")
    if not isinstance(observed_adapter_rows, dict):
        raise PipelineError(
            "Validated GWAS-to-VCF cleanup did not return per-chromosome "
            "adapter-input row counts; the final QC report cannot be trusted."
        )
    if set(observed_adapter_rows) != set(expected_adapter_rows):
        raise PipelineError(
            "Validated GWAS-to-VCF cleanup returned chromosome counts for %s; "
            "expected exactly %s."
            % (
                ", ".join(sorted(str(key) for key in observed_adapter_rows)),
                ", ".join(expected_adapter_rows),
            )
        )
    for chromosome, observed_rows in observed_adapter_rows.items():
        if (
            isinstance(observed_rows, bool)
            or not isinstance(observed_rows, int)
            or observed_rows < 0
        ):
            raise PipelineError(
                "Validated GWAS-to-VCF cleanup returned invalid adapter-input "
                "row count %r for chromosome %s; expected a non-negative "
                "integer." % (observed_rows, chromosome)
            )
    if (
        isinstance(total_variants_in_gwas2vcf_input, bool)
        or not isinstance(total_variants_in_gwas2vcf_input, int)
        or total_variants_in_gwas2vcf_input < 0
        or total_variants_in_gwas2vcf_input
        != sum(observed_adapter_rows.values())
    ):
        raise PipelineError(
            "Validated GWAS-to-VCF cleanup returned inconsistent total and "
            "per-chromosome adapter-input row counts."
        )
    for chromosome, expected_rows in expected_adapter_rows.items():
        observed_rows = observed_adapter_rows.get(chromosome)
        if expected_rows is not None and observed_rows != expected_rows:
            raise PipelineError(
                "GWAS-to-VCF adapter summary reports %r input rows for "
                "chromosome %s, but independent reconciliation reports %d."
                % (observed_rows, chromosome, expected_rows)
            )
    manifest["gwas2vcf_summary_audit"] = {
        "path": cleaned.get("summary"),
        "rows_by_chromosome": dict(observed_adapter_rows),
        "total_adapter_input_rows": total_variants_in_gwas2vcf_input,
    }

    outdir = Path(sample_column_dict["output_folder"])

    # =====================================================
    # Prepare the pre-VCF counts used by the combined QC report
    # =====================================================
    qc_input_build = grch_version

    raw_vcf_path = configured_output_path(
        outdir,
        output_layout["merged_build_vcf"],
        error_type=PipelineError,
        dataset_id=sample_id,
        build=qc_input_build,
    )

    missing_eaf_count = _metric_total(postgwas_qc_df2, "missing_eaf_count")
    invalid_beta_se_count = _metric_total(
        postgwas_qc_df2, "variants_removed_invalid_beta_se"
    )
    if missing_eaf_count is None:
        logger.warn(
            "The QC table has no 'missing_eaf_count' metric, so the number of "
            "variants with no frequency is reported as unavailable rather than "
            "as a fabricated zero."
        )
    if invalid_beta_se_count is None:
        logger.warn(
            "The QC table has no 'variants_removed_invalid_beta_se' metric, so "
            "that count is reported as unavailable."
        )

    dataset_reconciliation = dataset_extra.get("reconciliation") or {}
    pre_vcf_summary = {
        "total_variant_infile": qc_results.get("total_variant_infile"),
        "total_variant_read": qc_results.get("total_variant_read"),
        "total_variant_in_vcf_input": total_variants_in_gwas2vcf_input,
        "total_variant_removed_missing_values": qc_results.get(
            "total_variant_removed_missing_values"
        ),
        "total_variant_removed_duplicates": qc_results.get(
            "total_variant_removed_duplicates"
        ),
        "total_variant_removed_null_coords": qc_results.get(
            "total_variant_removed_null_coords"
        ),
        "total_variant_removed_unsupported_chromosomes": qc_results.get(
            "total_variant_removed_unsupported_chromosomes"
        ),
        "total_variant_removed_non_standard_alleles": qc_results.get(
            "total_variant_removed_non_standard_alleles"
        ),
        "total_variant_remaining_for_harmonisation": qc_results.get(
            "total_variant_remaining_for_harmonisation"
        ),
        "total_variant_ready_snps": qc_results.get(
            "total_variant_ready_snps"
        ),
        "total_variant_ready_indels_or_other": qc_results.get(
            "total_variant_ready_indels_or_other"
        ),
        "total_variant_removed_palindromic_ambiguous": (
            dataset_extra.get("reject_counts") or {}
        ).get("palindromic_ambiguous"),
        "total_variant_removed_palindromic_frequency_conflict": (
            dataset_extra.get("reject_counts") or {}
        ).get("palindromic_frequency_conflict"),
        "total_variant_removed_palindromic_orientation_unavailable": (
            dataset_extra.get("reject_counts") or {}
        ).get("palindromic_orientation_unavailable"),
        "total_variant_removed_palindromic_frequency_discordant": (
            dataset_extra.get("reject_counts") or {}
        ).get("palindromic_frequency_discordant"),
        "total_variant_removed_reference_unmatched": (
            dataset_extra.get("reject_counts") or {}
        ).get("reference_unmatched"),
        "total_variant_removed_reference_ambiguous": (
            dataset_extra.get("reject_counts") or {}
        ).get("reference_ambiguous"),
        "total_variant_removed_chromosome_harmonisation": (
            dataset_reconciliation.get("chromosome_stage_rejected")
        ),
        "total_variant_passed_chromosome_harmonisation": (
            dataset_reconciliation.get("rows_exported")
        ),
        "total_variant_rejected_end_to_end": dataset_reconciliation.get(
            "rejected"
        ),
        "total_variant_unprocessed_failed_chromosomes": (
            dataset_reconciliation.get(
                "unprocessed_failed_chromosome_rows"
            )
        ),
        "total_variant_accounted_end_to_end": dataset_reconciliation.get(
            "accounted_rows"
        ),
        "row_accounting_balanced": dataset_reconciliation.get("balanced"),
        "row_accounting_complete": dataset_reconciliation.get("complete"),
        "total_variant_with_missing_eaf": missing_eaf_count,
        "total_variant_with_invalid_beta_se": invalid_beta_se_count,
    }
    manifest["pre_vcf_summary"] = pre_vcf_summary

    # =====================================================
    # POST MERGE 04/05 — one VCF extraction, raw QC and virtual filters
    # =====================================================
    with logger.step(
        4, POST_MERGE_STEP_TOTAL, "Assess raw VCF quality",
        "qc_summary.service.run_qc_assessment",
    ) as ctx:
        ctx.info(
            "Using schema-validated modules.qc_summary rules for the "
            "input-build merged VCF and reference-frequency tag INFO/%s; "
            "the QC build is inferred from the VCF header."
            % comparison_cols
        )
        if (
            (dataset_extra.get("study_decisions") or {}).get(
                "info_score_type"
            ) == "mach_rsq"
            and float(qc_configuration.rules.info_max)
            < float(policies.get("info.mach_rsq_max"))
        ):
            ctx.warn(
                "The resolved imputation-quality type is mach_rsq, whose "
                "harmonisation maximum is {:g}, while the independently "
                "configured modules.qc_summary.rules.info_max is {:g}. Values "
                "above the QC maximum remain in the delivered raw VCF but are "
                "counted as failing the virtual QC assessment."
                .format(
                    float(policies.get("info.mach_rsq_max")),
                    float(qc_configuration.rules.info_max),
                )
            )
        assessment = run_qc_assessment(
            vcf_path=raw_vcf_path,
            output_directory=outdir,
            dataset_id=sample_id,
            external_af_name=comparison_cols,
            configuration=qc_configuration,
            bcftools_bin=bcftools_bin,
            genome_build_header=vcf_config["genome_build_header"],
            supported_genome_builds=tuple(vcf_config["target_builds"]),
            threads=threads,
            provenance_headers=dict(vcf_config["provenance"]["headers"]),
            logger=logger,
        )
        manifest["qc_genome_build"] = assessment["genome_build"]
        low_neff_summary = {
            "reference_quantile": assessment["sample_size_reference_quantile"],
            "reference_value": assessment["raw"][
                "effective_sample_size_reference_quantile_value"
            ],
            "minimum_fraction_of_reference": assessment[
                "sample_size_minimum_fraction_of_reference"
            ],
            "minimum_threshold": assessment["raw"][
                "effective_sample_size_minimum_threshold"
            ],
            "raw_below_threshold": assessment["raw"][
                "effective_sample_size_below_minimum_threshold"
            ],
            "raw_below_threshold_fraction": assessment["raw"][
                "effective_sample_size_below_minimum_threshold_fraction"
            ],
            "qc_passed_below_threshold": assessment["qc_passed"][
                "effective_sample_size_below_minimum_threshold"
            ],
            "qc_passed_below_threshold_fraction": assessment["qc_passed"][
                "effective_sample_size_below_minimum_threshold_fraction"
            ],
        }
        ctx.rows_in = assessment["raw"]["num_records"]
        ctx.set_rows(assessment["qc_passed"]["num_records"])
        ctx.extra["qc_assessment"] = {
            "raw_variants": assessment["raw"]["num_records"],
            "qc_passed_variants": assessment["qc_passed"]["num_records"],
            "excluded": assessment["excluded_total"],
            "low_neff": low_neff_summary,
            "reports": assessment["reports"],
        }
    manifest["qc_assessment"] = {
        "definition": assessment["definition"],
        "raw_variants": assessment["raw"]["num_records"],
        "qc_passed_variants": assessment["qc_passed"]["num_records"],
        "excluded": assessment["excluded_total"],
        "retained_fraction": assessment["retained_fraction"],
        "low_neff": low_neff_summary,
        "reports": assessment["reports"],
    }

    # =====================================================
    # POST MERGE 05/05 — final QC report, cleanup and logs
    # =====================================================
    keep_gwas2vcf_intermediate = bool(
        policies.get("vcf.keep_gwas2vcf_intermediate")
    )
    with logger.step(
        5, POST_MERGE_STEP_TOTAL, "Final QC report and logs",
        "service.run_harmonisation_pipeline",
        policy_keys=["vcf.keep_gwas2vcf_intermediate"],
    ) as ctx:
        ctx.info(
            "QC assessment reports: %s"
            % ", ".join(assessment["reports"].values())
        )
        ctx.info(
            "Raw GWAS-to-VCF intermediate policy after successful finalization: %s."
            % ("retain" if keep_gwas2vcf_intermediate else "delete")
        )

        _announce(
            logger,
            "\n" + "\n".join(
                harmonisation_qc_summary_lines(pre_vcf_summary, assessment)
            ) + "\n",
        )

        validated_outputs = _resolved_harmonisation_outputs(
            outdir, sample_id, grch_version, output_layout, vcf_config,
        )
        raw_gwas2vcf_intermediate = validated_outputs["gwas2vcf"]
        primary_outputs = dict(validated_outputs)
        if not keep_gwas2vcf_intermediate:
            primary_outputs.pop("gwas2vcf")

        status = _combined_dataset_status(
            dataset_extra.get("status", "OK"),
            concated_vcf_files.get("merge_status", "OK"),
        )
        manifest["status"] = status
        _announce(
            logger,
            "\n" + "\n".join(harmonisation_qc_takeaway_lines(
                pre_vcf_summary,
                assessment,
                genome_build_info=dataset_extra.get("genome_build") or {},
                study_decisions=dataset_extra.get("study_decisions") or {},
                population_frequency=population_frequency_result,
                dataset_summary=dataset_extra,
                merge_summary=concated_vcf_files,
                primary_outputs=primary_outputs,
                final_status=status,
                reference_source=default_comparison_af,
                reference_population=comparison_cols,
            )) + "\n",
        )
        combined_log = _write_combined_log(
            outdir,
            sample_id,
            dataset_extra.get("chromosomes") or {},
            output_layout,
        )
        manifest["combined_log"] = combined_log

        removed_gwas2vcf_artifacts = []
        if status == "OK" and not keep_gwas2vcf_intermediate:
            removed_gwas2vcf_artifacts = remove_merged_gwas2vcf_intermediate(
                raw_gwas2vcf_intermediate,
                output_directory=outdir,
                logger=logger,
            )
        retained_gwas2vcf_intermediate = Path(
            raw_gwas2vcf_intermediate
        ).is_file()
        manifest["gwas2vcf_intermediate"] = {
            "path": raw_gwas2vcf_intermediate,
            "keep_requested": keep_gwas2vcf_intermediate,
            "retained": retained_gwas2vcf_intermediate,
            "removed_artifacts": removed_gwas2vcf_artifacts,
            "retention_reason": (
                "explicit_keep_policy"
                if keep_gwas2vcf_intermediate
                else "non_ok_dataset_retained_for_diagnosis"
                if status != "OK"
                else "successful_run_default_cleanup"
            ),
        }
        _announce(
            logger,
            "\n" + "\n".join([
                screen_line(
                    "success" if status == "OK" else "warning",
                    "Harmonisation and GWAS-to-VCF finished",
                    indent=4,
                ),
                screen_field(
                    "info", "Dataset", sample_id,
                    indent=6, label_width=16,
                ),
                screen_field(
                    "success" if status == "OK" else "warning",
                    "Status", status,
                    indent=6, label_width=16,
                ),
            ]) + "\n",
        )

    return _PostMergeOutcome(
        primary_outputs=primary_outputs,
        assessment=assessment,
        population_frequency=population_frequency_result,
        concatenated_vcfs=concated_vcf_files,
        pre_vcf_summary=pre_vcf_summary,
        status=status,
        combined_log=combined_log,
    )


def run_harmonisation_pipeline(
    sample_column_dict: Dict[str, Any],
    default_cfg: Dict[str, Any],
    threads: int,
):
    """Full harmonisation + gwas2vcf + merge + cleanup pipeline for one dataset.

    Restructured to plan Part 2: the configuration and the file header are
    checked *before* the ten-million-row parse, so a mistyped column name costs
    a second rather than a full read.  Nothing here calls ``sys.exit``; a bad
    configuration raises :class:`ConfigError` and an incomplete run raises
    :class:`PipelineError`, so ``cli.py`` can carry on with the next dataset.
    Once the run manifest has been created, every ordinary exception records a
    terminal ``FAILED`` or ``PARTIAL`` state before it is re-raised; a user
    interruption records ``INTERRUPTED``.
    """
    started_at = time.time()

    # -----------------------------------------------------
    # Load DEFAULT config
    # -----------------------------------------------------
    try:
        default_cfg_obj = dict(default_cfg)
        if not default_cfg_obj:
            raise ValueError("engine configuration is empty")
    except Exception as e:
        raise ConfigError(
            [],
            message=(
                "Invalid harmonisation engine configuration '%s'.\n   Reason: %s"
                % (default_cfg, e)
            ),
        )

    # Extract all defaults safely
    default_eaf_reference_source = default_cfg_obj[
        "default_eaf_reference_source"
    ]
    default_eaf_reference_column = default_cfg_obj[
        "default_eaf_reference_column"
    ]
    default_comparison_af = default_cfg_obj["default_comparison_af_file"]
    comparison_cols = default_cfg_obj["default_comparison_af_column"]
    dbsnp = default_cfg_obj["default_dbsnp"]
    external_eaf_colmap = default_cfg_obj["external_eaf_colmap"]
    default_eaf_colmap = default_cfg_obj["default_eaf_colmap"]
    external_info_colmap = default_cfg_obj["external_info_colmap"]
    build_check_colmap = default_cfg_obj["build_check_colmap"]
    resource_layout = default_cfg_obj["resource_layout"]
    output_layout = default_cfg_obj["output_layout"]
    gwas2vcf_input = default_cfg_obj["gwas2vcf_input"]
    executables = default_cfg_obj["executables"]
    compression_executable = default_cfg_obj["compression_executable"]
    vcf_config = default_cfg_obj["vcf_processing"]
    qc_configuration = QCSummaryConfig.model_validate(
        default_cfg_obj["qc_summary"]
    )
    population_frequency_config = default_cfg_obj["population_frequency_qc"]
    external_reference_staging = default_cfg_obj[
        "external_reference_staging"
    ]
    terminal_progress = default_cfg_obj.get("terminal_progress")
    if not isinstance(terminal_progress, dict):
        raise ConfigError(
            [],
            message=(
                "Resolved harmonisation configuration is missing the canonical "
                "terminal-progress settings. Rebuild it through the PostGWAS "
                "configuration loader."
            ),
        )
    progress_enabled = terminal_progress.get("enabled")
    progress_label_width = terminal_progress.get("outcome_label_width")
    if not isinstance(progress_enabled, bool):
        raise ConfigError(
            [],
            message=(
                "Resolved logging.show_progress must be true or false; found %r."
                % progress_enabled
            ),
        )
    if (
        isinstance(progress_label_width, bool)
        or not isinstance(progress_label_width, int)
        or progress_label_width < 1
    ):
        raise ConfigError(
            [],
            message=(
                "Resolved logging.terminal_label_width must be a positive "
                "integer; found %r." % progress_label_width
            ),
        )
    executables_prevalidated = bool(
        default_cfg_obj.get("executables_prevalidated", False)
    )

    # Every value the run uses, resolved once, reported everywhere.
    policies = _resolve_policies(default_cfg_obj, sample_column_dict)

    # The sample-sheet boundary has already resolved the canonical dataset
    # harmonisation root. All configured result paths are relative to it.
    output_folder = Path(
        sample_column_dict["output_folder"]
    ).expanduser().resolve()
    gwas2vcf_script = str(
        Path(__file__).parent / "adapters" / "gwas2vcf" / "main.py"
    )

    # Save the normalized root back for every chromosome and post-merge caller.
    sample_column_dict["output_folder"] = str(output_folder)
    # Clean whitespace from keys AND values
    sample_column_dict = {
        k.strip(): (v.strip() if isinstance(v, str) else v)
        for k, v in sample_column_dict.items()
    }
    try:
        output_folder.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        raise PipelineError(
            "Failed to create the output directory.\n   Path: %s\n   Reason: %s"
            % (output_folder, e)
        )

    sample_id = _safe_glob_name(
        sample_column_dict.get("gwas_outputname"), "gwas_outputname"
    )
    log_dir = configured_output_path(
        output_folder, output_layout["logs_directory"], error_type=PipelineError,
    )
    progress = _HarmonisationProgress(
        enabled=progress_enabled,
        outcome_label_width=progress_label_width,
    )
    logger = PipelineLogger(
        sample_id=sample_id,
        scope="dataset",
        log_dir=str(log_dir),
        policies=policies,
        level=None,
        screen_level=None,
        log_path=str(configured_output_path(
            output_folder,
            output_layout["dataset_log"],
            error_type=PipelineError,
            dataset_id=sample_id,
        )),
        stage_progress=progress.dataset_stages,
    )
    manifest_path = configured_output_path(
        output_folder,
        output_layout["run_manifest"],
        error_type=PipelineError,
        dataset_id=sample_id,
    )
    manifest = None  # type: Optional[Dict[str, Any]]

    try:
        logger.blank()
        logger.info(
            "Harmonising %s from %s."
            % (sample_id, sample_column_dict.get("sumstat_file"))
        )
        if executables_prevalidated:
            logger.info(
                "External executables and the bcftools liftover plugin were validated "
                "during run preflight and will be reused by every chromosome."
            )
        else:
            executables = require_binaries(
                executables,
                plugins=(str(vcf_config["liftover_plugin"]),),
                logger=logger,
            )
        bcftools_bin = executables["bcftools"]
        _log_resolved_policies(logger, policies)

        manifest = {
            "sample_id": sample_id,
            "sumstat_file": sample_column_dict.get("sumstat_file"),
            "output_folder": str(output_folder),
            "defaults_file": str(default_cfg),
            "started": time.strftime("%Y-%m-%d %H:%M:%S"),
            "versions": {
                "python": platform.python_version(),
                "polars": getattr(pl, "__version__", "unknown"),
                "pandas": getattr(pd, "__version__", "unknown"),
                "platform": platform.platform(),
            },
            "config": dict(sample_column_dict),
            "policies": policies.resolved_dict(),
            "policies_changed_from_default": dict(
                (key, {"default": was, "value": now})
                for key, (was, now) in policies.changed_from_default().items()
            ),
            "status": "running",
        }
        _write_run_manifest(manifest_path, manifest)
        _announce(
            logger,
            "\n" + "\n".join([
                screen_line("info", "Resolved run settings", indent=4),
                screen_field(
                    "info", "Manifest", manifest_path,
                    indent=8, label_width=18,
                ),
                screen_field(
                    "count", "Settings", "%d total; %d changed from default" % (
                        len(policies.keys()), len(policies.changed_from_default()),
                    ), indent=8, label_width=18,
                ),
            ]),
        )
        _announce(
            logger,
            "\n" + "\n".join([
                screen_line("genetic", "Allele-frequency references", indent=4),
                screen_field(
                    "genetic", "VCF annotation / QC", "%s · %s" % (
                        default_comparison_af, comparison_cols,
                    ), indent=8, label_width=23,
                ),
                screen_field(
                    "genetic", "Strand / MAF–EAF", "%s · %s" % (
                        default_eaf_reference_source,
                        default_eaf_reference_column,
                    ), indent=8, label_width=23,
                ),
            ]),
        )

        # =====================================================
        # 01/08  validate_config    — no file is read
        # =====================================================
        with logger.step(
            1, DATASET_STEP_TOTAL, "Check the configuration", "validator.validate_config",
            policy_keys=[
            ],
        ) as ctx:
            ok, problems = validate_config(
                _validation_view(sample_column_dict, default_cfg_obj),
                policies=policies,
                logger=logger,
                external_reference_delimiters={
                    "eaffile": external_eaf_colmap["delimiter"],
                    "infofile": external_info_colmap["delimiter"],
                },
            )
            if not ok:
                report = format_problems(
                    problems,
                    title="Configuration check for %s" % sample_id,
                    footer="Nothing was read from the input file. Fix the config and re-run.",
                )
                ctx.failure_hint = (
                    "The configuration for this dataset is not usable as written."
                )
                raise ConfigError(problems, message=report)
            ctx.info("The configuration is complete and every named file exists.")

        # ---------------------------------------------------------
        # Resource folder and build check files
        # ---------------------------------------------------------
        resource_folder = Path(sample_column_dict["resource_folder"])
        if not resource_folder.exists():
            raise ConfigError(
                [], message="Resource folder does not exist:\n   %s" % resource_folder
            )
        if not resource_folder.is_dir():
            raise ConfigError(
                [], message="Resource folder is not a directory:\n   %s" % resource_folder
            )
        if not any(resource_folder.iterdir()):
            raise ConfigError(
                [], message="Resource folder is empty:\n   %s" % resource_folder
            )

        file_validator = validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
        )
        errors = []
        build_reference_files = {
            build: configured_output_path(
                sample_column_dict["resource_folder"],
                resource_layout["build_check"],
                error_type=PipelineError,
                build=build,
            )
            for build in vcf_config["target_builds"]
        }
        for path in [
            sample_column_dict["sumstat_file"],
            *build_reference_files.values(),
        ]:
            try:
                file_validator(path)
            except argparse.ArgumentTypeError as e:
                errors.append("File invalid: %s\n   → %s" % (path, e))
        if errors:
            raise ConfigError(
                [],
                message="One or more required files are invalid:\n" + "\n".join(errors),
            )

        # =====================================================
        # 02/08  validate_header    — reads one line
        # =====================================================
        with logger.step(
            2, DATASET_STEP_TOTAL, "Check the file header", "validator.validate_header",
            policy_keys=[
                "input.delimiter",
                "validation.declaration_mismatch_action",
                "validation.fuzzy_match_cutoff",
                "validation.on_ambiguous_column_mapping",
                "validation.on_duplicate_header",
            ],
        ) as ctx:
            ok, problems = validate_header(
                _validation_view(sample_column_dict, default_cfg_obj),
                sample_column_dict["sumstat_file"],
                policies=policies,
                logger=logger,
            )
            if not ok:
                report = format_problems(
                    problems,
                    title="Header check for %s" % sample_id,
                    footer="Nothing was read beyond the header. Fix the config and re-run.",
                )
                ctx.failure_hint = (
                    "The configuration names columns that are not in the input file."
                )
                raise ConfigError(problems, message=report)
            ctx.info("Every configured column name is present in the file header.")

        # --------------------------------------------------------
        # Truncation check, before a long run
        # ---------------------------------------------------------
        _truncation_warnings, input_line_count = inspect_summary_statistics_file(
            file_path=sample_column_dict["sumstat_file"],
            io_buffer_bytes=vcf_config["io_buffer_bytes"],
            policies=policies,
            logger=logger,
        )

        # ---------------------------------------------------------
        # STEPS 03-08 + the per-chromosome fan-out
        # ---------------------------------------------------------
        dataset_extra = {}  # type: Dict[str, Any]
        try:
            qc_results, grch_version = harmonise_chromosomes(
                sumstat_file=sample_column_dict["sumstat_file"],
                sample_column_dict=sample_column_dict,
                output_dir=sample_column_dict["output_folder"],
                resource_folder=sample_column_dict["resource_folder"],
                default_eaf_reference_source=default_eaf_reference_source,
                default_comparison_af_file=default_comparison_af,
                dbsnp=dbsnp,
                resource_layout=resource_layout,
                output_layout=output_layout,
                gwas2vcf_input=gwas2vcf_input,
                executables=executables,
                vcf_config=vcf_config,
                max_workers=threads,
                user_eaf_file=sample_column_dict["eaffile"],
                user_eaf_column=sample_column_dict["eafcolumn"],
                user_info_file=sample_column_dict["infofile"],
                user_info_column=sample_column_dict["infocolumn"],
                default_eaf_reference_column=default_eaf_reference_column,
                default_comparison_af_column=comparison_cols,
                build_reference_files=build_reference_files,
                build_check_colmap=build_check_colmap,
                default_eaf_colmap=default_eaf_colmap,
                gwas2vcf_main_script_path=gwas2vcf_script,
                external_eaf_colmap=external_eaf_colmap,
                external_info_colmap=external_info_colmap,
                external_reference_staging=external_reference_staging,
                policies=policies,
                logger=logger,
                input_line_count=input_line_count,
                progress=progress,
            )
        except PipelineError as exc:
            progress.fail_active_chromosomes()
            # The chromosomes that did complete still have QC worth keeping.
            if exc.results:
                partial = dict(exc.results)
                manifest["dataset"] = partial.pop("__dataset__", {})
                _save_qc_results(
                    qc_results=partial,
                    out_file=configured_output_path(
                        output_folder,
                        output_layout["qc_summary"],
                        error_type=PipelineError,
                        dataset_id=sample_id,
                    ),
                    allowed_chromosomes=policies.get("chromosome.allowed"),
                    delimiter=vcf_config["table_delimiter"],
                    null_output=vcf_config["table_null_output"],
                )
            manifest["status"] = getattr(exc, "status", "FAILED")
            manifest["failed_chromosomes"] = getattr(exc, "chromosomes", [])
            raise
        except BaseException:
            progress.fail_active_chromosomes()
            raise

        dataset_extra = qc_results.pop("__dataset__", {})
        manifest["dataset"] = dataset_extra

        qc_summary_path = configured_output_path(
            output_folder,
            output_layout["qc_summary"],
            error_type=PipelineError,
            dataset_id=sample_id,
        )
        postgwas_qc_df2 = _save_qc_results(
            qc_results=qc_results,
            out_file=qc_summary_path,
            allowed_chromosomes=policies.get("chromosome.allowed"),
            delimiter=vcf_config["table_delimiter"],
            null_output=vcf_config["table_null_output"],
        )

        logger.set_stage_progress(progress.start_post_merge())
        post_merge = _run_post_merge_stages(
            sample_column_dict=sample_column_dict,
            sample_id=sample_id,
            grch_version=grch_version,
            dataset_extra=dataset_extra,
            qc_results=qc_results,
            chromosome_qc_frame=postgwas_qc_df2,
            default_comparison_af=default_comparison_af,
            comparison_cols=comparison_cols,
            output_layout=output_layout,
            gwas2vcf_input=gwas2vcf_input,
            vcf_config=vcf_config,
            policies=policies,
            threads=threads,
            executables=executables,
            compression_executable=compression_executable,
            qc_configuration=qc_configuration,
            population_frequency_config=population_frequency_config,
            bcftools_bin=bcftools_bin,
            manifest=manifest,
            manifest_path=manifest_path,
            logger=logger,
        )
        primary_outputs = post_merge.primary_outputs
        assessment = post_merge.assessment
        population_frequency_result = post_merge.population_frequency
        status = post_merge.status
        manifest["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
        manifest["completed_at"] = datetime.now().astimezone().isoformat(
            timespec="seconds"
        )
        manifest["elapsed_seconds"] = round(time.time() - started_at, 2)
        _write_run_manifest(manifest_path, manifest)

        result = dict(primary_outputs)
        result.update(
            qc_assessment=assessment["reports"]["summary"],
            qc_filter_rules=assessment["reports"]["rules"],
            qc_assessment_json=assessment["reports"]["json"],
            population_frequency_qc=population_frequency_result.get("report"),
            status=status,
            manifest=str(manifest_path),
        )
        return result
    except KeyboardInterrupt as exc:
        _finalise_failed_run_manifest(
            manifest_path, manifest, exc, started_at, logger,
            status="INTERRUPTED",
        )
        raise
    except Exception as exc:
        _finalise_failed_run_manifest(
            manifest_path, manifest, exc, started_at, logger,
            status="FAILED",
        )
        raise
    finally:
        logger.close()
        progress.close()


def _write_combined_log(
    output_dir, sample_id, chromosome_statuses, output_layout,
):
    # type: (Any, str, Dict[str, Any]) -> Optional[str]
    """``logs/{sample}_combined.log`` — the dataset log then every chromosome log.

    Rewritten from scratch each time, unlike the per-chromosome logs, which are
    append-only so that every retry attempt is preserved.
    """
    dataset_log = configured_output_path(
        output_dir, output_layout["dataset_log"], dataset_id=sample_id,
    )
    combined = configured_output_path(
        output_dir, output_layout["combined_log"], dataset_id=sample_id,
    )
    order = sorted(chromosome_statuses, key=_chromosome_sort_key)
    try:
        with combined.open("w", encoding="utf-8") as out:
            if dataset_log.is_file():
                out.write("===== %s =====\n" % dataset_log.name)
                out.write(dataset_log.read_text(encoding="utf-8", errors="replace"))
                out.write("\n")
            for chrom in order:
                path = configured_output_path(
                    output_dir,
                    output_layout["chromosome_log"],
                    dataset_id=sample_id,
                    chromosome=chrom,
                )
                if not path.is_file():
                    continue
                out.write("===== %s =====\n" % path.name)
                out.write(path.read_text(encoding="utf-8", errors="replace"))
                out.write("\n")
    except OSError:
        return None
    return str(combined)
