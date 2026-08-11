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
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from multiprocessing import get_context

import polars as pl
import pandas as pd


# ----------------------------------------------------------------------
# PostGWAS imports
# ----------------------------------------------------------------------
from postgwas.core.execution.runtime import validate_path

from postgwas.modules.harmonisation.policies import (
    load_policies,
    default_policies,
    PolicyError,
)
from postgwas.core.pipeline_logging import (
    PipelineLogger,
    print_chromosome_summary,
)
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.paths import configured_output_matches, configured_output_path
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    RejectOutputError,
    concat_reject_files,
    reconcile,
    REASONS,
    REASON_STEPS,
)
from postgwas.modules.harmonisation.summary_statistics_io import (
    read_summary_statistics,
    inspect_summary_statistics_file,
)
from postgwas.modules.harmonisation.resource_paths import resolve_resource_file
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
from postgwas.modules.harmonisation.shared.variant_columns import canonical_variant_schema
from postgwas.modules.harmonisation.allele_frequency import harmonise_allele_frequency
from postgwas.modules.harmonisation.sample_size import (
    harmonise_sample_sizes,
    prepare_missing_sample_sizes,
)
from postgwas.modules.harmonisation.effect_type import harmonise_effect_estimates
from postgwas.modules.harmonisation.effect_from_z import derive_effect_and_standard_error_from_z
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
from postgwas.modules.harmonisation.imputation_quality import harmonise_imputation_quality
from postgwas.modules.harmonisation.variant_identifiers import harmonise_variant_identifiers
from postgwas.modules.harmonisation.vcf_processing import (
    annotate_and_liftover_vcf,
    concat_vcfs_by_build,
    require_binaries,
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
    remove_partial_chromosome_outputs,
)
from postgwas.modules.harmonisation.resource_preflight import (
    ResourcePreflightError,
    recheck_preflighted_resource_map,
    validate_harmonisation_resource_maps,
)
from postgwas.modules.harmonisation.qc_results import qc_results_to_dataframe
from postgwas.core.values import optional_text

from postgwas.modules.qc_summary.assessment import run_vcf_qc_assessment
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
    "extract_chromosome_from_filename",
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
        value = (
            "r=%s; mean |AF difference|=%s"
            % (
                "unavailable" if correlation is None else "%.6f" % correlation,
                "unavailable" if difference is None else "%.6f" % difference,
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
            "loss", "Missing required values", _count(missing_required),
            indent=8, label_width=30,
        ),
        screen_field(
            "loss", "Invalid coordinates", _count(invalid_coordinates),
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


def _chromosome_read_options(policies):
    # type: (Any) -> Tuple[List[str], int]
    """Null tokens and schema-inference depth for the per-chromosome read.

    Both values come from the resolved canonical policy registry.
    """
    resolved = policies if policies is not None else default_policies()
    return (
        list(resolved.get("input.chromosome_null_values")),
        int(resolved.get("input.chromosome_schema_inference_rows")),
    )


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
    """Return the three primary VCFs, failing if any promised file is absent."""
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
            "Harmonisation finished its analysis but cannot return all promised VCF "
            "outputs. Missing or empty: %s" % "; ".join(missing)
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


def _reject_counts_from_file(path, delimiter):
    # type: (Any) -> Dict[str, int]
    """Per-reason counts read back from a written reject file."""
    try:
        frame = pl.read_csv(str(path), separator=delimiter, infer_schema_length=0)
    except Exception:
        return {}
    if "reject_reason" not in frame.columns or frame.height == 0:
        return {}
    grouped = frame.group_by("reject_reason").len()
    counts = {}  # type: Dict[str, int]
    for reason, number in zip(
        grouped["reject_reason"].to_list(), grouped["len"].to_list()
    ):
        if reason is None:
            continue
        counts[str(reason)] = int(number)
    return counts


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


def _derive_parallelism(policies, requested_workers, n_chromosomes, logger=None):
    # type: (Any, int, int, Any) -> Tuple[int, int]
    """How many workers, and how many threads each may use.

    ``execution.threads_per_chromosome`` used to be the literal ``5`` passed to
    every worker regardless of ``--threads``, so the pipeline's real CPU use
    was about five times what the user asked for.  When
    ``execution.total_cpu_budget`` is set, workers x threads is held inside it:
    the worker count the user asked for is preserved for as long as possible
    and the per-chromosome thread count shrinks first, because a worker that
    cannot start does no work at all.
    """
    threads = int(policies.get("execution.threads_per_chromosome"))
    workers = max(1, min(int(requested_workers or 1), max(1, int(n_chromosomes or 1))))
    budget = policies.get("execution.total_cpu_budget")

    if budget in (None, ""):
        if logger is not None:
            logger.info(
                "Running %d chromosome%s at a time, %d thread%s each (no total CPU budget "
                "is set, so the product %d is not capped)."
                % (workers, "" if workers == 1 else "s",
                   threads, "" if threads == 1 else "s", workers * threads)
            )
        return workers, threads

    budget = max(1, int(budget))
    original = (workers, threads)
    workers = max(1, min(workers, budget))
    threads = max(1, budget // workers)
    while workers > 1 and workers * threads > budget:
        workers -= 1
        threads = max(1, budget // workers)

    if logger is not None:
        logger.decide(
            "How much CPU to use",
            {
                "total_cpu_budget": budget,
                "workers requested": original[0],
                "threads per chromosome requested": original[1],
            },
            "%d worker%s x %d thread%s = %d, inside the budget of %d"
            % (workers, "" if workers == 1 else "s",
               threads, "" if threads == 1 else "s",
               workers * threads, budget),
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
    logged.  ``status`` is ``'ok'`` or ``'failed'``.  This function never
    raises and never writes to stdout: the parent prints the returned block.
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
        null_values, infer_rows = _chromosome_read_options(pol)
        with logger.step(
            1, CHR_STEP_TOTAL, "Load the chromosome file", "polars.read_csv"
        ) as ctx:
            df = pl.read_csv(
                chr_file,
                separator=gwas2vcf_input["delimiter"],
                null_values=null_values,
                infer_schema_length=infer_rows,
                schema_overrides=canonical_variant_schema(sample_column_dict),
            )
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
            require_default_eaf = bool(pol.get("strand.enabled")) or (
                isinstance(study_decisions, dict)
                and study_decisions.get("eaf_is_maf") is True
            )
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
        effect_validation_state = {}
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
                rejects=rejects,
                validation_state=effect_validation_state,
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
            upstream_validation=effect_validation_state,
            step_number=10,
            step_total=CHR_STEP_TOTAL,
        )
        effect_validation_state.clear()

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
    except BaseException as exc:                      # noqa: BLE001 - deliberate
        status = "failed"
        qc_dict = {"error": "%s: %s" % (type(exc).__name__, exc)}
        # logger.step() already logged FAILED with the traceback when the
        # failure happened inside a step; this line names the chromosome so the
        # parent's block always ends with a clear verdict.
        logger.error(
            "Chromosome %s did not complete on attempt %d. %s: %s"
            % (chromosome, attempt, type(exc).__name__, exc)
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
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp_ctx) as executor:
            future_to_chr = {}
            for chrom in pending:
                chromosome_args = dict(worker_kwargs)
                chromosome_args.update(
                    chromosome=chrom,
                    chr_file=chr_file_by_chrom[chrom],
                    sample_column_dict=dict(worker_kwargs["sample_column_dict"]),
                    attempt=attempts.get(chrom, 0) + 1,
                    prevalidated_resource_map=dict(resource_maps[chrom]),
                )
                future = executor.submit(
                    process_one_chromosome,
                    **chromosome_args,
                )
                future_to_chr[future] = chrom
                started[chrom] = time.time()

            for future in as_completed(future_to_chr):
                chrom = future_to_chr[future]
                attempt = attempts.get(chrom, 0) + 1
                elapsed = time.time() - started.get(chrom, time.time())
                try:
                    chrom_res, qc, screen_text, status = future.result()
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
                    continue
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
    return results


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
    policies,
    logger,
    input_line_count=None,
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

    # ------------------------------------------------------------------
    # 03/08  read_summary_statistics — the expensive read
    # ------------------------------------------------------------------
    with logger.step(
        3, DATASET_STEP_TOTAL, "Read the summary statistics", "summary_statistics_io.read_summary_statistics",
        policy_keys=[
            "input.delimiter", "input.null_values", "input.strip_double_hash_lines",
            "columns.mandatory",
            "chromosome.allowed", "allele.pattern",
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
            non_standard_allele_count,
            removed_duplicates_count,
            removed_missing_count,
        ) = read_summary_statistics(
            sumstat_file=sumstat_file,
            output_dir=str(output_dir_path),
            sample_column_dict=sample_column_dict,
            output_layout=output_layout,
            table_delimiter=gwas2vcf_input["delimiter"],
            policies=pol,
            logger=logger,
            input_line_count=input_line_count,
        )
        ctx.set_rows(df.height)
        ctx.rows_in = polars_rows

    ready_variant_types = _ready_variant_type_counts(df, sample_column_dict)
    _announce(
        logger,
        "\n" + _input_validation_summary_block(
            input_variants=file_cvariant_count,
            variants_read=polars_rows,
            missing_required=removed_missing_count,
            invalid_coordinates=removed_coords,
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
            policies=pol,
            ctx=ctx,
            scope="dataset '%s'" % sample_id,
        )
        ok, problems = validate_content(sample_column_dict, df, policies=pol, logger=logger)
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
        dataset_extra["missing_sample_size"] = dict(missing_sample_size_plan)
        ctx.extra["missing_sample_size"] = dict(missing_sample_size_plan)
        ctx.set_rows(df.height, removed=0)

    # ------------------------------------------------------------------
    # 05/08  infer_genome_build
    # ------------------------------------------------------------------
    genome_build_info = infer_genome_build(
        df,
        build_reference_files,
        build_check_colmap,
        sample_column_dict,
        logger=logger,
        policies=pol,
        step_number=5,
        step_total=DATASET_STEP_TOTAL,
    )
    grch_version = genome_build_info["inferred_build"]
    build_label_width = 26
    build_match_lines = [
        screen_field(
            "genetic",
            "%s matches" % build,
            "%s (%s%%)" % (
                _count(genome_build_info["matches"][build]),
                genome_build_info["percentages"][build],
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
        ] + build_match_lines),
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
        6, DATASET_STEP_TOTAL, "Strand detection", "strand.resolve_strand_consensus",
        rows_in=df.height,
        policy_keys=[
            "strand.enabled", "strand.consensus_threshold",
            "strand.min_informative_variants",
        ],
    ) as ctx:
        strand_decision = resolve_strand_consensus(
            genome_build_info, grch_version, policies=pol,
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
    # 07/08  resolve_study_properties  ★ NEW
    # ------------------------------------------------------------------
    study_decisions = resolve_study_properties(
        df,
        sample_column_dict,
        policies=pol,
        logger=logger,
        step_number=7,
        step_total=DATASET_STEP_TOTAL,
    )
    study_decisions.update(strand_decision)
    dataset_extra["study_decisions"] = {
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
        "strand": study_decisions.get("strand"),
        "strand_consensus": {
            key: study_decisions.get(key)
            for key in (
                "forward", "reverse", "ambiguous", "informative",
                "dominant_fraction", "consensus_threshold", "reason",
            )
        },
    }
    _announce(
        logger,
        "\n" + _study_decisions_block(study_decisions),
    )
    _announce(
        logger,
        "\n" + _validated_column_mapping_block(sample_column_dict, study_decisions),
    )

    # The exact resources cannot be known at command preflight because both the
    # genome build and the chromosomes present in the study are data-derived.
    # Resolve and validate them now, before partition I/O or worker launch.
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
    require_default_eaf = bool(pol.get("strand.enabled")) or (
        study_decisions.get("eaf_is_maf") is True
    )
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
            policies=pol,
            require_default_eaf=require_default_eaf,
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

    dataset_extra["resource_preflight"] = resource_preflight
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
    # 08/08  split_by_chromosome
    # ------------------------------------------------------------------
    with logger.step(
        8, DATASET_STEP_TOTAL, "Split by chromosome", "chromosome_partition.write_chromosome_partitions",
        rows_in=df.height,
    ) as ctx:
        chr_file_by_chrom = write_chromosome_partitions(
            df,
            sample_gwas_dict=sample_column_dict,
            output_layout=output_layout,
            delimiter=gwas2vcf_input["delimiter"],
            source_snapshot=(
                source_snapshot if bool(pol.get("rejects.enabled")) else None
            ),
        )
        del source_snapshot
        ctx.info(
            "Wrote %d per-chromosome file%s: %s"
            % (len(chr_file_by_chrom), "" if len(chr_file_by_chrom) == 1 else "s",
               ", ".join(sorted(chr_file_by_chrom, key=_chromosome_sort_key)))
        )
        ctx.set_rows(df.height, removed=0)

    all_chromosomes = sorted(chr_file_by_chrom, key=_chromosome_sort_key)
    n_chr = len(all_chromosomes)
    if n_chr == 0:
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

    workers, threads_per_chromosome = _derive_parallelism(
        pol, max_workers, n_chr, logger=logger
    )
    dataset_extra["parallelism"] = {
        "workers": workers,
        "threads_per_chromosome": threads_per_chromosome,
        "requested_workers": max_workers,
        "total_cpu_budget": pol.get("execution.total_cpu_budget"),
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
                "count", "Parallel work", "%d at a time; %d thread%s each" % (
                    workers, threads_per_chromosome,
                    "" if threads_per_chromosome == 1 else "s",
                ), indent=8, label_width=20,
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
    # Rejected variants: one file per chromosome -> one dataset file, plus the
    # reason x chromosome table.
    # ------------------------------------------------------------------
    counts_by_source = {}  # type: Dict[str, Dict[str, int]]
    reject_paths = []      # type: List[str]
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

        alternatives = sorted(configured_output_matches(
            output_dir_path, output_layout["input_reject"] + "*",
            dataset_id=sample_id,
        ))
        input_reject = alternatives[0] if alternatives else None
        if input_reject is None or not Path(input_reject).is_file():
            raise PipelineError(
                "The required input-stage rejected-variants file is missing. "
                "The dataset cannot be finalised without rejection provenance.",
                results=per_chr_qc,
                status="FAILED",
            )
        reject_paths.append(str(input_reject))
        counts_by_source["input"] = _reject_counts_from_file(
            input_reject, gwas2vcf_input["delimiter"],
        )

        for chrom in all_chromosomes:
            alternatives = sorted(configured_output_matches(
                output_dir_path, output_layout["chromosome_reject"] + "*",
                dataset_id=sample_id, chromosome=chrom,
            ))
            if alternatives:
                reject_paths.append(str(alternatives[0]))
            elif chrom in completed:
                raise PipelineError(
                    "Chromosome %s completed without its required rejected-variants "
                    "file. The dataset cannot be finalised without rejection provenance."
                    % chrom,
                    results=per_chr_qc,
                    chromosomes=[chrom],
                    status="FAILED",
                )
            summary = summaries.get(chrom) or {}
            counts_by_source[chrom] = dict(summary.get("reject_counts") or {})

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

        try:
            combined = concat_reject_files(
                reject_paths,
                str(configured_output_path(
                    output_dir_path,
                    output_layout["dataset_reject"],
                    dataset_id=sample_id,
                )),
                delimiter=gwas2vcf_input["delimiter"],
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

    # ------------------------------------------------------------------
    # Dataset-level reconciliation over the chromosomes that completed.
    # ------------------------------------------------------------------
    rows_in_total = 0
    rows_out_total = 0
    rejected_total = 0
    for chrom in completed:
        summary = summaries.get(chrom) or {}
        rows_in_total += int(summary.get("rows_in") or 0)
        rows_out_total += int(summary.get("rows_out") or 0)
        rejected_total += int(summary.get("rejected") or 0)
    dataset_extra["reconciliation"] = {
        "chromosomes": completed,
        "rows_read": rows_in_total,
        "rejected": rejected_total,
        "rows_exported": rows_out_total,
        "balanced": (rows_in_total - rejected_total) == rows_out_total,
    }
    if not bool(pol.get("rejects.enabled")):
        # Without a collector there is no rejected count to balance against, so
        # reporting a total here would be arithmetic on a number nobody kept.
        dataset_extra["reconciliation"] = {
            "chromosomes": completed,
            "rows_read": rows_in_total,
            "rejected": None,
            "rows_exported": None,
            "balanced": None,
            "note": "rejects.enabled is false, so nothing could be reconciled",
        }
    elif completed:
        logger.info(
            "Across the %d chromosome%s that completed: %s read, %s rejected, "
            "%s exported. %s"
            % (
                len(completed), "" if len(completed) == 1 else "s",
                _count(rows_in_total), _count(rejected_total), _count(rows_out_total),
                "Balances." if dataset_extra["reconciliation"]["balanced"]
                else "DOES NOT BALANCE.",
            )
        )

    per_chr_qc["total_variant_infile"] = file_cvariant_count
    per_chr_qc["total_variant_read"] = polars_rows
    per_chr_qc["total_variant_removed_null_coords"] = removed_coords
    per_chr_qc["total_variant_removed_non_standard_alleles"] = non_standard_allele_count
    per_chr_qc["total_variant_remaining_for_harmonisation"] = df.height
    per_chr_qc["total_variant_ready_snps"] = ready_variant_types["snps"]
    per_chr_qc["total_variant_ready_indels_or_other"] = ready_variant_types[
        "indels_or_other_variants"
    ]
    per_chr_qc["total_variant_removed_missing_values"] = removed_missing_count
    per_chr_qc["total_variant_removed_duplicates"] = removed_duplicates_count

    dataset_extra["chromosomes"] = dict(
        (chrom, {"status": statuses[chrom], "attempts": attempts[chrom]})
        for chrom in all_chromosomes
    )
    dataset_extra["completed"] = completed
    dataset_extra["failed"] = failed
    dataset_extra["rounds"] = round_no
    dataset_extra["genome_build"] = genome_build_info
    dataset_extra["chromosome_summaries"] = summaries
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

def _save_qc_results(qc_results: Dict[str, Any], out_file: Path) -> None:
    """Save the QC dictionary to disk as JSON.

    ``default=str`` means a value the JSON encoder does not understand becomes
    its text form instead of demoting the whole file to ``str(dict)``, which
    ``qc_results_to_dataframe`` cannot read back.
    """
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w", encoding="utf-8") as f:
        json.dump(qc_results, f, indent=4, default=str)


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
    population_frequency_config = default_cfg_obj["population_frequency_qc"]
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
                policies=policies,
                logger=logger,
                input_line_count=input_line_count,
            )
        except PipelineError as exc:
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
                )
            manifest["status"] = getattr(exc, "status", "FAILED")
            manifest["failed_chromosomes"] = getattr(exc, "chromosomes", [])
            raise

        dataset_extra = qc_results.pop("__dataset__", {})
        manifest["dataset"] = dataset_extra

        qc_summary_path = configured_output_path(
            output_folder,
            output_layout["qc_summary"],
            error_type=PipelineError,
            dataset_id=sample_id,
        )
        _save_qc_results(qc_results=qc_results, out_file=qc_summary_path)

        try:
            postgwas_qc_df2 = qc_results_to_dataframe(
                data=qc_results,
                allowed_chromosomes=policies.get("chromosome.allowed"),
            )
            if postgwas_qc_df2 is None or getattr(postgwas_qc_df2, "empty", True):
                raise ValueError("qc_results_to_dataframe produced no rows.")
        except Exception as e:
            _announce(
                logger,
                screen_field(
                    "warning", "QC summary",
                    "per-chromosome table unavailable (%s); using dataset totals" % e,
                    indent=4, label_width=18,
                ),
                marker="WARNING",
            )
            postgwas_qc_df2 = _fallback_qc_frame(qc_results)

        # =====================================================
        # POST MERGE 01/05 — concatenate the per-chromosome VCFs
        # =====================================================
        with logger.step(
            1, POST_MERGE_STEP_TOTAL, "Concatenate the per-chromosome VCFs",
            "vcf_processing.concat_vcfs_by_build",
            policy_keys=["vcf.concat_max_workers", "vcf.concat_max_attempts",
                         "vcf.on_merge_failure", "vcf.min_valid_size_bytes"],
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
        with logger.step(
            3, POST_MERGE_STEP_TOTAL, "Merge the side files and clean up",
            "cleanup.finalise_harmonisation_outputs",
        ) as ctx:
            cleaned = finalise_harmonisation_outputs(
                output_dir=sample_column_dict["output_folder"],
                gwas_outputname=sample_id,
                output_layout=output_layout,
                threads=threads,
                compression_executable=compression_executable,
                policies=policies,
                logger=logger,
            )
            ctx.extra["cleaned"] = {
                "moved": len(cleaned.get("moved", [])),
                "removed": len(cleaned.get("removed", [])),
            }

        outdir = Path(sample_column_dict["output_folder"])

        # =====================================================
        # Prepare the pre-VCF counts used by the combined QC report
        # =====================================================
        target_build = str(policies.get("qc.target_build"))
        if target_build == "input":
            target_build = grch_version
        manifest["qc_target_build"] = target_build

        raw_vcf_path = configured_output_path(
            outdir,
            output_layout["merged_build_vcf"],
            error_type=PipelineError,
            dataset_id=sample_id,
            build=target_build,
        )

        # ---------------------------------------------------------
        # The gwas2vcf side summary supplies the number sent to conversion
        # ---------------------------------------------------------
        gwas2vcf_summary_path = configured_output_path(
            outdir,
            output_layout["gwas2vcf_summary"],
            error_type=PipelineError,
            dataset_id=sample_id,
        )
        try:
            gwas_to_vcf_qc_df = pd.read_csv(gwas2vcf_summary_path, sep="\t")
            gwas_to_vcf_qc_df.columns = [
                x.strip() for x in gwas_to_vcf_qc_df.columns
            ]
            gwas_to_vcf_qc_df = gwas_to_vcf_qc_df[
                gwas_to_vcf_qc_df["chromosome"] != "chromosome"
            ]
            snp_id_df = gwas_to_vcf_qc_df[
                gwas_to_vcf_qc_df["key"] == "snp_id_col"
            ]
            total_variants_in_gwas2vcf_input = int(
                snp_id_df["num_rows"].astype(int).sum()
            )
            gwas_to_vcf_qc_df.to_csv(
                gwas2vcf_summary_path, sep="\t", index=None
            )
        except (KeyError, OSError, ValueError, pd.errors.ParserError) as exc:
            message = (
                "Cannot prepare the pre-VCF QC counts from %s: %s. The raw VCF "
                "was not assessed; inspect the GWAS-to-VCF chromosome transcripts."
                % (gwas2vcf_summary_path, exc)
            )
            logger.error(message)
            raise PipelineError(message) from exc

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
            "total_variant_removed_reference_unmatched": (
                dataset_extra.get("reject_counts") or {}
            ).get("reference_unmatched"),
            "total_variant_removed_reference_ambiguous": (
                dataset_extra.get("reject_counts") or {}
            ).get("reference_ambiguous"),
            "total_variant_removed_chromosome_harmonisation": (
                dataset_extra.get("reconciliation") or {}
            ).get("rejected"),
            "total_variant_passed_chromosome_harmonisation": (
                dataset_extra.get("reconciliation") or {}
            ).get("rows_exported"),
            "total_variant_with_missing_eaf": missing_eaf_count,
            "total_variant_with_invalid_beta_se": invalid_beta_se_count,
        }
        manifest["pre_vcf_summary"] = pre_vcf_summary

        # =====================================================
        # POST MERGE 04/05 — one VCF extraction, raw QC and virtual filters
        # =====================================================
        with logger.step(
            4, POST_MERGE_STEP_TOTAL, "Assess raw VCF quality",
            "qc_summary.assessment.run_vcf_qc_assessment",
            policy_keys=[
                "qc.target_build",
                "qc.sample_size_outlier_standard_deviations",
                "filter.maf_cutoff", "filter.af_diff_cutoff", "filter.af_missing",
                "filter.info_cutoff", "filter.info_max", "filter.info_missing",
                "filter.lp_cutoff", "filter.lp_missing", "filter.include_indels",
                "filter.exclude_palindromic", "filter.palindromic_af_lower",
                "filter.palindromic_af_upper", "filter.remove_mhc",
                "filter.mhc_chrom", "filter.mhc_start", "filter.mhc_end",
            ],
        ) as ctx:
            assessment = run_vcf_qc_assessment(
                vcf_path=raw_vcf_path,
                output_directory=outdir,
                dataset_id=sample_id,
                genome_build=target_build,
                external_af_name=comparison_cols,
                vcf_fields=vcf_config["qc_fields"],
                output_layout=output_layout,
                table_delimiter=vcf_config["table_delimiter"],
                table_null_values=vcf_config["table_null_values"],
                table_null_output=vcf_config["table_null_output"],
                temporary_table_suffix=vcf_config["temporary_table_suffix"],
                io_buffer_bytes=vcf_config["io_buffer_bytes"],
                policies=policies,
                bcftools_bin=bcftools_bin,
                logger=logger,
            )
            ctx.rows_in = assessment["raw"]["num_records"]
            ctx.set_rows(assessment["qc_passed"]["num_records"])
            ctx.extra["qc_assessment"] = {
                "raw_variants": assessment["raw"]["num_records"],
                "qc_passed_variants": assessment["qc_passed"]["num_records"],
                "excluded": assessment["excluded_total"],
                "reports": assessment["reports"],
            }
        manifest["qc_assessment"] = {
            "definition": assessment["definition"],
            "raw_variants": assessment["raw"]["num_records"],
            "qc_passed_variants": assessment["qc_passed"]["num_records"],
            "excluded": assessment["excluded_total"],
            "retained_fraction": assessment["retained_fraction"],
            "reports": assessment["reports"],
        }

        # =====================================================
        # POST MERGE 05/05 — final QC report, cleanup and logs
        # =====================================================
        with logger.step(
            5, POST_MERGE_STEP_TOTAL, "Final QC report and logs",
            "service.run_harmonisation_pipeline",
        ) as ctx:
            ctx.info(
                "QC assessment reports: %s"
                % ", ".join(assessment["reports"].values())
            )

        _announce(
            logger,
            "\n" + "\n".join(
                harmonisation_qc_summary_lines(pre_vcf_summary, assessment)
            ) + "\n",
        )

        primary_outputs = _resolved_harmonisation_outputs(
            outdir, sample_id, grch_version, output_layout, vcf_config,
        )

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
        combined_log = _write_combined_log(
            outdir,
            sample_id,
            dataset_extra.get("chromosomes") or {},
            output_layout,
        )
        manifest["combined_log"] = combined_log
        manifest["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
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
