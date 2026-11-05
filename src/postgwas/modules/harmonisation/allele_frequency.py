"""Step 04 (per chromosome) — reference orientation and EAF harmonisation.

Plan v3 sections 5.2 and 5.3.

What this module does, in order:

1.  Reuse prepared study-reference matches, or join the raw chromosome reference
    when no study AF comparison was needed. Resolve
    forward/swapped/reverse-complement orientations, and transform allele-specific
    statistics consistently.
2.  If the study supplies its own frequency column, describe it, apply the configured
    out-of-range action without silently clipping under other policies, and follow
    the dataset AF acceptance decision. Before chromosome workers,
    a MAF-like study distribution is checked against the configured default EAF table:
    biallelic diagnostic matches use ALT frequency or its opposite-allele
    complement. This follows
    the GWAS-SSF requirement that frequency is aligned to the listed effect allele
    (https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics).
3.  A provisionally selected replacement first passes the separate whole-study
    A-versus-folded-B consistency and major-profile gate after dataset staging.
    A remains original study AF; B is aligned to the original listed effect
    allele. Only final passed, row/source-bound approval permits adoption.
    Otherwise validate the initially selected external source once after staging
    and before chromosome workers, retaining separate external-source evidence.
    Then take frequency from its approved cache, matching each variant in
    the two orientations named by ``eaf.match_orientations`` (direct and swap —
    never complement). For an ordinary sample-sheet external source, before swap
    alignment, screen the deduplicated raw
    frequency column for a MAF-like distribution. A non-MAF-like result retains
    the user's declared ALT/EAF contract. A MAF-like result must be confirmed
    against an independent default EAF panel; using the same file, confirmed
    MAF, or inconclusive evidence stops rather than guessing or auto-converting.
    Strand orientation has already been validated against reference REF/ALT. The
    external file's configured allele mapping determines whether its frequency
    is used as listed or reconstructed from the complete supported allele set;
    it never re-decides strand or changes study beta, Z, or ``strand_action``.
4.  Re-run the configured duplicate validation after REF/ALT and EAF alignment,
    catching opposite-strand representations that now identify the same
    physical variant.
5.  Apply policy ``eaf.degenerate`` to frequencies of exactly 0 or 1; the
    default removes them, and Z-only reconstruction always rejects any endpoint
    retained by an explicit alternative because its denominator would be zero.

Logging and rejected variants
-----------------------------
Every public function accepts ``logger=`` (a ``PipelineLogger``), ``ctx=`` (a
``StepContext``) and ``rejects=`` (a ``RejectCollector``) as optional keyword
arguments. Pipeline runs use one structured chromosome log and one unified
reject file.

Policies
--------
Every threshold is a policy. ``policies=None`` means "use the canonical YAML
defaults"; ``eaf.degenerate`` defaults to ``reject`` as documented there.
"""

import math
import os
from copy import deepcopy
from dataclasses import asdict
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Any, Dict, Optional, Tuple

import polars as pl

from postgwas.core.paths import configured_output_path
from postgwas.core.errors import MissingRequiredArgumentsError
from postgwas.core.required_arguments import RequiredArgument, require_resolved_arguments
from postgwas.core.completion import configuration_digest
from postgwas.core.checkpointing import software_identity
from postgwas.core.input_validation import capture_preflight_file_identities

from postgwas.core.values import (
    format_count,
    format_fraction_percentage,
    optional_text,
)

from .shared.runtime import (
    active_context,
    prepared_work_fingerprint,
    reject_rows,
    resolve_policies,
    step_context as shared_step_context,
)
from .shared.allele_join import allele_oriented_left_join, finish_allele_assignments, unused_column_name
from .shared.column_mapping import validate_scientific_column_mappings
from .shared.statistics import frequency_maf_screen, pool_frequency_maf_screens, valid_frequency_mask
from .shared.variant_columns import (
    canonicalize_variant_frame,
    reverse_complement_expression,
    palindromic_snp_expression,
    read_reference_variant_table,
)
from .summary_statistics_io import resolve_duplicate_variants
from .rejects import SOURCE_INPUT_ROW_COLUMN
from .strand import (
    POLICY_KEYS as STRAND_POLICY_KEYS,
    REFERENCE_AF_COLUMN,
    RESOLVED_STRAND_ACTIONS,
    STEP_LABEL as STRAND_STEP_LABEL,
    STRAND_ACTION_COLUMN,
    PREPARED_MATCHES_COLUMN,
    PREPARED_INPUT_COLUMN,
    harmonise_strand_orientation,
    prepare_study_reference_matches,
    frequency_diagnostic_candidates,
    validate_prepared_reference_matches,
    unpack_prepared_reference_matches,
)
from .chromosome_partition import _validate_source_row_ids, _write_parquet_atomically

__all__ = [
    "AlleleFrequencyError",
    "StudyAlleleFrequencyError",
    "validate_study_allele_frequency",
    "validate_replacement_frequency_pairs",
    "validate_replacement_allele_frequency",
    "validate_external_allele_frequency",
    "capture_external_frequency_sources",
    "STEP_LABEL",
    "STEP_TITLE",
    "POLICY_KEYS",
    "summarise_allele_frequency",
    "validate_allele_frequency",
    "confirm_eaf_with_reference",
    "merge_external_allele_frequencies",
    "harmonise_allele_frequency",
]


# -------------------------------------------------------
# 0. Module constants and the typed error
# -------------------------------------------------------
STEP_LABEL = "04 eaf_harmonisation"
STEP_TITLE = "Allele orientation and effect allele frequency"
FUNC_NAME = "allele_frequency.harmonise_allele_frequency"

#: The settings block printed at the top of the step.
POLICY_KEYS = (
    "eaf.out_of_range",
    "eaf.clip_tolerance",
    "eaf.degenerate",
    "eaf.maf_decision_cutoff",
    "eaf.maf_reference_min_overlap",
    "eaf.maf_reference_min_overlap_fraction",
    "eaf.maf_reference_correlation_method",
    "eaf.maf_reference_min_correlation",
    "eaf.maf_reference_max_mean_absolute_difference",
    "eaf.study_reference_min_minor_fraction",
    "eaf.study_reference_min_correlation",
    "eaf.study_reference_min_major_fraction",
    "eaf.reference_minor_fraction_cutoff",
    "eaf.maf_reference_error_margin",
    "eaf.external_min_match_fraction",
    "external_reference.exact_duplicate_action",
    "external_reference.non_identical_duplicate_action",
    "eaf.match_orientations",
    "chromosome.strip_chr_prefix",
    "chromosome.strip_leading_zero",
    "chromosome.rename_map",
    "eaf.invalid_fraction_cutoff",
    "eaf.missing_fraction_cutoff",
    "duplicates.key",
    "duplicates.consistency_fields",
    "duplicates.post_orientation_relative_tolerance",
    "duplicates.quality_fields",
    "duplicates.selection_order",
    "duplicates.conflicting_action",
)

class AlleleFrequencyError(RuntimeError):
    """A fatal problem in the EAF step.

    Derives from ``RuntimeError`` so the dataset orchestrator can catch it,
    record the chromosome failure, and apply its retry policy.
    """


class StudyAlleleFrequencyError(AlleleFrequencyError):
    """Stop the owning dataset's AF preflight without retrying it.

    The batch runner may continue with independent datasets. Optional measured
    study evidence is copied so failure reporting never needs another analysis.
    """

    def __init__(self, *args, study_decisions=None):
        super().__init__(*args)
        self.study_decisions = (
            deepcopy(study_decisions) if study_decisions is not None else None
        )


def validate_study_allele_frequency(
    df, sample_columns, study_decisions, resource_maps, reference_mapping,
    policies=None, logger=None, *, eaf_replacement_file=None, eaf_replacement_column=None,
):
    """Return the AF decision and study frame with reusable reference matches.

    The distribution screen is resolved once by study_properties. Only a
    MAF-like study needs reference concordance. Pool diagnostic pairs across
    the dataset, not chromosome averages, and compare unfolded frequencies
    aligned to the original effect allele. Thresholds are configurable QC
    heuristics, not proof of the semantic meaning of every individual value.
    """
    pol = resolve_policies(policies)
    emit = active_context(logger)
    decisions = dict(study_decisions)
    replacement_file = optional_text(eaf_replacement_file)
    replacement_column = optional_text(eaf_replacement_column)
    replacement_requirements = [
        RequiredArgument("--eaf-replacement-file", "modules.harmonisation.eaf_replacement_file", replacement_file),
        RequiredArgument("--eaf-replacement-column", "modules.harmonisation.eaf_replacement_column", replacement_column),
    ]
    if replacement_file is not None or replacement_column is not None:
        try:
            require_resolved_arguments(replacement_requirements)
        except MissingRequiredArgumentsError as exc:
            raise StudyAlleleFrequencyError(str(exc)) from exc
    if replacement_file is not None:
        decisions["eaf_replacement_specification"] = {"file": replacement_file, "column": replacement_column}
    if any(key in df.columns for key in (PREPARED_MATCHES_COLUMN, PREPARED_INPUT_COLUMN)):
        raise StudyAlleleFrequencyError("Study input contains reserved reference-preparation columns.")
    column = decisions.get("eaf_column")
    initial = decisions.get("eaf_is_maf")
    evidence = {"status": "not_applicable", "decision_reason": "no_study_frequency"}
    decisions["eaf_selected_source"] = "study_supplied" if column is not None else "external_eaf"
    decisions["eaf_input_interpretation"] = "eaf" if column is not None else None
    if column is None:
        decisions["eaf_external_validation_required"] = True
        decisions["eaf_external_study_variant_count"] = df.height
    if column is not None:
        screen = decisions.get("eaf_evidence") or frequency_maf_screen(
            df, column, float(pol.get("eaf.maf_decision_cutoff")))
        decisions["eaf_effective_frequency_screen"] = _effective_frequency_screen([screen], pol)
        if not decisions["eaf_effective_frequency_screen"]["usable"]:
            raise StudyAlleleFrequencyError("Selected GWAS effect-allele frequency has no usable values in [0,1].",
                study_decisions=decisions)
        if decisions.get("eaf_evidence"):
            decisions["eaf_evidence"] = {**decisions["eaf_evidence"],
                "scope": "merged_gwas_rows", "allele_direction": "original_gwas_effect_allele"}
        if initial is not True and initial is not False:
            raise StudyAlleleFrequencyError(
                "Study frequency preflight requires the whole-study MAF screen."
            )
        evidence = {
            "status": "not_required",
            "decision_reason": "study_frequency_not_predominantly_minor",
        }
        if initial:
            # A whole-genome reference used by several chromosomes is read once.
            references = {}
            for chromosome, resources in resource_maps.items():
                key = (
                    resources["default_eaf_file"],
                    resources["default_eaf_reference_column"],
                )
                references.setdefault(key, []).append(str(chromosome))
            pairs, prepared = [], []
            diagnostic_excluded = diagnostic_non_biallelic = 0
            diagnostic_exclusion_counts = {}
            prepared_references = {}
            projected = df.select([SOURCE_INPUT_ROW_COLUMN] + [
                sample_columns[key + "_col"] for key in ("chr", "pos", "ea", "oa")
            ] + [column])
            for (path, population), chromosomes in references.items():
                emit.info(
                    "Study AF preflight: reference %s; population %s; chromosomes %s."
                    % (path, population, ", ".join(chromosomes))
                )
                paired, matches, metadata = prepare_study_reference_matches(
                    projected.filter(
                        pl.col(sample_columns["chr_col"]).cast(pl.Utf8).is_in(chromosomes)
                    ),
                    sample_columns, column, path, population, reference_mapping,
                    pol, warn=emit.warn, diagnostic_require_biallelic=False,
                )
                diagnostic_stats = metadata["diagnostic_statistics"]
                diagnostic_non_biallelic += diagnostic_stats["excluded_non_biallelic_rows"]
                diagnostic_excluded += diagnostic_stats.get("excluded_orientation_rows",
                    diagnostic_stats["excluded_non_biallelic_rows"])
                for name, count in diagnostic_stats.get("excluded_counts", {}).items():
                    diagnostic_exclusion_counts[name] = diagnostic_exclusion_counts.get(name, 0) + count
                pairs.append(paired)
                prepared.append(matches)
                prepared_references.update({chromosome: metadata for chromosome in chromosomes})
            paired = pl.concat(pairs, rechunk=False) if pairs else pl.DataFrame(
                schema={
                    "__study_frequency": pl.Float64,
                    "__expected_effect_frequency": pl.Float64,
                }
            )
            interpretation, comparison = _compare_frequency_pairs(
                paired, "__study_frequency", "__expected_effect_frequency", pol, emit,
                accept_direct_eaf=True, study_variant_count=df.height,
            )
            reason = comparison["decision_reason"]
            evidence = {
                **comparison,
                "status": "passed" if interpretation == "eaf" else "failed",
                "correlation_method": comparison["maf_correlation_method"],
                "minimum_correlation": comparison["minimum_study_correlation"],
                "maximum_mean_absolute_difference": comparison["maximum_mean_absolute_maf_difference"],
                "correlation": comparison["direct_eaf_correlation"],
                "mean_absolute_difference": comparison["mean_eaf_error"],
                "diagnostic_excluded_non_biallelic_rows": diagnostic_non_biallelic,
                "diagnostic_excluded_orientation_rows": diagnostic_excluded,
                "diagnostic_excluded_counts": diagnostic_exclusion_counts,
                "reference_site_membership_assumption": "the selected reference represents all alleles at each listed site",
                "decision_reason": reason,
                "maf_reference_decision": interpretation,
            }
            decisions["eaf_reference_validation"] = evidence
            if interpretation == "maf":
                # Classification is established even if configuration prevents
                # activating a replacement source. Preserve it on failure.
                decisions["eaf_input_interpretation"] = "maf"
                _select_frequency_replacement(decisions, evidence)
                emit.info("Input MAF supported by the pooled reference comparison; selecting the supplied reference EAF. "
                          "This is a whole-column QC inference, not proof of every row's provenance.")
            elif interpretation == "eaf":
                emit.info("Study EAF source retained: allele-aligned reference agreement passed. "
                          "A MAF-like distribution does not invalidate EAF. Reference frequency of the original "
                          "listed effect allele at or below 0.5: %s of %s diagnostic pairs."
                          % (_fmt_pct(evidence["reference_effect_allele_minor_fraction"]), paired.height))
            emit.info("Frequency diagnostic exclusions: %s study rows with unsupported or unresolved orientation; "
                      "source rows and complete orientation candidates were retained." % diagnostic_excluded)
            emit.decide("Whole-study AF acceptance", evidence, evidence["status"])
            if evidence["status"] not in ("passed", "replacement_selected"):
                message = (
                    "Study frequency column '%s' cannot be accepted as EAF: %s; MAF interpretation is inconclusive. "
                    "Comparable non-palindromic variants: %s (required >= %s); usable coverage: %s "
                    "of %s study rows (required >= %s); %s correlation: %s "
                    "(required >= %s); reference EA AF <= 0.5 fraction: %s (required > %s); "
                    "mean absolute allele-aligned AF difference: %s (required <= %s). "
                    "Rule 1 failed requirements: %s; rule 2 failed requirements: %s. "
                    "Supplied AF versus reference MAF correlation: %s; mean absolute MAF difference: %s; "
                    "reference EA AF > 0.5 fraction: %s (required > %s). "
                    "No AF values were changed and no replacement was selected. "
                    "This dataset stopped before chromosome analysis; no dataset retry will be attempted."
                    % (column, reason, paired.height, evidence["minimum_overlap"],
                       _fmt_pct(evidence["observed_overlap_fraction"]), evidence["study_variant_count"],
                       _fmt_pct(evidence["minimum_overlap_fraction"]), evidence["correlation_method"],
                       _fmt_metric(evidence["correlation"]), evidence["minimum_correlation"],
                       _fmt_pct(evidence["reference_effect_allele_minor_fraction"]),
                       _fmt_pct(evidence["minimum_reference_effect_allele_minor_fraction"]),
                       _fmt_metric(evidence["mean_absolute_difference"]),
                       evidence["maximum_mean_absolute_difference"],
                       ", ".join(evidence["rule_1_failed_requirements"]) or "none",
                       ", ".join(evidence["rule_2_failed_requirements"]) or "none",
                       _fmt_metric(evidence["supplied_to_reference_maf_correlation"]),
                       _fmt_metric(evidence["mean_maf_error"]),
                       _fmt_pct(evidence["reference_effect_allele_major_fraction"]),
                       _fmt_pct(evidence["minimum_reference_effect_allele_major_fraction"]))
                )
                emit.error(message)
                raise StudyAlleleFrequencyError(message, study_decisions=decisions)
            # Preparation adds internal evidence only; preserve all source rows
            # and their order, including rows without diagnostic AF matches.
            matches = pl.concat(prepared, rechunk=False)
            if matches.height != df.height or matches[SOURCE_INPUT_ROW_COLUMN].n_unique() != df.height:
                raise StudyAlleleFrequencyError(
                    "Reference preparation did not cover every study row exactly once.",
                    study_decisions=decisions,
                )
            order_column = unused_column_name("__postgwas_preparation_order", df)
            df = df.with_row_index(order_column).join(
                matches, on=SOURCE_INPUT_ROW_COLUMN, how="left", coalesce=True,
            ).sort(order_column).drop(order_column)
            decisions["eaf_prepared_references"] = prepared_references
        else:
            emit.info(
                "Study AF accepted as declared EAF: the whole-study MAF screen "
                "did not exceed its cutoff; no additional frequency comparison was required."
            )
        decisions["eaf_is_maf_initial"] = initial
        decisions["eaf_is_maf_initial_source"] = decisions.get("eaf_is_maf_source")
        replacement_selected = decisions["eaf_selected_source"] == "replacement_reference"
        decisions["eaf_is_maf"] = replacement_selected
        decisions["eaf_is_maf_source"] = (
            "dataset_reference_maf_support" if replacement_selected else
            "dataset_reference_concordance" if initial else "study_level_statistic"
        )
    decisions["eaf_reference_validation"] = evidence
    return decisions, df


_fmt_int = format_count
_fmt_pct = format_fraction_percentage


def _fmt_metric(value: Any) -> str:
    """Render one optional finite QC metric without hiding missing evidence."""
    if value is None:
        return "not available"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return "%.6g" % number if math.isfinite(number) else "not available"


def _maf_reference_suitability_text(evidence: Dict[str, Any]) -> str:
    """Summarise the YAML-controlled folded-MAF suitability evidence."""
    return (
        "folded-MAF %s correlation: %s (minimum %s); mean absolute "
        "folded-MAF difference: %s (maximum %s)"
        % (
            evidence.get("maf_correlation_method", "unknown"),
            _fmt_metric(evidence.get("maf_correlation")),
            _fmt_metric(evidence.get("minimum_maf_correlation")),
            _fmt_metric(evidence.get("mean_absolute_maf_difference")),
            _fmt_metric(evidence.get(
                "maximum_mean_absolute_maf_difference"
            )),
        )
    )


def same_frequency_reference_source(path, column, mapping, reference_path, reference_column, reference_mapping):
    """Identify one declared AF source across aliases without comparing frequencies."""
    if (optional_text(path) is None or optional_text(reference_path) is None
            or optional_text(column) != optional_text(reference_column)
            or mapping is None or reference_mapping is None):
        return False
    try:
        same_file = os.path.samefile(path, reference_path)
    except OSError:
        same_file = os.path.realpath(path) == os.path.realpath(reference_path)
    delimiter, reference_delimiter = mapping.get("delimiter"), reference_mapping.get("delimiter")
    return (
        same_file
        and (delimiter == reference_delimiter or "auto" in (delimiter, reference_delimiter))
        and all(mapping.get(key) == reference_mapping.get(key) for key in ("chr", "pos", "a1", "a2"))
    )


def _frequency_validation_label(initial_external):
    return "Whole-study external AF validation" if initial_external else "EAF replacement validation"


def _effective_frequency_screen(screens, policies):
    """Pool GWAS-row counts once in original effect-allele direction."""
    return {**pool_frequency_maf_screens(
        ({"column_present": True, **screen} for screen in screens), float(policies.get("eaf.maf_decision_cutoff"))),
        "scope": "merged_gwas_rows", "allele_direction": "original_gwas_effect_allele"}


def _select_frequency_replacement(decisions, evidence):
    """Resolve the one explicit replacement after supported initial MAF."""
    specification = decisions.get("eaf_replacement_specification") or decisions.get("eaf_replacement") or {}
    try:
        require_resolved_arguments([
            RequiredArgument("--eaf-replacement-file", "modules.harmonisation.eaf_replacement_file", specification.get("file")),
            RequiredArgument("--eaf-replacement-column", "modules.harmonisation.eaf_replacement_column", specification.get("column")),
        ])
    except MissingRequiredArgumentsError as exc:
        evidence["failure_reason"] = "replacement_arguments_missing"
        raise StudyAlleleFrequencyError("Input MAF supported by the configured rules. " + str(exc),
            study_decisions=decisions) from exc
    evidence.update(status="replacement_selected", decision_reason="maf_interpretation_supported")
    decisions.update(eaf_selected_source="replacement_reference", eaf_input_interpretation="maf",
        eaf_is_maf=True, eaf_reference_validation=deepcopy(evidence),
        eaf_replacement={"file": specification["file"], "column": specification["column"]})


def _require_initial_frequency_inputs(decisions, *, verify_partitions=False):
    for item in (decisions.get("eaf_initial_frequency") or {}).get("inputs", []):
        if item["identity"] != _replacement_file_identity(item["path"], label="Initial selected AF input"):
            raise StudyAlleleFrequencyError("Initial selected AF inputs changed; rerun dataset preparation.")
    if verify_partitions:
        for expected in (decisions.get("eaf_initial_frequency") or {}).get("partitions", {}).values():
            if expected != prepared_work_fingerprint(expected["path"], error_type=StudyAlleleFrequencyError):
                raise StudyAlleleFrequencyError("Initial selected AF work rows changed; rerun dataset preparation.")


def _staged_frequency_bindings(resource_maps, source_files, staging_summary, label, source_identities=None):
    """Bind staged images to the producer that read the declared source."""
    bindings = {}
    for chromosome, resources in resource_maps.items():
        staged = optional_text(resources.get("user_eaf_file"))
        source = source_files[chromosome]
        if staged is None:
            continue
        produced = any("user_eaf_file" in item.get("roles", [])
            and item.get("partitions", {}).get(chromosome) == staged
            for item in (staging_summary or {}).get("resources", []))
        # A direct library caller may supply an unchanged alias, without a
        # staging producer. A claimed produced image must retain producer proof
        # even if its path was subsequently retargeted to the original file.
        if not produced and os.path.samefile(staged, source):
            continue
        identity = ((source_identities or {}).get(chromosome)
            or _replacement_file_identity(source, label=label + " original input"))
        producer = next((item for item in (staging_summary or {}).get("resources", [])
            if "user_eaf_file" in item.get("roles", []) and item.get("source_identity") == identity
            and item.get("partitions", {}).get(chromosome) == staged), None)
        expected = (producer or {}).get("partition_identities", {}).get(chromosome)
        if expected is None or expected != _replacement_file_identity(staged, label=label + " staged input"):
            raise StudyAlleleFrequencyError(label + " staged inputs changed or lack producer proof; rerun dataset preparation.")
        bindings[chromosome] = expected
    return bindings


def _effect_frequency_view(frame, matches, keys, policies, assignments=None, assignment=None):
    """Select one safely aligned AF per source row from an already prepared join.

    Default-unmatched rows with an unambiguous external assignment still enter
    the screen. Ambiguous orientations stay missing; AF never selects a strand.
    Only narrow cached-value joins are used here, without reference I/O.
    """
    original, _ = canonicalize_variant_frame(frame.select(SOURCE_INPUT_ROW_COLUMN, *keys.values()),
        keys, policies, error_type=StudyAlleleFrequencyError, label="effect-frequency screen keys")
    local_keys = {role: "__postgwas_effect_key_" + role for role in keys}
    original = original.rename({keys[role]: name for role, name in local_keys.items()})
    resolved = matches.filter(pl.col("__strand_candidates") == 1)
    ref_af = "__strand_ref_af" if "__strand_ref_af" in matches.columns else None
    default_value = pl.when(pl.col("__strand_swapped")).then(
        pl.col(ref_af) if ref_af else 1.0 - pl.col("__strand_af")
    ).otherwise(pl.col("__strand_af"))
    orientation = resolved.select(SOURCE_INPUT_ROW_COLUMN, STRAND_ACTION_COLUMN,
        default_value.alias("default_frequency"))
    known = matches.select(SOURCE_INPUT_ROW_COLUMN).unique().with_columns(pl.lit(True).alias("__known_reference"))
    selected = original.join(orientation, on=SOURCE_INPUT_ROW_COLUMN, how="left", coalesce=True).join(
        known, on=SOURCE_INPUT_ROW_COLUMN, how="left", coalesce=True)
    if assignments is None:
        return selected.select(SOURCE_INPUT_ROW_COLUMN,
            pl.col("default_frequency").alias("effect_frequency"), "default_frequency")
    reverse = pl.col(STRAND_ACTION_COLUMN).is_in(["reverse_complement", "reverse_complement_swapped"])
    selected = selected.with_columns([
        pl.when(reverse).then(reverse_complement_expression(pl.col(local_keys[role]))).otherwise(pl.col(local_keys[role])).alias(local_keys[role])
        for role in ("ea", "oa")])
    value = assignment["value"]
    narrow = assignments.select(SOURCE_INPUT_ROW_COLUMN,
        *[pl.col(keys[role]).alias(name) for role, name in local_keys.items()], value)
    selected = selected.join(narrow, on=[SOURCE_INPUT_ROW_COLUMN, *local_keys.values()], how="left", coalesce=True)
    fallback = original.with_columns([
        reverse_complement_expression(pl.col(local_keys[role])).alias(local_keys[role]) for role in ("ea", "oa")
    ]).join(narrow, on=[SOURCE_INPUT_ROW_COLUMN, *local_keys.values()], how="left", coalesce=True).select(
        SOURCE_INPUT_ROW_COLUMN, pl.col(value).alias("__reverse_frequency"))
    selected = selected.join(fallback, on=SOURCE_INPUT_ROW_COLUMN, how="left", coalesce=True)
    has_reference = pl.col("__known_reference").fill_null(False)
    unique_external = pl.col(value).is_not_null() ^ pl.col("__reverse_frequency").is_not_null()
    if str(policies.get("strand.mode")) == "reference_aligned":
        unmatched_value = pl.col(value)
    else:
        unmatched_value = pl.when(unique_external).then(pl.coalesce(value, "__reverse_frequency"))
    effect = pl.when(has_reference).then(
        pl.when(pl.col(STRAND_ACTION_COLUMN).is_not_null()).then(pl.col(value))
    ).otherwise(pl.when(~palindromic_snp_expression(pl.col(local_keys["ea"]), pl.col(local_keys["oa"]))).then(unmatched_value))
    return selected.select(SOURCE_INPUT_ROW_COLUMN, effect.alias("effect_frequency"), "default_frequency")


def _replacement_file_identity(path, *, label="EAF replacement validation input"):
    """Use the invocation-local identity contract for original and staged inputs."""
    identity, = capture_preflight_file_identities(
        [path], error_type=StudyAlleleFrequencyError, label=label,
    )
    return {**asdict(identity), "path": str(identity.path)}


def capture_external_frequency_sources(resource_maps, *, expected_identities=None):
    """Bind each declared alias before staging, including its resolved identity.

    When preflight identities are available, a retargeted alias must not bless a
    different source merely because the old resolved target remains unchanged.
    """
    bindings = {}
    label = _frequency_validation_label(True)
    for chromosome, resources in resource_maps.items():
        path = optional_text(resources.get("user_eaf_file"))
        if path is None:
            raise StudyAlleleFrequencyError(label + " lacks its selected original external source.")
        identity = _replacement_file_identity(path, label=label + " input")
        if expected_identities and identity not in expected_identities:
            raise StudyAlleleFrequencyError(label + " original source changed after resource preflight; rerun dataset preparation.")
        reference_path = optional_text(resources.get("default_eaf_file"))
        if reference_path is None:
            raise StudyAlleleFrequencyError(label + " lacks its selected default reference.")
        reference_identity = _replacement_file_identity(reference_path, label=label + " default reference")
        if expected_identities and reference_identity not in expected_identities:
            raise StudyAlleleFrequencyError(label + " default reference changed after resource preflight; rerun dataset preparation.")
        bindings[chromosome] = {"path": str(Path(path).expanduser().absolute()), "identity": identity,
            "reference_path": str(Path(reference_path).expanduser().absolute()), "reference_identity": reference_identity}
    return bindings



def _replacement_configuration(columns, decisions, mapping, policies):
    """Bind measured evidence to the declared source, matching and numeric rules."""
    return configuration_digest({
        "study_columns": {key: columns.get(key) for key in
                          ("chr_col", "pos_col", "ea_col", "oa_col", "eaf_col")},
        "replacement": decisions.get("eaf_replacement"),
        "external_source": decisions.get("eaf_external_source"),
        "external_validation_required": decisions.get("eaf_external_validation_required"),
        "external_column": columns.get("eafcolumn") if decisions.get("eaf_external_validation_required") else None,
        "external_study_variant_count": decisions.get("eaf_external_study_variant_count"),
        "replacement_specification": decisions.get("eaf_replacement_specification"),
        "initial_frequency": decisions.get("eaf_initial_frequency"),
        "original_evidence": decisions.get("eaf_reference_validation"),
        "original_interpretation": decisions.get("eaf_input_interpretation"),
        "mapping": mapping, "policies": policies.resolved_dict(),
        "software": software_identity(), "polars": pl.__version__,
    })


def _require_unchanged_replacement_binding(binding):
    """Check every source plus the checksummed generated row-evidence partition."""
    label = _frequency_validation_label(binding.get("source_mode") == "external_eaf")
    original_reference = binding["prepared_reference"]["contract"]["reference"]
    for pathname, expected in (
        (binding["source_path"], binding["source"]),
        (binding.get("staged_path"), binding.get("staged")),
        (binding["reference_path"], original_reference),
    ):
        if expected is not None and expected != _replacement_file_identity(pathname, label=label + " input"):
            raise StudyAlleleFrequencyError(label + " inputs changed; rerun dataset preparation.")
    if binding["partition"] != prepared_work_fingerprint(binding["partition_path"], error_type=StudyAlleleFrequencyError):
        raise StudyAlleleFrequencyError(label + " partition changed; rerun dataset preparation.")


def _replacement_candidate_keys(frame, matches, keys, policies, *, label="replacement assignment keys"):
    """Enumerate possible output keys without using replacement AF to orient rows."""
    columns = [SOURCE_INPUT_ROW_COLUMN, *keys.values()]
    original, _ = canonicalize_variant_frame(
        frame.select(columns), keys, policies, error_type=StudyAlleleFrequencyError,
        label=label,
    )
    candidates = original.join(matches.select(
        SOURCE_INPUT_ROW_COLUMN, "__strand_alt", "__strand_ref",
    ), on=SOURCE_INPUT_ROW_COLUMN, how="inner").with_columns(
        pl.col("__strand_alt").alias(keys["ea"]),
        pl.col("__strand_ref").alias(keys["oa"]),
    ).select(columns)
    reverse = original.filter(
        (pl.col(keys["ea"]).str.len_chars() == 1) & (pl.col(keys["oa"]).str.len_chars() == 1)
    ).with_columns(
        reverse_complement_expression(pl.col(keys["ea"])).alias(keys["ea"]),
        reverse_complement_expression(pl.col(keys["oa"])).alias(keys["oa"]),
    )
    return pl.concat([candidates, original, reverse]).unique(subset=columns, maintain_order=True)


def _reserve_assignment_columns(assignments, assignment, working, *, reserved=()):
    """Keep raw and planned assignment fields separate from scientific columns.

    Preparation matches a narrow key frame, so its private names must be rebound
    before joining original AF back in. The worker uses the same descriptor
    remapping against its wider frame; neither caller renames original data.
    """
    scalar_roles = ("value", "orientation", "row", "reason", "biallelic")
    list_roles = ("frequency_columns", "site_columns")
    names = [assignment[role] for role in scalar_roles if assignment[role] is not None]
    names.extend(name for role in list_roles for name in assignment[role])
    rename = {}
    for name in names:
        # Reserve untouched keys, other raw fields, and derived destinations,
        # including names that are not materialized until finalization.
        occupied = (set(assignments.columns) | set(names) | set(rename.values())) - {name}
        rename[name] = unused_column_name(name, working, reserved=(*reserved, *occupied))
    remapped = dict(assignment)
    remapped.update({role: rename.get(assignment[role]) for role in scalar_roles})
    remapped.update({role: [rename[name] for name in assignment[role]] for role in list_roles})
    return assignments.rename({name: target for name, target in rename.items()
        if name in assignments.columns}), remapped


def _consume_replacement_assignments(frame, cache, descriptor, keys, policies, *, chromosome, emit, logger, ctx, rejects):
    """Select only the key chosen by existing strand processing; never rematch a source."""
    label = "External AF" if descriptor.get("source_mode") == "external_eaf" else "EAF replacement"
    assignment = descriptor["assignment"]
    # Public destinations depend only on scientific input, as in a fresh merge.
    destination = unused_column_name(descriptor["statistics"]["source_column"], frame)
    packed_column = descriptor["column"]
    candidates = cache.explode(packed_column).unnest(packed_column)
    join_keys = [SOURCE_INPUT_ROW_COLUMN, *keys.values()]
    if candidates.select(pl.struct(join_keys).is_duplicated().any()).item():
        raise AlleleFrequencyError(label + " cache contains duplicate assignment keys.")
    candidates, assignment = _reserve_assignment_columns(candidates, assignment, frame, reserved=(destination,))
    marker = unused_column_name("__postgwas_cached_assignment", frame, candidates, reserved=(destination,))
    order = unused_column_name("__postgwas_cached_order", frame, candidates, reserved=(destination, marker))
    selected = frame.with_row_index(order).join(candidates.with_columns(pl.lit(True).alias(marker)),
        # Null keys identify an explicit cached unmatched assignment for this
        # source row; this does not permit null-key matches to reference alleles.
        on=join_keys, how="left", coalesce=True, join_nulls=True).sort(order)
    if selected.height != frame.height or selected[marker].null_count():
        raise AlleleFrequencyError(label + " cache lacks a known assignment for an oriented study row.")
    selected = selected.drop(marker)
    assignment["row"] = order
    stats = dict(descriptor["statistics"])
    reuse_marker = ("reused_prepared_external_assignments" if descriptor.get("source_mode") == "external_eaf"
                    else "reused_prepared_replacement_assignments")
    stats.update(external_rows_loaded=0, **{reuse_marker: True},
        min_match_fraction=float(policies.get("eaf.external_min_match_fraction")))
    assignment["reason"] = unused_column_name(assignment["reason"], selected, reserved=(destination,))
    assignment["biallelic"] = None
    stats["destination_column"] = destination
    return _finalize_external_allele_frequencies(selected, assignment, stats, keys, destination, policies,
        chromosome=chromosome, emit=emit, logger=logger, ctx=ctx, rejects=rejects)


def validate_replacement_allele_frequency(
    chromosome_partitions, sample_columns, study_decisions, resource_maps,
    external_eaf_colmap, *, policies=None, logger=None, staging_summary=None,
):
    """Validate one selected replacement before workers, without changing rows.

    Original prepared candidates supply independent orientation. One complete
    narrow chromosome join feeds both diagnostic pairs and a private assignment
    spool. After all gates and input bindings pass, append assignments atomically
    to existing work partitions. Workers select their independently oriented key
    without rereading the source. Temporary files are removed on success/failure;
    original scientific rows and the diagnostic denominator remain unchanged.
    """
    if study_decisions.get("eaf_selected_source") != "replacement_reference":
        return dict(study_decisions)
    replacement = study_decisions.get("eaf_replacement") or {}
    try:
        staged = _staged_frequency_bindings(resource_maps, replacement.get("resolved_files", {}),
            staging_summary, _frequency_validation_label(False))
    except Exception as exc:
        decisions = deepcopy(study_decisions)
        reason = "replacement_inputs_changed" if isinstance(exc, StudyAlleleFrequencyError) else "replacement_inputs_invalid"
        decisions["eaf_replacement_validation"] = {"status": "invalid",
            "failure_reason": reason, "decision_reason": reason,
            "source": {**replacement, "mapping": dict(external_eaf_colmap)}, "failure_message": str(exc)}
        raise StudyAlleleFrequencyError(str(exc), study_decisions=decisions) from exc
    return _prepare_selected_allele_frequency(chromosome_partitions, sample_columns,
        study_decisions, resource_maps, external_eaf_colmap, policies=policies, logger=logger, staged_bindings=staged)


def validate_external_allele_frequency(
    chromosome_partitions, sample_columns, study_decisions, resource_maps,
    external_eaf_colmap, default_eaf_colmap, *, original_resource_maps,
    policies=None, logger=None, source_bindings=None, staging_summary=None,
):
    """Approve initially selected external AF after staging, before any worker.

    Original preflight paths remain authoritative declarations after staging;
    no native study frequency is inferred or synthesized. Staged inputs require
    source/default alias bindings captured before staging and producer identities
    from the existing staging summary. Unstaged direct calls capture before reads.
    """
    if study_decisions.get("eaf_selected_source") != "external_eaf":
        return dict(study_decisions)
    decisions = deepcopy(study_decisions)
    decisions["eaf_external_validation_required"] = True
    decisions["eaf_external_source"] = {
        "column": sample_columns.get("eafcolumn"),
        "resolved_files": {c: r.get("user_eaf_file") for c, r in original_resource_maps.items()},
        "declared_files": {c: r.get("user_eaf_file") for c, r in original_resource_maps.items()},
    }
    try:
        if source_bindings is None:
            if any(Path(resources["user_eaf_file"]).expanduser().absolute() !=
                Path(original_resource_maps[c]["user_eaf_file"]).expanduser().absolute()
                for c, resources in resource_maps.items()):
                raise StudyAlleleFrequencyError(
                    _frequency_validation_label(True) + " requires original source proof captured before external staging.")
            source_bindings = capture_external_frequency_sources(original_resource_maps)
        if source_bindings != capture_external_frequency_sources(original_resource_maps):
            raise StudyAlleleFrequencyError(
                _frequency_validation_label(True) + " original inputs or aliases changed after staging; rerun dataset preparation.")
        staged_bindings = _staged_frequency_bindings(resource_maps,
            {c: b["path"] for c, b in source_bindings.items()}, staging_summary,
            _frequency_validation_label(True), {c: b["identity"] for c, b in source_bindings.items()})
    except Exception as exc:
        decisions["eaf_external_validation"] = {"status": "invalid", "failure_reason": "external_inputs_changed",
            "decision_reason": "external_inputs_changed", "study_variant_count": decisions.get("eaf_external_study_variant_count")}
        raise StudyAlleleFrequencyError(str(exc), study_decisions=decisions) from exc
    return _prepare_selected_allele_frequency(chromosome_partitions, sample_columns,
        decisions, resource_maps, external_eaf_colmap, policies=policies, logger=logger,
        initial_external=True, default_eaf_colmap=default_eaf_colmap, source_bindings=source_bindings, staged_bindings=staged_bindings)


def _validate_initial_external_pairs(paired_paths, raw_screens, bindings, source, mapping,
                                     exclusions, total, policies, emit, evidence, *, effective_screens):
    """Screen merged GWAS values, then compare same original-EA pairs once."""
    screen = _effective_frequency_screen(effective_screens.values(), policies)
    evidence.update(status="passed", decision_reason="external_frequency_not_predominantly_minor",
        raw_frequency_screen=pool_frequency_maf_screens(raw_screens.values(), float(policies.get("eaf.maf_decision_cutoff"))),
        effective_frequency_screen=screen, study_variant_count=total, excluded_counts=exclusions,
        source={**source, "mapping": dict(mapping),
            "identities": {c: b["source"] for c, b in bindings.items()}})
    if not screen["usable"]:
        evidence.update(status="inconclusive", decision_reason="no_usable_external_frequencies")
    elif screen["maf_like"]:
        if any(os.path.samefile(b["source_path"], b["reference_path"]) for b in bindings.values()):
            evidence.update(status="inconclusive", decision_reason="independent_reference_required")
        else:
            def usable_pairs():
                for paired_path in paired_paths:
                    batch = pl.read_parquet(paired_path)
                    yield batch.filter(pl.col("eligible").fill_null(False)
                        & valid_frequency_mask(pl.col("study_frequency"))
                        & valid_frequency_mask(pl.col("replacement_frequency")))
            count = sum(batch.height for batch in usable_pairs())
            interpretation, comparison = _compare_frequency_pairs(usable_pairs,
                "study_frequency", "replacement_frequency", policies, emit, accept_direct_eaf=True,
                study_variant_count=total, pair_count=count, comparison_evidence=evidence)
            evidence.update(comparison, status="passed" if interpretation == "eaf" else "inconclusive",
                maf_reference_decision=interpretation)
    evidence["failed_requirements"] = evidence.get("shared_evidence_failed_requirements") or (
        [] if evidence["status"] == "passed" else [evidence["decision_reason"]])
    emit.decide("Whole-study external AF validation", evidence, evidence["status"])
    return evidence


def _prepare_selected_allele_frequency(
    chromosome_partitions, sample_columns, study_decisions, resource_maps,
    external_eaf_colmap, *, policies=None, logger=None, initial_external=False,
    default_eaf_colmap=None, source_bindings=None, staged_bindings=None,
):
    """Share row-bound assignment, replay, publication and approval by AF mode."""
    pol = resolve_policies(policies)
    decisions = deepcopy(study_decisions)
    emit = active_context(logger)
    evidence_key = "eaf_external_validation" if initial_external else "eaf_replacement_validation"
    label = _frequency_validation_label(initial_external)
    replacement = decisions.get("eaf_external_source" if initial_external else "eaf_replacement") or {}
    if initial_external:
        decisions[evidence_key] = {"status": "pending",
            "study_variant_count": decisions.get("eaf_external_study_variant_count"),
            "source": {**replacement, "mapping": dict(external_eaf_colmap),
                "identities": {c: b["identity"] for c, b in source_bindings.items()}}}
    else:
        decisions[evidence_key] = {"status": "pending", "source": {**replacement, "mapping": dict(external_eaf_colmap)}}
    original = decisions.get("eaf_reference_validation") or {}
    try:
        if not initial_external and (original.get("status") != "replacement_selected"
                or decisions.get("eaf_input_interpretation") != "maf"
                or decisions.get("eaf_is_maf") is not True):
            raise StudyAlleleFrequencyError("EAF replacement validation requires the original whole-study MAF decision.")
        paths, rows = chromosome_partitions.paths, chromosome_partitions.rows
        total = decisions.get("eaf_external_study_variant_count") if initial_external else original.get("study_variant_count")
        if (isinstance(total, bool) or not isinstance(total, int) or total <= 0
                or set(paths) != set(rows) or set(paths) != set(resource_maps)
                or sum(rows.values()) != total):
            raise StudyAlleleFrequencyError(label + " partitions do not reconcile to the fixed original study variant count.")
        try:
            _require_initial_frequency_inputs(decisions, verify_partitions=not initial_external)
        except StudyAlleleFrequencyError:
            decisions[evidence_key].update(status="invalid", failure_reason="replacement_inputs_changed",
                decision_reason="replacement_inputs_changed")
            raise
        column = sample_columns.get("eaf_col") or (decisions.get("eaf_initial_frequency") or {}).get("column")
        keys = {key: sample_columns[key + "_col"] for key in ("chr", "pos", "ea", "oa")}
        projection = list(dict.fromkeys([SOURCE_INPUT_ROW_COLUMN, *keys.values()] +
            ([column, PREPARED_INPUT_COLUMN, PREPARED_MATCHES_COLUMN] if not initial_external else [])))
        configuration = _replacement_configuration(sample_columns, decisions, external_eaf_colmap, pol)
        raw_screens, effective_screens = {}, {}
        initial_column = (unused_column_name("__postgwas_initial_frequency",
            *(pl.scan_parquet(path) for path in paths.values())) if initial_external else None)
        external_exclusion = "external_unresolved_rows" if initial_external else "replacement_unresolved_rows"
        bindings, reused, exclusions = {}, [], {"original_orientation_ineligible_rows": 0,
            external_exclusion: 0}
        if initial_external:
            # These small observed counters remain explicitly per chromosome
            # until the complete whole-study preparation reaches its screen.
            decisions[evidence_key].update(raw_frequency_screens_by_chromosome=raw_screens,
                effective_frequency_screens_by_chromosome=effective_screens,
                excluded_counts=exclusions)
        # The spool is an internal temporary protocol, not a scientific output.
        # Native temporary names avoid a second output-layout filename policy.
        with TemporaryDirectory(dir=Path(next(iter(paths.values()))).parent) as temporary:
            paired_paths, cache_paths = [], {}
            for chromosome, path in paths.items():
                prepared = decisions.get("eaf_prepared_references", {}).get(chromosome)
                original_partition = prepared_work_fingerprint(path, error_type=StudyAlleleFrequencyError)
                if initial_external:
                    frame = pl.read_parquet(original_partition["path"], columns=projection)
                    if any(name in pl.scan_parquet(path).schema for name in (PREPARED_INPUT_COLUMN, PREPARED_MATCHES_COLUMN)):
                        raise StudyAlleleFrequencyError("Initial external AF input contains reserved preparation columns.")
                    _, packed, prepared = prepare_study_reference_matches(frame, sample_columns, None,
                        resource_maps[chromosome]["default_eaf_file"],
                        resource_maps[chromosome]["default_eaf_reference_column"],
                        default_eaf_colmap, pol, warn=emit.warn, diagnostic_require_biallelic=False)
                    frame = frame.join(packed, on=SOURCE_INPUT_ROW_COLUMN, how="left", coalesce=True)
                    decisions.setdefault("eaf_prepared_references", {})[chromosome] = prepared
                    if prepared["contract"]["reference"] != source_bindings[chromosome]["reference_identity"]:
                        raise StudyAlleleFrequencyError(label + " default reference changed during preparation; rerun dataset preparation.")
                if not prepared:
                    raise StudyAlleleFrequencyError((label + " lacks prepared original orientation for chromosome %s.") % chromosome)
                reference = prepared["contract"]
                source = replacement.get("resolved_files", {}).get(chromosome)
                if source is None or optional_text(replacement.get("column")) is None:
                    raise StudyAlleleFrequencyError(label + (" lacks the resolved external file or column." if initial_external else " lacks the resolved replacement file or column."))
                resources = resource_maps[chromosome]
                same_source = same_frequency_reference_source(source, replacement["column"], external_eaf_colmap,
                    reference["reference"]["path"], reference["population"], reference["mapping"])
                staged = resources.get("user_eaf_file") or (None if same_source else source)
                binding = {"partition": original_partition,
                    "partition_path": str(Path(path).expanduser().absolute()), "source": (source_bindings[chromosome]["identity"] if initial_external else _replacement_file_identity(source)),
                    "source_path": str(Path(replacement.get("declared_files", {}).get(chromosome, source)).expanduser().absolute()),
                    "staged": ((staged_bindings or {}).get(chromosome)
                        or _replacement_file_identity(staged, label=label + " input")) if staged else None,
                    "staged_path": str(Path(staged).expanduser().absolute()) if staged else None,
                    "reference_path": str(Path(resources["default_eaf_file"]).expanduser().absolute()),
                    "prepared_reference": deepcopy(prepared)}
                if initial_external:
                    binding["source_mode"] = "external_eaf"
                _require_unchanged_replacement_binding(binding)
                if not initial_external:
                    frame = pl.read_parquet(binding["partition"]["path"], columns=projection)
                if frame.height != rows[chromosome]:
                    raise StudyAlleleFrequencyError((label + " partition row count changed for chromosome %s.") % chromosome)
                validate_prepared_reference_matches(frame, sample_columns, prepared,
                    resources["default_eaf_file"], resources["default_eaf_reference_column"], reference["mapping"], pol)
                matches = unpack_prepared_reference_matches(frame, SOURCE_INPUT_ROW_COLUMN)
                eligible = frequency_diagnostic_candidates(matches, SOURCE_INPUT_ROW_COLUMN, require_biallelic=False)
                exclusions["original_orientation_ineligible_rows"] += frame.height - eligible.height
                working, _ = canonicalize_variant_frame(frame.drop([PREPARED_INPUT_COLUMN, PREPARED_MATCHES_COLUMN]),
                    keys, pol, error_type=StudyAlleleFrequencyError, warn=emit.warn,
                    label="replacement frequency diagnostic keys")
                working = working.join(eligible, on=SOURCE_INPUT_ROW_COLUMN, how="inner").with_columns(
                    pl.col("__strand_alt").alias(keys["ea"]), pl.col("__strand_ref").alias(keys["oa"]))
                if same_source:
                    effect_view = _effect_frequency_view(frame, matches, keys, pol)
                    reused.append(chromosome)
                    if initial_external:
                        raw_screens[chromosome] = prepared["duplicate_statistics"]["raw_frequency_screen"]
                        with NamedTemporaryFile(dir=temporary, delete=False) as handle:
                            cache_path = Path(handle.name)
                        packed.join(effect_view.select(SOURCE_INPUT_ROW_COLUMN,
                            pl.col("effect_frequency").alias(initial_column)), on=SOURCE_INPUT_ROW_COLUMN,
                            how="left", coalesce=True).write_parquet(cache_path,
                            compression=str(pol.get("input.chromosome_partition_compression")))
                        cache_paths[chromosome] = cache_path
                else:
                    candidates = _replacement_candidate_keys(frame, matches, keys, pol,
                        label="external AF assignment keys" if initial_external else "replacement assignment keys")
                    assignments, assignment, statistics = merge_external_allele_frequencies(
                        candidates, staged, replacement["column"], external_eaf_colmap, keys,
                        chromosome=chromosome, policies=pol, logger=logger, min_match_fraction=0.0,
                        study_columns_canonical=True, declared_alt_frequency=not initial_external,
                        retain_biallelic_eligibility=True, _prepare_only=True,
                        _screen_chromosomes=[chromosome] if initial_external else None)
                    if initial_external:
                        raw_screens[chromosome] = statistics.pop("raw_frequency_screen")
                    assignments, assignment = _reserve_assignment_columns(assignments, assignment, working)
                    raw_names = [assignment["value"], assignment["orientation"],
                        *assignment["frequency_columns"], *assignment["site_columns"]]
                    assignments = assignments.select(SOURCE_INPUT_ROW_COLUMN, *keys.values(), *raw_names)
                    effect_view = _effect_frequency_view(frame, matches, keys, pol, assignments, assignment)
                    cache_column = unused_column_name("__postgwas_replacement_assignments", pl.scan_parquet(binding["partition"]["path"]))
                    cache = assignments.group_by(SOURCE_INPUT_ROW_COLUMN).agg(
                        pl.struct([*keys.values(), *raw_names]).alias(cache_column))
                    if initial_external:
                        cache = packed.join(cache, on=SOURCE_INPUT_ROW_COLUMN, how="left", coalesce=True)
                        cache = cache.join(effect_view.select(SOURCE_INPUT_ROW_COLUMN,
                            pl.col("effect_frequency").alias(initial_column)), on=SOURCE_INPUT_ROW_COLUMN,
                            how="left", coalesce=True)
                    if cache.height != frame.height or cache[cache_column].null_count():
                        raise StudyAlleleFrequencyError(("External AF" if initial_external else "Replacement") + " assignment cache did not cover every source row.")
                    with NamedTemporaryFile(dir=temporary, delete=False) as handle:
                        cache_path = Path(handle.name)
                    cache.write_parquet(cache_path, compression=str(pol.get("input.chromosome_partition_compression")))
                    cache_paths[chromosome] = cache_path
                    binding["assignment_cache"] = {"column": cache_column, "assignment": assignment, "statistics": statistics}
                    if initial_external:
                        binding["assignment_cache"]["source_mode"] = "external_eaf"
                    binding["assignment_cache_digest"] = configuration_digest(binding["assignment_cache"])
                    del candidates, assignments, cache
                if initial_external:
                    # The proof has reached its private spool/cache. Do not keep
                    # the preceding shard alive while preparing the next one.
                    del packed
                effective_screens[chromosome] = frequency_maf_screen(effect_view, "effect_frequency",
                    float(pol.get("eaf.maf_decision_cutoff")))
                effect_name = unused_column_name("__postgwas_effect_frequency", working)
                default_name = unused_column_name("__postgwas_default_frequency", working, reserved=(effect_name,))
                working = working.join(effect_view.rename({"effect_frequency": effect_name,
                    "default_frequency": default_name}), on=SOURCE_INPUT_ROW_COLUMN, how="left", coalesce=True)
                exclusions[external_exclusion] += working.filter(~valid_frequency_mask(pl.col(effect_name))).height
                # Both selected A and reference f stay in the ORIGINAL GWAS EA
                # direction. Folding happens only in the pooled comparator.
                if initial_external:
                    paired = working.select(
                        pl.col(effect_name).alias("study_frequency"),
                        pl.col(default_name).alias("replacement_frequency"),
                        pl.lit(True).alias("eligible"))
                else:
                    paired = working.select(
                        pl.col(column).cast(pl.Float64, strict=False).alias("study_frequency"),
                        pl.col(effect_name).alias("replacement_frequency"),
                        pl.lit(True).alias("eligible"))
                with NamedTemporaryFile(dir=temporary, delete=False) as handle:
                    # Retain this native unique name for the private spool's
                    # lifetime, independently of configured partition basenames.
                    paired_path = Path(handle.name)
                paired.write_parquet(paired_path, compression=str(pol.get("input.chromosome_partition_compression")))
                paired_paths.append(paired_path)
                bindings[chromosome] = binding
                del frame, matches, eligible, working, paired, effect_view
            if initial_external:
                evidence = _validate_initial_external_pairs(paired_paths, raw_screens, bindings,
                    replacement, external_eaf_colmap, exclusions, total, pol, emit, decisions[evidence_key],
                    effective_screens=effective_screens)
            else:
                screen = _effective_frequency_screen(effective_screens.values(), pol)
                source_evidence = {**replacement, "mapping": dict(external_eaf_colmap),
                    "identities": {c: b["source"] for c, b in bindings.items()},
                    "reused_default_chromosomes": reused, "same_source_evidence_reused": bool(reused)}
                decisions[evidence_key] = {"status": "inconclusive", "effective_frequency_screen": screen,
                    "study_variant_count": total, "excluded_counts": exclusions, "source": source_evidence}
                if not screen["usable"] or screen["maf_like"]:
                    reason = "replacement_frequency_maf_like" if screen["maf_like"] else "no_usable_replacement_frequencies"
                    decisions[evidence_key].update(failure_reason=reason, decision_reason=reason,
                        failed_requirements=[reason])
                    raise StudyAlleleFrequencyError(
                        "EAF replacement rejected: %s. The merged GWAS effect-allele frequency failed its MAF-like screen. "
                        "Correct --eaf-replacement-file / --eaf-replacement-column; no comparison, fallback or retry will follow."
                        % reason, study_decisions=decisions)
                evidence = validate_replacement_frequency_pairs(
                    lambda: (pl.read_parquet(path) for path in paired_paths), "study_frequency", "replacement_frequency",
                policies=pol, logger=logger, study_variant_count=total, eligibility_column="eligible",
                excluded_counts=exclusions, source=source_evidence)
                evidence["effective_frequency_screen"] = screen
            decisions[evidence_key] = evidence
            if initial_external and evidence.get("maf_reference_decision") == "maf":
                if configuration != _replacement_configuration(sample_columns, decisions, external_eaf_colmap, pol):
                    raise StudyAlleleFrequencyError("AF validation configuration changed during preparation; rerun dataset preparation.")
                _select_frequency_replacement(decisions, evidence)
                inputs = {}
                for binding in bindings.values():
                    for pathname, identity in ((binding["source_path"], binding["source"]),
                            (binding["reference_path"], binding["prepared_reference"]["contract"]["reference"])):
                        inputs[pathname] = {"path": pathname, "identity": identity}
                decisions["eaf_initial_frequency"] = {"source_kind": "external_eaf", "column": initial_column,
                    "source": deepcopy(replacement), "inputs": list(inputs.values())}
                decisions.update(eaf_is_maf_initial=True, eaf_is_maf_initial_source="merged_external_frequency_screen",
                    eaf_is_maf_source="dataset_reference_maf_support")
                configuration = _replacement_configuration(sample_columns, decisions, external_eaf_colmap, pol)
            try:
                if configuration != _replacement_configuration(sample_columns, decisions, external_eaf_colmap, pol):
                    raise StudyAlleleFrequencyError("AF validation configuration changed during preparation; rerun dataset preparation.")
                for binding in bindings.values():
                    _require_unchanged_replacement_binding(binding)
            except Exception:
                evidence.update(status="invalid", failure_reason="external_inputs_changed" if initial_external else "replacement_inputs_changed")
                raise
            if evidence["status"] in ("passed", "replacement_selected"):
                try:
                    for chromosome, cache_path in cache_paths.items():
                        binding = bindings[chromosome]
                        _require_unchanged_replacement_binding(binding)
                        # A narrow in-memory build side and a streaming wide shard:
                        # never collect the full study simply to append private data.
                        cache = pl.read_parquet(cache_path)
                        if initial_external:
                            if evidence["status"] == "passed":
                                cache = cache.drop(initial_column)
                            elif binding.get("assignment_cache"):
                                cache = cache.drop(binding["assignment_cache"]["column"])
                        physical_path = binding["partition"]["path"]
                        joined = pl.scan_parquet(physical_path).join(cache.lazy(),
                            on=SOURCE_INPUT_ROW_COLUMN, how="left", coalesce=True)
                        _write_parquet_atomically(joined, Path(physical_path),
                            compression=str(pol.get("input.chromosome_partition_compression")),
                            expected_rows=rows[chromosome],
                            validate_before_replace=lambda binding=binding: _require_unchanged_replacement_binding(binding))
                        binding["partition"] = prepared_work_fingerprint(physical_path, error_type=StudyAlleleFrequencyError)
                        del cache, joined
                    for binding in bindings.values():
                        _require_unchanged_replacement_binding(binding)
                    if configuration != _replacement_configuration(sample_columns, decisions, external_eaf_colmap, pol):
                        raise StudyAlleleFrequencyError("AF validation configuration changed during publication; rerun dataset preparation.")
                except Exception:
                    evidence.update(status="invalid", failure_reason="external_cache_publication_failed" if initial_external else "replacement_cache_publication_failed")
                    raise
        if initial_external and evidence["status"] == "replacement_selected":
            decisions["eaf_initial_frequency"]["partitions"] = {
                chromosome: deepcopy(binding["partition"]) for chromosome, binding in bindings.items()}
            return decisions
        if initial_external and evidence["status"] != "passed":
            raise StudyAlleleFrequencyError(
                "Whole-study external AF validation stopped: %s. No chromosome analysis was started; "
                "no automatic replacement or dataset retry will be attempted." % evidence["decision_reason"],
                study_decisions=decisions)
        if evidence["status"] != "passed":
            raise StudyAlleleFrequencyError(
                "EAF replacement validation is inconclusive: %s. No chromosome analysis was started. "
                "Correct --eaf-replacement-file / --eaf-replacement-column or their declared allele mapping; "
                "no fallback or dataset retry will be attempted." % ", ".join(evidence["failed_requirements"]),
                study_decisions=decisions)
        # An approval contract is issued only for successful, unchanged inputs.
        evidence["contract"] = {"configuration": configuration, "chromosomes": bindings}
        return decisions
    except Exception as exc:
        # Invalid source/I/O errors retain their type; all failures retain the
        # original decision for canonical dataset failure finalization.
        evidence = decisions[evidence_key]
        if evidence.get("status") != "inconclusive":
            mode = "external" if initial_external else "replacement"
            reason = evidence.get("failure_reason") or (
                mode + "_preparation_failed" if isinstance(exc, StudyAlleleFrequencyError) else mode + "_inputs_invalid")
            evidence.update(status="invalid", failure_reason=reason, decision_reason=reason)
        evidence.pop("contract", None)
        evidence["failure_message"] = str(exc)
        if evidence.get("status") == "invalid":
            emit.decide(label, evidence, "invalid")
        if initial_external and not isinstance(exc, StudyAlleleFrequencyError):
            raise StudyAlleleFrequencyError(label + ": " + str(exc), study_decisions=decisions) from exc
        exc.study_decisions = deepcopy(decisions)
        raise


def _prepared_frequency_rows_equal(left, right, list_columns):
    """Compare row-bound frequency evidence exactly without nested row objects.

    Prepared matches and assignments are List<Struct<primitive>> protocols.
    Outer validity and per-row lengths preserve null/empty distinctions and all
    list boundaries; exploded parent validity and every struct field preserve
    entry order and values. This is lossless comparison, not a hash or tolerance.
    """
    if left.shape != right.shape or left.columns != right.columns:
        return False
    for name in list_columns:
        dtype = left.schema.get(name)
        if name not in left.columns:
            return False
        if (dtype != right.schema[name] or not isinstance(dtype, pl.List)
                or not isinstance(dtype.inner, pl.Struct)):
            # Preserve native equality semantics for other schema combinations,
            # including compatible numeric widths; the fast path is protocol-specific.
            if not left.select(name).equals(right.select(name), null_equal=True):
                return False
            continue
        lhs, rhs = left[name], right[name]
        # Polars 0.20 reports length zero for BOTH null and empty lists.
        if (not lhs.is_null().equals(rhs.is_null(), null_equal=True)
                or not lhs.list.len().equals(rhs.list.len(), null_equal=True)):
            return False
        lhs, rhs = lhs.explode(), rhs.explode()
        if not lhs.is_null().equals(rhs.is_null(), null_equal=True):
            return False
        if not lhs.struct.unnest().equals(rhs.struct.unnest(), null_equal=True):
            return False
    return left.drop(list_columns).equals(right.drop(list_columns), null_equal=True)


def _require_replacement_validation(chromosome, frame, columns, decisions,
                                    external_file, mapping, policies):
    """Allow only unchanged validated input rows, including legitimate subsets."""
    initial_external = decisions.get("eaf_selected_source") == "external_eaf"
    evidence = decisions.get("eaf_external_validation" if initial_external else "eaf_replacement_validation") or {}
    label = _frequency_validation_label(initial_external)
    _require_initial_frequency_inputs(decisions)
    contract = evidence.get("contract") or {}
    binding = contract.get("chromosomes", {}).get(str(chromosome))
    if evidence.get("status") != "passed" or not binding:
        raise AlleleFrequencyError(label + " requires passed, row-bound approval." if initial_external else
            "EAF replacement adoption requires passed, row-bound replacement validation.")
    if contract.get("configuration") != _replacement_configuration(columns, decisions, mapping, policies):
        raise AlleleFrequencyError(label + " configuration changed; rerun dataset preparation.")
    source = decisions["eaf_external_source" if initial_external else "eaf_replacement"]["resolved_files"].get(str(chromosome))
    current_staged = optional_text(external_file) or source
    expected_attachment = binding.get("staged") or binding["source"]
    _require_unchanged_replacement_binding(binding)
    if (expected_attachment != _replacement_file_identity(current_staged, label=label + " input")
            or binding["prepared_reference"] != decisions.get("eaf_prepared_references", {}).get(str(chromosome))):
        raise AlleleFrequencyError(label + " inputs changed; rerun dataset preparation.")
    _validate_source_row_ids(frame, label)
    projection = [SOURCE_INPUT_ROW_COLUMN, *[columns[key + "_col"] for key in ("chr", "pos", "ea", "oa")],
        *([columns["eaf_col"]] if columns.get("eaf_col") is not None else []), PREPARED_INPUT_COLUMN, PREPARED_MATCHES_COLUMN]
    initial_column = (decisions.get("eaf_initial_frequency") or {}).get("column")
    if initial_column is not None and initial_column not in projection:
        projection.append(initial_column)
    list_columns = [PREPARED_MATCHES_COLUMN]
    descriptor = binding.get("assignment_cache")
    if descriptor is None and binding.get("assignment_cache_digest") is not None:
        raise AlleleFrequencyError(label + " cache descriptor is missing; rerun dataset preparation.")
    if descriptor is not None:
        if configuration_digest(descriptor) != binding.get("assignment_cache_digest"):
            raise AlleleFrequencyError(label + " cache descriptor changed; rerun dataset preparation.")
        projection.append(descriptor["column"])
        list_columns.append(descriptor["column"])
    if any(name not in frame.columns for name in projection):
        raise AlleleFrequencyError(label + " input columns are incomplete; rerun dataset preparation.")
    expected = pl.read_parquet(binding["partition"]["path"], columns=projection)
    observed = frame.select(projection)
    if observed[SOURCE_INPUT_ROW_COLUMN].equals(expected[SOURCE_INPUT_ROW_COLUMN], check_dtypes=True):
        selected = expected
    else:
        selected = expected.join(observed.select(SOURCE_INPUT_ROW_COLUMN),
            on=SOURCE_INPUT_ROW_COLUMN, how="semi").sort(SOURCE_INPUT_ROW_COLUMN)
        observed = observed.sort(SOURCE_INPUT_ROW_COLUMN)
    if not _prepared_frequency_rows_equal(observed, selected, list_columns):
        raise AlleleFrequencyError(label + " rows changed or contain unknown source IDs; rerun dataset preparation.")
    # Attach through the verified resolved pathname, not a mutable alias that
    # could be retargeted while strand work runs between this guard and the read.
    return expected_attachment["path"]


def _study_decision_flag(study_decision) -> Optional[bool]:
    """Read 'is the study frequency column a MAF?' out of a study-level decision.

    Accepts a plain bool, or a dict carrying one of the documented keys.  Returns
    ``None`` when the decision says nothing about the frequency column.
    """
    if study_decision is None:
        return None
    if isinstance(study_decision, bool):
        return study_decision
    if isinstance(study_decision, dict):
        for key in ("eaf_is_maf", "is_maf", "eaf_is_minor_allele_frequency", "maf_like"):
            if key in study_decision and study_decision[key] is not None:
                return bool(study_decision[key])
        return None
    raise AlleleFrequencyError(
        "study_decision must be a bool or a dict such as {'eaf_is_maf': False}; got %r" % type(study_decision)
    )


# -------------------------------------------------------
# 5. Helper: Compute detailed AF/EAF statistics
# -------------------------------------------------------
def summarise_allele_frequency(
    df: pl.DataFrame,
    colname: str,
    clip_tolerance: Optional[float] = None,
    policies=None,
) -> Dict[str, Any]:
    """
    Comprehensive descriptive statistics for AF/EAF-like columns.
    Works on the DataFrame exactly as passed in.

    ``clip_tolerance`` (policy ``eaf.clip_tolerance``) is retained as a
    diagnostic boundary for overshoots. It never makes a value above 1 valid.
    """
    if colname not in df.columns:
        return {
            "column_present": False,
            "total_rows": df.height,
        }
    tolerance = float(
        clip_tolerance
        if clip_tolerance is not None
        else resolve_policies(policies).get("eaf.clip_tolerance")
    )
    total_rows = df.height
    if total_rows == 0:
        return {
            "column_present": True,
            "total_rows": 0,
            "non_null_rows": 0,
            "null_rows": 0,
            "clip_tolerance": tolerance,
        }
    value = pl.col(colname).cast(pl.Float64, strict=False)
    usable = valid_frequency_mask(value)
    stats_row = df.select([
        pl.len().alias("total_rows"),
        value.is_not_null().sum().alias("non_null_rows"),
        value.is_null().sum().alias("null_rows"),
        usable.sum().alias("usable_rows"),
        (value < 0.0).sum().alias("lt_0_rows"),
        (value > 1.0).sum().alias("gt_1_rows"),
        (value > tolerance).sum().alias("gt_1_05_rows"),
        (
            value.is_not_null() & ~value.is_finite()
        ).sum().alias("non_finite_rows"),
        (value < 0.0).sum().alias("invalid_negative_rows"),
        (
            value.is_not_null()
            & (
                ~value.is_finite()
                | (value < 0.0)
                | (value > 1.0)
            )
        ).sum().alias("invalid_total_rows"),
        (usable & (value == 0.0)).sum().alias("eq_0_rows"),
        (usable & (value == 1.0)).sum().alias("eq_1_rows"),
        (usable & (value <= 0.5)).sum().alias("le_0_5_rows"),
        pl.when(usable).then(value).min().alias("min"),
        pl.when(usable).then(value).max().alias("max"),
        pl.when(usable).then(value).mean().alias("mean"),
        pl.when(usable).then(value).median().alias("median"),
        pl.when(usable).then(value).quantile(0.25).alias("q25"),
        pl.when(usable).then(value).quantile(0.75).alias("q75"),
    ]).to_dicts()[0]
    usable_rows = stats_row["usable_rows"] or 0
    total_rows = stats_row["total_rows"] or 0

    def frac(n, d):
        return float(n) / float(d) if d else 0.0

    stats_row["null_fraction"] = frac(stats_row["null_rows"], total_rows)
    stats_row["invalid_total_fraction"] = frac(stats_row["invalid_total_rows"], total_rows)
    stats_row["unusable_rows"] = total_rows - usable_rows
    stats_row["usable_fraction"] = frac(usable_rows, total_rows)
    stats_row["le_0_5_fraction_among_total"] = frac(stats_row["le_0_5_rows"], total_rows)
    stats_row["le_0_5_fraction_among_usable"] = frac(stats_row["le_0_5_rows"], usable_rows)
    stats_row["degenerate_rows"] = (stats_row["eq_0_rows"] or 0) + (stats_row["eq_1_rows"] or 0)
    stats_row["clip_tolerance"] = tolerance
    stats_row["column_present"] = True
    return stats_row


def _summarise_frequency_after_changes(
    df: pl.DataFrame,
    colname: str,
    *,
    clip_tolerance: float,
    previous_values: Optional[pl.Series],
    previous_statistics: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """Reuse local scalar facts only for the same ordered numeric state."""
    if (
        previous_values is not None
        and previous_statistics is not None
        and previous_values.name == colname
        and previous_statistics.get("clip_tolerance") == clip_tolerance
        # Numeric equality treats +0 and -0 as equal, but descriptive extrema
        # can retain their different sign bits. Recalculate that boundary case.
        and not previous_statistics.get("eq_0_rows")
        and colname in df.columns
    ):
        current_values = df.get_column(colname)
        if (
            current_values.dtype == previous_values.dtype
            and current_values.chunk_lengths() == previous_values.chunk_lengths()
            and current_values.equals(
                previous_values, check_dtypes=True, null_equal=True,
            )
        ):
            # Before/after dictionaries are independently mutable public data.
            return dict(previous_statistics)
    return summarise_allele_frequency(
        df, colname, clip_tolerance=clip_tolerance,
    )


# -------------------------------------------------------
# 6. Helper: Validate + clean EAF stats
# -------------------------------------------------------
def validate_allele_frequency(
    df: pl.DataFrame,
    colname: str,
    maf_cutoff: Optional[float] = None,
    invalid_fraction_cutoff: Optional[float] = None,
    missing_fraction_cutoff: Optional[float] = None,
    policies=None,
    logger=None,
    ctx=None,
    rejects=None,
    out_of_range: Optional[str] = None,
    clip_tolerance: Optional[float] = None,
    study_eaf_accepted: bool = False,
) -> Tuple[bool, bool, int, int, pl.DataFrame, float, Dict[str, Any]]:
    """
    Validate and clean AF/EAF column.

    ``study_eaf_accepted`` reuses the whole-study decision without repeating
    MAF classification or its warning. Descriptive, range and missing QC remain.
    Omitted thresholds use the resolved EAF policies; explicit arguments take
    precedence, including zero-valued thresholds.

    Out-of-range values are counted *before* anything is modified, then handled
    according to ``eaf.out_of_range``:

    ``clip``   explicitly force them back into 0..1
    ``null``   blank them, keeping the variant
    ``reject`` remove the variant, reason ``eaf_out_of_range``
    ``fail``   raise :class:`AlleleFrequencyError`

    Returns
    -------
    (
        is_valid_range,
        is_suspicious_maf,
        out_of_range_count,
        missing_count,
        cleaned_df,
        low_freq_pct_after_cleaning,
        stats_dict
    )

    ``is_valid_range`` is reported honestly: it is ``False`` when the share of
    non-finite values or values outside the closed interval 0..1
    exceeds ``invalid_fraction_cutoff``.  ``missing_count`` is the null count of
    the column — it used to be a row-count difference across ``with_columns``,
    which never changes the row count and so was always 0.
    """
    emit = active_context(logger, ctx, df.height)
    resolved = resolve_policies(policies)
    maf_cutoff = (
        float(resolved.get("eaf.maf_decision_cutoff"))
        if maf_cutoff is None else maf_cutoff
    )
    invalid_fraction_cutoff = (
        float(resolved.get("eaf.invalid_fraction_cutoff"))
        if invalid_fraction_cutoff is None else invalid_fraction_cutoff
    )
    missing_fraction_cutoff = (
        float(resolved.get("eaf.missing_fraction_cutoff"))
        if missing_fraction_cutoff is None else missing_fraction_cutoff
    )
    tolerance = float(
        clip_tolerance
        if clip_tolerance is not None
        else resolved.get("eaf.clip_tolerance")
    )
    action = str(
        out_of_range
        if out_of_range is not None
        else resolved.get("eaf.out_of_range")
    )

    if colname not in df.columns:
        stats = {"column_present": False, "total_rows": df.height}
        return False, False, 0, df.height, df, 0.0, stats
    total_before = df.height
    if total_before == 0:
        stats = summarise_allele_frequency(df, colname, clip_tolerance=tolerance)
        return False, False, 0, 0, df, 0.0, stats

    initial_stats = summarise_allele_frequency(df, colname, clip_tolerance=tolerance)
    initial_values = df.get_column(colname).clone()

    # ``summarise_allele_frequency`` already calculated these values in one
    # aggregate pass over the untouched column; reuse them rather than scanning
    # the chromosome a second time before applying the configured action.
    out_of_range_count = int(initial_stats["invalid_total_rows"] or 0)
    outside_unit_count = out_of_range_count
    n_missing = int(initial_stats["null_rows"] or 0)

    out_of_range_fraction = out_of_range_count / total_before if total_before else 0.0
    is_valid_range = out_of_range_fraction <= invalid_fraction_cutoff

    if out_of_range_count and out_of_range_fraction > invalid_fraction_cutoff:
        emit.warn(
            "Column '%s' has %s non-finite values or values outside 0 to 1, which is %s of all %s variants - far more "
            "than the %s allowed by eaf.invalid_fraction_cutoff. The column may not be a frequency at all."
            % (
                colname,
                _fmt_int(out_of_range_count),
                _fmt_pct(out_of_range_fraction),
                _fmt_int(total_before),
                _fmt_pct(invalid_fraction_cutoff),
            )
        )

    # ---- apply the out-of-range policy ---------------------------------
    cleaned_df = df
    invalid_frequency = (
        pl.col(colname).is_not_null()
        & (
            ~pl.col(colname).is_finite()
            | (pl.col(colname) < 0.0)
            | (pl.col(colname) > 1.0)
        )
    ).fill_null(False)
    if out_of_range_count and action == "fail":
        raise AlleleFrequencyError(
            "Column '%s' has %s non-finite frequency values or values outside 0 to 1 and policy "
            "eaf.out_of_range is 'fail'." % (colname, _fmt_int(out_of_range_count))
        )
    if action == "reject":
        # fill_null(False): a null frequency is not out of range, it is missing,
        # and missingness is somebody else's decision
        cleaned_df, _ = reject_rows(
            cleaned_df,
            invalid_frequency,
            step_label=STEP_LABEL,
            reason="eaf_out_of_range",
            context=emit,
            collector=rejects,
            detail="%s non-finite or outside 0 to 1" % colname,
            check_name="frequency range",
            description=(
                "Non-finite frequencies and frequencies outside 0 to 1 cannot be probabilities; "
                "policy eaf.out_of_range is 'reject', so those variants were removed."
            ),
            warn_on_remove=True,
            warn_without_collector=True,
            record_empty=True,
        )
    elif action == "null":
        before = cleaned_df.height
        cleaned_df = cleaned_df.with_columns(
            pl.when(invalid_frequency)
            .then(None)
            .otherwise(pl.col(colname))
            .alias(colname)
        )
        emit.qc(
            "frequency range",
            "%s values are non-finite or outside 0 to 1; policy eaf.out_of_range is 'null', so those "
            "frequencies were blanked and the variants kept." % _fmt_int(out_of_range_count),
            before,
            cleaned_df.height,
            changed=out_of_range_count,
            warn=out_of_range_count > 0,
            step=STEP_LABEL,
        )
    elif action == "clip":  # explicit opt-in correction
        before = cleaned_df.height
        cleaned_df = cleaned_df.with_columns(
            pl.when(pl.col(colname).is_not_null() & pl.col(colname).is_finite())
            .then(pl.col(colname).clip(0.0, 1.0))
            .otherwise(None)
            .alias(colname)
        )
        emit.qc(
            "frequency range",
            "%s values are non-finite or outside 0 to 1; policy eaf.out_of_range is 'clip', so "
            "finite values were forced into range and non-finite values were blanked."
            % _fmt_int(outside_unit_count),
            before,
            before,
            changed=outside_unit_count,
            warn=outside_unit_count > 0,
            step=STEP_LABEL,
        )
    else:  # "fail" with no invalid values
        emit.qc(
            "frequency range",
            "No non-finite frequencies or values outside 0 to 1 were found; "
            "policy eaf.out_of_range is 'fail', so no action was needed.",
            cleaned_df.height,
            cleaned_df.height,
            changed=0,
            step=STEP_LABEL,
        )

    # ---- missingness (null_count, not a row-count difference) ----------
    final_stats = _summarise_frequency_after_changes(
        cleaned_df, colname, clip_tolerance=tolerance,
        previous_values=initial_values, previous_statistics=initial_stats,
    )
    final_missing = int(final_stats.get("null_rows") or 0)
    missing_fraction = final_missing / total_before if total_before else 0.0
    if missing_fraction > missing_fraction_cutoff:
        emit.warn(
            "Column '%s' is empty for %s of the %s variants (%s), more than the %s allowed by "
            "eaf.missing_fraction_cutoff."
            % (
                colname,
                _fmt_int(final_missing),
                _fmt_int(total_before),
                _fmt_pct(missing_fraction),
                _fmt_pct(missing_fraction_cutoff),
            )
        )

    # ---- MAF-likeness, measured over usable finite probabilities --------
    total_after = cleaned_df.height
    low_freq = int(final_stats.get("le_0_5_rows") or 0)
    low_freq_pct = float(final_stats.get("le_0_5_fraction_among_usable") or 0.0)
    usable_rows = int(final_stats.get("usable_rows") or 0)
    # Range/missing QC still runs, but accepted study AF is not reclassified.
    is_suspicious = not study_eaf_accepted and low_freq_pct > maf_cutoff
    if is_suspicious:
        emit.warn(
            "Column '%s': %s of its %s usable finite frequencies are at or below 0.5, more than the %s set by "
            "eaf.maf_decision_cutoff. This is a MAF-like distribution; EAF can also have this "
            "distribution when effect alleles are predominantly minor."
            % (colname, _fmt_pct(low_freq_pct), _fmt_int(usable_rows), _fmt_pct(maf_cutoff))
        )

    combined_stats = {
        "initial": initial_stats,
        "final": final_stats,
        "total_before_cleaning": total_before,
        "total_after_cleaning": total_after,
        "out_of_range_count": out_of_range_count,
        "out_of_range_fraction": out_of_range_fraction,
        "out_of_range_action": action,
        "outside_unit_interval_count": outside_unit_count,
        "is_valid_range": is_valid_range,
        # kept under the old key for the QC summary; it is now the honest null
        # count rather than a row-count difference that was always 0
        "initial_missing_count": n_missing,
        "missing_count_removed": final_missing,
        "missing_fraction_removed": missing_fraction,
        "missing_count": final_missing,
        "missing_fraction": missing_fraction,
        "low_freq_count_after_cleaning": low_freq,
        "low_freq_fraction_after_cleaning": low_freq_pct,
        "maf_like_flag": is_suspicious,
    }
    return is_valid_range, is_suspicious, out_of_range_count, final_missing, cleaned_df, low_freq_pct, combined_stats


# -------------------------------------------------------
# 7. Confirm a MAF-like frequency column against default EAF
# -------------------------------------------------------
def _measure_frequency_pairs(
    paired, study_column, reference_column, method, *, include_raw_maf=False, replacement_only=False,
    expected_pair_count=None,
):
    """Measure the same eligible pairs for all AF decisions.

    Original classifiers preserve the established Polars result. Replacement
    Pearson and replayable batches use centered cross-products and the parallel
    covariance identity
    (Chan et al., 1983, https://doi.org/10.1080/00031305.1983.10483115).
    This combines sufficient statistics, never correlations. Exact Spearman is
    evaluated on one global paired frame by the resource-checked caller.
    """
    study = pl.col(study_column).cast(pl.Float64, strict=False)
    reference = pl.col(reference_column).cast(pl.Float64, strict=False)
    study_maf = pl.min_horizontal(study, 1.0 - study)
    reference_maf = pl.min_horizontal(reference, 1.0 - reference)
    correlations = {
        "direct_eaf_correlation": (study, reference),
        "maf_correlation": (study_maf, reference_maf),
    }
    means = {
        "mean_absolute_maf_difference": (study_maf - reference_maf).abs(),
        "mean_eaf_error": (study - reference).abs(),
        "mean_maf_error": (study - reference_maf).abs(),
        "reference_effect_allele_minor_fraction": reference <= 0.5,
    }
    if include_raw_maf:
        correlations["supplied_to_reference_maf_correlation"] = (study, reference_maf)
        means["reference_effect_allele_major_fraction"] = reference > 0.5
    if replacement_only:
        correlations = {"supplied_to_reference_maf_correlation": (study, reference_maf)}
        means = {"mean_maf_error": (study - reference_maf).abs()}
        if isinstance(paired, pl.DataFrame) and method == "pearson":
            frame = paired
            paired = lambda: iter((frame,))
    if isinstance(paired, pl.DataFrame):
        expressions = {
            key: pl.when(((x.min() != x.max()) & (y.min() != y.max())).fill_null(False))
            .then(pl.corr(x, y, method=method)).otherwise(None)
            for key, (x, y) in correlations.items()
        }
        expressions.update({key: value.mean() for key, value in means.items()})
        measured = paired.select([value.alias(key) for key, value in expressions.items()]).row(0, named=True)
        return {key: float(value) if value is not None and math.isfinite(value) else None
                for key, value in measured.items()}
    if method != "pearson":
        raise StudyAlleleFrequencyError("Exact global Spearman requires a resource-checked paired frame.")
    # State per pair is n, mean(x), mean(y), centered SS(x), SS(y), cross-product.
    moments, ranges, sums, count = {}, {}, {key: 0.0 for key in means}, 0
    for batch in paired():
        n = batch.height
        if not n:
            continue
        for key, (x, y) in correlations.items():
            row = batch.select([
                x.mean().alias("mx"), y.mean().alias("my"),
                x.min().alias("xmin"), x.max().alias("xmax"),
                y.min().alias("ymin"), y.max().alias("ymax"),
                ((x - x.mean()) ** 2).sum().alias("xx"),
                ((y - y.mean()) ** 2).sum().alias("yy"),
                ((x - x.mean()) * (y - y.mean())).sum().alias("xy"),
            ]).row(0, named=True)
            old_range = ranges.get(key)
            ranges[key] = ((min(old_range[0], row["xmin"]), max(old_range[1], row["xmax"]),
                min(old_range[2], row["ymin"]), max(old_range[3], row["ymax"])) if old_range else
                (row["xmin"], row["xmax"], row["ymin"], row["ymax"]))
            if key not in moments:
                moments[key] = (n, row["mx"], row["my"], row["xx"], row["yy"], row["xy"])
                continue
            old_n, mx, my, xx, yy, xy = moments[key]
            total = old_n + n
            dx, dy = row["mx"] - mx, row["my"] - my
            weight = old_n * n / total
            moments[key] = (total, mx + dx * n / total, my + dy * n / total,
                math.fsum((xx, row["xx"], dx * dx * weight)),
                math.fsum((yy, row["yy"], dy * dy * weight)),
                math.fsum((xy, row["xy"], dx * dy * weight)))
        values = batch.select([value.sum().alias(key) for key, value in means.items()]).row(0, named=True)
        sums = {key: math.fsum((sums[key], float(values[key] or 0.0))) for key in means}
        count += n
    if expected_pair_count is not None and count != expected_pair_count:
        raise StudyAlleleFrequencyError("AF replay changed the pair count between counting and measurement.")
    result = {key: sums[key] / count if count else None for key in means}
    for key in correlations:
        moment = moments.get(key)
        bounds = ranges.get(key)
        value = (moment[5] / math.sqrt(moment[3]) / math.sqrt(moment[4])
                 if moment and moment[3] > 0.0 and moment[4] > 0.0
                 and bounds[0] < bounds[1] and bounds[2] < bounds[3] else None)
        # Pearson lies in [-1,1]; correct only final floating-point roundoff.
        result[key] = max(-1.0, min(1.0, value)) if value is not None and math.isfinite(value) else None
    return result


def _check_spearman_memory(n, pol, spearman_memory_budget_bytes=None, *, evidence=None):
    """Apply the shared conservative global-rank resource estimate before replay."""
    # Reuse the canonical conservative per-variant working estimate,
    # safety multiplier and parent/headroom reservation. This includes
    # ranks/sort workspace rather than budgeting only the two columns.
    estimate = n * float(pol.get("execution.memory_bytes_per_variant")) * float(pol.get("execution.memory_safety_factor"))
    budget = pol.get("execution.memory_budget_gb")
    available = None
    if budget is not None:
        import psutil
        budget_bytes = float(budget) * 1024**3  # canonical GiB units
        available = budget_bytes - max(psutil.Process().memory_info().rss,
            budget_bytes * float(pol.get("execution.memory_headroom_fraction")))
    if spearman_memory_budget_bytes is not None:
        if (isinstance(spearman_memory_budget_bytes, bool)
                or not math.isfinite(spearman_memory_budget_bytes) or spearman_memory_budget_bytes <= 0):
            raise StudyAlleleFrequencyError("Spearman memory budget must be positive finite bytes.")
        available = (float(spearman_memory_budget_bytes) if available is None
            else min(available, float(spearman_memory_budget_bytes)))
    resources = {"estimated_rank_memory_bytes": estimate, "available_rank_memory_bytes": available}
    if evidence is not None:
        evidence.update(resources)
    if available is not None and estimate > available:
        if evidence is not None:
            evidence["failure_reason"] = "spearman_memory_budget_exceeded"
        raise StudyAlleleFrequencyError(
            "Exact global Spearman AF comparison exceeds its memory budget: estimated %s bytes; available %s. "
            "Increase --memory-gb / execution.memory_budget_gb or adjust the canonical execution memory estimate; "
            "the correlation method was not changed." % (estimate, available))
    return resources


def validate_replacement_frequency_pairs(
    paired, study_column, replacement_column, policies=None, logger=None, ctx=None, *,
    study_variant_count, source=None, excluded_counts=None, eligibility_column=None,
    spearman_memory_budget_bytes=None,
):
    """Validate A against min(B, 1-B), keeping A and original-EA-aligned B intact.

    ``paired`` is a DataFrame or a replayable zero-argument batch factory over
    unchanged source rows (same values and row identity on every pass). The
    owning attachment/orientation code must prove
    unique non-palindromic biallelic eligibility, supplying only those rows or
    marking them in ``eligibility_column``. The fixed positive N includes every
    study row entering preflight. Numeric exclusions never change source rows.
    ``source`` retains the caller's immutable source identity and allele mapping;
    ``excluded_counts`` retains upstream diagnostic exclusions separately.

    Count/coverage precede measurement and ranking. Pearson is pooled using
    centered sufficient statistics. Spearman ranks all eligible pairs globally;
    its estimate uses canonical execution memory policies before a batch replay.
    An explicit byte budget may further constrain this estimate. Null canonical
    memory budget intentionally leaves direct library global ranks unbounded.
    Estimates are planning checks, not OS memory limits or capacity guarantees.
    """
    pol = resolve_policies(policies)
    emit = active_context(logger, ctx)
    if (isinstance(study_variant_count, bool) or not isinstance(study_variant_count, int)
            or study_variant_count <= 0):
        raise StudyAlleleFrequencyError("Replacement AF comparison requires a positive integer study variant count.")
    if not isinstance(paired, pl.DataFrame) and not callable(paired):
        raise StudyAlleleFrequencyError("Replacement AF pairs must be a DataFrame or replayable batch factory.")
    study = pl.col(study_column).cast(pl.Float64, strict=False)
    replacement = pl.col(replacement_column).cast(pl.Float64, strict=False)
    eligible = pl.col(eligibility_column).fill_null(False) if eligibility_column else pl.lit(True)
    missing = study.is_null() | replacement.is_null()
    finite = (study.is_finite() & replacement.is_finite()).fill_null(False)
    in_range = study.is_between(0.0, 1.0) & replacement.is_between(0.0, 1.0)
    usable = eligible & valid_frequency_mask(study) & valid_frequency_mask(replacement)
    counts = dict.fromkeys(("input_pair_rows", "comparable_variants", "major_count",
        "diagnostic_ineligible_rows", "missing_frequency_rows", "non_finite_frequency_rows",
        "out_of_range_frequency_rows"), 0)
    batches = (lambda: iter((paired,))) if isinstance(paired, pl.DataFrame) else paired
    for batch in batches():
        try:
            measured = batch.select([
                pl.len().alias("input_pair_rows"), usable.sum().alias("comparable_variants"),
                (usable & (replacement > 0.5)).sum().alias("major_count"),
                (~eligible).sum().alias("diagnostic_ineligible_rows"),
                (eligible & missing).sum().alias("missing_frequency_rows"),
                (eligible & ~missing & ~finite).sum().alias("non_finite_frequency_rows"),
                (eligible & ~missing & finite & ~in_range).sum().alias("out_of_range_frequency_rows"),
            ]).row(0, named=True)
        except pl.exceptions.ColumnNotFoundError as exc:
            raise StudyAlleleFrequencyError("Invalid replacement AF comparison column: %s" % exc) from exc
        counts = {key: counts[key] + int(measured[key] or 0) for key in counts}
    if counts["input_pair_rows"] > study_variant_count:
        raise StudyAlleleFrequencyError("Replacement AF comparison input pair count cannot exceed the study variant count.")
    n = counts["comparable_variants"]
    fraction = n / study_variant_count
    minimum = int(pol.get("eaf.maf_reference_min_overlap"))
    minimum_fraction = float(pol.get("eaf.maf_reference_min_overlap_fraction"))
    minimum_correlation = float(pol.get("eaf.study_reference_min_correlation"))
    maximum_difference = float(pol.get("eaf.maf_reference_max_mean_absolute_difference"))
    minimum_major = float(pol.get("eaf.study_reference_min_major_fraction"))
    method = str(pol.get("eaf.maf_reference_correlation_method"))
    failures = [reason for reason, failed in (
        ("insufficient_overlap", n < minimum),
        ("insufficient_overlap_fraction", fraction < minimum_fraction),
    ) if failed]
    evidence = {
        "source": deepcopy(source or {}), "correlation_method": method,
        "study_variant_count": study_variant_count, "comparable_variants": n,
        "input_pair_rows": counts["input_pair_rows"], "observed_overlap_fraction": fraction,
        "excluded_counts": {**deepcopy(excluded_counts or {}), **{key: counts[key] for key in (
            "diagnostic_ineligible_rows", "missing_frequency_rows", "non_finite_frequency_rows",
            "out_of_range_frequency_rows")}},
        "minimum_overlap": minimum, "minimum_overlap_fraction": minimum_fraction,
        "minimum_correlation": minimum_correlation, "maximum_mean_absolute_difference": maximum_difference,
        "minimum_major_fraction": minimum_major, "major_count": counts["major_count"],
        "major_fraction": counts["major_count"] / n if n else None,
        "correlation": None, "mean_absolute_difference": None,
        "shared_evidence_gate_passed": not failures,
    }
    if not failures:
        def usable_batches():
            for batch in batches():
                yield batch.filter(usable).select(study.alias(study_column), replacement.alias(replacement_column))
        selected = paired.filter(usable) if isinstance(paired, pl.DataFrame) else usable_batches
        if method == "spearman":
            evidence.update(_check_spearman_memory(n, pol, spearman_memory_budget_bytes))
            if not isinstance(selected, pl.DataFrame):
                frames = list(selected())
                selected = pl.concat(frames, rechunk=False) if frames else pl.DataFrame(
                    schema={study_column: pl.Float64, replacement_column: pl.Float64})
                if selected.height != n:
                    raise StudyAlleleFrequencyError("Replacement AF replay changed the pair count before global ranking.")
        metrics = _measure_frequency_pairs(selected, study_column, replacement_column, method,
            replacement_only=True, expected_pair_count=n)
        evidence.update(correlation=metrics["supplied_to_reference_maf_correlation"],
            mean_absolute_difference=metrics["mean_maf_error"])
        failures.extend(reason for reason, failed in (
            ("correlation_unavailable", evidence["correlation"] is None),
            ("correlation_below_minimum", evidence["correlation"] is not None and evidence["correlation"] < minimum_correlation),
            ("mean_absolute_difference_above_maximum", evidence["mean_absolute_difference"] is None
                or evidence["mean_absolute_difference"] > maximum_difference),
            ("major_fraction_not_above_minimum", evidence["major_fraction"] is None or evidence["major_fraction"] <= minimum_major),
        ) if failed)
    evidence.update(status="inconclusive" if failures else "passed", failed_requirements=failures,
        decision_reason=failures[0] if failures else "replacement_frequency_agreement")
    emit.decide("Whole-study replacement AF validation", evidence, evidence["status"])
    return evidence


def _compare_frequency_pairs(
    paired, study_column, reference_column, policies, emit, *, accept_direct_eaf=False,
    study_variant_count=None, pair_count=None, comparison_evidence=None,
):
    """Compare one eligible paired set once, retaining original EA direction.

    Folding as min(AF, 1-AF) is an allele-direction invariant for biallelic
    evidence. It cannot identify EAF. Whole-study callers evaluate both complete
    raw-source rules after the shared count/coverage gate, then use their truth
    table as the sole decision. Other callers retain the existing folded
    suitability and directional-margin rules.
    """
    n = paired.height if isinstance(paired, pl.DataFrame) else pair_count
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise StudyAlleleFrequencyError("Batched AF comparison requires a known nonnegative pair count.")
    minimum = int(policies.get("eaf.maf_reference_min_overlap"))
    minimum_fraction = None
    minimum_minor_fraction = None
    overlap_fraction = None
    if accept_direct_eaf:
        if (isinstance(study_variant_count, bool)
                or not isinstance(study_variant_count, int) or study_variant_count <= 0):
            raise StudyAlleleFrequencyError(
                "Whole-study AF comparison requires a positive integer study variant count."
            )
        if n > study_variant_count:
            raise StudyAlleleFrequencyError(
                "Whole-study AF comparison pair count cannot exceed the study variant count."
            )
        minimum_fraction = float(policies.get("eaf.maf_reference_min_overlap_fraction"))
        minimum_minor_fraction = float(policies.get("eaf.study_reference_min_minor_fraction"))
        overlap_fraction = n / study_variant_count
    method = str(policies.get("eaf.maf_reference_correlation_method"))
    minimum_correlation = float(policies.get("eaf.maf_reference_min_correlation"))
    maximum_difference = float(policies.get("eaf.maf_reference_max_mean_absolute_difference"))
    minor_cutoff = float(policies.get("eaf.reference_minor_fraction_cutoff"))
    margin = float(policies.get("eaf.maf_reference_error_margin"))
    shared_failures = []
    if accept_direct_eaf:
        shared_failures = [reason for reason, failed in (
            ("insufficient_overlap", n < minimum),
            ("insufficient_overlap_fraction", overlap_fraction < minimum_fraction),
        ) if failed]
    # Save established metadata before ranking/measurement can fail. Optional
    # caller-owned evidence does not change native/direct comparison semantics.
    stats = {
        "comparable_variants": n, "minimum_overlap": minimum,
        "maf_correlation_method": method, "minimum_maf_correlation": minimum_correlation,
        "maximum_mean_absolute_maf_difference": maximum_difference,
        "reference_minor_fraction_cutoff": minor_cutoff, "minimum_error_margin": margin,
    }
    if accept_direct_eaf:
        study_correlation = float(policies.get("eaf.study_reference_min_correlation"))
        major_minimum = float(policies.get("eaf.study_reference_min_major_fraction"))
        stats.update(study_variant_count=study_variant_count,
            observed_overlap_fraction=overlap_fraction, minimum_overlap_fraction=minimum_fraction,
            minimum_reference_effect_allele_minor_fraction=minimum_minor_fraction,
            minimum_study_correlation=study_correlation,
            minimum_reference_effect_allele_major_fraction=major_minimum,
            shared_evidence_gate_passed=not shared_failures,
            shared_evidence_failed_requirements=shared_failures,
            whole_study_rules_evaluated=False)
    measured = {key: None for key in (
        "direct_eaf_correlation", "maf_correlation", "mean_absolute_maf_difference",
        "mean_eaf_error", "mean_maf_error", "reference_effect_allele_minor_fraction",
        "supplied_to_reference_maf_correlation", "reference_effect_allele_major_fraction",
    )}
    if comparison_evidence is not None:
        comparison_evidence.update(stats, **measured)
    # Do not scan or rank when the count/coverage gate already failed.
    if not shared_failures:
        if method == "spearman" and not isinstance(paired, pl.DataFrame):
            stats.update(_check_spearman_memory(n, policies, evidence=comparison_evidence))
            paired = pl.concat(list(paired()), rechunk=False)
            if paired.height != n:
                raise StudyAlleleFrequencyError("AF replay changed the pair count before global ranking.")
        measured = _measure_frequency_pairs(paired, study_column, reference_column, method,
            include_raw_maf=accept_direct_eaf, expected_pair_count=n)
    stats.update({key: float(value) if value is not None and math.isfinite(value) else None
        for key, value in measured.items()})
    if comparison_evidence is not None:
        comparison_evidence.update(stats)
    eaf_error, maf_error = stats["mean_eaf_error"], stats["mean_maf_error"]
    stats["absolute_error_difference"] = (
        abs(eaf_error - maf_error) if eaf_error is not None and maf_error is not None else None
    )
    if accept_direct_eaf:
        # Whole-study thresholds are distinct from the ordinary external
        # classifier. Both raw-source rules share one measured paired set.
        stats["whole_study_rules_evaluated"] = not shared_failures
        for rule, correlation, error, profile, threshold, profile_failure in (
            (1, stats["direct_eaf_correlation"], eaf_error,
             stats["reference_effect_allele_minor_fraction"], minimum_minor_fraction,
             "reference_minor_fraction_not_above_minimum"),
            (2, stats["supplied_to_reference_maf_correlation"], maf_error,
             stats["reference_effect_allele_major_fraction"], major_minimum,
             "reference_major_fraction_not_above_minimum"),
        ):
            failures = list(shared_failures) if shared_failures else [
                reason for reason, failed in (
                    ("correlation_unavailable", correlation is None),
                    ("correlation_below_minimum", correlation is not None and correlation < study_correlation),
                    ("mean_absolute_difference_above_maximum", error is None or error > maximum_difference),
                    (profile_failure, profile is None or profile <= threshold),
                ) if failed
            ]
            stats["rule_%s_failed_requirements" % rule] = failures
            stats["rule_%s_passed" % rule] = not failures
        rule_one, rule_two = stats["rule_1_passed"], stats["rule_2_passed"]
        if shared_failures:
            decision, reason = "inconclusive", shared_failures[0]
        elif rule_one and not rule_two:
            decision, reason = "eaf", "allele_aligned_frequency_agreement"
        elif rule_two and not rule_one:
            decision, reason = "maf", "maf_interpretation_supported"
        else:
            decision, reason = "inconclusive", "both_rules_passed" if rule_one else "both_rules_failed"
        stats.update(
            # This established field means actual retention, including the
            # both-pass guard; rule_1_passed records the independent full rule.
            direct_eaf_accepted=decision == "eaf",
            direct_eaf_acceptance_reason=(stats["rule_1_failed_requirements"][0]
                if stats["rule_1_failed_requirements"] else "allele_aligned_frequency_agreement"),
            decision_reason=reason,
            # The ordinary folded-suitability classifier does not participate
            # here; its legacy metrics retain their meaning for audit only.
            reference_suitability_status="not_required",
            reference_suitability_reason="whole_study_rules_replace_folded_suitability",
        )
        return decision, stats

    # Ordinary external-source / chromosome classifier: preserve the existing
    # direct evidence, folded compatibility, profile veto and strict error gap.
    if n < minimum:
        direct_reason = "insufficient_overlap"
    elif stats["direct_eaf_correlation"] is None:
        direct_reason = "correlation_unavailable"
    elif stats["direct_eaf_correlation"] < minimum_correlation:
        direct_reason = "correlation_below_minimum"
    elif eaf_error is None or eaf_error > maximum_difference:
        direct_reason = "mean_absolute_difference_above_maximum"
    else:
        direct_reason = "allele_aligned_frequency_agreement"
    stats["direct_eaf_acceptance_reason"] = direct_reason
    stats["direct_eaf_accepted"] = direct_reason == "allele_aligned_frequency_agreement"
    if n < minimum:
        suitability_reason = "insufficient_overlap"
    elif stats["maf_correlation"] is None:
        suitability_reason = "maf_correlation_unavailable"
    elif stats["maf_correlation"] < minimum_correlation:
        suitability_reason = "maf_correlation_below_minimum"
    elif (stats["mean_absolute_maf_difference"] is None
          or stats["mean_absolute_maf_difference"] > maximum_difference):
        suitability_reason = "maf_difference_above_maximum"
    else:
        suitability_reason = "maf_correlation_and_difference_passed"
    suitable = suitability_reason == "maf_correlation_and_difference_passed"
    stats["reference_suitability_status"] = "suitable" if suitable else "inconclusive"
    stats["reference_suitability_reason"] = suitability_reason
    if not suitable:
        decision, reason = "inconclusive", suitability_reason
    elif stats["reference_effect_allele_minor_fraction"] > minor_cutoff:
        decision, reason = "inconclusive", "reference_effect_allele_mostly_minor"
    elif eaf_error + margin < maf_error:
        decision, reason = "eaf", "eaf_error_lower_by_required_margin"
    elif maf_error + margin < eaf_error:
        decision, reason = "maf", "maf_error_lower_by_required_margin"
    else:
        decision, reason = "inconclusive", "error_difference_below_required_margin"
    stats["decision_reason"] = reason
    emit.decide("Default EAF reference suitability", stats,
                "suitable for directional confirmation" if suitable else "inconclusive")
    return decision, stats


def confirm_eaf_with_reference(
    df: pl.DataFrame,
    input_af_col: str,
    reference_file: str,
    population_col: str,
    colmap: Dict[str, str],
    std_cols: Dict[str, str],
    policies=None,
    logger=None,
    ctx=None,
    aligned_reference_col: Optional[str] = None,
    study_columns_canonical: bool = False,
) -> Tuple[str, Dict[str, Any]]:
    """Return ``eaf``, ``maf`` or ``inconclusive`` for a MAF-like column.

    The reference population column is ALT frequency. Direct matches use it as
    supplied; allele-swapped matches use ``1-AF`` so every comparison is against
    the frequency of the study's listed effect allele. Palindromic SNPs are
    excluded because their strand cannot be inferred from the allele pair. Before
    testing direction, both frequencies are folded as ``min(AF, 1-AF)`` and the
    configured correlation plus mean-absolute-difference gates must show that
    the independent default panel meets the selected suitability criteria for
    these variants. Folding is a scientific invariant that removes allele
    direction; it cannot itself prove EAF. A reference set dominated by minor
    effect alleles is inconclusive because aligned EAF and folded MAF are then
    nearly identical; otherwise the two mean errors must differ by the
    configured minimum margin.
    """
    emit = active_context(logger, ctx, df.height)
    resolved = resolve_policies(policies)
    keys = [std_cols[name] for name in ("chr", "pos", "ea", "oa")]
    palindromic = palindromic_snp_expression(
        pl.col(std_cols["ea"]), pl.col(std_cols["oa"])
    )
    use_aligned = aligned_reference_col in df.columns if aligned_reference_col else False
    selected = keys + [input_af_col]
    if use_aligned:
        selected.append(aligned_reference_col)
        if STRAND_ACTION_COLUMN in df.columns:
            selected.append(STRAND_ACTION_COLUMN)
    study = df.select(selected).filter(
        valid_frequency_mask(pl.col(input_af_col))
        & ~palindromic
    ).rename({input_af_col: "_study_af"})
    if use_aligned:
        reference_col = "__maf_expected_study_eaf"
        if STRAND_ACTION_COLUMN in study.columns:
            swapped = pl.col(STRAND_ACTION_COLUMN).is_in([
                "forward_swapped", "reverse_complement_swapped"
            ])
            merged = study.with_columns(
                pl.when(swapped)
                .then(1.0 - pl.col(aligned_reference_col))
                .otherwise(pl.col(aligned_reference_col))
                .alias(reference_col)
            )
        else:
            merged = study.with_columns(
                pl.col(aligned_reference_col).alias(reference_col)
            )
        matched_rows = merged.filter(pl.col(reference_col).is_not_null())
        matched_count = matched_rows.height
        if STRAND_ACTION_COLUMN in matched_rows.columns:
            actions = matched_rows.select([
                pl.col(STRAND_ACTION_COLUMN).is_in([
                    "forward", "reverse_complement"
                ]).sum().alias("direct"),
                pl.col(STRAND_ACTION_COLUMN).is_in([
                    "forward_swapped", "reverse_complement_swapped"
                ]).sum().alias("swapped"),
            ]).to_dicts()[0]
        else:
            actions = {"direct": matched_count, "swapped": 0}
        merge_stats = {
            "direct_match_rows": int(actions["direct"] or 0),
            "flip_match_rows": int(actions["swapped"] or 0),
        }
        emit.info(
            "Reused the population AF already aligned during reference orientation; "
            "the chromosome reference was not read a second time."
        )
    else:
        merged, reference_col, merge_stats = merge_external_allele_frequencies(
            study,
            reference_file,
            population_col,
            colmap,
            std_cols,
            policies=policies,
            logger=logger,
            ctx=ctx,
            min_match_fraction=0.0,
            study_columns_canonical=study_columns_canonical,
        )
    matched = merged.filter(pl.col(reference_col).is_not_null())

    reference_value = pl.col(reference_col).cast(pl.Float64, strict=False)
    invalid_reference = int(matched.select(
        (
            ~reference_value.is_finite()
            | (reference_value < 0.0)
            | (reference_value > 1.0)
        ).sum().alias("invalid_reference")
    ).item() or 0)
    if invalid_reference:
        raise AlleleFrequencyError(
            "The default EAF table '%s' contains %s matched non-finite values or values outside "
            "0 to 1 in column '%s'."
            % (reference_file, _fmt_int(invalid_reference), population_col)
        )
    decision, stats = _compare_frequency_pairs(
        matched, "_study_af", reference_col, resolved, emit,
    )
    stats.update(
        direct_matches=merge_stats["direct_match_rows"],
        swapped_matches=merge_stats["flip_match_rows"],
    )
    emit.decide("MAF-like frequency column '%s'" % input_af_col, stats, decision)
    return decision, stats


# -------------------------------------------------------
# 9. Helper: Load Generic External
# ------------------------------------------------------
def merge_external_allele_frequencies(
    df: pl.DataFrame,
    path: str,
    eaf_col: str,
    colmap: dict,
    std_cols: dict,
    chromosome: str = "",
    policies=None,
    logger=None,
    ctx=None,
    rejects=None,
    min_match_fraction: Optional[float] = None,
    study_columns_canonical: bool = False,
    validation_reference_col: Optional[str] = None,
    validation_reference_file: Optional[str] = None,
    validation_reference_column: Optional[str] = None,
    declared_alt_frequency: bool = False,
    retain_biallelic_eligibility: bool = False,
    _prepare_only: bool = False,
    _screen_chromosomes=None,
) -> Tuple[pl.DataFrame, str | Dict[str, Any], Dict[str, Any]]:
    """Take the frequency from an external panel.

    Variants are matched in the orientations listed by ``eaf.match_orientations``
    (``direct`` then ``swap``); complement passes are deliberately not offered,
    see plan section 5.3. Unmatched variants keep a null frequency. The coverage
    gate counts only allele-key matches whose frequency is finite and between 0
    and 1; if that usable share falls below
    ``eaf.external_min_match_fraction`` the step fails instead of quietly
    reporting success. The raw external column is screened before swap-based
    alignment. When that screen is MAF-like and ``validation_reference_col``
    is available on the study frame, the existing EAF-versus-folded-MAF
    classifier validates the external column against that independent,
    allele-aligned reference before the caller adopts it.

    ``retain_biallelic_eligibility`` exposes the shared join's complete-site
    diagnostic flag in the returned frame; its collision-safe column name is
    returned as ``statistics['biallelic_eligibility_column']``. It does not
    filter attachment rows or change the ordinary coverage gate.

    Private ``_prepare_only`` returns raw assignments, their column descriptor
    (instead of the final column name), and source statistics. Initial-external
    and explicit-replacement preparation both use it; finalization and coverage
    run later on independently oriented worker rows. Ordinary calls retain the
    public tuple.
    """
    resolved_policies = resolve_policies(policies)
    emit = active_context(logger, ctx, df.height)
    min_match_fraction = float(
        min_match_fraction
        if min_match_fraction is not None
        else resolved_policies.get("eaf.external_min_match_fraction")
    )
    duplicate_exact_action = str(
        resolved_policies.get("external_reference.exact_duplicate_action")
    )
    duplicate_non_identical_action = str(
        resolved_policies.get(
            "external_reference.non_identical_duplicate_action"
        )
    )
    orientations = [
        str(value).lower()
        for value in resolved_policies.get("eaf.match_orientations")
    ]
    unsupported = [o for o in orientations if o not in ("direct", "swap")]
    if unsupported:
        raise AlleleFrequencyError(
            "eaf.match_orientations asks for %s, but only 'direct' and 'swap' are implemented. "
            "Study rows have already received a validated per-variant strand_action; the external "
            "frequency join must not perform a second complement pass or re-decide effect direction."
            % unsupported
        )
    stats: Dict[str, Any] = {
        "external_rows_loaded": 0,
        "input_rows_before_merge": df.height,
        "direct_match_rows": 0,
        "flip_match_rows": 0,
        "key_match_rows": 0,
        "key_match_fraction": 0.0,
        "usable_direct_match_rows": 0,
        "usable_flip_match_rows": 0,
        "usable_match_rows": 0,
        "usable_match_fraction": 0.0,
        "matched_missing_frequency_rows": 0,
        "matched_non_finite_frequency_rows": 0,
        "matched_out_of_range_frequency_rows": 0,
        "matched_unusable_frequency_rows": 0,
        "unusable_frequency_rows": 0,
        "unmatched_rows": 0,
        "total_missing": 0,
        "missing_fraction": 0.0,
        "match_fraction": 0.0,
        "min_match_fraction": min_match_fraction,
        "orientations": orientations,
        "palindromic_evaluated_rows": 0,
        "palindromic_orientation_resolved_rows": 0,
        "palindromic_orientation_rejected_rows": 0,
        "palindromic_orientation_state_counts": {},
    }
    # --------------------------------------------------
    # Load file
    # --------------------------------------------------
    load_path = path
    if not os.path.exists(load_path):
        raise AlleleFrequencyError(
            "The external effect-allele-frequency file was not found: '%s'. Check the 'eaffile' setting "
            "for this dataset." % load_path
        )
    ext_df, detected = read_reference_variant_table(
        load_path,
        colmap,
        resolved_policies,
        value_columns=[eaf_col],
        error_type=AlleleFrequencyError,
        description="external frequency file",
    )
    input_format = (
        detected.method
        if detected.kind == "unknown"
        else "%s separated" % detected.method
    )
    emit.info(
        "Read the external frequency file %s: %s rows, %d required columns, %s."
        % (
            os.path.basename(load_path),
            _fmt_int(ext_df.height),
            ext_df.width,
            input_format,
        )
    )
    stats["external_rows_loaded"] = ext_df.height
    source_column = eaf_col
    eaf_col = unused_column_name(source_column, df)
    stats.update(
        source_file=load_path,
        source_column=source_column,
        destination_column=eaf_col,
    )
    row_col = unused_column_name(
        "__postgwas_external_eaf_row", df, reserved=(eaf_col,),
    )
    study_for_join = df.with_row_index(row_col)
    assignment_reason_col = unused_column_name(
        "__postgwas_external_eaf_assignment_reason", study_for_join, ext_df,
        reserved=(eaf_col,),
    )
    biallelic_col = (
        unused_column_name(
            "__postgwas_external_eaf_biallelic", study_for_join, ext_df,
            reserved=(eaf_col, assignment_reason_col),
        )
        if retain_biallelic_eligibility or (not declared_alt_frequency and validation_reference_col is not None)
        else None
    )
    matched_all, orientation_col, join_stats = allele_oriented_left_join(
        study_for_join,
        ext_df,
        study_columns=std_cols,
        reference_columns={
            "chr": colmap["chr"],
            "pos": colmap["pos"],
            "ea": colmap["a1"],
            "oa": colmap["a2"],
        },
        policies=resolved_policies,
        value_column=source_column,
        output_column=eaf_col,
        duplicate_exact_action=duplicate_exact_action,
        duplicate_non_identical_action=duplicate_non_identical_action,
        orientations=orientations,
        swapped_value="one_minus",
        prefer_non_null_value=True,
        study_columns_canonical=study_columns_canonical,
        error_type=AlleleFrequencyError,
        warn=emit.warn,
        reference_label="external frequency file '%s'" % load_path,
        screen_frequency=not declared_alt_frequency,
        screen_chromosomes=_screen_chromosomes,
        require_complete_frequency_site=True,
        frequency_reason_column=assignment_reason_col,
        frequency_biallelic_column=biallelic_col,
        retain_assignment_metadata=True,
    )
    if join_stats["reference_duplicate_groups"]:
        emit.info(
            "External frequency reference duplicates: %s exact key/value groups "
            "(action '%s'); %s non-identical key groups covering %s rows "
            "(action '%s'); %s reference rows removed in total."
            % (
                _fmt_int(join_stats["reference_exact_duplicate_groups"]),
                duplicate_exact_action,
                _fmt_int(
                    join_stats["reference_non_identical_duplicate_groups"]
                ),
                _fmt_int(join_stats["reference_non_identical_duplicate_rows"]),
                duplicate_non_identical_action,
                _fmt_int(join_stats["reference_duplicate_rows_removed"]),
            )
        )
    stats.update(join_stats)

    assignment = stats.pop("assignment_metadata")
    assignment.update(row=row_col, reason=assignment_reason_col, biallelic=biallelic_col)
    if _prepare_only:
        return matched_all, assignment, stats
    return _finalize_external_allele_frequencies(
        matched_all, assignment, stats, std_cols, eaf_col, resolved_policies,
        chromosome=chromosome, emit=emit, logger=logger, ctx=ctx, rejects=rejects,
        validation_reference_col=validation_reference_col,
        validation_reference_file=validation_reference_file,
        validation_reference_column=validation_reference_column, colmap=colmap,
        declared_alt_frequency=declared_alt_frequency,
        retain_biallelic_eligibility=retain_biallelic_eligibility,
    )


def _finalize_external_allele_frequencies(
    matched_all, assignment, statistics, std_cols, eaf_col, resolved_policies, *,
    chromosome, emit, logger=None, ctx=None, rejects=None,
    validation_reference_col=None, validation_reference_file=None,
    validation_reference_column=None, colmap=None, declared_alt_frequency=True,
    retain_biallelic_eligibility=False,
):
    """Apply the same scientific finalization to fresh and approved cached joins."""
    stats = dict(statistics)
    input_rows = matched_all.height
    stats["input_rows_before_merge"] = input_rows
    min_match_fraction = stats["min_match_fraction"]
    load_path, source_column = stats["source_file"], stats["source_column"]
    row_col, biallelic_col = assignment["row"], assignment["biallelic"]
    matched_all, orientation_col, join_stats = finish_allele_assignments(
        matched_all, assignment["orientation"], assignment["value"], stats,
        output_column=eaf_col, orientations=stats["orientations"],
        frequency_columns=assignment["frequency_columns"],
        frequency_reason_column=assignment["reason"], frequency_biallelic_column=biallelic_col,
        site_columns=assignment["site_columns"], site_keys=[std_cols["chr"], std_cols["pos"]],
        warn=emit.warn, reference_label="external frequency file '%s'" % load_path,
    )
    stats = join_stats
    raw_frequency_screen = dict(
        join_stats.get("raw_frequency_screen") or {}
    )
    stats["raw_frequency_screen"] = raw_frequency_screen
    stats["maf_reference_decision"] = "not_required"
    stats["maf_reference_validation"] = {}
    if (
        not declared_alt_frequency
        and raw_frequency_screen.get("maf_like") is True
        and validation_reference_col is not None
    ):
        if validation_reference_col not in matched_all.columns:
            raise AlleleFrequencyError(
                "Internal PostGWAS execution-order error: external AF was "
                "MAF-like, but the aligned default reference column '%s' was "
                "not available for confirmation."
                % validation_reference_col
            )
        candidate_col = unused_column_name("__postgwas_external_raw_af", matched_all)
        expected_col = unused_column_name(
            "__postgwas_external_validation_eaf", matched_all,
            reserved=(candidate_col,),
        )
        swapped = pl.col(orientation_col) == 1
        # The existing EAF/MAF diagnostic complements the opposite allele.
        # This comparison is defined only for two-allele sites; at a known
        # multiallelic site 1-EAF includes all remaining alleles, not just OA.
        semantic_rows = matched_all.filter(pl.col(biallelic_col))
        stats["maf_reference_multiallelic_or_unresolved_rows_excluded"] = (
            matched_all.height - semantic_rows.height
        )
        semantic_frame = semantic_rows.select([
            *[pl.col(std_cols[name]) for name in ("chr", "pos", "ea", "oa")],
            pl.when(swapped)
            .then(1.0 - pl.col(eaf_col).cast(pl.Float64, strict=False))
            .otherwise(pl.col(eaf_col).cast(pl.Float64, strict=False))
            .alias(candidate_col),
            pl.when(swapped)
            .then(
                1.0
                - pl.col(validation_reference_col).cast(
                    pl.Float64, strict=False,
                )
            )
            .otherwise(
                pl.col(validation_reference_col).cast(
                    pl.Float64, strict=False,
                )
            )
            .alias(expected_col),
        ])
        decision, validation = confirm_eaf_with_reference(
            semantic_frame,
            candidate_col,
            validation_reference_file or "aligned default EAF reference",
            validation_reference_column or validation_reference_col,
            colmap,
            std_cols,
            policies=resolved_policies,
            logger=logger,
            ctx=ctx,
            aligned_reference_col=expected_col,
            study_columns_canonical=True,
        )
        stats["maf_reference_decision"] = decision
        stats["maf_reference_validation"] = validation

    # The study row is already reference-oriented. The external mapping declares
    # which columns are REF and ALT/effect allele, so the shared direct/swap join
    # has already put its frequency on the final study EA.
    # A palindrome is retained only when upstream strand_action proves that this
    # orientation was resolved; a reference-panel frequency must never be used
    # to infer the missing study strand.
    is_palindromic = palindromic_snp_expression(
        pl.col(std_cols["ea"]), pl.col(std_cols["oa"])
    )
    has_key_match = pl.col(orientation_col).is_not_null()
    matched_palindromic = is_palindromic & has_key_match
    orientation_resolved = (
        pl.col(STRAND_ACTION_COLUMN).is_in(list(RESOLVED_STRAND_ACTIONS))
        if STRAND_ACTION_COLUMN in matched_all.columns
        else pl.lit(False)
    ).fill_null(False)
    state_col = unused_column_name("palindromic_orientation_state", matched_all)
    palindromic_orientation = (
        matched_all.filter(matched_palindromic)
        .select(row_col, orientation_resolved.alias("__orientation_resolved"))
        .with_columns(
            pl.when(pl.col("__orientation_resolved"))
            .then(pl.lit("resolved_by_study_orientation"))
            .otherwise(pl.lit("orientation_unavailable"))
            .alias(state_col)
        )
        .drop("__orientation_resolved")
    )
    matched_all = matched_all.join(
        palindromic_orientation, on=row_col, how="left", coalesce=True
    ).with_columns(
        pl.when(
            matched_palindromic
            & (
                pl.col(state_col)
                != "resolved_by_study_orientation"
            ).fill_null(True)
        )
        .then(None)
        .otherwise(pl.col(eaf_col))
        .alias(eaf_col)
    )
    palindromic_state_counts = {
        str(row[state_col]): int(row["len"])
        for row in palindromic_orientation.group_by(
            state_col
        ).len().to_dicts()
    }
    stats["palindromic_evaluated_rows"] = palindromic_orientation.height
    stats["palindromic_orientation_resolved_rows"] = int(
        palindromic_state_counts.get("resolved_by_study_orientation", 0)
    )
    stats["palindromic_orientation_rejected_rows"] = (
        palindromic_orientation.height
        - stats["palindromic_orientation_resolved_rows"]
    )
    stats["palindromic_orientation_state_counts"] = palindromic_state_counts

    # An allele-key match is only evidence that the panel contains the variant;
    # it does not supply an EAF when the matched value is null, non-finite or
    # outside the probability interval. Count all value states in one Polars
    # aggregation so the configured coverage gate measures usable frequencies.
    reference_value = pl.col(eaf_col).cast(pl.Float64, strict=False)
    finite_value = reference_value.is_finite().fill_null(False)
    in_frequency_range = (
        (reference_value >= 0.0) & (reference_value <= 1.0)
    ).fill_null(False)
    usable_value = has_key_match & finite_value & in_frequency_range
    value_counts = matched_all.select([
        ((pl.col(orientation_col) == 0) & usable_value).sum().alias("direct_usable"),
        ((pl.col(orientation_col) == 1) & usable_value).sum().alias("swapped_usable"),
        usable_value.sum().alias("usable"),
        (has_key_match & reference_value.is_null()).sum().alias("matched_missing"),
        (
            has_key_match & reference_value.is_not_null() & ~finite_value
        ).sum().alias("matched_non_finite"),
        (
            has_key_match & finite_value & ~in_frequency_range
        ).sum().alias("matched_out_of_range"),
    ]).row(0, named=True)

    stats["direct_match_rows"] = join_stats["direct_key_matches"]
    stats["flip_match_rows"] = join_stats["swapped_key_matches"]
    stats["key_match_rows"] = stats["direct_match_rows"] + stats["flip_match_rows"]
    stats["key_match_fraction"] = (
        stats["key_match_rows"] / input_rows if input_rows else 0.0
    )
    stats["usable_direct_match_rows"] = int(value_counts["direct_usable"] or 0)
    stats["usable_flip_match_rows"] = int(value_counts["swapped_usable"] or 0)
    stats["usable_match_rows"] = int(value_counts["usable"] or 0)
    stats["usable_match_fraction"] = (
        stats["usable_match_rows"] / input_rows if input_rows else 0.0
    )
    stats["matched_missing_frequency_rows"] = int(
        value_counts["matched_missing"] or 0
    )
    stats["matched_non_finite_frequency_rows"] = int(
        value_counts["matched_non_finite"] or 0
    )
    stats["matched_out_of_range_frequency_rows"] = int(
        value_counts["matched_out_of_range"] or 0
    )
    stats["matched_unusable_frequency_rows"] = (
        stats["matched_missing_frequency_rows"]
        + stats["matched_non_finite_frequency_rows"]
        + stats["matched_out_of_range_frequency_rows"]
    )

    # ==================================================
    # STEP 3: UNMATCHED -> NULL EAF from the left join
    # ==================================================
    stats["unmatched_rows"] = join_stats["unmatched_rows"]
    stats["total_missing"] = (
        stats["unmatched_rows"] + stats["matched_missing_frequency_rows"]
    )
    stats["missing_fraction"] = (
        stats["total_missing"] / input_rows if input_rows else 0.0
    )
    stats["unusable_frequency_rows"] = input_rows - stats["usable_match_rows"]
    # Keep the established public name, but make its meaning scientifically
    # correct: it is now the usable EAF share, not merely the allele-key share.
    stats["match_fraction"] = stats["usable_match_fraction"]
    final_df = matched_all
    if stats["palindromic_evaluated_rows"]:
        emit.info(
            "External EAF for palindromic variants: %s evaluated after study "
            "orientation, %s retained, and %s rejected because orientation "
            "evidence was unavailable."
            % (
                _fmt_int(stats["palindromic_evaluated_rows"]),
                _fmt_int(stats["palindromic_orientation_resolved_rows"]),
                _fmt_int(stats["palindromic_orientation_rejected_rows"]),
            )
        )
    final_df, _ = reject_rows(
        final_df,
        matched_palindromic
        & (
            pl.col(state_col)
            != "resolved_by_study_orientation"
        ).fill_null(True),
        step_label=STEP_LABEL,
        reason="palindromic_orientation_unavailable",
        context=emit,
        collector=rejects,
        detail=(
            "external reference AF was available, but no resolved upstream "
            "strand_action proved the study palindrome's orientation"
        ),
        check_name="external-EAF palindromic orientation",
        description=(
            "A matched A/T or C/G external frequency is used only after automatic "
            "consensus/internal-study-EAF resolution or an explicit "
            "reference_aligned declaration has oriented the study row."
        ),
        warn_on_remove=True,
        warn_without_collector=True,
        record_empty=True,
    )

    final_df = final_df.drop([
        orientation_col,
        row_col,
        state_col,
        *([biallelic_col] if biallelic_col is not None and not retain_biallelic_eligibility else []),
    ])
    if retain_biallelic_eligibility:
        stats["biallelic_eligibility_column"] = biallelic_col

    emit.info("Direct allele-key matches  : %s" % _fmt_int(stats["direct_match_rows"]))
    emit.info("Swapped allele-key matches : %s" % _fmt_int(stats["flip_match_rows"]))
    emit.info("No allele-key match        : %s" % _fmt_int(stats["unmatched_rows"]))
    emit.info(
        "External frequency file %s, source column %r, was assigned to study "
        "column %r: %s of %s variants have usable frequencies (%s)."
        % (
            load_path,
            source_column,
            eaf_col,
            _fmt_int(stats["usable_match_rows"]),
            _fmt_int(input_rows),
            _fmt_pct(stats["usable_match_fraction"]),
        )
    )
    if stats["matched_unusable_frequency_rows"]:
        emit.warn(
            "%s allele-key matches did not supply a usable EAF: %s missing, %s "
            "non-finite and %s outside 0 to 1."
            % (
                _fmt_int(stats["matched_unusable_frequency_rows"]),
                _fmt_int(stats["matched_missing_frequency_rows"]),
                _fmt_int(stats["matched_non_finite_frequency_rows"]),
                _fmt_int(stats["matched_out_of_range_frequency_rows"]),
            )
        )

    if input_rows and stats["match_fraction"] < min_match_fraction:
        raise AlleleFrequencyError(
            "The external frequency file '%s' matched chromosome, position and alleles for %s of %s "
            "variants on chromosome %s (%s), but supplied a usable EAF for only %s variants (%s). "
            "A usable EAF is finite and between 0 and 1. This is below the %s required by "
            "eaf.external_min_match_fraction. Among allele-key matches, %s frequencies were missing, "
            "%s were non-finite and %s were outside 0 to 1; %s variants had no allele-key match. "
            "Check the reference panel's genome build, chromosome and allele columns, frequency "
            "column, and completeness, or explicitly lower the threshold for an intentionally "
            "limited panel."
            % (
                load_path,
                _fmt_int(stats["key_match_rows"]),
                _fmt_int(input_rows),
                chromosome,
                _fmt_pct(stats["key_match_fraction"]),
                _fmt_int(stats["usable_match_rows"]),
                _fmt_pct(stats["usable_match_fraction"]),
                _fmt_pct(min_match_fraction),
                _fmt_int(stats["matched_missing_frequency_rows"]),
                _fmt_int(stats["matched_non_finite_frequency_rows"]),
                _fmt_int(stats["matched_out_of_range_frequency_rows"]),
                _fmt_int(stats["unmatched_rows"]),
            )
        )
    return final_df, eaf_col, stats


# -------------------------------------------------------
# 10. Helper: Write structured statistics into qc_info
# -------------------------------------------------------
def _add_prefixed_statistics(
    qc_info: Dict[str, Any],
    prefix: str,
    stats: Dict[str, Any]
) -> None:
    for key, value in stats.items():
        if isinstance(value, dict):
            for sub_key, sub_val in value.items():
                qc_info["%s_%s_%s" % (prefix, key, sub_key)] = sub_val
        else:
            qc_info["%s_%s" % (prefix, key)] = value


def _out_of_bound_frame(frame: pl.DataFrame, colname: str) -> pl.DataFrame:
    """Rows whose frequency is non-finite or outside 0..1 before correction."""
    if colname not in frame.columns:
        return frame.head(0)
    # Cast defensively: a frequency column read with schema inference disabled
    # arrives as Utf8, and comparing that against a float raises ComputeError.
    # Every sibling module casts with strict=False first; match them.
    value = pl.col(colname).cast(pl.Float64, strict=False)
    return frame.filter(
        (
            value.is_not_null()
            & (~value.is_finite() | (value < 0.0) | (value > 1.0))
        ).fill_null(False)
    )


def _resolve_reference_aligned_duplicates(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    output_layout: Dict[str, str],
    output_delimiter: str,
    policies,
    logger=None,
    ctx=None,
    rejects=None,
) -> Tuple[pl.DataFrame, Dict[str, Any]]:
    """Resolve rows that become identical only after reference orientation.

    Strand orientation has already made the study effect allele reference ALT
    and the other allele reference REF, and the final effect-frequency source
    has been aligned to ALT. Running the shared duplicate resolver at this
    boundary catches opposite-strand representations without applying a
    complement-based key to unnormalised indels. With the canonical key,
    distinct alternate alleles at one coordinate remain distinct because both
    REF and ALT participate in the key.
    """
    retained, report, statistics = resolve_duplicate_variants(
        df,
        sample_column_dict,
        policies=policies,
        logger=logger,
        rejects=rejects,
        step_label=STRAND_STEP_LABEL,
        reference_aligned=True,
    )
    report_path = configured_output_path(
        sample_column_dict["output_folder"],
        output_layout["post_orientation_duplicates"],
        error_type=AlleleFrequencyError,
        dataset_id=sample_column_dict["gwas_outputname"],
        chromosome=chromosome,
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report.write_csv(report_path, separator=output_delimiter)
    statistics["report"] = str(report_path)
    if ctx is not None:
        ctx.extra["post_orientation_duplicates"] = dict(statistics)
        if report.height:
            ctx.info(
                "Reference-aligned duplicate evidence for %s row%s was written "
                "to %s."
                % (
                    "{:,}".format(report.height),
                    "" if report.height == 1 else "s",
                    report_path,
                )
            )

    unsafe_groups = (
        int(statistics["conflicting_groups"])
        + int(statistics["swapped_orientation_groups"])
    )
    if (
        unsafe_groups
        and str(policies.get("duplicates.conflicting_action"))
        == "fail_dataset"
    ):
        raise AlleleFrequencyError(
            "Chromosome %s contains %s unsafe reference-aligned duplicate "
            "group%s. Their rows were recorded in %s and the rejected-variant "
            "output; processing stopped because duplicates.conflicting_action "
            "is 'fail_dataset'."
            % (
                chromosome,
                "{:,}".format(unsafe_groups),
                "" if unsafe_groups == 1 else "s",
                report_path,
            )
        )
    return retained, statistics


# -------------------------------------------------------
# 11. MAIN FUNCTION
# -------------------------------------------------------
def harmonise_allele_frequency(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    output_layout: Dict[str, str],
    output_delimiter: str,
    eaffile: Optional[str] = None,
    external_eaf_colmap: Optional[Dict] = None,
    default_eaf_file: Optional[str] = None,
    default_eaf_column: Optional[str] = None,
    default_eaf_colmap: Optional[Dict] = None,
    policies=None,
    logger=None,
    rejects=None,
    step_context=None,
    step_number: int = 4,
    step_total: int = 16,
    study_decision: Optional[Any] = None,
) -> Tuple[pl.DataFrame, Dict, dict]:
    """Give every variant an effect allele frequency.

    Parameters added by the plan v3 migration, all optional:

    policies       a ``Policies`` object; ``None`` means the registry defaults.
    logger         a ``PipelineLogger``.  When given, this function opens its own
                   ``step()`` and nothing is written to stdout.
    rejects        a ``RejectCollector``.  Removed variants are written to the
                   per-chromosome rejected-variants file with their reason.
    step_context   an already-open ``StepContext``, when the caller prefers to
                   own the step.
    study_decision the dataset-level decision from step 07, e.g.
                   ``{"eaf_is_maf": False}``. The decision is made once from
                   the whole study and then applied consistently to each chromosome.
    default_eaf_*  the raw chromosome REF/ALT/population-AF reference used for
                   strand orientation and independent MAF/EAF confirmation.

    A selected replacement is a whole-column source decision, supported by the
    dataset's MAF evidence. Its declared ALT frequencies govern numerical QC;
    missing/invalid original MAF values remain unchanged as provenance and do
    not independently reject a row with a valid replacement EAF. An invalid
    individual input cell never activates replacement by itself.
    """
    gwas_outputname = sample_column_dict["gwas_outputname"]
    output_dir = sample_column_dict["output_folder"]
    merge_stats = {}

    pol = resolve_policies(policies)
    degenerate_action = str(pol.get("eaf.degenerate"))

    # Create subdirectory for EAF QC outputs
    eaf_qc_dir = configured_output_path(
        output_dir, output_layout["frequency_qc_directory"],
    )
    eaf_qc_dir.mkdir(parents=True, exist_ok=True)

    qc_info: Dict[str, Any] = {
        "chromosome": chromosome,
        "gwas_outputname": gwas_outputname,
        "Total_number_of_variants": df.height,
        "final_status": "started",
        "decision_source": None,
        "eaf_provenance": None,
        "final_eaf_col": None,
    }

    final_df = None
    with shared_step_context(
        logger,
        step_context,
        number=step_number,
        total=step_total,
        title=STEP_TITLE,
        operation=FUNC_NAME,
        rows_in=df.height,
        policy_keys=list(dict.fromkeys(STRAND_POLICY_KEYS + POLICY_KEYS)),
    ) as ctx:
        emit = ctx
        emit.info(
            "Giving every variant on chromosome %s an effect allele frequency. Starting with %s variants."
            % (chromosome, _fmt_int(df.height))
        )

        try:
            validate_scientific_column_mappings(
                sample_column_dict, stage=STEP_LABEL, error_type=AlleleFrequencyError,
            )
            std_cols = {
                "chr": sample_column_dict["chr_col"],
                "pos": sample_column_dict["pos_col"],
                "ea": sample_column_dict["ea_col"],
                "oa": sample_column_dict["oa_col"],
            }

            replacement_cache = replacement_cache_descriptor = None
            external_approved = (isinstance(study_decision, dict)
                and study_decision.get("eaf_selected_source") != "replacement_reference" and (
                study_decision.get("eaf_external_validation_required") is True
                or "eaf_external_validation" in study_decision))
            if (isinstance(study_decision, dict)
                    and (study_decision.get("eaf_selected_source") == "replacement_reference" or external_approved)):
                eaffile = _require_replacement_validation(chromosome, df, sample_column_dict,
                    study_decision, eaffile, external_eaf_colmap, pol)
                binding = study_decision["eaf_external_validation" if external_approved else "eaf_replacement_validation"]["contract"]["chromosomes"][str(chromosome)]
                replacement_cache_descriptor = binding.get("assignment_cache")
                if replacement_cache_descriptor is not None:
                    cache_column = replacement_cache_descriptor["column"]
                    replacement_cache = df.select(SOURCE_INPUT_ROW_COLUMN, cache_column)
                    df = df.drop(cache_column)
                initial_column = (study_decision.get("eaf_initial_frequency") or {}).get("column")
                if initial_column is not None and initial_column in df.columns:
                    df = df.drop(initial_column)

            df, strand_qc, sample_column_dict = harmonise_strand_orientation(
                chromosome=chromosome,
                df=df,
                sample_column_dict=sample_column_dict,
                reference_file=default_eaf_file,
                population_col=default_eaf_column,
                reference_colmap=default_eaf_colmap,
                study_decision=study_decision,
                policies=pol,
                logger=logger,
                ctx=ctx,
                rejects=rejects,
            )
            qc_info["strand_orientation"] = strand_qc

            internal_eaf = optional_text(sample_column_dict.get("eaf_col"))
            replacement_selected = (
                isinstance(study_decision, dict)
                and study_decision.get("eaf_selected_source") == "replacement_reference"
            )
            replacement = study_decision.get("eaf_replacement", {}) if replacement_selected else {}
            if replacement_selected and (
                study_decision.get("eaf_input_interpretation") != "maf"
                or study_decision.get("eaf_reference_validation", {}).get("status") != "replacement_selected"
                or study_decision.get("eaf_is_maf") is not True
                or optional_text(replacement.get("column")) is None
            ):
                raise AlleleFrequencyError("Selected reference EAF replacement lacks a validated whole-study MAF decision.")
            qc_info.update(
                eaf_selected_source="replacement_reference" if replacement_selected else
                "study_supplied" if internal_eaf else "external_eaf",
                eaf_input_interpretation=(study_decision.get("eaf_input_interpretation")
                    if isinstance(study_decision, dict) else None),
            )
            if external_approved:
                qc_info["eaf_external_validation"] = deepcopy(study_decision["eaf_external_validation"])
            if replacement_selected:
                qc_info["eaf_replacement"] = dict(replacement)
                qc_info["eaf_replacement_validation"] = deepcopy(study_decision["eaf_replacement_validation"])

            if external_eaf_colmap is None:
                raise AlleleFrequencyError(
                    "External EAF column mapping is missing from the resolved YAML configuration."
                )

            final_eaf_col = None
            validated_values = None
            validated_statistics = None
            out_of_bound_df = df.head(0)
            study_says_maf = _study_decision_flag(study_decision)
            qc_info["study_decision_eaf_is_maf"] = study_says_maf

            # ------------------------------
            # STEP 1: Internal EAF
            # ------------------------------
            if internal_eaf and internal_eaf in df.columns and not replacement_selected:
                emit.info("The study supplies its own frequency column '%s'." % internal_eaf)
                qc_info["eaf_provenance"] = "study_supplied"

                # Capture the original invalid values before the configured
                # reject, null, clip or fail action changes the working frame.
                out_of_bound_df = _out_of_bound_frame(df, internal_eaf)

                valid_range, is_suspicious, out_of_range, n_missing, cleaned_df, low_freq_pct, stats = \
                    validate_allele_frequency(
                        df,
                        internal_eaf,
                        policies=pol,
                        logger=logger,
                        ctx=ctx,
                        rejects=rejects,
                        study_eaf_accepted=(
                            isinstance(study_decision, dict)
                            and study_decision.get("eaf_reference_validation", {}).get("status")
                            in ("passed", "not_required")
                        ),
                    )

                validated_statistics = stats.get("final")
                validated_values = (
                    cleaned_df.get_column(internal_eaf).clone()
                    if validated_statistics is not None else None
                )
                qc_info["Internal_AF"] = "present"
                _add_prefixed_statistics(qc_info, "internal_af", stats)
                qc_info["internal_af_valid_range"] = valid_range

                # NOTE: `valid_range` is reported, not used as the gate for
                # adopting the column.  eaf.out_of_range has already decided
                # what happens to out-of-range values; diverting an otherwise
                # usable column to the reference panel because a few values were
                # bad would change behaviour silently.
                if cleaned_df.height == 0:
                    raise AlleleFrequencyError(
                        "No variants remain after validating effect-allele-frequency column "
                        "'%s'." % internal_eaf
                    )
                if study_says_maf is not None:
                    emit.decide(
                        "Frequency column '%s'" % internal_eaf,
                        "study-level decision, made once on the whole file before the split",
                        "a minor allele frequency" if study_says_maf else "an effect allele frequency",
                    )
                    if study_says_maf:
                        reference_file = optional_text(default_eaf_file)
                        reference_column = optional_text(default_eaf_column)
                        if reference_file is None or reference_column is None or default_eaf_colmap is None:
                            raise AlleleFrequencyError(
                                "Column '%s' is MAF-like, but the resolved default EAF file, population "
                                "column, or column mapping is missing." % internal_eaf
                            )
                        reference_decision, reference_stats = confirm_eaf_with_reference(
                            cleaned_df,
                            internal_eaf,
                            reference_file,
                            reference_column,
                            default_eaf_colmap,
                            std_cols,
                            policies=pol,
                            logger=logger,
                            ctx=ctx,
                            aligned_reference_col=REFERENCE_AF_COLUMN,
                            study_columns_canonical=True,
                        )
                        _add_prefixed_statistics(qc_info, "maf_reference", reference_stats)
                        qc_info["maf_reference_decision"] = reference_decision
                        if reference_decision == "maf":
                            raise AlleleFrequencyError(
                                "Column '%s' was confirmed as minor allele frequency rather than frequency "
                                "aligned to the listed effect allele. Provide an effect-allele-frequency "
                                "column or a dataset-supplied external EAF file." % internal_eaf
                            )
                        if reference_decision == "inconclusive":
                            minor_fraction = reference_stats.get(
                                "reference_effect_allele_minor_fraction"
                            )
                            minor_text = (
                                "not available"
                                if minor_fraction is None
                                else _fmt_pct(minor_fraction)
                            )
                            raise AlleleFrequencyError(
                                "Column '%s' is MAF-like, but the reference comparison was "
                                "inconclusive (%s). Comparable non-palindromic variants: %s "
                                "(minimum %s); %s; mean error as EAF: %s; mean error as MAF: %s; "
                                "required error separation: %s; reference-aligned EAF values at "
                                "or below 0.5: %s (uninformative above %s). PostGWAS will not "
                                "guess or apply the deferred allele-swap frequency flip. Provide "
                                "a frequency aligned to the listed effect allele, or configure a "
                                "dataset-supplied external EAF file."
                                % (
                                    internal_eaf,
                                    reference_stats.get("decision_reason", "unknown reason"),
                                    _fmt_int(reference_stats.get("comparable_variants", 0)),
                                    _fmt_int(reference_stats.get("minimum_overlap", 0)),
                                    _maf_reference_suitability_text(
                                        reference_stats
                                    ),
                                    reference_stats.get("mean_eaf_error"),
                                    reference_stats.get("mean_maf_error"),
                                    reference_stats.get("minimum_error_margin"),
                                    minor_text,
                                    _fmt_pct(reference_stats.get(
                                        "reference_minor_fraction_cutoff", 0.0
                                    )),
                                )
                            )
                        else:
                            if STRAND_ACTION_COLUMN in cleaned_df.columns:
                                swapped = pl.col(STRAND_ACTION_COLUMN).is_in([
                                    "forward_swapped", "reverse_complement_swapped"
                                ])
                                cleaned_df = cleaned_df.with_columns(
                                    pl.when(swapped)
                                    .then(1.0 - pl.col(internal_eaf))
                                    .otherwise(pl.col(internal_eaf))
                                    .alias(internal_eaf)
                                )
                            qc_info["decision_source"] = "internal_eaf_reference_validation"
                    else:
                        qc_info["decision_source"] = "internal_eaf_study_decision"
                    final_df = cleaned_df
                    final_eaf_col = internal_eaf
                elif not is_suspicious:
                    final_df = cleaned_df
                    final_eaf_col = internal_eaf
                    qc_info["decision_source"] = "internal_eaf_direct"
                else:
                    qc_info["final_status"] = "failed"
                    raise AlleleFrequencyError(
                        "Column '%s' has a frequency distribution consistent with a minor "
                        "allele frequency, but no study-wide MAF/EAF decision is available. "
                        "PostGWAS will not use an unoriented frequency column. Run the normal "
                        "dataset-level study-property step, or provide a dataset-supplied "
                        "external EAF file."
                        % internal_eaf
                    )
            else:
                qc_info["Internal_AF"] = "MAF retained as input provenance" if replacement_selected else "Absent"
                emit.info("Original input MAF retained; assigning the selected reference EAF after independent orientation."
                          if replacement_selected else "The study has no frequency column of its own.")

            # ------------------------------
            # STEP 2: Dataset-supplied external source
            # ------------------------------
            if final_df is None:
                emit.info(
                    "Taking the effect allele frequency from the dataset's configured "
                    "external frequency file."
                )

                original_replacement_file = optional_text(
                    replacement.get("resolved_files", {}).get(str(chromosome))
                ) if replacement_selected else None
                external_file = optional_text(eaffile) or original_replacement_file
                if external_file is None or not os.path.exists(external_file):
                    raise AlleleFrequencyError(
                        "No usable effect-allele-frequency source remains. Provide either "
                        "effect_allele_frequency_column or external_eaf_file plus "
                        "external_eaf_column in the sample sheet."
                    )
                target_file = external_file
                target_col = optional_text(replacement.get("column")) if replacement_selected else optional_text(sample_column_dict.get("eafcolumn"))
                if target_col is None:
                    raise AlleleFrequencyError("external_eaf_column is required with external_eaf_file.")
                target_map = external_eaf_colmap
                qc_info["external_eaf_file"] = os.path.basename(target_file)
                qc_info["external_eaf_column"] = target_col

                prepared_reference = (
                    study_decision.get("eaf_prepared_references", {}).get(str(chromosome), {})
                    if isinstance(study_decision, dict) else {}
                )
                default_original = prepared_reference.get("contract", {}).get("reference", {}).get("path") or default_eaf_file
                source_identity_path = binding["source"]["path"] if replacement_selected or external_approved else target_file
                same_reference_file = same_frequency_reference_source(
                    source_identity_path, target_col, target_map,
                    default_original, target_col, target_map,
                )
                same_reference_source = same_frequency_reference_source(
                    source_identity_path, target_col, target_map,
                    default_original, default_eaf_column, default_eaf_colmap,
                )
                if same_reference_source:
                    if REFERENCE_AF_COLUMN not in df.columns:
                        raise AlleleFrequencyError(
                            "Internal PostGWAS execution-order error: the explicitly selected "
                            "external EAF source is identical to the strand reference, but its "
                            "validated aligned frequency column is unavailable."
                        )
                    raw_frequency_screen = {} if replacement_selected or external_approved else strand_qc.get("raw_frequency_screen")
                    if not isinstance(raw_frequency_screen, dict):
                        raise AlleleFrequencyError(
                            "Internal PostGWAS execution-order error: the raw "
                            "external-frequency screen was not retained from "
                            "the identical strand reference."
                        )
                    reference_value = pl.col(REFERENCE_AF_COLUMN).cast(
                        pl.Float64, strict=False
                    )
                    finite_value = reference_value.is_finite().fill_null(False)
                    in_range = (
                        (reference_value >= 0.0) & (reference_value <= 1.0)
                    ).fill_null(False)
                    usable_value = finite_value & in_range
                    value_counts = df.select([
                        usable_value.sum().alias("usable"),
                        reference_value.is_null().sum().alias("missing"),
                        (
                            reference_value.is_not_null() & ~finite_value
                        ).sum().alias("non_finite"),
                        (
                            reference_value.is_not_null()
                            & finite_value
                            & ~in_range
                        ).sum().alias("out_of_range"),
                        palindromic_snp_expression(
                            pl.col(std_cols["ea"]), pl.col(std_cols["oa"])
                        ).sum().alias("palindromic"),
                    ]).to_dicts()[0]
                    usable_rows = int(value_counts["usable"] or 0)
                    usable_fraction = usable_rows / df.height if df.height else 0.0
                    minimum_fraction = float(
                        pol.get("eaf.external_min_match_fraction")
                    )
                    if df.height and usable_fraction < minimum_fraction:
                        raise AlleleFrequencyError(
                            "The external frequency file '%s' is identical to the validated "
                            "strand reference, but its aligned column '%s' supplied a usable EAF "
                            "for only %s of %s retained variants on chromosome %s (%s). This is "
                            "below the %s required by eaf.external_min_match_fraction."
                            % (
                                target_file,
                                target_col,
                                _fmt_int(usable_rows),
                                _fmt_int(df.height),
                                chromosome,
                                _fmt_pct(usable_fraction),
                                _fmt_pct(minimum_fraction),
                            )
                        )
                    actions = strand_qc.get("actions", {})
                    direct_rows = int(actions.get("forward", 0) or 0) + int(
                        actions.get("reverse_complement", 0) or 0
                    )
                    swapped_rows = int(
                        actions.get("forward_swapped", 0) or 0
                    ) + int(actions.get("reverse_complement_swapped", 0) or 0)
                    palindromic_rows = int(value_counts["palindromic"] or 0)
                    missing_rows = int(value_counts["missing"] or 0)
                    non_finite_rows = int(value_counts["non_finite"] or 0)
                    out_of_range_rows = int(value_counts["out_of_range"] or 0)
                    merge_stats = {
                        "external_rows_loaded": 0,
                        "input_rows_before_merge": df.height,
                        "direct_match_rows": direct_rows,
                        "flip_match_rows": swapped_rows,
                        "key_match_rows": df.height,
                        "key_match_fraction": 1.0 if df.height else 0.0,
                        "usable_match_rows": usable_rows,
                        "usable_match_fraction": usable_fraction,
                        "match_fraction": usable_fraction,
                        "matched_missing_frequency_rows": missing_rows,
                        "matched_non_finite_frequency_rows": non_finite_rows,
                        "matched_out_of_range_frequency_rows": out_of_range_rows,
                        "matched_unusable_frequency_rows": (
                            missing_rows + non_finite_rows + out_of_range_rows
                        ),
                        "unusable_frequency_rows": df.height - usable_rows,
                        "unmatched_rows": 0,
                        "total_missing": missing_rows,
                        "missing_fraction": (
                            missing_rows / df.height if df.height else 0.0
                        ),
                        "min_match_fraction": minimum_fraction,
                        "palindromic_evaluated_rows": palindromic_rows,
                        "palindromic_orientation_resolved_rows": palindromic_rows,
                        "palindromic_orientation_rejected_rows": 0,
                        "palindromic_orientation_state_counts": (
                            {"resolved_by_study_orientation": palindromic_rows}
                            if palindromic_rows
                            else {}
                        ),
                        "reused_aligned_strand_reference": True,
                        "raw_frequency_screen": raw_frequency_screen,
                        "maf_reference_decision": "not_required",
                        "maf_reference_validation": {},
                    }
                    merged_df = df
                    new_col_name = REFERENCE_AF_COLUMN
                    if external_approved and not same_frequency_reference_source(
                        target_file, target_col, target_map, default_original, default_eaf_column, default_eaf_colmap,
                    ):
                        # Shared staging used to attach this exact source via a
                        # separate file. Preserve that public destination schema
                        # while reusing the already aligned value, without I/O.
                        new_col_name = unused_column_name(target_col, merged_df)
                        merged_df = merged_df.with_columns(pl.col(REFERENCE_AF_COLUMN).alias(new_col_name))
                    emit.info(
                        "Reused the allele-aligned population AF already retained from "
                        "strand orientation: external EAF file %s, source column %r, "
                        "study column %r, %s usable frequencies among %s variants. "
                        "The explicitly selected file, column and allele mapping are "
                        "identical; the reference table was not read a second time."
                        % (
                            target_file, target_col, new_col_name,
                            _fmt_int(usable_rows), _fmt_int(merged_df.height),
                        )
                    )
                elif replacement_selected or external_approved:
                    if replacement_cache_descriptor is None:
                        raise AlleleFrequencyError(_frequency_validation_label(external_approved) + " lacks its approved assignment cache.")
                    merged_df, new_col_name, merge_stats = _consume_replacement_assignments(
                        df, replacement_cache, replacement_cache_descriptor, std_cols, pol,
                        chromosome=chromosome, emit=emit, logger=logger, ctx=ctx, rejects=rejects)
                else:
                    merged_df, new_col_name, merge_stats = merge_external_allele_frequencies(
                        df,
                        target_file,
                        target_col,
                        target_map,
                        std_cols,
                        chromosome,
                        policies=pol,
                        logger=logger,
                        ctx=ctx,
                        rejects=rejects,
                        study_columns_canonical=True,
                        validation_reference_col=(
                            None
                            if same_reference_file or replacement_selected
                            else REFERENCE_AF_COLUMN
                        ),
                        validation_reference_file=default_eaf_file,
                        validation_reference_column=default_eaf_column,
                        declared_alt_frequency=replacement_selected,
                    )

                merge_stats.update(
                    source_file=target_file,
                    source_column=target_col,
                    destination_column=new_col_name,
                )
                raw_frequency_screen = dict(
                    merge_stats.get("raw_frequency_screen") or {}
                )
                _add_prefixed_statistics(
                    qc_info,
                    "external_raw_af",
                    raw_frequency_screen,
                )
                if raw_frequency_screen and not external_approved:
                    emit.decide(
                        "Raw external frequency column '%s'" % target_col,
                        {
                            "usable frequencies": _fmt_int(
                                raw_frequency_screen.get("usable", 0)
                            ),
                            "at or below 0.5": "%s (%s)"
                            % (
                                _fmt_int(raw_frequency_screen.get(
                                    "at_or_below_0_5", 0
                                )),
                                _fmt_pct(raw_frequency_screen.get(
                                    "low_fraction_of_usable", 0.0
                                )),
                            ),
                            "decision cutoff": raw_frequency_screen.get(
                                "decision_cutoff"
                            ),
                        },
                        (
                            "MAF-like and requiring independent reference confirmation"
                            if raw_frequency_screen.get("maf_like") is True
                            else (
                                "not MAF-like; retaining the declared external "
                                "ALT/EAF contract"
                            )
                        ),
                    )
                if raw_frequency_screen.get("maf_like") is True and not external_approved:
                    if same_reference_file:
                        raise AlleleFrequencyError(
                            "The external AF column '%s' is MAF-like: %s of %s "
                            "usable values are at or below 0.5 (%s), exceeding "
                            "eaf.maf_decision_cutoff=%s. external_eaf_file and "
                            "the resource selected by "
                            "modules.harmonisation.default_eaf.source resolve "
                            "to the same file '%s', so PostGWAS cannot independently "
                            "determine whether the column is ALT/EAF or MAF. "
                            "Configure a different modules.harmonisation.default_eaf.source "
                            "and rerun."
                            % (
                                target_col,
                                _fmt_int(raw_frequency_screen.get(
                                    "at_or_below_0_5", 0
                                )),
                                _fmt_int(raw_frequency_screen.get("usable", 0)),
                                _fmt_pct(raw_frequency_screen.get(
                                    "low_fraction_of_usable", 0.0
                                )),
                                raw_frequency_screen.get("decision_cutoff"),
                                target_file,
                            )
                        )
                    external_decision = merge_stats.get(
                        "maf_reference_decision"
                    )
                    external_validation = dict(
                        merge_stats.get("maf_reference_validation") or {}
                    )
                    _add_prefixed_statistics(
                        qc_info,
                        "external_maf_reference",
                        external_validation,
                    )
                    qc_info["external_maf_reference_decision"] = (
                        external_decision
                    )
                    if external_decision == "maf":
                        raise AlleleFrequencyError(
                            "The external AF column '%s' was confirmed as minor "
                            "allele frequency rather than ALT/effect-allele "
                            "frequency. Configure an external_eaf_column containing "
                            "the frequency of the file's declared ALT allele."
                            % target_col
                        )
                    if external_decision != "eaf":
                        raise AlleleFrequencyError(
                            "The external AF column '%s' is MAF-like, but comparison "
                            "with the independent default EAF reference was "
                            "inconclusive (%s). Comparable non-palindromic variants: "
                            "%s (minimum %s); %s; mean error as EAF: %s; mean error as "
                            "MAF: %s; required error separation: %s. PostGWAS will "
                            "not guess. Configure a different external ALT/EAF "
                            "column or a more suitable independent "
                            "modules.harmonisation.default_eaf reference."
                            % (
                                target_col,
                                external_validation.get(
                                    "decision_reason", "unknown reason"
                                ),
                                _fmt_int(external_validation.get(
                                    "comparable_variants", 0
                                )),
                                _fmt_int(external_validation.get(
                                    "minimum_overlap", 0
                                )),
                                _maf_reference_suitability_text(
                                    external_validation
                                ),
                                external_validation.get("mean_eaf_error"),
                                external_validation.get("mean_maf_error"),
                                external_validation.get("minimum_error_margin"),
                            )
                        )

                _add_prefixed_statistics(qc_info, "external_merge", merge_stats)

                # Capture invalid reference values before the configured
                # frequency-range action changes the working frame.
                out_of_bound_df = _out_of_bound_frame(merged_df, new_col_name)

                _, _, _, _, cleaned_df, _, stats = validate_allele_frequency(
                    merged_df,
                    new_col_name,
                    policies=pol,
                    logger=logger,
                    ctx=ctx,
                    rejects=rejects,
                    study_eaf_accepted=replacement_selected or external_approved,
                )
                if replacement_selected:
                    # Selected output frequency has a separate column; original
                    # MAF remains unchanged for provenance and checkpoint binding.
                    destination = unused_column_name(target_col, cleaned_df)
                    cleaned_df = cleaned_df.with_columns(pl.col(new_col_name).alias(destination))
                    new_col_name = destination
                    qc_info["eaf_replacement_source_column"] = target_col
                    qc_info["eaf_replacement_source_file"] = original_replacement_file or target_file

                validated_statistics = stats.get("final")
                validated_values = (
                    cleaned_df.get_column(new_col_name).clone()
                    if validated_statistics is not None else None
                )
                _add_prefixed_statistics(qc_info, "external_af", stats)

                final_df = cleaned_df
                final_eaf_col = new_col_name
                if qc_info["decision_source"] != "internal_eaf_rejected_inconclusive":
                    qc_info["decision_source"] = "replacement_reference" if replacement_selected else "external_eaf"
                qc_info["eaf_provenance"] = "reference_imputed"

            # ------------------------------
            # FINALIZE
            # ------------------------------
            if final_df is None or final_eaf_col is None:
                raise AlleleFrequencyError("Failed to determine which column holds the effect allele frequency.")

            # The reference has now fixed ALT/REF orientation, and any deferred
            # MAF-like-column proof has also applied the required 1-EAF change
            # to swapped rows. Comparing duplicates before this boundary could
            # falsely call an EAF disagreement that step 04 was still required
            # to reconcile.
            sample_column_dict["eaf_col"] = final_eaf_col
            validate_scientific_column_mappings(
                sample_column_dict, stage=STEP_LABEL, error_type=AlleleFrequencyError,
            )
            final_df, post_orientation_duplicate_qc = (
                _resolve_reference_aligned_duplicates(
                    chromosome=chromosome,
                    df=final_df,
                    sample_column_dict=sample_column_dict,
                    output_layout=output_layout,
                    output_delimiter=output_delimiter,
                    policies=pol,
                    logger=logger,
                    ctx=ctx,
                    rejects=rejects,
                )
            )
            qc_info["post_orientation_duplicates"] = (
                post_orientation_duplicate_qc
            )

            missing_eaf_df = final_df.filter(pl.col(final_eaf_col).is_null())
            missing_eaf_count = missing_eaf_df.height
            out_of_bound_count = out_of_bound_df.height

            # Save QC files inside subdirectory
            missing_file = configured_output_path(
                output_dir,
                output_layout["missing_eaf"],
                dataset_id=gwas_outputname,
                chromosome=chromosome,
            )
            outofbound_file = configured_output_path(
                output_dir,
                output_layout["out_of_range_eaf"],
                dataset_id=gwas_outputname,
                chromosome=chromosome,
            )

            if missing_eaf_count > 0:
                missing_eaf_df.write_csv(missing_file, separator=output_delimiter)
                emit.info("Variants with no frequency were written to %s" % missing_file)

            if out_of_bound_count > 0:
                out_of_bound_df.write_csv(outofbound_file, separator=output_delimiter)
                emit.info("Variants whose frequency was outside 0 to 1 were written to %s" % outofbound_file)

            qc_info["missing_eaf_count"] = missing_eaf_count
            qc_info["out_of_bound_eaf_count"] = out_of_bound_count

            # Reference AF is supporting concordance evidence after alleles and
            # the final EAF source are aligned. It must not arbitrate a column
            # that remains only provisionally accepted as possible MAF.
            frequency_is_aligned = (
                replacement_selected or not study_says_maf
                or qc_info.get("maf_reference_decision") == "eaf"
            )
            if REFERENCE_AF_COLUMN in final_df.columns and frequency_is_aligned:
                af_difference_col = "strand_af_difference"
                tolerance = float(pol.get("strand.af_tolerance"))
                non_palindromic_action = str(
                    pol.get("strand.af_discordance_action")
                )
                palindromic_action = str(
                    pol.get("strand.palindromic_af_discordance_action")
                )
                final_df = final_df.with_columns(
                    (
                        pl.col(final_eaf_col).cast(pl.Float64, strict=False)
                        - pl.col(REFERENCE_AF_COLUMN).cast(pl.Float64, strict=False)
                    ).abs().alias(af_difference_col)
                )
                comparable = (
                    pl.col(final_eaf_col).is_not_null()
                    & pl.col(REFERENCE_AF_COLUMN).is_not_null()
                )
                palindromic = palindromic_snp_expression(
                    pl.col(std_cols["ea"]), pl.col(std_cols["oa"])
                )
                discordant = (
                    comparable & (pl.col(af_difference_col) > tolerance)
                ).fill_null(False)
                palindromic_discordant = discordant & palindromic
                non_palindromic_discordant = discordant & ~palindromic
                counts = final_df.select([
                    comparable.sum().alias("comparable"),
                    discordant.sum().alias("discordant"),
                    (comparable & palindromic).sum().alias("pal_comparable"),
                    palindromic_discordant.sum().alias("pal_discordant"),
                    (comparable & ~palindromic).sum().alias("non_pal_comparable"),
                    non_palindromic_discordant.sum().alias("non_pal_discordant"),
                ]).to_dicts()[0]
                n_comparable = int(counts["comparable"] or 0)
                n_discordant = int(counts["discordant"] or 0)
                n_pal_comparable = int(counts["pal_comparable"] or 0)
                n_pal_discordant = int(counts["pal_discordant"] or 0)
                n_non_pal_comparable = int(counts["non_pal_comparable"] or 0)
                n_non_pal_discordant = int(
                    counts["non_pal_discordant"] or 0
                )
                qc_info["strand_af_comparable"] = n_comparable
                qc_info["strand_af_discordant"] = n_discordant
                qc_info["strand_af_tolerance"] = tolerance
                qc_info["strand_af_discordance_action"] = non_palindromic_action
                qc_info["strand_palindromic_af_comparable"] = n_pal_comparable
                qc_info["strand_palindromic_af_discordant"] = n_pal_discordant
                qc_info["strand_palindromic_af_discordance_action"] = (
                    palindromic_action
                )
                qc_info["strand_non_palindromic_af_comparable"] = (
                    n_non_pal_comparable
                )
                qc_info["strand_non_palindromic_af_discordant"] = (
                    n_non_pal_discordant
                )
                if n_pal_discordant and palindromic_action == "fail":
                    raise AlleleFrequencyError(
                        "%s palindromic variants differ from the aligned population "
                        "reference AF by more than %s after resolved study orientation, "
                        "and strand.palindromic_af_discordance_action is 'fail'."
                        % (_fmt_int(n_pal_discordant), tolerance)
                    )
                if n_non_pal_discordant and non_palindromic_action == "fail":
                    raise AlleleFrequencyError(
                        "%s non-palindromic variants differ from the aligned population "
                        "reference AF by more than %s and "
                        "strand.af_discordance_action is 'fail'."
                        % (_fmt_int(n_non_pal_discordant), tolerance)
                    )
                if palindromic_action == "reject":
                    final_df, _ = reject_rows(
                        final_df,
                        palindromic_discordant,
                        step_label=STEP_LABEL,
                        reason="af_discordant",
                        context=emit,
                        collector=rejects,
                        detail=(
                            "palindromic variant has absolute EAF difference above %s "
                            "after resolved study orientation" % tolerance
                        ),
                        check_name="palindromic study/reference frequency concordance",
                        description=(
                            "After automatic consensus/internal study EAF or an explicit "
                            "reference-aligned declaration selected the orientation, "
                            "the palindromic variant's aligned study EAF differs from the "
                            "configured population reference by more than %s." % tolerance
                        ),
                        warn_on_remove=True,
                        warn_without_collector=True,
                        record_empty=True,
                    )
                else:
                    emit.qc(
                        "palindromic study/reference frequency concordance",
                        "%s of %s comparable variants differ from the aligned population "
                        "reference AF by more than %s; policy action is '%s'."
                        % (
                            _fmt_int(n_pal_discordant), _fmt_int(n_pal_comparable),
                            tolerance, palindromic_action,
                        ),
                        final_df.height,
                        final_df.height,
                        matched=n_pal_discordant,
                        outcome="retained_without_frequency_change",
                        warn=n_pal_discordant > 0,
                        step=STEP_LABEL,
                    )
                if non_palindromic_action == "reject":
                    final_df, _ = reject_rows(
                        final_df,
                        non_palindromic_discordant,
                        step_label=STEP_LABEL,
                        reason="af_discordant",
                        context=emit,
                        collector=rejects,
                        detail="absolute EAF difference above %s" % tolerance,
                        check_name="non-palindromic study/reference frequency concordance",
                        description=(
                            "The aligned study EAF differs from the configured population "
                            "reference by more than %s." % tolerance
                        ),
                        warn_on_remove=True,
                        warn_without_collector=True,
                        record_empty=True,
                    )
                else:
                    emit.qc(
                        "non-palindromic study/reference frequency concordance",
                        "%s of %s comparable variants differ from the aligned population "
                        "reference AF by more than %s; policy action is '%s'."
                        % (
                            _fmt_int(n_non_pal_discordant),
                            _fmt_int(n_non_pal_comparable), tolerance,
                            non_palindromic_action,
                        ),
                        final_df.height,
                        final_df.height,
                        matched=n_non_pal_discordant,
                        outcome="retained_without_frequency_change",
                        warn=n_non_pal_discordant > 0,
                        step=STEP_LABEL,
                    )
            else:
                qc_info["strand_af_comparable"] = 0
                qc_info["strand_af_discordant"] = 0
                qc_info["strand_palindromic_af_comparable"] = 0
                qc_info["strand_palindromic_af_discordant"] = 0
                qc_info["strand_non_palindromic_af_comparable"] = 0
                qc_info["strand_non_palindromic_af_discordant"] = 0

            # ---- frequencies of exactly 0 or 1 ----------------------------
            # These make either configured step-06 reconstruction denominator
            # exactly zero, which yields inf rather than null without its
            # reconstruction-specific guard.
            # fill_null(False): a missing frequency is not a degenerate one
            degenerate_expr = (
                (pl.col(final_eaf_col) == 0.0) | (pl.col(final_eaf_col) == 1.0)
            ).fill_null(False)
            n_degenerate = int(final_df.select(degenerate_expr.sum()).item() or 0)
            qc_info["eaf_degenerate_count"] = n_degenerate
            qc_info["eaf_degenerate_action"] = degenerate_action
            if degenerate_action == "fail" and n_degenerate:
                raise AlleleFrequencyError(
                    "%s variants have a frequency of exactly 0 or 1 and policy eaf.degenerate is 'fail'."
                    % _fmt_int(n_degenerate)
                )
            if degenerate_action == "reject":
                before_degenerate = final_df.height
                final_df, _ = reject_rows(
                    final_df,
                    degenerate_expr,
                    step_label=STEP_LABEL,
                    reason="eaf_degenerate",
                    context=emit,
                    collector=rejects,
                    detail="%s is exactly 0 or 1" % final_eaf_col,
                    check_name="degenerate frequencies",
                    description=(
                        "A frequency of exactly 0 or 1 makes the standard-error formula divide by zero, "
                        "so those variants were removed."
                    ),
                    warn_on_remove=True,
                    warn_without_collector=True,
                    record_empty=True,
                )
                qc_info["eaf_degenerate_removed"] = before_degenerate - final_df.height
            elif degenerate_action == "null":
                before_degenerate = final_df.height
                final_df = final_df.with_columns(
                    pl.when(degenerate_expr)
                    .then(None)
                    .otherwise(pl.col(final_eaf_col))
                    .alias(final_eaf_col)
                )
                emit.qc(
                    "degenerate frequencies",
                    "%s variants have a frequency of exactly 0 or 1, which makes the standard-error "
                    "formula divide by zero; policy eaf.degenerate is 'null', so the frequency was "
                    "blanked and the variants kept." % _fmt_int(n_degenerate),
                    before_degenerate,
                    final_df.height,
                    changed=n_degenerate,
                    warn=n_degenerate > 0,
                    step=STEP_LABEL,
                )
                qc_info["eaf_degenerate_removed"] = 0
            elif degenerate_action == "keep":
                emit.qc(
                    "degenerate frequencies",
                    "%s variants have a frequency of exactly 0 or 1; policy eaf.degenerate is 'keep', so "
                    "they were left in place." % _fmt_int(n_degenerate),
                    final_df.height,
                    final_df.height,
                    warn=n_degenerate > 0,
                    step=STEP_LABEL,
                )
                qc_info["eaf_degenerate_removed"] = 0
            else:  # "fail" with no degenerate values
                emit.qc(
                    "degenerate frequencies",
                    "No frequencies of exactly 0 or 1 were found; policy "
                    "eaf.degenerate is 'fail', so no action was needed.",
                    final_df.height,
                    final_df.height,
                    changed=0,
                    step=STEP_LABEL,
                )
                qc_info["eaf_degenerate_removed"] = 0

            final_stats = _summarise_frequency_after_changes(
                final_df, final_eaf_col,
                clip_tolerance=float(pol.get("eaf.clip_tolerance")),
                previous_values=validated_values,
                previous_statistics=validated_statistics,
            )
            _add_prefixed_statistics(qc_info, "final_eaf", final_stats)
            if replacement_selected:
                qc_info["eaf_replacement_rows"] = int(final_stats["non_null_rows"])

            sample_column_dict["eaf_col"] = final_eaf_col
            qc_info["final_eaf_col"] = final_eaf_col
            qc_info["final_variants"] = final_df.height
            qc_info["variants_removed"] = qc_info["Total_number_of_variants"] - final_df.height
            qc_info["final_status"] = "success"

            emit.info("Frequency column used          : %s" % final_eaf_col)
            emit.info("Initial variants               : %s" % _fmt_int(qc_info["Total_number_of_variants"]))
            emit.info("Variants with no frequency     : %s" % _fmt_int(missing_eaf_count))
            emit.info("Frequencies outside 0 to 1     : %s" % _fmt_int(out_of_bound_count))
            emit.info("Final variants                 : %s" % _fmt_int(qc_info["final_variants"]))
            emit.info("Variants removed               : %s" % _fmt_int(qc_info["variants_removed"]))
            if ctx is not None:
                ctx.set_rows(final_df.height)

        except Exception as exc:
            qc_info["final_status"] = "failed"
            qc_info["error_message"] = str(exc)
            raise

    return final_df, qc_info, sample_column_dict
