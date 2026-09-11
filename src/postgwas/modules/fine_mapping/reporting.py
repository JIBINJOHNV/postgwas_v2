"""Scientific HTML reporting for authoritative fine-mapping results.

The report is presentation-only. It reads the already validated, post-overlap
credible-set handoff and the run's QC tables; it never refits a model, changes
set membership, filters a scientific result, or recalculates PIP.
"""

from __future__ import annotations

from html import escape
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd

from postgwas.core.io.reports import collect_html_reports, read_static_report_content, write_html_report
from postgwas.modules.fine_mapping.overlap_resolution import _read_index


class FineMappingReportError(ValueError):
    """A published fine-mapping artifact cannot support a valid report."""


def _read_table(
    source: str | Path | None,
    *,
    delimiter: str,
    required: bool = False,
) -> pd.DataFrame:
    if source is None:
        if required:
            raise FineMappingReportError("Required fine-mapping report input is missing")
        return pd.DataFrame()
    path = Path(source)
    if not path.is_file():
        if required:
            raise FineMappingReportError(
                "Required fine-mapping report input does not exist: %s" % path
            )
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep=delimiter, dtype=object)
    except (OSError, UnicodeError, ValueError, pd.errors.ParserError) as exc:
        raise FineMappingReportError(
            "Cannot parse fine-mapping report input %s: %s" % (path, exc)
        ) from exc


def _clean(value: Any) -> Any:
    if value is None or value is pd.NA:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    return value


def _records(frame: pd.DataFrame, configured_columns: Sequence[str]) -> dict[str, Any]:
    columns = [column for column in configured_columns if column in frame.columns]
    if len(frame.columns) and not columns:
        raise FineMappingReportError(
            "Configured fine-mapping report columns do not match the retained "
            "table columns: configured=%s; available=%s"
            % (
                ", ".join(configured_columns),
                ", ".join(str(column) for column in frame.columns),
            )
        )
    return {
        "columns": columns,
        "rows": [
            {column: _clean(row.get(column)) for column in columns}
            for row in frame.to_dict(orient="records")
        ],
        "total_rows": int(len(frame)),
    }


def _reason_counts(values: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    return [
        {"reason": str(reason), "loci": int(count)}
        for reason, count in sorted((values or {}).items())
        if int(count)
    ]


def _probability(value: Any, label: str, source: Path) -> float:
    try:
        converted = float(value)
    except (OverflowError, TypeError, ValueError) as exc:
        raise FineMappingReportError(
            "%s contains a non-numeric %s: %r" % (source, label, value)
        ) from exc
    if not math.isfinite(converted) or not 0 <= converted <= 1:
        raise FineMappingReportError(
            "%s contains %s outside [0, 1]: %r" % (source, label, value)
        )
    return converted


def _integer(value: Any, label: str, source: Path) -> int:
    try:
        converted = int(float(value))
    except (OverflowError, TypeError, ValueError) as exc:
        raise FineMappingReportError(
            "%s contains an invalid %s: %r" % (source, label, value)
        ) from exc
    if converted < 1 or float(value) != converted:
        raise FineMappingReportError(
            "%s contains an invalid %s: %r" % (source, label, value)
        )
    return converted


def _credible_set_evidence(
    final_sets: pd.DataFrame,
    *,
    report_config: Mapping[str, Any],
    target_coverage: float,
    coverage_tolerance: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Read exact final set files and return metadata plus all member variants."""
    if final_sets.empty:
        return [], []
    required_metadata = {
        "final_credible_set_file",
        "final_genomic_locus",
        "analysis_round",
        "credible_set",
    }
    missing_metadata = sorted(required_metadata - set(final_sets.columns))
    if missing_metadata:
        raise FineMappingReportError(
            "Final combined credible-set table is missing report columns: %s"
            % ", ".join(missing_metadata)
        )
    if final_sets["final_credible_set_file"].duplicated().any():
        raise FineMappingReportError(
            "Final combined credible-set table repeats a final credible-set file"
        )

    file_columns = report_config["credible_set_file_columns"]
    rank_column = file_columns["rank"]
    variant_column = file_columns["variant_id"]
    probability_column = file_columns["posterior_probability"]
    required_file_columns = {rank_column, variant_column, probability_column}
    set_records: list[dict[str, Any]] = []
    variant_records: list[dict[str, Any]] = []

    for set_number, metadata in enumerate(
        final_sets.to_dict(orient="records"), start=1,
    ):
        source = Path(str(metadata["final_credible_set_file"]))
        frame = _read_table(source, delimiter=r"\s+", required=True)
        missing = sorted(required_file_columns - set(frame.columns))
        if missing:
            raise FineMappingReportError(
                "%s is missing configured credible-set columns: %s"
                % (source, ", ".join(missing))
            )
        if frame.empty:
            raise FineMappingReportError(
                "Authoritative credible-set file contains no variants: %s" % source
            )
        if frame[variant_column].astype(str).duplicated().any():
            raise FineMappingReportError(
                "Authoritative credible-set file repeats a variant: %s" % source
            )

        ranks = [
            _integer(value, rank_column, source) for value in frame[rank_column]
        ]
        if len(set(ranks)) != len(ranks):
            raise FineMappingReportError(
                "Authoritative credible-set file repeats a rank: %s" % source
            )
        if sorted(ranks) != list(range(1, len(frame) + 1)):
            raise FineMappingReportError(
                "Authoritative credible-set file ranks must be consecutive from 1: %s"
                % source
            )
        probabilities = [
            _probability(value, probability_column, source)
            for value in frame[probability_column]
        ]
        probability_mass = float(sum(probabilities))
        if probability_mass < target_coverage - coverage_tolerance:
            raise FineMappingReportError(
                "%s has posterior probability mass %.6g below configured coverage %.6g"
                % (source, probability_mass, target_coverage)
            )
        reported_variants = _clean(metadata.get("n_variants"))
        if reported_variants is not None:
            manifest_variants = _integer(reported_variants, "n_variants", source)
            if manifest_variants != len(frame):
                raise FineMappingReportError(
                    "%s contains %d variants but its manifest reports %s"
                    % (source, len(frame), reported_variants)
                )
        reported_target = _clean(metadata.get("target_coverage"))
        if reported_target is not None:
            manifest_target = _probability(
                reported_target, "target_coverage", source,
            )
            if abs(manifest_target - target_coverage) > coverage_tolerance:
                raise FineMappingReportError(
                    "%s reports target coverage %.12g, inconsistent with the "
                    "resolved value %.12g"
                    % (source, manifest_target, target_coverage)
                )
        model_probability = _clean(metadata.get("model_posterior_probability"))
        if model_probability is not None:
            _probability(
                model_probability, "model_posterior_probability", source,
            )

        locus = str(metadata["final_genomic_locus"])
        analysis_round = str(metadata["analysis_round"])
        credible_set = str(metadata["credible_set"])
        members = [
            {
                "final_genomic_locus": locus,
                "analysis_round": analysis_round,
                "credible_set": credible_set,
                "rank": rank,
                "variant_id": str(variant),
                "posterior_probability": probability,
                "warning_reason": _clean(metadata.get("warning_reason")),
                "final_credible_set_file": str(source),
            }
            for rank, variant, probability in zip(
                ranks,
                frame[variant_column].astype(str),
                probabilities,
            )
        ]
        members.sort(
            key=lambda row: (-row["posterior_probability"], row["rank"])
        )
        variant_records.extend(members)
        set_records.append({
            "report_set_number": set_number,
            "final_genomic_locus": locus,
            "analysis_round": analysis_round,
            "credible_set": credible_set,
            "n_variants": len(frame),
            "posterior_probability_mass": probability_mass,
            "top_variant": members[0]["variant_id"],
            "top_variant_probability": members[0]["posterior_probability"],
            "warning_reason": _clean(metadata.get("warning_reason")),
            "source_primary_genomic_loci": _clean(
                metadata.get("source_primary_genomic_loci")
            ),
            "final_selection_reason": _clean(metadata.get("final_selection_reason")),
            "target_coverage": _clean(metadata.get("target_coverage")),
            "achieved_component_coverage": _clean(
                metadata.get("achieved_component_coverage")
            ),
            "achieved_coverage": _clean(metadata.get("achieved_coverage")),
            "global_pip_sum": _clean(metadata.get("global_pip_sum")),
            "component_log10bf": _clean(metadata.get("component_log10bf")),
            "min_abs_ld": _clean(metadata.get("min_abs_ld")),
            "mean_abs_ld": _clean(metadata.get("mean_abs_ld")),
            "median_abs_ld": _clean(metadata.get("median_abs_ld")),
            "model_k": _clean(metadata.get("model_k")),
            "model_posterior_probability": _clean(
                metadata.get("model_posterior_probability")
            ),
            "probability_definition": _clean(metadata.get("probability_definition")),
            "credible_set_membership_definition": _clean(
                metadata.get("credible_set_membership_definition")
            ),
            "final_credible_set_file": str(source),
        })
    return set_records, variant_records


def _parameter_rows(args) -> list[dict[str, Any]]:
    configuration = args.resolved_fine_mapping_configuration
    rows = [
        {"Parameter": "Engine", "Value": configuration["engine"]},
        {"Parameter": "Genome build", "Value": configuration["genome_build"]},
        {"Parameter": "Locus definition", "Value": configuration["locus_type"]},
        {"Parameter": "Locus flank (kb)", "Value": configuration["locus_window_kb"]},
        {"Parameter": "Locus LP threshold", "Value": configuration["locus_lp_threshold"]},
        {
            "Parameter": "Credible-set coverage target",
            "Value": configuration["credible_set_coverage"],
        },
        {
            "Parameter": "MHC policy",
            "Value": (
                "exclude chr%s:%s-%s"
                % (
                    configuration["mhc_chromosome"],
                    configuration["mhc_start"],
                    configuration["mhc_end"],
                )
                if configuration["skip_mhc"] else "include"
            ),
        },
        {
            "Parameter": "Per-locus sample size",
            "Value": "%s (%s policy)"
            % (
                configuration["sample_size"]["summary_statistic"],
                configuration["sample_size"]["policy"],
            ),
        },
        {
            "Parameter": "NEF relative-range warning threshold",
            "Value": configuration["sample_size"][
                "relative_range_warning_threshold"
            ],
        },
        {
            "Parameter": "Maximum variants per locus",
            "Value": configuration["ld_resource_guard"][
                "maximum_variants_per_locus"
            ],
        },
        {
            "Parameter": "Post-primary overlap policy",
            "Value": configuration["overlap_resolution"]["policy"],
        },
    ]
    if configuration["engine"] == "susie":
        engine = configuration["engines"]["susie"]
        rows.extend((
            {
                "Parameter": "Maximum causal components (L)",
                "Value": engine["max_causal_components"],
            },
            {"Parameter": "Minimum credible-set purity", "Value": engine["minimum_purity"]},
            {"Parameter": "Primary maximum iterations", "Value": engine["fitting"]["main_max_iter"]},
            {"Parameter": "SuSiE fitting timeout (seconds)", "Value": engine["execution"]["susie_timeout_seconds"]},
        ))
    else:
        engine = configuration["engines"]["finemap"]
        rows.extend((
            {"Parameter": "FINEMAP algorithm", "Value": engine["algorithm"]},
            {"Parameter": "Maximum causal SNPs", "Value": engine["n_causal_snps"]},
            {"Parameter": "Prior effect-size SD", "Value": engine["prior_std"]},
            {"Parameter": "Iterations", "Value": engine["n_iter"]},
            {"Parameter": "FINEMAP timeout (seconds)", "Value": engine["finemap_timeout_seconds"]},
        ))
    return rows


def _initial_result_summary(result: Mapping[str, Any]) -> dict[str, int]:
    """Count initial sets from the retained index, before any overlap reruns."""
    # Overlap resolution preserves flames_index and n_credible_sets from the
    # engine's initial pass; flames_input and n_final_credible_sets refer to the
    # later output and must not be used to reconstruct initial counts.
    expected = result.get("n_credible_sets")
    if expected is None:
        raise FineMappingReportError("Initial credible-set count is missing")
    initial_index = _read_index(dict(result))
    if len(initial_index) != int(expected):
        raise FineMappingReportError(
            "Initial credible-set index does not match the validated result count"
        )
    if initial_index["GenomicLocus"].isna().any() or initial_index["Filename"].duplicated().any():
        raise FineMappingReportError("Initial credible-set index has missing loci or duplicate files")
    loci_with_sets = int(initial_index["GenomicLocus"].nunique())
    completed = int(result["n_successful"])
    if loci_with_sets > completed:
        raise FineMappingReportError("Initial loci with credible sets exceed completed loci")
    return {
        "initial_credible_sets": len(initial_index),
        "initial_loci_with_credible_sets": loci_with_sets,
        "initial_loci_without_credible_sets": completed - loci_with_sets,
    }


def _qc_range(
    frame: pd.DataFrame, columns: Sequence[str], digits: int, *, percent: bool = False,
) -> str:
    """Describe observed QC values without treating absent diagnostics as zero."""
    values = pd.to_numeric(
        pd.Series(frame.reindex(columns=columns).to_numpy().ravel()), errors="coerce",
    )
    if values.empty or not values.map(math.isfinite).all():
        return "Not fully recorded in locus QC"
    scale = 100 if percent else 1
    low, high = float(values.min()) * scale, float(values.max()) * scale
    lower = _value(int(low) if low.is_integer() else low, significant_digits=digits)
    upper = _value(int(high) if high.is_integer() else high, significant_digits=digits)
    return (lower if low == high else lower + "–" + upper) + ("%" if percent else "")


def _warning_explanations(
    locus_status: pd.DataFrame, configuration: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Explain recorded initial-locus warnings; do not assign new QC outcomes."""
    if "warning_reason" not in locus_status:
        return []
    # Semicolons separate engine warning codes in the retained QC contract.
    tokens = locus_status["warning_reason"].fillna("").astype(str).str.split(";")
    digits = int(configuration["html_report"]["probability_significant_digits"])
    explanations = []
    for code in ("nef_relative_range_exceeds_threshold", "ld_z_mismatch_warning"):
        affected = locus_status.loc[tokens.map(lambda values: code in values)]
        if affected.empty:
            continue
        if code == "nef_relative_range_exceeds_threshold":
            settings = configuration["sample_size"]
            explanation = {
                "title": "Effective sample size differs between variants",
                "meaning": (
                    "NEF is the effective sample size supplied for each variant. "
                    "This warning means that (maximum NEF − minimum NEF) / median NEF "
                    "within a locus exceeds the configured threshold. Different "
                    "variants may have been measured in different numbers of people "
                    "or contributing studies; this diagnostic alone does not establish why."
                ),
                "impact": (
                    "The fit uses one sample-size value for the whole locus. When "
                    "variants have different sample sizes or contributing samples, "
                    "that approximation can affect posterior probabilities and "
                    "credible-set coverage. The relative range does not reveal how "
                    "many variants have low NEF or which participants overlap."
                ),
                "action": (
                    "Configured policy: %s; locus sample-size summary: %s. This "
                    "policy permits fitting with the warning recorded; it does not "
                    "correct sample-size heterogeneity. Check the per-variant NEF "
                    "distribution and study contributions before interpreting the sets."
                ) % (settings["policy"], settings["summary_statistic"]),
                "measurements": [
                    {"Measure": "Within-locus relative range (across affected loci)", "Value": _qc_range(affected, ["nef_relative_range"], digits, percent=True)},
                    {"Measure": "Configured warning threshold", "Value": _value(float(settings["relative_range_warning_threshold"]) * 100, significant_digits=digits) + "%"},
                    {"Measure": "Per-variant NEF range across affected loci", "Value": _qc_range(affected, ["nef_min", "nef_max"], digits)},
                    {"Measure": "Selected locus sample sizes", "Value": _qc_range(affected, ["selected_nef"], digits)},
                ],
                "reference": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9337707/",
                "reference_label": "SuSiE-RSS paper: sample consistency and reference LD",
            }
        elif configuration["engine"] == "susie":
            settings = configuration["engines"]["susie"]["ld_validation"]
            explanation = {
                "title": "GWAS statistics and reference LD do not fully agree",
                "meaning": (
                    "LD describes correlation between variants. SuSiE compares the "
                    "GWAS z-scores with the reference LD using estimate_s_rss. A "
                    "larger diagnostic λ indicates greater disagreement. It is not "
                    "genomic inflation λ, an error rate or the percentage of incorrect variants."
                ),
                "impact": (
                    "Reference-panel sampling variation, population differences, "
                    "allele-coding errors or inconsistent GWAS samples can contribute. "
                    "The warning does not identify the cause. Inaccurate LD can change "
                    "variant probabilities and credible sets, including producing false signals."
                ),
                "action": (
                    "PostGWAS permits fitting above the warning threshold but fails "
                    "the locus above the failure threshold. These are configured "
                    "decision limits, not guarantees of accuracy. The warning itself "
                    "does not correct the LD or GWAS data. Check allele alignment, "
                    "reference suitability and variant-level z-score/LD diagnostics."
                ),
                "measurements": [
                    {"Measure": "Diagnostic λ across affected loci", "Value": _qc_range(affected, ["ld_z_mismatch_lambda"], digits)},
                    {"Measure": "Configured warning threshold (λ >)", "Value": settings["mismatch_warning_threshold"]},
                    {"Measure": "Configured failure threshold (λ >)", "Value": settings["mismatch_failure_threshold"]},
                ],
                "reference": "https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html",
                "reference_label": "Official SuSiE LD/statistics diagnostic",
            }
        else:
            continue
        explanation.update(code=code, loci=int(affected["genomic_locus"].nunique()))
        explanations.append(explanation)
    return explanations


def build_fine_mapping_report(
    args,
    preflight,
    result: Mapping[str, Any],
    output_paths: Mapping[str, Path],
    *,
    upstream_results: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build one report model from final result and QC artifacts."""
    report_config = dict(args.fine_mapping_html_report)
    final_path = result.get("final_combined_credible_sets")
    final_sets = _read_table(final_path, delimiter="\t", required=True)
    expected_sets = int(result.get("n_final_credible_sets", 0) or 0)
    if len(final_sets) != expected_sets:
        raise FineMappingReportError(
            "Final combined credible-set row count does not match the validated "
            "result summary: %d != %d" % (len(final_sets), expected_sets)
        )
    set_records, variant_records = _credible_set_evidence(
        final_sets,
        report_config=report_config,
        target_coverage=float(args.resolved_fine_mapping_configuration[
            "credible_set_coverage"
        ]),
        coverage_tolerance=float(args.fine_mapping_validation[
            "credible_set_coverage_tolerance"
        ]),
    )
    locus_status = _read_table(
        result.get("locus_status"), delimiter="\t", required=True,
    )
    preflight_table = _read_table(
        result.get("preflight_validation"), delimiter="\t", required=True,
    )
    overlap = _read_table(
        result.get("overlap_resolution"), delimiter="\t", required=True,
    )
    if args.finemap_method == "susie":
        diagnostics_path = result.get("recovery_audit")
        diagnostics = _read_table(
            diagnostics_path, delimiter="\t", required=True,
        )
        diagnostic_columns = report_config["susie_recovery_columns"]
        diagnostic_label = "SuSiE fitting and recovery audit"
    else:
        diagnostics_path = output_paths["finemap_model_summary_file"]
        diagnostics = _read_table(
            diagnostics_path, delimiter=",", required=True,
        )
        diagnostic_columns = report_config["finemap_model_columns"]
        diagnostic_label = "FINEMAP causal-count model selection"

    authoritative_loci = sorted({
        record["final_genomic_locus"] for record in set_records
    })
    probability_values = [
        record["posterior_probability"] for record in variant_records
    ]
    input_matches = int(preflight.reference.locus_variant_id_matches)
    input_requested = int(preflight.reference.locus_variant_ids_requested)
    concordance = (
        100.0 * input_matches / input_requested if input_requested else None
    )
    report_path = Path(output_paths["html_report_file"])
    candidate_output_links = [
        {"label": "HTML report", "path": str(output_paths["html_report_file"])},
        {"label": "Final combined credible sets", "path": final_path},
        {"label": "Initial credible-set files", "path": result.get("flames_index")},
        {"label": "Final files for FLAMES", "path": result.get("flames_input")},
        {"label": "Input and resource validation", "path": result.get("preflight_validation")},
        {"label": "Primary locus QC", "path": result.get("locus_status")},
        {"label": "Overlap-resolution audit", "path": result.get("overlap_resolution")},
        {"label": "Overlap-resolution summary", "path": result.get("overlap_resolution_summary")},
        {"label": diagnostic_label, "path": diagnostics_path},
        {"label": "Resolved run configuration", "path": output_paths["run_configuration_file"]},
        {"label": "Detailed pipeline log", "path": output_paths["pipeline_log_file"]},
        {"label": "Diagnostic plots", "path": output_paths["diagnostic_plots_directory"]},
    ]
    output_links = []
    for row in candidate_output_links:
        raw_path = row["path"]
        if not raw_path:
            continue
        candidate = Path(raw_path)
        if candidate.is_dir() and not any(candidate.iterdir()):
            continue
        if candidate == report_path or candidate.exists():
            output_links.append({"label": row["label"], "path": str(candidate)})

    return {
        "schema_version": 1,
        "dataset_id": str(args.dataset_id),
        "engine": str(args.finemap_method),
        "engine_label": "SuSiE-RSS" if args.finemap_method == "susie" else "FINEMAP",
        "status": str(result.get("status", "unknown")),
        "genome_build": str(args.genome_build),
        "report_path": str(output_paths["html_report_file"]),
        "config": report_config,
        "summary": {
            **_initial_result_summary(result),
            "input_variants": int(preflight.inputs.summary_statistics_rows),
            "variants_in_eligible_loci": int(preflight.inputs.locus_summary_variants),
            "input_loci": int(preflight.inputs.input_loci),
            "eligible_loci": int(preflight.inputs.eligible_loci),
            "primary_loci_attempted": int(result.get("n_attempted", 0) or 0),
            "primary_loci_successful": int(result.get("n_successful", 0) or 0),
            "primary_loci_failed": int(result.get("n_failed", 0) or 0),
            "loci_with_warnings": int(result.get("n_warnings", 0) or 0),
            "authoritative_loci": len(authoritative_loci),
            "final_credible_sets": len(set_records),
            "credible_set_variant_memberships": len(variant_records),
            "maximum_posterior_probability": (
                max(probability_values) if probability_values else None
            ),
            "overlap_groups": int(result.get("n_overlap_groups", 0) or 0),
            "successful_joint_reruns": int(
                result.get("n_joint_rerun_successful", 0) or 0
            ),
            "excluded_or_failed_joint_reruns": int(
                result.get("n_joint_rerun_failed", 0) or 0
            ),
            "joint_reruns_without_sets": int(
                result.get("n_joint_rerun_without_credible_sets", 0) or 0
            ),
            "reference_variant_matches": input_matches,
            "coordinate_concordant_variant_matches": int(
                preflight.reference.coordinate_concordant_variant_matches
            ),
            "locus_variants_requested": input_requested,
            "loci_with_reference_variants": int(
                preflight.reference.loci_with_reference_variants
            ),
            "reference_match_percent": concordance,
        },
        "warning_reasons": _reason_counts(result.get("warning_reason_counts")),
        "warning_explanations": _warning_explanations(
            locus_status, args.resolved_fine_mapping_configuration,
        ),
        "failure_reasons": _reason_counts(result.get("failure_reason_counts")),
        "credible_sets": set_records,
        "credible_set_table": _records(
            final_sets, report_config["credible_set_columns"],
        ),
        "variant_table": {
            "columns": [
                "final_genomic_locus",
                "analysis_round",
                "credible_set",
                "rank",
                "variant_id",
                "posterior_probability",
                "warning_reason",
            ],
            "rows": variant_records,
            "total_rows": len(variant_records),
        },
        "locus_status": _records(
            locus_status, report_config["locus_status_columns"],
        ),
        "overlap": _records(overlap, report_config["overlap_columns"]),
        "preflight": _records(
            preflight_table, report_config["preflight_columns"],
        ),
        "diagnostics": _records(diagnostics, diagnostic_columns),
        "diagnostic_label": diagnostic_label,
        "parameters": _parameter_rows(args),
        "inputs": [
            {"Input": "Summary statistics", "Path": str(preflight.inputs.summary_statistics_file)},
            {"Input": "Locus definitions", "Path": str(preflight.inputs.locus_file)},
            {"Input": "PLINK BED", "Path": str(preflight.reference.bed_file)},
            {"Input": "PLINK BIM", "Path": str(preflight.reference.bim_file)},
            {"Input": "PLINK FAM", "Path": str(preflight.reference.fam_file)},
        ],
        "tools": [
            {"Runtime": tool.name, "Version": tool.version, "Path": tool.path}
            for tool in preflight.tools
        ],
        "outputs": output_links,
        "upstream_reports": collect_html_reports(upstream_results),
    }


def _value(value: Any, *, significant_digits: int = 4) -> str:
    value = _clean(value)
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return format(value, ".%dg" % significant_digits)
    return str(value)


def _title(column: str) -> str:
    replacements = {
        "ld": "LD",
        "nef": "NEF",
        "pip": "PIP",
        "bf": "BF",
        "chr": "chromosome",
    }
    return " ".join(
        replacements.get(part.lower(), part).upper()
        if part.lower() in {"ld", "nef", "pip", "bf"}
        else replacements.get(part.lower(), part).capitalize()
        for part in str(column).split("_")
    )


def _static_table(
    rows: Sequence[Mapping[str, Any]],
    columns: Sequence[str],
    *,
    path_columns: Sequence[str] = (),
    report_path: Path | None = None,
    significant_digits: int = 4,
) -> str:
    if not rows or not columns:
        return '<p class="empty">No records were produced for this section.</p>'
    header = "".join("<th>%s</th>" % escape(_title(column)) for column in columns)
    body_rows = []
    for row in rows:
        cells = []
        for column in columns:
            raw = _clean(row.get(column))
            rendered = escape(_value(raw, significant_digits=significant_digits))
            if raw is not None and column in path_columns and report_path is not None:
                href = os.path.relpath(str(raw), report_path.parent)
                rendered = '<a href="%s">%s</a>' % (
                    escape(href, quote=True), rendered,
                )
            cells.append("<td>%s</td>" % rendered)
        body_rows.append("<tr>%s</tr>" % "".join(cells))
    return (
        '<div class="table-wrap"><table><thead><tr>%s</tr></thead>'
        '<tbody>%s</tbody></table></div>' % (header, "".join(body_rows))
    )


def _interactive_table(
    table: Mapping[str, Any],
    *,
    table_id: str,
    page_size: int,
    significant_digits: int,
) -> str:
    columns = table.get("columns") or []
    rows = table.get("rows") or []
    if not columns or not rows:
        return '<p class="empty">No records were produced for this section.</p>'
    payload_rows = [
        [_value(row.get(column), significant_digits=significant_digits) for column in columns]
        for row in rows
    ]
    payload = json.dumps(
        payload_rows, ensure_ascii=False, separators=(",", ":")
    ).replace("</", "<\\/")
    header = "".join(
        '<th><button type="button" data-column="%d">%s</button></th>'
        % (index, escape(_title(column)))
        for index, column in enumerate(columns)
    )
    body = "".join(
        "<tr>%s</tr>" % "".join("<td>%s</td>" % escape(value) for value in row)
        for row in payload_rows
    )
    return """
<div class="data-browser" id="%s" data-page-size="%d" data-source="%s-data">
  <div class="table-tools" hidden>
    <label>Search all rows <input type="search" class="table-search" autocomplete="off" placeholder="locus, variant, status, warning…"></label>
    <div class="pager"><button type="button" class="previous">Previous</button><span class="page-state" aria-live="polite"></span><button type="button" class="next">Next</button></div>
  </div>
  <div class="table-wrap"><table><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>
  <noscript><p class="notice">All rows are shown. JavaScript adds search, sorting and pagination.</p></noscript>
</div>
<script type="application/json" id="%s-data">%s</script>""" % (
        escape(table_id, quote=True),
        page_size,
        escape(table_id, quote=True),
        header,
        body,
        escape(table_id, quote=True),
        payload,
    )


def _metric(label: str, value: Any, tone: str, significant_digits: int) -> str:
    return (
        '<div class="metric %s"><div class="metric-label">%s</div>'
        '<div class="metric-value">%s</div></div>'
        % (
            escape(tone),
            escape(label),
            escape(_value(value, significant_digits=significant_digits)),
        )
    )


def _reason_section(
    warnings: Sequence[Mapping[str, Any]],
    failures: Sequence[Mapping[str, Any]],
    explanations: Sequence[Mapping[str, Any]],
    significant_digits: int,
) -> str:
    rows = [
        {"Severity": "Warning", "Reason": row["reason"], "Loci": row["loci"]}
        for row in warnings
    ] + [
        {"Severity": "Failure", "Reason": row["reason"], "Loci": row["loci"]}
        for row in failures
    ]
    if not rows and not explanations:
        return (
            '<div class="notice good"><strong>No primary-locus warnings or failures '
            'were recorded.</strong> Continue to inspect the per-locus table and '
            'method assumptions before biological interpretation.</div>'
        )
    cards = []
    for item in explanations:
        cards.append(
            '<article class="concept"><h3>%s</h3><p><strong>%s initial loci affected.</strong></p>'
            '<p>%s</p>%s<h3>Why it matters</h3><p>%s</p>'
            '<h3>What the analysis did and what to review</h3><p>%s</p>'
            '<p><a href="%s">%s</a></p><p class="card-note">Log code: <code>%s</code></p></article>'
            % (
                escape(item["title"]), escape(_value(item["loci"])),
                escape(item["meaning"]), _static_table(
                    item["measurements"], ("Measure", "Value"),
                    significant_digits=significant_digits,
                ),
                escape(item["impact"]), escape(item["action"]),
                escape(item["reference"], quote=True), escape(item["reference_label"]),
                escape(item["code"]),
            )
        )
    return (
        '<div class="notice warning"><strong>Review required.</strong> A completed '
        'run can contain warned or failed loci. The explanations below describe '
        'recorded warnings from the initial analysis. A locus can have more than '
        'one warning; affected-locus counts should not be added together.</div>'
        + '<div class="concept-grid">' + "".join(cards) + '</div>'
        + '<details><summary>All recorded warning and failure codes</summary>'
        + _static_table(rows, ("Severity", "Reason", "Loci"))
        + '</details>'
    )


def _probability_bar(value: float) -> str:
    width = max(0.0, min(100.0, 100.0 * float(value)))
    return '<span class="prob-fill" style="width:%.6g%%"></span>' % width


def _credible_set_cards(
    report: Mapping[str, Any],
    *,
    significant_digits: int,
) -> str:
    records = report["credible_sets"]
    if not records:
        return (
            '<div class="empty-state"><strong>No credible set was retained.</strong>'
            '<p>This is a valid fine-mapping outcome when the fitted model has no set '
            'meeting the configured reporting criteria. It is not evidence that the '
            'locus contains no causal variant.</p></div>'
        )
    # A source file uniquely identifies a validated set, even when loci reuse
    # a component label or a variant belongs to more than one set.
    members_by_file: dict[str, list[dict[str, Any]]] = {}
    for member in report["variant_table"]["rows"]:
        members_by_file.setdefault(member["final_credible_set_file"], []).append(member)
    cards = []
    probability_label = "PIP" if report["engine"] == "susie" else "SNP posterior probability"
    for record, metadata in zip(records, report["credible_set_table"]["rows"]):
        members = members_by_file[record["final_credible_set_file"]]
        source_link = _static_table(
            [{"Source file": record["final_credible_set_file"]}], ("Source file",),
            path_columns=("Source file",), report_path=Path(report["report_path"]),
            significant_digits=significant_digits,
        )
        variants = []
        for member in members:
            details = _static_table(
                [{"Field": probability_label if column == "posterior_probability" else _title(column), "Value": member[column]}
                 for column in report["variant_table"]["columns"]],
                ("Field", "Value"), significant_digits=significant_digits,
            )
            variants.append(
                '<details class="variant-detail"><summary class="prob-row">'
                '<span class="variant">%s</span><span class="prob-track" aria-hidden="true">%s</span>'
                '<strong>%s <span>%s</span></strong></summary>%s<div class="output-table">%s</div></details>'
                % (escape(member["variant_id"]),
                   _probability_bar(member["posterior_probability"]),
                   escape(probability_label),
                   escape(_value(member["posterior_probability"], significant_digits=significant_digits)),
                   details, source_link)
            )
        round_label = {"primary": "Initial fine-mapping", "joint": "Additional fine-mapping"}.get(
            record["analysis_round"], record["analysis_round"],
        )
        warning = (
            '<span class="mini-alert">Warning: %s</span>' % escape(str(record["warning_reason"]))
            if record["warning_reason"] else ""
        )
        set_details = _static_table(
            [{"Field": _title(column), "Value": metadata.get(column)}
             for column in report["credible_set_table"]["columns"]],
            ("Field", "Value"), significant_digits=significant_digits,
        )
        cards.append("""
<details class="cs-card">
  <summary class="cs-summary"><span class="cs-heading"><span><span class="kicker">%s</span><strong class="set-title">%s · set %s</strong></span><span class="set-size">%s variants</span></span>%s<span class="card-note">Expand to explore variants and set details</span></summary>
  <div class="cs-content"><p class="provenance">Initially analysed loci: %s</p>
  <p class="section-note">All %s variants are listed, highest posterior probability first. Expand a variant for its recorded details. The source file below contains this set's variants.</p>
  <div class="prob-list">%s</div>
  <details class="set-details"><summary>Set statistics and source file</summary>%s<div class="output-table">%s</div></details></div>
</details>""" % (
            escape(round_label), escape(record["final_genomic_locus"]),
            escape(record["credible_set"]), escape(_value(record["n_variants"])), warning,
            escape(str(record.get("source_primary_genomic_loci") or "not reported")),
            escape(_value(record["n_variants"])), "".join(variants), set_details, source_link,
        ))
    return (
        '<div class="credible-set-browser"><div class="table-tools" hidden>'
        '<label>Search sets and variants <input type="search" class="set-search" '
        'placeholder="locus, set ID, variant or warning" autocomplete="off"></label>'
        '<div class="pager"><button type="button" class="expand-sets">Expand matching sets</button>'
        '<button type="button" class="collapse-sets">Collapse all</button></div></div>'
        '<p class="set-search-state" role="status" aria-live="polite"></p>'
        '<div class="cs-grid">%s</div></div>' % "".join(cards)
    )


_STYLES = """
:root{--ink:#152238;--muted:#617188;--canvas:#f3f7fa;--panel:#fff;--line:#dbe5ea;--navy:#102a43;--teal:#0f766e;--teal2:#14b8a6;--blue:#0369a1;--good:#166534;--goodbg:#ecfdf3;--amber:#92400e;--amberbg:#fff8e7;--red:#991b1b;--redbg:#fff1f2;--violet:#6d28d9}*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--canvas);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.hero{background:radial-gradient(circle at 86% 20%,rgba(45,212,191,.34),transparent 28%),linear-gradient(128deg,#0f2740,#075985 56%,#0f766e);color:#fff}.hero-inner{max-width:1440px;margin:auto;padding:42px 26px 40px}.eyebrow{margin:0 0 7px;font-size:12px;font-weight:850;letter-spacing:.16em;text-transform:uppercase;color:#99f6e4}.hero h1{margin:0;font-size:clamp(31px,4.5vw,52px);line-height:1.04;letter-spacing:-.035em}.hero-copy{max-width:920px;margin:14px 0 20px;color:#dbeafe;font-size:17px}.hero-meta{display:flex;gap:9px;flex-wrap:wrap}.badge{display:inline-flex;align-items:center;padding:6px 10px;border-radius:999px;font-size:12px;font-weight:850;letter-spacing:.035em;text-transform:uppercase}.badge.light{color:#083344;background:#ccfbf1}.badge.good{color:var(--good);background:#dcfce7}.badge.warning{color:var(--amber);background:#fef3c7}.badge.error{color:var(--red);background:#fee2e2}.report-nav{position:sticky;top:0;z-index:9;background:rgba(255,255,255,.94);backdrop-filter:blur(12px);border-bottom:1px solid var(--line)}.report-nav div{max-width:1440px;margin:auto;padding:10px 26px;display:flex;gap:18px;overflow:auto;white-space:nowrap}.report-nav a{color:#33536e;text-decoration:none;font-weight:750;font-size:13px}.report-nav a:hover{color:var(--teal)}main{max-width:1440px;margin:auto;padding:26px 24px 64px}.notice,section{background:var(--panel);border:1px solid var(--line);border-radius:14px;box-shadow:0 7px 28px rgba(15,23,42,.045)}.notice{margin:0 0 20px;padding:15px 17px;border-left:5px solid var(--blue)}.notice.good{border-left-color:#16a34a;background:var(--goodbg)}.notice.warning{border-left-color:#d97706;background:var(--amberbg)}.notice.error{border-left-color:#dc2626;background:var(--redbg)}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(185px,1fr));gap:12px;margin:0 0 22px}.metric{min-height:94px;padding:15px 16px;border:1px solid var(--line);border-top:4px solid var(--blue);border-radius:11px;background:#fff}.metric.good{border-top-color:#16a34a}.metric.warning{border-top-color:#d97706}.metric.analysis{border-top-color:var(--violet)}.metric-label{font-size:11px;font-weight:850;color:var(--muted);letter-spacing:.065em;text-transform:uppercase}.metric-value{margin-top:5px;font-size:23px;font-weight:800;overflow-wrap:anywhere}section{margin:19px 0;padding:23px}section h2{margin:0 0 5px;font-size:22px;letter-spacing:-.012em}section h3{margin:20px 0 7px;font-size:17px}.section-note{max-width:105ch;margin:0 0 17px;color:var(--muted)}.concept-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(285px,1fr));gap:13px}.concept{padding:17px;border:1px solid var(--line);border-radius:11px;background:#f8fafc}.concept h3{margin:0 0 6px}.concept p{margin:0;color:#475569}.cs-grid{display:grid;grid-template-columns:1fr;gap:14px}.cs-card{border:1px solid var(--line);border-radius:12px;padding:0;overflow:hidden;background:#fbfdff}.cs-heading{display:flex;justify-content:space-between;gap:13px;align-items:start}.kicker{color:var(--teal);font-size:10px;font-weight:900;letter-spacing:.1em;text-transform:uppercase}.set-size{white-space:nowrap;padding:4px 8px;border-radius:999px;background:#e0f2fe;color:#075985;font-size:11px;font-weight:800}.provenance,.card-note{color:var(--muted);font-size:12px}.mini-alert{margin:8px 0;padding:7px 9px;border-radius:7px;background:var(--amberbg);color:var(--amber);font-size:12px;overflow-wrap:anywhere}.prob-list{display:grid;gap:7px;margin:13px 0}.prob-row{display:grid;grid-template-columns:minmax(115px,1.25fr) minmax(110px,2fr) 62px;gap:8px;align-items:center;font-size:12px}.variant{white-space:normal;overflow-wrap:anywhere;font-family:ui-monospace,SFMono-Regular,Menlo,monospace}.prob-track{height:8px;border-radius:999px;background:#e2e8f0;overflow:hidden}.prob-fill{display:block;height:100%;border-radius:999px;background:linear-gradient(90deg,var(--teal),var(--teal2))}.prob-row strong{text-align:right;font-variant-numeric:tabular-nums}.table-tools{display:flex;align-items:end;justify-content:space-between;gap:12px;flex-wrap:wrap;margin:0 0 12px}.table-tools label{display:flex;flex-direction:column;gap:5px;font-size:11px;font-weight:850;color:var(--muted);letter-spacing:.05em;text-transform:uppercase}.table-tools input{min-width:min(440px,85vw);padding:9px 11px;border:1px solid #b8c6d2;border-radius:8px;background:#fff;color:var(--ink);font:inherit;text-transform:none;letter-spacing:normal}.pager{display:flex;align-items:center;gap:9px}.pager button{padding:8px 11px;border:1px solid var(--line);border-radius:8px;background:#fff;color:var(--blue);font-weight:750;cursor:pointer}.pager button:disabled{color:#94a3b8;cursor:not-allowed}.page-state{min-width:210px;text-align:center;color:var(--muted);font-size:13px}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:9px;max-height:68vh}table{width:100%;border-collapse:collapse;background:#fff;font-size:13px}th,td{padding:9px 11px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top;white-space:nowrap}thead th{position:sticky;top:0;z-index:1;background:#eaf5f6;color:#334155;font-size:10px;text-transform:uppercase;letter-spacing:.045em}th button{border:0;background:transparent;color:inherit;font:inherit;font-weight:850;text-transform:inherit;letter-spacing:inherit;padding:0;cursor:pointer}.empty,.card-note{margin-bottom:0}.empty-state{padding:22px;border:1px dashed #94a3b8;border-radius:11px;background:#f8fafc}.empty-state p{margin-bottom:0}.output-table td:last-child,.path{white-space:normal;overflow-wrap:anywhere;word-break:break-word}.output-table a,a{color:#0369a1}.callout-list{margin:9px 0 0;padding-left:20px}.callout-list li{margin:6px 0}.footer{max-width:1440px;margin:0 auto;padding:0 24px 34px;color:var(--muted);font-size:13px}@media(max-width:720px){.hero-inner{padding:32px 18px}.report-nav div{padding:9px 15px}main{padding:17px 10px 45px}section{padding:16px}.page-state{min-width:auto}.prob-row{grid-template-columns:minmax(90px,1fr) minmax(75px,1.6fr) 52px}}@media print{body{background:#fff}.hero{background:#fff;color:var(--ink)}.hero-copy{color:var(--muted)}.eyebrow{color:var(--teal)}.report-nav,.table-tools{display:none}main{max-width:none;padding:0}.notice,section{box-shadow:none;break-inside:avoid}.table-wrap{max-height:none;overflow:visible}a{color:inherit;text-decoration:none}}
[hidden]{display:none!important}.cs-summary{position:relative;padding:17px 17px 17px 40px;cursor:pointer;list-style:none}.cs-summary::-webkit-details-marker,.variant-detail>summary::-webkit-details-marker{display:none}.cs-summary:before,.variant-detail>summary:before{content:"+";position:absolute;left:16px;font-weight:800;color:var(--teal)}.cs-card[open]>.cs-summary:before,.variant-detail[open]>summary:before{content:"−"}.cs-summary .card-note,.cs-summary .mini-alert,.set-title{display:block}.set-title{font-size:17px;overflow-wrap:anywhere}.cs-content{padding:0 17px 17px}.cs-card[open]>.cs-summary{border-bottom:1px solid var(--line)}.cs-content .prob-list{max-height:65vh;overflow:auto;border:1px solid var(--line);border-radius:8px;gap:0}.variant-detail{border-bottom:1px solid var(--line);background:white}.variant-detail:last-child{border-bottom:0}.variant-detail>summary{position:relative;padding:12px 14px 12px 36px;cursor:pointer;list-style:none;grid-template-columns:minmax(150px,1.3fr) minmax(80px,2fr) minmax(90px,1fr)}.variant-detail>summary strong{font-size:11px}.variant-detail>summary strong span{display:block;font-size:14px}.variant-detail>.table-wrap{margin:0 14px 14px;border:0;max-height:none}.cs-content td{white-space:normal;overflow-wrap:anywhere}.cs-content .table-wrap{max-height:none}.set-details{margin-top:15px}.set-details>summary{cursor:pointer;font-weight:750;padding:8px 0}.set-search-state{color:var(--muted)}summary:focus-visible,button:focus-visible,input:focus-visible{outline:3px solid var(--blue);outline-offset:3px}@media(max-width:720px){.variant-detail>summary{grid-template-columns:minmax(95px,1.4fr) minmax(45px,1fr) minmax(70px,1fr);gap:7px}.cs-summary .cs-heading{flex-wrap:wrap}.cs-content{padding:0 10px 12px}}@media print{.cs-content .prob-list{max-height:none;overflow:visible}.cs-card{break-inside:auto}}

.upstream-detail>summary{cursor:pointer;font-weight:800;font-size:18px;padding:12px 0}.embedded-report{margin-top:15px}.embedded-report section{box-shadow:none}.embedded-report td{white-space:normal;overflow-wrap:anywhere}.embedded-report h1{font-size:23px}
"""


_TABLE_SCRIPT = """
<script>
document.querySelectorAll('.data-browser').forEach(function(container) {
  const payload = document.getElementById(container.dataset.source);
  const body = container.querySelector('tbody');
  const search = container.querySelector('.table-search');
  const previous = container.querySelector('.previous');
  const next = container.querySelector('.next');
  const state = container.querySelector('.page-state');
  const pageSize = Number(container.dataset.pageSize);
  container.querySelector('.table-tools').hidden = false;
  const rows = JSON.parse(payload.textContent);
  let filtered = rows.slice();
  let page = 0;
  let sortColumn = null;
  let ascending = true;
  function sortable(value) {
    const text = String(value).replaceAll(',', '').trim();
    if (text !== '' && Number.isFinite(Number(text))) return Number(text);
    return text.toLocaleLowerCase();
  }
  function render() {
    const pages = Math.max(1, Math.ceil(filtered.length / pageSize));
    page = Math.max(0, Math.min(page, pages - 1));
    const start = page * pageSize;
    const end = Math.min(start + pageSize, filtered.length);
    body.replaceChildren();
    filtered.slice(start, end).forEach(function(row) {
      const tr = document.createElement('tr');
      row.forEach(function(value) {
        const td = document.createElement('td');
        td.textContent = value;
        tr.appendChild(td);
      });
      body.appendChild(tr);
    });
    state.textContent = filtered.length
      ? `${start + 1}–${end} of ${filtered.length.toLocaleString()} · page ${page + 1} of ${pages}`
      : '0 matching rows';
    previous.disabled = page === 0;
    next.disabled = page + 1 >= pages;
  }
  search.addEventListener('input', function() {
    const term = search.value.trim().toLocaleLowerCase();
    filtered = term ? rows.filter(function(row) {
      return row.some(function(value) {
        return String(value).toLocaleLowerCase().includes(term);
      });
    }) : rows.slice();
    if (sortColumn !== null) {
      filtered.sort(function(a, b) {
        const left = sortable(a[sortColumn]);
        const right = sortable(b[sortColumn]);
        const order = left < right ? -1 : left > right ? 1 : 0;
        return ascending ? order : -order;
      });
    }
    page = 0;
    render();
  });
  previous.addEventListener('click', function() { page -= 1; render(); });
  next.addEventListener('click', function() { page += 1; render(); });
  container.querySelectorAll('th button').forEach(function(button) {
    button.addEventListener('click', function() {
      const column = Number(button.dataset.column);
      ascending = sortColumn === column ? !ascending : true;
      sortColumn = column;
      filtered.sort(function(a, b) {
        const left = sortable(a[column]);
        const right = sortable(b[column]);
        const order = left < right ? -1 : left > right ? 1 : 0;
        return ascending ? order : -order;
      });
      page = 0;
      render();
    });
  });
  render();
});
document.querySelectorAll('.credible-set-browser').forEach(function(browser) {
  const search = browser.querySelector('.set-search');
  const status = browser.querySelector('.set-search-state');
  const cards = Array.from(browser.querySelectorAll('.cs-card')).map(function(card) {
    return {
      element: card,
      text: [card.querySelector('.cs-summary'), card.querySelector('.provenance'),
             card.querySelector('.set-details')].map(function(el) {
        return el.textContent.toLocaleLowerCase();
      }).join(' '),
      variants: Array.from(card.querySelectorAll('.variant-detail')).map(function(variant) {
        return {element: variant, text: variant.textContent.toLocaleLowerCase()};
      })
    };
  });
  const totalVariants = cards.reduce(function(total, card) { return total + card.variants.length; }, 0);
  let savedOpen = null;
  function filter() {
    const term = search.value.trim().toLocaleLowerCase();
    if (term && savedOpen === null) {
      savedOpen = Array.from(browser.querySelectorAll('details')).map(function(el) {
        return {element: el, open: el.open};
      });
    }
    let matchedSets = 0;
    let matchedVariants = 0;
    cards.forEach(function(card) {
      const setMatches = !term || card.text.includes(term);
      let count = 0;
      card.variants.forEach(function(variant) {
        variant.element.hidden = !(setMatches || variant.text.includes(term));
        if (!variant.element.hidden) count += 1;
      });
      card.element.hidden = count === 0;
      if (count) {
        matchedSets += 1;
        matchedVariants += count;
      }
      if (term) card.element.open = count > 0;
    });
    if (!term && savedOpen !== null) {
      savedOpen.forEach(function(saved) { saved.element.open = saved.open; });
      savedOpen = null;
    }
    status.textContent = matchedSets
      ? `${matchedSets} of ${cards.length} sets · ${matchedVariants} of ${totalVariants} variant entries`
      : 'No matching sets or variants.';
  }
  search.addEventListener('input', filter);
  browser.querySelector('.expand-sets').addEventListener('click', function() {
    cards.forEach(function(card) { if (!card.element.hidden) card.element.open = true; });
  });
  browser.querySelector('.collapse-sets').addEventListener('click', function() {
    browser.querySelectorAll('details').forEach(function(el) { el.open = false; });
  });
  browser.querySelector('.table-tools').hidden = false;
  filter();
});

</script>
"""


def _secondary_analysis_section(report: Mapping[str, Any]) -> str:
    """Explain the recorded overlap decision separately from initial results."""
    summary = report["summary"]
    digits = int(report["config"]["probability_significant_digits"])
    groups = summary["overlap_groups"]
    if not groups:
        outcome = (
            "No overlapping loci were found among the completed initial analyses. "
            "No additional fine-mapping was needed. The initial credible sets "
            "were kept unchanged in the final results."
        )
        rerun_metrics = ""
    else:
        outcome = (
            "%s groups of overlapping loci were found. Each group eligible under "
            "the configured limits was analysed again as one larger region, "
            "using the combined boundaries without adding another flank. "
            "The final results use the new credible sets for those groups and "
            "keep the initial sets from non-overlapping loci. Groups that failed, "
            "exceeded a configured limit or produced no credible set contribute "
            "no sets to the final results; their initial sets are not restored."
        ) % _value(groups)
        rerun_metrics = "".join(
            _metric(label, summary[key], tone, digits)
            for label, key, tone in (
                ("Groups producing credible sets", "successful_joint_reruns", "analysis"),
                ("Groups completed without a credible set", "joint_reruns_without_sets", "analysis"),
                ("Groups failed or excluded", "excluded_or_failed_joint_reruns", "warning"),
            )
        )
    final_metrics = "".join((
        _metric("Loci with credible sets in final results", summary["authoritative_loci"], "analysis", digits),
        _metric("Credible sets in final results", summary["final_credible_sets"], "analysis", digits),
    ))
    return (
        '<section id="overlap"><h2>Additional analysis: overlapping loci</h2>'
        '<p class="section-note">After the initial fine-mapping, PostGWAS checked '
        'whether the analysed regions overlap. Overlapping regions can contain '
        'the same association signal. Analysing them together avoids carrying '
        'separate regional results for the same stretch of genome into the final '
        'result files.</p><h3>What happened in this run</h3><p>%s</p>'
        '<div class="metrics">%s%s</div><details><summary>Details of the overlap check</summary>%s'
        '</details></section>'
    ) % (
        escape(outcome), rerun_metrics, final_metrics,
        _interactive_table(
            report["overlap"], table_id="overlap-table",
            page_size=int(report["config"]["page_size"]), significant_digits=digits,
        ),
    )


def _final_credible_set_description(summary: Mapping[str, Any]) -> str:
    """Explain which analysis produced the sets displayed in this run."""
    if not summary["final_credible_sets"]:
        return (
            "No credible sets are available in the final results. Review the initial "
            "locus outcomes and any additional analysis above for the reasons."
        )
    if not summary["overlap_groups"]:
        return (
            "The initial analysis found %s credible %s. No overlapping regions "
            "required another analysis. The results shown below are unchanged "
            "from the initial analysis."
        ) % (
            _value(summary["initial_credible_sets"]),
            "set" if summary["initial_credible_sets"] == 1 else "sets",
        )
    return (
        "The initial analysis examined each region separately. Some regions covered "
        "the same part of the genome, so they were assessed for a second analysis "
        "together. This section combines the unchanged sets from regions that did "
        "not overlap with new sets from the second analysis. For overlapping "
        "regions, the new sets replace the initial sets. Groups without a valid "
        "new set are excluded, as explained above."
    )


def _upstream_report_section(report: Mapping[str, Any]) -> str:
    """Link recorded earlier reports and expose missing files without guessing."""
    rows = report["upstream_reports"]
    if not rows:
        return ""
    tables = []
    for row in rows:
        available = Path(row["path"]).is_file()
        label = {
            "ld_clump": "LD clumping",
            "formatter": "Summary-statistics formatting",
        }.get(row["module"], _title(row["module"]))
        link = _static_table(
            [{"Module": label, "Report": row["path"],
              "Availability": "Available" if available else "Recorded report file not found"}],
            ("Module", "Report", "Availability"),
            path_columns=("Report",) if available else (),
            report_path=Path(report["report_path"]),
            significant_digits=int(report["config"]["probability_significant_digits"]),
        )
        content = read_static_report_content(row["path"]) if available else ""
        if available and not content.strip():
            content = '<p>No static report content was available to embed. Open the original report above.</p>'
        tables.append(
            '<details class="upstream-detail"%s><summary>%s</summary>'
            '<div class="output-table">%s</div><div class="embedded-report">%s</div></details>'
            % (' open' if row["module"] == "ld_clump" else '', escape(label), link, content)
        )
    return (
        '<section id="upstream-reports"><h2>Locus selection and earlier pipeline details</h2>'
        '<p class="section-note">The recorded reports are included below so you '
        'can inspect locus boundaries, lead SNPs, clumping settings, exclusions '
        'and input preparation here. Their table values and explanations are '
        'copied from the original reports, without repeating the analysis. '
        'LD-clumping locus boundaries describe the selection step; the regions '
        'analysed by fine-mapping are shown in the next section.</p>%s</section>'
        % "".join(tables)
    )


def render_fine_mapping_html(report: Mapping[str, Any]) -> str:
    """Render one self-contained, accessible fine-mapping results report."""
    summary = report["summary"]
    config = report["config"]
    digits = int(config["probability_significant_digits"])
    page_size = int(config["page_size"])
    status = report["status"]
    if status == "success":
        needs_review = bool(summary["loci_with_warnings"] or summary["primary_loci_failed"])
        status_label = "Analysis completed — review warnings or exclusions" if needs_review else "Analysis completed"
        status_tone = "warning" if needs_review else "good"
        status_notice = (
            "Start with the initial fine-mapping results below. Warnings and failed "
            "or skipped loci remain part of the interpretation, even when the run "
            "has finished. Any additional analysis is explained separately."
        )
    elif status == "completed_no_credible_sets":
        status_label, status_tone = "complete · no credible sets", "warning"
        status_notice = (
            "The fitting workflow completed, but no credible set met the configured "
            "reporting criteria. No downstream FLAMES handoff is available."
        )
    elif status.startswith("completed"):
        status_label, status_tone = "Analysis completed with excluded regions", "warning"
        status_notice = (
            "Initial results are shown first. Some overlapping regions were excluded "
            "from the final results; the additional-analysis section explains why."
        )
    else:
        status_label, status_tone = "Analysis incomplete", "error"
        status_notice = "The run did not produce a complete final result. Inspect the recorded failures."

    probability_name = (
        "model-wide posterior inclusion probability (PIP)"
        if report["engine"] == "susie"
        else "SNP posterior probability within the selected FINEMAP model"
    )
    if report["engine"] == "susie":
        method_explanation = (
            "SuSiE-RSS models the regional signal as a sum of single effects. Each "
            "reported set is defined from one component's posterior weights (alpha) "
            "at the configured coverage target and must pass the configured LD-purity "
            "filter. The variant value exported here is model-wide PIP, not alpha; "
            "therefore its sum within a set is not the set's component coverage."
        )
    else:
        method_explanation = (
            "FINEMAP evaluates multi-variant causal configurations. PostGWAS first "
            "selects the validated causal-count model by its model posterior "
            "probability (Post-Pr), then reports that model's credible-set members. "
            "The SNP posterior shown here and the causal-count model probability are "
            "different quantities and are kept in separate columns."
        )

    metrics = "".join((
        _metric("Loci analysed", summary["primary_loci_attempted"], "analysis", digits),
        _metric("Loci with credible sets", summary["initial_loci_with_credible_sets"], "analysis", digits),
        _metric("Credible sets found", summary["initial_credible_sets"], "analysis", digits),
        _metric("Completed without a credible set", summary["initial_loci_without_credible_sets"], "analysis", digits),
        _metric("Failed or skipped loci", summary["primary_loci_failed"], "warning", digits),
        _metric("Loci with warnings", summary["loci_with_warnings"], "warning", digits),
    ))
    report_path = Path(report["report_path"])
    output_table = _static_table(
        report["outputs"],
        ("label", "path"),
        path_columns=("path",),
        report_path=report_path,
        significant_digits=digits,
    )
    return """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s · PostGWAS fine-mapping report</title><style>%s</style></head>
<body>
<header class="hero"><div class="hero-inner"><p class="eyebrow">PostGWAS scientific results</p><h1>Fine-mapping report</h1><p class="hero-copy">%s · %s · %s. Posterior evidence ranks variants under the fitted model; it does not by itself establish biological causality.</p><div class="hero-meta"><span class="badge light">%s</span><span class="badge %s">%s</span></div></div></header>
<nav class="report-nav" aria-label="Report sections"><div><a href="#overview">Initial results</a>%s<a href="#loci">Initial locus details</a><a href="#overlap">Additional analysis</a><a href="#credible-sets">Credible sets and variants</a><a href="#interpretation">How to interpret</a><a href="#validation">Input checks</a><a href="#provenance">Methods and files</a></div></nav>
<main>
<div class="notice %s"><strong>%s.</strong> %s</div>
<section id="overview"><h2>Results at a glance</h2><h3>Initial fine-mapping results</h3><p class="section-note">Each locus (a region of the genome) was analysed separately. These counts describe that initial analysis, before checking for overlapping regions. A credible set is a group of candidate variants for an association signal; one locus can produce more than one set.</p><div class="metrics">%s</div><p class="section-note">Completed loci without a credible set are counted separately from failed or skipped loci. Warnings can apply to either outcome and must be reviewed before interpreting the results.</p>%s</section>
%s<section id="loci"><h2>Initial fine-mapping: locus details</h2><p class="section-note">These are the outcomes and quality checks from analysing each locus separately. They include loci with warnings and loci that failed or were skipped, as well as completed analyses without a credible set.</p>%s</section>
%s
<section id="credible-sets"><h2>Credible sets and variant details</h2><p class="section-note">%s</p><p class="section-note">Each bar shows the probability assigned to one variant by the fitted model. A longer bar means stronger support for that variant within this model. These values do not measure the coverage of the whole credible set.</p>%s</section>
<section id="interpretation"><h2>How to interpret this report</h2><div class="concept-grid"><div class="concept"><h3>Posterior inclusion probability</h3><p>The displayed variant statistic is %s. It is conditional on the summary statistics, LD reference, locus boundary, priors, model settings and candidate variants supplied to this run.</p></div><div class="concept"><h3>Credible set</h3><p>A credible set is a model-based collection of variants intended to contain the relevant effect at the configured posterior coverage. Coverage is not a frequentist confidence guarantee, and every member is not necessarily causal.</p></div><div class="concept"><h3>Engine-specific meaning</h3><p>%s</p></div><div class="concept"><h3>Warnings and failed loci</h3><p>A completed run does not erase locus-level warnings. LD/Z mismatch, heterogeneous effective sample size, low purity, recovery, resource limits and failed fits can change coverage or interpretation.</p></div></div></section>
<section id="validation"><h2>Input and resource validation</h2><p class="section-note">These checks were completed before expensive modelling. Genome build and ancestry remain declared properties that must match the GWAS and LD panel; PLINK file contents cannot independently prove population ancestry.</p>%s%s<h3>%s</h3><p class="section-note">This audit records fitting attempts within the initial analysis. It is separate from the later analysis of overlapping regions. For SuSiE, retries use unchanged LD for fitting non-convergence; invalid LD follows the separately recorded validation policy.</p>%s</section>
<section id="provenance"><h2>Resolved model, inputs and reproducibility</h2><p class="section-note">All values below came from the schema-validated effective configuration or observed preflight evidence. This report introduces no independent analysis defaults.</p><h3>Scientific parameters</h3>%s<h3>Input files</h3>%s<h3>Validated runtimes</h3>%s<h3>Output guide</h3><div class="output-table">%s</div></section>
<section><h2>Interpretation guardrails</h2><ul class="callout-list"><li>High PIP prioritises a variant under this model; it does not prove molecular mechanism, target gene, tissue, direction of intervention or clinical relevance.</li><li>Missing or ancestry-mismatched LD, allele errors, summary-statistic/LD inconsistency, variable per-SNP sample size and incomplete candidate-variant coverage can distort posterior results.</li><li>Locus boundaries and the maximum number of effects constrain the candidate model. Compare sensitivity analyses when conclusions depend on these choices.</li><li>Functional annotation and experimental evidence are downstream evidence layers, not replacements for the statistical QC reported here.</li></ul><h3>Method references</h3><ul class="callout-list"><li><a href="https://doi.org/10.1111/rssb.12388">Wang et al. (2020), SuSiE</a> and the <a href="https://stephenslab.github.io/susieR/reference/susie_rss.html">official SuSiE-RSS documentation</a>.</li><li><a href="https://doi.org/10.1093/bioinformatics/btw018">Benner et al. (2016), FINEMAP</a> and the <a href="https://www.christianbenner.com/">official FINEMAP documentation</a>.</li><li><a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8837575/">Schaid et al., statistical fine-mapping review</a>.</li></ul></section>
</main><footer class="footer">PostGWAS fine-mapping report · initial results, additional analysis and final credible sets; consult the linked TSV, model and log files for complete audit evidence.</footer>%s</body></html>""" % (
        escape(report["dataset_id"]),
        _STYLES,
        escape(report["dataset_id"]),
        escape(report["engine_label"]),
        escape(report["genome_build"]),
        escape(report["engine_label"]),
        status_tone,
        escape(status_label),
        '<a href="#upstream-reports">Clumping and input details</a>' if report["upstream_reports"] else "",
        status_tone,
        escape(status_label),
        escape(status_notice),
        metrics,
        _reason_section(
            report["warning_reasons"], report["failure_reasons"],
            report["warning_explanations"], digits,
        ),
        _upstream_report_section(report),
        _interactive_table(
            report["locus_status"], table_id="locus-table",
            page_size=page_size, significant_digits=digits,
        ),
        _secondary_analysis_section(report),
        escape(_final_credible_set_description(summary)),
        _credible_set_cards(report, significant_digits=digits),
        escape(probability_name),
        escape(method_explanation),
        _static_table(
            [
                {"Check": label, "Count": summary[key]}
                for label, key in (
                    ("Input variants", "input_variants"),
                    ("Variants in eligible loci", "variants_in_eligible_loci"),
                    ("Eligible loci", "eligible_loci"),
                    ("Reference ID matches", "reference_variant_matches"),
                    ("Reference coordinate matches", "coordinate_concordant_variant_matches"),
                    ("Locus variants checked against the reference", "locus_variants_requested"),
                    ("Loci with reference variants", "loci_with_reference_variants"),
                )
            ],
            ("Check", "Count"), significant_digits=digits,
        ),
        _interactive_table(
            report["preflight"], table_id="preflight-table",
            page_size=page_size, significant_digits=digits,
        ),
        escape(report["diagnostic_label"]),
        _interactive_table(
            report["diagnostics"], table_id="diagnostic-table",
            page_size=page_size, significant_digits=digits,
        ),
        _static_table(report["parameters"], ("Parameter", "Value"), significant_digits=digits),
        _static_table(
            report["inputs"], ("Input", "Path"), path_columns=("Path",),
            report_path=report_path, significant_digits=digits,
        ),
        _static_table(
            report["tools"], ("Runtime", "Version", "Path"),
            path_columns=("Path",), report_path=report_path,
            significant_digits=digits,
        ),
        output_table,
        _TABLE_SCRIPT,
    )


def write_fine_mapping_html_report(
    args,
    preflight,
    result: Mapping[str, Any],
    output_paths: Mapping[str, Path],
    *,
    upstream_results: Mapping[str, Any] | None = None,
) -> Path:
    """Build, render, and atomically write the final fine-mapping report."""
    report = build_fine_mapping_report(
        args, preflight, result, output_paths, upstream_results=upstream_results,
    )
    destination = output_paths["html_report_file"]
    return write_html_report(render_fine_mapping_html(report), destination)


__all__ = [
    "FineMappingReportError",
    "build_fine_mapping_report",
    "render_fine_mapping_html",
    "write_fine_mapping_html_report",
]
