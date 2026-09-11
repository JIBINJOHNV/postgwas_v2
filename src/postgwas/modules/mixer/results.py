"""Compact reports for official single-trait MiXeR and GSA-MiXeR output."""

from __future__ import annotations

import csv
import heapq
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence, Type

from postgwas.core.io.delimiters import open_text
from postgwas.core.io.reports import write_delimited_report, write_yaml_report
from postgwas.core.reference_resources import validate_table_header
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.values import (
    format_count,
    format_fraction_percentage,
    format_number,
)


# These expressions describe MiXeR v2's upstream log protocol. They are not
# scientific thresholds or user policy; unavailable values remain null.
_LOG_PATTERNS = {
    "version": re.compile(r"MiXeR v([^\s]+) software"),
    "coordinate_fallback_matches": re.compile(r"(\d+) lines matched via CHR:BP:A1:A2"),
    "reference_unmatched": re.compile(r"(\d+) lines were ignored .* did not match reference"),
    "strand_ambiguous": re.compile(r"(\d+) variants .* strand-ambiguous"),
    "alleles_flipped": re.compile(r"(\d+) variants had flipped A1/A2 alleles"),
    "well_defined_variants": re.compile(r"Found (\d+) variants with well-defined Z and N"),
}


def inspect_mixer_input(
    input_file: str | Path,
    *,
    delimiter: str,
    required_columns: Sequence[str],
    chromosome_column: str,
    configured_chromosomes: Sequence[str],
    error_type: Type[Exception] = RuntimeError,
) -> dict[str, Any]:
    """Stream the formatter artifact once and verify its exact MiXeR schema."""
    counts: dict[str, int] = {}
    row_count = 0
    try:
        with open_text(input_file) as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            header = next(reader, None)
            if not header:
                raise error_type("MiXeR input has no header: %s" % input_file)
            if len(header) != len(set(header)):
                raise error_type("MiXeR input contains duplicate column names: %s" % input_file)
            missing = [name for name in required_columns if name not in header]
            if missing:
                raise error_type(
                    "MiXeR input is missing configured formatter columns: %s"
                    % ", ".join(missing)
                )
            chromosome_index = header.index(chromosome_column)
            for line_number, row in enumerate(reader, 2):
                if not row or not any(value.strip() for value in row):
                    continue
                if len(row) != len(header):
                    raise error_type(
                        "MiXeR input row %d has %d fields; the header has %d"
                        % (line_number, len(row), len(header))
                    )
                chromosome = row[chromosome_index].strip()
                counts[chromosome] = counts.get(chromosome, 0) + 1
                row_count += 1
    except error_type:
        raise
    except (OSError, UnicodeError, ValueError, csv.Error) as exc:
        raise error_type("Could not validate MiXeR input %s: %s" % (input_file, exc)) from exc
    if row_count == 0:
        raise error_type("MiXeR input contains no variant rows: %s" % input_file)

    configured = list(configured_chromosomes)
    observed = set(counts)
    return {
        "formatted_variants": row_count,
        "columns": header,
        "chromosome_counts": {
            chromosome: counts[chromosome]
            for chromosome in configured if chromosome in counts
        },
        "chromosomes_present": [value for value in configured if value in observed],
        "chromosomes_missing": [value for value in configured if value not in observed],
        "chromosomes_unexpected": sorted(observed.difference(configured)),
    }


def _load_result(
    path: str | Path,
    label: str,
    analysis: str,
    error_type: Type[Exception],
) -> dict[str, Any]:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise error_type("Cannot read %s MiXeR JSON %s: %s" % (label, path, exc)) from exc
    if not isinstance(value, dict) or value.get("analysis") != analysis:
        raise error_type(
            "%s is not an official MiXeR %s result: %s" % (label, analysis, path)
        )
    return value


def _estimate(
    result: Mapping[str, Any],
    metric: str,
    statistic: str,
    label: str,
    error_type: Type[Exception],
) -> float:
    values = (result.get("ci") or {}).get(metric) or {}
    value = values.get(statistic)
    if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise error_type(
            "%s MiXeR result has no finite ci.%s.%s value"
            % (label, metric, statistic)
        )
    return float(value)


def _final_optimisation(result: Mapping[str, Any], field: str) -> dict[str, Any]:
    stages = result.get(field) or []
    if not stages or not isinstance(stages[-1], list) or len(stages[-1]) != 2:
        return {}
    values = stages[-1][1]
    return values if isinstance(values, dict) else {}


def _model_assessment(
    fit: Mapping[str, Any], error_type: Type[Exception],
) -> dict[str, Any]:
    mixture = _final_optimisation(fit, "optimize")
    infinitesimal = _final_optimisation(fit, "inft_optimize")
    if mixture.get("success") is not True or infinitesimal.get("success") is not True:
        raise error_type("MiXeR mixture or infinitesimal optimisation did not converge")
    mixture_aic, mixture_bic = mixture.get("AIC"), mixture.get("BIC")
    infinitesimal_aic = infinitesimal.get("AIC")
    infinitesimal_bic = infinitesimal.get("BIC")
    values = (mixture_aic, mixture_bic, infinitesimal_aic, infinitesimal_bic)
    if not all(
        isinstance(value, (int, float)) and math.isfinite(float(value))
        for value in values
    ):
        raise error_type("MiXeR fit result does not contain finite fit-model AIC/BIC values")
    return {
        "converged": True,
        "mixture_aic": float(mixture_aic),
        "mixture_bic": float(mixture_bic),
        "infinitesimal_aic": float(infinitesimal_aic),
        "infinitesimal_bic": float(infinitesimal_bic),
        "mixture_vs_infinitesimal_delta_aic": float(infinitesimal_aic) - float(mixture_aic),
        "mixture_vs_infinitesimal_delta_bic": float(infinitesimal_bic) - float(mixture_bic),
    }


def _log_metrics(path: str | Path) -> dict[str, Any]:
    values: dict[str, Any] = {}
    source = Path(path)
    if not source.is_file():
        return values
    with source.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            for name, pattern in _LOG_PATTERNS.items():
                match = pattern.search(line)
                if match:
                    values[name] = (
                        match.group(1) if name == "version" else int(match.group(1))
                    )
    return values


def _nonfinite_count(value: Any) -> int:
    if isinstance(value, float):
        return int(not math.isfinite(value))
    if isinstance(value, dict):
        return sum(_nonfinite_count(item) for item in value.values())
    if isinstance(value, list):
        return sum(_nonfinite_count(item) for item in value)
    return 0


def _has_uncertainty(result: Mapping[str, Any]) -> bool:
    uncertainty_names = {"std", "lower", "upper", "low", "high", "l", "u"}
    return any(
        uncertainty_names.intersection(values)
        for values in (result.get("ci") or {}).values()
        if isinstance(values, dict)
    )


def _quality_issue(
    warnings: list[str], action: str, message: str, error_type: Type[Exception],
) -> None:
    if action == "error":
        raise error_type(message)
    if action == "warning":
        warnings.append(message)


def build_mixer_summary(
    *,
    dataset_id: str,
    run_id: str,
    result: Mapping[str, Any],
    input_metrics: Mapping[str, Any],
    configuration,
    input_warnings: Sequence[str] = (),
    formatter_metrics: Mapping[str, Any] | None = None,
    error_type: Type[Exception] = RuntimeError,
) -> dict[str, Any]:
    """Extract compact estimates and quality evidence from fit1/test1 outputs."""
    reporting = configuration.modules.mixer.reporting
    statistic = reporting.figure_statistics[0]
    fit = _load_result(result["fit"], "fit", "univariate", error_type)
    test = _load_result(result["test"], "test", "univariate", error_type)
    fit_options = fit.get("options") or {}
    test_options = test.get("options") or {}
    fit_seed, test_seed = fit_options.get("seed"), test_options.get("seed")
    log_metrics = _log_metrics(result["fit_log"])
    qq_nonfinite = _nonfinite_count(test.get("qqplot") or {})
    model = _model_assessment(fit, error_type)

    warnings = list(input_warnings)
    quality = reporting.quality
    if fit_seed != test_seed:
        _quality_issue(
            warnings, quality.seed_mismatch_action,
            "MiXeR fit1 and test1 did not use the same configured seed", error_type,
        )
    if qq_nonfinite:
        _quality_issue(
            warnings, quality.nonfinite_qq_action,
            "MiXeR QQ diagnostics contain %d nonfinite values" % qq_nonfinite,
            error_type,
        )
    required_log_metrics = ("reference_unmatched", "strand_ambiguous")
    unavailable_log_metrics = [
        name for name in required_log_metrics if log_metrics.get(name) is None
    ]
    if unavailable_log_metrics:
        _quality_issue(
            warnings, quality.missing_log_metrics_action,
            "MiXeR fit log did not expose reference-matching metrics: %s"
            % ", ".join(unavailable_log_metrics), error_type,
        )
    if not _has_uncertainty(fit):
        _quality_issue(
            warnings, quality.missing_uncertainty_action,
            "MiXeR estimates do not include confidence intervals", error_type,
        )
    if model["mixture_vs_infinitesimal_delta_aic"] <= 0:
        _quality_issue(
            warnings, quality.insufficient_model_power_action,
            "Fit AIC does not support the finite-mixture model over the infinitesimal model",
            error_type,
        )
    elif model["mixture_vs_infinitesimal_delta_bic"] <= 0:
        _quality_issue(
            warnings, quality.insufficient_model_power_action,
            "Fit AIC supports the finite-mixture model but conservative BIC does not",
            error_type,
        )

    architecture = {
        "estimate_statistic": statistic,
        "polygenicity_pi": _estimate(fit, "pi", statistic, "fit", error_type),
        "estimated_causal_variants": _estimate(fit, "nc", statistic, "fit", error_type),
        "variants_explaining_90_percent_h2": _estimate(
            fit, "nc@p9", statistic, "fit", error_type,
        ),
        "discoverability_sig2_beta": _estimate(
            fit, "sig2_beta", statistic, "fit", error_type,
        ),
        "native_snp_heritability": _estimate(fit, "h2", statistic, "fit", error_type),
        "residual_inflation_fit": _estimate(
            fit, "sig2_zero", statistic, "fit", error_type,
        ),
        "residual_inflation_test": _estimate(
            test, "sig2_zero", statistic, "test", error_type,
        ),
        "uncertainty": fit.get("ci"),
    }
    variants_analysed = int(
        fit_options.get("num_tag") or log_metrics.get("well_defined_variants") or 0
    )
    formatted = int(input_metrics["formatted_variants"])
    if variants_analysed <= 0 or variants_analysed > formatted:
        raise error_type(
            "MiXeR reports %d analysed variants from %d formatted variants"
            % (variants_analysed, formatted)
        )
    if not unavailable_log_metrics:
        accounted = (
            variants_analysed
            + int(log_metrics["reference_unmatched"])
            + int(log_metrics["strand_ambiguous"])
        )
        if accounted != formatted:
            raise error_type(
                "MiXeR input accounting failed: %d analysed + unmatched + "
                "strand-ambiguous variants does not equal %d formatted variants"
                % (accounted, formatted)
            )
    input_qc = {
        **dict(input_metrics),
        "formatter_input_variants": (
            formatter_metrics.get("rows_in") if formatter_metrics else None
        ),
        "formatter_excluded_variants": (
            formatter_metrics.get("rows_excluded") if formatter_metrics else None
        ),
        "coordinate_fallback_matches": log_metrics.get("coordinate_fallback_matches"),
        "reference_unmatched": log_metrics.get("reference_unmatched"),
        "strand_ambiguous": log_metrics.get("strand_ambiguous"),
        "alleles_flipped": log_metrics.get("alleles_flipped"),
        "variants_analysed": variants_analysed,
        "analysed_fraction": variants_analysed / formatted if formatted else None,
    }
    return {
        "run": {
            "status": "completed",
            "analysis": "univariate",
            "dataset_id": dataset_id,
            "run_id": run_id,
            "execution_backend": result["execution_backend"],
            "genome_build": result["genome_build"],
            "mixer_version": log_metrics.get("version"),
            "container_image": (
                configuration.resources.containers.mixer.image
                if result["execution_backend"] == "docker" else None
            ),
            "fit_seed": fit_seed,
            "test_seed": test_seed,
            "threads": fit_options.get("threads", configuration.execution.threads),
        },
        "input_qc": input_qc,
        "architecture": architecture,
        "model_fit": {**model, "model_selection_source": "fit"},
        "diagnostics": {
            "qq_nonfinite_values": qq_nonfinite,
            "qq_variants": (test.get("qqplot") or {}).get("n_snps"),
        },
        "quality": {
            "status": "warning" if warnings else "pass",
            "warnings": warnings,
        },
        "artifacts": {
            "mixer_input": result["mixer_input"],
            "fit_json": result["fit"],
            "test_json": result["test"],
            "fit_log": result["fit_log"],
            "test_log": result["test_log"],
            "canonical_log": result["log_file"],
        },
    }


def _summary_records(summary: Mapping[str, Any]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []

    def add(section, metric, value, unit, status, source, description):
        if isinstance(value, list):
            value = ",".join(str(item) for item in value)
        records.append({
            "section": section, "metric": metric, "value": value, "unit": unit,
            "status": status, "source": source, "description": description,
        })

    qc = summary["input_qc"]
    for metric, unit, source, description in (
        ("formatted_variants", "variants", "formatter input", "Variants supplied to MiXeR"),
        ("reference_unmatched", "variants", "fit log", "Variants not matched to the reference"),
        ("strand_ambiguous", "variants", "fit log", "Strand-ambiguous variants excluded by MiXeR"),
        ("coordinate_fallback_matches", "variants", "fit log", "Variants matched by chromosome, position and alleles"),
        ("alleles_flipped", "variants", "fit log", "Variants whose A1/A2 orientation and Z sign were flipped"),
        ("variants_analysed", "variants", "fit JSON", "Reference-compatible tag variants analysed"),
        ("analysed_fraction", "proportion", "fit JSON and input", "Fraction of formatted variants analysed"),
    ):
        add("input_qc", metric, qc.get(metric), unit, "reported", source, description)
    add("input_qc", "chromosomes_missing", qc["chromosomes_missing"], "chromosomes", "warning" if qc["chromosomes_missing"] else "pass", "formatter input", "Configured chromosomes absent from input")
    add("input_qc", "chromosomes_unexpected", qc["chromosomes_unexpected"], "chromosomes", "warning" if qc["chromosomes_unexpected"] else "pass", "formatter input", "Input chromosomes outside the configured range")

    architecture = summary["architecture"]
    for metric, unit, description in (
        ("polygenicity_pi", "proportion", "Estimated non-null proportion of reference variants"),
        ("estimated_causal_variants", "variants", "Estimated number of non-null variants"),
        ("variants_explaining_90_percent_h2", "variants", "Estimated variants explaining 90% of model SNP heritability"),
        ("discoverability_sig2_beta", "effect-size variance", "Variance of effects among non-null variants"),
        ("native_snp_heritability", "proportion", "MiXeR native-scale SNP heritability"),
        ("residual_inflation_fit", "factor", "Residual inflation estimated during fit1"),
        ("residual_inflation_test", "factor", "Residual inflation estimated during test1"),
    ):
        add("architecture", metric, architecture[metric], unit, "reported", "fit JSON" if not metric.endswith("_test") else "test JSON", description)

    model = summary["model_fit"]
    add("model_fit", "converged", model["converged"], "boolean", "pass", "fit JSON", "Mixture and infinitesimal optimisations converged")
    add("model_fit", "mixture_vs_infinitesimal_delta_aic", model["mixture_vs_infinitesimal_delta_aic"], "AIC units", "pass" if model["mixture_vs_infinitesimal_delta_aic"] > 0 else "warning", "fit JSON", "Positive values favour the finite-mixture model")
    add("model_fit", "mixture_vs_infinitesimal_delta_bic", model["mixture_vs_infinitesimal_delta_bic"], "BIC units", "pass" if model["mixture_vs_infinitesimal_delta_bic"] > 0 else "warning", "fit JSON", "Positive values favour the finite-mixture model")
    add("diagnostics", "qq_nonfinite_values", summary["diagnostics"]["qq_nonfinite_values"], "values", "warning" if summary["diagnostics"]["qq_nonfinite_values"] else "pass", "test JSON", "Nonfinite values in QQ diagnostics")
    for number, warning in enumerate(summary["quality"]["warnings"], 1):
        add("quality", "warning_%d" % number, warning, "message", "warning", "PostGWAS assessment", "Scientifically relevant limitation")
    return records


def write_mixer_summary(
    summary: Mapping[str, Any],
    yaml_path: str | Path,
    tsv_path: str | Path,
    reporting,
) -> tuple[str, str]:
    written_tsv = write_delimited_report(
        _summary_records(summary),
        tsv_path,
        fieldnames=[
            "section", "metric", "value", "unit", "status", "source", "description",
        ],
        delimiter=reporting.table_delimiter,
        null_value=reporting.null_value,
    )
    written_yaml = write_yaml_report(summary, yaml_path)
    return str(written_yaml), str(written_tsv)


def _float_or_none(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def validate_gsa_go_file(path: str | Path, settings, error_type=RuntimeError) -> None:
    """Check ontology header and first-record presence through the shared reader."""
    validate_table_header(
        path, delimiter=settings.go_file_delimiter,
        required_columns=settings.go_file_required_columns,
        label="GSA-MiXeR ontology file", error_type=error_type,
    )


def _gsa_test_identifiers(path: str | Path, settings, error_type) -> set[str]:
    try:
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter=settings.go_file_delimiter)
            header = reader.fieldnames or []
            missing = [name for name in settings.go_file_required_columns if name not in header]
            if missing:
                raise error_type(
                    "GSA-MiXeR test ontology file is missing configured columns: %s"
                    % ", ".join(missing)
                )
            identifiers = {
                str(row.get(settings.go_identifier_column) or "").strip()
                for row in reader
            }
    except error_type:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise error_type("Cannot read GSA-MiXeR test ontology file %s: %s" % (path, exc)) from exc
    identifiers.discard("")
    if not identifiers:
        raise error_type("GSA-MiXeR test ontology file contains no identifiers: %s" % path)
    return identifiers


def _gsa_enrichment_metrics(
    result_path: str | Path,
    test_go_file: str | Path,
    settings,
    error_type,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    test_identifiers = _gsa_test_identifiers(test_go_file, settings, error_type)
    columns = settings.result_columns.model_dump()
    required = list(columns.values())
    total_rows = 0
    evidence_rows = 0
    missing_aic = 0
    seen: set[str] = set()
    top: list[tuple[float, int, dict[str, Any]]] = []
    try:
        with open_text(result_path) as handle:
            reader = csv.DictReader(handle, delimiter=settings.result_delimiter)
            header = reader.fieldnames or []
            missing = [name for name in required if name not in header]
            if missing:
                raise error_type(
                    "GSA-MiXeR enrichment result is missing configured columns: %s"
                    % ", ".join(missing)
                )
            for row_number, row in enumerate(reader, 1):
                total_rows += 1
                identifier = str(row.get(columns["identifier"]) or "").strip()
                if identifier not in test_identifiers:
                    continue
                seen.add(identifier)
                record = {
                    semantic: (
                        identifier if semantic == "identifier"
                        else _float_or_none(row.get(source))
                    )
                    for semantic, source in columns.items()
                }
                model_aic = record["model_aic"]
                if model_aic is None:
                    missing_aic += 1
                    continue
                if model_aic <= settings.evidence_aic_threshold:
                    continue
                evidence_rows += 1
                ranking = record["enrichment"]
                if ranking is None:
                    continue
                candidate = (ranking, row_number, record)
                if len(top) < settings.top_results:
                    heapq.heappush(top, candidate)
                elif candidate[:2] > top[0][:2]:
                    heapq.heapreplace(top, candidate)
    except error_type:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise error_type("Cannot read GSA-MiXeR enrichment result %s: %s" % (result_path, exc)) from exc
    missing_identifiers = sorted(test_identifiers.difference(seen))
    if missing_identifiers:
        raise error_type(
            "GSA-MiXeR enrichment output is missing %d configured test identifiers; "
            "examples: %s"
            % (len(missing_identifiers), ", ".join(missing_identifiers[:5]))
        )
    records = [item[2] for item in sorted(top, reverse=True)]
    return ({
        "result_rows": total_rows,
        "tested_identifiers": len(test_identifiers),
        "identifiers_with_positive_aic_evidence": evidence_rows,
        "identifiers_missing_aic": missing_aic,
        "evidence_aic_threshold": settings.evidence_aic_threshold,
        "ranking_metric": "enrichment",
        "top_results_written": len(records),
    }, records)


def build_gsa_summary(
    *,
    dataset_id: str,
    run_id: str,
    result: Mapping[str, Any],
    input_metrics: Mapping[str, Any],
    configuration,
    error_type: Type[Exception] = RuntimeError,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Validate GSA base/full models and summarize configured test gene sets."""
    settings = configuration.modules.mixer.gsa
    baseline = _load_result(result["baseline_json"], "GSA baseline", "plsa", error_type)
    full = _load_result(result["full_json"], "GSA full", "plsa", error_type)
    baseline_options = baseline.get("options") or {}
    full_options = full.get("options") or {}
    if baseline_options.get("gsa_base") is not True:
        raise error_type("GSA baseline JSON was not produced with --gsa-base")
    if full_options.get("gsa_full") is not True:
        raise error_type("GSA full JSON was not produced with --gsa-full")
    if str(baseline_options.get("seed")) != str(full_options.get("seed")):
        raise error_type("GSA baseline and full models did not use the same configured seed")
    metrics, top_records = _gsa_enrichment_metrics(
        result["enrichment_results"], result["test_go_file"], settings, error_type,
    )
    log_metrics = _log_metrics(result["full_log"])
    summary = {
        "run": {
            "status": "completed",
            "analysis": "single_trait_gsa_mixer",
            "dataset_id": dataset_id,
            "run_id": run_id,
            "execution_backend": result["execution_backend"],
            "genome_build": result["genome_build"],
            "mixer_version": log_metrics.get("version"),
            "seed": full_options.get("seed"),
            "threads": full_options.get("threads", configuration.execution.threads),
        },
        "input_qc": dict(input_metrics),
        "models": {
            "baseline": {
                "variants": baseline_options.get("num_snp"),
                "tag_variants": baseline_options.get("num_tag"),
                "sample_size": baseline_options.get("trait1_nval"),
            },
            "full": {
                "variants": full_options.get("num_snp"),
                "tag_variants": full_options.get("num_tag"),
                "sample_size": full_options.get("trait1_nval"),
            },
        },
        "gene_set_assessment": metrics,
        "interpretation": {
            "evidence_rule": (
                "loglike_aic greater than the configured threshold supports the "
                "gene-set-specific model over its null comparison"
            ),
            "multiple_testing_p_values_computed": False,
        },
        "artifacts": {
            "mixer_input": result["mixer_input"],
            "split_sumstats_pattern": result["split_sumstats_pattern"],
            "baseline_json": result["baseline_json"],
            "baseline_log": result["baseline_log"],
            "full_json": result["full_json"],
            "full_log": result["full_log"],
            "enrichment_results": result["enrichment_results"],
            "canonical_log": result["log_file"],
        },
    }
    return summary, top_records


def write_gsa_summary(
    summary: Mapping[str, Any],
    top_records: Sequence[Mapping[str, Any]],
    yaml_path: str | Path,
    tsv_path: str | Path,
    configuration,
) -> tuple[str, str]:
    reporting = configuration.modules.mixer.reporting
    fieldnames = list(configuration.modules.mixer.gsa.result_columns.model_dump())
    written_tsv = write_delimited_report(
        top_records,
        tsv_path,
        fieldnames=fieldnames,
        delimiter=reporting.table_delimiter,
        null_value=reporting.null_value,
    )
    written_yaml = write_yaml_report(summary, yaml_path)
    return str(written_yaml), str(written_tsv)


def render_mixer_summary(summary: Mapping[str, Any]) -> str:
    run = summary["run"]
    qc = summary["input_qc"]
    architecture = summary["architecture"]
    model = summary["model_fit"]
    artifacts = summary["artifacts"]
    lines = ["", screen_line("analysis", "MiXeR single-trait analysis", indent=2)]
    lines.extend([
        screen_field("info", "Dataset", run["dataset_id"], indent=6, label_width=30),
        screen_field("success", "Analysis status", run["status"], indent=6, label_width=30),
        screen_field("info", "Execution backend", run["execution_backend"], indent=6, label_width=30),
        "",
        screen_line("genetic", "Input and reference matching", indent=6),
        screen_field("count", "Formatted variants", format_count(qc["formatted_variants"]), indent=10, label_width=30),
        screen_field("loss", "Unmatched to reference", format_count(qc.get("reference_unmatched"), missing="not available"), indent=10, label_width=30),
        screen_field("loss", "Strand-ambiguous variants", format_count(qc.get("strand_ambiguous"), missing="not available"), indent=10, label_width=30),
        screen_field("count", "Variants analysed", "%s  (%s)" % (format_count(qc["variants_analysed"]), format_fraction_percentage(qc.get("analysed_fraction"))), indent=10, label_width=30),
        "",
        screen_line("decision", "Estimated genetic architecture", indent=6),
        screen_field("info", "Polygenicity", format_number(architecture["polygenicity_pi"], "%.4f"), indent=10, label_width=30),
        screen_field("info", "Estimated causal variants", format_count(architecture["estimated_causal_variants"]), indent=10, label_width=30),
        screen_field("info", "Variants explaining 90% h²", format_count(architecture["variants_explaining_90_percent_h2"]), indent=10, label_width=30),
        screen_field("info", "Discoverability", format_number(architecture["discoverability_sig2_beta"], "%.4e"), indent=10, label_width=30),
        screen_field("info", "Native SNP heritability", format_number(architecture["native_snp_heritability"], "%.4f"), indent=10, label_width=30),
        screen_field("info", "Residual inflation", format_number(architecture["residual_inflation_test"], "%.4f"), indent=10, label_width=30),
        "",
        screen_line("analysis", "Model assessment", indent=6),
        screen_field("success", "Optimisation", "converged", indent=10, label_width=30),
        screen_field("info", "Mixture vs infinitesimal ΔAIC", format_number(model["mixture_vs_infinitesimal_delta_aic"], "%.4f"), indent=10, label_width=30),
        screen_field("info", "Mixture vs infinitesimal ΔBIC", format_number(model["mixture_vs_infinitesimal_delta_bic"], "%.4f"), indent=10, label_width=30),
    ])
    warnings = summary["quality"]["warnings"]
    if warnings:
        lines.extend(["", screen_line("warning", "Quality warnings", indent=6)])
        lines.extend(screen_line("loss", warning, indent=10) for warning in warnings)
    lines.extend([
        "",
        screen_field("info", "Compact YAML summary", artifacts["summary_yaml"], indent=6, label_width=30),
        screen_field("info", "Analysis-ready TSV", artifacts["summary_tsv"], indent=6, label_width=30),
        screen_field("info", "Diagnostic files", _figure_location(artifacts.get("figures")), indent=6, label_width=30),
        screen_field("info", "Full canonical log", artifacts["canonical_log"], indent=6, label_width=30),
        "",
    ])
    return "\n".join(lines)


def render_gsa_summary(summary: Mapping[str, Any]) -> str:
    run = summary["run"]
    models = summary["models"]
    assessment = summary["gene_set_assessment"]
    artifacts = summary["artifacts"]
    return "\n".join([
        "",
        screen_line("analysis", "GSA-MiXeR single-trait gene-set analysis", indent=2),
        screen_field("info", "Dataset", run["dataset_id"], indent=6, label_width=32),
        screen_field("success", "Analysis status", run["status"], indent=6, label_width=32),
        "",
        screen_line("genetic", "Model and gene-set assessment", indent=6),
        screen_field("count", "Baseline tag variants", format_count(models["baseline"]["tag_variants"]), indent=10, label_width=32),
        screen_field("count", "Full-model tag variants", format_count(models["full"]["tag_variants"]), indent=10, label_width=32),
        screen_field("count", "Tested gene-set identifiers", format_count(assessment["tested_identifiers"]), indent=10, label_width=32),
        screen_field("success", "Positive AIC evidence", format_count(assessment["identifiers_with_positive_aic_evidence"]), indent=10, label_width=32),
        screen_field("info", "Evidence threshold", "loglike_aic > %s" % format_number(assessment["evidence_aic_threshold"]), indent=10, label_width=32),
        "",
        screen_field("info", "Complete enrichment table", artifacts["enrichment_results"], indent=6, label_width=32),
        screen_field("info", "Compact YAML summary", artifacts["summary_yaml"], indent=6, label_width=32),
        screen_field("info", "Ranked evidence TSV", artifacts["top_results_tsv"], indent=6, label_width=32),
        screen_field("info", "Full canonical log", artifacts["canonical_log"], indent=6, label_width=32),
        "",
    ])


def _figure_location(values: Sequence[str] | None) -> str:
    if not values:
        return "not generated"
    return "%d files in %s" % (len(values), Path(values[0]).parent)


__all__ = [
    "build_gsa_summary",
    "build_mixer_summary",
    "inspect_mixer_input",
    "render_gsa_summary",
    "render_mixer_summary",
    "validate_gsa_go_file",
    "write_gsa_summary",
    "write_mixer_summary",
]
