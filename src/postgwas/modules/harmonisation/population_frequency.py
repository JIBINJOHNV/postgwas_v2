"""Population-frequency compatibility QC for an unfiltered merged GWAS-VCF.

This is a descriptive dataset check, not genetic-ancestry inference. It compares
the study ALT frequency with configured population ALT frequencies already
aligned on each VCF record. The comparison uses the same complete records for
every population, preventing unequal INFO missingness from changing the
variant set behind each correlation. It never filters variants or changes a
declared population.

The population labels follow the reference panel configured by the user. The
packaged fields use the 1000 Genomes Phase 3 super-population frequencies
described by The 1000 Genomes Project Consortium (Nature 2015,
doi:10.1038/nature15393).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
import tempfile
from typing import Any, Mapping

import polars as pl

from postgwas.core.paths import configured_output_path
from postgwas.core.vcf import extract_vcf_table

from .shared.statistics import valid_frequency_mask


class PopulationFrequencyQCError(RuntimeError):
    """The merged-VCF population-frequency check could not be completed."""


def _finite_number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def assess_population_frequency_table(
    table_path: str | Path,
    settings: Mapping[str, Any],
    *,
    delimiter: str,
    null_values: list[str],
) -> dict[str, Any]:
    """Calculate missingness and like-for-like population AF similarity."""
    populations = list(settings["population_fields"])
    columns = ["STUDY_AF", *populations]
    schema = {column: pl.Float64 for column in columns}
    try:
        table = pl.scan_csv(
            table_path,
            separator=delimiter,
            null_values=null_values,
            schema_overrides=schema,
        )
        complete = pl.all_horizontal(
            [valid_frequency_mask(pl.col(column)) for column in columns]
        )
        study_complete = pl.when(complete).then(pl.col("STUDY_AF"))
        expressions: list[pl.Expr] = [
            pl.len().alias("total_records"),
            complete.sum().alias("comparable_variants"),
        ]
        for column in columns:
            present = pl.col(column).is_not_null()
            expressions.extend([
                pl.col(column).is_null().sum().alias("%s_missing" % column),
                (present & ~valid_frequency_mask(pl.col(column))).sum().alias(
                    "%s_invalid" % column
                ),
            ])
        for population in populations:
            population_complete = pl.when(complete).then(pl.col(population))
            expressions.extend([
                pl.corr(
                    study_complete, population_complete, method="pearson"
                ).alias("%s_correlation" % population),
                pl.when(complete)
                .then((pl.col("STUDY_AF") - pl.col(population)).abs())
                .mean()
                .alias("%s_mae" % population),
                pl.when(complete)
                .then(
                    (
                        (1.0 - pl.col("STUDY_AF"))
                        - pl.col(population)
                    ).abs()
                )
                .mean()
                .alias("%s_inverted_mae" % population),
            ])
        values = table.select(expressions).collect().row(0, named=True)
    except Exception as exc:
        raise PopulationFrequencyQCError(
            "Cannot calculate population-frequency statistics from %s: %s"
            % (table_path, exc)
        ) from exc

    total = int(values["total_records"])
    comparable = int(values["comparable_variants"] or 0)
    fields = {}
    for column in columns:
        missing = int(values["%s_missing" % column] or 0)
        invalid = int(values["%s_invalid" % column] or 0)
        fields[column] = {
            "total_records": total,
            "missing": missing,
            "missing_fraction": (missing / total) if total else None,
            "invalid": invalid,
            "invalid_fraction": (invalid / total) if total else None,
        }

    metrics = {}
    for population in populations:
        metrics[population] = {
            "comparable_variants": comparable,
            "pearson_correlation": _finite_number(
                values["%s_correlation" % population]
            ),
            "mean_absolute_difference": _finite_number(
                values["%s_mae" % population]
            ),
            "inverted_mean_absolute_difference": _finite_number(
                values["%s_inverted_mae" % population]
            ),
        }

    minimum_count = int(settings["minimum_comparable_variants"])
    minimum_correlation = float(settings["minimum_correlation"])
    minimum_gap = float(settings["minimum_correlation_gap"])
    inversion_minimum_correlation = float(
        settings["inversion_minimum_absolute_correlation"]
    )
    inversion_maximum_error = float(
        settings["inversion_maximum_mean_absolute_difference"]
    )
    inversion_minimum_improvement = float(
        settings[
            "inversion_minimum_mean_absolute_difference_improvement"
        ]
    )
    candidates = [
        population
        for population, metric in metrics.items()
        if comparable >= minimum_count
        and metric["pearson_correlation"] is not None
        and metric["mean_absolute_difference"] is not None
        and metric["inverted_mean_absolute_difference"] is not None
    ]

    closest = None
    reason = "insufficient comparable population-frequency evidence"
    correlation_gap = None
    inversion_population = None
    negative_candidates = [
        population
        for population in candidates
        if metrics[population]["pearson_correlation"]
        <= -inversion_minimum_correlation
    ]
    if negative_candidates:
        strongest_inversion = min(
            negative_candidates,
            key=lambda population: (
                metrics[population]["inverted_mean_absolute_difference"],
                metrics[population]["pearson_correlation"],
            ),
        )
        inversion_metric = metrics[strongest_inversion]
        if (
            inversion_metric["inverted_mean_absolute_difference"]
            <= inversion_maximum_error
            and inversion_metric["inverted_mean_absolute_difference"]
            + inversion_minimum_improvement
            <= inversion_metric["mean_absolute_difference"]
        ):
            inversion_population = strongest_inversion
            reason = (
                "the study AF appears inverted relative to %s: Pearson correlation "
                "is %.6g, mean absolute difference is %.6g as written and %.6g "
                "after one-minus inversion"
                % (
                    strongest_inversion,
                    inversion_metric["pearson_correlation"],
                    inversion_metric["mean_absolute_difference"],
                    inversion_metric["inverted_mean_absolute_difference"],
                )
            )
    if inversion_population is None and len(candidates) >= 2:
        ranked = sorted(
            candidates,
            key=lambda population: metrics[population]["pearson_correlation"],
            reverse=True,
        )
        correlation_gap = (
            metrics[ranked[0]]["pearson_correlation"]
            - metrics[ranked[1]]["pearson_correlation"]
        )
        if metrics[ranked[0]]["pearson_correlation"] < minimum_correlation:
            reason = (
                "best Pearson correlation %.6g is below %.6g"
                % (
                    metrics[ranked[0]]["pearson_correlation"],
                    minimum_correlation,
                )
            )
        elif correlation_gap <= minimum_gap:
            reason = (
                "best-versus-second Pearson correlation gap %.6g does not exceed %.6g"
                % (correlation_gap, minimum_gap)
            )
        else:
            mae_best = min(
                candidates,
                key=lambda population: metrics[population][
                    "mean_absolute_difference"
                ],
            )
            if bool(settings["require_mae_agreement"]) and mae_best != ranked[0]:
                reason = (
                    "highest correlation (%s) and smallest absolute difference (%s) disagree"
                    % (ranked[0], mae_best)
                )
            else:
                closest = ranked[0]
                reason = (
                    "highest Pearson correlation with supporting absolute-frequency difference"
                )

    status = (
        "frequency_inversion_suspected"
        if inversion_population is not None
        else "complete"
        if closest is not None
        else "inconclusive"
    )
    warnings = (
        [
            "The study frequency appears to belong to the non-effect allele "
            "rather than the effect allele. %s" % reason
        ]
        if inversion_population is not None
        else []
    )
    return {
        "status": status,
        "total_records": total,
        "comparable_variants": comparable,
        "fields": fields,
        "populations": metrics,
        "closest_population": closest,
        "inversion_candidate_population": inversion_population,
        "correlation_gap": correlation_gap,
        "decision_reason": reason,
        "warnings": warnings,
        "decision_thresholds": {
            "minimum_comparable_variants": minimum_count,
            "minimum_correlation": minimum_correlation,
            "minimum_correlation_gap": minimum_gap,
            "inversion_minimum_absolute_correlation": (
                inversion_minimum_correlation
            ),
            "inversion_maximum_mean_absolute_difference": (
                inversion_maximum_error
            ),
            "inversion_minimum_mean_absolute_difference_improvement": (
                inversion_minimum_improvement
            ),
            "require_mae_agreement": bool(settings["require_mae_agreement"]),
        },
    }


def _filename_population_tokens(
    filename: str | Path,
    populations: list[str],
) -> tuple[str, list[str]]:
    """Normalize a basename and find exact configured population tokens."""
    basename = Path(str(filename)).name.upper()
    normalized = "_".join(re.findall(r"[A-Z0-9]+", basename))
    filename_tokens = set(normalized.split("_"))
    return normalized, [
        population
        for population in populations
        if population in filename_tokens
    ]


def check_selected_population(
    closest_population: str | None,
    selected_population: str,
    populations: list[str],
) -> tuple[dict[str, Any], list[str]]:
    """Compare the resolved comparison-AF column with the similarity result."""
    normalized = str(selected_population).strip().upper()
    if normalized not in populations:
        status = "not_comparable"
    elif closest_population is None:
        status = "not_compared"
    elif normalized == closest_population:
        status = "match"
    else:
        status = "mismatch"
    warnings = []
    if status == "mismatch":
        warnings.append(
            "Study frequencies are closest to %s, but the selected comparison-AF "
            "column is %s. Processing will continue."
            % (closest_population, selected_population)
        )
    return {
        "selected_column": str(selected_population),
        "normalized_population": normalized,
        "status": status,
    }, warnings


def check_external_population_filenames(
    closest_population: str | None,
    external_files: Mapping[str, str | Path | None],
    populations: list[str],
) -> tuple[list[dict[str, Any]], list[str]]:
    """Compare only explicit population tokens in sample-sheet basenames."""
    checks = []
    warnings = []
    for label, path in external_files.items():
        if path is None or not str(path).strip():
            continue
        basename = Path(str(path)).name
        normalized, tokens = _filename_population_tokens(path, populations)
        if len(tokens) == 1:
            indicated = tokens[0]
            status = (
                "match"
                if closest_population is not None and indicated == closest_population
                else "mismatch"
                if closest_population is not None
                else "not_compared"
            )
            if status == "mismatch":
                warnings.append(
                    "Study frequencies are closest to %s, but %s filename %s "
                    "indicates %s after filename normalization. Processing will "
                    "continue."
                    % (closest_population, label, basename, indicated)
                )
        elif len(tokens) > 1:
            indicated = None
            status = "ambiguous"
            warnings.append(
                "%s filename %s contains multiple population tokens (%s), so it "
                "cannot be checked against the study-frequency result. Processing "
                "will continue."
                % (label, basename, ", ".join(tokens))
            )
        else:
            indicated = None
            status = "no_population_token"
        checks.append({
            "source": label,
            "filename": basename,
            "normalized_filename": normalized,
            "population_tokens": tokens,
            "indicated_population": indicated,
            "status": status,
        })
    return checks, warnings


def run_population_frequency_qc(
    *,
    vcf_path: str | Path,
    output_directory: str | Path,
    dataset_id: str,
    genome_build: str,
    settings: Mapping[str, Any],
    output_layout: Mapping[str, str],
    table_delimiter: str,
    table_null_values: list[str],
    temporary_table_suffix: str,
    io_buffer_bytes: int,
    bcftools_bin: str,
    selected_population: str,
    external_files: Mapping[str, str | Path | None],
    logger=None,
) -> dict[str, Any]:
    """Extract the configured INFO AF fields and persist the dataset QC."""
    if not bool(settings["enabled"]):
        return {
            "status": "disabled",
            "closest_population": None,
            "warnings": [],
            "external_file_checks": [],
        }

    population_fields = dict(settings["population_fields"])
    columns = {"STUDY_AF": str(settings["study_field"])}
    columns.update(population_fields)
    report_path = configured_output_path(
        output_directory,
        output_layout["population_frequency_qc"],
        error_type=PopulationFrequencyQCError,
        dataset_id=dataset_id,
        build=genome_build,
    )
    temporary_prefix = configured_output_path(
        output_directory,
        output_layout["population_frequency_temporary"],
        error_type=PopulationFrequencyQCError,
        dataset_id=dataset_id,
        build=genome_build,
    )
    temporary_prefix.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix=temporary_prefix.name,
        suffix=temporary_table_suffix,
        dir=temporary_prefix.parent,
        delete=False,
    )
    table_path = Path(handle.name)
    handle.close()
    try:
        extract_vcf_table(
            vcf_path,
            table_path,
            dataset_id,
            columns,
            bcftools_bin,
            delimiter=table_delimiter,
            io_buffer_bytes=io_buffer_bytes,
            allow_undefined_tags=True,
            logger=logger,
            error_type=PopulationFrequencyQCError,
            purpose="Extracting population-frequency fields from the raw merged VCF",
        )
        result = assess_population_frequency_table(
            table_path,
            settings,
            delimiter=table_delimiter,
            null_values=table_null_values,
        )
        selected_check, selected_warnings = check_selected_population(
            result["closest_population"],
            selected_population,
            list(population_fields),
        )
        checks, filename_warnings = check_external_population_filenames(
            result["closest_population"],
            external_files,
            list(population_fields),
        )
        result.update({
            "raw_merged_vcf": str(Path(vcf_path).expanduser().resolve()),
            "vcf_fields": columns,
            "selected_population_check": selected_check,
            "external_file_checks": checks,
            "warnings": list(result.get("warnings") or []) + selected_warnings + (
                filename_warnings
                if bool(settings["warn_on_filename_mismatch"])
                else []
            ),
        })
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        result["report"] = str(report_path)
        if logger is not None:
            logger.record(
                "RESULT",
                "population_frequency_similarity",
                raw_merged_vcf=result["raw_merged_vcf"],
                records=result["total_records"],
                comparable=result["comparable_variants"],
                closest=result["closest_population"] or "inconclusive",
                report=report_path,
            )
        return result
    except OSError as exc:
        raise PopulationFrequencyQCError(
            "Cannot write population-frequency QC report %s: %s"
            % (report_path, exc)
        ) from exc
    finally:
        table_path.unlink(missing_ok=True)


__all__ = [
    "PopulationFrequencyQCError",
    "assess_population_frequency_table",
    "check_external_population_filenames",
    "check_selected_population",
    "run_population_frequency_qc",
]
