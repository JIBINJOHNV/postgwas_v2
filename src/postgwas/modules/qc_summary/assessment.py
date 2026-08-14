"""Fast QC assessment of a merged GWAS-VCF without creating a filtered VCF.

The VCF is extracted once with ``bcftools query``. Polars evaluates every
configured QC rule independently against the raw records, combines all active
failure masks, and calculates the final virtual QC-passed metrics in one query.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
from typing import Any

import polars as pl

from postgwas.config.models.modules.qc_summary import QCSummaryConfig
from postgwas.core.vcf import VCF_TAG, extract_vcf_table
from postgwas.core.paths import configured_output_path


class VcfAssessmentError(RuntimeError):
    """The raw VCF could not be converted to a trustworthy QC assessment."""


_TABLE_COLUMNS = (
    "CHROM", "POS", "REF", "ALT",
    "INFO_AF", "INFO_EXTERNAL_AF", "FORMAT_AF", "FORMAT_SI", "FORMAT_LP",
    "FORMAT_NEF",
)


def extract_vcf_assessment_table(
    vcf_path: str | Path,
    table_path: str | Path,
    dataset_id: str,
    external_af_name: str,
    vcf_fields: dict[str, str],
    table_delimiter: str,
    io_buffer_bytes: int,
    bcftools_bin: str,
    logger=None,
) -> str:
    """Extract the ten fields needed by all QC calculations in one VCF pass."""
    if not VCF_TAG.fullmatch(str(external_af_name)):
        raise VcfAssessmentError(
            "The comparison allele-frequency column must be a valid VCF INFO tag "
            "beginning with a letter: %r" % external_af_name
        )
    vcf = Path(vcf_path).expanduser().resolve()
    table = Path(table_path)
    if not vcf.is_file() or vcf.stat().st_size <= 0:
        raise VcfAssessmentError("Merged VCF does not exist or is empty: %s" % vcf)
    columns = _resolved_vcf_fields(vcf_fields, external_af_name)
    if logger is not None:
        logger.record(
            "INPUT", "vcf_qc_extraction",
            vcf=str(vcf), external_af="INFO/%s" % external_af_name,
        )
    return extract_vcf_table(
        vcf, table, dataset_id, columns, bcftools_bin,
        delimiter=table_delimiter,
        io_buffer_bytes=io_buffer_bytes,
        allow_undefined_tags=True,
        logger=logger,
        error_type=VcfAssessmentError,
        purpose="Extracting VCF fields for QC assessment",
    )


def _resolved_vcf_fields(
    configured: dict[str, str],
    external_af_name: str,
) -> dict[str, str]:
    """Map configured bcftools queries onto the stable internal QC schema."""
    keys = (
        "chromosome", "position", "reference_allele", "alternate_allele",
        "study_info_af", "external_info_af", "study_format_af",
        "imputation_format", "log_pvalue_format", "effective_sample_size_format",
    )
    missing = [key for key in keys if not str(configured.get(key, "")).strip()]
    if missing:
        raise VcfAssessmentError(
            "QC VCF field mapping is missing: %s" % ", ".join(missing)
        )
    values = [
        str(configured[key]).format(external_af=external_af_name)
        for key in keys
    ]
    return dict(zip(_TABLE_COLUMNS, values))


def _field_label(query: str) -> str:
    """Render a bcftools query expression as a readable VCF field name."""
    query = query.strip()
    if query.startswith("[") and query.endswith("]"):
        return "FORMAT/" + query[1:-1].lstrip("%")
    return query.lstrip("%")


def _qc_field_labels(vcf_fields: dict[str, str], external_af_name: str) -> dict[str, str]:
    resolved = _resolved_vcf_fields(vcf_fields, external_af_name)
    return {
        "study_info_af": _field_label(resolved["INFO_AF"]),
        "external_info_af": _field_label(resolved["INFO_EXTERNAL_AF"]),
        "study_format_af": _field_label(resolved["FORMAT_AF"]),
        "imputation_format": _field_label(resolved["FORMAT_SI"]),
        "log_pvalue_format": _field_label(resolved["FORMAT_LP"]),
        "effective_sample_size_format": _field_label(resolved["FORMAT_NEF"]),
    }


def _text(column: str) -> pl.Expr:
    return pl.col(column).cast(pl.String, strict=False)


def _number(column: str) -> pl.Expr:
    return pl.col(column).cast(pl.Float64, strict=False)


def _missing(column: str) -> pl.Expr:
    return _text(column).is_null()


def _present(column: str) -> pl.Expr:
    return ~_missing(column)


def _usable_sample_size() -> pl.Expr:
    sample_size = _number("FORMAT_NEF")
    return _present("FORMAT_NEF") & sample_size.is_finite() & (sample_size > 0.0)


def _safe(mask: pl.Expr) -> pl.Expr:
    return mask.fill_null(False)


def _missing_outcome(policy_value: str) -> str:
    return "fail QC" if policy_value == "remove" else "pass QC"


def _count(mask: pl.Expr, name: str) -> pl.Expr:
    return _safe(mask).cast(pl.UInt64).sum().alias(name)


def _all_rows() -> pl.Expr:
    """A true expression with the table's row cardinality, including null CHROM."""
    chromosome = pl.col("CHROM")
    return chromosome.is_null() | chromosome.is_not_null()


def _variant_masks() -> dict[str, pl.Expr]:
    ref = _text("REF").str.to_uppercase()
    alt = _text("ALT").str.to_uppercase()
    bases = ["A", "C", "G", "T"]
    snp = (
        (ref.str.len_chars() == 1)
        & (alt.str.len_chars() == 1)
        & ref.is_in(bases)
        & alt.is_in(bases)
    )
    transition = snp & (
        ((ref == "A") & (alt == "G"))
        | ((ref == "G") & (alt == "A"))
        | ((ref == "C") & (alt == "T"))
        | ((ref == "T") & (alt == "C"))
    )
    palindromic = snp & (
        ((ref == "A") & (alt == "T"))
        | ((ref == "T") & (alt == "A"))
        | ((ref == "C") & (alt == "G"))
        | ((ref == "G") & (alt == "C"))
    )
    return {
        "snp": _safe(snp),
        "transition": _safe(transition),
        "transversion": _safe(snp & ~transition),
        "palindromic": _safe(palindromic),
    }


def build_qc_rules(
    configuration: QCSummaryConfig,
    genome_build: str,
    external_af_name: str,
    field_labels: dict[str, str],
) -> list[dict[str, Any]]:
    """Build independent rules from the resolved QC configuration."""
    if not VCF_TAG.fullmatch(str(external_af_name)):
        raise VcfAssessmentError(
            "Invalid comparison allele-frequency INFO tag: %r" % external_af_name
        )
    masks = _variant_masks()
    af = _number("FORMAT_AF")
    study_info_af = _number("INFO_AF")
    external_af = _number("INFO_EXTERNAL_AF")
    info = _number("FORMAT_SI")
    lp = _number("FORMAT_LP")
    rules: list[dict[str, Any]] = []
    format_af_label = field_labels["study_format_af"]
    format_info_label = field_labels["imputation_format"]
    format_lp_label = field_labels["log_pvalue_format"]
    study_info_af_label = field_labels["study_info_af"]
    external_info_af_label = field_labels["external_info_af"]

    rules_config = configuration.rules
    lp_cutoff = rules_config.minimum_neglog10_p
    if lp_cutoff is not None:
        lp_missing = rules_config.missing_pvalue_action
        missing_lp = _missing("FORMAT_LP")
        below = _present("FORMAT_LP") & (lp < float(lp_cutoff))
        rules.append({
            "key": "significance",
            "label": "Association significance",
            "criterion": "%s ≥ %.6g; missing values %s"
            % (format_lp_label, float(lp_cutoff), _missing_outcome(lp_missing)),
            "failure": below | (missing_lp if lp_missing == "remove" else pl.lit(False)),
            "details": [
                {"label": "Missing %s" % format_lp_label, "mask": missing_lp,
                 "removes": lp_missing == "remove"},
                {"label": "Below the LP threshold", "mask": below, "removes": True},
            ],
        })

    maf = float(rules_config.maf_min)
    af_missing_action = rules_config.missing_af_action
    missing_format_af = _missing("FORMAT_AF")
    outside_maf = _present("FORMAT_AF") & ((af < maf) | (af > 1.0 - maf))
    rules.append({
        "key": "allele_frequency",
        "label": "Study allele frequency",
        "criterion": "%.6g ≤ %s ≤ %.6g; missing values %s"
        % (maf, format_af_label, 1.0 - maf, _missing_outcome(af_missing_action)),
        "failure": outside_maf | (
            missing_format_af if af_missing_action == "remove" else pl.lit(False)
        ),
        "details": [
            {"label": "Missing %s" % format_af_label, "mask": missing_format_af,
             "removes": af_missing_action == "remove"},
            {"label": "Outside the allele-frequency range", "mask": outside_maf,
             "removes": True},
        ],
    })

    info_min = float(rules_config.info_min)
    info_max = float(rules_config.info_max)
    info_missing_action = rules_config.missing_info_action
    missing_info = _missing("FORMAT_SI")
    outside_info = _present("FORMAT_SI") & (
        (info < info_min) | (info > info_max)
    )
    rules.append({
        "key": "imputation_quality",
        "label": "Imputation quality",
        "criterion": "%.6g ≤ %s ≤ %.6g; missing values %s"
        % (info_min, format_info_label, info_max, _missing_outcome(info_missing_action)),
        "failure": outside_info | (
            missing_info if info_missing_action == "remove" else pl.lit(False)
        ),
        "details": [
            {"label": "Missing %s" % format_info_label, "mask": missing_info,
             "removes": info_missing_action == "remove"},
            {"label": "Outside the imputation-quality range", "mask": outside_info,
             "removes": True},
        ],
    })

    difference_cutoff = float(rules_config.maximum_af_difference)
    missing_study_info_af = _missing("INFO_AF")
    missing_external_af = _missing("INFO_EXTERNAL_AF")
    discordant = (
        _present("INFO_AF")
        & _present("INFO_EXTERNAL_AF")
        & ((study_info_af - external_af).abs() > difference_cutoff)
    )
    missing_for_comparison = missing_study_info_af | missing_external_af
    rules.append({
        "key": "external_af_concordance",
        "label": "Reference-frequency concordance",
        "criterion": "|%s − %s| ≤ %.6g; missing values %s"
        % (
            study_info_af_label, external_info_af_label, difference_cutoff,
            _missing_outcome(af_missing_action),
        ),
        "failure": discordant | (
            missing_for_comparison
            if af_missing_action == "remove" else pl.lit(False)
        ),
        "details": [
            {"label": "Missing %s" % study_info_af_label, "mask": missing_study_info_af,
             "removes": af_missing_action == "remove"},
            {"label": "Missing %s" % external_info_af_label,
             "mask": missing_external_af, "removes": af_missing_action == "remove"},
            {"label": "Absolute frequency difference above %.6g" % difference_cutoff,
             "mask": discordant, "removes": True},
        ],
    })

    if not rules_config.include_indels:
        non_snp = ~masks["snp"]
        rules.append({
            "key": "variant_type",
            "label": "Variant type",
            "criterion": "single-nucleotide variants only",
            "failure": non_snp,
            "details": [
                {"label": "Indels and other non-SNP variants", "mask": non_snp,
                 "removes": True},
            ],
        })

    if rules_config.remove_palindromic:
        lower = float(rules_config.palindromic_lower)
        upper = float(rules_config.palindromic_upper)
        ambiguous = (
            masks["palindromic"]
            & _present("FORMAT_AF")
            & (af >= lower)
            & (af <= upper)
        )
        rules.append({
            "key": "palindromic_variants",
            "label": "Strand-ambiguous variants",
            "criterion": "palindromic SNP with %.6g ≤ %s ≤ %.6g"
            % (lower, format_af_label, upper),
            "failure": ambiguous,
            "details": [
                {"label": "Frequency-ambiguous palindromic SNPs", "mask": ambiguous,
                 "removes": True},
            ],
        })

    if rules_config.remove_mhc:
        mhc_region = configuration.mhc_region(genome_build)
        configured_chromosome = str(mhc_region.chromosome)
        chromosome = configured_chromosome.lower()
        if chromosome.startswith("chr"):
            chromosome = chromosome[3:]
        start = int(mhc_region.start)
        end = int(mhc_region.end)
        normalized_chromosome = _text("CHROM").str.to_lowercase().str.replace(
            r"^chr", ""
        )
        in_mhc = (
            (normalized_chromosome == chromosome)
            & (pl.col("POS").cast(pl.Int64, strict=False) >= start)
            & (pl.col("POS").cast(pl.Int64, strict=False) <= end)
        )
        rules.append({
            "key": "mhc_region",
            "label": "MHC region",
            "criterion": "%s:%s-%s inclusive" % (
                configured_chromosome, format(start, ","), format(end, ",")
            ),
            "failure": in_mhc,
            "details": [
                {"label": "Variants inside the configured MHC region", "mask": in_mhc,
                 "removes": True},
            ],
        })

    return rules


def _metric_expressions(
    prefix: str,
    subset: pl.Expr,
    external_difference_cutoff: float,
    sample_size_sd_multiplier: float,
) -> list[pl.Expr]:
    masks = _variant_masks()
    study_af_missing = _missing("INFO_AF")
    external_af_missing = _missing("INFO_EXTERNAL_AF")
    comparable = ~study_af_missing & ~external_af_missing
    discordant = comparable & (
        (_number("INFO_AF") - _number("INFO_EXTERNAL_AF")).abs()
        > external_difference_cutoff
    )
    sample_size = _number("FORMAT_NEF")
    usable_sample_size = _usable_sample_size()
    stage_sample_size = _safe(subset) & usable_sample_size
    mean_column = "__%s_sample_size_mean" % prefix
    sd_column = "__%s_sample_size_sd" % prefix
    outlier_threshold = (
        pl.col(mean_column) + sample_size_sd_multiplier * pl.col(sd_column)
    )
    return [
        _count(subset, prefix + "__num_records"),
        _count(subset & masks["snp"], prefix + "__num_snps"),
        _count(subset & ~masks["snp"], prefix + "__num_non_snps"),
        _count(subset & masks["transition"], prefix + "__transitions"),
        _count(subset & masks["transversion"], prefix + "__transversions"),
        _count(subset & _missing("FORMAT_AF"), prefix + "__format_af_missing"),
        _count(subset & _missing("FORMAT_SI"), prefix + "__format_si_missing"),
        _count(subset & study_af_missing, prefix + "__study_af_missing"),
        _count(subset & external_af_missing, prefix + "__external_af_missing"),
        _count(subset & comparable, prefix + "__af_comparable"),
        _count(subset & discordant, prefix + "__af_difference_above_cutoff"),
        _count(stage_sample_size, prefix + "__effective_sample_size_available"),
        _count(
            _safe(subset) & ~usable_sample_size,
            prefix + "__effective_sample_size_missing_or_invalid",
        ),
        pl.when(stage_sample_size).then(sample_size).min().alias(
            prefix + "__effective_sample_size_minimum"
        ),
        pl.when(stage_sample_size).then(sample_size).max().alias(
            prefix + "__effective_sample_size_maximum"
        ),
        pl.col(mean_column).first().alias(prefix + "__effective_sample_size_mean"),
        pl.col(sd_column).first().alias(
            prefix + "__effective_sample_size_standard_deviation"
        ),
        outlier_threshold.first().alias(
            prefix + "__effective_sample_size_outlier_threshold"
        ),
        _count(
            stage_sample_size & (sample_size > outlier_threshold),
            prefix + "__effective_sample_size_above_outlier_threshold",
        ),
    ]


def _with_sample_size_statistics(
    frame: pl.LazyFrame,
    stage_subsets: dict[str, pl.Expr],
) -> pl.LazyFrame:
    """Add stage-specific Neff mean and sample SD without another table read."""
    sample_size = _number("FORMAT_NEF")
    usable = _usable_sample_size()
    expressions = []
    for prefix, subset in stage_subsets.items():
        values = pl.when(_safe(subset) & usable).then(sample_size)
        expressions.extend([
            values.mean().over(pl.lit(1)).alias(
                "__%s_sample_size_mean" % prefix
            ),
            values.std(ddof=1).over(pl.lit(1)).alias(
                "__%s_sample_size_sd" % prefix
            ),
        ])
    return frame.with_columns(expressions)


def assess_variant_table(
    table_path: str | Path,
    configuration: QCSummaryConfig,
    genome_build: str,
    external_af_name: str,
    vcf_fields: dict[str, str],
    *,
    delimiter: str,
    null_values: list[str],
) -> dict[str, Any]:
    """Calculate raw-rule and combined final QC metrics in one table scan."""
    table = Path(table_path)
    schema = dict((column, pl.String) for column in _TABLE_COLUMNS)
    frame = pl.scan_csv(
        table,
        separator=delimiter,
        null_values=null_values,
        schema_overrides=schema,
        infer_schema_length=0,
    ).with_columns(pl.col("POS").cast(pl.Int64, strict=False))
    field_labels = _qc_field_labels(vcf_fields, external_af_name)
    rules = build_qc_rules(
        configuration, genome_build, external_af_name, field_labels,
    )
    difference_cutoff = float(configuration.rules.maximum_af_difference)
    sample_size_sd_multiplier = float(
        configuration.rules.sample_size_outlier_standard_deviations
    )
    rule_selections = []

    combined_failure = pl.lit(False)
    failure_count = pl.lit(0, dtype=pl.UInt16)
    for index, rule in enumerate(rules, 1):
        failure = _safe(rule["failure"])
        rule_selections.extend([
            _count(failure, "rule_%02d__failed_raw" % index),
        ])
        for detail_index, detail in enumerate(rule["details"], 1):
            detail_mask = _safe(detail["mask"])
            rule_selections.extend([
                _count(
                    detail_mask,
                    "rule_%02d__detail_%02d__matched_raw" % (index, detail_index),
                ),
            ])
        combined_failure = combined_failure | failure
        failure_count = failure_count + failure.cast(pl.UInt16)

    combined_failure_column = "__combined_qc_failure"
    frame = frame.with_columns(
        _safe(combined_failure).alias(combined_failure_column)
    )
    qc_passed_subset = ~pl.col(combined_failure_column)
    frame = _with_sample_size_statistics(
        frame,
        {"raw": _all_rows(), "qc_passed": qc_passed_subset},
    )
    selections = _metric_expressions(
        "raw", _all_rows(), difference_cutoff, sample_size_sd_multiplier,
    )
    selections.extend(rule_selections)
    selections.extend([
        _count(failure_count > 1, "assessment__overlap_variants"),
        (failure_count.cast(pl.Int64) - 1)
        .clip(lower_bound=0)
        .sum()
        .alias("assessment__extra_rule_matches"),
    ])
    selections.extend(_metric_expressions(
        "qc_passed", qc_passed_subset, difference_cutoff,
        sample_size_sd_multiplier,
    ))
    values = frame.select(selections).collect(engine="streaming").to_dicts()[0]

    def stage_metrics(prefix: str) -> dict[str, Any]:
        integer_metrics = {
            "num_records", "num_snps", "num_non_snps", "transitions",
            "transversions", "format_af_missing", "format_si_missing",
            "study_af_missing", "external_af_missing", "af_comparable",
            "af_difference_above_cutoff", "effective_sample_size_available",
            "effective_sample_size_missing_or_invalid",
            "effective_sample_size_above_outlier_threshold",
        }
        metrics = {}
        for key, value in values.items():
            if not key.startswith(prefix + "__"):
                continue
            metric = key.split("__", 1)[1]
            metrics[metric] = (
                int(value or 0) if metric in integer_metrics
                else (None if value is None else float(value))
            )
        transversions = metrics.get("transversions", 0)
        metrics["ts_tv_ratio"] = (
            metrics.get("transitions", 0) / transversions
            if transversions else None
        )
        return metrics

    rule_results = []
    for index, rule in enumerate(rules, 1):
        prefix = "rule_%02d__" % index
        details = []
        for detail_index, detail in enumerate(rule["details"], 1):
            detail_prefix = prefix + "detail_%02d__" % detail_index
            details.append({
                "label": detail["label"],
                "action": "fail_qc" if detail.get("removes") else "pass_qc",
                "matched_raw": int(values.get(detail_prefix + "matched_raw") or 0),
            })
        rule_results.append({
            "number": index,
            "key": rule["key"],
            "label": rule["label"],
            "criterion": rule["criterion"],
            "failed_raw": int(values.get(prefix + "failed_raw") or 0),
            "details": details,
        })

    raw = stage_metrics("raw")
    qc_passed = stage_metrics("qc_passed")
    rule_match_total = sum(rule["failed_raw"] for rule in rule_results)
    excluded_total = raw["num_records"] - qc_passed["num_records"]
    extra_rule_matches = int(values.get("assessment__extra_rule_matches") or 0)
    if raw["num_records"] != excluded_total + qc_passed["num_records"]:
        raise VcfAssessmentError(
            "QC accounting failed: %s raw records do not equal %s excluded plus "
            "%s retained. No assessment report was accepted."
            % (raw["num_records"], excluded_total, qc_passed["num_records"])
        )
    if rule_match_total != excluded_total + extra_rule_matches:
        raise VcfAssessmentError(
            "QC rule accounting failed: %s raw-rule matches do not equal %s unique "
            "failed variants plus %s overlapping matches."
            % (rule_match_total, excluded_total, extra_rule_matches)
        )
    return {
        "definition": (
            "all active configured QC conditions are applied together to every raw "
            "merged-VCF record; a variant is QC-passed only when it passes every "
            "condition. The raw VCF is unchanged and no filtered VCF is created"
        ),
        "raw": raw,
        "rules": rule_results,
        "active_rule_count": len(rule_results),
        "qc_passed": qc_passed,
        "excluded_total": excluded_total,
        "rule_match_total": rule_match_total,
        "accounting_balanced": True,
        "retained_fraction": (
            qc_passed["num_records"] / raw["num_records"]
            if raw["num_records"] else 0.0
        ),
        "overlap_variants": int(values.get("assessment__overlap_variants") or 0),
        "extra_rule_matches": extra_rule_matches,
        "external_af_name": external_af_name,
        "field_labels": field_labels,
        "af_difference_cutoff": difference_cutoff,
        "sample_size_outlier_standard_deviations": sample_size_sd_multiplier,
    }


def _write_assessment_reports(
    assessment: dict[str, Any],
    report_paths: dict[str, Path],
    *,
    delimiter: str,
    null_output: str,
) -> dict[str, str]:
    for path in report_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    staged_paths = {}
    for name, destination in report_paths.items():
        handle = tempfile.NamedTemporaryFile(
            prefix=".%s." % destination.name,
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        )
        handle.close()
        staged_paths[name] = Path(handle.name)

    try:
        with staged_paths["summary"].open(
            "w", encoding="utf-8", newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["stage", "metric", "value"],
                delimiter=delimiter,
            )
            writer.writeheader()
            for stage in ("raw", "qc_passed"):
                for metric, value in assessment[stage].items():
                    writer.writerow({
                        "stage": stage,
                        "metric": metric,
                        "value": null_output if value is None else value,
                    })
            for metric in (
                "active_rule_count", "excluded_total", "rule_match_total",
                "accounting_balanced", "retained_fraction", "overlap_variants",
                "extra_rule_matches", "sample_size_outlier_standard_deviations",
                "definition",
            ):
                value = assessment[metric]
                writer.writerow({
                    "stage": "assessment", "metric": metric,
                    "value": null_output if value is None else value,
                })

        with staged_paths["rules"].open(
            "w", encoding="utf-8", newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "record_type", "rule_number", "rule", "criterion", "detail",
                    "action", "variants_matching_in_raw_vcf",
                ],
                delimiter=delimiter,
            )
            writer.writeheader()
            for rule in assessment["rules"]:
                writer.writerow({
                    "record_type": "rule",
                    "rule_number": rule["number"],
                    "rule": rule["label"],
                    "criterion": rule["criterion"],
                    "action": "fail_qc",
                    "variants_matching_in_raw_vcf": rule["failed_raw"],
                })
                for detail in rule["details"]:
                    writer.writerow({
                        "record_type": "detail",
                        "rule_number": rule["number"],
                        "rule": rule["label"],
                        "criterion": rule["criterion"],
                        "detail": detail["label"],
                        "action": detail["action"],
                        "variants_matching_in_raw_vcf": detail["matched_raw"],
                    })

        staged_paths["json"].write_text(
            json.dumps(assessment, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        for name, destination in report_paths.items():
            staged_paths[name].replace(destination)
    finally:
        for path in staged_paths.values():
            path.unlink(missing_ok=True)
    return {
        "summary": str(report_paths["summary"]),
        "rules": str(report_paths["rules"]),
        "json": str(report_paths["json"]),
    }


def run_vcf_qc_assessment(
    *,
    vcf_path: str | Path,
    output_directory: str | Path,
    dataset_id: str,
    genome_build: str,
    external_af_name: str,
    configuration: QCSummaryConfig,
    bcftools_bin: str,
    logger=None,
) -> dict[str, Any]:
    """Extract once, assess raw rules, persist summaries, and remove the TSV."""
    vcf_fields = configuration.vcf_fields.model_dump()
    output_layout = configuration.output_layout
    table_configuration = configuration.table
    report_paths = {
        "summary": configured_output_path(
            output_directory, output_layout.metric_report,
            error_type=VcfAssessmentError, dataset_id=dataset_id, build=genome_build,
        ),
        "rules": configured_output_path(
            output_directory, output_layout.rule_report,
            error_type=VcfAssessmentError, dataset_id=dataset_id, build=genome_build,
        ),
        "json": configured_output_path(
            output_directory, output_layout.assessment_json,
            error_type=VcfAssessmentError, dataset_id=dataset_id, build=genome_build,
        ),
    }
    temporary_prefix = configured_output_path(
        output_directory, output_layout.temporary_table_prefix,
        error_type=VcfAssessmentError, dataset_id=dataset_id, build=genome_build,
    )
    temporary_prefix.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix=temporary_prefix.name,
        suffix=table_configuration.temporary_suffix,
        dir=temporary_prefix.parent,
        delete=False,
    )
    table_path = Path(handle.name)
    handle.close()
    try:
        extract_vcf_assessment_table(
            vcf_path=vcf_path,
            table_path=table_path,
            dataset_id=dataset_id,
            external_af_name=external_af_name,
            vcf_fields=vcf_fields,
            table_delimiter=table_configuration.delimiter,
            io_buffer_bytes=table_configuration.io_buffer_bytes,
            bcftools_bin=bcftools_bin,
            logger=logger,
        )
        assessment = assess_variant_table(
            table_path, configuration, genome_build, external_af_name, vcf_fields,
            delimiter=table_configuration.delimiter,
            null_values=table_configuration.null_values,
        )
        assessment["raw_vcf"] = str(Path(vcf_path).expanduser().resolve())
        assessment["reports"] = _write_assessment_reports(
            assessment,
            report_paths,
            delimiter=table_configuration.delimiter,
            null_output=table_configuration.null_output,
        )
        if logger is not None:
            logger.record(
                "RESULT", "vcf_qc_assessment",
                raw=assessment["raw"]["num_records"],
                qc_passed=assessment["qc_passed"]["num_records"],
                excluded=assessment["excluded_total"],
                sample_size_outlier_standard_deviations=assessment[
                    "sample_size_outlier_standard_deviations"
                ],
                raw_sample_size_outliers=assessment["raw"][
                    "effective_sample_size_above_outlier_threshold"
                ],
                qc_passed_sample_size_outliers=assessment["qc_passed"][
                    "effective_sample_size_above_outlier_threshold"
                ],
                vcf_output="raw_merged_only",
            )
        return assessment
    finally:
        try:
            table_path.unlink()
        except OSError:
            pass
