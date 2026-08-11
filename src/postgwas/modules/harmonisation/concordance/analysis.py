"""Pure Polars analysis for input-versus-GWAS-VCF concordance.

The comparison is intentionally independent of the GWAS-to-VCF adapter.  It
uses the sample-sheet mapping, transforms every statistic to the VCF ALT-allele
orientation, and reports both retention and per-value agreement.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np
import polars as pl
from scipy.stats import norm

from ..coordinates import position_text_expression


INPUT_SOURCE_ROW_COLUMN = "__concordance_input_source_row"
VCF_SOURCE_ROW_COLUMN = "__concordance_vcf_source_row"
PARTITION_COLUMN = "__concordance_chromosome_partition"

VARIANT_TYPES = (
    ("snps", "snp"),
    ("indels", "indel"),
    ("other_variants", "other"),
)
MATCH_TYPE_FIELDS = {
    "direct": "direct_matches",
    "allele_swapped": "swapped_matches",
    "strand_complement": "complement_matches",
    "strand_complement_swapped": "complement_swapped_matches",
    "palindromic_as_listed": "palindromic_as_listed_matches",
    "palindromic_swapped": "palindromic_swapped_matches",
    "palindromic_forward": "palindromic_forward_matches",
    "palindromic_forward_swapped": "palindromic_forward_swapped_matches",
    "palindromic_reverse_complement": "palindromic_reverse_complement_matches",
    "palindromic_reverse_complement_swapped": (
        "palindromic_reverse_complement_swapped_matches"
    ),
}
POSITION_COUNT_FIELDS = (
    "input_unmatched_positions",
    "vcf_unmatched_positions",
    "unmatched_position_union",
    "shared_unmatched_positions",
    "unambiguous_position_pairs",
    "ambiguous_shared_positions",
    "same_position_input_variants",
    "same_position_vcf_variants",
    "input_specific_variants",
    "vcf_specific_variants",
    "input_specific_positions",
    "vcf_specific_positions",
    "position_value_checked_pairs",
    "position_value_concordant_pairs",
    "position_value_mismatch_pairs",
)
METRIC_OBSERVED_COLUMNS = {
    "effect": "vcf_effect",
    "standard_error": "vcf_se",
    "allele_frequency": "observed_frequency",
    "z_score": "vcf_z",
}


@dataclass
class ConcordanceAnalysis:
    status: str
    summary: dict[str, Any]
    metric_summary: dict[str, dict[str, Any]]
    metric_summary_by_variant_type: dict[str, dict[str, dict[str, Any]]]
    position_metric_summary_by_variant_type: dict[str, dict[str, dict[str, Any]]]
    matched: pl.DataFrame | pl.LazyFrame
    position_matches: pl.DataFrame | pl.LazyFrame
    mismatches: pl.DataFrame | pl.LazyFrame
    input_only: pl.DataFrame | pl.LazyFrame
    vcf_only: pl.DataFrame | pl.LazyFrame
    vcf_duplicates: pl.DataFrame | pl.LazyFrame


def _chromosome(expr: pl.Expr) -> pl.Expr:
    value = (
        expr.cast(pl.String, strict=False)
        .str.strip_chars()
        .str.to_uppercase()
        .str.replace(r"^CHR", "")
    )
    return value.replace({"23": "X", "24": "Y", "25": "XY", "26": "MT", "M": "MT"})


def _allele(expr: pl.Expr) -> pl.Expr:
    return expr.cast(pl.String, strict=False).str.strip_chars().str.to_uppercase()


def _complement(expr: pl.Expr) -> pl.Expr:
    return (
        pl.when(expr == "A").then(pl.lit("T"))
        .when(expr == "T").then(pl.lit("A"))
        .when(expr == "C").then(pl.lit("G"))
        .when(expr == "G").then(pl.lit("C"))
        .otherwise(None)
    )


def _canonical_key(chrom: str, pos: str, first: str, second: str) -> pl.Expr:
    left = pl.col(first)
    right = pl.col(second)
    ordered = pl.when(left <= right).then(
        pl.concat_str([left, right], separator=":")
    ).otherwise(pl.concat_str([right, left], separator=":"))
    return pl.concat_str([
        pl.col(chrom).cast(pl.String),
        pl.col(pos).cast(pl.String),
        ordered,
    ], separator=":")


def _numeric(frame: pl.DataFrame, name: str | None) -> pl.Expr:
    if not name or name not in frame.columns:
        return pl.lit(None, dtype=pl.Float64)
    return pl.col(name).cast(pl.Float64, strict=False)


def _input_coordinates(row, policies) -> tuple[pl.Expr, pl.Expr]:
    if row.chromosome_column and row.position_column:
        return pl.col(row.chromosome_column), pl.col(row.position_column)
    combined = pl.col(row.chromosome_position_column).cast(pl.String, strict=False)
    pattern = r"(?i)^(?:chr)?([^:_\-]+)[:_\-]+(.+)$"
    return (
        combined.str.extract(pattern, 1),
        position_text_expression(
            combined.str.extract(pattern, 2),
            str(policies.get("position.extraction")),
        ),
    )


def input_chromosome_expression(row, policies) -> pl.Expr:
    """Return the canonical chromosome expression shared by staging and analysis."""
    chromosome, _ = _input_coordinates(row, policies)
    return _chromosome(chromosome)


def reference_chromosome_expression(column: str) -> pl.Expr:
    """Return the canonical chromosome expression for a configured reference column."""
    return _chromosome(pl.col(column))


def vcf_chromosome_expression() -> pl.Expr:
    """Return the canonical chromosome expression for the extracted VCF table."""
    return _chromosome(pl.col("CHROM"))


def _p_to_raw(values: np.ndarray, p_value_type: str) -> np.ndarray:
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        if p_value_type == "raw":
            return values.astype(float, copy=True)
        if p_value_type == "neglog10":
            return np.power(10.0, -values, dtype=float)
        if p_value_type == "negln":
            return np.exp(-values, dtype=float)
    return np.full(values.shape, np.nan, dtype=float)


def prepare_input_table(
    frame: pl.DataFrame,
    row,
    *,
    effect_type: str,
    se_scale: str | None = None,
    p_value_type: str,
    eaf_is_maf: bool | None,
    settings,
    policies,
    external_eaf_frame: pl.DataFrame | None = None,
    external_eaf_mapping=None,
) -> pl.DataFrame:
    """Normalize only representation; indels are never reference-normalized."""
    chrom, pos = _input_coordinates(row, policies)
    effect = _numeric(frame, row.effect_column)
    beta = (
        pl.when(effect > 0).then(effect.log()).otherwise(None)
        if effect_type == "odds_ratio" else effect
    )
    supplied_z = _numeric(frame, row.z_score_column)
    supplied_standard_error = _numeric(frame, row.standard_error_column)
    if (
        effect_type == "odds_ratio"
        and row.standard_error_column
        and se_scale not in ("log_odds", "as_given")
    ):
        raise ValueError(
            "Odds-ratio concordance requires the study-level SE-scale decision "
            "('log_odds' or 'as_given')."
        )
    standard_error = (
        pl.when(effect > 0)
        .then(supplied_standard_error / effect)
        .otherwise(None)
        if effect_type == "odds_ratio" and se_scale == "as_given"
        else supplied_standard_error
    )
    input_af = _numeric(frame, row.effect_allele_frequency_column)
    p_value = _numeric(frame, row.p_value_column)
    duplicate_action = (
        [pl.col("duplicate_action").cast(pl.String, strict=False)]
        if "duplicate_action" in frame.columns else []
    )
    duplicate_input_row = (
        [pl.col("duplicate_input_row").cast(pl.UInt32, strict=False)]
        if "duplicate_input_row" in frame.columns else []
    )

    indexed = (
        frame
        if INPUT_SOURCE_ROW_COLUMN in frame.columns
        else frame.with_row_index(INPUT_SOURCE_ROW_COLUMN, offset=1)
    )
    result = (
        indexed
        .select(
            pl.col(INPUT_SOURCE_ROW_COLUMN).cast(pl.UInt32).alias("input_row"),
            _chromosome(chrom).alias("input_chrom"),
            pos.cast(pl.Int64, strict=False).alias("input_pos"),
            _allele(pl.col(row.effect_allele_column)).alias("input_effect_allele"),
            _allele(pl.col(row.other_allele_column)).alias("input_other_allele"),
            beta.alias("input_beta"),
            standard_error.alias("input_se"),
            supplied_z.alias("input_z_supplied"),
            input_af.alias("input_af"),
            p_value.alias("input_p"),
            *duplicate_action,
            *duplicate_input_row,
        )
        .with_columns(
            (
                pl.col("input_z_supplied").is_null()
                & pl.col("input_beta").is_not_null()
                & pl.col("input_se").is_not_null()
                & (pl.col("input_se") > 0)
            ).alias("_can_z_from_se")
        )
        .with_columns(
            pl.when(pl.col("input_z_supplied").is_not_null())
            .then(pl.col("input_z_supplied"))
            .when(pl.col("_can_z_from_se"))
            .then(pl.col("input_beta") / pl.col("input_se"))
            .otherwise(None)
            .alias("input_z"),
            pl.when(pl.col("input_z_supplied").is_not_null())
            .then(pl.lit("supplied"))
            .when(pl.col("_can_z_from_se"))
            .then(pl.lit("calculated_beta_div_se"))
            .otherwise(pl.lit("unavailable"))
            .alias("input_z_source"),
        )
        .drop("_can_z_from_se")
    )

    # If neither supplied Z nor BETA/SE can provide Z, use the configured
    # P-value scale and tail. scipy's survival function remains stable at
    # genome-wide p-values where `1 - p` loses floating-point precision.
    needs_p = result.get_column("input_z").is_null().to_numpy()
    beta_values = result.get_column("input_beta").to_numpy()
    p_values = result.get_column("input_p").to_numpy()
    raw_p = _p_to_raw(p_values, p_value_type)
    usable_p = needs_p & np.isfinite(beta_values) & np.isfinite(raw_p)
    clipped = np.clip(
        raw_p,
        float(policies.get("pvalue.clip_low")),
        float(policies.get("pvalue.clip_high")),
    )
    divisor = float(policies.get("pvalue.se_tail"))
    calculated = np.sign(beta_values) * norm.isf(clipped / divisor)
    calculated[~usable_p] = np.nan
    result = result.with_columns(pl.Series("_z_from_p", calculated, dtype=pl.Float64))
    result = result.with_columns(
        pl.when(pl.col("input_z").is_null() & pl.col("_z_from_p").is_not_null())
        .then(pl.col("_z_from_p"))
        .otherwise(pl.col("input_z"))
        .alias("input_z"),
        pl.when(pl.col("input_z").is_null() & pl.col("_z_from_p").is_not_null())
        .then(pl.lit("calculated_effect_and_p"))
        .otherwise(pl.col("input_z_source"))
        .alias("input_z_source"),
    ).drop("_z_from_p")
    result = result.with_columns(
        pl.when(pl.col("input_se").is_not_null())
        .then(pl.col("input_se"))
        .when(
            pl.col("input_beta").is_not_null()
            & pl.col("input_z").is_not_null()
            & (pl.col("input_z") != 0)
        )
        .then((pl.col("input_beta") / pl.col("input_z")).abs())
        .otherwise(None)
        .alias("input_se"),
        pl.when(pl.col("input_se").is_not_null())
        .then(pl.lit("supplied_canonical_scale"))
        .when(
            pl.col("input_beta").is_not_null()
            & pl.col("input_z").is_not_null()
            & (pl.col("input_z") != 0)
        )
        .then(pl.lit("calculated_abs_beta_div_z"))
        .otherwise(pl.lit("unavailable"))
        .alias("input_se_source"),
    )

    result = result.with_columns(
        (
            pl.col("input_chrom").is_not_null()
            & pl.col("input_pos").is_not_null()
            & (pl.col("input_pos") > 0)
            & pl.col("input_effect_allele").str.contains(r"^[ACGT]+$")
            & pl.col("input_other_allele").str.contains(r"^[ACGT]+$")
            & (pl.col("input_effect_allele") != pl.col("input_other_allele"))
        ).fill_null(False).alias("valid_input_variant"),
        (
            (pl.col("input_effect_allele").str.len_chars() == 1)
            & (pl.col("input_other_allele").str.len_chars() == 1)
        ).fill_null(False).alias("input_is_snp"),
        (
            pl.col("input_effect_allele").str.len_chars()
            != pl.col("input_other_allele").str.len_chars()
        ).fill_null(False).alias("input_is_indel"),
        pl.lit(eaf_is_maf, dtype=pl.Boolean).alias("input_frequency_is_maf"),
    ).with_columns(
        pl.when(~pl.col("valid_input_variant"))
        .then(pl.lit("invalid"))
        .when(pl.col("input_is_snp"))
        .then(pl.lit("snp"))
        .when(pl.col("input_is_indel"))
        .then(pl.lit("indel"))
        .otherwise(pl.lit("other"))
        .alias("input_variant_type"),
        (
            pl.col("input_is_snp")
            & (
                ((pl.col("input_effect_allele") == "A") & (pl.col("input_other_allele") == "T"))
                | ((pl.col("input_effect_allele") == "T") & (pl.col("input_other_allele") == "A"))
                | ((pl.col("input_effect_allele") == "C") & (pl.col("input_other_allele") == "G"))
                | ((pl.col("input_effect_allele") == "G") & (pl.col("input_other_allele") == "C"))
            )
        ).alias("input_is_palindromic"),
        _complement(pl.col("input_effect_allele")).alias("input_effect_complement"),
        _complement(pl.col("input_other_allele")).alias("input_other_complement"),
    ).with_columns(
        _canonical_key(
            "input_chrom", "input_pos", "input_effect_allele", "input_other_allele"
        ).alias("direct_key"),
        pl.when(pl.col("input_is_snp"))
        .then(_canonical_key(
            "input_chrom", "input_pos", "input_effect_complement", "input_other_complement"
        ))
        .otherwise(None)
        .alias("complement_key"),
    )
    if external_eaf_frame is not None:
        if external_eaf_mapping is None or not row.external_eaf_column:
            raise ValueError("External EAF data require a column mapping and frequency column.")
        external = (
            external_eaf_frame.select(
                _chromosome(pl.col(external_eaf_mapping.chromosome)).alias("external_chrom"),
                pl.col(external_eaf_mapping.position).cast(pl.Int64, strict=False).alias("external_pos"),
                _allele(pl.col(external_eaf_mapping.effect_allele)).alias("external_effect_allele"),
                _allele(pl.col(external_eaf_mapping.other_allele)).alias("external_other_allele"),
                pl.col(row.external_eaf_column).cast(pl.Float64, strict=False).alias("external_af"),
            )
            .with_columns(
                _canonical_key(
                    "external_chrom", "external_pos",
                    "external_effect_allele", "external_other_allele",
                ).alias("external_key")
            )
            .unique(subset=["external_key"], keep="first", maintain_order=True)
        )
        result = (
            result.join(external, left_on="direct_key", right_on="external_key", how="left")
            .with_columns(
                pl.when(pl.col("input_effect_allele") == pl.col("external_effect_allele"))
                .then(pl.col("external_af"))
                .when(pl.col("input_effect_allele") == pl.col("external_other_allele"))
                .then(1.0 - pl.col("external_af"))
                .otherwise(None)
                .alias("input_af")
            )
            .drop(
                "external_chrom", "external_pos", "external_effect_allele",
                "external_other_allele", "external_af",
            )
        )
    return result


def _apply_duplicate_report(
    input_table: pl.DataFrame,
    duplicate_report: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Apply harmonisation's recorded duplicate decisions to concordance input."""
    valid_report = duplicate_report.filter(pl.col("valid_input_variant"))
    invalid_report = duplicate_report.filter(~pl.col("valid_input_variant"))
    if invalid_report.height:
        raise ValueError(
            "The harmonisation duplicate report contains invalid coordinates or alleles."
        )

    # These labels are the complete schema-validated duplicate protocol:
    # consistent groups have one kept row and removed alternatives; conflicting
    # groups have no kept row and carry the configured dataset action.
    actions = {value for value in valid_report["duplicate_action"].unique() if value}
    unknown_actions = actions - {"kept", "removed", "remove_all", "fail_dataset"}
    if unknown_actions or valid_report["duplicate_action"].null_count():
        raise ValueError(
            "The harmonisation duplicate report contains missing or unrecognized "
            "duplicate_action values."
        )

    report_rows = valid_report.select("duplicate_input_row")
    if (
        valid_report["duplicate_input_row"].null_count()
        or valid_report["duplicate_input_row"].n_unique() != valid_report.height
    ):
        raise ValueError(
            "The harmonisation duplicate report contains missing or repeated "
            "duplicate_input_row values."
        )
    kept = valid_report.filter(pl.col("duplicate_action") == "kept")
    if kept.get_column("direct_key").n_unique() != kept.height:
        raise ValueError(
            "The harmonisation duplicate report retains more than one row for a variant."
        )
    group_actions = valid_report.group_by("direct_key").agg(
        pl.col("duplicate_action").unique().alias("_actions")
    )
    invalid_groups = group_actions.filter(
        (
            pl.col("_actions").list.contains("kept")
            & (
                pl.col("_actions").list.set_difference(pl.lit(["kept", "removed"]))
                .list.len() > 0
            )
        )
        | (
            ~pl.col("_actions").list.contains("kept")
            & (
                pl.col("_actions").list.set_difference(
                    pl.lit(["remove_all", "fail_dataset"])
                ).list.len() > 0
            )
        )
    )
    if invalid_groups.height:
        raise ValueError(
            "The harmonisation duplicate report contains an inconsistent set of "
            "actions for one or more variants."
        )

    match_columns = [
        "direct_key", "input_beta", "input_se", "input_z_supplied",
        "input_af", "input_p",
    ]
    verified_report = input_table.join(
        valid_report.select(["duplicate_input_row"] + match_columns),
        left_on=["input_row"] + match_columns,
        right_on=["duplicate_input_row"] + match_columns,
        how="inner",
        nulls_equal=True,
    )
    if verified_report.height != valid_report.height:
        raise ValueError(
            "The duplicate_input_row values in the harmonisation duplicate report "
            "do not match the source summary statistics."
        )

    selected = (
        input_table.join(
            kept.select("duplicate_input_row"),
            left_on="input_row",
            right_on="duplicate_input_row",
            how="inner",
        )
        .sort("input_row")
    )
    if selected.height != kept.height:
        raise ValueError(
            "The rows marked 'kept' in the harmonisation duplicate report do not "
            "match the source summary statistics."
        )

    duplicate_input = (
        input_table.join(
            report_rows, left_on="input_row", right_on="duplicate_input_row",
            how="inner",
        )
        .join(selected.select("input_row"), on="input_row", how="anti")
        .with_columns(pl.lit("duplicate_input_variant").alias("not_retained_reason"))
    )
    nonduplicates = input_table.join(
        report_rows, left_on="input_row", right_on="duplicate_input_row", how="anti",
    )
    retained = pl.concat([nonduplicates, selected], how="diagonal_relaxed").sort("input_row")
    return retained, duplicate_input


def prepare_vcf_table(frame: pl.DataFrame) -> pl.DataFrame:
    indexed = (
        frame
        if VCF_SOURCE_ROW_COLUMN in frame.columns
        else frame.with_row_index(VCF_SOURCE_ROW_COLUMN, offset=1)
    )
    result = (
        indexed
        .select(
            pl.col(VCF_SOURCE_ROW_COLUMN).cast(pl.UInt32).alias("vcf_row"),
            _chromosome(pl.col("CHROM")).alias("vcf_chrom"),
            pl.col("POS").cast(pl.Int64, strict=False).alias("vcf_pos"),
            pl.col("ID").cast(pl.String, strict=False).alias("vcf_id"),
            _allele(pl.col("REF")).alias("vcf_ref"),
            _allele(pl.col("ALT")).alias("vcf_alt"),
            pl.col("ES").cast(pl.Float64, strict=False).alias("vcf_effect"),
            pl.col("SE").cast(pl.Float64, strict=False).alias("vcf_se"),
            pl.col("EZ").cast(pl.Float64, strict=False).alias("vcf_z"),
            pl.col("AF").cast(pl.Float64, strict=False).alias("vcf_af"),
        )
        .with_columns(
            (
                pl.col("vcf_chrom").is_not_null()
                & pl.col("vcf_pos").is_not_null()
                & (pl.col("vcf_pos") > 0)
                & pl.col("vcf_ref").str.contains(r"^[ACGT]+$")
                & pl.col("vcf_alt").str.contains(r"^[ACGT]+$")
                & (pl.col("vcf_ref") != pl.col("vcf_alt"))
            ).fill_null(False).alias("valid_vcf_variant")
        )
        .with_columns(
            (
                (pl.col("vcf_ref").str.len_chars() == 1)
                & (pl.col("vcf_alt").str.len_chars() == 1)
            ).fill_null(False).alias("vcf_is_snp"),
            (
                pl.col("vcf_ref").str.len_chars()
                != pl.col("vcf_alt").str.len_chars()
            ).fill_null(False).alias("vcf_is_indel"),
            _canonical_key("vcf_chrom", "vcf_pos", "vcf_ref", "vcf_alt").alias("vcf_key")
        )
        .with_columns(
            pl.when(~pl.col("valid_vcf_variant"))
            .then(pl.lit("invalid"))
            .when(pl.col("vcf_is_snp"))
            .then(pl.lit("snp"))
            .when(pl.col("vcf_is_indel"))
            .then(pl.lit("indel"))
            .otherwise(pl.lit("other"))
            .alias("vcf_variant_type")
        )
    )
    return result


def _join_matches(input_unique: pl.DataFrame, vcf_unique: pl.DataFrame, allow_complement: bool):
    vcf_join = vcf_unique.with_columns(pl.col("vcf_key").alias("matched_vcf_key"))
    direct = input_unique.join(vcf_join, left_on="direct_key", right_on="vcf_key", how="left")
    direct_match = direct.filter(pl.col("vcf_row").is_not_null()).with_columns(
        pl.lit("listed").alias("_match_basis"), pl.lit(0).alias("_match_priority")
    )
    unmatched = direct.filter(pl.col("vcf_row").is_null()).select(input_unique.columns)

    complement_match = pl.DataFrame()
    if allow_complement and unmatched.height:
        complement = unmatched.filter(
            pl.col("input_is_snp") & ~pl.col("input_is_palindromic")
        ).join(vcf_join, left_on="complement_key", right_on="vcf_key", how="left")
        complement_match = complement.filter(pl.col("vcf_row").is_not_null()).with_columns(
            pl.lit("complement").alias("_match_basis"),
            pl.lit(1).alias("_match_priority"),
        )
        complemented_rows = complement_match.select("input_row")
        if complemented_rows.height:
            unmatched = unmatched.join(complemented_rows, on="input_row", how="anti")

    matched = direct_match
    if complement_match.height:
        matched = pl.concat([direct_match, complement_match], how="diagonal_relaxed")
    if matched.height:
        matched = matched.sort(["_match_priority", "input_row"])
        collisions = matched.filter(pl.col("matched_vcf_key").is_duplicated())
        matched = matched.unique(subset=["matched_vcf_key"], keep="first", maintain_order=True)
        if collisions.height:
            retained_rows = matched.select("input_row")
            collisions = collisions.join(retained_rows, on="input_row", how="anti")
            if collisions.height:
                unmatched = pl.concat(
                    [unmatched, collisions.select(input_unique.columns)],
                    how="diagonal_relaxed",
                ).unique(subset=["input_row"], keep="first")
    return matched, unmatched


def _metric_columns(frame: pl.DataFrame, name: str, expected: pl.Expr, observed: pl.Expr,
                    orientation: pl.Expr, tolerance) -> pl.DataFrame:
    checked = orientation & expected.is_not_null()
    difference = (expected - observed).abs()
    threshold = float(tolerance.absolute) + float(tolerance.relative) * expected.abs()
    return frame.with_columns(
        expected.alias("expected_%s" % name),
        difference.alias("%s_absolute_difference" % name),
        orientation.fill_null(False).alias("%s_orientation_eligible" % name),
        checked.fill_null(False).alias("%s_checked" % name),
        (checked & observed.is_not_null()).fill_null(False).alias("%s_observed" % name),
        (checked & observed.is_not_null() & (difference <= threshold))
        .fill_null(False).alias("%s_concordant" % name),
    )


def _safe_number(value):
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _metric_aggregate_expressions(name: str, observed: str) -> list[pl.Expr]:
    """Build one metric's expressions for a shared single-pass aggregation."""
    expected = "expected_%s" % name
    difference = "%s_absolute_difference" % name
    orientation_column = "%s_orientation_eligible" % name
    checked_column = "%s_checked" % name
    observed_column = "%s_observed" % name
    concordant_column = "%s_concordant" % name
    present = pl.col(observed_column)
    expected_present = pl.when(present).then(pl.col(expected)).otherwise(None)
    observed_present = pl.when(present).then(pl.col(observed)).otherwise(None)
    difference_present = pl.when(present).then(pl.col(difference)).otherwise(None)
    prefix = "%s__" % name
    expressions = [
        pl.len().alias(prefix + "rows"),
        pl.col(orientation_column).sum().alias(prefix + "orientation_eligible"),
        pl.col(checked_column).sum().alias(prefix + "checked"),
        pl.col(observed_column).sum().alias(prefix + "observed"),
        pl.col(concordant_column).sum().alias(prefix + "concordant"),
        difference_present.mean().alias(prefix + "mean"),
        difference_present.median().alias(prefix + "median"),
        difference_present.max().alias(prefix + "maximum"),
        difference_present.pow(2).mean().sqrt().alias(prefix + "rmse"),
        pl.corr(expected_present, observed_present, method="pearson")
        .alias(prefix + "pearson"),
        pl.corr(expected_present, observed_present, method="spearman")
        .alias(prefix + "spearman"),
    ]
    if name in ("effect", "z_score"):
        expressions.append(
            pl.when(present)
            .then(pl.col(expected).sign() == pl.col(observed).sign())
            .otherwise(None)
            .mean()
            .alias(prefix + "sign")
        )
    return expressions


def _metric_summary_from_stats(stats: dict[str, Any], name: str) -> dict[str, Any]:
    prefix = "%s__" % name
    rows = int(stats.get(prefix + "rows") or 0)
    orientation_eligible = int(stats.get(prefix + "orientation_eligible") or 0)
    checked = int(stats.get(prefix + "checked") or 0)
    observed_count = int(stats.get(prefix + "observed") or 0)
    concordant = int(stats.get(prefix + "concordant") or 0)
    return {
        "orientation_eligible": orientation_eligible,
        "unavailable_in_input": orientation_eligible - checked,
        "excluded_orientation": rows - orientation_eligible,
        "checked": checked,
        "observed": observed_count,
        "missing_in_vcf": checked - observed_count,
        "concordant": concordant,
        "mismatches": checked - concordant,
        "concordance_fraction": concordant / checked if checked else None,
        "pearson": (
            _safe_number(stats.get(prefix + "pearson"))
            if observed_count >= 2 else None
        ),
        "spearman": (
            _safe_number(stats.get(prefix + "spearman"))
            if observed_count >= 2 else None
        ),
        "sign_concordance": _safe_number(stats.get(prefix + "sign")),
        "mean_absolute_difference": _safe_number(stats.get(prefix + "mean")),
        "median_absolute_difference": _safe_number(stats.get(prefix + "median")),
        "maximum_absolute_difference": _safe_number(stats.get(prefix + "maximum")),
        "rmse": _safe_number(stats.get(prefix + "rmse")),
    }


def _metric_summaries(
    frame: pl.DataFrame | pl.LazyFrame,
    *,
    eaf_is_maf: bool | None,
    group_column: str | None = None,
) -> dict[str, dict[str, Any]] | dict[str, dict[str, dict[str, Any]]]:
    """Aggregate every configured value metric in one eager or lazy scan."""
    lazy = frame.lazy() if isinstance(frame, pl.DataFrame) else frame
    expressions = [
        expression
        for name, observed in METRIC_OBSERVED_COLUMNS.items()
        for expression in _metric_aggregate_expressions(name, observed)
    ]
    if group_column is None:
        records = lazy.select(expressions).collect().to_dicts()
        stats = records[0] if records else {}
        summaries = {
            name: _metric_summary_from_stats(stats, name)
            for name in METRIC_OBSERVED_COLUMNS
        }
        summaries["allele_frequency"]["frequency_type"] = (
            "minor_allele_frequency"
            if eaf_is_maf is True
            else "effect_allele_frequency"
            if eaf_is_maf is False
            else "unresolved"
        )
        return summaries

    records = {
        str(record[group_column]): record
        for record in lazy.group_by(group_column).agg(expressions).collect().to_dicts()
    }
    by_type: dict[str, dict[str, dict[str, Any]]] = {}
    for label, variant_type in VARIANT_TYPES:
        stats = records.get(variant_type, {})
        summaries = {
            name: _metric_summary_from_stats(stats, name)
            for name in METRIC_OBSERVED_COLUMNS
        }
        summaries["allele_frequency"]["frequency_type"] = (
            "minor_allele_frequency"
            if eaf_is_maf is True
            else "effect_allele_frequency"
            if eaf_is_maf is False
            else "unresolved"
        )
        by_type[label] = summaries
    return by_type


def _counts_by_type(frame: pl.DataFrame, type_column: str) -> dict[str, int]:
    """Count valid rows by the already-computed variant-type label."""
    if not frame.height or type_column not in frame.columns:
        return {}
    return {
        str(variant_type): int(count)
        for variant_type, count in (
            frame.filter(pl.col(type_column).is_not_null())
            .group_by(type_column)
            .len()
            .iter_rows()
        )
    }


def _position_diagnostics(
    unmatched_input: pl.DataFrame,
    unmatched_vcf: pl.DataFrame,
    *,
    settings,
) -> tuple[
    pl.DataFrame,
    pl.DataFrame,
    pl.DataFrame,
    dict[str, dict[str, int]],
]:
    """Compare unmatched one-to-one records at the same position diagnostically.

    Allele correspondence is unknown once exact allele matching has failed.
    Position-level effect and Z comparisons therefore use magnitudes, and AF is
    folded to minor frequency. These diagnostics never establish orientation
    and never determine pass/fail status.
    """
    input_keys = ["input_chrom", "input_pos"]
    vcf_keys = ["vcf_chrom", "vcf_pos"]
    input_positions = unmatched_input.group_by(input_keys).agg(
        pl.len().alias("input_records_at_position"),
        pl.col("input_variant_type").first().alias("input_position_variant_type"),
    )
    vcf_positions = unmatched_vcf.group_by(vcf_keys).agg(
        pl.len().alias("vcf_records_at_position"),
        pl.col("vcf_variant_type").first().alias("vcf_position_variant_type"),
    )
    shared_positions = input_positions.join(
        vcf_positions,
        left_on=input_keys,
        right_on=vcf_keys,
        how="inner",
    )

    annotated_input = (
        unmatched_input.join(input_positions, on=input_keys, how="left")
        .join(vcf_positions, left_on=input_keys, right_on=vcf_keys, how="left")
        .with_columns(pl.col("vcf_records_at_position").fill_null(0))
        .with_columns(
            pl.when(pl.col("vcf_records_at_position") == 0)
            .then(pl.lit("input_position_only"))
            .when(
                (pl.col("input_records_at_position") == 1)
                & (pl.col("vcf_records_at_position") == 1)
                & (
                    pl.col("input_variant_type")
                    == pl.col("vcf_position_variant_type")
                )
            )
            .then(pl.lit("one_to_one_same_position"))
            .when(
                (pl.col("input_records_at_position") == 1)
                & (pl.col("vcf_records_at_position") == 1)
            )
            .then(pl.lit("variant_type_mismatch"))
            .otherwise(pl.lit("multiallelic_position_ambiguous"))
            .alias("position_comparison")
        )
    )
    annotated_vcf = (
        unmatched_vcf.join(vcf_positions, on=vcf_keys, how="left")
        .join(input_positions, left_on=vcf_keys, right_on=input_keys, how="left")
        .with_columns(pl.col("input_records_at_position").fill_null(0))
        .with_columns(
            pl.when(pl.col("input_records_at_position") == 0)
            .then(pl.lit("vcf_position_only"))
            .when(
                (pl.col("input_records_at_position") == 1)
                & (pl.col("vcf_records_at_position") == 1)
                & (
                    pl.col("vcf_variant_type")
                    == pl.col("input_position_variant_type")
                )
            )
            .then(pl.lit("one_to_one_same_position"))
            .when(
                (pl.col("input_records_at_position") == 1)
                & (pl.col("vcf_records_at_position") == 1)
            )
            .then(pl.lit("variant_type_mismatch"))
            .otherwise(pl.lit("multiallelic_position_ambiguous"))
            .alias("position_comparison")
        )
    )

    position_pairs = annotated_input.filter(
        pl.col("position_comparison") == "one_to_one_same_position"
    ).join(
        annotated_vcf.filter(
            pl.col("position_comparison") == "one_to_one_same_position"
        ).select(unmatched_vcf.columns),
        left_on=input_keys,
        right_on=vcf_keys,
        how="inner",
    )
    orientation_unknown = pl.lit(True)
    position_pairs = _metric_columns(
        position_pairs,
        "effect",
        pl.col("input_beta").abs(),
        pl.col("vcf_effect").abs(),
        orientation_unknown,
        settings.effect,
    )
    position_pairs = _metric_columns(
        position_pairs,
        "standard_error",
        pl.col("input_se"),
        pl.col("vcf_se"),
        orientation_unknown,
        settings.standard_error,
    )
    position_pairs = position_pairs.with_columns(
        pl.min_horizontal(pl.col("input_af"), 1.0 - pl.col("input_af"))
        .alias("_position_expected_af"),
        pl.min_horizontal(pl.col("vcf_af"), 1.0 - pl.col("vcf_af"))
        .alias("observed_frequency"),
    )
    position_pairs = _metric_columns(
        position_pairs,
        "allele_frequency",
        pl.col("_position_expected_af"),
        pl.col("observed_frequency"),
        orientation_unknown,
        settings.allele_frequency,
    ).drop("_position_expected_af")
    position_pairs = _metric_columns(
        position_pairs,
        "z_score",
        pl.col("input_z").abs(),
        pl.col("vcf_z").abs(),
        orientation_unknown,
        settings.z_score,
    )
    position_pairs = position_pairs.with_columns(
        pl.concat_str([
            pl.when(pl.col("effect_checked") & ~pl.col("effect_concordant"))
            .then(pl.lit("effect_mismatch")).otherwise(pl.lit("")),
            pl.when(
                pl.col("standard_error_checked")
                & ~pl.col("standard_error_concordant")
            ).then(pl.lit("standard_error_mismatch")).otherwise(pl.lit("")),
            pl.when(
                pl.col("allele_frequency_checked")
                & ~pl.col("allele_frequency_concordant")
            ).then(pl.lit("allele_frequency_mismatch")).otherwise(pl.lit("")),
            pl.when(pl.col("z_score_checked") & ~pl.col("z_score_concordant"))
            .then(pl.lit("z_score_mismatch")).otherwise(pl.lit("")),
        ], separator=";")
        .str.replace_all(r";+", ";")
        .str.strip_chars(";")
        .alias("position_value_reasons")
    )

    counts: dict[str, dict[str, int]] = {}
    checked_any = (
        pl.col("effect_checked")
        | pl.col("standard_error_checked")
        | pl.col("allele_frequency_checked")
        | pl.col("z_score_checked")
    )
    checked_pairs = position_pairs.filter(checked_any)
    mismatch_pairs = checked_pairs.filter(pl.col("position_value_reasons") != "")
    global_pair_positions = position_pairs.select(input_keys).unique().height
    global_input_specific_positions = annotated_input.filter(
        pl.col("vcf_records_at_position") == 0
    ).select(input_keys).unique().height
    global_vcf_specific_positions = annotated_vcf.filter(
        pl.col("input_records_at_position") == 0
    ).select(vcf_keys).unique().height
    counts["all_variants"] = {
        "input_unmatched_positions": input_positions.height,
        "vcf_unmatched_positions": vcf_positions.height,
        "unmatched_position_union": (
            input_positions.height + vcf_positions.height - shared_positions.height
        ),
        "shared_unmatched_positions": shared_positions.height,
        "unambiguous_position_pairs": global_pair_positions,
        "ambiguous_shared_positions": shared_positions.height - global_pair_positions,
        "variant_type_mismatch_positions": shared_positions.filter(
            (pl.col("input_records_at_position") == 1)
            & (pl.col("vcf_records_at_position") == 1)
            & (
                pl.col("input_position_variant_type")
                != pl.col("vcf_position_variant_type")
            )
        ).height,
        "multiallelic_ambiguous_positions": shared_positions.filter(
            (pl.col("input_records_at_position") != 1)
            | (pl.col("vcf_records_at_position") != 1)
        ).height,
        "same_position_input_variants": annotated_input.filter(
            pl.col("vcf_records_at_position") > 0
        ).height,
        "same_position_vcf_variants": annotated_vcf.filter(
            pl.col("input_records_at_position") > 0
        ).height,
        "input_specific_variants": annotated_input.filter(
            pl.col("vcf_records_at_position") == 0
        ).height,
        "vcf_specific_variants": annotated_vcf.filter(
            pl.col("input_records_at_position") == 0
        ).height,
        "input_specific_positions": global_input_specific_positions,
        "vcf_specific_positions": global_vcf_specific_positions,
        "position_value_checked_pairs": checked_pairs.height,
        "position_value_concordant_pairs": (
            checked_pairs.height - mismatch_pairs.height
        ),
        "position_value_mismatch_pairs": mismatch_pairs.height,
    }

    for label, variant_type in VARIANT_TYPES:
        input_type = annotated_input.filter(
            pl.col("input_variant_type") == variant_type
        )
        vcf_type = annotated_vcf.filter(pl.col("vcf_variant_type") == variant_type)
        input_type_positions = input_type.select(input_keys).unique()
        vcf_type_positions = vcf_type.select(vcf_keys).unique()
        shared_type = input_type_positions.join(
            vcf_type_positions,
            left_on=input_keys,
            right_on=vcf_keys,
            how="inner",
        )
        pairs_type = position_pairs.filter(
            pl.col("input_variant_type") == variant_type
        )
        checked_pairs = pairs_type.filter(checked_any)
        mismatch_pairs = checked_pairs.filter(pl.col("position_value_reasons") != "")
        counts[label] = {
            "input_unmatched_positions": input_type_positions.height,
            "vcf_unmatched_positions": vcf_type_positions.height,
            "unmatched_position_union": (
                input_type_positions.height
                + vcf_type_positions.height
                - shared_type.height
            ),
            "shared_unmatched_positions": shared_type.height,
            "unambiguous_position_pairs": pairs_type.height,
            "ambiguous_shared_positions": shared_type.height - pairs_type.height,
            "same_position_input_variants": input_type.filter(
                pl.col("vcf_records_at_position") > 0
            ).height,
            "same_position_vcf_variants": vcf_type.filter(
                pl.col("input_records_at_position") > 0
            ).height,
            "input_specific_variants": input_type.filter(
                pl.col("vcf_records_at_position") == 0
            ).height,
            "vcf_specific_variants": vcf_type.filter(
                pl.col("input_records_at_position") == 0
            ).height,
            "input_specific_positions": input_type.filter(
                pl.col("vcf_records_at_position") == 0
            ).select("input_chrom", "input_pos").unique().height,
            "vcf_specific_positions": vcf_type.filter(
                pl.col("input_records_at_position") == 0
            ).select("vcf_chrom", "vcf_pos").unique().height,
            "position_value_checked_pairs": checked_pairs.height,
            "position_value_concordant_pairs": (
                checked_pairs.height - mismatch_pairs.height
            ),
            "position_value_mismatch_pairs": mismatch_pairs.height,
        }
    return annotated_input, annotated_vcf, position_pairs, counts


def _variant_type_summary(
    *,
    input_count: int,
    vcf_count: int,
    matched_count: int,
    input_only_count: int,
    vcf_only_count: int,
    duplicate_count: int,
    checked_count: int,
    mismatch_count: int,
    match_type_counts: dict[str, int] | None = None,
    position_counts: dict[str, int] | None = None,
) -> dict[str, Any]:
    union_count = input_count + vcf_count - matched_count
    summary = {
        "input_unique_variants": input_count,
        "vcf_unique_variants": vcf_count,
        "exact_matched_variants": matched_count,
        "input_only_variants": input_only_count,
        "vcf_only_variants": vcf_only_count,
        "vcf_duplicate_records": duplicate_count,
        "exact_match_fraction": matched_count / input_count if input_count else None,
        "variant_union": union_count,
        "common_union_fraction": matched_count / union_count if union_count else None,
        "value_checked_variants": checked_count,
        "value_mismatch_variants": mismatch_count,
        "value_mismatch_fraction": (
            mismatch_count / checked_count if checked_count else 0.0
        ),
    }
    supplied_matches = match_type_counts or {}
    summary.update({
        field: int(supplied_matches.get(match_type, 0))
        for match_type, field in MATCH_TYPE_FIELDS.items()
    })
    supplied_positions = position_counts or {}
    summary.update({
        field: int(supplied_positions.get(field, 0))
        for field in POSITION_COUNT_FIELDS
    })
    return summary


def _variant_type_counts(
    input_unique: pl.DataFrame,
    unmatched_input: pl.DataFrame,
    vcf_unique: pl.DataFrame,
    vcf_only: pl.DataFrame,
    vcf_duplicates: pl.DataFrame,
    matched: pl.DataFrame,
    position_counts: dict[str, dict[str, int]],
) -> dict[str, dict[str, Any]]:
    """Build exact-match statistics independently for SNPs and indels."""
    input_counts = _counts_by_type(input_unique, "input_variant_type")
    input_only_counts = _counts_by_type(unmatched_input, "input_variant_type")
    vcf_counts = _counts_by_type(vcf_unique, "vcf_variant_type")
    valid_vcf_only = vcf_only.filter(pl.col("valid_vcf_variant"))
    vcf_only_counts = _counts_by_type(valid_vcf_only, "vcf_variant_type")
    duplicate_counts = _counts_by_type(vcf_duplicates, "vcf_variant_type")
    matched_counts = _counts_by_type(matched, "input_variant_type")
    checked_counts: dict[str, int] = {}
    mismatch_counts: dict[str, int] = {}
    match_counts: dict[str, dict[str, int]] = {}
    if matched.height and "effect_checked" in matched.columns:
        checked = matched.filter(
            pl.col("effect_checked")
            | pl.col("standard_error_checked")
            | pl.col("allele_frequency_checked")
            | pl.col("z_score_checked")
        )
        checked_counts = _counts_by_type(checked, "input_variant_type")
        mismatch_counts = _counts_by_type(
            matched.filter(pl.col("failure_reasons") != ""), "input_variant_type"
        )
        for variant_type, match_type, count in matched.group_by(
            "input_variant_type", "match_type"
        ).len().iter_rows():
            match_counts.setdefault(str(variant_type), {})[str(match_type)] = int(count)

    summaries: dict[str, dict[str, Any]] = {}
    for label, variant_type in VARIANT_TYPES:
        input_count = input_counts.get(variant_type, 0)
        matched_count = matched_counts.get(variant_type, 0)
        input_only_count = input_only_counts.get(variant_type, 0)
        vcf_only_count = vcf_only_counts.get(variant_type, 0)
        checked_count = checked_counts.get(variant_type, 0)
        mismatch_count = mismatch_counts.get(variant_type, 0)
        summaries[label] = _variant_type_summary(
            input_count=input_count,
            vcf_count=vcf_counts.get(variant_type, 0),
            matched_count=matched_count,
            input_only_count=input_only_count,
            vcf_only_count=vcf_only_count,
            duplicate_count=duplicate_counts.get(variant_type, 0),
            checked_count=checked_count,
            mismatch_count=mismatch_count,
            match_type_counts=match_counts.get(variant_type),
            position_counts=position_counts.get(label),
        )
    return summaries


def _evaluate_variant_type(
    label: str,
    summary: dict[str, Any],
    settings,
) -> tuple[str, list[str]]:
    """Apply value-integrity thresholds; unmatched variants are report-only."""
    if not summary["input_unique_variants"] and not summary["vcf_unique_variants"]:
        return "NOT_APPLICABLE", ["no %s were available for comparison" % label]

    failure_reasons: list[str] = []
    warning_reasons: list[str] = []
    if summary["vcf_duplicate_records"] > int(
        settings.failure.maximum_vcf_duplicate_records
    ):
        failure_reasons.append(
            "duplicate %s VCF records exceed the configured maximum" % label
        )
    if summary["value_mismatch_fraction"] > float(
        settings.failure.maximum_value_mismatch_fraction
    ):
        failure_reasons.append(
            "%s value mismatches exceed the configured maximum fraction" % label
        )
    if summary["input_only_variants"]:
        warning_reasons.append(
            "some valid input %s have no allele-aware VCF match (report-only)" % label
        )
    if summary["vcf_only_variants"]:
        warning_reasons.append(
            "some VCF %s have no allele-aware input match (report-only)" % label
        )
    if summary["position_value_mismatch_pairs"]:
        warning_reasons.append(
            "same-position unmatched %s contain diagnostic value differences"
            % label
        )

    if failure_reasons:
        return "FAIL", failure_reasons
    if warning_reasons:
        return "WARNING", warning_reasons
    return "PASS", ["all configured %s checks passed" % label]


def _finalise_variant_types(
    variant_types: dict[str, dict[str, Any]],
    settings,
) -> tuple[str, list[str], str, list[str], str, list[str]]:
    snp_status, snp_reasons = _evaluate_variant_type(
        "SNPs", variant_types["snps"], settings,
    )
    indel_status, indel_reasons = _evaluate_variant_type(
        "indels",
        variant_types["indels"],
        settings,
    )
    other_status, other_reasons = _evaluate_variant_type(
        "other variants", variant_types["other_variants"], settings,
    )
    variant_types["snps"].update(status=snp_status, status_reasons=snp_reasons)
    variant_types["indels"].update(status=indel_status, status_reasons=indel_reasons)
    variant_types["other_variants"].update(
        status=other_status,
        status_reasons=other_reasons,
    )
    return (
        snp_status, snp_reasons,
        indel_status, indel_reasons,
        other_status, other_reasons,
    )


def _finalise_summary_status(
    summary: dict[str, Any],
    metric_summary: dict[str, dict[str, Any]],
    *,
    eaf_is_maf: bool | None,
    settings,
) -> str:
    statuses = _finalise_variant_types(summary["variant_types"], settings)
    snp_status, snp_reasons, indel_status, indel_reasons, other_status, other_reasons = statuses
    unavailable_values = sum(
        values["unavailable_in_input"] for values in metric_summary.values()
    )
    failure_reasons = [
        "%s: %s" % (label, reason)
        for label, status, reasons in (
            ("SNP concordance", snp_status, snp_reasons),
            ("indel concordance", indel_status, indel_reasons),
            ("other-variant concordance", other_status, other_reasons),
        )
        if status == "FAIL"
        for reason in reasons
    ]
    if summary["vcf_invalid_rows"] > int(
        settings.failure.maximum_invalid_vcf_records
    ):
        failure_reasons.append(
            "invalid VCF records exceed the configured maximum"
        )
    warning_reasons = [
        "%s: %s" % (label, reason)
        for label, status, reasons in (
            ("SNP concordance", snp_status, snp_reasons),
            ("indel concordance", indel_status, indel_reasons),
            ("other-variant concordance", other_status, other_reasons),
        )
        if status == "WARNING"
        for reason in reasons
    ]
    if summary["input_invalid_rows"]:
        warning_reasons.append("some input rows have invalid coordinates or alleles")
    if summary["input_duplicate_rows"]:
        warning_reasons.append("the input contains duplicate variants")
    if unavailable_values:
        warning_reasons.append("some requested value comparisons lack an input value")
    if summary["palindromic_variants_excluded_from_orientation"]:
        warning_reasons.append(
            "some palindromic variants lacked a safely resolved orientation for "
            "effect, EAF, and Z comparison"
        )
    if eaf_is_maf is None:
        warning_reasons.append(
            "allele-frequency concordance was skipped because the final frequency "
            "type is unresolved"
        )
    if summary["vcf_invalid_rows"]:
        warning_reasons.append("some VCF rows have invalid coordinates or alleles")

    if failure_reasons:
        status, status_reasons = "FAIL", failure_reasons
    elif warning_reasons:
        status, status_reasons = "WARNING", warning_reasons
    else:
        status, status_reasons = "PASS", ["all configured checks passed"]
    summary["status"] = status
    summary["status_reasons"] = status_reasons
    return status


def compare_input_to_vcf(
    input_frame: pl.DataFrame,
    vcf_frame: pl.DataFrame,
    row,
    *,
    effect_type: str,
    se_scale: str | None = None,
    p_value_type: str,
    eaf_is_maf: bool | None,
    strand_consensus: str | None = None,
    settings,
    policies,
    external_eaf_frame: pl.DataFrame | None = None,
    external_eaf_mapping=None,
    duplicate_report_frame: pl.DataFrame | None = None,
) -> ConcordanceAnalysis:
    input_table = prepare_input_table(
        input_frame, row, effect_type=effect_type, se_scale=se_scale,
        p_value_type=p_value_type,
        eaf_is_maf=eaf_is_maf, settings=settings, policies=policies,
        external_eaf_frame=external_eaf_frame,
        external_eaf_mapping=external_eaf_mapping,
    )
    vcf_table = prepare_vcf_table(vcf_frame)

    input_invalid = input_table.filter(~pl.col("valid_input_variant")).with_columns(
        pl.lit("invalid_coordinate_or_allele").alias("not_retained_reason")
    )
    input_eligible = input_table.filter(pl.col("valid_input_variant")).with_columns(
        pl.len().over("direct_key").alias("input_duplicate_count")
    )
    if duplicate_report_frame is not None:
        prepared_report = prepare_input_table(
            duplicate_report_frame, row,
            effect_type=effect_type, se_scale=se_scale,
            p_value_type=p_value_type,
            eaf_is_maf=eaf_is_maf, settings=settings, policies=policies,
            external_eaf_frame=external_eaf_frame,
            external_eaf_mapping=external_eaf_mapping,
        )
        input_unique, duplicate_input = _apply_duplicate_report(
            input_eligible, prepared_report,
        )
    else:
        input_unique = input_eligible.unique(
            subset=["direct_key"], keep="first", maintain_order=True
        )
        retained_first = input_unique.select(
            "direct_key", pl.col("input_row").alias("_retained_input_row")
        )
        duplicate_input = (
            input_eligible.join(retained_first, on="direct_key", how="left")
            .filter(pl.col("input_row") != pl.col("_retained_input_row"))
            .drop("_retained_input_row")
            .with_columns(pl.lit("duplicate_input_variant").alias("not_retained_reason"))
        )

    vcf_invalid = vcf_table.filter(~pl.col("valid_vcf_variant")).with_columns(
        pl.lit("invalid_vcf_coordinate_or_allele").alias("vcf_only_reason")
    )
    vcf_eligible = vcf_table.filter(pl.col("valid_vcf_variant")).with_columns(
        pl.len().over("vcf_key").alias("vcf_duplicate_count")
    )
    vcf_unique = vcf_eligible.unique(subset=["vcf_key"], keep="first", maintain_order=True)
    retained_vcf_rows = vcf_unique.select(
        "vcf_key", pl.col("vcf_row").alias("_retained_vcf_row")
    )
    vcf_duplicates = (
        vcf_eligible.join(retained_vcf_rows, on="vcf_key", how="left")
        .filter(pl.col("vcf_row") != pl.col("_retained_vcf_row"))
        .drop("_retained_vcf_row")
        .with_columns(pl.lit("duplicate_vcf_variant").alias("vcf_duplicate_reason"))
    )

    matched, unmatched_unique = _join_matches(
        input_unique, vcf_unique, bool(settings.allow_strand_complement)
    )
    listed_effect_is_alt = (
        pl.when(pl.col("_match_basis") == "listed")
        .then(pl.col("input_effect_allele") == pl.col("vcf_alt"))
        .otherwise(pl.col("input_effect_complement") == pl.col("vcf_alt"))
    )
    if (
        settings.palindromic_action == "compare_resolved"
        and strand_consensus == "reverse"
    ):
        effect_is_alt = pl.when(pl.col("input_is_palindromic")).then(
            pl.col("input_effect_complement") == pl.col("vcf_alt")
        ).otherwise(listed_effect_is_alt)
    else:
        effect_is_alt = listed_effect_is_alt
    matched = matched.with_columns(effect_is_alt.alias("_effect_is_alt"))

    palindromic_resolved = (
        settings.palindromic_action == "compare_resolved"
        and strand_consensus in ("forward", "reverse")
    )
    if palindromic_resolved:
        palindromic_aligned = "palindromic_%s" % (
            "forward" if strand_consensus == "forward" else "reverse_complement"
        )
        palindromic_swapped = palindromic_aligned + "_swapped"
    else:
        palindromic_aligned = "palindromic_as_listed"
        palindromic_swapped = "palindromic_swapped"

    matched = matched.with_columns(
        pl.when(pl.col("input_is_palindromic"))
        .then(
            pl.when(pl.col("_effect_is_alt"))
            .then(pl.lit(palindromic_aligned))
            .otherwise(pl.lit(palindromic_swapped))
        )
        .when(pl.col("_match_basis") == "complement")
        .then(
            pl.when(pl.col("_effect_is_alt"))
            .then(pl.lit("strand_complement"))
            .otherwise(pl.lit("strand_complement_swapped"))
        )
        .otherwise(
            pl.when(pl.col("_effect_is_alt"))
            .then(pl.lit("direct"))
            .otherwise(pl.lit("allele_swapped"))
        )
        .alias("match_type"),
        pl.when(pl.col("_effect_is_alt"))
        .then(pl.lit(1.0))
        .otherwise(pl.lit(-1.0))
        .alias("orientation_factor"),
    )
    include_palindromic = (
        settings.palindromic_action == "compare_as_listed"
        or palindromic_resolved
    )
    orientation = (
        pl.lit(True) if include_palindromic else ~pl.col("input_is_palindromic")
    )
    matched = matched.with_columns(orientation.alias("orientation_comparable"))
    matched = _metric_columns(
        matched, "effect",
        pl.col("input_beta") * pl.col("orientation_factor"),
        pl.col("vcf_effect"), pl.col("orientation_comparable"), settings.effect,
    )
    matched = _metric_columns(
        matched, "standard_error",
        pl.col("input_se"), pl.col("vcf_se"), pl.lit(True),
        settings.standard_error,
    )
    if eaf_is_maf is True:
        expected_af = pl.min_horizontal(pl.col("input_af"), 1.0 - pl.col("input_af"))
        observed_af = pl.min_horizontal(pl.col("vcf_af"), 1.0 - pl.col("vcf_af"))
        af_orientation = pl.lit(True)
    elif eaf_is_maf is False:
        expected_af = pl.when(pl.col("orientation_factor") > 0).then(
            pl.col("input_af")
        ).otherwise(1.0 - pl.col("input_af"))
        observed_af = pl.col("vcf_af")
        af_orientation = pl.col("orientation_comparable")
    else:
        expected_af = pl.lit(None, dtype=pl.Float64)
        observed_af = pl.col("vcf_af")
        af_orientation = pl.lit(False)
    matched = matched.with_columns(observed_af.alias("observed_frequency"))
    matched = _metric_columns(
        matched, "allele_frequency", expected_af,
        pl.col("observed_frequency"), af_orientation, settings.allele_frequency,
    )
    matched = _metric_columns(
        matched, "z_score",
        pl.col("input_z") * pl.col("orientation_factor"),
        pl.col("vcf_z"), pl.col("orientation_comparable"), settings.z_score,
    )
    matched = matched.with_columns(
        pl.concat_str([
            pl.when(pl.col("effect_checked") & ~pl.col("effect_concordant"))
            .then(pl.lit("effect_mismatch")).otherwise(pl.lit("")),
            pl.when(
                pl.col("standard_error_checked")
                & ~pl.col("standard_error_concordant")
            ).then(pl.lit("standard_error_mismatch")).otherwise(pl.lit("")),
            pl.when(pl.col("allele_frequency_checked") & ~pl.col("allele_frequency_concordant"))
            .then(pl.lit("allele_frequency_mismatch")).otherwise(pl.lit("")),
            pl.when(pl.col("z_score_checked") & ~pl.col("z_score_concordant"))
            .then(pl.lit("z_score_mismatch")).otherwise(pl.lit("")),
        ], separator=";")
        .str.replace_all(r";+", ";")
        .str.strip_chars(";")
        .alias("failure_reasons")
    )

    matched_vcf = matched.select("matched_vcf_key")
    valid_vcf_only = vcf_unique.join(
        matched_vcf, left_on="vcf_key", right_on="matched_vcf_key", how="anti"
    ).with_columns(pl.lit("no_input_coordinate_and_allele_match").alias("vcf_only_reason"))
    (
        unmatched_unique,
        valid_vcf_only,
        position_matches,
        position_counts,
    ) = _position_diagnostics(
        unmatched_unique,
        valid_vcf_only,
        settings=settings,
    )
    unmatched_unique = unmatched_unique.with_columns(
        pl.lit("no_coordinate_and_allele_match").alias("not_retained_reason")
    )
    input_only = pl.concat(
        [input_invalid, duplicate_input, unmatched_unique], how="diagonal_relaxed"
    ).sort("input_row")
    vcf_only = valid_vcf_only
    if vcf_invalid.height:
        vcf_only = pl.concat([vcf_only, vcf_invalid], how="diagonal_relaxed")

    mismatches = matched.filter(pl.col("failure_reasons") != "")
    metric_summary = _metric_summaries(matched, eaf_is_maf=eaf_is_maf)
    metric_summary_by_variant_type = _metric_summaries(
        matched,
        eaf_is_maf=eaf_is_maf,
        group_column="input_variant_type",
    )
    position_metric_summary_by_variant_type = _metric_summaries(
        position_matches,
        eaf_is_maf=eaf_is_maf,
        group_column="input_variant_type",
    )

    checked_any = 0
    if matched.height:
        checked_any = int(matched.select(
            (
                pl.col("effect_checked")
                | pl.col("standard_error_checked")
                | pl.col("allele_frequency_checked")
                | pl.col("z_score_checked")
            ).sum()
        ).item() or 0)
    mismatch_fraction = mismatches.height / checked_any if checked_any else 0.0
    retained_fraction = matched.height / input_unique.height if input_unique.height else 0.0
    match_counts = (
        dict(matched.group_by("match_type").len().iter_rows()) if matched.height else {}
    )
    variant_types = _variant_type_counts(
        input_unique,
        unmatched_unique,
        vcf_unique,
        valid_vcf_only,
        vcf_duplicates,
        matched,
        position_counts,
    )
    palindromic_total = sum(
        int(count)
        for match_type, count in match_counts.items()
        if str(match_type).startswith("palindromic_")
    )
    palindromic_compared = palindromic_total if include_palindromic else 0
    palindromic_excluded = palindromic_total - palindromic_compared
    if palindromic_resolved:
        palindromic_basis = "study_wide_%s_strand_consensus" % strand_consensus
    elif settings.palindromic_action == "compare_as_listed":
        palindromic_basis = "listed_allele_order"
    elif settings.palindromic_action == "exclude":
        palindromic_basis = "excluded_by_policy"
    else:
        palindromic_basis = "strand_consensus_unavailable"

    summary = {
        "input_rows": input_table.height,
        "input_valid_rows": input_eligible.height,
        "input_unique_variants": input_unique.height,
        "input_duplicate_rows": duplicate_input.height,
        "input_invalid_rows": input_invalid.height,
        "vcf_records": vcf_table.height,
        "vcf_unique_variants": vcf_unique.height,
        "vcf_duplicate_records": vcf_duplicates.height,
        "vcf_invalid_rows": vcf_invalid.height,
        "matched_variants": matched.height,
        "variant_union": input_unique.height + vcf_unique.height - matched.height,
        "input_only_variants": unmatched_unique.height,
        "vcf_only_variants": valid_vcf_only.height,
        "retained_fraction": retained_fraction,
        "value_mismatch_variants": mismatches.height,
        "value_mismatch_fraction": mismatch_fraction,
        "palindromic_variants_compared": palindromic_compared,
        "palindromic_variants_excluded_from_orientation": palindromic_excluded,
        "palindromic_comparison_basis": palindromic_basis,
        "match_types": match_counts,
        "variant_types": variant_types,
        "position_diagnostics": position_counts["all_variants"],
    }
    summary["common_union_fraction"] = (
        summary["matched_variants"] / summary["variant_union"]
        if summary["variant_union"] else None
    )
    status = _finalise_summary_status(
        summary,
        metric_summary,
        eaf_is_maf=eaf_is_maf,
        settings=settings,
    )
    return ConcordanceAnalysis(
        status=status,
        summary=summary,
        metric_summary=metric_summary,
        metric_summary_by_variant_type=metric_summary_by_variant_type,
        position_metric_summary_by_variant_type=(
            position_metric_summary_by_variant_type
        ),
        matched=matched,
        position_matches=position_matches,
        mismatches=mismatches,
        input_only=input_only,
        vcf_only=vcf_only,
        vcf_duplicates=vcf_duplicates,
    )


def combine_concordance_partitions(
    partition_summaries: list[dict[str, Any]],
    *,
    matched: pl.LazyFrame,
    position_matches: pl.LazyFrame,
    input_only: pl.LazyFrame,
    vcf_only: pl.LazyFrame,
    vcf_duplicates: pl.LazyFrame,
    eaf_is_maf: bool | None,
    settings,
) -> ConcordanceAnalysis:
    """Combine chromosome analyses without materializing their report tables."""
    count_fields = (
        "input_rows",
        "input_valid_rows",
        "input_unique_variants",
        "input_duplicate_rows",
        "input_invalid_rows",
        "vcf_records",
        "vcf_unique_variants",
        "vcf_duplicate_records",
        "vcf_invalid_rows",
        "matched_variants",
        "input_only_variants",
        "vcf_only_variants",
        "value_mismatch_variants",
        "palindromic_variants_compared",
        "palindromic_variants_excluded_from_orientation",
    )
    summary = {
        field: sum(int(partition[field]) for partition in partition_summaries)
        for field in count_fields
    }
    summary["retained_fraction"] = (
        summary["matched_variants"] / summary["input_unique_variants"]
        if summary["input_unique_variants"] else 0.0
    )
    summary["variant_union"] = (
        summary["input_unique_variants"]
        + summary["vcf_unique_variants"]
        - summary["matched_variants"]
    )
    summary["common_union_fraction"] = (
        summary["matched_variants"] / summary["variant_union"]
        if summary["variant_union"] else None
    )
    checked_any = int(
        matched.select(
            (
                pl.col("effect_checked")
                | pl.col("standard_error_checked")
                | pl.col("allele_frequency_checked")
                | pl.col("z_score_checked")
            ).sum().alias("checked_any")
        ).collect().item()
        or 0
    )
    summary["value_mismatch_fraction"] = (
        summary["value_mismatch_variants"] / checked_any if checked_any else 0.0
    )
    match_types: dict[str, int] = {}
    for partition in partition_summaries:
        for match_type, count in partition["match_types"].items():
            match_types[match_type] = match_types.get(match_type, 0) + int(count)
    summary["match_types"] = match_types
    bases = {
        str(partition["palindromic_comparison_basis"])
        for partition in partition_summaries
    }
    summary["palindromic_comparison_basis"] = (
        bases.pop() if len(bases) == 1 else "mixed"
    )

    variant_types: dict[str, dict[str, Any]] = {}
    for label, _variant_type in VARIANT_TYPES:
        values = [partition["variant_types"][label] for partition in partition_summaries]
        combined_match_counts = {
            match_type: sum(int(value[field]) for value in values)
            for match_type, field in MATCH_TYPE_FIELDS.items()
        }
        combined_position_counts = {
            field: sum(int(value[field]) for value in values)
            for field in POSITION_COUNT_FIELDS
        }
        variant_types[label] = _variant_type_summary(
            input_count=sum(int(value["input_unique_variants"]) for value in values),
            vcf_count=sum(int(value["vcf_unique_variants"]) for value in values),
            matched_count=sum(int(value["exact_matched_variants"]) for value in values),
            input_only_count=sum(int(value["input_only_variants"]) for value in values),
            vcf_only_count=sum(int(value["vcf_only_variants"]) for value in values),
            duplicate_count=sum(int(value["vcf_duplicate_records"]) for value in values),
            checked_count=sum(int(value["value_checked_variants"]) for value in values),
            mismatch_count=sum(int(value["value_mismatch_variants"]) for value in values),
            match_type_counts=combined_match_counts,
            position_counts=combined_position_counts,
        )
    summary["variant_types"] = variant_types
    position_fields = partition_summaries[0]["position_diagnostics"].keys()
    summary["position_diagnostics"] = {
        field: sum(
            int(partition["position_diagnostics"][field])
            for partition in partition_summaries
        )
        for field in position_fields
    }

    metric_summary = _metric_summaries(matched, eaf_is_maf=eaf_is_maf)
    metric_summary_by_variant_type = _metric_summaries(
        matched,
        eaf_is_maf=eaf_is_maf,
        group_column="input_variant_type",
    )
    position_metric_summary_by_variant_type = _metric_summaries(
        position_matches,
        eaf_is_maf=eaf_is_maf,
        group_column="input_variant_type",
    )
    status = _finalise_summary_status(
        summary,
        metric_summary,
        eaf_is_maf=eaf_is_maf,
        settings=settings,
    )
    ordered_matches = matched.sort(["_match_priority", "input_row"])
    return ConcordanceAnalysis(
        status=status,
        summary=summary,
        metric_summary=metric_summary,
        metric_summary_by_variant_type=metric_summary_by_variant_type,
        position_metric_summary_by_variant_type=(
            position_metric_summary_by_variant_type
        ),
        matched=ordered_matches,
        position_matches=position_matches.sort(["input_chrom", "input_pos"]),
        mismatches=ordered_matches.filter(pl.col("failure_reasons") != ""),
        input_only=input_only.sort("input_row"),
        vcf_only=vcf_only.sort(
            ["valid_vcf_variant", "vcf_row"], descending=[True, False],
        ),
        vcf_duplicates=vcf_duplicates.sort("vcf_row"),
    )
