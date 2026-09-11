"""Reference-based allele orientation for one chromosome.

The supplied reference remains a normal CHROM/POS/REF/ALT/population-AF table.
One coordinate join is followed by vectorised tests of direct/swapped allele
order and, for SNVs, the two reverse-complement forms; no fourfold reference
table is materialised. Study-wide strand consensus and, when necessary, a true
internal study effect-allele frequency arbitrate palindromic SNPs, whose allele
letters cannot identify the strand on their own.

This is an internal PostGWAS step. Effect estimates must already be on the beta
scale: the preceding effect-normalization step converts OR to log(OR), rejects
non-positive OR under the resolved policy, and records that provenance. A swap
therefore negates beta; raw odds ratios are never reciprocated here.

Scientific methods: NHGRI-EBI GWAS Catalog palindrome orientation and the
frequency-informed TwoSampleMR harmonisation method:
https://ebispot.github.io/gwas-sumstats-harmoniser-documentation/Introduction/Orientation-of-palindromic-variants/
https://mrcieu.github.io/TwoSampleMR/articles/harmonise.html

The study-frequency guard follows the GWAS-SSF and GWAS-VCF requirement that
association frequency belongs to the listed effect/ALT allele. It compares the
declared EAF interpretation with its one-minus alternative on uniquely oriented
non-palindromic variants and fails only when all configured evidence thresholds
support the non-effect-allele interpretation:
https://www.ebi.ac.uk/gwas/docs/summary-statistics-format
https://github.com/MRCIEU/gwas-vcf-specification
"""

import math
from typing import Any, Dict, Optional, Tuple

import polars as pl

from postgwas.core.values import optional_text

from .shared.allele_join import deduplicate_reference_values
from .shared.runtime import reject_rows, resolve_policies
from .shared.statistics import frequency_maf_screen, valid_frequency_mask
from .shared.variant_columns import (
    canonicalize_variant_frame,
    has_canonical_variant_columns,
    mark_canonical_variant_columns,
    palindromic_snp_expression,
    read_reference_variant_table,
    reverse_complement_expression,
)

__all__ = [
    "STRAND_ACTION_COLUMN",
    "REFERENCE_AF_COLUMN",
    "RESOLVED_STRAND_ACTIONS",
    "REFERENCE_UNMATCHED_RETAINED_ACTIONS",
    "EXPORTABLE_STRAND_ACTIONS",
    "StrandOrientationError",
    "resolve_strand_consensus",
    "assess_study_frequency_allele_alignment",
    "resolve_palindromic_frequency_evidence",
]


STRAND_ACTION_COLUMN = "strand_action"
REFERENCE_AF_COLUMN = "strand_reference_af"
STEP_LABEL = "04 strand_orientation"
RESOLVED_STRAND_ACTIONS = (
    "forward",
    "forward_swapped",
    "reverse_complement",
    "reverse_complement_swapped",
)
REFERENCE_UNMATCHED_RETAINED_ACTIONS = (
    "reference_unmatched_retained",
    "reference_unmatched_retained_reverse_complement",
)
EXPORTABLE_STRAND_ACTIONS = (
    RESOLVED_STRAND_ACTIONS + REFERENCE_UNMATCHED_RETAINED_ACTIONS
)
POLICY_KEYS = (
    "strand.mode",
    "strand.consensus_threshold",
    "strand.min_informative_variants",
    "strand.palindromic_af_ambiguity_lower",
    "strand.palindromic_af_ambiguity_upper",
    "strand.palindromic_af_max_difference",
    "strand.palindromic_af_min_error_margin",
    "strand.unmatched_action",
    "strand.ambiguous_action",
    "strand.af_tolerance",
    "strand.af_discordance_action",
    "strand.palindromic_af_discordance_action",
    "eaf.non_effect_frequency_min_overlap",
    "eaf.non_effect_frequency_min_correlation",
    "eaf.non_effect_frequency_max_error",
    "eaf.non_effect_frequency_error_margin",
    "external_reference.exact_duplicate_action",
    "external_reference.non_identical_duplicate_action",
    "chromosome.strip_chr_prefix",
    "chromosome.strip_leading_zero",
    "chromosome.rename_map",
)


class StrandOrientationError(RuntimeError):
    """A chromosome cannot be oriented safely against its reference."""


def resolve_strand_consensus(genome_build_info, build, policies=None):
    """Resolve one study-wide forward/reverse decision from build evidence."""
    pol = resolve_policies(policies)
    mode = str(pol.get("strand.mode"))
    if mode == "reference_aligned":
        return {
            "strand": "reference_aligned",
            "strand_mode": mode,
            "forward": 0,
            "reverse": 0,
            "ambiguous": 0,
            "informative": 0,
            "dominant_fraction": None,
            "minimum_informative_variants": int(
                pol.get("strand.min_informative_variants")
            ),
            "consensus_threshold": float(pol.get("strand.consensus_threshold")),
            "reason": (
                "the configured strand.mode explicitly declares reference-aligned "
                "study alleles; chromosome rows will still be validated against REF/ALT"
            ),
        }
    evidence = dict((genome_build_info.get("strand_evidence") or {}).get(build) or {})
    forward = int(evidence.get("forward", 0) or 0)
    reverse = int(evidence.get("reverse", 0) or 0)
    ambiguous = int(evidence.get("ambiguous", 0) or 0)
    informative = forward + reverse
    minimum = int(pol.get("strand.min_informative_variants"))
    threshold = float(pol.get("strand.consensus_threshold"))
    if informative == 0:
        decision = "unresolved"
        dominant_fraction = None
        reason = "no non-palindromic reference matches"
    else:
        dominant_fraction = max(forward, reverse) / informative
        if informative < minimum:
            decision = "unresolved"
            reason = "%s informative variants are below the configured minimum %s" % (
                "{:,}".format(informative), "{:,}".format(minimum),
            )
        elif dominant_fraction < threshold:
            decision = "mixed"
            reason = "dominant strand fraction %.4f is below %.4f" % (
                dominant_fraction, threshold,
            )
        else:
            decision = "forward" if forward >= reverse else "reverse"
            reason = "dominant strand fraction %.4f" % dominant_fraction
    return {
        "strand": decision,
        "strand_mode": mode,
        "forward": forward,
        "reverse": reverse,
        "ambiguous": ambiguous,
        "informative": informative,
        "dominant_fraction": dominant_fraction,
        "minimum_informative_variants": minimum,
        "consensus_threshold": threshold,
        "reason": reason,
    }


def _load_reference(
    path, population_col, colmap, policies, chromosome, *, warn=None,
):
    raw, _detected = read_reference_variant_table(
        path,
        colmap,
        policies,
        value_columns=[population_col],
        error_type=StrandOrientationError,
        description="chromosome %s strand reference" % chromosome,
    )
    selected = raw.rename({
        colmap["chr"]: "__strand_chr",
        colmap["pos"]: "__strand_pos",
        colmap["a1"]: "__strand_alt",
        colmap["a2"]: "__strand_ref",
        population_col: "__strand_af",
    })
    reference, _normalization = canonicalize_variant_frame(
        selected,
        {
            "chr": "__strand_chr",
            "pos": "__strand_pos",
            "ea": "__strand_alt",
            "oa": "__strand_ref",
        },
        policies,
        error_type=StrandOrientationError,
        warn=warn,
        label="strand reference",
    )
    reference = reference.with_columns(
        pl.col("__strand_af").cast(pl.Float64, strict=False)
    )
    invalid = reference.filter(
        pl.col("__strand_af").is_not_null()
        & ((pl.col("__strand_af") < 0.0) | (pl.col("__strand_af") > 1.0))
    ).height
    if invalid:
        raise StrandOrientationError(
            "The strand reference '%s' contains %s values outside 0 to 1 in column '%s'."
            % (path, "{:,}".format(invalid), population_col)
        )
    reference, duplicate_stats = deduplicate_reference_values(
        reference,
        ["__strand_chr", "__strand_pos", "__strand_ref", "__strand_alt"],
        "__strand_af",
        exact_action=str(
            policies.get("external_reference.exact_duplicate_action")
        ),
        non_identical_action=str(
            policies.get("external_reference.non_identical_duplicate_action")
        ),
        reference_label="strand reference '%s'" % path,
        error_type=StrandOrientationError,
    )
    duplicate_stats["raw_frequency_screen"] = frequency_maf_screen(
        reference,
        "__strand_af",
        float(policies.get("eaf.maf_decision_cutoff")),
    )
    return reference, duplicate_stats


def _study_value(decision, key, default=None):
    if key == "eaf_is_maf" and isinstance(decision, bool):
        return decision
    if isinstance(decision, dict):
        return decision.get(key, default)
    return default


def _finite_float(value):
    """Return one finite diagnostic scalar, otherwise ``None``."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def assess_study_frequency_allele_alignment(
    candidates: pl.DataFrame,
    *,
    row_column: str,
    frequency_column: str,
    policies=None,
) -> Dict[str, Any]:
    """Test whether a declared study EAF instead follows the non-effect allele.

    ``candidates`` is the already materialised result of the chromosome's
    coordinate/reference join. Only non-palindromic study rows with exactly one
    valid allele orientation contribute. For an unswapped row, the expected
    frequency of the study's listed effect allele is reference ALT AF; for a
    swapped row it is ``1 - ALT_AF``. Reverse complementation does not change
    frequency. No row is transformed here.

    A non-effect diagnosis requires sufficient overlap, strong positive
    correlation after one-minus inversion, a small inverted mean absolute
    error, and the configured improvement over the declared-EAF error. An
    ancestry-mismatched or uninformative comparison therefore remains
    inconclusive instead of being corrected or rejected.
    """
    pol = resolve_policies(policies)
    minimum_overlap = int(pol.get("eaf.non_effect_frequency_min_overlap"))
    minimum_correlation = float(
        pol.get("eaf.non_effect_frequency_min_correlation")
    )
    maximum_inverted_error = float(
        pol.get("eaf.non_effect_frequency_max_error")
    )
    minimum_error_margin = float(
        pol.get("eaf.non_effect_frequency_error_margin")
    )
    evidence: Dict[str, Any] = {
        "status": "inconclusive",
        "decision": "inconclusive",
        "comparable_variants": 0,
        "minimum_overlap": minimum_overlap,
        "declared_eaf_correlation": None,
        "inverted_eaf_correlation": None,
        "declared_eaf_mean_absolute_error": None,
        "inverted_eaf_mean_absolute_error": None,
        "minimum_inverted_correlation": minimum_correlation,
        "maximum_inverted_mean_absolute_error": maximum_inverted_error,
        "minimum_error_margin": minimum_error_margin,
        "decision_reason": "no uniquely oriented non-palindromic evidence",
    }
    if candidates.is_empty() or frequency_column not in candidates.columns:
        return evidence

    counted_candidates = (
        candidates
        if "__strand_candidates" in candidates.columns
        else candidates.with_columns(
            pl.len().over(row_column).alias("__strand_candidates")
        )
    )
    unique_non_palindromic = counted_candidates.filter(
        ~pl.col("__strand_palindromic")
        & (pl.col("__strand_candidates") == 1)
    )
    study_frequency = pl.col(frequency_column).cast(pl.Float64, strict=False)
    reference_alt_frequency = pl.col("__strand_af").cast(
        pl.Float64, strict=False
    )
    expected_effect_frequency = (
        pl.when(pl.col("__strand_swapped"))
        .then(1.0 - reference_alt_frequency)
        .otherwise(reference_alt_frequency)
    )
    comparable = (
        unique_non_palindromic.select([
            study_frequency.alias("__study_frequency"),
            expected_effect_frequency.alias("__expected_effect_frequency"),
        ])
        .filter(
            valid_frequency_mask(pl.col("__study_frequency"))
            & valid_frequency_mask(pl.col("__expected_effect_frequency"))
        )
    )
    if comparable.is_empty():
        return evidence

    study = pl.col("__study_frequency")
    expected = pl.col("__expected_effect_frequency")
    inverted = 1.0 - study
    statistics = comparable.select([
        pl.len().alias("comparable_variants"),
        pl.corr(study, expected, method="pearson").alias(
            "declared_correlation"
        ),
        pl.corr(inverted, expected, method="pearson").alias(
            "inverted_correlation"
        ),
        (study - expected).abs().mean().alias("declared_error"),
        (inverted - expected).abs().mean().alias("inverted_error"),
    ]).to_dicts()[0]
    count = int(statistics["comparable_variants"] or 0)
    declared_correlation = _finite_float(statistics["declared_correlation"])
    inverted_correlation = _finite_float(statistics["inverted_correlation"])
    declared_error = _finite_float(statistics["declared_error"])
    inverted_error = _finite_float(statistics["inverted_error"])
    evidence.update({
        "comparable_variants": count,
        "declared_eaf_correlation": declared_correlation,
        "inverted_eaf_correlation": inverted_correlation,
        "declared_eaf_mean_absolute_error": declared_error,
        "inverted_eaf_mean_absolute_error": inverted_error,
    })

    if count < minimum_overlap:
        evidence["decision_reason"] = "insufficient_overlap"
        return evidence
    if None in (
        declared_correlation,
        inverted_correlation,
        declared_error,
        inverted_error,
    ):
        evidence["decision_reason"] = "non_finite_or_constant_frequency_evidence"
        return evidence

    inverted_supported = (
        inverted_correlation >= minimum_correlation
        and declared_correlation <= -minimum_correlation
        and inverted_error <= maximum_inverted_error
        and inverted_error + minimum_error_margin <= declared_error
    )
    declared_supported = (
        declared_correlation >= minimum_correlation
        and declared_error <= maximum_inverted_error
        and declared_error + minimum_error_margin <= inverted_error
    )
    if inverted_supported:
        evidence.update({
            "status": "failed",
            "decision": "non_effect_allele_frequency",
            "decision_reason": "inverted_interpretation_passed_all_thresholds",
        })
    elif declared_supported:
        evidence.update({
            "status": "passed",
            "decision": "effect_allele_frequency",
            "decision_reason": "declared_interpretation_passed_all_thresholds",
        })
    else:
        evidence["decision_reason"] = "evidence_did_not_support_one_interpretation"
    return evidence


def resolve_palindromic_frequency_evidence(
    candidates: pl.DataFrame,
    *,
    row_column: str,
    candidate_column: str,
    frequency_column: str,
    reference_frequency_column: str,
    policies=None,
) -> pl.DataFrame:
    """Select one palindromic orientation only when AF evidence is decisive.

    ``candidates`` must contain every competing orientation for each row.  The
    returned frame has one row per input row, a ``palindromic_frequency_state``
    value, and the winning candidate/frequency only for a scientifically unique
    result.  The ambiguity interval, maximum reference difference, and required
    separation from the competing orientation all come from canonical YAML.
    """
    pol = resolve_policies(policies)
    lower = float(pol.get("strand.palindromic_af_ambiguity_lower"))
    upper = float(pol.get("strand.palindromic_af_ambiguity_upper"))
    maximum_difference = float(
        pol.get("strand.palindromic_af_max_difference")
    )
    minimum_margin = float(
        pol.get("strand.palindromic_af_min_error_margin")
    )
    output_columns = [
        row_column,
        "palindromic_frequency_state",
        "palindromic_frequency_candidate",
        "palindromic_frequency_value",
        "palindromic_frequency_error",
    ]
    if candidates.is_empty():
        return pl.DataFrame(schema={
            row_column: candidates.schema[row_column],
            "palindromic_frequency_state": pl.Utf8,
            "palindromic_frequency_candidate": pl.Utf8,
            "palindromic_frequency_value": pl.Float64,
            "palindromic_frequency_error": pl.Float64,
        }).select(output_columns)

    frequency = pl.col(frequency_column).cast(pl.Float64, strict=False)
    reference = pl.col(reference_frequency_column).cast(
        pl.Float64, strict=False
    )
    frequency_valid = (
        frequency.is_not_null()
        & frequency.is_finite()
        & (frequency >= 0.0)
        & (frequency <= 1.0)
    ).fill_null(False)
    reference_valid = (
        reference.is_not_null()
        & reference.is_finite()
        & (reference >= 0.0)
        & (reference <= 1.0)
    ).fill_null(False)
    frequency_outside = ((frequency < lower) | (frequency > upper)).fill_null(
        False
    )
    reference_outside = ((reference < lower) | (reference > upper)).fill_null(
        False
    )
    informative = frequency_valid & reference_valid & frequency_outside & reference_outside

    scored = candidates.with_columns([
        frequency.alias("__pal_frequency"),
        reference.alias("__pal_reference_frequency"),
        frequency_valid.alias("__pal_frequency_valid"),
        reference_valid.alias("__pal_reference_valid"),
        informative.alias("__pal_frequency_informative"),
        pl.when(informative)
        .then((frequency - reference).abs())
        .otherwise(None)
        .alias("__pal_frequency_error"),
    ]).with_columns([
        pl.len().over(row_column).alias("__pal_candidate_count"),
        pl.col("__pal_frequency_informative")
        .sum()
        .over(row_column)
        .alias("__pal_informative_count"),
        pl.col("__pal_frequency_error")
        .min()
        .over(row_column)
        .alias("__pal_min_error"),
    ]).with_columns(
        (
            (
                pl.col("__pal_frequency_error")
                < pl.col("__pal_min_error") + minimum_margin
            )
            .fill_null(False)
            .sum()
            .over(row_column)
        ).alias("__pal_competing_best_count")
    )

    all_informative = (
        pl.col("__pal_informative_count") == pl.col("__pal_candidate_count")
    )
    winner = (
        all_informative
        & (pl.col("__pal_frequency_error") == pl.col("__pal_min_error"))
        & (pl.col("__pal_min_error") <= maximum_difference)
        & (pl.col("__pal_competing_best_count") == 1)
    ).fill_null(False)
    winners = scored.filter(winner).select([
        row_column,
        pl.col(candidate_column)
        .cast(pl.Utf8)
        .alias("palindromic_frequency_candidate"),
        pl.col("__pal_frequency").alias("palindromic_frequency_value"),
        pl.col("__pal_frequency_error").alias("palindromic_frequency_error"),
    ])

    states = scored.group_by(row_column, maintain_order=True).agg([
        pl.col("__pal_candidate_count").first().alias("__pal_candidate_count"),
        pl.col("__pal_informative_count").first().alias("__pal_informative_count"),
        pl.col("__pal_min_error").first().alias("__pal_min_error"),
        pl.col("__pal_competing_best_count")
        .first()
        .alias("__pal_competing_best_count"),
        (~pl.col("__pal_frequency_valid") | ~pl.col("__pal_reference_valid"))
        .any()
        .alias("__pal_has_unusable"),
        (
            pl.col("__pal_frequency_valid")
            & pl.col("__pal_reference_valid")
            & ~pl.col("__pal_frequency_informative")
        )
        .any()
        .alias("__pal_near_half"),
    ]).with_columns(
        pl.when(
            (pl.col("__pal_informative_count") == pl.col("__pal_candidate_count"))
            & (pl.col("__pal_min_error") <= maximum_difference)
            & (pl.col("__pal_competing_best_count") == 1)
        )
        .then(pl.lit("resolved"))
        .when(pl.col("__pal_has_unusable"))
        .then(pl.lit("unusable_frequency"))
        .when(pl.col("__pal_near_half"))
        .then(pl.lit("ambiguous_frequency_band"))
        .when(pl.col("__pal_min_error") > maximum_difference)
        .then(pl.lit("maximum_difference_exceeded"))
        .otherwise(pl.lit("insufficient_error_margin"))
        .alias("palindromic_frequency_state")
    )
    return states.join(
        winners, on=row_column, how="left", coalesce=True
    ).select(output_columns)


def harmonise_strand_orientation(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    reference_file: str,
    population_col: str,
    reference_colmap: Dict[str, str],
    study_decision: Optional[Any] = None,
    policies=None,
    logger=None,
    ctx=None,
    rejects=None,
) -> Tuple[pl.DataFrame, Dict[str, Any], dict]:
    """Orient one chromosome after effect estimates have reached beta scale."""
    pol = resolve_policies(policies)
    mode = str(pol.get("strand.mode"))
    if not reference_file or not population_col or not reference_colmap:
        raise StrandOrientationError(
            "Strand orientation requires a raw chromosome reference file, its population AF "
            "column, and the configured CHROM/POS/REF/ALT mapping."
        )

    chr_col = sample_column_dict["chr_col"]
    pos_col = sample_column_dict["pos_col"]
    ea_col = sample_column_dict["ea_col"]
    oa_col = sample_column_dict["oa_col"]
    eaf_col = optional_text(sample_column_dict.get("eaf_col"))
    eaf_decision = _study_value(study_decision, "eaf_is_maf")
    if eaf_col is not None and eaf_col in df.columns and eaf_decision is None:
        raise StrandOrientationError(
            "The study-level EAF/MAF decision is missing; strand orientation cannot "
            "safely decide whether swapped frequencies should be inverted."
        )
    eaf_is_maf = eaf_decision is True
    effect_col = optional_text(sample_column_dict.get("beta_or_col"))
    effect_type = _study_value(study_decision, "effect_type")
    if effect_col and effect_col in df.columns and effect_type not in (
        "beta", "odds_ratio"
    ):
        raise StrandOrientationError(
            "The study-level effect decision is missing; strand orientation cannot "
            "verify that the effect column has already been normalized to beta."
        )
    if effect_col and effect_col in df.columns and effect_type == "odds_ratio":
        raise StrandOrientationError(
            "Internal PostGWAS execution-order error: strand orientation received "
            "raw odds ratios on chromosome %s. Call harmonise_effect_estimates() "
            "first so non-positive OR values follow effect.or_non_positive with "
            "effect_non_positive_or provenance and retained effects are converted "
            "to beta=log(OR) before allele orientation."
            % chromosome
        )
    row_col = "__postgwas_strand_row"
    while row_col in df.columns:
        row_col += "_"

    if ctx is not None:
        ctx.input(
            "Chromosome strand reference",
            chromosome=chromosome,
            file=str(reference_file),
            population_frequency_column=population_col,
            delimiter=reference_colmap["delimiter"],
            chromosome_column=reference_colmap["chr"],
            position_column=reference_colmap["pos"],
            reference_allele_column=reference_colmap["a2"],
            alternate_allele_column=reference_colmap["a1"],
        )
        ctx.observed(
            "Applied study strand consensus",
            value=str(_study_value(study_decision, "strand", "unresolved") or "unresolved"),
        )

    reference, reference_duplicate_stats = _load_reference(
        reference_file,
        population_col,
        reference_colmap,
        pol,
        chromosome,
        warn=ctx.warn if ctx is not None else None,
    )
    if ctx is not None and reference_duplicate_stats["reference_duplicate_groups"]:
        ctx.info(
            "Strand reference duplicate groups: %s exact and %s non-identical; "
            "%s groups discarded by the configured shared reference policy."
            % (
                reference_duplicate_stats["reference_exact_duplicate_groups"],
                reference_duplicate_stats[
                    "reference_non_identical_duplicate_groups"
                ],
                reference_duplicate_stats[
                    "reference_duplicate_groups_discarded"
                ],
            )
        )
    study = df.with_row_index(row_col)
    if not has_canonical_variant_columns(df, sample_column_dict):
        study, _normalization = canonicalize_variant_frame(
            study,
            {"chr": chr_col, "pos": pos_col, "ea": ea_col, "oa": oa_col},
            pol,
            error_type=StrandOrientationError,
            warn=ctx.warn if ctx is not None else None,
            label="study before strand-reference matching",
        )
    joined = study.join(
        reference,
        left_on=[chr_col, pos_col],
        right_on=["__strand_chr", "__strand_pos"],
        how="inner",
    ).with_columns([
        reverse_complement_expression(
            pl.col(ea_col)
        ).alias("__strand_ea_rc"),
        reverse_complement_expression(
            pl.col(oa_col)
        ).alias("__strand_oa_rc"),
    ])

    ea = pl.col(ea_col)
    oa = pl.col(oa_col)
    ref = pl.col("__strand_ref")
    alt = pl.col("__strand_alt")
    snv = (
        (ea.str.len_chars() == 1)
        & (oa.str.len_chars() == 1)
        & (ref.str.len_chars() == 1)
        & (alt.str.len_chars() == 1)
    )
    orientations = [
        ("forward", (ea == alt) & (oa == ref), False),
        ("forward_swapped", (ea == ref) & (oa == alt), True),
        (
            "reverse_complement",
            snv
            & (pl.col("__strand_ea_rc") == alt)
            & (pl.col("__strand_oa_rc") == ref),
            False,
        ),
        (
            "reverse_complement_swapped",
            snv
            & (pl.col("__strand_ea_rc") == ref)
            & (pl.col("__strand_oa_rc") == alt),
            True,
        ),
    ]
    if mode == "reference_aligned":
        orientations = orientations[:2]
    candidate_frames = [
        joined.filter(condition).with_columns([
            pl.lit(action).alias(STRAND_ACTION_COLUMN),
            pl.lit(swapped).alias("__strand_swapped"),
        ])
        for action, condition, swapped in orientations
    ]
    candidates = pl.concat(candidate_frames, how="vertical_relaxed").unique(
        subset=[
            row_col, STRAND_ACTION_COLUMN, "__strand_ref", "__strand_alt", "__strand_af"
        ],
        maintain_order=True,
    )

    candidates = candidates.with_columns([
        palindromic_snp_expression(
            pl.col(ea_col), pl.col(oa_col)
        ).alias("__strand_palindromic"),
        pl.len().over(row_col).alias("__strand_candidates"),
    ])
    candidate_counts = candidates.select([
        row_col, "__strand_candidates",
    ]).unique(subset=[row_col], maintain_order=True)
    frequency_alignment: Dict[str, Any] = {
        "status": "not_applicable",
        "decision": "not_evaluated",
        "decision_reason": (
            "no study frequency column"
            if eaf_col is None or eaf_col not in candidates.columns
            else "study frequency was screened as MAF-like"
            if eaf_is_maf
            else "frequency alignment check was not applicable"
        ),
    }
    if (
        eaf_col is not None
        and eaf_col in candidates.columns
        and not eaf_is_maf
    ):
        frequency_alignment = assess_study_frequency_allele_alignment(
            candidates,
            row_column=row_col,
            frequency_column=eaf_col,
            policies=pol,
        )
        if ctx is not None:
            ctx.decide(
                "Study frequency allele alignment",
                {
                    "column": eaf_col,
                    "uniquely oriented non-palindromic variants": (
                        frequency_alignment["comparable_variants"]
                    ),
                    "minimum required": frequency_alignment["minimum_overlap"],
                    "correlation as declared EAF": frequency_alignment[
                        "declared_eaf_correlation"
                    ],
                    "correlation after one-minus inversion": frequency_alignment[
                        "inverted_eaf_correlation"
                    ],
                    "mean error as declared EAF": frequency_alignment[
                        "declared_eaf_mean_absolute_error"
                    ],
                    "mean error after one-minus inversion": frequency_alignment[
                        "inverted_eaf_mean_absolute_error"
                    ],
                    "decision reason": frequency_alignment["decision_reason"],
                },
                frequency_alignment["decision"],
            )
        if frequency_alignment["decision"] == "non_effect_allele_frequency":
            raise StrandOrientationError(
                "The supplied frequency column '%s' appears to contain the frequency "
                "of the non-effect allele rather than the listed effect allele on "
                "chromosome %s. Across %s uniquely oriented non-palindromic variants, "
                "mean absolute error was %.6g as declared and %.6g after one-minus "
                "inversion; the corresponding correlations were %.6g and %.6g. "
                "PostGWAS did not modify the data. Supply the correct effect-allele-"
                "frequency column or a dataset-supplied external EAF file."
                % (
                    eaf_col,
                    chromosome,
                    frequency_alignment["comparable_variants"],
                    frequency_alignment["declared_eaf_mean_absolute_error"],
                    frequency_alignment["inverted_eaf_mean_absolute_error"],
                    frequency_alignment["declared_eaf_correlation"],
                    frequency_alignment["inverted_eaf_correlation"],
                )
            )
    consensus = str(_study_value(study_decision, "strand", "unresolved") or "unresolved")

    candidate_strand = (
        pl.when(pl.col(STRAND_ACTION_COLUMN).str.starts_with("forward"))
        .then(pl.lit("forward"))
        .otherwise(pl.lit("reverse"))
    )
    candidates = candidates.with_columns(candidate_strand.alias("__strand_basis"))

    pal_candidates = candidates.filter(pl.col("__strand_palindromic"))
    if (
        mode == "auto"
        and eaf_col is not None
        and eaf_col in candidates.columns
        and not eaf_is_maf
    ):
        candidate_eaf_col = "__strand_candidate_eaf"
        candidate_eaf = pl.col(eaf_col).cast(pl.Float64, strict=False)
        pal_candidates = pal_candidates.with_columns(
            pl.when(pl.col("__strand_swapped"))
            .then(1.0 - candidate_eaf)
            .otherwise(candidate_eaf)
            .alias(candidate_eaf_col)
        )
        pal_frequency = resolve_palindromic_frequency_evidence(
            pal_candidates,
            row_column=row_col,
            candidate_column=STRAND_ACTION_COLUMN,
            frequency_column=candidate_eaf_col,
            reference_frequency_column="__strand_af",
            policies=pol,
        )
    else:
        pal_frequency = pal_candidates.select(row_col).unique().with_columns([
            pl.lit("not_available").alias("palindromic_frequency_state"),
            pl.lit(None, dtype=pl.Utf8).alias(
                "palindromic_frequency_candidate"
            ),
            pl.lit(None, dtype=pl.Float64).alias(
                "palindromic_frequency_value"
            ),
            pl.lit(None, dtype=pl.Float64).alias(
                "palindromic_frequency_error"
            ),
        ])

    consensus_known = consensus in ("forward", "reverse")
    if mode == "reference_aligned":
        declared_counts = pal_candidates.group_by(row_col).len().rename({
            "len": "__strand_declared_candidates"
        })
        selected_pal = (
            pal_candidates.join(
                declared_counts, on=row_col, how="left", coalesce=True,
            )
            .filter(pl.col("__strand_declared_candidates") == 1)
            .with_columns(
                pl.lit("reference_aligned_declaration")
                .alias("__strand_resolution_basis")
            )
        )
        frequency_conflicts = candidates.head(0).select(row_col).with_columns(
            pl.lit(False).alias("__strand_frequency_conflict")
        )
    else:
        consensus_pool = (
            candidates.filter(
                pl.col("__strand_palindromic")
                & (pl.col("__strand_basis") == pl.lit(consensus))
            )
            if consensus_known
            else candidates.head(0)
        )
        consensus_counts = consensus_pool.group_by(row_col).len().rename({
            "len": "__strand_consensus_candidates"
        })
        consensus_selected = consensus_pool.join(
            consensus_counts, on=row_col, how="left", coalesce=True
        ).filter(pl.col("__strand_consensus_candidates") == 1).join(
            pal_frequency, on=row_col, how="left", coalesce=True
        )
        frequency_conflicts = consensus_selected.filter(
            (pl.col("palindromic_frequency_state") == "resolved")
            & (
                pl.col("palindromic_frequency_candidate")
                != pl.col(STRAND_ACTION_COLUMN)
            )
        ).select(row_col).unique().with_columns(
            pl.lit(True).alias("__strand_frequency_conflict")
        )
        consensus_selected = consensus_selected.join(
            frequency_conflicts, on=row_col, how="left", coalesce=True
        ).filter(
            ~pl.col("__strand_frequency_conflict").fill_null(False)
        ).with_columns(
            pl.when(pl.col("palindromic_frequency_state") == "resolved")
            .then(pl.lit("consensus_and_internal_eaf"))
            .otherwise(pl.lit("study_wide_consensus"))
            .alias("__strand_resolution_basis")
        )

        frequency_selected = pal_candidates.join(
            pal_frequency.filter(
                pl.col("palindromic_frequency_state") == "resolved"
            ),
            left_on=[row_col, STRAND_ACTION_COLUMN],
            right_on=[row_col, "palindromic_frequency_candidate"],
            how="inner",
        ).with_columns(
            pl.lit("internal_eaf").alias("__strand_resolution_basis")
        )
        selected_pal = (
            consensus_selected if consensus_known else frequency_selected
        )
    selected_columns = list(candidates.columns) + ["__strand_resolution_basis"]
    eligible = pl.concat([
        candidates.filter(~pl.col("__strand_palindromic")).with_columns(
            pl.lit("allele_match").alias("__strand_resolution_basis")
        ).select(selected_columns),
        selected_pal.select(selected_columns),
    ], how="vertical_relaxed").unique(
        subset=[row_col, STRAND_ACTION_COLUMN, "__strand_ref", "__strand_alt"],
        maintain_order=True,
    )
    eligible_counts = eligible.group_by(row_col).len().rename({"len": "__strand_eligible"})
    resolved = eligible.join(
        eligible_counts, on=row_col, how="left", coalesce=True
    ).filter(
        pl.col("__strand_eligible") == 1
    ).select([
        row_col,
        STRAND_ACTION_COLUMN,
        "__strand_swapped",
        pl.col("__strand_ref").alias("__strand_output_oa"),
        pl.col("__strand_alt").alias("__strand_output_ea"),
        pl.col("__strand_af").alias(REFERENCE_AF_COLUMN),
        "__strand_resolution_basis",
    ])
    pal_rows = (
        study.filter(palindromic_snp_expression(
            pl.col(ea_col), pl.col(oa_col)
        ))
        .select(row_col)
        .unique()
        .with_columns(pl.lit(True).alias("__strand_is_palindromic"))
    )
    tagged = (
        study.join(candidate_counts, on=row_col, how="left", coalesce=True)
        .join(pal_rows, on=row_col, how="left", coalesce=True)
        .join(frequency_conflicts, on=row_col, how="left", coalesce=True)
        .join(pal_frequency, on=row_col, how="left", coalesce=True)
        .join(resolved, on=row_col, how="left", coalesce=True)
        .with_columns([
            pl.col("__strand_candidates").fill_null(0),
            pl.col("__strand_is_palindromic").fill_null(False),
            pl.col("__strand_frequency_conflict").fill_null(False),
        ])
    )

    unmatched = pl.col("__strand_candidates") == 0
    unmatched_action = str(pol.get("strand.unmatched_action"))
    unmatched_palindromic_orientation_unavailable = (
        unmatched
        & pl.col("__strand_is_palindromic")
        & pl.lit(
            unmatched_action == "retain"
            and mode == "auto"
            and not consensus_known
        )
    ).fill_null(False)
    ambiguous_pal = (
        pl.col("__strand_is_palindromic")
        & pl.col(STRAND_ACTION_COLUMN).is_null()
        & ~unmatched
        & ~pl.col("__strand_frequency_conflict")
    )
    palindromic_orientation_unavailable = (
        (
            ambiguous_pal
            & pl.col("palindromic_frequency_state")
            .is_in(["not_available", "unusable_frequency"])
            .fill_null(True)
        )
        | unmatched_palindromic_orientation_unavailable
    )
    palindromic_frequency_discordant = ambiguous_pal & (
        pl.col("palindromic_frequency_state")
        == "maximum_difference_exceeded"
    ).fill_null(False)
    palindromic_ambiguous = (
        ambiguous_pal
        & ~palindromic_orientation_unavailable
        & ~palindromic_frequency_discordant
    )
    ambiguous_reference = (
        ~pl.col("__strand_is_palindromic")
        & pl.col(STRAND_ACTION_COLUMN).is_null()
        & ~unmatched
    )
    initial = tagged.height
    n_unmatched_detected = tagged.filter(unmatched).height
    if unmatched_action == "fail" and n_unmatched_detected:
        raise StrandOrientationError(
            "Chromosome %s contains reference-unmatched variants and strand.unmatched_action is 'fail'."
            % chromosome
        )
    if unmatched_action == "reject":
        tagged, n_unmatched = reject_rows(
            tagged, unmatched, step_label=STEP_LABEL,
            reason="reference_unmatched", context=ctx, collector=rejects,
            detail=(
                "no REF/ALT orientation matched the supplied chromosome "
                "reference"
            ),
        )
        n_unmatched_retained = 0
    elif unmatched_action == "retain":
        # The population-frequency panel could not orient these rows. Preserve
        # their supplied effect-allele order as explicit unverified provenance;
        # a strong reverse study consensus can still establish the strand of an
        # SNV, so apply only that complement here. The GWAS-to-VCF adapter later
        # requires the resulting REF allele (or its swap) to match the genome
        # FASTA, and record-count validation fails the chromosome if it skips a
        # retained row. This opt-in therefore relaxes panel coverage, not the
        # final reference-genome contract.
        retained_unmatched = (
            unmatched & ~unmatched_palindromic_orientation_unavailable
        ).fill_null(False)
        study_snv = (
            (pl.col(ea_col).str.len_chars() == 1)
            & (pl.col(oa_col).str.len_chars() == 1)
        ).fill_null(False)
        reverse_unmatched = (
            retained_unmatched
            & study_snv
            & pl.lit(consensus == "reverse")
        ).fill_null(False)
        tagged = tagged.with_columns([
            pl.when(retained_unmatched)
            .then(
                pl.when(reverse_unmatched)
                .then(reverse_complement_expression(pl.col(ea_col)))
                .otherwise(pl.col(ea_col))
            )
            .otherwise(pl.col("__strand_output_ea"))
            .alias("__strand_output_ea"),
            pl.when(retained_unmatched)
            .then(
                pl.when(reverse_unmatched)
                .then(reverse_complement_expression(pl.col(oa_col)))
                .otherwise(pl.col(oa_col))
            )
            .otherwise(pl.col("__strand_output_oa"))
            .alias("__strand_output_oa"),
            pl.when(retained_unmatched)
            .then(False)
            .otherwise(pl.col("__strand_swapped"))
            .alias("__strand_swapped"),
            pl.when(retained_unmatched)
            .then(
                pl.when(reverse_unmatched)
                .then(pl.lit(
                    "reference_unmatched_retained_reverse_complement"
                ))
                .otherwise(pl.lit("reference_unmatched_retained"))
            )
            .otherwise(pl.col(STRAND_ACTION_COLUMN))
            .alias(STRAND_ACTION_COLUMN),
            pl.when(retained_unmatched)
            .then(pl.lit("user_requested_reference_unmatched_retention"))
            .otherwise(pl.col("__strand_resolution_basis"))
            .alias("__strand_resolution_basis"),
        ])
        n_unmatched = 0
        n_unmatched_retained = tagged.filter(retained_unmatched).height
        if ctx is not None:
            ctx.qc(
                "reference-unmatched variants",
                "%s variants had no compatible allele row in the configured "
                "population-frequency reference. strand.unmatched_action is "
                "'retain', so they remain with no aligned panel AF and must "
                "still pass genome-FASTA validation during GWAS-to-VCF."
                % n_unmatched_retained,
                tagged.height,
                tagged.height,
                changed=n_unmatched_retained,
                warn=n_unmatched_retained > 0,
                step=STEP_LABEL,
            )
            if n_unmatched_retained:
                ctx.warn(
                    "Retained %s reference-unmatched variants by explicit "
                    "policy. They have no population-reference AF comparison; "
                    "GWAS-to-VCF will stop the chromosome if any supplied "
                    "allele pair cannot be validated against the genome FASTA."
                    % n_unmatched_retained
                )
    elif unmatched_action == "fail":
        # The non-zero case stopped above. Keeping the zero-count path explicit
        # avoids treating the valid configured policy itself as unsupported.
        n_unmatched = 0
        n_unmatched_retained = 0
    else:
        raise StrandOrientationError(
            "Unsupported strand.unmatched_action value %r." % unmatched_action
        )
    if str(pol.get("strand.ambiguous_action")) == "fail" and tagged.filter(
        ambiguous_pal
        | ambiguous_reference
        | unmatched_palindromic_orientation_unavailable
        | pl.col("__strand_frequency_conflict")
    ).height:
        raise StrandOrientationError(
            "Chromosome %s contains ambiguously oriented variants and strand.ambiguous_action is 'fail'."
            % chromosome
        )
    tagged, n_frequency_conflicts = reject_rows(
        tagged,
        pl.col("__strand_frequency_conflict"),
        step_label=STEP_LABEL,
        reason="palindromic_frequency_conflict",
        context=ctx,
        collector=rejects,
        detail=(
            "study-wide strand consensus and decisive internal effect-allele "
            "frequency evidence select different palindromic orientations"
        ),
    )
    tagged, n_orientation_unavailable = reject_rows(
        tagged,
        palindromic_orientation_unavailable,
        step_label=STEP_LABEL,
        reason="palindromic_orientation_unavailable",
        context=ctx,
        collector=rejects,
        detail=(
            "automatic mode had neither a strong study-wide strand consensus "
            "nor usable internal study EAF evidence for this palindrome"
        ),
    )
    tagged, n_frequency_discordant = reject_rows(
        tagged,
        palindromic_frequency_discordant,
        step_label=STEP_LABEL,
        reason="palindromic_frequency_discordant",
        context=ctx,
        collector=rejects,
        detail=(
            "both candidate orientations exceeded the configured maximum "
            "difference from the ancestry-matched reference AF"
        ),
    )
    tagged, n_palindromic = reject_rows(
        tagged,
        palindromic_ambiguous,
        step_label=STEP_LABEL,
        reason="palindromic_ambiguous",
        context=ctx, collector=rejects,
        detail=(
            "internal study EAF was near 0.5 or did not uniquely distinguish "
            "the competing reference orientations under the configured rules"
        ),
    )
    tagged, n_ambiguous = reject_rows(
        tagged, ambiguous_reference, step_label=STEP_LABEL, reason="reference_ambiguous",
        context=ctx, collector=rejects,
        detail="more than one reference orientation remained valid",
    )

    swap = pl.col("__strand_swapped").fill_null(False)
    updates = [
        pl.col("__strand_output_ea").alias(ea_col),
        pl.col("__strand_output_oa").alias(oa_col),
    ]
    if effect_col and effect_col in tagged.columns:
        effect = pl.col(effect_col).cast(pl.Float64, strict=False)
        transformed = pl.when(swap).then(-effect).otherwise(effect)
        updates.append(transformed.alias(effect_col))
    z_col = optional_text(sample_column_dict.get("imp_z_col"))
    if z_col and z_col in tagged.columns:
        z = pl.col(z_col).cast(pl.Float64, strict=False)
        updates.append(pl.when(swap).then(-z).otherwise(z).alias(z_col))
    if eaf_col and eaf_col in tagged.columns and not eaf_is_maf:
        eaf = pl.col(eaf_col).cast(pl.Float64, strict=False)
        updates.append(pl.when(swap).then(1.0 - eaf).otherwise(eaf).alias(eaf_col))

    output = tagged.with_columns(updates).sort(row_col)
    resolution_counts = {
        str(row["__strand_resolution_basis"]): int(row["len"])
        for row in output.filter(pl.col("__strand_is_palindromic"))
        .group_by("__strand_resolution_basis")
        .len()
        .to_dicts()
    }
    frequency_state_counts = {
        str(row["palindromic_frequency_state"]): int(row["len"])
        for row in pal_frequency.group_by("palindromic_frequency_state")
        .len()
        .to_dicts()
    }
    internal = [
        row_col, "__strand_candidates", "__strand_is_palindromic",
        "__strand_swapped", "__strand_output_oa", "__strand_output_ea",
        "__strand_frequency_conflict", "__strand_resolution_basis",
        "palindromic_frequency_state", "palindromic_frequency_candidate",
        "palindromic_frequency_value", "palindromic_frequency_error",
    ]
    output = output.drop([column for column in internal if column in output.columns])
    action_counts = {
        str(row[STRAND_ACTION_COLUMN]): int(row["len"])
        for row in output.group_by(STRAND_ACTION_COLUMN).len().to_dicts()
    }
    qc = {
        "status": "success",
        "initial_variants": initial,
        "final_variants": output.height,
        "reference_file": reference_file,
        "reference_population_column": population_col,
        "strand_mode": mode,
        "study_strand_consensus": consensus,
        "study_frequency_allele_alignment": frequency_alignment,
        "palindromic_orientation_basis": (
            "reference_aligned_declaration"
            if mode == "reference_aligned"
            else "study_wide_non_palindromic_consensus_or_internal_eaf"
        ),
        "palindromic_resolution_counts": resolution_counts,
        "palindromic_frequency_state_counts": frequency_state_counts,
        "actions": action_counts,
        "reference_unmatched_action": unmatched_action,
        "reference_unmatched_detected": n_unmatched_detected,
        "reference_unmatched": n_unmatched,
        "reference_unmatched_retained": n_unmatched_retained,
        "palindromic_frequency_conflict": n_frequency_conflicts,
        "palindromic_orientation_unavailable": n_orientation_unavailable,
        "palindromic_frequency_discordant": n_frequency_discordant,
        "palindromic_ambiguous": n_palindromic,
        "reference_ambiguous": n_ambiguous,
    }
    qc.update(reference_duplicate_stats)
    if ctx is not None:
        ctx.info(
            "Reference orientation retained %s of %s variants; actions: %s."
            % ("{:,}".format(output.height), "{:,}".format(initial), action_counts)
        )
        ctx.extra["strand_orientation"] = dict(qc)
    mark_canonical_variant_columns(sample_column_dict)
    return output, qc, sample_column_dict
