"""Reference-based allele orientation for one chromosome.

The supplied reference remains a normal CHROM/POS/REF/ALT/population-AF table.
One coordinate join is followed by vectorised tests of direct/swapped allele
order and, for SNVs, the two reverse-complement forms; no fourfold reference
table is materialised. Study-wide strand consensus is used only to arbitrate
palindromic SNPs, whose allele letters cannot identify the strand on their own.

Scientific method: NHGRI-EBI GWAS Catalog summary-statistics harmonisation,
https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics.
"""

from typing import Any, Dict, Optional, Tuple

import polars as pl

from postgwas.core.dataframes import chromosome_expression, position_expression
from postgwas.core.io.tables import read_delimited_table
from postgwas.core.values import optional_text

from .shared.allele_join import deduplicate_reference_values
from .shared.runtime import reject_rows, resolve_policies
from .shared.variant_columns import (
    has_canonical_variant_columns,
    mark_canonical_variant_columns,
)

__all__ = [
    "STRAND_ACTION_COLUMN",
    "REFERENCE_AF_COLUMN",
    "StrandOrientationError",
    "resolve_strand_consensus",
    "harmonise_strand_orientation",
]


STRAND_ACTION_COLUMN = "strand_action"
REFERENCE_AF_COLUMN = "strand_reference_af"
STEP_LABEL = "04 strand_orientation"
PALINDROMIC_PAIRS = ("AT", "TA", "CG", "GC")
POLICY_KEYS = (
    "strand.enabled",
    "strand.consensus_threshold",
    "strand.min_informative_variants",
    "strand.unmatched_action",
    "strand.ambiguous_action",
    "strand.af_tolerance",
    "strand.af_discordance_action",
    "strand.palindromic_af_discordance_action",
    "external_reference.exact_duplicate_action",
    "external_reference.non_identical_duplicate_action",
)


class StrandOrientationError(RuntimeError):
    """A chromosome cannot be oriented safely against its reference."""


def resolve_strand_consensus(genome_build_info, build, policies=None):
    """Resolve one study-wide forward/reverse decision from build evidence."""
    pol = resolve_policies(policies)
    if not bool(pol.get("strand.enabled")):
        return {
            "strand": "disabled",
            "forward": 0,
            "reverse": 0,
            "ambiguous": 0,
            "informative": 0,
            "dominant_fraction": None,
            "minimum_informative_variants": int(
                pol.get("strand.min_informative_variants")
            ),
            "consensus_threshold": float(pol.get("strand.consensus_threshold")),
            "reason": "strand.enabled is false",
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
        "forward": forward,
        "reverse": reverse,
        "ambiguous": ambiguous,
        "informative": informative,
        "dominant_fraction": dominant_fraction,
        "minimum_informative_variants": minimum,
        "consensus_threshold": threshold,
        "reason": reason,
    }


def _reverse_complement(column: str) -> pl.Expr:
    return (
        pl.col(column)
        .cast(pl.Utf8)
        .str.to_uppercase()
        .str.replace_many(["A", "C", "G", "T"], ["T", "G", "C", "A"])
        .str.reverse()
    )


def _load_reference(path, population_col, colmap, policies, chromosome):
    raw, _detected = read_delimited_table(
        path,
        colmap["delimiter"],
        candidates=list(policies.get("input.delimiter_candidates")),
        minimum_columns=int(policies.get("input.delimiter_min_columns")),
        maximum_columns=int(policies.get("input.delimiter_max_columns")),
        sample_lines=int(policies.get("input.delimiter_sample_rows")),
        null_values=list(policies.get("input.null_values")),
        infer_schema_length=int(policies.get("input.schema_inference_rows")),
        error_type=StrandOrientationError,
        description="chromosome %s strand reference" % chromosome,
    )
    required = [colmap[key] for key in ("chr", "pos", "a1", "a2")] + [population_col]
    missing = [column for column in required if column not in raw.columns]
    if missing:
        raise StrandOrientationError(
            "The strand reference '%s' is missing configured columns: %s"
            % (path, ", ".join(missing))
        )
    selected = raw.select(required).rename({
        colmap["chr"]: "__strand_chr",
        colmap["pos"]: "__strand_pos",
        colmap["a1"]: "__strand_alt",
        colmap["a2"]: "__strand_ref",
        population_col: "__strand_af",
    })
    schema = dict(selected.schema)
    reference = selected.with_columns([
        chromosome_expression(
            "__strand_chr", schema["__strand_chr"],
            strip_chr_prefix=bool(policies.get("chromosome.strip_chr_prefix")),
            strip_leading_zero=True,
        ),
        position_expression("__strand_pos", schema["__strand_pos"]),
        pl.col("__strand_ref").cast(pl.Utf8).str.to_uppercase(),
        pl.col("__strand_alt").cast(pl.Utf8).str.to_uppercase(),
        pl.col("__strand_af").cast(pl.Float64, strict=False),
    ])
    invalid = reference.filter(
        pl.col("__strand_af").is_not_null()
        & ((pl.col("__strand_af") < 0.0) | (pl.col("__strand_af") > 1.0))
    ).height
    if invalid:
        raise StrandOrientationError(
            "The strand reference '%s' contains %s values outside 0 to 1 in column '%s'."
            % (path, "{:,}".format(invalid), population_col)
        )
    return deduplicate_reference_values(
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


def _study_value(decision, key, default=None):
    if key == "eaf_is_maf" and isinstance(decision, bool):
        return decision
    if isinstance(decision, dict):
        return decision.get(key, default)
    return default


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
    """Orient one chromosome with a single raw-reference coordinate join."""
    pol = resolve_policies(policies)
    if not bool(pol.get("strand.enabled")):
        output = df.with_columns(
            pl.lit("disabled").alias(STRAND_ACTION_COLUMN)
        )
        return output, {
            "status": "disabled",
            "initial_variants": df.height,
            "final_variants": output.height,
            "actions": {"disabled": output.height},
        }, sample_column_dict
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
            "safely choose between negating beta and reciprocating an odds ratio."
        )
    se_col = optional_text(sample_column_dict.get("se_col"))
    if effect_type == "odds_ratio" and se_col and se_col in df.columns:
        se_scale = _study_value(study_decision, "se_scale", pol.get("effect.se_scale"))
        if se_scale != "log_odds":
            raise StrandOrientationError(
                "Odds ratios with a raw or unresolved standard-error scale cannot be "
                "oriented safely. Call harmonise_effect_estimates() first so OR and SE "
                "are both converted to log-odds units; strand orientation can then "
                "negate beta while leaving SE unchanged."
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
        reference_file, population_col, reference_colmap, pol, chromosome,
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
        study_schema = dict(df.schema)
        study = study.with_columns([
            chromosome_expression(
                chr_col,
                study_schema[chr_col],
                strip_chr_prefix=bool(pol.get("chromosome.strip_chr_prefix")),
                strip_leading_zero=True,
            ),
            position_expression(pos_col, study_schema[pos_col]),
            pl.col(ea_col).cast(pl.Utf8).str.strip_chars().str.to_uppercase(),
            pl.col(oa_col).cast(pl.Utf8).str.strip_chars().str.to_uppercase(),
        ])
    joined = study.join(
        reference,
        left_on=[chr_col, pos_col],
        right_on=["__strand_chr", "__strand_pos"],
        how="inner",
    ).with_columns([
        _reverse_complement(ea_col).alias("__strand_ea_rc"),
        _reverse_complement(oa_col).alias("__strand_oa_rc"),
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

    pair = pl.concat_str([pl.col(ea_col), pl.col(oa_col)])
    candidates = candidates.with_columns(
        pair.is_in(list(PALINDROMIC_PAIRS)).alias("__strand_palindromic")
    )
    consensus = str(_study_value(study_decision, "strand", "unresolved") or "unresolved")

    candidate_strand = (
        pl.when(pl.col(STRAND_ACTION_COLUMN).str.starts_with("forward"))
        .then(pl.lit("forward"))
        .otherwise(pl.lit("reverse"))
    )
    candidates = candidates.with_columns(candidate_strand.alias("__strand_basis"))

    consensus_known = consensus in ("forward", "reverse")
    pal_pool = (
        candidates.filter(
            pl.col("__strand_palindromic")
            & (pl.col("__strand_basis") == pl.lit(consensus))
        )
        if consensus_known
        else candidates.head(0)
    )

    eligible = pl.concat([
        candidates.filter(~pl.col("__strand_palindromic")),
        pal_pool,
    ], how="vertical_relaxed").unique(
        subset=[row_col, STRAND_ACTION_COLUMN, "__strand_ref", "__strand_alt"],
        maintain_order=True,
    )
    eligible_counts = eligible.group_by(row_col).len().rename({"len": "__strand_eligible"})
    resolved = eligible.join(eligible_counts, on=row_col, how="left").filter(
        pl.col("__strand_eligible") == 1
    ).select([
        row_col,
        STRAND_ACTION_COLUMN,
        "__strand_swapped",
        pl.col("__strand_ref").alias("__strand_output_oa"),
        pl.col("__strand_alt").alias("__strand_output_ea"),
        pl.col("__strand_af").alias(REFERENCE_AF_COLUMN),
    ])
    candidate_counts = candidates.group_by(row_col).len().rename({"len": "__strand_candidates"})
    pal_rows = candidates.filter(pl.col("__strand_palindromic")).select(row_col).unique().with_columns(
        pl.lit(True).alias("__strand_is_palindromic")
    )
    tagged = (
        study.join(candidate_counts, on=row_col, how="left")
        .join(pal_rows, on=row_col, how="left")
        .join(resolved, on=row_col, how="left")
        .with_columns([
            pl.col("__strand_candidates").fill_null(0),
            pl.col("__strand_is_palindromic").fill_null(False),
        ])
    )

    unmatched = pl.col("__strand_candidates") == 0
    ambiguous_pal = (
        pl.col("__strand_is_palindromic")
        & pl.col(STRAND_ACTION_COLUMN).is_null()
        & ~unmatched
    )
    ambiguous_reference = (
        ~pl.col("__strand_is_palindromic")
        & pl.col(STRAND_ACTION_COLUMN).is_null()
        & ~unmatched
    )
    initial = tagged.height
    if str(pol.get("strand.unmatched_action")) == "fail" and tagged.filter(unmatched).height:
        raise StrandOrientationError(
            "Chromosome %s contains reference-unmatched variants and strand.unmatched_action is 'fail'."
            % chromosome
        )
    tagged, n_unmatched = reject_rows(
        tagged, unmatched, step_label=STEP_LABEL, reason="reference_unmatched",
        context=ctx, collector=rejects,
        detail="no REF/ALT orientation matched the supplied chromosome reference",
    )
    if str(pol.get("strand.ambiguous_action")) == "fail" and tagged.filter(
        ambiguous_pal | ambiguous_reference
    ).height:
        raise StrandOrientationError(
            "Chromosome %s contains ambiguously oriented variants and strand.ambiguous_action is 'fail'."
            % chromosome
        )
    tagged, n_palindromic = reject_rows(
        tagged, ambiguous_pal, step_label=STEP_LABEL, reason="palindromic_ambiguous",
        context=ctx, collector=rejects,
        detail=(
            "study-wide non-palindromic strand consensus was mixed, unresolved, "
            "or did not leave exactly one safe reference orientation"
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
        transformed = (
            pl.when(swap & (effect > 0.0)).then(1.0 / effect).otherwise(effect)
            if effect_type == "odds_ratio"
            else pl.when(swap).then(-effect).otherwise(effect)
        )
        updates.append(transformed.alias(effect_col))
    z_col = optional_text(sample_column_dict.get("imp_z_col"))
    if z_col and z_col in tagged.columns:
        z = pl.col(z_col).cast(pl.Float64, strict=False)
        updates.append(pl.when(swap).then(-z).otherwise(z).alias(z_col))
    if eaf_col and eaf_col in tagged.columns and not eaf_is_maf:
        eaf = pl.col(eaf_col).cast(pl.Float64, strict=False)
        updates.append(pl.when(swap).then(1.0 - eaf).otherwise(eaf).alias(eaf_col))

    output = tagged.with_columns(updates).sort(row_col)
    internal = [
        row_col, "__strand_candidates", "__strand_is_palindromic",
        "__strand_swapped", "__strand_output_oa", "__strand_output_ea",
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
        "study_strand_consensus": consensus,
        "palindromic_orientation_basis": "study_wide_non_palindromic_consensus",
        "actions": action_counts,
        "reference_unmatched": n_unmatched,
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
