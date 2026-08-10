"""Step 02 - chromosome and position harmonisation.

Splits a combined chromosome:position column when one is configured, normalises
the chromosome label, turns the position into an integer, and removes variants
whose coordinates cannot be used.

Everything this step does is driven by the policy registry
(``harmonisation/policies.py``) and recorded through the pipeline logger
(``core/pipeline_logging.py``). Removed variants go through the reject
collector (``harmonisation/rejects.py``) so that every dropped row is written to
the reject file with a reason, and rows in minus rejected equals rows out.

Policies read here
------------------
position.extraction              how the position is pulled out of a combined field
position.min_value               smallest position accepted
chromosome.strip_chr_prefix      treat 'chr7' and '7' as the same
chromosome.rename_map            how non-standard labels are renamed
chromosome.allowed_after_split   which chromosomes survive this step
chromosome.drop_mt               drop mitochondrial variants
"""

from typing import Any, Dict, Optional, Tuple

import polars as pl

from postgwas.core.dataframes import chromosome_expression
from postgwas.core.values import optional_text

from .policies import default_policies
from .shared.runtime import reject_rows, step_context


__all__ = ["harmonise_coordinates_and_alleles"]


#: The label written into the reject file's ``reject_step`` column.
STEP_LABEL = "02 fix_chr_pos_allele"

#: Resolved values recorded as PARAM entries at the top of the step.
POLICY_KEYS = [
    "position.extraction",
    "position.min_value",
    "chromosome.strip_chr_prefix",
    "chromosome.rename_map",
    "chromosome.allowed_after_split",
    "chromosome.drop_mt",
]


# =============================================================================
# small helpers
# =============================================================================

_clean_name = optional_text


def _resolve_drop_mt(drop_mt, policies, ctx):
    # type: (Optional[bool], Any, Any) -> bool
    """The argument still wins, unless the configuration set the policy itself.

    ``drop_mt`` used to default to True; it now defaults to None, meaning "use
    chromosome.drop_mt", whose own default is True - so nothing changes for a
    caller that says nothing.  A caller that passes the flag explicitly keeps
    winning, except when the run's configuration set chromosome.drop_mt, which
    is the only way a user can express an opinion.
    """
    policy_value = bool(policies.chromosome.drop_mt)
    try:
        policy_is_default = policies.is_default("chromosome.drop_mt")
    except Exception:
        policy_is_default = True

    if not policy_is_default:
        if drop_mt is not None and bool(drop_mt) != policy_value:
            ctx.info(
                "The caller asked for drop_mt=%s, but chromosome.drop_mt was set to %s in the "
                "configuration, so the configuration wins."
                % (bool(drop_mt), policy_value)
            )
        return policy_value
    if drop_mt is None:
        return policy_value
    return bool(drop_mt)


def _allowed_chromosomes(policies, drop_mt):
    # type: (Any, bool) -> list
    """The chromosomes kept after this step, as upper-case strings."""
    allowed = [str(c).strip().upper() for c in (policies.chromosome.allowed_after_split or [])]
    if drop_mt:
        allowed = [c for c in allowed if c != "MT"]
    elif "MT" not in allowed:
        allowed.append("MT")
    return allowed


# =============================================================================
# the step
# =============================================================================

def harmonise_coordinates_and_alleles(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    drop_mt: Optional[bool] = None,
    logger: Any = None,
    policies: Any = None,
    rejects: Any = None,
    ctx: Any = None,
    step_number: int = 2,
    step_total: int = 8,
    return_qc_info: bool = False,
):
    """Harmonise the chromosome and position columns.

    Returns ``(df, sample_column_dict)``, or ``(df, sample_column_dict, qc_info)``
    when ``return_qc_info=True`` - the extra return value is opt-in so that
    existing two-value call sites keep working unchanged.  The same dictionary is
    always attached to the step context, so it reaches ``logger.summary()``
    whether or not the caller asks for it.

    Parameters
    ----------
    chromosome : the chromosome being processed, used as the log scope.
    df, sample_column_dict : as before.
    drop_mt : True/False as before; None (the new default) means "use the policy
        chromosome.drop_mt", whose default is True - identical behaviour.
    logger : the pipeline logger; direct library calls may omit it.
    policies : a Policies object; the defaults are used when it is None.
    rejects : a RejectCollector.  Without one the variants are still removed and
        still counted, they are just not written to the reject file.
    ctx : an open StepContext, when the caller has already opened the step.
    """
    policies = default_policies() if policies is None else policies

    rows_in = df.height
    counters = {}       # type: Dict[str, int]
    qc_info = {
        "step": STEP_LABEL,
        "chromosome": str(chromosome),
        "initial_variants": rows_in,
    }               # type: Dict[str, Any]

    with step_context(
        logger, ctx,
        number=step_number, total=step_total,
        title="Chromosome and position",
        operation="coordinates.harmonise_coordinates_and_alleles",
        rows_in=rows_in, policy_keys=POLICY_KEYS,
    ) as step_ctx:
        df, sample_column_dict = _harmonise(
            df=df,
            sample_column_dict=sample_column_dict,
            policies=policies,
            ctx=step_ctx,
            rejects=rejects,
            drop_mt=drop_mt,
            counters=counters,
            qc_info=qc_info,
        )
        qc_info["removed_by_reason"] = dict(counters)
        qc_info["final_variants"] = df.height
        qc_info["removed_total"] = rows_in - df.height
        step_ctx.set_rows(df.height)
        step_ctx.extra.update(qc_info)

    if return_qc_info:
        return df, sample_column_dict, qc_info
    return df, sample_column_dict


def _harmonise(df, sample_column_dict, policies, ctx, rejects,
               drop_mt, counters, qc_info):
    # type: (pl.DataFrame, dict, Any, Any, Any, Optional[bool], Dict[str, int], Dict[str, Any]) -> Tuple[pl.DataFrame, dict]
    chr_col = _clean_name(sample_column_dict.get("chr_col"))
    pos_col = _clean_name(sample_column_dict.get("pos_col"))
    chr_pos_col = _clean_name(sample_column_dict.get("chr_pos_col"))

    effective_drop_mt = _resolve_drop_mt(drop_mt, policies, ctx)
    allowed = _allowed_chromosomes(policies, effective_drop_mt)
    qc_info["drop_mt"] = effective_drop_mt
    qc_info["allowed_chromosomes"] = list(allowed)

    # ------------------------------------------------------------------
    # 1. split a combined chromosome:position column
    # ------------------------------------------------------------------
    qc_info["split_performed"] = False
    if chr_pos_col and chr_pos_col in df.columns and (not chr_col or not pos_col):
        df, chr_col, pos_col = _split_chromosome_position(
            df, chr_pos_col, policies, ctx, qc_info
        )
        sample_column_dict.update({
            "chr_col": "Chr",
            "pos_col": "Pos",
            "chr_pos_col": None,
        })

    if not chr_col or chr_col not in df.columns:
        raise ValueError(
            "No chromosome column to work with. The configuration points at %r and the file "
            "has %s. Set chr_col, or set chr_pos_col to a combined chromosome:position column."
            % (chr_col, ", ".join(df.columns[:20]) or "no columns")
        )

    # ------------------------------------------------------------------
    # 2. normalise the chromosome label and the position
    # ------------------------------------------------------------------
    schema = dict(df.schema)
    updates = [
        chromosome_expression(
            chr_col,
            schema[chr_col],
            strip_chr_prefix=bool(policies.chromosome.strip_chr_prefix),
            strip_leading_zero=True,
        ).str.replace(r"\.0$", "", literal=False).str.to_uppercase(),
    ]

    has_pos = bool(pos_col) and pos_col in df.columns
    if has_pos:
        unreadable = df.select(
            (
                pl.col(pos_col).is_not_null()
                & pl.col(pos_col).cast(pl.Float64, strict=False).is_null()
            ).sum()
        ).item()
        qc_info["positions_unreadable"] = int(unreadable or 0)
        if unreadable:
            ctx.warn(
                "%d positions are present but could not be read as a number, so they become "
                "empty and are removed below as invalid positions." % int(unreadable)
            )
        # Decimal-looking integer coordinates such as ``100.0`` are common in
        # exported tables and were historically accepted by this input step.
        updates.append(
            pl.col(pos_col)
            .cast(pl.Float64, strict=False)
            .cast(pl.Int64, strict=False)
            .alias(pos_col)
        )
    else:
        qc_info["positions_unreadable"] = 0
        ctx.warn(
            "No position column is present, so positions cannot be checked at this step."
        )

    for allele_key in ("ea_col", "oa_col"):
        allele_col = _clean_name(sample_column_dict.get(allele_key))
        if allele_col and allele_col in df.columns:
            updates.append(
                pl.col(allele_col)
                .cast(pl.String, strict=False)
                .str.strip_chars()
                .replace("", None)
                .str.to_uppercase()
                .alias(allele_col)
            )
    df = df.with_columns(updates)

    # ------------------------------------------------------------------
    # 3. rename non-standard chromosome labels
    # ------------------------------------------------------------------
    df = _rename_chromosomes(df, chr_col, policies, ctx, effective_drop_mt, qc_info)

    # ------------------------------------------------------------------
    # 4. remove the variants that cannot be used
    # ------------------------------------------------------------------
    df, _ = reject_rows(
        df, pl.col(chr_col).is_null(), step_label=STEP_LABEL,
        reason="invalid_chromosome", context=ctx, collector=rejects,
        counters=counters, detail="chromosome is empty",
    )

    if has_pos:
        df, _ = reject_rows(
            df, pl.col(pos_col).is_null(), step_label=STEP_LABEL,
            reason="invalid_position", context=ctx, collector=rejects,
            counters=counters,
            detail="position is empty or could not be read as a number",
        )
        min_value = policies.position.min_value
        if min_value is not None:
            df, _ = reject_rows(
                df,
                pl.col(pos_col).is_not_null() & (pl.col(pos_col) < int(min_value)),
                step_label=STEP_LABEL, reason="invalid_position", context=ctx,
                collector=rejects, counters=counters,
                detail="position below position.min_value (%d)" % int(min_value),
            )
        qc_info["position_min_value"] = None if min_value is None else int(min_value)

    _warn_about_allowed_sets(df, chr_col, policies, allowed, ctx)
    df, _ = reject_rows(
        df, ~pl.col(chr_col).is_in(allowed), step_label=STEP_LABEL,
        reason="invalid_chromosome", context=ctx, collector=rejects,
        counters=counters,
        detail="not one of the chromosomes being kept (%s)" % ", ".join(allowed),
    )

    return df, sample_column_dict


# =============================================================================
# the pieces
# =============================================================================

def _split_chromosome_position(df, chr_pos_col, policies, ctx, qc_info):
    # type: (pl.DataFrame, str, Any, Any, Dict[str, Any]) -> Tuple[pl.DataFrame, str, str]
    """Split '1:123456' (or '1:123456:A:G') into 'Chr' and 'Pos'."""
    ctx.info(
        "The chromosome and position columns are not both configured, so the combined column "
        "'%s' is being split into 'Chr' and 'Pos'." % chr_pos_col
    )

    non_null_before = df.select(pl.col(chr_pos_col).is_not_null().sum()).item()
    dup_before = df.select(pl.col(chr_pos_col).is_duplicated().sum()).item()

    used_sep = None
    for sep in [":", "-", "_"]:
        try:
            temp = df.with_columns([
                pl.col(chr_pos_col).cast(pl.Utf8, strict=False)
                  .str.split_exact(sep, 1).struct.field("field_0").alias("Chr"),
                pl.col(chr_pos_col).cast(pl.Utf8, strict=False)
                  .str.split_exact(sep, 1).struct.field("field_1").alias("Pos"),
            ])
            non_null_after = temp.select(pl.col("Pos").is_not_null().sum()).item()
            if non_null_after == non_null_before:
                df = temp
                used_sep = sep
                break
        except Exception as exc:                       # pragma: no cover - defensive
            ctx.info("Separator '%s' could not be used: %s" % (sep, exc))

    if used_sep is None:
        raise RuntimeError(
            "Could not split '%s' with ':', '-' or '_'. %s values are present, and none of the "
            "three separators produced a position for all of them. Check the column really does "
            "hold a combined chromosome and position."
            % (chr_pos_col, "{:,}".format(int(non_null_before)))
        )

    ctx.decide(
        "Combined chromosome/position column '%s'" % chr_pos_col,
        {"separator that worked": used_sep, "values split": int(non_null_before)},
        "split into 'Chr' and 'Pos'",
    )

    # A non-null combined value must produce a non-null chromosome.  This is the
    # honest version of the old "duplicate mismatch" guard, which compared
    # duplicates of the full identifier against duplicates of the (Chr, Pos)
    # pair - two alleles at one position are one coordinate but two identifiers,
    # so those two counts differ for perfectly good data and the guard raised on
    # the most common identifier convention.
    lost = df.select(
        (pl.col(chr_pos_col).is_not_null() & pl.col("Chr").is_null()).sum()
    ).item()
    if lost:
        raise RuntimeError(
            "%s rows have a value in '%s' but produced no chromosome after splitting on '%s'."
            % ("{:,}".format(int(lost)), chr_pos_col, used_sep)
        )

    df, changed, empty_after = _extract_position(df, "Pos", chr_pos_col, policies, ctx)

    dup_after = df.select(pl.struct(["Chr", "Pos"]).is_duplicated().sum()).item()
    ctx.info(
        "After the split, %s rows share an identical '%s' value and %s rows share an identical "
        "chromosome and position. These two numbers are not expected to match when the combined "
        "column also carries the alleles (for example '1:100:A:G'): two alleles at one position "
        "are one coordinate but two identifiers."
        % ("{:,}".format(int(dup_before)), chr_pos_col, "{:,}".format(int(dup_after)))
    )

    qc_info.update({
        "split_performed": True,
        "split_separator": used_sep,
        "split_values": int(non_null_before),
        "split_duplicate_ids": int(dup_before),
        "split_duplicate_coordinates": int(dup_after),
        "positions_rewritten": int(changed),
        "positions_empty_after_extraction": int(empty_after),
    })
    return df, "Chr", "Pos"


def _extract_position(df, pos_col, source_col, policies, ctx):
    # type: (pl.DataFrame, str, str, Any, Any) -> Tuple[pl.DataFrame, int, int]
    """Pull the position out of the split field, the way the policy asks.

    'strip_non_digits' deletes every
    non-digit and glues the rest together, so '1:123456.0' becomes 1234560 - ten
    times too large - and '1:123456-123457' becomes 123456123457.  Because the
    later cast is non-strict, those never become null; they become plausible
    wrong integers.  'leading_digits' is the safe default and takes only the
    digits at the front. 'none' leaves the text alone so the numeric cast decides.
    """
    extraction = str(policies.position.extraction)

    changed = df.select(
        (pl.col(pos_col).is_not_null() & pl.col(pos_col).str.contains(r"\D")).sum()
    ).item() or 0

    expr = position_text_expression(pl.col(pos_col), extraction)
    if extraction == "none":
        changed = 0

    df = df.with_columns(expr.alias(pos_col))

    empty_after = df.select(
        (
            pl.col(source_col).is_not_null()
            & (pl.col(pos_col).is_null() | (pl.col(pos_col) == ""))
        ).sum()
    ).item() or 0

    plain = (
        "Positions were taken out of '%s' using position.extraction = '%s'. %s of them contained "
        "characters other than digits."
        % (source_col, extraction, "{:,}".format(int(changed)))
    )
    if extraction == "strip_non_digits" and changed:
        plain += (
            " Every non-digit was deleted and the remaining digits glued together, so a value "
            "like '123456.0' becomes 1234560 - ten times too large - and '123456-123457' becomes "
            "123456123457. 'leading_digits' is the recommended setting."
        )
    if empty_after:
        plain += (
            " %s values produced no position at all and are removed below."
            % "{:,}".format(int(empty_after))
        )

    before = df.height
    ctx.qc(
        "position extraction", plain, before, before,
        changed=int(changed),
        warn=bool(empty_after) or (extraction == "strip_non_digits" and bool(changed)),
    )
    return df, int(changed), int(empty_after)


def position_text_expression(expression, extraction):
    """Apply the configured position-text policy to a Polars expression."""
    if extraction == "leading_digits":
        return expression.str.extract(r"^(\d+)", 1)
    if extraction == "none":
        return expression
    # Explicit ``strip_non_digits`` compatibility choice.
    return expression.str.replace_all(r"\D", "")


def _rename_chromosomes(df, chr_col, policies, ctx, drop_mt, qc_info):
    # type: (pl.DataFrame, str, Any, Any, bool, Dict[str, Any]) -> pl.DataFrame
    """Rename non-standard chromosome labels, and say what that costs.

    The default map is the one the code has always used, 23->X, 24->Y, 25->MT,
    M->MT.  In PLINK, 25 is the pseudoautosomal region (XY) and 26 is
    mitochondrial, so with the default map every PAR variant is relabelled MT
    and then deleted by chromosome.drop_mt, while true mitochondrial variants
    labelled 26 are not recognised at all.  The recommended map is
    23->X, 24->Y, 25->XY, 26->MT, M->MT.
    """
    mapping = dict(policies.chromosome.rename_map or {})
    mapping = {str(k).strip().upper(): str(v).strip().upper() for k, v in mapping.items()}
    qc_info["rename_map"] = dict(mapping)

    if not mapping:
        qc_info["chromosomes_renamed"] = 0
        return df

    keys = list(mapping.keys())
    counted = df.select([
        pl.col(chr_col).is_in(keys).sum().alias("renamed"),
        (pl.col(chr_col) == "25").sum().alias("plink_25"),
        (pl.col(chr_col) == "26").sum().alias("plink_26"),
    ]).to_dicts()[0]
    renamed = int(counted["renamed"] or 0)
    n_25 = int(counted["plink_25"] or 0)
    n_26 = int(counted["plink_26"] or 0)

    before = df.height
    df = df.with_columns(pl.col(chr_col).replace(mapping).alias(chr_col))

    ctx.qc(
        "chromosome renaming",
        "Non-standard chromosome labels were renamed using chromosome.rename_map (%s)."
        % ", ".join("%s->%s" % (k, v) for k, v in mapping.items()),
        before, before, changed=renamed,
    )

    if n_25 and mapping.get("25") == "MT" and drop_mt:
        ctx.warn(
            "%s variants are on chromosome 25. In PLINK, 25 is the pseudoautosomal region (XY), "
            "not mitochondrial, but chromosome.rename_map sends 25 to MT and chromosome.drop_mt "
            "is on, so all %s of them are about to be deleted. Set chromosome.rename_map to "
            "23->X, 24->Y, 25->XY, 26->MT, M->MT to keep them."
            % ("{:,}".format(n_25), "{:,}".format(n_25))
        )
    if n_26 and "26" not in mapping:
        ctx.warn(
            "%s variants are on chromosome 26, which is mitochondrial in PLINK, but "
            "chromosome.rename_map has no entry for 26, so they are not recognised as MT and are "
            "about to be removed as an invalid chromosome. Set chromosome.rename_map to "
            "23->X, 24->Y, 25->XY, 26->MT, M->MT." % "{:,}".format(n_26)
        )

    qc_info["chromosomes_renamed"] = renamed
    qc_info["plink_25_variants"] = n_25
    qc_info["plink_26_variants"] = n_26
    return df


def _warn_about_allowed_sets(df, chr_col, policies, allowed, ctx):
    # type: (pl.DataFrame, str, Any, list, Any) -> None
    """Say so when a chromosome that survived the read is dropped here.

    chromosome.allowed (used while reading, io.py) keeps 1-22, X, Y and MT,
    while chromosome.allowed_after_split keeps only 1-22, X and Y.  The two
    settings disagree, so mitochondrial variants survive the read and are then
    silently dropped at this step.  The two should be made the same.
    """
    at_read = [str(c).strip().upper() for c in (policies.chromosome.allowed or [])]
    only_at_read = [c for c in at_read if c not in allowed]
    if not only_at_read:
        return
    affected = df.select(pl.col(chr_col).is_in(only_at_read).sum()).item() or 0
    if not affected:
        return
    ctx.warn(
        "%s variants are on %s, which chromosome.allowed keeps while the input file is read but "
        "chromosome.allowed_after_split does not keep here. They are removed below and recorded "
        "as an invalid chromosome. These two settings disagree by default and should be made the "
        "same." % ("{:,}".format(int(affected)), ", ".join(only_at_read))
    )
