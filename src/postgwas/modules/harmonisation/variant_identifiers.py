"""Step 12 - the SNP identifier column.

Makes sure every variant has an identifier: creates one from chromosome,
position and the two alleles when the study has none, replaces ':' with '_' so
the identifier survives the VCF round trip, and fills in the identifiers that
are missing.

This step never removes a variant, so it takes a reject collector only for
symmetry with the other steps and does not use it.

Policies read here
------------------
input.null_values   configured reader tokens added to shared identifier sentinels
"""

from typing import Any, Dict, List, Optional, Tuple

import polars as pl

from postgwas.core.values import missing_tokens, optional_text

from .policies import default_policies
from .shared.runtime import step_context


__all__ = ["harmonise_variant_identifiers"]


STEP_LABEL = "12 snp_column"

POLICY_KEYS = ["input.null_values"]


# =============================================================================
# small helpers
# =============================================================================

_clean_name = optional_text


def _missing_sentinels(policies):
    # type: (Any) -> List[str]
    """The upper-case strings that mean "this identifier is missing".

    Identifier matching combines the configured input.null_values with the
    shared optional-text placeholders. This deliberately recognizes '-' as a
    missing identifier without asking the table reader to erase '-' from allele
    columns, where older indel representations may use it for a deletion.
    """
    return list(missing_tokens(policies.input.null_values))


def _missing_expr(snp_col, sentinels):
    # type: (str, List[str]) -> pl.Expr
    return (
        pl.col(snp_col).is_null()
        | pl.col(snp_col).cast(pl.Utf8, strict=False).str.strip_chars()
                         .str.to_uppercase().is_in(sentinels)
    )


def _identifier_expr(chr_col, pos_col, ea_col, oa_col):
    # type: (str, str, str, str) -> pl.Expr
    """chr_pos_ea_oa, built from whatever those columns hold."""
    return (
        pl.col(chr_col).cast(pl.Utf8, strict=False)
        + pl.lit("_")
        + pl.col(pos_col).cast(pl.Utf8, strict=False)
        + pl.lit("_")
        + pl.col(ea_col).cast(pl.Utf8, strict=False)
        + pl.lit("_")
        + pl.col(oa_col).cast(pl.Utf8, strict=False)
    )


def _coordinate_columns(df, sample_column_dict):
    # type: (pl.DataFrame, dict) -> Tuple[List[str], List[str]]
    """The four columns an identifier is built from, and the ones that are missing."""
    names = []
    missing = []
    for key in ("chr_col", "pos_col", "ea_col", "oa_col"):
        name = _clean_name(sample_column_dict.get(key))
        if name is None or name not in df.columns:
            missing.append("%s (%s)" % (key, name if name else "not configured"))
        names.append(name)
    return names, missing


# =============================================================================
# the step
# =============================================================================

def harmonise_variant_identifiers(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    logger: Any = None,
    policies: Any = None,
    rejects: Any = None,
    ctx: Any = None,
    step_number: int = 12,
    step_total: int = 16,
) -> Tuple[pl.DataFrame, dict]:
    """Ensure the dataframe has a usable variant identifier column.

    - creates one from chromosome, position and the two alleles when the study
      has none, or names one that is not in the file;
    - replaces ':' with '_';
    - fills in identifiers that are missing, empty, or one of the configured or
      shared optional-text placeholders.

    Returns ``(df, sample_column_dict)``.

    Parameters
    ----------
    logger : the pipeline logger; direct library calls may omit it.
    policies : a Policies object; the defaults are used when it is None.
    rejects : accepted for symmetry with the other steps; this step removes no
        variants, so it is not used.
    ctx : an open StepContext, when the caller has already opened the step.
    """
    policies = default_policies() if policies is None else policies

    rows_in = df.height
    with step_context(
        logger, ctx,
        number=step_number, total=step_total,
        title="SNP identifier",
        operation="variant_identifiers.harmonise_variant_identifiers",
        rows_in=rows_in, policy_keys=POLICY_KEYS,
    ) as step_ctx:
        df, sample_column_dict, info = _harmonise(
            df=df,
            sample_column_dict=sample_column_dict,
            policies=policies,
            ctx=step_ctx,
        )
        step_ctx.set_rows(df.height)
        step_ctx.extra.update(info)

    return df, sample_column_dict


def _harmonise(df, sample_column_dict, policies, ctx):
    # type: (pl.DataFrame, dict, Any, Any) -> Tuple[pl.DataFrame, dict, Dict[str, Any]]
    info = {"step": STEP_LABEL}     # type: Dict[str, Any]

    snp_col = _clean_name(sample_column_dict.get("snp_id_col"))
    (chr_col, pos_col, ea_col, oa_col), missing_coords = _coordinate_columns(
        df, sample_column_dict
    )
    can_build = not missing_coords

    # ------------------------------------------------------------------
    # 1. create the column when there is none
    # ------------------------------------------------------------------
    if snp_col is None or snp_col not in df.columns:
        if not can_build:
            raise ValueError(
                "There is no SNP identifier column (snp_id_col = %r) and one cannot be built, "
                "because these columns are not in the data: %s. Set snp_id_col, or make sure the "
                "chromosome, position and allele columns are configured."
                % (sample_column_dict.get("snp_id_col"), "; ".join(missing_coords))
            )
        ctx.info(
            "No SNP identifier column is present, so one is being built from %s, %s, %s and %s "
            "in the form chromosome_position_effect_other."
            % (chr_col, pos_col, ea_col, oa_col)
        )
        snp_col = "snp_col"
        df = df.with_columns(
            _identifier_expr(chr_col, pos_col, ea_col, oa_col).alias(snp_col)
        )
        sample_column_dict["snp_id_col"] = snp_col
        info["source"] = "built_from_coordinates"
    else:
        ctx.info(
            "The study supplies the SNP identifier column '%s'. It is being cleaned and checked "
            "for missing values." % snp_col
        )
        info["source"] = "study_column:%s" % snp_col

    # ------------------------------------------------------------------
    # 2. clean: ':' becomes '_'
    # ------------------------------------------------------------------
    # The cast to text has to come first.  A numeric identifier column - which
    # is what a file with a plain integer ID gives you - is not a string, and
    # calling a string method on it raises before anything is written.
    before = df.height
    changed = df.select(
        pl.col(snp_col).cast(pl.Utf8, strict=False).str.contains(":", literal=True).sum()
    ).item() or 0
    df = df.with_columns(
        pl.col(snp_col).cast(pl.Utf8, strict=False).str.replace_all(":", "_").alias(snp_col)
    )
    ctx.qc(
        "identifier separator",
        "Colons in the identifier were replaced with underscores, so the identifier survives the "
        "VCF round trip.",
        before, before, changed=int(changed),
    )
    info["identifiers_rewritten"] = int(changed)

    # ------------------------------------------------------------------
    # 3. fill in the identifiers that are missing
    # ------------------------------------------------------------------
    sentinels = _missing_sentinels(policies)
    sentinel_text = ", ".join("'%s'" % s for s in sentinels if s) or "the empty string"
    missing_expr = _missing_expr(snp_col, sentinels)
    n_missing = int(df.select(missing_expr.sum()).item() or 0)
    info["missing_identifiers"] = n_missing
    info["missing_sentinels"] = list(sentinels)

    if n_missing == 0:
        ctx.qc(
            "missing identifiers",
            "Every variant has an identifier. A missing one is one that is empty or one of %s."
            % sentinel_text,
            before, before, changed=0,
        )
        return df, sample_column_dict, info

    if not can_build:
        ctx.warn(
            "%s variants have no usable identifier, and one cannot be built because these columns "
            "are not in the data: %s. They keep the value they have."
            % ("{:,}".format(n_missing), "; ".join(missing_coords))
        )
        ctx.qc(
            "missing identifiers",
            "Identifiers that are empty or one of %s could not be filled in." % sentinel_text,
            before, before, changed=0, warn=True,
        )
        return df, sample_column_dict, info

    df = df.with_columns(
        pl.when(missing_expr)
          .then(_identifier_expr(chr_col, pos_col, ea_col, oa_col))
          .otherwise(pl.col(snp_col))
          .alias(snp_col)
    )
    still_missing = int(df.select(_missing_expr(snp_col, sentinels).sum()).item() or 0)
    info["identifiers_filled"] = n_missing - still_missing
    info["identifiers_still_missing"] = still_missing

    ctx.qc(
        "missing identifiers",
        "%s variants had no usable identifier - empty, or one of %s - and were given "
        "chromosome_position_effect_other instead."
        % ("{:,}".format(n_missing), sentinel_text),
        before, before, changed=n_missing - still_missing,
        warn=bool(still_missing),
    )
    if still_missing:
        ctx.warn(
            "%s variants still have no identifier, because their chromosome, position or alleles "
            "are themselves missing." % "{:,}".format(still_missing)
        )
    return df, sample_column_dict, info
