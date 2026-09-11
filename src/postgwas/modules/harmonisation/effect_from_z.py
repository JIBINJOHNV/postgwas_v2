"""Step 06 (per chromosome) — derive BETA and SE from an imputed Z score.

Public entry point: ``derive_effect_and_standard_error_from_z``.

When neither statistic was supplied, the default ``metal_large_n`` method uses
the large-sample standardized-trait approximation::

    BETA = sigma_y * z / sqrt(2 p (1 - p) Neff)
    SE   = sigma_y     / sqrt(2 p (1 - p) Neff)

This is compatible with sample-size-weighted Z statistics such as those
produced by METAL.  The optional ``zhu_2016`` method uses::

    BETA = sigma_y * z / sqrt(2 p (1 - p) (Neff + z^2))
    SE   = sigma_y     / sqrt(2 p (1 - p) (Neff + z^2))

where ``p`` is the effect allele frequency.  The denominator is *exactly* zero
when ``p`` is 0 or 1, and polars then yields ``inf`` — not null, not an error.
Those endpoints are therefore unusable for this reconstruction even when a
general EAF policy permits them for supplied effect estimates.  The calculation
also requires a finite positive denominator and enough expected additive
genotype variance, measured here by the configurable proxy
``2*p*(1-p)*Neff``.

Neither method reconstructs a study's original regression coefficient without
the assumptions above.  With the default null phenotype standard deviation,
``sigma_y`` is 1 and the result is a standardized additive effect.  A configured
positive phenotype standard deviation scales both BETA and SE into those
phenotype units without changing Z.  Both methods assume genotype variance
``2p(1-p)``, that ``p`` represents the analysed samples and that the supplied
sample size describes the same association test.
For a binary trait, ``Neff`` is the balanced-design effective sample size from
the preceding sample-size step; the result is an approximation, not a direct
recovery of the logistic-regression log odds ratio or a liability-scale effect.

That diploid genotype-variance assumption is not generally valid for
chromosome X. Male non-PAR dosages may be coded 0/1 or 0/2, PAR and non-PAR
regions have different ploidy, and pooled EAF/Neff do not recover the missing
sex-specific information. Consequently the default
``effect_from_z.x_chromosome_z_only_action=fail`` refuses Combination 4 on the
canonical X chromosome. Exact identities involving a supplied BETA or SE are
unaffected. An explicit ``allow_autosomal_assumption`` override preserves the
legacy approximation but records that it is an assumption, not an X-specific
reconstruction.

When BETA and a signed Z score are supplied but SE is missing, a positive SE
exists only when BETA and Z have the same non-zero sign. PostGWAS therefore
uses ``SE = abs(BETA / Z)`` only for sign-consistent finite pairs. Opposite
signs, or zero BETA with non-zero Z, cannot satisfy ``Z = BETA / SE`` for any
positive SE and follow ``effect_from_z.beta_z_sign_mismatch`` instead of being
converted into a negative or fabricated standard error. Before chromosome
fan-out, the complete study is also checked against
``effect_from_z.max_beta_z_sign_mismatch_fraction``. Missing and unusable
statistics do not enter that denominator. A fraction equal to the configured
limit continues to row-level rejection; a larger fraction stops the study.

References:

* Willer CJ et al. METAL. Bioinformatics 2010; PMID 20616382;
  https://pubmed.ncbi.nlm.nih.gov/20616382/
* Rietveld CA et al. Science 2013; PMID 23722424;
  https://pubmed.ncbi.nlm.nih.gov/23722424/
* Okbay A et al. Nature 2016; PMID 27225129;
  https://pubmed.ncbi.nlm.nih.gov/27225129/
* Lee JJ et al. Nature Genetics 2018; PMID 30038396;
  https://pubmed.ncbi.nlm.nih.gov/30038396/
* Zhu Z et al. Nature Genetics 2016; PMID 27019110;
  https://pubmed.ncbi.nlm.nih.gov/27019110/
* Clayton D. Biostatistics 2008; PMID 18441336;
  https://pubmed.ncbi.nlm.nih.gov/18441336/
* PLINK 2 chromosome-X association models;
  https://www.cog-genomics.org/plink/2.0/assoc

Logging and rejected variants
-----------------------------
Pass ``logger`` (a ``PipelineLogger``), ``policies`` (a ``Policies``) and
``rejects`` (a ``RejectCollector``). The pipeline supplies the collector so
every removed variant is written once to the unified reject file.

Give the logger the *same* ``Policies`` object you pass here — the settings
block is read from ``logger.policies``.  ``POLICY_KEYS`` lists every policy this
step reads, so a caller opening the step itself can pass
``policy_keys=POLICY_KEYS`` to ``logger.step(...)``.
"""

from functools import partial
import math
from typing import Any, Dict, Optional, Sequence, Tuple

import polars as pl

from postgwas.core.dataframes import count_matching_rows
from postgwas.core.values import format_number, optional_text

from .shared.runtime import (
    log_info,
    log_warning,
    reject_rows,
    resolve_policies,
)


__all__ = [
    "derive_effect_and_standard_error_from_z",
    "assess_beta_z_sign_discordance",
    "enforce_x_chromosome_z_only_policy",
    "SE_UNAVAILABLE_FROM_Z_COLUMN",
    "BetaZSignMismatchError",
    "XChromosomeZOnlyReconstructionError",
    "describe_effect_from_z_reconstruction",
    "POLICY_KEYS",
    "STEP_LABEL",
]


#: Internal row-level provenance used by the following standard-error step.
#: True means the study supplied no SE and BETA / Z could not recover one for
#: this row because Z or BETA was unusable. It is dropped before export with
#: the other ``__`` working columns.
SE_UNAVAILABLE_FROM_Z_COLUMN = "__se_unavailable_from_z"


#: Every policy this step reads.  Pass to ``logger.step(policy_keys=...)``.
POLICY_KEYS = (
    "eaf.degenerate",
    "effect_from_z.beta_z_sign_mismatch",
    "effect_from_z.x_chromosome_z_only_action",
    "effect_from_z.method",
    "effect_from_z.phenotype_standard_deviation",
    "effect_from_z.minimum_effective_variance",
    "effect_from_z.low_effective_variance_action",
    "validation.z_invalid",
    "sample_size.trait_type",
    "sample_size.min_value",
)

#: The free-text step label recorded in the reject file.
STEP_LABEL = "06 beta_se_from_z"


# Shared mechanics stay outside the scientific implementation.
_info = log_info
_warn = log_warning
_fmt = partial(format_number, pattern="%10.6f", missing="       n/a")
_count = partial(count_matching_rows, null_is_match=True)


# These URLs are part of the scientific provenance, not mutable runtime
# configuration.  The canonical YAML selects the method; this registry keeps
# its equation and cited implementation metadata identical in calculation,
# logging and VCF provenance.
_METHOD_METADATA = {
    "metal_large_n": {
        "label": "METAL-compatible large-N standardized approximation",
        "denominator": "2*EAF*(1-EAF)*Neff",
        "citation": (
            "Willer et al. 2010 (METAL; PMID 20616382); "
            "Lee et al. 2018 (PMID 30038396)"
        ),
        "pubmed_urls": (
            "https://pubmed.ncbi.nlm.nih.gov/20616382/",
            "https://pubmed.ncbi.nlm.nih.gov/23722424/",
            "https://pubmed.ncbi.nlm.nih.gov/27225129/",
            "https://pubmed.ncbi.nlm.nih.gov/30038396/",
        ),
    },
    "zhu_2016": {
        "label": "Zhu et al. standardized approximation",
        "denominator": "2*EAF*(1-EAF)*(Neff+Z^2)",
        "citation": "Zhu et al. 2016 (PMID 27019110)",
        "pubmed_urls": ("https://pubmed.ncbi.nlm.nih.gov/27019110/",),
    },
}


class XChromosomeZOnlyReconstructionError(ValueError):
    """The configured policy forbids assumption-based X reconstruction."""


class BetaZSignMismatchError(ValueError):
    """A missing SE cannot be recovered from incompatible signed statistics."""


def _beta_z_sign_statistic_expressions(
    effect_col: str,
    z_col: str,
    se_col: Optional[str] = None,
    prefix: str = "",
) -> Sequence[pl.Expr]:
    """Return one-pass study-wide BETA/Z sign evidence expressions.

    The source effect may still be an odds ratio at dataset step 07. Both
    possible interpretations are therefore counted in the same Polars
    aggregation, and the already-resolved study effect type selects the
    applicable pair of counts afterwards. For a positive odds ratio, the sign
    of ``log(OR)`` is exactly the sign of ``OR - 1``; no approximation or
    duplicate effect transformation is needed.
    """
    effect = pl.col(effect_col).cast(pl.Float64, strict=False)
    z_score = pl.col(z_col).cast(pl.Float64, strict=False)
    needs_se = (
        pl.col(se_col).cast(pl.Float64, strict=False).is_null()
        if se_col is not None
        else pl.lit(True)
    )
    usable_z = (
        z_score.is_not_null()
        & z_score.is_finite()
        & (z_score != 0.0)
    )
    usable_effect = effect.is_not_null() & effect.is_finite()

    beta_comparable = (needs_se & usable_effect & usable_z).fill_null(False)
    beta_same_sign = (
        ((effect > 0.0) & (z_score > 0.0))
        | ((effect < 0.0) & (z_score < 0.0))
    ).fill_null(False)

    # A non-positive odds ratio is invalid input and follows the existing
    # effect.or_non_positive policy. It cannot provide BETA/Z sign evidence.
    or_comparable = (
        needs_se & usable_effect & (effect > 0.0) & usable_z
    ).fill_null(False)
    or_same_sign = (
        ((effect > 1.0) & (z_score > 0.0))
        | ((effect > 0.0) & (effect < 1.0) & (z_score < 0.0))
    ).fill_null(False)

    return (
        beta_comparable.sum().alias(prefix + "beta_comparable"),
        (beta_comparable & ~beta_same_sign).sum().alias(
            prefix + "beta_discordant"
        ),
        or_comparable.sum().alias(prefix + "odds_ratio_comparable"),
        (or_comparable & ~or_same_sign).sum().alias(
            prefix + "odds_ratio_discordant"
        ),
    )


def assess_beta_z_sign_discordance(
    statistics: Dict[str, Any],
    *,
    effect_type: str,
    effect_col: str,
    z_col: str,
    se_col: Optional[str] = None,
    policies: Optional[Any] = None,
) -> Dict[str, Any]:
    """Apply the configured complete-study BETA/Z discordance guard.

    The denominator contains only rows that need SE recovery and have a finite
    valid effect plus a finite non-zero signed Z. Missing, unparseable and
    non-finite cells cannot dilute evidence of a wrong mapping or sign
    convention. At the exact configured fraction the study continues and the
    chromosome step rejects affected rows; only a strictly larger fraction
    triggers the dataset failure.
    """
    pol = resolve_policies(policies)
    resolved_type = str(effect_type)
    if resolved_type not in ("beta", "odds_ratio"):
        raise ValueError(
            "BETA/Z sign assessment requires effect_type 'beta' or "
            "'odds_ratio', got {!r}.".format(effect_type)
        )

    count_prefix = "beta" if resolved_type == "beta" else "odds_ratio"
    comparable = int(statistics.get(count_prefix + "_comparable") or 0)
    discordant = int(statistics.get(count_prefix + "_discordant") or 0)
    if discordant > comparable:
        raise ValueError(
            "BETA/Z sign evidence is internally inconsistent: {:,} discordant "
            "rows exceed {:,} comparable rows.".format(discordant, comparable)
        )

    maximum_fraction = float(
        pol.get("effect_from_z.max_beta_z_sign_mismatch_fraction")
    )
    action = str(pol.get("effect_from_z.beta_z_sign_mismatch"))
    fraction = discordant / comparable if comparable else 0.0
    result = {
        "applies": True,
        "effect_type": resolved_type,
        "effect_column": effect_col,
        "z_column": z_col,
        "se_column": se_col,
        "comparable_missing_se_rows": comparable,
        "discordant_rows": discordant,
        "discordant_fraction": fraction,
        "maximum_fraction": maximum_fraction,
        "row_action": action,
        "decision": "continue_and_reject" if discordant else "continue",
    }

    fail_on_any = bool(discordant and action == "fail")
    exceeds_fraction = bool(
        comparable and fraction > maximum_fraction
    )
    if fail_on_any or exceeds_fraction:
        if fail_on_any:
            trigger = (
                "effect_from_z.beta_z_sign_mismatch='fail' stops on any "
                "discordant row"
            )
        else:
            trigger = (
                "the fraction exceeds "
                "effect_from_z.max_beta_z_sign_mismatch_fraction={:g} "
                "({:.2%})".format(maximum_fraction, maximum_fraction)
            )
        raise BetaZSignMismatchError(
            "Study-wide BETA/Z sign validation found {:,} discordant row(s) "
            "among {:,} comparable row(s) needing SE recovery ({:.4%}); {}. "
            "A positive SE cannot satisfy Z = BETA / SE for opposite signs, "
            "or for zero BETA with non-zero Z. Check effect column {!r}, Z-score "
            "column {!r}, their allele orientation and the declared effect type. "
            "No chromosome processing was started."
            .format(
                discordant,
                comparable,
                fraction,
                trigger,
                effect_col,
                z_col,
            )
        )

    return result


def _recover_positive_se_from_signed_z(
    chromosome: str,
    df: pl.DataFrame,
    beta_col: str,
    z_col: str,
    output_se_col: str,
    *,
    policies: Any,
    logger: Optional[Any],
    rejects: Optional[Any],
) -> Tuple[pl.DataFrame, Dict[str, int]]:
    """Recover missing SE cells without hiding BETA/Z direction conflicts.

    A standard error is strictly positive and ``Z = BETA / SE``. Therefore a
    finite, non-zero signed Z can recover SE only when BETA has the same
    non-zero sign. This shared implementation serves both a wholly absent SE
    column and null cells in a supplied SE column.
    """
    action = str(policies.get("effect_from_z.beta_z_sign_mismatch"))
    has_output = output_se_col in df.columns
    needs_se = pl.col(output_se_col).is_null() if has_output else pl.lit(True)
    beta = pl.col(beta_col)
    z_score = pl.col(z_col)
    finite_nonzero_z_pair = (
        needs_se
        & beta.is_not_null()
        & beta.is_finite()
        & z_score.is_not_null()
        & z_score.is_finite()
        & (z_score != 0.0)
    ).fill_null(False)
    same_nonzero_sign = (
        ((beta > 0.0) & (z_score > 0.0))
        | ((beta < 0.0) & (z_score < 0.0))
    ).fill_null(False)
    sign_mismatch = finite_nonzero_z_pair & ~same_nonzero_sign
    n_sign_mismatch = _count(df, sign_mismatch)

    if n_sign_mismatch and action == "fail":
        raise BetaZSignMismatchError(
            "Chromosome {}: {:,} variant(s) need a standard error but their "
            "supplied BETA and signed Z score are directionally incompatible: "
            "the signs are opposite, or BETA is zero while Z is non-zero. No "
            "positive SE can satisfy Z = BETA / SE. Policy "
            "effect_from_z.beta_z_sign_mismatch is 'fail', so no affected SE "
            "was calculated. Correct the BETA/Z columns or set the policy to "
            "'reject' to remove and record only the affected variants."
            .format(chromosome, n_sign_mismatch)
        )

    n_rejected = 0
    if n_sign_mismatch:
        df, n_rejected = reject_rows(
            df,
            sign_mismatch,
            step_label=STEP_LABEL,
            reason="beta_z_sign_discordant",
            context=logger,
            collector=rejects,
            detail=(
                "columns '{}' and '{}'; policy "
                "effect_from_z.beta_z_sign_mismatch='reject'"
            ).format(beta_col, z_col),
            check_name="BETA/Z sign agreement for SE recovery",
            description=(
                "BETA and the signed Z score have incompatible directions while "
                "SE is missing, so no positive standard error can satisfy "
                "Z = BETA / SE."
            ),
            warn_on_remove=True,
            warn_without_collector=True,
        )

    # Rebuild expressions against the surviving frame. The absolute value is
    # a numerical-domain safeguard; sign agreement is what makes this an exact
    # recovery instead of a silent correction of contradictory inputs.
    needs_se = pl.col(output_se_col).is_null() if has_output else pl.lit(True)
    beta = pl.col(beta_col)
    z_score = pl.col(z_col)
    usable_for_se = (
        needs_se
        & beta.is_not_null()
        & beta.is_finite()
        & z_score.is_not_null()
        & z_score.is_finite()
        & (z_score != 0.0)
        & (
            ((beta > 0.0) & (z_score > 0.0))
            | ((beta < 0.0) & (z_score < 0.0))
        )
    ).fill_null(False)
    n_recovered = _count(df, usable_for_se)
    existing_se = pl.col(output_se_col) if has_output else pl.lit(None)
    df = df.with_columns(
        pl.when(usable_for_se)
        .then((beta / z_score).abs())
        .otherwise(existing_se)
        .cast(pl.Float64, strict=False)
        .alias(output_se_col)
    )
    unavailable = pl.col(output_se_col).is_null()
    n_unavailable = _count(df, unavailable)
    df = df.with_columns(unavailable.alias(SE_UNAVAILABLE_FROM_Z_COLUMN))

    return df, {
        "beta_z_sign_mismatches": n_sign_mismatch,
        "beta_z_sign_mismatches_rejected": n_rejected,
        "se_recovered_from_z": n_recovered,
        "se_unavailable_from_z": n_unavailable,
    }


def _requires_z_only_reconstruction(
    columns: Sequence[str], sample_column_dict: dict,
) -> bool:
    """Return whether the current schema reaches Combination 4.

    This schema check matches the four-way dispatch in the public function. It
    performs no row scan and can therefore be reused in dataset preflight
    before chromosome partitions are written.
    """
    available = set(columns)
    z_column = optional_text(sample_column_dict.get("imp_z_col"))
    effect_column = optional_text(sample_column_dict.get("beta_or_col"))
    se_column = optional_text(sample_column_dict.get("se_col"))
    return bool(
        z_column is not None
        and z_column in available
        and (effect_column is None or effect_column not in available)
        and (se_column is None or se_column not in available)
    )


def enforce_x_chromosome_z_only_policy(
    chromosome: str,
    columns: Sequence[str],
    sample_column_dict: dict,
    *,
    policies: Optional[Any] = None,
    logger: Optional[Any] = None,
    scope: Optional[str] = None,
) -> Dict[str, Any]:
    """Enforce policy before assumption-based chromosome-X reconstruction.

    Chromosome labels are canonical by the time this step runs. Dataset
    preflight and the chromosome worker share this function so the early
    failure and row-processing invariant cannot drift apart.
    """
    pol = resolve_policies(policies)
    action = str(pol.get("effect_from_z.x_chromosome_z_only_action"))
    applies = (
        str(chromosome).strip().upper() == "X"
        and _requires_z_only_reconstruction(columns, sample_column_dict)
    )
    decision = {
        "applies": applies,
        "chromosome": "X" if applies else str(chromosome),
        "action": action,
        "assumption": (
            "diploid_autosomal_genotype_variance_2p1mp"
            if applies and action == "allow_autosomal_assumption"
            else None
        ),
    }
    if not applies:
        return decision

    location = scope or "The study"
    explanation = (
        "{} contains chromosome X and supplies a Z score but neither BETA nor "
        "SE. PostGWAS cannot scientifically reconstruct chromosome-X BETA and "
        "SE from pooled EAF and Neff with the diploid autosomal variance "
        "2*EAF*(1-EAF): male dosage coding (0/1 versus 0/2), sex composition, "
        "sex-specific allele frequency and PAR status are not available."
    ).format(location)
    if action == "fail":
        raise XChromosomeZOnlyReconstructionError(
            explanation
            + " No chromosome-X BETA or SE was calculated. Provide a BETA "
            "column, provide an SE column, or exclude X by setting both "
            "chromosome.allowed and chromosome.allowed_after_split to lists "
            "containing chromosomes 1 through 22 only. "
            "Only when this autosomal assumption is knowingly acceptable, set "
            "effect_from_z.x_chromosome_z_only_action to "
            "'allow_autosomal_assumption'. The default is 'fail'."
        )

    _warn(
        logger,
        "SCIENTIFIC OVERRIDE: {} The explicit policy "
        "effect_from_z.x_chromosome_z_only_action='allow_autosomal_assumption' "
        "permits the autosomal approximation. The resulting chromosome-X "
        "BETA and SE are assumption-dependent and are not an X-specific "
        "reconstruction.".format(explanation),
    )
    return decision


def describe_effect_from_z_reconstruction(
    method: str,
    phenotype_standard_deviation: Optional[float] = None,
) -> Dict[str, Any]:
    """Return one authoritative description of a configured reconstruction.

    ``phenotype_standard_deviation=None`` denotes the standardized-trait
    assumption ``sigma_y=1``.  The returned strings are used in the canonical
    log and VCF header, so provenance cannot drift from the calculation.
    """
    try:
        registered = _METHOD_METADATA[str(method)]
    except KeyError as exc:
        raise ValueError(
            "Unknown effect_from_z.method {!r}; expected one of {}.".format(
                method, ", ".join(sorted(_METHOD_METADATA))
            )
        ) from exc

    if phenotype_standard_deviation is None:
        scale = 1.0
        beta_formula = "BETA = Z / sqrt({})".format(registered["denominator"])
        se_formula = "SE = 1 / sqrt({})".format(registered["denominator"])
        output_scale = "standardized additive effect approximation"
        scale_name = "standardized_effect_estimate"
        scale_source = "unit_phenotype_standard_deviation_assumption"
    else:
        scale = float(phenotype_standard_deviation)
        if not math.isfinite(scale) or scale <= 0.0:
            raise ValueError(
                "effect_from_z.phenotype_standard_deviation must be a finite "
                "number greater than zero when provided."
            )
        beta_formula = (
            "BETA = phenotype_SD * Z / sqrt({}); phenotype_SD={}".format(
                registered["denominator"], scale
            )
        )
        se_formula = (
            "SE = phenotype_SD / sqrt({}); phenotype_SD={}".format(
                registered["denominator"], scale
            )
        )
        output_scale = "additive effect approximation in configured phenotype units"
        scale_name = "phenotype_scale_effect_estimate"
        scale_source = "configured_phenotype_standard_deviation"

    return {
        "method": str(method),
        "label": registered["label"],
        "denominator": registered["denominator"],
        "beta_formula": beta_formula,
        "se_formula": se_formula,
        "output_scale": output_scale,
        "scale_name": scale_name,
        "scale_source": scale_source,
        "scale_multiplier": scale,
        "citation": registered["citation"],
        "pubmed_urls": tuple(registered["pubmed_urls"]),
    }


# public entry point
# ---------------------------------------------------------------------------
def derive_effect_and_standard_error_from_z(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    logger: Optional[Any] = None,
    policies: Optional[Any] = None,
    rejects: Optional[Any] = None,
) -> Tuple[pl.DataFrame, Dict, dict]:
    """Fill in whichever of BETA and SE the study did not supply.

    Four combinations of BETA/SE presence are handled explicitly:

    ==============  =========  ==================================================
    BETA present    SE present  action
    ==============  =========  ==================================================
    yes             yes        preserve populated SE; fill null cells with
                               abs(BETA / Z) only after signed agreement
    yes             no         SE = abs(BETA / Z) after signed agreement, otherwise
                               left to ``derive_standard_error_from_effect_and_p_value``
    no              yes        BETA = Z * SE  (needs a Z score)
    no              no         both derived from Z, EAF and Neff
    ==============  =========  ==================================================
    """
    sample_column_dict = sample_column_dict.copy()  # SAFE
    pol = resolve_policies(policies)

    qc_info = {"initial_variants": df.height}  # type: Dict[str, Any]

    if logger is not None:
        logger.settings(list(POLICY_KEYS))

    # -------------------------------------------------------
    # Column mappings
    # -------------------------------------------------------
    imp_z_col = optional_text(sample_column_dict.get("imp_z_col"))
    beta_col = optional_text(sample_column_dict.get("beta_or_col"))
    se_col = optional_text(sample_column_dict.get("se_col"))
    eaf_col = optional_text(sample_column_dict.get("eaf_col"))

    has_z = imp_z_col is not None and imp_z_col in df.columns
    has_beta = beta_col is not None and beta_col in df.columns
    has_se = se_col is not None and se_col in df.columns

    for name, column in (
        ("beta_or_col", beta_col),
        ("se_col", se_col),
        ("imp_z_col", imp_z_col),
    ):
        if column is not None and column not in df.columns:
            _warn(
                logger,
                "The config sets {} = '{}' but that column is not in the data for "
                "chromosome {}, so it is treated as absent.".format(name, column, chromosome),
            )

    qc_info.update({
        "has_beta": has_beta,
        "has_se": has_se,
        "has_zscore": has_z,
        "effect_column": beta_col if has_beta else None,
        "se_column": se_col if has_se else None,
        "z_column": imp_z_col if has_z else None,
    })

    # -------------------------------------------------------
    # Combination 1 of 4 — BETA and SE both supplied.
    # -------------------------------------------------------
    if has_beta and has_se:
        df = df.with_columns(
            [
                pl.col(beta_col).cast(pl.Float64, strict=False),
                pl.col(se_col).cast(pl.Float64, strict=False),
            ]
            + (
                [pl.col(imp_z_col).cast(pl.Float64, strict=False)]
                if has_z else []
            )
        )
        missing_se = pl.col(se_col).is_null()
        n_missing_se = int(
            df.select(missing_se.sum().alias("count")).item() or 0
        )
        n_recovered_from_z = 0
        n_unavailable_from_z = n_missing_se
        n_sign_mismatch = 0
        n_sign_mismatch_rejected = 0
        if has_z:
            df, recovery = _recover_positive_se_from_signed_z(
                chromosome,
                df,
                beta_col,
                imp_z_col,
                se_col,
                policies=pol,
                logger=logger,
                rejects=rejects,
            )
            n_recovered_from_z = recovery["se_recovered_from_z"]
            n_unavailable_from_z = recovery["se_unavailable_from_z"]
            n_sign_mismatch = recovery["beta_z_sign_mismatches"]
            n_sign_mismatch_rejected = recovery[
                "beta_z_sign_mismatches_rejected"
            ]
            _info(
                logger,
                "The study supplies both an effect size ('{}') and a standard error "
                "column ('{}'). Populated SE values were preserved; {} null cell(s) "
                "were recovered as positive abs(BETA / Z) values after signed "
                "BETA/Z agreement against '{}'; {} incompatible cell(s) followed "
                "effect_from_z.beta_z_sign_mismatch='{}', and {} remain for the "
                "configured p-value fallback.".format(
                    beta_col,
                    se_col,
                    n_recovered_from_z,
                    imp_z_col,
                    n_sign_mismatch,
                    pol.get("effect_from_z.beta_z_sign_mismatch"),
                    n_unavailable_from_z,
                ),
            )
        else:
            _info(
                logger,
                "The study supplies both an effect size ('{}') and a standard error "
                "column ('{}'). Populated SE values were preserved; {} null cell(s) "
                "remain for the configured p-value fallback because no Z score "
                "column is available.".format(beta_col, se_col, n_missing_se),
            )
        qc_info.update(
            {
                "status": (
                    "Partial SE recovered from BETA and Z"
                    if n_recovered_from_z else "Beta and SE columns are present"
                ),
                "beta_computed": False,
                "se_computed": n_recovered_from_z > 0,
                "se_missing_on_entry": n_missing_se,
                "se_recovered_from_z": n_recovered_from_z,
                "se_unavailable_from_z": n_unavailable_from_z,
                "beta_z_sign_mismatch_action": pol.get(
                    "effect_from_z.beta_z_sign_mismatch"
                ),
                "beta_z_sign_mismatches": n_sign_mismatch,
                "removed_beta_z_sign_mismatches": n_sign_mismatch_rejected,
                "variants_removed_total": n_sign_mismatch_rejected,
                "variants_remaining": df.height,
            }
        )
        return df, qc_info, sample_column_dict

    # -------------------------------------------------------
    # Nothing to work from: no effect size and no Z score.
    # -------------------------------------------------------
    if not has_beta and not has_z:
        qc_info["status"] = "skipped_no_zscore_or_effectsize_column"
        raise KeyError(
            "Chromosome {}: the summary statistics have neither an effect size "
            "column nor a Z score column, so BETA and SE cannot be produced. "
            "Provide beta_or_col (with se_col) or imp_z_col.".format(chromosome)
        )

    # -------------------------------------------------------
    # Combination 2 of 4 — BETA supplied, SE missing.
    # -------------------------------------------------------
    if has_beta and not has_se:
        if not has_z:
            _info(
                logger,
                "The study supplies an effect size ('{}') but no standard error, and "
                "there is no Z score column to recover one from. SE is calculated in a "
                "later step from BETA and the p-value.".format(beta_col),
            )
            qc_info.update(
                {
                    "status": "Beta column is present, but SE column is missing",
                    "beta_computed": False,
                    "se_computed": False,
                    "variants_removed_total": 0,
                    "variants_remaining": df.height,
                }
            )
            return df, qc_info, sample_column_dict

        df = df.with_columns(
            [
                pl.col(beta_col).cast(pl.Float64, strict=False),
                pl.col(imp_z_col).cast(pl.Float64, strict=False),
            ]
        )
        # A positive SE can satisfy Z = BETA / SE only for same-sign,
        # non-zero BETA/Z pairs. The shared helper handles incompatible signs
        # before calculating abs(BETA / Z), and leaves other unusable inputs
        # null for the configured p-value fallback.
        df, recovery = _recover_positive_se_from_signed_z(
            chromosome,
            df,
            beta_col,
            imp_z_col,
            "SE",
            policies=pol,
            logger=logger,
            rejects=rejects,
        )
        n_recovered = recovery["se_recovered_from_z"]
        n_undefined = recovery["se_unavailable_from_z"]
        n_sign_mismatch = recovery["beta_z_sign_mismatches"]
        n_sign_mismatch_rejected = recovery[
            "beta_z_sign_mismatches_rejected"
        ]

        sample_column_dict["se_col"] = "SE"
        _info(
            logger,
            "The study supplies an effect size ('{}') but no standard error, so SE was "
            "recovered as positive abs(BETA / Z) after requiring signed agreement "
            "with Z score column '{}'. {} incompatible variant(s) followed "
            "effect_from_z.beta_z_sign_mismatch='{}'.".format(
                beta_col,
                imp_z_col,
                n_sign_mismatch,
                pol.get("effect_from_z.beta_z_sign_mismatch"),
            ),
        )
        if logger is not None:
            logger.qc(
                "standard error recovered from Z",
                "SE was derived as abs(BETA / Z) only for finite, non-zero, "
                "same-sign BETA/Z pairs.",
                df.height,
                df.height,
                changed=n_recovered,
                step=STEP_LABEL,
                warn=n_undefined > 0,
            )
        if n_undefined > 0:
            _warn(
                logger,
                "{:,} variants have an unusable BETA or a zero, missing or non-finite "
                "Z score, so no standard error could be derived for them. The "
                "following standard-error step applies "
                "pvalue.zero_missing_se to any row whose raw p-value was zero and, "
                "when pvalue.derive_partial_missing_se is enabled, attempts to fill "
                "other null SE cells from BETA and p. Values that remain unresolved "
                "reach the configured standard-error validation gate."
                .format(n_undefined),
            )
        qc_info.update(
            {
                "status": "SE derived from BETA and Z",
                "beta_computed": False,
                "se_computed": n_recovered > 0,
                "se_recovered_from_z": n_recovered,
                "se_undefined_from_z": n_undefined,
                "beta_z_sign_mismatch_action": pol.get(
                    "effect_from_z.beta_z_sign_mismatch"
                ),
                "beta_z_sign_mismatches": n_sign_mismatch,
                "removed_beta_z_sign_mismatches": n_sign_mismatch_rejected,
                "variants_removed_total": n_sign_mismatch_rejected,
                "variants_remaining": df.height,
            }
        )
        return df, qc_info, sample_column_dict

    # -------------------------------------------------------
    # Combination 3 of 4 — SE supplied, BETA missing.
    # (has_z is guaranteed here: the no-BETA-no-Z case raised above.)
    # -------------------------------------------------------
    if has_se and not has_beta:
        df = df.with_columns(
            [
                pl.col(se_col).cast(pl.Float64, strict=False),
                pl.col(imp_z_col).cast(pl.Float64, strict=False),
            ]
        )
        # BETA = Z * SE is exact and needs no allele frequency.
        df = df.with_columns(
            pl.when(pl.col(imp_z_col).is_finite() & pl.col(se_col).is_finite())
            .then(pl.col(imp_z_col) * pl.col(se_col))
            .otherwise(None)
            .alias("BETA")
        )
        n_undefined = int(df.select(pl.col("BETA").is_null().sum()).item() or 0)

        sample_column_dict["beta_or_col"] = "BETA"
        sample_column_dict["beta_col"] = "BETA"
        _info(
            logger,
            "The study supplies a standard error ('{}') but no effect size, so BETA was "
            "recovered exactly as Z * SE from the Z score column '{}'. No allele "
            "frequency is needed for this recovery.".format(se_col, imp_z_col),
        )
        if logger is not None:
            logger.qc(
                "effect size recovered from Z",
                "BETA was derived as Z * SE for every variant with a usable Z score "
                "and standard error.",
                df.height,
                df.height,
                changed=df.height - n_undefined,
                step=STEP_LABEL,
                warn=n_undefined > 0,
            )
        if n_undefined > 0:
            _warn(
                logger,
                "{:,} variants have a non-finite Z score or standard error, so no effect "
                "size could be derived for them; their BETA is empty and the effect size "
                "check removes them.".format(n_undefined),
            )
        qc_info.update(
            {
                "status": "BETA derived from Z and SE",
                "beta_computed": True,
                "se_computed": False,
                "beta_undefined_from_z": n_undefined,
                "variants_removed_total": 0,
                "variants_remaining": df.height,
            }
        )
        return df, qc_info, sample_column_dict

    # -------------------------------------------------------
    # Combination 4 of 4 — neither BETA nor SE; derive both from Z.
    # -------------------------------------------------------
    x_reconstruction_policy = enforce_x_chromosome_z_only_policy(
        chromosome,
        df.columns,
        sample_column_dict,
        policies=pol,
        logger=logger,
    )
    _info(
        logger,
        "Neither an effect size nor a standard error was supplied, so both are derived "
        "from the Z score column '{}' using the effect allele frequency '{}' and the "
        "effective sample size.".format(imp_z_col, eaf_col),
    )

    if eaf_col is None or eaf_col not in df.columns:
        raise KeyError(
            "Chromosome {}: BETA and SE must be derived from the Z score column '{}', "
            "which needs an effect allele frequency, but the configured column '{}' is "
            "not present.".format(chromosome, imp_z_col, eaf_col)
        )

    if "Neff" not in df.columns:
        raise KeyError(
            "Chromosome {}: BETA and SE must be derived from the Z score column '{}', "
            "which needs the effective sample size, but the 'Neff' column is not "
            "present. It is produced by the sample size step.".format(chromosome, imp_z_col)
        )

    # -------------------------------------------------------
    # Cast numeric
    # -------------------------------------------------------
    df = df.with_columns(
        [
            pl.col(imp_z_col).cast(pl.Float64, strict=False),
            pl.col(eaf_col).cast(pl.Float64, strict=False),
            pl.col("Neff").cast(pl.Float64, strict=False),
        ]
    )

    z_invalid_action = pol.get("validation.z_invalid")
    degenerate_action = pol.get("eaf.degenerate")
    reconstruction = describe_effect_from_z_reconstruction(
        str(pol.get("effect_from_z.method")),
        pol.get("effect_from_z.phenotype_standard_deviation"),
    )
    minimum_effective_variance = float(
        pol.get("effect_from_z.minimum_effective_variance")
    )
    low_effective_variance_action = str(
        pol.get("effect_from_z.low_effective_variance_action")
    )
    trait_type = str(pol.get("sample_size.trait_type"))
    min_sample_size = pol.get("sample_size.min_value")

    n_initial = qc_info["initial_variants"]

    # --- 1. Z score ---------------------------------------------------------
    cond_z_invalid = pl.col(imp_z_col).is_null() | ~pl.col(imp_z_col).is_finite()
    n_z_null = 0
    if z_invalid_action == "fail":
        n_z_null = _count(df, cond_z_invalid)
        if n_z_null:
            raise ValueError(
                "Chromosome {}: {:,} variants have a missing or non-finite Z score in "
                "column '{}' and policy validation.z_invalid is 'fail'.".format(
                    chromosome, n_z_null, imp_z_col
                )
            )
    else:
        if z_invalid_action == "keep":
            n_z_null = _count(df, cond_z_invalid)
            if n_z_null:
                _warn(
                    logger,
                    "Policy validation.z_invalid is 'keep', which may retain a row "
                    "when effect statistics already exist. It cannot retain a "
                    "missing or non-finite Z score when both BETA and SE must be "
                    "reconstructed from Z; {:,} affected variants are rejected here "
                    "with reason z_invalid.".format(n_z_null),
                )
        df, n_z_null = reject_rows(
            df,
            cond_z_invalid,
            step_label=STEP_LABEL, reason="z_invalid",
            check_name="Z score usable",
            description=(
                "The Z score is missing or not finite, so neither BETA nor SE "
                "can be reconstructed from it."
            ),
            collector=rejects, context=logger,
            detail="column '{}'".format(imp_z_col),
        )

    # --- 2. Effect allele frequency: missing --------------------------------
    cond_eaf_null = pl.col(eaf_col).is_null() | ~pl.col(eaf_col).is_finite()
    df, n_eaf_null = reject_rows(
        df,
        cond_eaf_null,
        step_label=STEP_LABEL, reason="eaf_null",
        check_name="effect allele frequency present",
        description="The effect allele frequency is missing or not a finite number, so the denominator of the effect size formula cannot be evaluated.",
        collector=rejects, context=logger,
        detail="column '{}'".format(eaf_col),
    )

    # --- 3. Effect allele frequency: degenerate -----------------------------
    # At exactly 0 or 1 both configured denominators are exactly zero and
    # polars yields inf, so the predicate must be <= 0 / >= 1, not < 0 / > 1.
    cond_eaf_lt_zero = pl.col(eaf_col) <= 0
    cond_eaf_gt_one = pl.col(eaf_col) >= 1
    cond_eaf_degenerate = cond_eaf_lt_zero | cond_eaf_gt_one

    edge = df.select(
        [
            cond_eaf_lt_zero.fill_null(False).sum().alias("le_zero"),
            cond_eaf_gt_one.fill_null(False).sum().alias("ge_one"),
        ]
    ).to_dicts()[0]
    n_eaf_lt_zero = int(edge["le_zero"] or 0)
    n_eaf_gt_one = int(edge["ge_one"] or 0)
    n_eaf_invalid = 0

    if degenerate_action == "fail":
        if n_eaf_lt_zero + n_eaf_gt_one:
            raise ValueError(
                "Chromosome {}: {:,} variants have an effect allele frequency of exactly "
                "0 or 1 (or outside that range) in column '{}' and policy eaf.degenerate "
                "is 'fail'.".format(chromosome, n_eaf_lt_zero + n_eaf_gt_one, eaf_col)
            )
    else:
        if (
            degenerate_action in ("keep", "null")
            and n_eaf_lt_zero + n_eaf_gt_one
        ):
            _warn(
                logger,
                "Policy eaf.degenerate is '{}', which may retain an endpoint EAF when "
                "BETA and SE were supplied. It cannot retain one for Z-only effect "
                "reconstruction because the denominator is zero; {:,} affected variants "
                "are rejected here with reason eaf_degenerate.".format(
                    degenerate_action, n_eaf_lt_zero + n_eaf_gt_one
                ),
            )
        df, n_eaf_invalid = reject_rows(
            df,
            cond_eaf_degenerate,
            step_label=STEP_LABEL, reason="eaf_degenerate",
            check_name="degenerate effect allele frequency",
            description="The effect allele frequency is exactly 0 or 1, which makes the standard error formula divide by zero and yield an infinite value.",
            collector=rejects, context=logger,
            detail="{:,} at or below 0, {:,} at or above 1".format(
                n_eaf_lt_zero, n_eaf_gt_one
            ),
        )

    # --- 4. Effective sample size -------------------------------------------
    cond_neff_null = pl.col("Neff").is_null() | ~pl.col("Neff").is_finite()
    df, n_neff_null = reject_rows(
        df,
        cond_neff_null,
        step_label=STEP_LABEL, reason="neff_invalid",
        check_name="effective sample size present",
        description="The effective sample size is missing or not a finite number, so the denominator of the effect size formula cannot be evaluated.",
        collector=rejects, context=logger,
        detail="Neff missing or non-finite",
    )

    # sample_size.min_value sets the minimum usable effective sample size.
    cond_neff_invalid = (pl.col("Neff") <= 0) | (pl.col("Neff") < min_sample_size)
    df, n_neff_invalid = reject_rows(
        df,
        cond_neff_invalid,
        step_label=STEP_LABEL, reason="neff_invalid",
        check_name="effective sample size usable",
        description="The effective sample size is not a positive number, so the effect size formula would divide by zero or return a complex value.",
        collector=rejects, context=logger,
        detail="Neff at or below {}".format(min_sample_size)
        if min_sample_size > 0
        else "Neff at or below 0",
    )

    # --- 5. Effective genotype-variance information ------------------------
    # This is a reconstruction-specific reliability proxy, not an exact minor
    # allele count: Neff may differ from total N and EAF may come from a
    # reference panel.  It is deliberately applied only when both BETA and SE
    # must be approximated from Z.
    effective_variance_col = "__effect_from_z_effective_variance"
    df = df.with_columns(
        (
            2.0
            * pl.col(eaf_col)
            * (1.0 - pl.col(eaf_col))
            * pl.col("Neff")
        ).alias(effective_variance_col)
    )
    effective_variance_summary = df.select([
        pl.col(effective_variance_col).min().alias("minimum"),
        pl.col(effective_variance_col).max().alias("maximum"),
        (
            pl.col(effective_variance_col) < minimum_effective_variance
        ).fill_null(True).sum().alias("below_floor"),
    ]).to_dicts()[0]
    n_low_effective_variance = int(
        effective_variance_summary["below_floor"] or 0
    )
    low_effective_variance = (
        pl.col(effective_variance_col) < minimum_effective_variance
    ).fill_null(True)
    if low_effective_variance_action == "fail" and n_low_effective_variance:
        raise ValueError(
            "Chromosome {}: {:,} variants have 2*EAF*(1-EAF)*Neff below the "
            "configured minimum {} required for Z-only BETA/SE reconstruction. "
            "Policy effect_from_z.low_effective_variance_action is 'fail'. This "
            "quantity is an effective genotype-variance proxy, not an exact minor "
            "allele count. Check EAF and Neff, or change "
            "effect_from_z.minimum_effective_variance explicitly.".format(
                chromosome,
                n_low_effective_variance,
                minimum_effective_variance,
            )
        )
    df, n_low_effective_variance_removed = reject_rows(
        df,
        low_effective_variance,
        step_label=STEP_LABEL,
        reason="effect_from_z_low_effective_variance",
        check_name="effective variance for Z-only effect reconstruction",
        description=(
            "The proxy 2*EAF*(1-EAF)*Neff is below the configured minimum, so "
            "the Z-only BETA and SE approximation is not reliable."
        ),
        collector=rejects,
        context=logger,
        detail="minimum {}".format(minimum_effective_variance),
    )

    # --- 6. Full denominator ------------------------------------------------
    # The row-level guard is an invariant even when an upstream policy has
    # retained a problematic value. For zhu_2016 it also catches overflow from
    # an extreme but finite Z score. No division may use a non-positive or
    # non-finite denominator.
    denominator_col = "__effect_from_z_denominator"
    sample_information = pl.col("Neff")
    if reconstruction["method"] == "zhu_2016":
        sample_information = sample_information + pl.col(imp_z_col) ** 2
    df = df.with_columns(
        (
            2.0
            * pl.col(eaf_col)
            * (1.0 - pl.col(eaf_col))
            * sample_information
        ).alias(denominator_col)
    )
    denominator_invalid = (
        pl.col(denominator_col).is_null()
        | ~pl.col(denominator_col).is_finite()
        | (pl.col(denominator_col) <= 0.0)
    )
    df, n_denominator_invalid = reject_rows(
        df,
        denominator_invalid,
        step_label=STEP_LABEL,
        reason="effect_from_z_invalid_denominator",
        check_name="finite positive Z-reconstruction denominator",
        description=(
            "The Z-only BETA/SE denominator is missing, non-finite or not "
            "positive, so division would produce an invalid statistic."
        ),
        collector=rejects,
        context=logger,
        detail=reconstruction["denominator"],
    )

    # Every removal above came out of `df`, so the difference is the truth and
    # cannot drift from the per-reason counts.
    n_final = df.height
    n_removed = n_initial - n_final
    pct_removed = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

    qc_info.update(
        {
            "removed_z_null": n_z_null,
            "removed_eaf_null": n_eaf_null,
            "removed_eaf_invalid": n_eaf_invalid,
            "removed_eaf_lt_zero": n_eaf_lt_zero,
            "removed_eaf_gt_one": n_eaf_gt_one,
            "removed_neff_null": n_neff_null,
            "removed_neff_invalid": n_neff_invalid,
            "minimum_effective_variance": minimum_effective_variance,
            "low_effective_variance_action": low_effective_variance_action,
            "effective_variance_min_before_filter": effective_variance_summary["minimum"],
            "effective_variance_max_before_filter": effective_variance_summary["maximum"],
            "low_effective_variance_detected": n_low_effective_variance,
            "removed_low_effective_variance": n_low_effective_variance_removed,
            "removed_invalid_denominator": n_denominator_invalid,
            "variants_removed_total": n_removed,
            "variants_remaining": n_final,
        }
    )

    if n_final == 0:
        raise ValueError(
            "Chromosome {}: zero usable variants remain, so BETA and SE cannot be "
            "derived from the Z score. Started with {:,} variants and removed all of "
            "them ({:,} unusable Z score, {:,} missing frequency, {:,} degenerate "
            "frequency, {:,} unusable effective sample size, {:,} below the effective-"
            "variance floor, {:,} invalid denominator). Check the Z score, frequency "
            "and sample size columns and the effect_from_z policies for this chromosome."
            .format(
                chromosome,
                n_initial,
                n_z_null,
                n_eaf_null,
                n_eaf_invalid,
                n_neff_null + n_neff_invalid,
                n_low_effective_variance_removed,
                n_denominator_invalid,
            )
        )

    # -------------------------------------------------------
    # COMPUTE BETA AND SE
    # -------------------------------------------------------
    usable_denominator = (
        pl.col(denominator_col).is_finite()
        & (pl.col(denominator_col) > 0.0)
    ).fill_null(False)
    scale_multiplier = float(reconstruction["scale_multiplier"])
    df = df.with_columns(
        [
            pl.when(usable_denominator)
            .then(
                scale_multiplier
                * pl.col(imp_z_col)
                / pl.col(denominator_col).sqrt()
            )
            .otherwise(None)
            .alias("BETA"),
            pl.when(usable_denominator)
            .then(scale_multiplier / pl.col(denominator_col).sqrt())
            .otherwise(None)
            .alias("SE"),
        ]
    ).drop([denominator_col, effective_variance_col])

    invalid_reconstruction = (
        pl.col("BETA").is_null()
        | ~pl.col("BETA").is_finite()
        | pl.col("SE").is_null()
        | ~pl.col("SE").is_finite()
        | (pl.col("SE") <= 0.0)
    )
    n_invalid_reconstruction = _count(df, invalid_reconstruction)
    if n_invalid_reconstruction:
        raise ValueError(
            "Chromosome {}: the configured Z-only reconstruction produced a "
            "non-finite BETA or non-positive/non-finite SE for {:,} variant(s) "
            "despite a valid denominator. Check Z, EAF, Neff and "
            "effect_from_z.phenotype_standard_deviation; no reconstructed "
            "statistics were accepted.".format(
                chromosome, n_invalid_reconstruction
            )
        )

    sample_column_dict["beta_or_col"] = "BETA"
    sample_column_dict["beta_col"] = "BETA"
    sample_column_dict["se_col"] = "SE"
    beta_col = "BETA"
    se_col = "SE"
    qc_info["beta_computed"] = True
    qc_info["se_computed"] = True
    qc_info["effect_estimate_method"] = reconstruction["method"]
    qc_info["effect_estimate_method_label"] = reconstruction["label"]
    qc_info["effect_estimate_scale"] = reconstruction["scale_name"]
    qc_info["effect_estimate_scale_source"] = reconstruction["scale_source"]
    qc_info["effect_estimate_phenotype_standard_deviation"] = scale_multiplier
    qc_info["effect_estimate_beta_formula"] = reconstruction["beta_formula"]
    qc_info["effect_estimate_se_formula"] = reconstruction["se_formula"]
    qc_info["effect_estimate_citation"] = reconstruction["citation"]
    qc_info["effect_estimate_pubmed_urls"] = list(reconstruction["pubmed_urls"])
    qc_info["invalid_reconstructed_statistics"] = 0
    qc_info["effect_estimate_trait_type"] = trait_type
    qc_info["x_chromosome_z_only_policy"] = x_reconstruction_policy
    qc_info["effect_estimate_uses_neff"] = True
    qc_info["effect_estimate_eaf_assumption"] = "same_analysed_samples_as_z"
    qc_info["effect_estimate_binary_interpretation"] = (
        (
            "approximate_standardized_not_log_odds_or_liability_scale"
            if reconstruction["scale_name"] == "standardized_effect_estimate"
            else "approximate_phenotype_scaled_not_log_odds_or_liability_scale"
        )
        if trait_type == "binary"
        else None
    )

    _info(
        logger,
        "BETA and SE were reconstructed with effect_from_z.method='{}' "
        "({}): {} and {}. Reference: {}; PubMed: {}.".format(
            reconstruction["method"],
            reconstruction["label"],
            reconstruction["beta_formula"],
            reconstruction["se_formula"],
            reconstruction["citation"],
            ", ".join(reconstruction["pubmed_urls"]),
        ),
    )
    _warn(
        logger,
        "The BETA and SE approximation assumes that EAF was measured in "
        "the same analysed samples as the Z score. If EAF comes from an external "
        "reference panel, ancestry or sample-frequency differences can make the "
        "derived estimates inaccurate.",
    )
    if trait_type == "binary":
        _warn(
            logger,
            "For this binary trait, the approximation uses the case-control effective "
            "sample size. The derived BETA is not a direct reconstruction "
            "of the study's logistic-regression log odds ratio and is not a "
            "liability-scale effect.",
        )

    _info(
        logger,
        "Z to BETA/SE conversion:\n"
        "Variants in: {:,}\n"
        "Removed: {:,} ({:.2f}%)\n"
        "  unusable Z score: {:,}\n"
        "  missing frequency: {:,}\n"
        "  degenerate frequency (exactly 0 or 1): {:,} "
        "(at or below 0: {:,}; at or above 1: {:,})\n"
        "  missing effective sample size: {:,}\n"
        "  non-positive effective sample size: {:,}\n"
        "  below effective-variance floor ({}): {:,}\n"
        "  invalid reconstruction denominator: {:,}\n"
        "Variants out: {:,}".format(
            n_initial,
            n_removed,
            pct_removed,
            n_z_null,
            n_eaf_null,
            n_eaf_invalid,
            n_eaf_lt_zero,
            n_eaf_gt_one,
            n_neff_null,
            n_neff_invalid,
            minimum_effective_variance,
            n_low_effective_variance_removed,
            n_denominator_invalid,
            n_final,
        ),
    )

    # -------------------------------------------------------
    # FINAL SUMMARY
    # -------------------------------------------------------
    summary = df.select(
        [
            pl.col(beta_col).min().alias("beta_min"),
            pl.col(beta_col).max().alias("beta_max"),
            pl.col(beta_col).mean().alias("beta_mean"),
            pl.col(beta_col).std().alias("beta_std"),
            pl.col(se_col).min().alias("se_min"),
            pl.col(se_col).max().alias("se_max"),
            pl.col(se_col).mean().alias("se_mean"),
            pl.col(se_col).std().alias("se_std"),
            pl.len().alias("total"),
        ]
    ).to_dicts()[0]

    qc_info.update(summary)

    _info(
        logger,
        "BETA and SE summary:\n"
        "BETA -> Min: {} | Max: {} | Mean: {} | Std: {}\n"
        "SE   -> Min: {} | Max: {} | Mean: {} | Std: {} | Total: {:,}".format(
            _fmt(summary["beta_min"]),
            _fmt(summary["beta_max"]),
            _fmt(summary["beta_mean"]),
            _fmt(summary["beta_std"]),
            _fmt(summary["se_min"]),
            _fmt(summary["se_max"]),
            _fmt(summary["se_mean"]),
            _fmt(summary["se_std"]),
            summary["total"],
        ),
    )

    return df, qc_info, sample_column_dict
