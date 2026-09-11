"""Shared scientific-policy contract for QC assessment and VCF filtering."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import isfinite
from typing import Any, Mapping


_MISSING_ACTIONS = frozenset(("keep", "remove"))

# MaCH Rsq is conventionally supported through 2; both QC backends share this
# schema/runtime ceiling while retaining independently configured defaults.
MAXIMUM_SUPPORTED_INFO_SCORE = 2.0


def variant_qc_rule_display_groups(
    fields: Mapping[str, str],
) -> dict[str, dict[str, str]]:
    """Build common QC-rule headings from resolved VCF field names."""
    required = (
        "chromosome", "position", "reference_allele", "alternate_allele",
        "study_af", "imputation_quality", "log_pvalue", "study_info_af",
        "external_info_af",
    )
    missing = [key for key in required if not str(fields.get(key, "")).strip()]
    if missing:
        raise ValueError(
            "Variant-QC display fields are missing: %s" % ", ".join(missing)
        )

    def group(label: str, *values: str, kind: str) -> dict[str, str]:
        context = " and ".join(str(value) for value in values if value)
        return {
            "display_group": (
                "%s · %s" % (label, context) if context else label
            ),
            "display_group_kind": kind,
        }

    allele_fields = "%s/%s" % (
        fields["reference_allele"], fields["alternate_allele"],
    )
    variant_type = str(fields.get("variant_type") or allele_fields)
    return {
        "filtering_policy": group(
            "Filtering policy", kind="decision",
        ),
        "significance": group(
            "Statistical significance", fields["log_pvalue"], kind="analysis",
        ),
        "allele_frequency": group(
            "Allele frequency", fields["study_af"], kind="genetic",
        ),
        "imputation_quality": group(
            "Imputation quality", fields["imputation_quality"], kind="analysis",
        ),
        "external_af_concordance": group(
            "Study/reference frequency concordance",
            fields["study_info_af"],
            fields["external_info_af"],
            kind="genetic",
        ),
        "variant_type": group(
            "Variant type", variant_type, kind="genetic",
        ),
        "palindromic_variants": group(
            "Palindromic allele ambiguity",
            allele_fields,
            fields["study_af"],
            kind="genetic",
        ),
        "mhc_region": group(
            "Genomic region",
            "%s/%s" % (fields["chromosome"], fields["position"]),
            kind="genetic",
        ),
    }


@dataclass(frozen=True)
class VariantQCPolicy:
    """Normalised values consumed by the Polars and bcftools QC backends."""

    minimum_neglog10_p: float | None
    missing_pvalue_action: str
    maf_min: float | None
    missing_af_action: str
    info_min: float | None
    info_max: float | None
    missing_info_action: str
    maximum_af_difference: float | None
    include_indels: bool
    remove_palindromic: bool
    palindromic_lower: float
    palindromic_upper: float
    remove_mhc: bool
    mhc_chromosome: str | None = None
    mhc_start: int | None = None
    mhc_end: int | None = None

    def __post_init__(self) -> None:
        for name in (
            "missing_pvalue_action", "missing_af_action", "missing_info_action",
        ):
            if getattr(self, name) not in _MISSING_ACTIONS:
                raise ValueError("%s must be either 'keep' or 'remove'" % name)
        self._bounded("minimum_neglog10_p", self.minimum_neglog10_p, 0.0, None)
        self._bounded("maf_min", self.maf_min, 0.0, 0.5)
        self._bounded("info_min", self.info_min, 0.0, 1.0)
        self._bounded(
            "info_max", self.info_max, 0.0, MAXIMUM_SUPPORTED_INFO_SCORE,
        )
        self._bounded(
            "maximum_af_difference", self.maximum_af_difference, 0.0, 1.0,
        )
        self._bounded("palindromic_lower", self.palindromic_lower, 0.0, 0.5)
        self._bounded("palindromic_upper", self.palindromic_upper, 0.5, 1.0)
        if (
            self.info_min is not None
            and self.info_max is not None
            and self.info_min > self.info_max
        ):
            raise ValueError("info_max must be greater than or equal to info_min")
        if self.palindromic_upper <= self.palindromic_lower:
            raise ValueError(
                "palindromic_upper must be greater than palindromic_lower"
            )
        if self.remove_mhc:
            if not str(self.mhc_chromosome or "").strip():
                raise ValueError("MHC removal requires a chromosome")
            if self.mhc_start is None or self.mhc_end is None:
                raise ValueError("MHC removal requires start and end coordinates")
            if int(self.mhc_start) < 1 or int(self.mhc_end) < int(self.mhc_start):
                raise ValueError("MHC coordinates must be a valid inclusive interval")

    @staticmethod
    def _bounded(
        name: str,
        value: float | None,
        minimum: float,
        maximum: float | None,
    ) -> None:
        if value is None:
            return
        numeric = float(value)
        if not isfinite(numeric) or numeric < minimum:
            raise ValueError("%s is outside its supported range" % name)
        if maximum is not None and numeric > maximum:
            raise ValueError("%s is outside its supported range" % name)

    def active_rule_keys(self) -> tuple[str, ...]:
        """Return the common scientific rules active in either backend."""
        keys = []
        if self.minimum_neglog10_p is not None:
            keys.append("significance")
        if self.maf_min is not None:
            keys.append("allele_frequency")
        if self.info_min is not None or self.info_max is not None:
            keys.append("imputation_quality")
        if self.maximum_af_difference is not None:
            keys.append("external_af_concordance")
        if not self.include_indels:
            keys.append("variant_type")
        if self.remove_palindromic:
            keys.append("palindromic_variants")
        if self.remove_mhc:
            keys.append("mhc_region")
        return tuple(keys)

    def as_dict(self) -> dict[str, Any]:
        values = asdict(self)
        values["active_rule_keys"] = list(self.active_rule_keys())
        return values
