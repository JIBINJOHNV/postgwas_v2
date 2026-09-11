"""Scientifically reviewed contracts for downstream formatter outputs."""

from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True)
class FormatContract:
    """Document the columns and scientific meaning of one exported format."""

    name: str
    frequency: str
    source: str


FORMAT_CONTRACTS = {
    "gcta_gene": FormatContract(
        name="GCTA COJO / fastBAT / mBAT-combo",
        frequency="freq is effect-allele frequency for the shared GCTA .ma schema",
        source="GCTA documentation: https://yanglab.westlake.edu.cn/software/gcta/",
    ),
    "magma": FormatContract(
        name="MAGMA",
        frequency="not required",
        source=(
            "MAGMA manual v1.09a: https://ibg.colorado.edu/cdrom2021/"
            "Day10-posthuma/magma_session/manual_v1.09a.pdf"
        ),
    ),
    "finemap": FormatContract(
        name="FINEMAP",
        frequency="maf = min(effect-allele frequency, 1 − effect-allele frequency)",
        source=(
            "FINEMAP v1.4 documentation: https://christianbenner.com/; "
            "PostGWAS FINEMAP adapter"
        ),
    ),
    "susie": FormatContract(
        name="SuSiE-RSS",
        frequency="not used by the configured PostGWAS SuSiE-RSS call",
        source=(
            "susieR susie_rss documentation: "
            "https://stephenslab.github.io/susieR/reference/susie_rss.html; "
            "PostGWAS SuSiE engine"
        ),
    ),
    "pred_ld": FormatContract(
        name="PRED-LD",
        frequency="not consumed by PRED-LD; AF is retained for re-harmonisation",
        source=(
            "PRED-LD documentation: https://github.com/pbagos/PRED-LD; "
            "PostGWAS imputation handoff"
        ),
    ),
    "ldsc": FormatContract(
        name="LDSC",
        frequency="FRQ is effect-allele frequency; munge_sumstats derives MAF for QC",
        source=(
            "CBIIT LDSC munge_sumstats.py: "
            "https://github.com/CBIIT/ldsc/blob/ldsc39/munge_sumstats.py"
        ),
    ),
    "mixer": FormatContract(
        name="MiXeR (univariate)",
        frequency="not required by fit1/test1",
        source="MiXeR user documentation: https://github.com/precimed/mixer",
    ),
}

CUSTOM_OUTPUT_TARGET = "custom"

DUPLICATE_POLICY_DESCRIPTIONS = {
    "exclude_all": (
        "exclude every conflicting record in a duplicated identifier group"
    ),
    "error": (
        "collapse exact repeats and stop if a conflicting duplicated-ID group remains"
    ),
    "lowest_p": (
        "retain the lowest-p-value record in each duplicated identifier group"
    ),
    "most_significant": (
        "retain only a unique most-significant record; exclude tied or missing ranks"
    ),
    "highest_maf": (
        "retain only a unique highest-MAF record; exclude tied or missing ranks"
    ),
    "highest_info": (
        "retain only a unique highest-INFO record; exclude tied or missing ranks"
    ),
}

# Stable Python result keys are internal handoff contracts, not configurable
# scientific values. Keeping them here prevents exporters and resume recovery
# from defining parallel interfaces.
NAMED_OUTPUT_RESULT_KEYS = {
    "gcta_gene": {
        "summary_statistics": "summary_statistics_input_file",
    },
    "magma": {
        "snp_location": "snp_loc_file",
        "p_values": "pval_file",
    },
}
SINGLE_OUTPUT_RESULT_KEYS = {
    "susie": "susie_input",
    "finemap": "finemap_input",
    "ldsc": "ldsc_file",
    "mixer": "mixer_input",
    CUSTOM_OUTPUT_TARGET: "custom_output_file",
}
PARTITIONED_OUTPUT_RESULT_KEYS = {
    "pred_ld": {
        "directory": "pred_ld_folder",
        "files": "files",
    },
}
GCTA_SAMPLE_SIZE_MODE = "total_sample_size_from_FORMAT_SS"


def formatter_result_targets(config, selected: list[str]) -> list[str]:
    """Return built-in targets plus the optional additive custom export."""
    targets = list(selected)
    if config.custom_output.active:
        targets.append(CUSTOM_OUTPUT_TARGET)
    return targets


def formatter_target_display_name(
    target: str,
    variant_id_observations: Mapping[str, Mapping] | None = None,
) -> str:
    """Return the resolved downstream consumer label for one formatter target."""
    if target == CUSTOM_OUTPUT_TARGET:
        return "Custom CLI table"
    observation = (variant_id_observations or {}).get(target, {})
    consumers = [
        str(consumer).strip()
        for consumer in observation.get("consumers", ())
        if str(consumer).strip()
    ]
    return " / ".join(consumers) if consumers else FORMAT_CONTRACTS[target].name


def required_formats(config, modules=(), requested=()):
    """Return downstream formats once, in stable execution order."""
    selected = set(requested or ())
    for module in modules or ():
        selected.update(config.module_formats.get(module, ()))
    return [name for name in config.format_order if name in selected]


__all__ = [
    "FORMAT_CONTRACTS",
    "CUSTOM_OUTPUT_TARGET",
    "DUPLICATE_POLICY_DESCRIPTIONS",
    "GCTA_SAMPLE_SIZE_MODE",
    "NAMED_OUTPUT_RESULT_KEYS",
    "PARTITIONED_OUTPUT_RESULT_KEYS",
    "SINGLE_OUTPUT_RESULT_KEYS",
    "FormatContract",
    "formatter_result_targets",
    "formatter_target_display_name",
    "required_formats",
]
