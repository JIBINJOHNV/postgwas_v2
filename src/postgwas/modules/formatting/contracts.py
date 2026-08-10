"""Scientifically reviewed contracts for downstream formatter outputs."""

from dataclasses import dataclass


@dataclass(frozen=True)
class FormatContract:
    """Document the columns and scientific meaning of one exported format."""

    name: str
    frequency: str
    sample_size: str
    source: str


FORMAT_CONTRACTS = {
    "gcta_gene": FormatContract(
        name="GCTA COJO / fastBAT / mBAT-combo",
        frequency="freq is effect-allele frequency for the shared GCTA .ma schema",
        sample_size=(
            "N is total sample size in the shared GCTA .ma schema used by "
            "COJO, fastBAT, and mBAT-combo"
        ),
        source="GCTA documentation: https://yanglab.westlake.edu.cn/software/gcta/",
    ),
    "magma": FormatContract(
        name="MAGMA",
        frequency="not required",
        sample_size="N_COL is total sample size, including case-control studies",
        source=(
            "MAGMA manual v1.09a: https://ibg.colorado.edu/cdrom2021/"
            "Day10-posthuma/magma_session/manual_v1.09a.pdf"
        ),
    ),
    "finemap": FormatContract(
        name="FINEMAP",
        frequency="maf = min(effect-allele frequency, 1 − effect-allele frequency)",
        sample_size="NEF supplies PostGWAS n_samples for each locus",
        source=(
            "FINEMAP v1.4 documentation: https://christianbenner.com/; "
            "PostGWAS FINEMAP adapter"
        ),
    ),
    "susie": FormatContract(
        name="SuSiE-RSS",
        frequency="not used by the configured PostGWAS SuSiE-RSS call",
        sample_size="NEF supplies the n argument",
        source=(
            "susieR susie_rss documentation: "
            "https://stephenslab.github.io/susieR/reference/susie_rss.html; "
            "PostGWAS SuSiE engine"
        ),
    ),
    "pred_ld": FormatContract(
        name="PRED-LD",
        frequency="not consumed by PRED-LD; AF is retained for re-harmonisation",
        sample_size="not consumed by PRED-LD; NC/SS are retained for re-harmonisation",
        source=(
            "PRED-LD documentation: https://github.com/pbagos/PRED-LD; "
            "PostGWAS imputation handoff"
        ),
    ),
    "ldsc": FormatContract(
        name="LDSC",
        frequency="FRQ is effect-allele frequency; munge_sumstats derives MAF for QC",
        sample_size=(
            "binary: NC/NCO case/control counts; quantitative: N is taken from NCO"
        ),
        source=(
            "CBIIT LDSC munge_sumstats.py: "
            "https://github.com/CBIIT/ldsc/blob/ldsc39/munge_sumstats.py"
        ),
    ),
    "mixer": FormatContract(
        name="MiXeR (univariate)",
        frequency="not required by fit1/test1",
        sample_size=(
            "binary: effective N from NEF = 4/(1/Ncase + 1/Ncontrol); "
            "quantitative: total N from NCO"
        ),
        source="MiXeR user documentation: https://github.com/precimed/mixer",
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
GCTA_SAMPLE_SIZE_MODE = "total_sample_size_from_FORMAT_SS"


def required_formats(config, modules=(), requested=()):
    """Return downstream formats once, in stable execution order."""
    selected = set(requested or ())
    for module in modules or ():
        selected.update(config.module_formats.get(module, ()))
    return [name for name in config.format_order if name in selected]


__all__ = [
    "FORMAT_CONTRACTS",
    "GCTA_SAMPLE_SIZE_MODE",
    "NAMED_OUTPUT_RESULT_KEYS",
    "FormatContract",
    "required_formats",
]
