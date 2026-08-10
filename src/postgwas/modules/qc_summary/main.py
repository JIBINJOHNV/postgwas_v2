"""Build the tabular QC summary for one harmonised GWAS-VCF."""

from pathlib import Path

from postgwas.modules.qc_summary.qc_summary import (
    bcftools_essential_summary,
    parse_bcftools_stats,
    run_bcftools_stats,
)


def run_qc_summary(
    vcf_path: str,
    output_directory: str,
    dataset_id: str,
    external_af_name: str = "EUR",
    allelefreq_diff_cutoff: float = 0.2,
    threads: int = 5,
    bcftools_bin: str = "bcftools",
):
    """Run VCF statistics and return a one-column QC table."""
    af_qc = run_bcftools_stats(
        vcf_path=vcf_path,
        external_af_name=external_af_name,
        allelefreq_diff_cutoff=allelefreq_diff_cutoff,
        threads=threads,
        bcftools_bin=bcftools_bin,
    )

    data_df = bcftools_essential_summary(
        parse_bcftools_stats("%s.stats" % vcf_path)
    ).T
    data_df.columns = ["raw_variants"]

    for metric in (
        "external_af_missing",
        "study_af_missing",
        "af_comparable",
        "af_difference_above_cutoff",
    ):
        data_df.loc[metric] = af_qc[metric]
    data_df = data_df.reset_index()

    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    data_df.to_csv(
        output_directory / ("%s_qc_summary.tsv" % dataset_id),
        sep="\t",
        index=False,
    )
    return data_df
