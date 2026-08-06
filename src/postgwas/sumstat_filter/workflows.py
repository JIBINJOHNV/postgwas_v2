
from postgwas.sumstat_filter.sumstat_filter import filter_gwas_vcf_bcftools
from pathlib import Path
from postgwas.utils.main import require_executable


def run_sumstat_filter_direct(args,ctx=None):
    # Your logic — YOU said do not enforce checks here
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    require_executable("bcftools")
    require_executable("tabix")

    # -------------------------------------------------
    # Validation
    # -------------------------------------------------
    if not args.vcf :
        raise ValueError("run_sumstat_filter_direct: --vcf arguments are required.")

    if not args.outdir:
        raise ValueError("run_sumstat_filter_direct: --outdir is required for execution.")

    outputs=filter_gwas_vcf_bcftools(
        vcf_path=args.vcf,
        output_folder=str(outdir),
        output_prefix=args.sample_id,
        pval_cutoff=args.pval_cutoff,
        maf_cutoff=args.maf_cutoff,
        allelefreq_diff_cutoff=args.allelefreq_diff_cutoff,
        info_cutoff=args.info_cutoff,
        info_missing=args.info_missing,
        info_min=args.info_min,
        info_max=args.info_max,
        external_af_name=args.external_af_name,
        include_indels=args.include_indels,
        exclude_palindromic=args.exclude_palindromic,
        remove_mhc=args.remove_mhc,
        mhc_chrom=args.mhc_chrom,
        mhc_start=args.mhc_start,
        mhc_end=args.mhc_end,
        threads=args.nthreads,
        max_mem=args.max_mem,

    )
    # -------------------------------------------------
    # Pipeline mode: register outputs
    # -------------------------------------------------
    if ctx is not None:
        ctx["sumstat_filter"] = outputs
    return outputs
