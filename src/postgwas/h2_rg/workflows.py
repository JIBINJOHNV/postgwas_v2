from pathlib import Path
from postgwas.h2_rg.ldsc_runner import run_ldsc

def run_ldsc_direct(args):
    """
    Wrap your existing run_ldsc_pipeline(sumstats_tsv, out_prefix, ...) function
    so it fits the CLI 'args' style.
    """
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    out_prefix = Path(args.outdir) / args.sample_id  # depending on get_common_out_parser impl
    if args.heritability_tool == "ldsc":
        run_ldsc(
            sumstats_tsv=str(args.ldsc_inut),
            out_prefix=str(out_prefix),
            hm3_snplist=str(args.merge_alleles),
            ldscore_dir=str(args.ref_ld_chr),
            info_min=args.heritability_info_min,
            maf_min=args.heritability_maf_min,
            samp_prev=args.samp_prev,
            pop_prev=args.pop_prev,
        )
    else:
        print(f"The requested heritability_tool {args.heritability_tool} not implimented yet")

