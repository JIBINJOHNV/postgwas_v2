from pathlib import Path
from postgwas.modules.ldsc.ldsc_runner import run_ldsc

def run_ldsc_direct(args):
    """
    Wrap your existing run_ldsc_pipeline(sumstats_tsv, out_prefix, ...) function
    so it fits the CLI 'args' style.
    """
    output_directory = Path(args.output_directory).expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    out_prefix = output_directory / args.dataset_id
    run_ldsc(
        sumstats_tsv=str(args.ldsc_input),
        out_prefix=str(out_prefix),
        hm3_snplist=str(args.merge_alleles),
        ldscore_dir=str(args.ref_ld_chr),
        weight_ldscore_dir=str(args.w_ld_chr),
        info_min=args.ldsc_minimum_info,
        maf_min=args.ldsc_minimum_maf,
        samp_prev=args.samp_prev,
        pop_prev=args.pop_prev,
    )
