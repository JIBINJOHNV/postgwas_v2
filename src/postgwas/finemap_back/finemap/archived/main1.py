
from postgwas.finemap.finemap.main import run_finemap_pipeline


def main(args):
    run_finemap_pipeline(
        sig_loci_path=args.sig_loci_path,
        sumstats_path=args.sumstats_path,
        finemap_ld_ref_bim=args.finemap_ld_ref_bim,
        outdir=args.outdir,
        max_workers=args.max_workers,
        n_threads_ldstore=args.n_threads_ldstore,
    )

