from postgwas.qc_summary.main import run_qc_summary



def run_qc_summary_direct(args):
    run_qc_summary(
        vcf_path=args.vcf,
        qc_outdir=args.outdir,
        sample_id=args.sample_id,
        external_af_name=args.qc_summary_external_af_name,
        allelefreq_diff_cutoff=args.qc_summary_allelefreq_diff_cutoff,
        n_threads=args.nthreads,
        bcftools_bin=args.bcftools,
    )



