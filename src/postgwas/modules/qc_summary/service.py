from postgwas.modules.qc_summary.main import run_qc_summary
from postgwas.cli.compute import resolve_compute_args



def run_qc_summary_direct(args):
    resolve_compute_args(args)
    run_qc_summary(
        vcf_path=args.vcf,
        output_directory=args.output_directory,
        dataset_id=args.dataset_id,
        external_af_name=args.reference_af_column,
        allelefreq_diff_cutoff=args.maximum_af_difference,
        threads=args.threads,
        bcftools_bin=args.bcftools,
    )
