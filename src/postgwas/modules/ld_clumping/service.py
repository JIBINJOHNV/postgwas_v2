

from postgwas.modules.ld_clumping.ld_prune_region import ld_clump_by_regions
from postgwas.modules.ld_clumping.ld_prune_standard import ld_clump_standard
import sys
import logging
from postgwas.cli.compute import resolve_compute_args

def run_ld_clump_direct(args, ctx=None):
    """
    Wrapper that dispatches LD clumping functions with robust error handling.
    Exits if both region-based and standard clumping fail.
    """
    resolve_compute_args(args)
    if ctx is None:
        ctx = {}

    region_outputs = None
    region_standard = None

    # 1. Attempt Region-based Clumping
    try:
        print("\t\t\t\t [*] Attempting Region-based LD clumping...")
        region_outputs = ld_clump_by_regions(
            sumstat_vcf=args.vcf,
            output_directory=args.output_directory,
            dataset_id=args.dataset_id,
            population=args.population,
            threads=args.threads,
            bcftools=args.bcftools
        )
    except Exception as e:
        print(f"\t\t\t\t [!] Warning: ld_clump_by_regions failed: {e}")

    # 2. Attempt Standard LD Clumping
    try:
        print("\t\t\t\t [*] Attempting Standard LD clumping...")
        region_standard = ld_clump_standard(
            vcf_path=args.vcf,
            output_directory=args.output_directory,
            dataset_id=args.dataset_id,
            r2_clump=args.r2_clump,
            r2_lead=args.r2_lead,
            lead_p=args.lead_p,
            merge_dist=args.merge_dist,
            threads=args.threads,
            ld_folder=args.ld_folder,
            pop=args.population,
            bcftools_bin=args.bcftools
        )
    except Exception as e:
        print(f"\t\t\t\t [!] Warning: ld_clump_standard failed: {e}")

    # 3. Final Validation: Did both fail?
    if region_outputs is None and region_standard is None:
        print("\n" + "="*50)
        print("CRITICAL FAILURE: Both LD clumping methods failed.")
        print("Please check your input VCF, LD folder paths, and logs.")
        print("="*50)
        sys.exit(1) # Exit with error code 1

    # Update context with whatever results were obtained
    ctx["ld_clump_region"] = region_outputs
    ctx["ld_clump_standard"] = region_standard

    return ctx
