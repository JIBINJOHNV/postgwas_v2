

from postgwas.ld_clump.ld_prune_region import ld_clump_by_regions
from postgwas.ld_clump.ld_prune_standard import ld_clump_standard 
import sys
import logging

def run_ld_clump_direct(args, ctx=None):
    """
    Wrapper that dispatches LD clumping functions with robust error handling.
    Exits if both region-based and standard clumping fail.
    """
    if ctx is None:
        ctx = {}

    region_outputs = None
    region_standard = None
    
    # 1. Attempt Region-based Clumping
    try:
        print("\t\t\t\t [*] Attempting Region-based LD clumping...")
        region_outputs = ld_clump_by_regions(
            sumstat_vcf=args.vcf,
            output_folder=args.outdir,
            sample_name=args.sample_id,
            population=args.ld_clump_population,
            nthreads=args.nthreads,
            bcftools=args.bcftools
        )
    except Exception as e:
        print(f"\t\t\t\t [!] Warning: ld_clump_by_regions failed: {e}")

    # 2. Attempt Standard LD Clumping
    try:
        print("\t\t\t\t [*] Attempting Standard LD clumping...")
        region_standard = ld_clump_standard(
            vcf_path=args.vcf,
            outdir=args.outdir,
            sample_id=args.sample_id,
            r2_clump=args.r2_clump,
            r2_lead=args.r2_lead,
            lead_p=args.lead_p,
            merge_dist=args.merge_dist,
            n_threads=args.nthreads,
            ld_folder=args.ld_folder,
            pop=args.ld_clump_population,
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