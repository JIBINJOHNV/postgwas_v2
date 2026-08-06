import sys,os
from pathlib import Path
import shutil
import traceback

from postgwas.formatter.to_ldpred import vcf_to_ldpred


def create_magma_inputs(args):
    """Handler for MAGMA conversion."""
    try:
        from postgwas.formatter.to_magma import vcf_to_magma
        magma_inputs=vcf_to_magma(
            sumstat_vcf=args.vcf,
            output_folder=args.outdir,
            sample_name=args.sample_id
        )
    except ImportError:
        print("❌ Error: postgwas.formatter.to_magma module not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error during MAGMA conversion: {e}", file=sys.stderr)
        sys.exit(1)
    return magma_inputs

def create_finemap_susie_inputs(args):
    """Handler for FINEMAP conversion (Placeholder)."""
    try:
        from postgwas.formatter.to_finemap_susie import vcf_to_finemap_susie
        finemap_inputs=vcf_to_finemap_susie(
            sumstat_vcf=args.vcf, 
            output_folder=args.outdir,
            sample_name=args.sample_id
            )
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
    return(finemap_inputs)

def create_finemap_finemap_inputs(args):
    """Handler for FINEMAP conversion (Placeholder)."""
    try:
        from postgwas.formatter.to_finemap_finemap import vcf_to_finemap_finemap
        finemap_inputs=vcf_to_finemap_finemap(
            sumstat_vcf=args.vcf, 
            output_folder=args.outdir,
            sample_name=args.sample_id
            )
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
    return(finemap_inputs)


def create_ldpred_inputs(args):
    """
    Wrapper to call vcf_to_ldpred from CLI args.
    """
    # 1. Resolve BCFtools path
    bcf_cmd = getattr(args, "bcftools", None)
    
    # Fallback/Safety Check
    if not bcf_cmd:
        bcf_cmd = shutil.which("bcftools")
        if not bcf_cmd:
             # Last resort fallback if shutil fails but we assume it's in path
             bcf_cmd = "bcftools"

    try:
        # 2. Call the logic function
        return vcf_to_ldpred(
            sumstat_vcf=args.vcf, 
            output_folder=args.outdir,
            sample_name=args.sample_id,
            bcftools_path=bcf_cmd,
            nthreads=args.nthreads
        )

    except Exception as e:
        print("\n❌ Error in create_ldpred_inputs:")
        traceback.print_exc() # This will show you EXACTLY where it failed
        raise e



def create_ldsc_inputs(args):
    """Handler for ldsc conversion (Placeholder)."""
    try:
        from postgwas.formatter.to_ldsc import vcf_to_ldsc
        ldsc_inputs=vcf_to_ldsc(
            sumstat_vcf=args.vcf, 
            output_folder=args.outdir,
            sample_name=args.sample_id
            )
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
    return(ldsc_inputs)