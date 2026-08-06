import argparse
import sys
from rich_argparse import RichHelpFormatter 
# Assuming utility/logic imports
from postgwas.utils.main import validate_path 
from postgwas.annot_ldblock.annot_ldblock import annotate_ldblocks 
from pathlib import Path
from postgwas.utils.main import require_executable




# ---------------------------------------------------------
# 2. EXECUTABLE LOGIC (Engines)
# ---------------------------------------------------------
def run_annot_ldblock(args, ctx=None):
    """
    Executes LD-block annotation.
    Works in both direct mode and pipeline mode.
    """
    
    # -------------------------------------------------
    # Validation
    # -------------------------------------------------
    if not args.vcf or not args.ld_region_dir:
        raise ValueError("annot_ldblock: --vcf and --ld-dir arguments are required.")

    if not args.outdir:
        raise ValueError("annot_ldblock: --outdir is required for pipeline execution.")

    # -------------------------------------------------
    # Run engine
    # -------------------------------------------------
    outputs = annotate_ldblocks(
        vcf_path=args.vcf,
        outdir=args.outdir,
        genome_version=args.genome_version,
        ld_dir=args.ld_region_dir,
        populations=tuple(args.ld_block_population),
        n_threads=args.nthreads,
        max_memory_gb=args.max_mem,
        smple_id=args.sample_id
    )

    # -------------------------------------------------
    # Pipeline mode: register outputs
    # -------------------------------------------------
    if ctx is not None:
        ctx["annot_ldblock"] = outputs
    # -------------------------------------------------
    # Direct mode: still succeed cleanly
    # -------------------------------------------------
    return outputs