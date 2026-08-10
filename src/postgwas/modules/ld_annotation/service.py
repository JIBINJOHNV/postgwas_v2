import argparse
import sys
from rich_argparse import RichHelpFormatter
# Assuming utility/logic imports
from postgwas.core.execution.runtime import validate_path
from postgwas.modules.ld_annotation.annot_ldblock import annotate_ldblocks
from pathlib import Path
from postgwas.core.execution.runtime import require_executable
from postgwas.cli.compute import memory_limit, resolve_compute_args




# ---------------------------------------------------------
# 2. EXECUTABLE LOGIC (Engines)
# ---------------------------------------------------------
def run_annot_ldblock(args, ctx=None):
    """
    Executes LD-block annotation.
    Works in both direct mode and pipeline mode.
    """

    resolve_compute_args(args)
    # -------------------------------------------------
    # Validation
    # -------------------------------------------------
    if not args.vcf or not args.ld_region_dir:
        raise ValueError("annot_ldblock: --vcf and --ld-dir arguments are required.")

    if not args.output_directory:
        raise ValueError(
            "annot_ldblock: --output-directory is required for pipeline execution."
        )

    # -------------------------------------------------
    # Run engine
    # -------------------------------------------------
    outputs = annotate_ldblocks(
        vcf_path=args.vcf,
        output_directory=args.output_directory,
        genome_build=args.genome_build,
        ld_dir=args.ld_region_dir,
        populations=tuple(args.ld_block_populations),
        threads=args.threads,
        max_memory_gb=memory_limit(args.memory_gb),
        dataset_id=args.dataset_id,
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
