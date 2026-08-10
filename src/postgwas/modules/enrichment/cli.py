#!/usr/bin/env python3
import argparse
from postgwas.modules.enrichment.arguments import get_geneset_common_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples

# Shared parsers
from postgwas.cli.common import (
    get_common_out_parser,
)

# =========================================================
# SHARED LDSC INPUT ARGUMENTS
# =========================================================

def get_geneset_parser(add_help: bool = False) -> argparse.ArgumentParser:
    """
    Defines input arguments for gene set over representation analysis.
    """
    # 1. Capture the parent parsers
    # These functions likely return parser objects containing the arguments
    geneset_p = get_geneset_common_parser()
    out_p = get_common_out_parser()
    resource_p = get_compute_parser()

    # 2. Initialize your main parser with 'parents'
    parser = argparse.ArgumentParser(
        prog="postgwas pathway_enrichment",
        usage="postgwas pathway_enrichment --gene-input-file PATH [options]",
        description=(
            "Run single-list pathway and interaction enrichment across the configured "
            "external providers."
        ),
        epilog=format_cli_examples(
            (
                "Run pathway enrichment for a gene list:",
                "postgwas pathway_enrichment",
                (
                    "--gene-input-file genes.txt",
                    "--biogrid-key YOUR_BIOGRID_KEY",
                    "--david-email you@example.org",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
        ),
        add_help=add_help,
        formatter_class=AlignedRichHelpFormatter,
        parents=[geneset_p, out_p, resource_p],
    )

    return parser
