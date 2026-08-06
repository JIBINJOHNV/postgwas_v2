#!/usr/bin/env python3
import argparse
from postgwas.utils.main import validate_path
from rich_argparse import RichHelpFormatter
from postgwas.clis.geneset_cli import get_geneset_common_parser

# Shared parsers
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
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
    resource_p = get_defaultresourse_parser()

    # 2. Initialize your main parser with 'parents'
    parser = argparse.ArgumentParser(
        add_help=add_help,
        formatter_class=RichHelpFormatter,
        parents=[geneset_p, out_p, resource_p] 
    )

    return parser
