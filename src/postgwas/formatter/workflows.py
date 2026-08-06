#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from rich_argparse import RichHelpFormatter
from postgwas.utils.main import run_cmd,require_executable


# ---- Formatter engines ----
from postgwas.formatter.main import (
    create_magma_inputs,
    create_finemap_susie_inputs,
    create_finemap_finemap_inputs,
    create_ldpred_inputs,
    create_ldsc_inputs
)



# ============================================================
#  DIRECT MODE ENGINE
# ============================================================
from pathlib import Path
import sys

def run_formatter_direct(args, ctx=None):
    """
    Formatter engine.
    Works in both direct mode and pipeline mode.
    Returns structured downstream inputs.
    """
    # -------------------------------------------------
    # Dependency checks (FAIL FAST)
    # -------------------------------------------------
    selected_formats = args.format  # nargs="+"

    format_map = {
        "magma": create_magma_inputs,
        "finemap_susie": create_finemap_susie_inputs,
        "finemap_finemap": create_finemap_finemap_inputs,
        "ldpred": create_ldpred_inputs,
        "ldsc": create_ldsc_inputs,
    }

    outputs = {
        "magma": {},
        "finemap_susie": {},
        "finemap_finemap": {},
        "ldsc": {},
        "ldpred": {},
    }

    try:
        for fmt in selected_formats:
            if fmt not in format_map:
                raise ValueError(f"Unknown format: {fmt}")

            # -------------------------------------------------
            # Each formatter returns a dict → capture it
            # -------------------------------------------------
            result = format_map[fmt](args)

            if not isinstance(result, dict):
                raise ValueError(
                    f"Formatter '{fmt}' must return a dict, got {type(result)}"
                )

            outputs[fmt] = result

    except Exception as e:
        print(f"\n❌ Formatter failed: {e}", file=sys.stderr)
        sys.exit(1)

    # -------------------------------------------------
    # Pipeline mode: register outputs
    # -------------------------------------------------
    if ctx is not None:
        ctx["formatter"] = outputs
    return outputs


