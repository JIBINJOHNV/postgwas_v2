#!/usr/bin/env python3
"""
Unified CLI for the postgwas toolkit.
"""
import argparse
import sys
from rich_argparse import RichHelpFormatter


def main():
    # -----------------------------------------------------------
    # 1. Top-Level Parser
    # -----------------------------------------------------------
    parser = argparse.ArgumentParser(
        prog="postgwas",
        description="PostGWAS — Unified Post-GWAS Analysis Toolkit",
        formatter_class=RichHelpFormatter,
        epilog=(
            "Examples:\n"
            "  postgwas finemap pipeline --help\n"
            "  postgwas harmonise --help\n"
            "  postgwas qc --help\n"
        ),
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="COMMAND",
        title="Available subcommands",
        help="Run `postgwas COMMAND --help` for module-specific options.",
    )

    # -----------------------------------------------------------
    # Subcommands (alphabetical by COMMAND)
    # -----------------------------------------------------------

    # 1. ANNOT LD BLOCK (src/postgwas/annot_ldblock/cli.py)
    try:
        from postgwas.annot_ldblock import cli as annot_cli
        sub = subparsers.add_parser(
            "annot-ld",
            help="Annotate variants with LD-block information",
            add_help=False,
        )
        sub.set_defaults(module=annot_cli)
    except ImportError:
        pass

    # 2. FINEMAP (src/postgwas/finemap/cli.py)
    try:
        from postgwas.finemap import cli as finemap_cli
        sub = subparsers.add_parser(
            "finemap",
            help="Fine-mapping (SuSiE / FINEMAP)",
            add_help=False,
        )
        sub.set_defaults(module=finemap_cli)
    except ImportError:
        pass

    # 3. HARMONISATION (src/postgwas/harmonisation/cli.py)
    try:
        from postgwas.harmonisation import cli as harmonise_cli
        sub = subparsers.add_parser(
            "harmonise",
            help="Harmonise GWAS summary statistics",
            add_help=False,
        )
        sub.set_defaults(module=harmonise_cli)
    except ImportError:
        pass

    # 4. IMPUTATION (src/postgwas/imputation/cli.py)
    try:
        from postgwas.imputation import cli as impute_cli
        sub = subparsers.add_parser(
            "impute",
            help="Impute z-scores or missing summary statistics",
            add_help=False,
        )
        sub.set_defaults(module=impute_cli)
    except ImportError:
        pass

    # 5. MANHATTAN (src/postgwas/manhattan/cli.py)
    try:
        from postgwas.manhattan import cli as manhattan_cli
        sub = subparsers.add_parser(
            "manhattan",
            help="Generate Manhattan plots",
            add_help=False,
        )
        sub.set_defaults(module=manhattan_cli)
    except ImportError:
        pass

    # 6. LD PRUNE (src/postgwas/ld_prune/cli.py)
    try:
        from postgwas.ld_prune import cli as prune_cli
        sub = subparsers.add_parser(
            "prune",
            help="LD pruning / clumping utilities",
            add_help=False,
        )
        sub.set_defaults(module=prune_cli)
    except ImportError:
        pass

    # 7. QC SUMMARY (src/postgwas/qc_summary/cli.py)
    try:
        from postgwas.qc_summary import cli as qc_cli
        sub = subparsers.add_parser(
            "qc",
            help="GWAS QC summary and diagnostics",
            add_help=False,
        )
        sub.set_defaults(module=qc_cli)
    except ImportError:
        pass

    # 8. SUMSTAT FILTER (src/postgwas/sumstat_qc/cli.py)
    try:
        from postgwas.sumstat_qc import cli as sumstat_cli
        sub = subparsers.add_parser(
            "sumstat-filter",
            help="Filter GWAS summary statistics",
            add_help=False,
        )
        sub.set_defaults(module=sumstat_cli)
    except ImportError:
        pass

    # 9. TRANSFMT / FORMATTER (src/postgwas/formatter/cli.py)
    try:
        from postgwas.formatter import cli as transfmt_cli
        sub = subparsers.add_parser(
            "transfmt",
            help="Transform / format files for downstream tools",
            add_help=False,
        )
        sub.set_defaults(module=transfmt_cli)
    except ImportError:
        pass

    # -----------------------------------------------------------
    # Execution Logic
    # -----------------------------------------------------------
    # parse_known_args stops parsing when it hits the subcommand,
    # leaving all subsequent flags in 'remaining'.
    args, remaining = parser.parse_known_args()

    # If no subcommand is provided, print general help
    if args.command is None:
        parser.print_help()
        sys.exit(0)

    # Reconstruct sys.argv to simulate running the submodule directly
    # e.g., ["postgwas finemap", "--locus_file", "x", "--verbose"]
    sys.argv = [f"postgwas {args.command}"] + remaining

    # Pass control to the specific module's main() function
    try:
        return args.module.main()
    except AttributeError:
        print(
            f"Error: The module for '{args.command}' does not have a main() function.",
            file=sys.stderr,
        )
        sys.exit(1)
    except Exception as e:
        print(f"Error executing {args.command}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
