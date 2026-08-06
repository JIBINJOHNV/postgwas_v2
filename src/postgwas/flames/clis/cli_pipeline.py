# postgwas/flames/clis/cli_pipeline.py

from __future__ import annotations
import argparse
from rich.console import Console
from rich_argparse import ArgumentDefaultsRichHelpFormatter

# import submodules
from .pipeline.cli_gwas import add_gwas_arguments
from .pipeline.cli_magma import add_magma_arguments
from .pipeline.cli_susie import add_susie_arguments
from .pipeline.cli_pops import add_pops_arguments
from .pipeline.cli_runtime import add_runtime_arguments

from postgwas.flames.clis.pipeline.cli_magmacovar import add_magma_covar_arguments


def run_flames_pipeline(args):
    console = Console()
    console.rule("[bold magenta]FLAMES — Pipeline Mode[/bold magenta]")
    console.print(vars(args))


def add_pipeline_mode(subparsers, parent):
    p = subparsers.add_parser(
        "pipeline",
        help="Full FLAMES workflow: VCF → MAGMA → SuSiE → PoPS → ML.",
        description=(
            "[bold cyan]PIPELINE MODE[/bold cyan]\n\n"
            "Runs:\n"
            "  1) MAGMA gene-level statistics\n"
            "  2) SuSiE Bayesian fine-mapping\n"
            "  3) PoPS gene prioritisation\n"
            "  4) FLAMES ML effector-gene modelling"
        ),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
        add_help=False,
        parents=[parent],
    )

    # custom help
    p.add_argument("-h", "--help", action="store_true", help="Show help and exit.")

    # add modular groups
    add_gwas_arguments(p)
    add_magma_arguments(p)
    add_susie_arguments(p)
    add_pops_arguments(p)
    add_runtime_arguments(p)
    add_magma_covar_arguments(p)

    p.set_defaults(func=run_flames_pipeline)
    return p
