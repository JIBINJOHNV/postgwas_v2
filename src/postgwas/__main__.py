#!/usr/bin/env python3
"""
postgwas — unified top-level CLI for all PostGWAS modules.

Example usage:
  postgwas finemap --help
  postgwas harmonisation pipeline ...
  postgwas heritability direct --help
"""
from postgwas._mp_fix import *
import sys
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text

# Import module entrypoints
from postgwas.harmonisation.cli import main as harmonisation_main
from postgwas.finemap.cli import main as finemap_main
from postgwas.flames.cli import main as flames_main
from postgwas.pops.cli import main as pops_main
from postgwas.qc_summary.cli import main as qc_main
from postgwas.gene_assoc.cli import main as magma_main
from postgwas.magmacovar.cli import main as magmacovar_main
from postgwas.annot_ldblock.cli import main as annot_ldblock_main
from postgwas.ld_clump.cli import main as ld_clump_main
from postgwas.formatter.cli import main as formatter_main
from postgwas.sumstat_filter.cli import main as sumstat_filter_main
from postgwas.manhattan.cli import main as manhattan_main
from postgwas.h2_rg.cli import main as heritability_main
from postgwas.imputation.cli import main as imputation_main
from postgwas.pipeline.cli import main as pipeline_main
from postgwas.enrichment.main import main as enrichment_main
# -------------------------------------------------------------------
# MODULE REGISTRY (alphabetical by key)
# -------------------------------------------------------------------
MODULES = {
    "annot_ldblock": (annot_ldblock_main, "Run LD-block annotation module"),
    "finemap":       (finemap_main,       "Run fine-mapping pipeline (SuSiE/FINEMAP)"),
    "flames":        (flames_main,        "Run FLAMES effector-gene pipeline"),
    "formatter":     (formatter_main,     "Format sumstats → MAGMA/PoPS inputs"),
    "harmonisation": (harmonisation_main, "Run harmonisation pipeline"),
    "heritability":  (heritability_main,  "Run heritability estimation (LDSC)"),
    "imputation":    (imputation_main,    "Run summary-statistic imputation module"),
    "ld_clump":      (ld_clump_main,      "Run LD clumping and pruning module"),
    "magma":         (magma_main,         "Run MAGMA gene/pathway analysis"),
    "magmacovar":    (magmacovar_main,    "Run MAGMA gene-property (covariate) model"),
    "manhattan":     (manhattan_main,     "Generate Manhattan/QQ plots"),
    "pops":          (pops_main,          "Run PoPS gene-prioritisation module"),
    "qc":            (qc_main,            "Run QC summary module"),
    "pathway_enrichmnet":(enrichment_main,"Run pathway enrichment module"),
    "sumstat_filter":(sumstat_filter_main,"Summary-statistics filtering module"),
    "pipeline":(pipeline_main,"Summary-statistics pipeline module"),
}

# Sort alphabetically
MODULES = dict(sorted(MODULES.items(), key=lambda x: x[0].lower()))

console = Console()

# -------------------------------------------------------------------
# Rich Global Help
# -------------------------------------------------------------------
def print_global_help(prog: str = "postgwas") -> None:
    """Render Rich-formatted top-level help."""

    # Header banner
    header = Text("PostGWAS — Unified Toolkit for Post-GWAS Analyses", style="bold cyan")
    console.print(Panel(header, expand=False))

    # Module table
    table = Table(title="Available Modules", title_style="bold magenta", padding=(0,1))
    table.add_column("Command", style="cyan", no_wrap=True)
    table.add_column("Description", style="green")

    for name, (_, desc) in MODULES.items():
        table.add_row(name, desc)

    console.print(table)

    console.print("\n[bold yellow]Usage:[/bold yellow]  postgwas <module> --help")
    console.print("Example: [green]postgwas finemap --help[/green]\n")


# -------------------------------------------------------------------
# Dispatcher
# -------------------------------------------------------------------
def main():
    argv = sys.argv
    prog = argv[0]

    # If no command or help requested
    if len(argv) == 1 or argv[1] in ("-h", "--help"):
        print_global_help(prog)
        sys.exit(0)

    cmd = argv[1]

    # Invalid module name
    if cmd not in MODULES:
        console.print(f"[red]{prog}: error: unknown module '{cmd}'[/red]\n")
        print_global_help(prog)
        sys.exit(1)

    module_main, _ = MODULES[cmd]

    # Replace argv so the submodule sees correct prog name
    sys.argv = [f"{prog}-{cmd}"] + argv[2:]

    return module_main()


if __name__ == "__main__":
    main()
