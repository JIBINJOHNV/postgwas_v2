#!/usr/bin/env python3
"""
FLAMES Command-Line Interface
Fine-mapping → Gene Prioritization

Supports:
  1) direct     — FLAMES using precomputed SuSiE, MAGMA, PoPS outputs
  2) pipeline   — Full workflow from GWAS VCF → FineMap → MAGMA → PoPS → FLAMES
"""

import argparse
import sys
from pathlib import Path
from postgwas.utils.main import auto_detect_workers, create_default_log, validate_path
import multiprocessing
from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich.rule import Rule
from rich.markdown import Markdown

# ✨ Rich help formatter
from rich_argparse import RichHelpFormatter


# -------------------------------------------------------------------
# VALIDATORS
# -------------------------------------------------------------------

DEFAULT_OUTDIR = Path.cwd() / "flame_analysis"
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_MODEL_DIR = SCRIPT_DIR / "model"

# -------------------------------------------------------------------
# DISPATCH IMPLEMENTATION
# -------------------------------------------------------------------

def run_flames_direct(args):
    print("\n🔥 Running FLAMES Direct Mode with:")
    print(args)
    print("\n➡ Plug in: FLAMES annotation + FLAMES scoring here.\n")


def run_flames_pipeline(args):
    print("\n🚀 Running FLAMES Full Pipeline with:")
    print(args)
    print("\n➡ Plug in: GWAS → FineMap → MAGMA → PoPS → FLAMES pipeline here.\n")


# ===================================================================
# MAIN CLI
# ===================================================================
def main():
    # ---------------------------------------------------------------
    # Banner when no arguments
    # ---------------------------------------------------------------
    if len(sys.argv) == 1:
        console = Console()

        title_text = Text.assemble(
            ("🔥 FLAMES", "bold magenta"),
            (" — Fine-mapped Locus Assessment Model of Effector geneS", "bold white"),
        )

        # Top Title Panel (Gradient-like using two colors)
        console.print(
            Panel(
                title_text,
                border_style="bright_magenta",
                padding=(1, 2),
                expand=False
            )
        )

        console.print(Rule("[bold cyan]✨ Overview[/bold cyan]"))

        console.print(
            "[bold white]FLAMES integrates fine-mapping + functional annotation + machine learning\n"
            "to prioritize the most plausible causal gene(s) at each locus.[/bold white]\n"
        )

        console.print(Rule("[bold cyan]🎯 Modes of Operation[/bold cyan]", style="cyan"))

        console.print(
            "[bold violet]1. flames direct[/bold violet]\n"
            "    • Use FLAMES with precomputed outputs (SuSiE, MAGMA, PoPS)\n"
            "    • Fastest mode; assumes all upstream components already run\n\n"
            "[bold violet]2. flames pipeline[/bold violet]\n"
            "    • Complete workflow: VCF → Harmonisation → MAGMA → SuSiE → PoPS → FLAMES\n"
            "    • End-to-end reproducible causal gene prioritization\n"
        )

        console.print(Rule("[bold cyan]📘 Quick Start[/bold cyan]"))

        console.print(
            "[green]# Run direct mode (precomputed inputs)[/green]\n"
            "  [bold white]flames direct[/bold white] \\\n"
            "      --sample_id PGC3_SCZ \\\n"
            "      --magma_out results.magma.genes.out \\\n"
            "      --susie susie_output_dir \\\n"
            "      --pops pops_scores.preds \\\n"
            "      --covariates covars.tsv\n\n"
            "[green]# Run full pipeline[/green]\n"
            "  [bold white]flames pipeline[/bold white] \\\n"
            "      --sample_id PGC3_SCZ \\\n"
            "      --vcf gwas.vcf.gz \\\n"
            "      --ref_panel 1000G.EUR.bed \\\n"
            "      --gene_loc magma.genes.loc \\\n"
            "      --gtex gtex.txt \\\n"
            "      --out output_folder\n"
        )

        console.print(Rule("[bold cyan]📑 Help Commands[/bold cyan]"))

        console.print(
            "  [yellow]flames direct --help[/yellow]   • Detailed help for direct mode\n"
            "  [yellow]flames pipeline --help[/yellow] • Detailed help for full pipeline mode\n"
        )

        console.print()
        sys.exit(0)


    parser = argparse.ArgumentParser(
        prog="flames",
        description="FLAMES (Fine-Mapping + Gene Prioritization)",
        formatter_class=RichHelpFormatter,
    )

    # ---------------------------------------------------------------
    # COMMON SHARED OPTIONS
    # ---------------------------------------------------------------
    common_parser = argparse.ArgumentParser(
        add_help=False, formatter_class=RichHelpFormatter
    )

    common_group = common_parser.add_argument_group("🔧 Common Settings")

    common_group.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging."
    )

    common_group.add_argument(
        "--num_workers",
        type=int,
        default=auto_detect_workers(),
        help="Number of CPU workers (auto-detected, default: %(default)s)."
    )

    common_group.add_argument(
        "--out_dir",
        type=validate_path(create_if_missing=True, must_be_dir=True),
        default=str(DEFAULT_OUTDIR),
        help="Output directory for FLAMES results (auto-created if missing). Default: %(default)s"
    )

    common_group.add_argument(
        "--sample_id",
        type=str,
        required=True,
        help="Sample identifier / GWAS phenotype name (e.g., 'PGC3_SCZ_european')."
    )

    # ---------------------------------------------------------------
    # FLAMES SHARED SETTINGS (annotation, tools, and model config)
    # ---------------------------------------------------------------
    common_flame_related = common_parser.add_argument_group(
        "🔥 FLAMES Annotation, Tools, and Model Configuration"
    )

    # ---- Genome build -----------------------------------------------------
    common_flame_related.add_argument(
        "-b", "--build",
        type=str,
        choices=["GRCh37", "GRCh38"],
        default="GRCh37",
        help=(
            "Genome build (GRCh37 or GRCh38). Must be consistent with your GWAS VCF / sumstat file, "
            "VEP cache, CADD scores, LD reference panel, and FLAMES annotation assets. "
            "(default: %(default)s)"
        )
    )

    # ---- Annotation directory ---------------------------------------------
    common_flame_related.add_argument(
        "-a", "--annotation_dir",
        type=validate_path(
            must_exist=True,
            must_be_dir=True,
            dir_must_have_files=True,
            dir_require_nonempty_files=True
        ),
        required=True,
        help=(
            "Directory containing required FLAMES annotation files "
            "(must contain the official FLAMES annotation dataset)."
        )
    )

    # ---- VEP executable ----------------------------------------------------
    common_flame_related.add_argument(
        "-cv", "--cmd_vep",
        type=str,
        default=None,
        help="Path to the VEP executable. Must match the --build version. (default: %(default)s)"
    )

    # ---- VEP cache ---------------------------------------------------------
    common_flame_related.add_argument(
        "-vc", "--vep_cache",
        type=validate_path(must_exist=True, must_be_dir=True),
        default=None,
        help=(
            "Path to the VEP cache directory. Must match --build (GRCh37 or GRCh38). "
            "(default: %(default)s)"
        )
    )

    # ---- tabix -------------------------------------------------------------
    common_flame_related.add_argument(
        "-t", "--tabix",
        type=str,
        default=None,
        help=(
            "Path to the tabix executable (used only when local CADD scores are provided). "
            "(default: %(default)s)"
        )
    )

    # ---- Local CADD file ---------------------------------------------------
    common_flame_related.add_argument(
        "-cf", "--CADD_file",
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
            allowed_suffixes=[".tsv", ".tsv.gz", ".txt", ".txt.gz"]
        ),
        default=None,
        help=(
            "Local CADD score file (optional). Must match genome build (--build). "
            "Supports .tsv/.txt (gzipped allowed). (default: %(default)s)"
        )
    )

    # ---- FLAMES ML Model Settings ------------------------------------------
    common_flame_related.add_argument(
        "--weight",
        type=float,
        default=0.725,
        help="XGBoost weight used in FLAMES (PoPS = 1 - weight). (default: %(default)s)"
    )

    common_flame_related.add_argument(
        "--distance",
        type=int,
        default=750000,
        help="Gene inclusion distance threshold around credible variants (bp). (default: %(default)s)"
    )

    common_flame_related.add_argument(
        "--modelpath",
        type=validate_path(must_exist=True, must_be_dir=True),
        default=DEFAULT_MODEL_DIR,
        help=(
            "Directory containing the pretrained FLAMES XGBoost model "
            "(auto-detected relative to installation). (default: %(default)s)"
        )
    )






    # ---------------------------------------------------------------
    # SUBCOMMAND ROOT
    # ---------------------------------------------------------------
    subparsers = parser.add_subparsers(
        dest="cmd",
        metavar="",
        help="Choose a FLAMES mode",
    )

    # ===============================================================
    # MODE 1 — DIRECT
    # ===============================================================
    parser_direct = subparsers.add_parser(
        "direct",
        parents=[common_parser],
        help="Run FLAMES using precomputed MAGMA, SuSiE, and PoPS outputs.",
        formatter_class=RichHelpFormatter,
    )

    direct_io = parser_direct.add_argument_group("1. Input / Output Files")

    direct_io.add_argument(
        "--magma_out",
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
            allowed_suffixes=[".genes.out"]
        ),
        required=True,
        help=(
            "MAGMA gene-level results file (.genes.out). "
            "This file must be generated by the [bold green]postgwas gene_assoc[/bold green] module."
        )
    )

    direct_io.add_argument(
        "--magma_covar",
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
            allowed_suffixes=[".gsa.out"]
        ),
        required=True,
        help=(
            "MAGMA gene-covariate results file (.gsa.out). "
            "This file must be generated by the [bold green]postgwas magma_covar[/bold green] module "
            "using FLAMES-compatible settings."
        )
    )

    direct_io.add_argument(
        "--susie",
        type=validate_path(
            must_exist=True,
            must_be_dir=True,
            dir_must_have_files=True,
            dir_require_nonempty_files=True
        ),
        required=True,
        help=(
            "Directory containing SuSiE credible-set outputs in FLAMES-compatible format. "
            "These files must be generated by the [bold green]postgwas finemap[/bold green] module."
        )
    )

    direct_io.add_argument(
        "--pops",
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
            allowed_suffixes=[".preds"]
        ),
        required=True,
        help=(
            "PoPS gene-score file (.preds). "
            "This file must be generated by the [bold green]postgwas pops[/bold green] module."
        )
    )



    # ===============================================================
    # MODE 2 — PIPELINE
    # ===============================================================
    parser_pipeline = subparsers.add_parser(
        "pipeline",
        parents=[common_parser],
        help="Run full workflow: VCF → MAGMA → SuSiE → PoPS → FLAMES.",
        formatter_class=RichHelpFormatter,
    )

    pipeline_gwas = parser_pipeline.add_argument_group("1. GWAS Inputs")

    pipeline_gwas.add_argument(
        "--vcf",
        type=validate_path(
            must_exist=True, must_be_file=True, must_not_be_empty=True,
            allowed_suffixes=[".vcf", ".vcf.gz"]
        ),
        required=True,
        help="GWAS VCF from postgwas harmonisation."
    )

    pipeline_gwas.add_argument(
        "--ref_panel",
        type=validate_path(
            must_exist=True, must_be_file=True,
            allowed_suffixes=[".bed"]
        ),
        required=True,
        help="PLINK .bed file from LD reference panel."
    )

    pipeline_magma = parser_pipeline.add_argument_group("2. MAGMA Settings")

    pipeline_magma.add_argument(
        "--gene_loc",
        type=validate_path(
            must_exist=True, must_be_file=True,
            allowed_suffixes=[".loc", ".txt"]
        ),
        required=True,
        help="MAGMA gene location file."
    )

    pipeline_magma.add_argument(
        "--gtex",
        type=validate_path(must_exist=True, must_be_file=True),
        required=True,
        help="GTEx tissue expression matrix (MAGMA tissue analysis)."
    )

    pipeline_susie = parser_pipeline.add_argument_group("3. SuSiE Settings")

    pipeline_susie.add_argument(
        "--susie_L",
        type=int,
        default=10,
        help="Max number of SuSiE components."
    )

    pipeline_susie.add_argument(
        "--susie_skip_mhc",
        action="store_true",
        help="Skip MHC region during fine-mapping."
    )

    pipeline_pops = parser_pipeline.add_argument_group("4. PoPS Settings")

    pipeline_pops.add_argument(
        "--feature_mat_prefix",
        type=str,
        help="Prefix for PoPS feature matrix splits."
    )

    pipeline_pops.add_argument(
        "--num_feature_chunks",
        type=int,
        default=1,
        help="Number of PoPS feature matrix chunks."
    )

    pipeline_runtime = parser_pipeline.add_argument_group("5. Runtime Settings")

    pipeline_runtime.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of compute threads."
    )

    pipeline_runtime.add_argument(
        "--out",
        type=validate_path(create_if_missing=True, must_be_dir=True),
        required=True,
        help="Output directory for FLAMES pipeline."
    )

    # ===============================================================
    # DISPATCH
    # ===============================================================
    args = parser.parse_args()

    if args.cmd == "direct":
        return run_flames_direct(args)

    if args.cmd == "pipeline":
        return run_flames_pipeline(args)


# -------------------------------------------------------------------
# ENTRY POINT
# -------------------------------------------------------------------
if __name__ == "__main__":
    main()
