# postgwas/flames/clis/cli_banner.py

from rich.console import Console
from rich.panel import Panel
from rich.text import Text

console = Console()

def short_rule(title: str = "", style: str = "cyan", total: int = 50):
    """Modern professional separators (works in all Rich versions)."""
    if title:
        title_text = f" {title} "
        side = (total - len(title_text)) // 2
        line = f"{'─' * side}{title_text}{'─' * side}"
        if len(line) < total:
            line += "─"
    else:
        line = "─" * total
    console.print(f"[{style}]{line}[/{style}]")


def print_flames_banner():
    """Compact modern FLAMES banner."""

    # Compact banner width
    banner_width = 62

    title = Text.assemble(
        ("🔥  FLAMES", "bold magenta"),
        (" — For GWAS gene prioritization", "bold white"),
    )

    console.print(
        Panel(
            title,
            padding=(1, 2),
            border_style="magenta",
            width=banner_width,   # <<<<< FIXED WIDTH, NOT FULL SCREEN
        )
    )

    console.print("[green]Fine-mapping ➜ MAGMA ➜ PoPS ➜ Effector gene prioritisation.[/green]\n")

    short_rule("Available Subcommands")
    console.print(
        "  [bold violet]direct[/bold violet]   – Run FLAMES using pre-computed SuSiE/MAGMA/PoPS.\n"
        "  [bold violet]pipeline[/bold violet] – Full workflow from VCF → SuSiE → MAGMA → PoPS → FLAMES.\n"
    )

    short_rule("Usage")
    console.print("  [yellow]flames direct --help[/yellow]   Show help for direct mode.")
    console.print("  [yellow]flames pipeline --help[/yellow] Show help for pipeline mode.\n")
