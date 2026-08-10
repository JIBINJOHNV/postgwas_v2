# postgwas/flames/clis/cli_utils.py
from __future__ import annotations

import argparse
import re
from typing import Any

from rich_argparse import RichHelpFormatter

# ------------------------------------------------------------------
# Central formatter – subclass once, reuse everywhere
# ------------------------------------------------------------------
class FlamesHelpFormatter(RichHelpFormatter):
    """Rich formatter that inserts a visually empty line after every optional argument."""

    # rich_argparse ≥ 1.4  exposes this hook – safest & cleanest
    def add_argument(self, action: argparse.Action) -> None:
        super().add_argument(action)
        if action.option_strings and action.help is not argparse.SUPPRESS:
            self.add_text("")  # blank line


# ------------------------------------------------------------------
# One-line helper used by every parser / group
# ------------------------------------------------------------------
def install_formatter(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Install FlamesHelpFormatter on parser *and* all its argument groups."""
    parser.formatter_class = FlamesHelpFormatter
    for group in parser._action_groups:
        group.formatter_class = FlamesHelpFormatter
    return parser