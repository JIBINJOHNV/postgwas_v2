"""Shared user-facing terminal formatting.

File logs deliberately remain plain and machine-searchable. These helpers are
only for terminal output, where semantic symbols and stable indentation make
results easier to scan across every PostGWAS module.
"""

import shutil
import textwrap

from rich.cells import cell_len


SYMBOLS = {
    "genetic": "🧬",
    "count": "🧮",
    "analysis": "🔬",
    "decision": "🧠",
    "info": "🔹",
    "loss": "🔻",
    "success": "✅",
    "warning": "⚠️",
    "error": "❌",
    "run": "▶️",
    "retry": "🔁",
}


def screen_line(kind, text, indent=4):
    """Return one indented terminal line with a meaning-specific symbol."""
    return "%s%s  %s" % (" " * int(indent), SYMBOLS[kind], str(text))


def screen_field(
    kind,
    label,
    value,
    width=None,
    indent=6,
    label_width=20,
    break_long_values=False,
):
    """Return one field with fixed label and value columns.

    Labels wider than ``label_width`` wrap inside the label column. Values then
    start at the same terminal cell for short and long labels alike.
    """
    terminal_width = (
        shutil.get_terminal_size(fallback=(128, 24)).columns
        if width is None
        else int(width)
    )
    indent = int(indent)
    label_width = int(label_width)
    label_lines = textwrap.wrap(
        str(label),
        width=label_width,
        break_long_words=True,
        break_on_hyphens=False,
    ) or [""]
    label_indent = " " * (indent + cell_len(SYMBOLS[kind]) + 2)
    leading_lines = [
        "%s%s  %s" % (
            " " * indent if index == 0 else label_indent,
            SYMBOLS[kind] if index == 0 else "",
            fragment,
        )
        for index, fragment in enumerate(label_lines[:-1])
    ]
    final_label_prefix = (
        "%s%s  " % (" " * indent, SYMBOLS[kind])
        if len(label_lines) == 1
        else label_indent
    )
    prefix = "%s%-*s : " % (
        final_label_prefix, label_width, label_lines[-1],
    )
    continuation = " " * cell_len(prefix)
    unicode_width_adjustment = max(0, cell_len(prefix) - len(prefix))
    value_lines = textwrap.wrap(
        prefix + str(value),
        width=max(72, terminal_width - unicode_width_adjustment),
        subsequent_indent=continuation,
        break_long_words=bool(break_long_values),
        break_on_hyphens=False,
    ) or [prefix.rstrip()]
    return "\n".join([*leading_lines, *value_lines])


__all__ = ["SYMBOLS", "screen_field", "screen_line"]
