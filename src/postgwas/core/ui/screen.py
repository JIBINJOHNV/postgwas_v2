"""Shared user-facing terminal formatting.

File logs deliberately remain plain and machine-searchable. These helpers are
only for terminal output, where semantic symbols and stable indentation make
results easier to scan across every PostGWAS module.
"""

import textwrap


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


def screen_field(kind, label, value, width=128, indent=6, label_width=20):
    """Return one wrapped field whose continuation lines align with its value."""
    prefix = "%s%s  %-*s : " % (
        " " * int(indent), SYMBOLS[kind], int(label_width), str(label),
    )
    continuation = " " * (int(indent) + 4 + int(label_width) + 3)
    return "\n".join(textwrap.wrap(
        prefix + str(value),
        width=max(72, int(width)),
        subsequent_indent=continuation,
        break_long_words=False,
        break_on_hyphens=False,
    ) or [prefix.rstrip()])


__all__ = ["SYMBOLS", "screen_field", "screen_line"]
