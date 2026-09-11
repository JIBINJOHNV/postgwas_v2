"""Shared user-facing terminal formatting.

File logs deliberately remain plain and machine-searchable. These helpers are
only for terminal output, where semantic symbols and stable indentation make
results easier to scan across every PostGWAS module.
"""

import shutil
import textwrap
from contextlib import contextmanager
from contextvars import ContextVar
from functools import lru_cache

from rich.cells import cell_len
from rich.text import Text

from postgwas.core.ui.terminal_text import TerminalTextDecoder


SYMBOLS = {
    "genetic": "🧬",
    "count": "🧮",
    "analysis": "🔬",
    "decision": "🧠",
    "attention": "🚨",
    "info": "🔹",
    "loss": "🔻",
    "success": "✅",
    # A single-codepoint emoji avoids the platform-dependent column width of
    # the warning-sign plus variation-selector sequence ("⚠️").
    "warning": "❗",
    "error": "❌",
    "run": "▶️",
    "retry": "🔁",
}

_TERMINAL_SETTINGS = ContextVar("postgwas_terminal_settings", default=None)


@lru_cache(maxsize=1)
def default_terminal_settings():
    """Load only the canonical logging defaults, lazily to avoid import cycles."""
    from postgwas.config.loader import _packaged_yaml
    from postgwas.config.models.logging import LoggingConfig

    return LoggingConfig.model_validate(
        _packaged_yaml("defaults/application.yaml")["logging"],
    )


def terminal_settings():
    """Return the resolved run presentation, or canonical defaults outside a run."""
    configured = _TERMINAL_SETTINGS.get()
    return default_terminal_settings() if configured is None else configured


@contextmanager
def terminal_presentation(logging_config):
    """Scope one validated theme to a run without mutating module/global palettes."""
    token = _TERMINAL_SETTINGS.set(logging_config)
    try:
        yield
    finally:
        _TERMINAL_SETTINGS.reset(token)


def terminal_style(kind):
    settings = terminal_settings().terminal_style
    return "none" if settings.color == "never" else settings.styles[kind]


def _screen_indent(indent):
    configured = _TERMINAL_SETTINGS.get()
    if configured is None:
        return int(indent)
    layout = configured.terminal_style
    # Existing callers use the shared 2/6/10/... cell grid (some older cards
    # are offset by two cells). Translate that interface into the configured
    # heading/field levels; these numbers describe the caller protocol, not
    # independent runtime defaults.
    depth = max(0, (int(indent) - 2 + 3) // 4)
    return (layout.heading_indent if depth == 0 else
            layout.field_indent + (depth - 1) * layout.indent_step)

# The shared QC-card vocabulary includes labels up to 28 cells wide
# (for example, ``Reference unmatched retained``). Keeping that complete
# vocabulary on one line gives standalone QC and harmonisation one stable
# evidence-value column in the terminal and both durable logs.
SUMMARY_CARD_LABEL_WIDTH = 28


def normalize_stage_outcome_fields(fields) -> tuple[tuple, ...]:
    """Validate and normalize semantic fields shared by progress and loggers."""
    normalized = []
    for field in fields or ():
        if len(field) not in {2, 3}:
            raise ValueError(
                "Each stage outcome entry must contain kind and label, with "
                "an optional value"
            )
        kind, label = str(field[0]).strip(), str(field[1]).strip()
        if kind not in SYMBOLS or not label:
            raise ValueError(
                "Stage outcome fields require a known symbol kind and label"
            )
        normalized.append(
            (kind, label) if len(field) == 2 else (kind, label, field[2])
        )
    return tuple(normalized)


def screen_line(kind, text, indent=4):
    """Return one indented terminal line with a meaning-specific symbol."""
    return "%s%s  %s" % (" " * _screen_indent(indent), SYMBOLS[kind], str(text))


def screen_section_lines(kind, title, indent=6):
    """Return the shared blank-line and heading pattern for result sections."""
    return ["", "", screen_line(kind, title, indent=indent), ""]


def _wrap_path(value: str, width: int) -> list[str]:
    """Wrap a POSIX path at directory boundaries without changing its text."""
    segments = [part + "/" for part in value.split("/")[:-1]]
    segments.append(value.split("/")[-1])
    lines: list[str] = []
    current = ""
    for segment in segments:
        if cell_len(current + segment) <= width:
            current += segment
            continue
        if current:
            lines.append(current)
        while cell_len(segment) > width:
            split_at = 1
            while (
                split_at < len(segment)
                and cell_len(segment[:split_at + 1]) <= width
            ):
                split_at += 1
            lines.append(segment[:split_at])
            segment = segment[split_at:]
        current = segment
    if current or not lines:
        lines.append(current)
    return lines


def screen_field(
    kind,
    label,
    value,
    width=None,
    indent=6,
    label_width=20,
    break_long_values=False,
    path_value=False,
    separator_column=None,
):
    """Return one field with fixed label and value columns.

    Labels wider than ``label_width`` wrap inside the label column. Values then
    start at the same terminal cell for short and long labels alike. When
    ``separator_column`` is provided, the label width is adjusted so `` : ``
    starts at that terminal-cell column regardless of indentation. Local paths
    may request directory-aware wrapping with ``path_value``.
    """
    terminal_width = (
        shutil.get_terminal_size(fallback=(128, 24)).columns
        if width is None
        else int(width)
    )
    indent = _screen_indent(indent)
    label_width = int(label_width)
    configured = _TERMINAL_SETTINGS.get()
    if configured is not None:
        # One absolute value column even when nested sections use different
        # indentation. Reserve half a narrow terminal for the value itself.
        separator_column = min(
            configured.terminal_style.field_indent
            + cell_len(SYMBOLS[kind]) + 2 + configured.terminal_label_width,
            terminal_width // 2,
        )
        indent = min(indent, max(0, separator_column - cell_len(SYMBOLS[kind]) - 3))
        break_long_values = True
    if separator_column is not None:
        label_width = (
            int(separator_column)
            - indent
            - cell_len(SYMBOLS[kind])
            - 2
        )
        if label_width < 1:
            raise ValueError(
                "separator_column leaves no space for the field label"
            )
    label_lines = textwrap.wrap(
        str(label),
        width=label_width,
        break_long_words=True,
        break_on_hyphens=False,
    ) or [""]
    label_indent = " " * (indent + cell_len(SYMBOLS[kind]) + 2)
    leading_lines = [
        ((" " * indent + SYMBOLS[kind] + "  ") if index == 0 else label_indent)
        + fragment
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
    if path_value:
        path_lines = _wrap_path(
            str(value), max(1, terminal_width - cell_len(prefix)),
        )
        value_lines = [
            prefix + path_lines[0],
            *(continuation + line for line in path_lines[1:]),
        ]
    elif break_long_values:
        # Wrap the value independently: a long first token must not push the
        # complete value onto the next row after an otherwise empty label.
        fragments = textwrap.wrap(
            str(value),
            width=max(1, terminal_width - cell_len(prefix)),
            break_long_words=True,
            break_on_hyphens=False,
        ) or [""]
        value_lines = [
            prefix + fragments[0],
            *(continuation + fragment for fragment in fragments[1:]),
        ]
    else:
        value_lines = textwrap.wrap(
            prefix + str(value),
            width=max(72, terminal_width - unicode_width_adjustment),
            subsequent_indent=continuation,
            break_long_words=bool(break_long_values),
            break_on_hyphens=False,
        ) or [prefix.rstrip()]
    return "\n".join([*leading_lines, *value_lines])


def screen_summary_card(
    kind,
    title,
    fields,
    *,
    title_indent=8,
    field_indent=12,
    label_width=SUMMARY_CARD_LABEL_WIDTH,
):
    """Return one shared summary card with aligned evidence fields."""
    normalized = normalize_stage_outcome_fields(fields)
    incomplete = [field for field in normalized if len(field) != 3]
    if incomplete:
        raise ValueError(
            "Summary-card fields must contain kind, label, and value"
        )
    lines = ["", screen_line(kind, title, indent=title_indent)]
    for field_kind, label, value in normalized:
        lines.extend(screen_field(
            field_kind,
            label,
            value,
            indent=field_indent,
            label_width=label_width,
        ).splitlines())
    return lines


def render_action_plan(
    title,
    items,
    *,
    empty_message,
    title_indent=4,
    item_indent=8,
    footer=(),
):
    """Render one ordered action plan with shared filtering-style spacing."""
    normalized = []
    for item in items or ():
        action = str(item.get("action") or "").strip().upper()
        label = str(item.get("label") or "").strip()
        if not action or not label:
            raise ValueError("Action-plan items require an action and label")
        normalized.append((action, label))

    footer_fields = normalize_stage_outcome_fields(footer)
    if any(len(field) != 2 for field in footer_fields):
        raise ValueError("Action-plan footer entries require kind and text")

    lines = [
        "",
        screen_line("analysis", title, indent=title_indent),
        "",
    ]
    if normalized:
        lines.extend(
            "%s%d. %s · %s" % (" " * item_indent, number, action, label)
            for number, (action, label) in enumerate(normalized, 1)
        )
    else:
        lines.append(screen_line(
            "info", empty_message, indent=item_indent,
        ))
    if footer_fields:
        lines.append("")
        lines.extend(
            screen_line(kind, text, indent=item_indent)
            for kind, text in footer_fields
        )
    lines.append("")
    return "\n".join(lines)


def style_screen_block(value: str | Text) -> Text:
    """Apply shared semantic colors while preserving a block's printable text.

    Section headings and standalone status lines receive their complete
    semantic style. Ordinary aligned fields color only the symbol-and-label
    prefix so explanatory values remain neutral and easy to read. Attention
    fields color their complete wrapped content. Module-specific styles are
    deliberately not accepted: all roles use logging.terminal_style.styles.
    Incoming Rich styles are discarded, while literal text (including square
    brackets in paths and diagnostics) is preserved. Terminal controls are
    removed before display and logging. Rich omits colour from
    non-terminal logs. No severity is inferred from words in a message.
    """
    value = value.plain if isinstance(value, Text) else str(value)
    value = TerminalTextDecoder("utf-8").clean(value)
    rendered = Text(style=terminal_style("text"))
    active_kind = None
    active_style = None
    active_indent = 0
    for line in value.splitlines(keepends=True):
        body = line.rstrip("\r\n")
        ending = line[len(body):]
        stripped = body.lstrip()
        field_separator = body.find(" : ")
        kind = next(
            (
                candidate
                for candidate, symbol in SYMBOLS.items()
                if stripped.startswith(symbol)
                and stripped[len(symbol):len(symbol) + 1].isspace()
            ),
            None,
        )
        if kind is not None:
            active_kind = kind
            active_style = terminal_style(kind)
            active_indent = len(body) - len(stripped)
        elif not stripped or len(body) - len(stripped) <= active_indent:
            active_kind = None
            active_style = None

        if active_style is None:
            rendered.append(body)
        elif active_kind in {"attention", "warning", "error"}:
            # High-attention fields deliberately style their complete wrapped
            # value. Ordinary semantic fields retain neutral values so dense
            # scientific summaries remain readable.
            rendered.append(body, style=active_style)
        elif field_separator >= 0:
            value_start = field_separator + len(" : ")
            rendered.append(body[:value_start], style=active_style)
            rendered.append(body[value_start:])
        elif kind is not None:
            rendered.append(body, style=active_style)
        else:
            rendered.append(body)
        rendered.append(ending)
    return rendered


__all__ = [
    "default_terminal_settings",
    "terminal_settings",
    "terminal_presentation",
    "terminal_style",
    "SUMMARY_CARD_LABEL_WIDTH",
    "SYMBOLS",
    "screen_field",
    "screen_line",
    "screen_section_lines",
    "screen_summary_card",
    "render_action_plan",
    "style_screen_block",
]
