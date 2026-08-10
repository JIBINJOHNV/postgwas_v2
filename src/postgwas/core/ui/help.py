"""Consistent, readable command-line help formatting."""

from collections.abc import Sequence
import re
import shutil

from rich.markup import escape
from rich_argparse import RichHelpFormatter, RawTextRichHelpFormatter


_PARENTHETICAL_DEFAULT = re.compile(
    r"\((?P<label>(?:(?:configured(?: YAML)?|effective|packaged) )?default):\s*"
    r"(?P<value>[^)\n]+)\)",
    flags=re.IGNORECASE,
)
_STANDALONE_DEFAULT = re.compile(
    r"(?<!\])\b(?P<label>(?:(?:Configured(?: YAML)?|Effective|Packaged) )?Default):"
    r"\s*(?P<value>[^\n]+)",
    flags=re.IGNORECASE,
)
_DEFAULT_SEPARATOR = re.compile(
    r"((?:(?:configured(?: YAML)?|effective|packaged) )?default:) ",
    flags=re.IGNORECASE,
)
_NONBREAKING_DISPLAY_SPACE = "\ue000"


class AlignedRichHelpFormatter(RawTextRichHelpFormatter):
    """Align option descriptions while preserving deliberate example layouts."""

    @staticmethod
    def group_name_formatter(value: str) -> str:
        return value.capitalize() if value == value.lower() else value

    def __init__(self, prog: str, **kwargs) -> None:
        width = kwargs.get("width") or shutil.get_terminal_size().columns
        max_help_position = min(42, max(24, width // 3))
        super().__init__(prog, max_help_position=max_help_position, **kwargs)

    def _expand_help(self, action):
        return style_cli_defaults(super()._expand_help(action))

    def _rich_split_lines(self, text, width):
        """Wrap option help while retaining deliberate description/epilog layout."""
        protected = text.copy()
        protected.plain = _DEFAULT_SEPARATOR.sub(
            r"\1%s" % _NONBREAKING_DISPLAY_SPACE,
            protected.plain,
        )
        lines = RichHelpFormatter._rich_split_lines(self, protected, width)
        for line in lines:
            line.plain = line.plain.replace(_NONBREAKING_DISPLAY_SPACE, " ")
        return lines


def format_cli_default(value, *, label: str = "Default") -> str:
    """Format one configured CLI default with distinct label and value styles."""
    if value is None:
        rendered = "unset"
    elif isinstance(value, bool):
        rendered = str(value).lower()
    else:
        rendered = str(value)
    return "[bold green]%s:[/bold green] [cyan]%s[/cyan]" % (
        escape(label),
        escape(rendered),
    )


def help_with_default(description: str, value, *, label: str = "Default") -> str:
    """Append a centrally styled, configuration-derived default to CLI help."""
    return "%s (%s)." % (
        description.rstrip("."),
        format_cli_default(value, label=label),
    )


def style_cli_defaults(help_text: str) -> str:
    """Color legacy plain-text defaults when shared module help is rendered."""
    def replace(match: re.Match) -> str:
        label = match.group("label")
        label = label[:1].upper() + label[1:]
        return "(%s)" % format_cli_default(
            match.group("value"), label=label,
        )

    styled = _PARENTHETICAL_DEFAULT.sub(replace, help_text)

    def replace_standalone(match: re.Match) -> str:
        return format_cli_default(
            match.group("value"), label=match.group("label"),
        )

    return _STANDALONE_DEFAULT.sub(replace_standalone, styled)


def format_cli_examples(
    *examples: tuple[str, str, Sequence[str]],
    notes: Sequence[str] = (),
    title: str = "Examples",
) -> str:
    """Render consistently indented, copyable multi-line CLI examples."""
    lines = [title]
    for label, command, arguments in examples:
        lines.extend(("", label))
        arguments = tuple(arguments)
        lines.append("  %s%s" % (command, " \\" if arguments else ""))
        for index, argument in enumerate(arguments):
            continuation = " \\" if index < len(arguments) - 1 else ""
            lines.append("    %s%s" % (argument, continuation))
    if notes:
        lines.extend(("", "Notes"))
        lines.extend("  %s" % note for note in notes)
    return "\n".join(lines)


__all__ = [
    "AlignedRichHelpFormatter",
    "format_cli_default",
    "format_cli_examples",
    "help_with_default",
    "style_cli_defaults",
]
