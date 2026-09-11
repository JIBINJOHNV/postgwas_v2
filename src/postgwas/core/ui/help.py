"""Consistent, readable command-line help formatting."""

import argparse
from collections.abc import Collection, Iterable, Sequence
import re
import shutil

from rich.markup import escape
from rich.text import Text
from rich_argparse import RichHelpFormatter, RawTextRichHelpFormatter


_PARENTHETICAL_DEFAULT = re.compile(
    r"\((?P<label>(?:(?:configured(?: YAML)?|effective|packaged|upstream) )?default):\s*"
    r"(?P<value>[^)\n]+)\)",
    flags=re.IGNORECASE,
)
_STANDALONE_DEFAULT = re.compile(
    r"(?<!\])\b(?P<label>(?:(?:Configured(?: YAML)?|Effective|Packaged|Upstream) )?Default):"
    r"\s*(?P<value>[^\n]+)",
    flags=re.IGNORECASE,
)
_STYLED_DEFAULT = re.compile(
    r"\[bold green\](?P<label>[^\]]*default):\[/bold green\]\s*"
    r"\[(?P<style>cyan|green)\](?P<value>.*?)\[/(?P=style)\]",
    flags=re.IGNORECASE,
)
_COLORED_VALUE_DEFAULT = re.compile(
    r"(?<!\])\b(?P<label>(?:(?:Configured(?: YAML)?|Effective|Packaged|Upstream) )?Default):"
    r"\s*\[(?:cyan|green)\](?P<value>.*?)\[/(?:cyan|green)\]",
    flags=re.IGNORECASE,
)
_DEFAULT_SEPARATOR = re.compile(
    r"((?:(?:configured(?: YAML)?|effective|packaged|upstream) )?default:) ",
    flags=re.IGNORECASE,
)
_LEADING_REQUIRED = re.compile(
    r"^(?:"
    r"\[bold(?: bright_red| red)?\]Required:?\[/bold(?: bright_red| red)?\]"
    r"|Required"
    r")(?:\s+in\s+standalone\s+mode)?\s*[:.]?\s*",
    flags=re.IGNORECASE,
)
_REQUIRED_PARENTHETICAL = re.compile(r"\s*\(required\)", flags=re.IGNORECASE)
_HELP_REQUIRED_ATTRIBUTE = "_postgwas_required_help"
_NONBREAKING_DISPLAY_SPACE = "\ue000"
_STYLED_DEFAULT_PLACEHOLDER = "\ue001%d\ue002"
_WRAPPED_HELP_BLOCK_PREFIX = "\ue003"
_STYLED_CHOICES_LABEL = (
    "[bold bright_yellow]Available options:[/bold bright_yellow]"
)


class AlignedRichHelpFormatter(RawTextRichHelpFormatter):
    """Align option descriptions while preserving deliberate example layouts."""

    @staticmethod
    def group_name_formatter(value: str) -> str:
        return value.capitalize() if value == value.lower() else value

    def __init__(self, prog: str, **kwargs) -> None:
        width = kwargs.get("width") or shutil.get_terminal_size().columns
        max_help_position = min(42, max(24, width // 3))
        super().__init__(prog, max_help_position=max_help_position, **kwargs)

    def _get_help_string(self, action):
        help_text = super()._get_help_string(action)
        if help_text is None:
            return None
        required = action.required or getattr(
            action, _HELP_REQUIRED_ATTRIBUTE, False,
        )
        return style_cli_requirement(
            style_cli_defaults(help_text), required=required,
        )

    def _expand_help(self, action):
        help_text = super()._expand_help(action)
        if _is_finite_choice_action(action):
            return help_with_choices(help_text, action.choices)
        return help_text

    def _rich_expand_help(self, action):
        help_text = super()._rich_expand_help(action)
        if (
            not _is_finite_choice_action(action)
            or _STYLED_CHOICES_LABEL in str(action.help)
        ):
            return help_text
        help_text.append(" ")
        help_text.append(
            Text.from_markup(
                "%s." % format_cli_choices(action.choices),
                style="argparse.help",
            )
        )
        return help_text

    def _rich_split_lines(self, text, width):
        """Wrap option help while retaining deliberate description/epilog layout."""
        protected = text.copy()
        protected.plain = _DEFAULT_SEPARATOR.sub(
            r"\1%s" % _NONBREAKING_DISPLAY_SPACE,
            protected.plain,
        )
        lines = []
        for logical_line in RawTextRichHelpFormatter._rich_split_lines(
            self, protected, width,
        ):
            wrapped = RichHelpFormatter._rich_split_lines(
                self, logical_line, width,
            )
            lines.extend(wrapped or [logical_line])
        for line in lines:
            line.plain = line.plain.replace(_NONBREAKING_DISPLAY_SPACE, " ")
        return lines

    def _rich_fill_text(self, text, width, indent):
        """Wrap explicitly marked prose blocks while preserving their line breaks."""
        if not text.plain.startswith(_WRAPPED_HELP_BLOCK_PREFIX):
            return super()._rich_fill_text(text, width, indent)
        protected = text[len(_WRAPPED_HELP_BLOCK_PREFIX):]
        lines = []
        for logical_line in protected.split():
            leading_spaces = len(logical_line.plain) - len(
                logical_line.plain.lstrip(" ")
            )
            content = logical_line[leading_spaces:]
            wrapped = RichHelpFormatter._rich_split_lines(
                self, content, max(1, width - leading_spaces),
            )
            prefix = Text(" " * leading_spaces)
            lines.extend(prefix + line for line in (wrapped or [content]))
        return Text("\n").join(indent + line for line in lines) + "\n\n"


def format_cli_default(value) -> str:
    """Format one CLI default with the single public ``Default:`` label."""
    if value is None:
        rendered = "unset"
    elif isinstance(value, bool):
        rendered = str(value).lower()
    else:
        rendered = str(value)
    return "[bold green]Default:[/bold green] [green]%s[/green]" % escape(
        rendered
    )


def help_with_default(description: str, value) -> str:
    """Append a centrally styled default with the exact public label."""
    return "%s (%s)." % (
        description.rstrip("."),
        format_cli_default(value),
    )


def help_with_conditional_requirement(description: str, condition: str) -> str:
    """Prefix conditionally mandatory help without changing parser validation."""
    rendered_condition = escape(condition.strip().rstrip("."))
    return (
        "[bold red]Required:[/bold red] [bold]%s[/bold].\n%s"
        % (rendered_condition, description.lstrip())
    )


def format_cli_choices(choices: Iterable[object]) -> str:
    """Format the accepted values for one choice-constrained CLI option."""
    rendered = ", ".join(escape(str(choice)) for choice in choices)
    return "%s [bright_yellow]%s[/bright_yellow]" % (
        _STYLED_CHOICES_LABEL,
        rendered,
    )


def help_with_choices(description: str, choices: Iterable[object]) -> str:
    """Append a consistently styled list of accepted CLI values once."""
    if _STYLED_CHOICES_LABEL in description:
        return description
    choice_values = tuple(choices)
    if not choice_values:
        return description
    rendered = format_cli_choices(choice_values)
    return "%s %s." % (description.rstrip(), rendered)


def _is_finite_choice_action(action: argparse.Action) -> bool:
    """Return whether the action validates values against finite choices."""
    return (
        isinstance(action.choices, Collection)
        and bool(action.choices)
    )


def style_cli_defaults(help_text: str) -> str:
    """Color legacy plain-text defaults when shared module help is rendered."""
    protected_defaults: list[str] = []

    def protect_styled(match: re.Match) -> str:
        protected_defaults.append(
            format_cli_default(match.group("value"))
        )
        return _STYLED_DEFAULT_PLACEHOLDER % (len(protected_defaults) - 1)

    styled = _STYLED_DEFAULT.sub(protect_styled, help_text)

    def protect_colored_value(match: re.Match) -> str:
        protected_defaults.append(
            format_cli_default(match.group("value"))
        )
        return _STYLED_DEFAULT_PLACEHOLDER % (len(protected_defaults) - 1)

    styled = _COLORED_VALUE_DEFAULT.sub(protect_colored_value, styled)

    def replace(match: re.Match) -> str:
        return "(%s)" % format_cli_default(match.group("value"))

    styled = _PARENTHETICAL_DEFAULT.sub(replace, styled)

    def replace_standalone(match: re.Match) -> str:
        return format_cli_default(match.group("value"))

    styled = _STANDALONE_DEFAULT.sub(replace_standalone, styled)
    for index, protected in enumerate(protected_defaults):
        styled = styled.replace(_STYLED_DEFAULT_PLACEHOLDER % index, protected)
    return styled


def style_cli_requirement(help_text: str, *, required: bool = False) -> str:
    """Put one red ``Required:`` label at the start of required option help."""
    legacy_required = _LEADING_REQUIRED.match(help_text) is not None
    if not required and not legacy_required:
        return help_text
    description = _LEADING_REQUIRED.sub("", help_text, count=1)
    description = _REQUIRED_PARENTHETICAL.sub("", description)
    return "[bold red]Required:[/bold red] %s" % description.lstrip()


def mark_cli_required_help(
    parser: argparse.ArgumentParser,
    destinations: Iterable[str],
) -> argparse.ArgumentParser:
    """Mark config-validated requirements without changing argparse semantics."""
    required_destinations = set(destinations)
    for action in parser._actions:
        if (
            action.dest in required_destinations
            and action.help is not argparse.SUPPRESS
        ):
            setattr(action, _HELP_REQUIRED_ATTRIBUTE, True)
    return parser


def move_cli_help_actions(
    parser: argparse.ArgumentParser,
    group: argparse._ArgumentGroup,
    destinations: Iterable[str],
) -> None:
    """Move existing actions between help groups without changing parsing."""
    actions = [
        action
        for destination in destinations
        for action in parser._actions
        if action.dest == destination
    ]
    for source_group in parser._action_groups:
        source_group._group_actions[:] = [
            action
            for action in source_group._group_actions
            if action not in actions
        ]
    group._group_actions.extend(actions)


def cli_option_values(
    arguments: Sequence[str], option: str,
) -> tuple[str, ...]:
    """Return values explicitly supplied to one long CLI option.

    This deliberately performs no default resolution.  It is used only to
    route context-sensitive help before the complete parser validates the
    command.  Values following an ``nargs`` option are collected until the
    next option token; ``--option=value`` is supported as well.
    """
    for index, argument in enumerate(arguments):
        if argument == option:
            values = []
            for value in arguments[index + 1:]:
                if value.startswith("-"):
                    break
                values.append(value)
            return tuple(values)
        if argument.startswith(option + "="):
            first = argument.partition("=")[2]
            values = [first] if first else []
            for value in arguments[index + 1:]:
                if value.startswith("-"):
                    break
                values.append(value)
            return tuple(values)
    return ()


def cli_option_value(arguments: Sequence[str], option: str) -> str | None:
    """Return one explicitly supplied long-option value, without a default."""
    values = cli_option_values(arguments, option)
    return values[0] if values else None


def format_cli_examples(
    *examples: tuple[str, str, Sequence[str]],
    notes: Sequence[str] = (),
    title: str = "Examples",
) -> str:
    """Render consistently indented, copyable, colour-coded CLI examples."""
    lines = ["[bold cyan]%s[/bold cyan]" % escape(title)]
    for label, command, arguments in examples:
        lines.extend((
            "",
            "[bold bright_yellow]%s[/bold bright_yellow]" % escape(label),
        ))
        arguments = tuple(arguments)
        lines.append(
            "  [bold magenta]%s[/bold magenta]%s"
            % (escape(command), " \\" if arguments else "")
        )
        for index, argument in enumerate(arguments):
            continuation = " \\" if index < len(arguments) - 1 else ""
            lines.append(
                "    [cyan]%s[/cyan]%s" % (escape(argument), continuation)
            )
    if notes:
        lines.extend(("", "[bold cyan]Notes[/bold cyan]"))
        lines.extend("  [dim]%s[/dim]" % escape(note) for note in notes)
    return "\n".join(lines)


def format_cli_help_block(text: str) -> str:
    """Mark deliberate multi-line help prose for width-aware shared wrapping."""
    return _WRAPPED_HELP_BLOCK_PREFIX + text


__all__ = [
    "AlignedRichHelpFormatter",
    "cli_option_value",
    "cli_option_values",
    "format_cli_choices",
    "format_cli_default",
    "format_cli_examples",
    "format_cli_help_block",
    "help_with_choices",
    "help_with_conditional_requirement",
    "help_with_default",
    "mark_cli_required_help",
    "move_cli_help_actions",
    "style_cli_defaults",
    "style_cli_requirement",
]
