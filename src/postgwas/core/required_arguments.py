"""Shared validation for requirements resolved from CLI overrides and YAML."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from postgwas.core.errors import MissingRequiredArgumentsError


@dataclass(frozen=True)
class RequiredArgument:
    """One public option and its canonical resolved configuration value."""

    option: str
    configuration_path: str | None
    value: object


@dataclass(frozen=True)
class RequiredAlternative:
    """A requirement satisfied by any one of several public options."""

    arguments: tuple[RequiredArgument, ...]


def _missing(value: object) -> bool:
    return value is None or (isinstance(value, str) and not value.strip())


def require_resolved_arguments(
    requirements: Iterable[RequiredArgument],
    *,
    alternatives: Sequence[RequiredAlternative] = (),
) -> None:
    """Report all independently missing resolved requirements in one error."""
    messages = []
    for requirement in requirements:
        if not _missing(requirement.value):
            continue
        if requirement.configuration_path:
            message = (
                "Required argument not provided: %s. Provide %s VALUE or set %s "
                "in the run configuration."
                % (
                    requirement.option,
                    requirement.option,
                    requirement.configuration_path,
                )
            )
        else:
            message = (
                "Required argument not provided: %s. Provide %s VALUE."
                % (requirement.option, requirement.option)
            )
        if message not in messages:
            messages.append(message)
    for alternative in alternatives:
        if not alternative.arguments:
            raise ValueError("RequiredAlternative must contain at least one argument")
        if all(_missing(argument.value) for argument in alternative.arguments):
            options = " or ".join(argument.option for argument in alternative.arguments)
            configuration_paths = " or ".join(
                argument.configuration_path
                for argument in alternative.arguments
                if argument.configuration_path
            )
            message = (
                "Required argument not provided: %s. Provide one of these options "
                "with a VALUE%s."
                % (
                    options,
                    (
                        " or set %s in the run configuration"
                        % configuration_paths
                        if configuration_paths else ""
                    ),
                )
            )
            if message not in messages:
                messages.append(message)
    if messages:
        raise MissingRequiredArgumentsError("\n".join(messages))


__all__ = [
    "RequiredAlternative",
    "RequiredArgument",
    "require_resolved_arguments",
]
