"""Validation error formatting shared by every execution mode."""

from pydantic import ValidationError

from postgwas.core.errors import ConfigurationError


def configuration_error(error: ValidationError, source: str) -> ConfigurationError:
    lines = ["Invalid PostGWAS configuration (%s):" % source]
    for item in error.errors(include_url=False):
        location = ".".join(str(part) for part in item["loc"])
        lines.append("  - %s: %s" % (location or "configuration", item["msg"]))
    return ConfigurationError("\n".join(lines))
