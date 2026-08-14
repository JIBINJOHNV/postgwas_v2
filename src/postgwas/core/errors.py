"""Typed errors shared by CLIs, modules, and pipeline orchestration."""


class PostGWASError(Exception):
    """Base class for expected PostGWAS failures."""


class ConfigurationError(PostGWASError):
    """The resolved user configuration is invalid."""


class MissingRequiredArgumentsError(ConfigurationError):
    """Resolved CLI and YAML settings omit one or more required values."""


class PipelinePlanningError(ConfigurationError):
    """The requested module graph cannot be constructed safely."""


class ModuleExecutionError(PostGWASError):
    """A module failed after pipeline planning completed."""

    def __init__(self, module_name, message):
        self.module_name = str(module_name)
        super().__init__("%s: %s" % (self.module_name, message))
