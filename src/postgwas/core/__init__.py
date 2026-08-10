"""Shared architectural contracts for PostGWAS.

Scientific modules should depend on this package, while this package must not
import any scientific module.  Keeping that dependency direction makes module
code usable from both the CLI and the pipeline executor.
"""

from .contracts import Artifact, ModuleResult, RunContext
from .errors import (
    ConfigurationError,
    ModuleExecutionError,
    PipelinePlanningError,
    PostGWASError,
)

__all__ = [
    "Artifact",
    "ConfigurationError",
    "ModuleExecutionError",
    "ModuleResult",
    "PipelinePlanningError",
    "PostGWASError",
    "RunContext",
]
