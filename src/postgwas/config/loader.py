"""Load, layer, validate, and record PostGWAS configuration."""

from __future__ import annotations

from copy import deepcopy
from importlib.resources import files
import math
import os
from pathlib import Path
from typing import Any, Iterable, Mapping

import yaml
from pydantic import ValidationError

from postgwas.config.merger import apply_dotted_overrides, deep_merge
from postgwas.config.models import PostGWASConfig
from postgwas.config.models.application import ModulesConfig
from postgwas.config.validator import configuration_error
from postgwas.core.errors import ConfigurationError
from postgwas.core.io.reports import write_yaml_report


MODULE_NAMES = tuple(ModulesConfig.model_fields)
CONFIG_ALIASES = {
    "sumstat_filter": "filtering",
    "post_imputation_filter": "filtering",
    "annot_ldblock": "ld_annotation",
    "formatter": "formatting",
    "finemap": "fine_mapping",
    "ld_clump": "ld_clumping",
    "heritability": "ldsc",
    "qc": "qc_summary",
    "pathway_enrichment": "enrichment",
}
_ROOT_OVERRIDE_PREFIXES = (
    "run.", "execution.", "logging.", "resources.", "pipeline.", "modules.",
)


def canonical_module_name(name: str) -> str:
    canonical = CONFIG_ALIASES.get(name, name)
    if canonical not in MODULE_NAMES:
        raise ConfigurationError(
            "Unknown configuration module %r. Available modules: %s"
            % (name, ", ".join(MODULE_NAMES))
        )
    return canonical


def _module_cli_overrides(
    name: str, overrides: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Qualify module-local keys while preserving explicit root dotted paths."""
    return {
        (
            key
            if key.startswith(_ROOT_OVERRIDE_PREFIXES)
            else "modules.%s.%s" % (name, key)
        ): value
        for key, value in (overrides or {}).items()
    }


def _read_yaml(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = yaml.safe_load(handle) or {}
    except OSError as exc:
        raise ConfigurationError("Cannot read configuration %s: %s" % (path, exc)) from exc
    except yaml.YAMLError as exc:
        raise ConfigurationError("Invalid YAML in %s: %s" % (path, exc)) from exc
    if not isinstance(value, dict):
        raise ConfigurationError("Configuration root must be a mapping: %s" % path)
    return value


def _packaged_yaml(relative_path: str) -> dict[str, Any]:
    resource = files("postgwas.config").joinpath(relative_path)
    try:
        value = yaml.safe_load(resource.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigurationError(
            "Cannot load packaged configuration %s: %s" % (relative_path, exc)
        ) from exc
    if not isinstance(value, dict):
        raise ConfigurationError("Packaged configuration must be a mapping: %s" % relative_path)
    return value


def _packaged_module_defaults(name: str) -> dict[str, Any]:
    """Load a module's defaults from its canonical configuration document."""
    document = _packaged_yaml("defaults/modules/%s.yaml" % name)
    module = document.get("module", document)
    if not isinstance(module, dict):
        raise ConfigurationError(
            "Packaged module configuration must be a mapping: %s" % name
        )
    return module


def packaged_defaults(profile: str | None = None) -> dict[str, Any]:
    configuration = _packaged_yaml("defaults/application.yaml")
    configuration = deep_merge(configuration, {"pipeline": _packaged_yaml("defaults/pipeline.yaml")})
    configuration = deep_merge(configuration, {"resources": _packaged_yaml("defaults/resources.yaml")})
    modules = {
        name: _packaged_module_defaults(name)
        for name in MODULE_NAMES
    }
    configuration = deep_merge(configuration, {"modules": modules})
    if profile:
        configuration = deep_merge(
            configuration, _packaged_yaml("profiles/%s.yaml" % profile)
        )
    return configuration


def detected_logical_threads() -> int:
    return max(1, int(os.cpu_count() or 1))


def detected_physical_memory_gb() -> float:
    try:
        import psutil

        total_bytes = int(psutil.virtual_memory().total)
    except ImportError:
        try:
            total_bytes = int(os.sysconf("SC_PAGE_SIZE")) * int(os.sysconf("SC_PHYS_PAGES"))
        except (AttributeError, OSError, ValueError):
            total_bytes = 1024**3
    return total_bytes / float(1024**3)


def resolve_execution_defaults(configuration: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve hardware-derived values after all configuration layers merge."""
    resolved = deep_merge({}, configuration)
    execution = resolved.get("execution")
    if not isinstance(execution, dict):
        raise ConfigurationError("execution configuration must be a mapping")
    try:
        reserve_threads = int(execution["reserve_threads"])
        usable_fraction = float(execution["usable_memory_fraction"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ConfigurationError(
            "execution.reserve_threads and execution.usable_memory_fraction are required"
        ) from exc
    threads = execution.get("threads")
    if isinstance(threads, str) and threads.strip().lower() == "auto":
        execution["threads"] = max(1, detected_logical_threads() - reserve_threads)
    memory = execution.get("memory_gb")
    if isinstance(memory, str) and memory.strip().lower() == "auto":
        execution["memory_gb"] = max(
            1, math.floor(detected_physical_memory_gb() * usable_fraction)
        )
    return resolved


def _load_user_layers(path: Path, seen: set[Path] | None = None) -> dict[str, Any]:
    """Load optional ``include`` files before applying the current file."""
    path = path.expanduser().resolve()
    seen = set() if seen is None else seen
    if path in seen:
        raise ConfigurationError("Configuration include cycle detected at %s" % path)
    seen.add(path)
    current = _read_yaml(path)
    include_value = current.pop("include", [])
    includes = [include_value] if isinstance(include_value, str) else include_value
    if not isinstance(includes, list) or not all(isinstance(item, str) for item in includes):
        raise ConfigurationError("'include' must be a path or list of paths in %s" % path)
    merged: dict[str, Any] = {}
    for included in includes:
        included_path = (path.parent / included).resolve()
        merged = deep_merge(merged, _load_user_layers(included_path, seen))
    seen.remove(path)
    return deep_merge(merged, current)


def load_configuration(
    config_file: str | Path | None = None,
    *,
    profile: str | None = None,
    cli_overrides: Mapping[str, Any] | None = None,
) -> PostGWASConfig:
    """Resolve ``defaults < profile < user YAML < explicit CLI``."""
    resolved = packaged_defaults(profile)
    source = "packaged defaults"
    if config_file is not None:
        path = Path(config_file)
        resolved = deep_merge(resolved, _load_user_layers(path))
        source = str(path)
    resolved = apply_dotted_overrides(resolved, cli_overrides)
    resolved = resolve_execution_defaults(resolved)
    try:
        return PostGWASConfig.model_validate(resolved)
    except ValidationError as exc:
        raise configuration_error(exc, source) from exc


def load_module_configuration(
    module: str,
    config_file: str | Path | None = None,
    *,
    profile: str | None = None,
    cli_overrides: Mapping[str, Any] | None = None,
):
    """Resolve one module from either a full pipeline or module-only YAML file."""
    name = canonical_module_name(module)
    module_overrides = _module_cli_overrides(name, cli_overrides)
    if config_file is None:
        root = load_configuration(profile=profile, cli_overrides=module_overrides)
        return getattr(root.modules, name)

    path = Path(config_file)
    supplied = _load_user_layers(path)
    if not any(key in supplied for key in ("config_version", "modules", "pipeline")):
        supplied = {"modules": {name: supplied}}
    defaults = packaged_defaults(profile)
    resolved = deep_merge(defaults, supplied)
    resolved = apply_dotted_overrides(resolved, module_overrides)
    resolved = resolve_execution_defaults(resolved)
    try:
        root = PostGWASConfig.model_validate(resolved)
    except ValidationError as exc:
        raise configuration_error(exc, str(path)) from exc
    return getattr(root.modules, name)


def load_run_configuration_for_module(
    module: str,
    config_file: str | Path | None = None,
    *,
    module_overrides: Mapping[str, Any] | None = None,
    global_overrides: Mapping[str, Any] | None = None,
) -> PostGWASConfig:
    """Resolve a full run from either full-run or module-only YAML.

    Module services share this boundary so module-only files cannot cause
    duplicated fallback logic or hide an invalid full-run configuration.
    """
    name = canonical_module_name(module)
    local = _module_cli_overrides(name, module_overrides)
    global_values = dict(global_overrides or {})
    if config_file is None:
        overrides = dict(global_values)
        overrides.update(local)
        return load_configuration(cli_overrides=overrides)

    path = Path(config_file)
    supplied = _load_user_layers(path)
    module_only = not any(
        key in supplied for key in ("config_version", "modules", "pipeline")
    )
    if not module_only:
        overrides = dict(global_values)
        overrides.update(local)
        return load_configuration(path, cli_overrides=overrides)

    root_overrides = {
        key: value for key, value in local.items()
        if not key.startswith("modules.")
    }
    root_overrides.update(global_values)
    module_values = {
        key: value for key, value in local.items()
        if key.startswith("modules.")
    }
    root = load_configuration(cli_overrides=root_overrides)
    setattr(
        root.modules,
        name,
        load_module_configuration(name, path, cli_overrides=module_values),
    )
    return root


def select_configuration_values(
    source: Mapping[str, Any], dotted_paths: str | Iterable[str],
) -> dict[str, Any]:
    """Copy selected configuration leaves while retaining their YAML hierarchy."""
    selected: dict[str, Any] = {}
    paths = (dotted_paths,) if isinstance(dotted_paths, str) else dotted_paths
    for dotted_path in dict.fromkeys(paths):
        parts = dotted_path.split(".")
        if not parts or any(not part for part in parts):
            raise ConfigurationError(
                "Invalid resolved-configuration path: %r" % dotted_path
            )
        source_cursor: Any = source
        for part in parts:
            if not isinstance(source_cursor, Mapping) or part not in source_cursor:
                raise ConfigurationError(
                    "Unknown resolved-configuration path: %s" % dotted_path
                )
            source_cursor = source_cursor[part]
        destination_cursor = selected
        for part in parts[:-1]:
            destination_cursor = destination_cursor.setdefault(part, {})
        destination_cursor[parts[-1]] = deepcopy(source_cursor)
    return selected


def resolved_configuration_values(
    config: PostGWASConfig,
    *,
    modules: str | Iterable[str] | None = None,
    resource_paths: str | Iterable[str] = (),
    module_paths: Mapping[str, str | Iterable[str]] | None = None,
) -> dict[str, Any]:
    """Return a reloadable full or module-scoped resolved configuration.

    Run metadata is scoped to the modules that actually consumed the values.
    A file without ``modules`` remains the complete application configuration
    used by ``postgwas config show --output``.
    """
    values = config.model_dump(mode="json", exclude_none=False)
    if modules is None:
        if module_paths:
            raise ConfigurationError(
                "module_paths requires a module-scoped resolved configuration"
            )
        return values

    module_names = (modules,) if isinstance(modules, str) else modules
    selected_modules = tuple(dict.fromkeys(
        canonical_module_name(name) for name in module_names
    ))
    if not selected_modules:
        raise ConfigurationError(
            "At least one module is required for scoped run metadata"
        )
    selected_paths = {
        canonical_module_name(name): paths
        for name, paths in (module_paths or {}).items()
    }
    unknown_path_modules = sorted(set(selected_paths) - set(selected_modules))
    if unknown_path_modules:
        raise ConfigurationError(
            "Resolved-configuration paths were supplied for unselected modules: %s"
            % ", ".join(unknown_path_modules)
        )
    scoped = {
        "config_version": values["config_version"],
        "run": values["run"],
        "execution": values["execution"],
        "logging": values["logging"],
    }
    resources = select_configuration_values(values["resources"], resource_paths)
    if resources:
        scoped["resources"] = resources
    scoped["modules"] = {}
    for name in selected_modules:
        paths = selected_paths.get(name)
        scoped["modules"][name] = (
            deepcopy(values["modules"][name])
            if paths is None
            else select_configuration_values(values["modules"][name], paths)
        )
    return scoped


def write_resolved_configuration(
    config: PostGWASConfig,
    output_file: str | Path,
    *,
    modules: str | Iterable[str] | None = None,
    resource_paths: str | Iterable[str] = (),
    module_paths: Mapping[str, str | Iterable[str]] | None = None,
) -> Path:
    """Atomically write reloadable resolved configuration or scoped run metadata."""
    serialised = resolved_configuration_values(
        config,
        modules=modules,
        resource_paths=resource_paths,
        module_paths=module_paths,
    )
    return write_yaml_report(serialised, output_file)
