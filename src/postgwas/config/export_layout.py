"""Validate presentation metadata without adding analysis configuration keys."""

from __future__ import annotations

from typing import Any

from pydantic import Field, ValidationError, field_validator, model_validator

from postgwas.config.loader import _packaged_yaml, select_configuration_values
from postgwas.config.models.common import StrictModel
from postgwas.core.errors import ConfigurationError


class ExportSection(StrictModel):
    title: str = Field(min_length=1)
    fields: list[str] = Field(min_length=1)


class ExportLayout(StrictModel):
    introduction: list[str] = Field(min_length=1)
    common_introduction: list[str] = Field(min_length=1)
    sections: list[ExportSection] = Field(min_length=1)
    common_fields: list[str] = Field(min_length=1)
    common_run_fields: list[str] = Field(min_length=1)
    comments: dict[str, list[str]]
    run_comments: dict[str, list[str]]
    policy_summaries: dict[str, str]

    @field_validator("common_fields", "common_run_fields")
    @classmethod
    def unique_paths(cls, values: list[str]) -> list[str]:
        if len(values) != len(set(values)) or any(
            not part.strip() for value in values for part in value.split(".")
        ):
            raise ValueError("export paths must be unique non-empty dotted paths")
        return values

    @model_validator(mode="after")
    def unique_sections(self):
        fields = [field for section in self.sections for field in section.fields]
        if len(fields) != len(set(fields)):
            raise ValueError("each module field must appear in exactly one export section")
        return self


def load_export_layout(module: str, values: dict[str, Any]) -> ExportLayout:
    """Check the canonical export layout against the resolved module fields."""
    document = _packaged_yaml("defaults/modules/%s.yaml" % module)
    try:
        layout = ExportLayout.model_validate(document.get("export_layout"))
    except ValidationError as exc:
        raise ConfigurationError("Invalid %s export layout: %s" % (module, exc)) from exc
    fields = {field for section in layout.sections for field in section.fields}
    if fields != set(values):
        raise ConfigurationError(
            "%s export sections must cover every module field exactly once; "
            "missing=%s, unknown=%s"
            % (module, sorted(set(values) - fields), sorted(fields - set(values)))
        )
    select_configuration_values(values, [*layout.common_fields, *layout.comments])
    if set(layout.policy_summaries) != set(document["policies"]):
        raise ConfigurationError("%s export summaries must cover every policy group" % module)
    common = document["common"]
    known = {"%s.%s" % (group, key) for group, entries in document["policies"].items() for key in entries}
    if len(common) != len(set(common)) or set(common) - known:
        raise ConfigurationError("%s common policies must be unique registered keys" % module)
    return layout


def changed_configuration_paths(values: dict, defaults: dict, prefix: str = "") -> list[str]:
    """Retain non-default values when exporting a smaller configuration.

    Lists and empty mappings are atomic. Unknown children retain their parent
    mapping, which also handles literal dots in user-defined mapping keys.
    """
    changed = []
    for key, value in values.items():
        path = "%s.%s" % (prefix, key) if prefix else key
        if key in defaults and value == defaults[key]:
            continue
        previous = defaults.get(key)
        if (
            isinstance(value, dict) and value
            and isinstance(previous, dict)
            and set(value) <= set(previous)
            and not any("." in str(child) for child in value)
        ):
            changed.extend(changed_configuration_paths(value, previous, path))
        else:
            changed.append(path)
    return changed
