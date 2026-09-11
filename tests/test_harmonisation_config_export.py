"""User-facing exports must not change resolved analysis settings."""

from copy import deepcopy
import sys

import pytest
import yaml

from postgwas.config import load_configuration, load_module_configuration, load_run_configuration_for_module
from postgwas.config import export_layout
from postgwas.config.cli import main
from postgwas.config.exporter import render_module_configuration
from postgwas.config.loader import _packaged_yaml
from postgwas.core.errors import ConfigurationError
from postgwas.modules.harmonisation.policies import _yaml_scalar, default_policies, load_policies


def _resolved(config):
    values = config.model_dump(mode="json")
    module = values["modules"]["harmonisation"]
    module["policies"] = load_policies(module["policies"]).resolved_dict()
    return values


@pytest.mark.parametrize("style", ["full", "minimal", "values"])
def test_common_export_reloads_identical_defaults(tmp_path, style):
    text = render_module_configuration("harmonisation", scope="common", style=style)
    destination = tmp_path / "common.yaml"
    destination.write_text(text)
    document = yaml.safe_load(text)
    module = document["modules"]["harmonisation"]
    assert "resource_layout" not in module
    assert "output_layout" not in module
    assert len(module["policies"]) < len(default_policies().nested_dict())
    assert _resolved(load_configuration(destination)) == _resolved(load_configuration())
    assert document["run"]["resume"] is True
    assert document["run"]["overwrite"] is False
    if style == "values":
        assert not any(line.lstrip().startswith("#") for line in text.splitlines())
    else:
        assert "without removing additional variants" in text
        assert "shared policies appear once" in text


@pytest.mark.parametrize("style", ["full", "minimal", "values"])
def test_common_preserves_included_advanced_and_global_overrides(tmp_path, style):
    included = tmp_path / "included.yaml"
    included.write_text(yaml.safe_dump({
        "modules": {"harmonisation": {
            "concordance_validation": {"staging_batch_rows": 1234},
            "policies": {"eaf": {"missing_fraction_cutoff": 0.025}},
        }},
    }))
    source = tmp_path / "input.yaml"
    supplied = {
        "include": [included.name],
        "run": {"resume": False, "overwrite": True, "output_directory": "results with spaces"},
        "execution": {"threads": 2, "memory_gb": 4},
        "logging": {"show_screen": False},
        "resources": {"root": "/references/custom panel"},
        "modules": {
            "harmonisation": {
                "default_eaf": {"source": "1000G"},
                "policies": {
                    "strand": {"unmatched_action": "fail"},
                    "chromosome": {"allowed_after_split": ["1", "2", "X"]},
                    "effect_from_z": {"phenotype_standard_deviation": None},
                },
            },
            "qc_summary": {"rules": {"maf_min": 0.02, "remove_mhc": True}},
            "formatting": {"formats": ["magma", "ldsc"]},
        },
    }
    source.write_text(yaml.safe_dump(supplied))
    original = source.read_bytes()
    expected = _resolved(load_run_configuration_for_module("harmonisation", source))
    text = render_module_configuration("harmonisation", config_file=source, scope="common", style=style)
    destination = tmp_path / "export.yaml"
    destination.write_text(text)
    assert source.read_bytes() == original
    assert _resolved(load_configuration(destination)) == expected
    assert "batch_rows: 1234" in text
    assert "missing_fraction_cutoff:" in text
    assert "include:" not in text


def test_module_only_input_and_cli_override_still_work(tmp_path):
    source = tmp_path / "module.yaml"
    source.write_text("policies:\n  strand:\n    unmatched_action: retain\n")
    destination = tmp_path / "common.yaml"
    destination.write_text(render_module_configuration("harmonisation", config_file=source, scope="common"))
    config = load_run_configuration_for_module(
        "harmonisation", destination,
        module_overrides={"policies.strand.unmatched_action": "fail"},
        global_overrides={"run.output_directory": "cli-results"},
    )
    assert load_policies(config.modules.harmonisation.policies).strand.unmatched_action == "fail"
    assert str(config.run.output_directory) == "cli-results"


def test_all_export_keeps_every_setting():
    expected = load_module_configuration("harmonisation").model_dump(mode="json")
    expected["policies"] = default_policies().nested_dict()
    for style in ("full", "minimal", "values"):
        document = yaml.safe_load(render_module_configuration("harmonisation", style=style))
        assert document == expected
        keys = list(document)
        assert keys.index("policies") < keys.index("resource_layout")
        policy_groups = list(document["policies"])
        assert policy_groups.index("info") < policy_groups.index("duplicates")
        assert policy_groups.index("sample_size") < policy_groups.index("build")


@pytest.mark.parametrize("shorthand", [False, True])
def test_cli_common_export_and_validation(tmp_path, monkeypatch, capsys, shorthand):
    destination = tmp_path / "common.yaml"
    monkeypatch.setattr(sys, "argv", ["postgwas --config" if shorthand else "postgwas config"]
                        + ([] if shorthand else ["export"])
                        + ["--module", "harmonisation", "--scope", "common", "--output", str(destination)])
    assert main() == 0
    monkeypatch.setattr(sys, "argv", ["postgwas config", "validate", "--config", str(destination)])
    assert main() == 0
    assert "Configuration is valid" in capsys.readouterr().out


@pytest.mark.parametrize("target", [["--module", "magma"], ["--pipeline", "magma"]])
def test_common_scope_rejects_unsupported_targets(tmp_path, monkeypatch, capsys, target):
    destination = tmp_path / "must-not-exist.yaml"
    monkeypatch.setattr(sys, "argv", ["postgwas config", "export", *target,
                                     "--scope", "common", "--output", str(destination)])
    assert main() == 2
    assert "requires --module harmonisation" in capsys.readouterr().err
    assert not destination.exists()


def test_invalid_source_cannot_produce_short_config(tmp_path):
    source = tmp_path / "invalid.yaml"
    source.write_text("policies:\n  strand:\n    unmatched_action: silently_guess\n")
    with pytest.raises(ConfigurationError):
        render_module_configuration("harmonisation", config_file=source, scope="common")
    with pytest.raises(ConfigurationError):
        render_module_configuration("harmonisation", config_file=tmp_path / "missing.yaml", scope="common")


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "unknown", "common", "summary", "run_comment"])
def test_layout_rejects_stale_or_incomplete_paths(monkeypatch, mutation):
    document = deepcopy(_packaged_yaml("defaults/modules/harmonisation.yaml"))
    layout = document["export_layout"]
    if mutation == "missing":
        layout["sections"][0]["fields"].remove("enabled")
    elif mutation == "duplicate":
        layout["sections"][1]["fields"].append("enabled")
    elif mutation == "unknown":
        layout["common_fields"].append("nonexistent.path")
    elif mutation == "common":
        document["common"].append("strand.nonexistent")
    elif mutation == "summary":
        layout["policy_summaries"].pop("strand")
    else:
        layout["run_comments"]["run.nonexistent"] = ["Invalid path"]
    monkeypatch.setattr(export_layout, "_packaged_yaml", lambda _: document)
    with pytest.raises(ConfigurationError):
        render_module_configuration("harmonisation", scope="common")


def test_changed_paths_preserve_empty_lists_nulls_and_literal_dotted_keys():
    current = {"a": {"items": [], "nullable": None}, "b": {"literal.key": "new"}}
    baseline = {"a": {"items": ["1"], "nullable": 1}, "b": {"literal.key": "old"}}
    assert export_layout.changed_configuration_paths(current, baseline) == ["a.items", "a.nullable", "b"]


@pytest.mark.parametrize("value", ['a "quoted" name', r'(?P<chr>\d+)', 'a\nb', '2026-09-07', 'off', 1e-300, 1e20,
                                  {"quoted:key": ["X", "01", None, False]}])
def test_policy_export_preserves_yaml_sensitive_values(value):
    assert yaml.safe_load(_yaml_scalar(value)) == value
