"""Accelerated configuration parsing retains the safe, fresh-read contract."""

from datetime import date
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from postgwas.config import loader
from postgwas.core.errors import ConfigurationError


def _read(tmp_path, monkeypatch, text, packaged):
    path = tmp_path / "fixture.yaml"
    path.write_text(text, encoding="utf-8")
    if packaged:
        monkeypatch.setattr(loader, "files", lambda package: tmp_path)
        return loader._packaged_yaml(path.name)
    return loader._read_yaml(path)


@pytest.mark.skipif(not hasattr(yaml, "CSafeLoader"), reason="PyYAML has no optional LibYAML extension")
@pytest.mark.parametrize("packaged", [False, True])
def test_configuration_uses_available_c_safe_parser_without_python_scan(tmp_path, monkeypatch, packaged):
    def python_scan(*args, **kwargs):
        pytest.fail("available safe C parser was bypassed by the pure Python scan")
    monkeypatch.setattr(yaml.SafeLoader, "get_single_data", python_scan)
    assert _read(tmp_path, monkeypatch, "execution: {threads: 3}\n", packaged) == {
        "execution": {"threads": 3}}


@pytest.mark.parametrize("packaged", [False, True])
def test_safe_configuration_parser_preserves_yaml_types_aliases_order_and_unicode(tmp_path, monkeypatch, packaged):
    text = ("first: &base {enabled: true, n: 0x10, missing: null}\n"
            "second: {<<: *base, n: 2}\n"
            "label: 'off'\n"
            "unicode: 'βeta — 遺伝子'\n"
            "day: 2026-11-05\n"
            "block: |\n  first line\n  second line\n"
            "values: [0, 1.0, .inf, -.inf, .nan]\n")
    actual = _read(tmp_path, monkeypatch, text, packaged)
    assert list(actual) == ["first", "second", "label", "unicode", "day", "block", "values"]
    assert actual["first"] == {"enabled": True, "n": 16, "missing": None}
    assert actual["second"] == {"enabled": True, "n": 2, "missing": None}
    assert actual["label"] == "off" and actual["unicode"] == "βeta — 遺伝子"
    assert actual["day"] == date(2026, 11, 5)
    assert actual["block"] == "first line\nsecond line\n"
    assert actual["values"][:4] == [0, 1.0, math.inf, -math.inf]
    assert math.isnan(actual["values"][4])


@pytest.mark.parametrize("packaged", [False, True])
@pytest.mark.parametrize("text", ["", "# an empty document\n", "null\n"])
def test_empty_configuration_documents_remain_empty_mappings(tmp_path, monkeypatch, packaged, text):
    assert _read(tmp_path, monkeypatch, text, packaged) == {}


@pytest.mark.parametrize("packaged", [False, True])
@pytest.mark.parametrize("text", [
    "execution: [unterminated\n",
    "key: !!python/object/new:builtins.list [1, 2]\n",
    "key: !!python/name:builtins.len ''\n",
    "key: !unsupported value\n",
    "---\na: 1\n---\nb: 2\n",
])
def test_invalid_and_nonstandard_tagged_configuration_is_rejected_with_source(tmp_path, monkeypatch, packaged, text):
    with pytest.raises(ConfigurationError) as caught:
        _read(tmp_path, monkeypatch, text, packaged)
    assert "fixture.yaml" in str(caught.value)
    assert isinstance(caught.value.__cause__, yaml.YAMLError)


@pytest.mark.parametrize("packaged", [False, True])
@pytest.mark.parametrize("text", ["[one, two]\n", "a scalar\n"])
def test_nonmapping_configuration_roots_remain_invalid(tmp_path, monkeypatch, packaged, text):
    with pytest.raises(ConfigurationError, match="must be a mapping"):
        _read(tmp_path, monkeypatch, text, packaged)


def test_all_packaged_default_documents_match_standard_safe_parser():
    paths = ["defaults/application.yaml", "defaults/pipeline.yaml", "defaults/resources.yaml"]
    paths.extend("defaults/modules/%s.yaml" % name for name in loader.MODULE_NAMES)
    assert len(paths) == 25
    for path in paths:
        text = loader.files("postgwas.config").joinpath(path).read_text(encoding="utf-8")
        expected = yaml.safe_load(text) or {}
        actual = loader._packaged_yaml(path)
        assert actual == expected, path
        assert list(actual) == list(expected), path


def test_complete_validated_defaults_match_standard_safe_parser(monkeypatch):
    original = yaml.load
    with monkeypatch.context() as patch:
        patch.setattr(yaml, "load", lambda stream, Loader: original(stream, Loader=yaml.SafeLoader))
        expected = loader.load_configuration().model_dump(mode="json")
    assert loader.load_configuration().model_dump(mode="json") == expected


def test_configuration_still_loads_without_optional_c_extension(tmp_path):
    path = tmp_path / "run.yaml"
    path.write_text("execution: {threads: 3}\n", encoding="utf-8")
    unsafe = tmp_path / "unsafe.yaml"
    unsafe.write_text("value: !!python/object/new:builtins.list [1, 2]\n", encoding="utf-8")
    script = ("import json, sys, yaml\n"
              "if hasattr(yaml, 'CSafeLoader'): del yaml.CSafeLoader\n"
              "from postgwas.config import load_configuration\n"
              "from postgwas.core.errors import ConfigurationError\n"
              "value = load_configuration(sys.argv[1])\n"
              "try: load_configuration(sys.argv[2])\n"
              "except ConfigurationError: unsafe_rejected = True\n"
              "else: unsafe_rejected = False\n"
              "print(json.dumps({'threads': value.execution.threads, 'resume': value.run.resume, 'unsafe_rejected': unsafe_rejected}))\n")
    completed = subprocess.run([sys.executable, "-B", "-c", script, str(path), str(unsafe)],
        capture_output=True, text=True, check=False,
        env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"), timeout=60)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout) == {"threads": 3, "resume": True, "unsafe_rejected": True}


def test_loaded_default_documents_and_models_do_not_share_mutable_values():
    first = loader.packaged_defaults()
    expected_template = first["modules"]["harmonisation"]["output_layout"]["chromosome_raw_vcf"]
    first["modules"]["harmonisation"]["output_layout"]["chromosome_raw_vcf"] = "mutated.vcf"
    assert loader.packaged_defaults()["modules"]["harmonisation"]["output_layout"]["chromosome_raw_vcf"] == expected_template
    configuration = loader.load_configuration()
    expected_modules = list(configuration.pipeline.modules)
    configuration.pipeline.modules.append("not_a_real_module")
    assert loader.load_configuration().pipeline.modules == expected_modules


def test_changed_run_and_include_files_are_reread_even_with_same_size_and_mtime(tmp_path):
    included = tmp_path / "base.yaml"
    included.write_text("execution: {threads: 3}\n", encoding="utf-8")
    run = tmp_path / "run.yaml"
    run.write_text("include: base.yaml\nlogging: {show_screen: true}\n", encoding="utf-8")
    before = included.stat()
    first = loader.load_configuration(run)
    included.write_text("execution: {threads: 7}\n", encoding="utf-8")
    os.utime(included, ns=(before.st_atime_ns, before.st_mtime_ns))
    run.write_text("include: base.yaml\nlogging: {show_screen: false}\n", encoding="utf-8")
    second = loader.load_configuration(run)
    assert first.execution.threads == 3 and first.logging.show_screen is True
    assert second.execution.threads == 7 and second.logging.show_screen is False
