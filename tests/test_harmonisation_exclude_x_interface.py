"""Explicit chromosome-X scope uses canonical configuration and checkpoints."""

import argparse
import re
import shlex
import sys
from types import SimpleNamespace

import pytest
import yaml

from postgwas.config import load_run_configuration_for_module
from postgwas.config.exporter import render_module_configuration
from postgwas.core.completion import configuration_digest
from postgwas.core.direct_execution import _checkpoint_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.modules.harmonisation import cli
from postgwas.modules.harmonisation.policies import POLICIES, PolicyError, load_policies


def test_exclude_x_default_preserves_x_and_argparse_owns_no_default():
    args = cli.get_harmonisation_parser().parse_args([])
    config = cli.resolve_harmonisation_configuration([])
    assert not hasattr(args, "exclude_chr_x")
    policy = POLICIES["chromosome.exclude_x"]
    assert policy.type == "bool"
    assert policy.default is False
    assert load_policies(config.modules.harmonisation.policies).chromosome.exclude_x is False
    assert cli._engine_defaults(config)["policies"]["chromosome"]["exclude_x"] is False


def test_exclude_x_flag_is_presence_only_and_overrides_yaml_false(tmp_path):
    run = tmp_path / "run.yaml"
    run.write_text("policies:\n  chromosome:\n    exclude_x: false\n")
    parser = cli.get_harmonisation_parser()
    action = next(action for action in parser._actions if action.dest == "exclude_chr_x")
    assert isinstance(action, argparse._StoreTrueAction)
    assert action.default == argparse.SUPPRESS
    config = cli.resolve_harmonisation_configuration([
        "--run-config", str(run), "--exclude-chr-x",
    ])
    assert load_policies(config.modules.harmonisation.policies).chromosome.exclude_x is True
    assert cli._engine_defaults(config)["policies"]["chromosome"]["exclude_x"] is True
    with pytest.raises(ConfigurationError, match="unrecognized arguments"):
        cli.resolve_harmonisation_configuration(["--exclude-chr-x", "false"])


@pytest.mark.parametrize("nested", [False, True])
def test_exclude_x_yaml_true_survives_omitted_flag(tmp_path, nested):
    run = tmp_path / "run.yaml"
    module = {"policies": {"chromosome": {"exclude_x": True}}}
    run.write_text(yaml.safe_dump({"modules": {"harmonisation": module}} if nested else module))
    config = cli.resolve_harmonisation_configuration(["--run-config", str(run)])
    assert load_policies(config.modules.harmonisation.policies).chromosome.exclude_x is True


@pytest.mark.parametrize("invalid", [None, 1, 0, "autosomes", [], {}])
def test_exclude_x_rejects_non_boolean_policy_values(tmp_path, invalid):
    run = tmp_path / "invalid.yaml"
    run.write_text(yaml.safe_dump({"policies": {"chromosome": {"exclude_x": invalid}}}))
    with pytest.raises(PolicyError, match="chromosome.exclude_x.*must be true or false"):
        cli._engine_defaults(load_run_configuration_for_module("harmonisation", run))


@pytest.mark.parametrize("scope", ["all", "common"])
@pytest.mark.parametrize("style", ["full", "minimal", "values"])
@pytest.mark.parametrize("enabled", [False, True])
def test_exclude_x_export_round_trip_preserves_scope(tmp_path, scope, style, enabled):
    run = tmp_path / "run.yaml"
    run.write_text(yaml.safe_dump({"policies": {"chromosome": {"exclude_x": enabled}}}))
    rendered = render_module_configuration("harmonisation", config_file=run, scope=scope, style=style)
    document = yaml.safe_load(rendered)
    module = document["modules"]["harmonisation"] if scope == "common" else document
    assert module["policies"]["chromosome"]["exclude_x"] is enabled
    exported = tmp_path / "export.yaml"
    exported.write_text(rendered)
    reloaded = load_run_configuration_for_module("harmonisation", exported)
    assert load_policies(reloaded.modules.harmonisation.policies).chromosome.exclude_x is enabled


@pytest.mark.parametrize("width", [80, 120])
def test_exclude_x_help_uses_shared_default_and_alignment(monkeypatch, width):
    monkeypatch.setenv("COLUMNS", str(width))
    monkeypatch.setenv("NO_COLOR", "1")
    parser = cli.get_harmonisation_parser(add_help=True)
    lines = parser.format_help().splitlines()
    start = next(index for index, line in enumerate(lines) if line.lstrip().startswith("--exclude-chr-x"))
    block = []
    for line in lines[start:]:
        if block and (not line.strip() or line.lstrip().startswith("--")):
            break
        block.append(line)
    assert "Default: false" in " ".join(" ".join(block).split())
    assert "Configured default:" not in " ".join(block)
    action = next(action for action in parser._actions if action.dest == "exclude_chr_x")
    surrounding = next(line for line in lines if "Folder where results" in line)
    description_column = surrounding.index("Folder")
    assert block[0].index("Exclude") == description_column
    assert all(len(line) - len(line.lstrip()) == description_column for line in block[1:])
    assert action.default == argparse.SUPPRESS


def test_exclude_x_public_help_has_complete_parseable_example(monkeypatch, capsys):
    monkeypatch.setenv("NO_COLOR", "1")
    monkeypatch.setattr(sys, "argv", ["postgwas harmonisation", "--help"])
    with pytest.raises(SystemExit) as error:
        cli.main()
    assert error.value.code == 0
    text = re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", capsys.readouterr().out)
    match = re.search(r"Exclude chromosome X[^\n]*:\s*(postgwas harmonisation(?:[^\n]*\\\n)*[^\n]*)", text)
    assert match is not None
    arguments = shlex.split(match.group(1).replace("\\\n", " "))[2:]
    args = cli.get_harmonisation_parser().parse_args(arguments)
    assert args.exclude_chr_x is True
    assert args.sample_sheet and args.resource_directory and args.output_directory


def test_exclude_x_changes_resolved_scientific_checkpoint_digest():
    retained = cli.resolve_harmonisation_configuration([])
    excluded = cli.resolve_harmonisation_configuration(["--exclude-chr-x"])
    # Omit raw CLI arguments here so this proves the resolved scientific policy
    # itself affects identity, including when configured only through YAML.
    first = _checkpoint_configuration(retained, "harmonisation", "harmonisation", ())
    second = _checkpoint_configuration(excluded, "harmonisation", "harmonisation", ())
    assert configuration_digest(first) != configuration_digest(second)


def test_run_summary_preserves_explicit_exclusion_count_without_changing_rejection_total(tmp_path):
    row = SimpleNamespace(dataset_id="study", input_file=tmp_path / "study.tsv")
    record = cli._run_summary_record(
        row, status="FAILED", manifest_path=tmp_path / "manifest.json",
        manifest={
            "pre_vcf_summary": {"total_variant_removed_excluded_chromosomes": 5},
            "dataset": {"reconciliation": {"rows_read": 20, "rejected": 8}},
        },
        strand_reference_panel="1000G", strand_reference_population="EUR",
        comparison_af_panel="ALFA", comparison_af_population="EUR",
    )
    assert "excluded_chromosome_variants" in cli._RUN_SUMMARY_FIELDS
    assert record["excluded_chromosome_variants"] == 5
    assert record["parsed_variants"] == 20
    assert record["rejected_variants"] == 8


@pytest.mark.parametrize("literal", [False, True])
def test_exclude_x_checkpoint_inputs_preserve_shared_files_and_omit_x_templates(tmp_path, literal):
    panel_x = tmp_path / "panel_GRCh37_X.tsv"
    panel_one = tmp_path / "panel_GRCh37_1.tsv"
    for panel in (panel_x, panel_one):
        panel.write_text("CHROM\tPOS\tREF\tALT\tEUR\n1\t10\tA\tG\t0.2\n")
    specification = str(panel_x) if literal else str(tmp_path / "panel_{build}_{chromosome}.tsv")
    arguments = ["--eaf-replacement-file", specification, "--eaf-replacement-column", "EUR"]
    retained_config = cli.resolve_harmonisation_configuration(arguments)
    retained = cli.harmonisation_checkpoint_inputs(arguments, retained_config)
    excluded_arguments = [*arguments, "--exclude-chr-x"]
    excluded_config = cli.resolve_harmonisation_configuration(excluded_arguments)
    excluded = cli.harmonisation_checkpoint_inputs(excluded_arguments, excluded_config)
    assert set(retained.values()) == ({panel_x} if literal else {panel_one, panel_x})
    assert set(excluded.values()) == ({panel_x} if literal else {panel_one})


@pytest.mark.parametrize("initially_excluded", [False, True])
def test_public_direct_checkpoint_restarts_changed_scope_and_resumes_identical_scope(
    tmp_path, monkeypatch, initially_excluded,
):
    from postgwas import __main__ as entry

    output = tmp_path / "results"
    resource = tmp_path / "resources"
    resource.mkdir()
    arguments = ["--sample-sheet", str(tmp_path / "unused.csv"),
                 "--output-directory", str(output), "--resource-directory", str(resource)]
    calls = []

    def operation():
        config = cli.resolve_harmonisation_configuration(arguments)
        scope = load_policies(config.modules.harmonisation.policies).chromosome.exclude_x
        calls.append(scope)
        output.mkdir(exist_ok=True)
        (output / "scientific_output.txt").write_text(str(scope))
        return 0

    real_resolve = entry.resolve_reference
    monkeypatch.setattr(entry, "resolve_reference", lambda reference: operation
                        if reference == "postgwas.modules.harmonisation.cli:main"
                        else real_resolve(reference))
    monkeypatch.setattr(entry, "_run_recorded_command", lambda operation, settings: operation())
    monkeypatch.setattr(entry, "run_with_progress", lambda operation, **kwargs: operation())

    def invoke():
        monkeypatch.setattr(sys, "argv", ["postgwas", "harmonisation", *arguments])
        return entry.main()

    if initially_excluded:
        arguments.append("--exclude-chr-x")
    assert invoke() == 0
    assert invoke() == 0
    assert calls == [initially_excluded]
    if initially_excluded:
        arguments.remove("--exclude-chr-x")
    else:
        arguments.append("--exclude-chr-x")
    assert invoke() == 0
    assert invoke() == 0
    assert calls == [initially_excluded, not initially_excluded]
    assert (output / "scientific_output.txt").read_text() == str(not initially_excluded)
    (output / "scientific_output.txt").write_text("externally modified")
    if initially_excluded:
        arguments.append("--exclude-chr-x")
    else:
        arguments.remove("--exclude-chr-x")
    with pytest.raises(RuntimeError, match="Refusing automatic restart.*was modified"):
        invoke()
    assert calls == [initially_excluded, not initially_excluded]
    assert (output / "scientific_output.txt").read_text() == "externally modified"
