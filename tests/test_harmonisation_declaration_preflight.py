"""Invalid scientific declarations stop at the real public CLI input boundary."""

from contextlib import ExitStack
import json
import os
from pathlib import Path
import subprocess
from types import SimpleNamespace
import sys
from unittest.mock import patch

import pytest

from postgwas.__main__ import main as public_main
from postgwas.config import write_resolved_configuration
from postgwas.core import direct_execution
from postgwas.core.paths import configured_output_path
from postgwas.modules.harmonisation import cli, service, vcf_processing

from test_harmonisation_batch_preflight import config, report, sheet, write_sheet


INVALID_DECLARATIONS = (
    ("trait_type", "Case_Contol", ("auto", "quantitative", "case_control")),
    ("effect_type", "betaa", ("auto", "beta", "odds_ratio")),
    ("p_value_type", "negln", ("auto", "raw", "neglog10")),
    ("delimiter", "tabb", ("auto", "tab", "comma", "semicolon", "space", "whitespace")),
)


def _dispatch_public_cli(arguments, study_paths):
    study_paths = {Path(path).resolve() for path in study_paths}
    events = {"preparation": [], "content_read": [], "workers": [], "external_tools": []}
    original_open = Path.open

    def guarded_open(path, *args, **kwargs):
        if path.resolve() in study_paths:
            events["content_read"].append(str(path))
            raise AssertionError("Input validation must not open invalid study content")
        return original_open(path, *args, **kwargs)

    def prepare(**kwargs):
        dataset = kwargs["sample_column_dict"]["gwas_outputname"]
        events["preparation"].append(dataset)
        raise AssertionError("CSV errors must stop before scientific preparation: " + dataset)

    def read_content(*args, **kwargs):
        events["content_read"].append("summary statistics reader")
        raise AssertionError("Invalid declarations must not reach the study reader")

    def dispatch_valid_rows(*args, **kwargs):
        # Orchestration now receives invalid-row evidence for reporting, but
        # scientific preparation receives only independently valid models.
        rejected = kwargs.get("preflight_failures") or {}
        assert all(row.dataset_id not in rejected for row in args[2])
        raise cli.HarmonisationPreflightError("Boundary checked; fixture stops before valid-row execution")

    def external_command(*args, **kwargs):
        events["external_tools"].append("native command")
        raise AssertionError("Invalid declarations must not execute external tools")

    # Shared version/resource checks are independent, expensive preflight
    # boundaries. Keep the actual parser, CSV model, issue collection,
    # public exit handling, canonical reporting and screen recorder in use.
    with ExitStack() as guards:
        for owner, name, replacement in (
            (Path, "open", guarded_open),
            (direct_execution, "software_identity", lambda *args, **kwargs: {}),
            (cli, "require_binaries", lambda *args, **kwargs: {}),
            (cli, "preflight_configured_reference_scope", lambda *args, **kwargs: None),
            (cli, "run_harmonisation_pipeline", prepare),
            (service, "read_summary_statistics", read_content),
            (cli, "_run_validated_rows", dispatch_valid_rows),
            (vcf_processing, "run_checked_command", external_command),
            (sys, "argv", arguments),
        ):
            guards.enter_context(patch.object(owner, name, replacement))
        exit_code = public_main()
    return exit_code, events


def _invoke_public_cli(tmp_path, capfd, records, *, hidden=False):
    cfg = config(tmp_path)
    cfg.logging.terminal_style.color = "never"
    configuration = tmp_path / "run.yaml"
    write_resolved_configuration(cfg, configuration)
    manifest = write_sheet(tmp_path / "studies.csv", records)
    study_paths = {Path(record["input_file"]).resolve() for record in records}
    original_sources = {path: path.read_bytes() for path in study_paths}
    arguments = [
        "postgwas", "harmonisation", "--sample-sheet", str(manifest),
        "--run-config", str(configuration),
        "--resource-directory", str(cfg.resources.root),
        "--output-directory", str(cfg.run.output_directory),
    ]
    capfd.readouterr()
    if hidden:
        arguments.append("--hide-screen")
        request = tmp_path / "cli_request.json"
        observed = tmp_path / "cli_observed.json"
        request.write_text(json.dumps({
            "arguments": arguments, "study_paths": [str(path) for path in study_paths],
        }))
        # A real child process keeps stdout/stderr attached to descriptors 1/2;
        # pytest's substituted streams bypass the process-level screen recorder.
        driver = (
            "import json, sys\nfrom pathlib import Path\n"
            "from test_harmonisation_declaration_preflight import _dispatch_public_cli\n"
            "request = json.loads(Path(sys.argv[1]).read_text())\n"
            "code, events = _dispatch_public_cli(**request)\n"
            "Path(sys.argv[2]).write_text(json.dumps(events))\n"
            "raise SystemExit(code)\n"
        )
        repository = Path(__file__).resolve().parents[1]
        environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
        environment["PYTHONPATH"] = os.pathsep.join([
            str(repository / "tests"), str(repository / "src"),
            environment.get("PYTHONPATH", ""),
        ])
        result = subprocess.run(
            [sys.executable, "-B", "-c", driver, str(request), str(observed)],
            capture_output=True, text=True, check=False, env=environment,
        )
        exit_code = result.returncode
        events = json.loads(observed.read_text())
        captured = SimpleNamespace(out=result.stdout, err=result.stderr)
    else:
        exit_code, events = _dispatch_public_cli(arguments, study_paths)
        captured = capfd.readouterr()
    assert exit_code != 0
    assert events["preparation"] == []
    assert events["content_read"] == events["workers"] == events["external_tools"] == []
    assert all(path.read_bytes() == original for path, original in original_sources.items())
    runtime = cfg.modules.harmonisation.runtime
    metadata = configured_output_path(cfg.run.output_directory, runtime.top_metadata_directory)
    canonical_log = configured_output_path(metadata, cfg.logging.filename).read_text(encoding="utf-8")
    screen_log = configured_output_path(cfg.run.output_directory, cfg.logging.screen_log_file).read_text(encoding="utf-8")
    observations = report(cfg)
    assert "preflight" in (canonical_log + screen_log).lower()
    assert not list(cfg.run.output_directory.rglob("*.parquet"))
    assert not list(cfg.run.output_directory.rglob("*.vcf*"))
    if hidden:
        assert captured.out == captured.err == ""
    return canonical_log, observations, " ".join((captured.out + captured.err).split())


@pytest.mark.parametrize("field,original,choices", INVALID_DECLARATIONS)
def test_invalid_declaration_is_reported_before_scientific_preparation(
    tmp_path, capfd, field, original, choices,
):
    records = sheet(tmp_path, 1)
    records[0][field] = original

    canonical_log, observations, terminal = _invoke_public_cli(tmp_path, capfd, records)

    failures = [item for item in observations if item["stage"] == "sample sheet" and item["status"] == "failed"]
    assert failures
    message = "\n".join(item["problem"] for item in failures)
    for expected in ("row 2", "study0", field, original, *choices):
        assert expected in message
        assert expected in canonical_log
        assert expected in terminal
    assert "normalized to 'auto'" not in message


def test_all_invalid_rows_and_declarations_are_preserved_in_hidden_preflight_reports(
    tmp_path, capfd,
):
    records = sheet(tmp_path, 2)
    records[0].update(trait_type="Case_Contol", effect_type="betaa")
    records[1].update(p_value_type="negln", delimiter="tabb")

    canonical_log, observations, _ = _invoke_public_cli(
        tmp_path, capfd, records, hidden=True,
    )

    failures = [item for item in observations if item["stage"] == "sample sheet" and item["status"] == "failed"]
    assert {item["row"] for item in failures} == {"2", "3"}
    assert {item["dataset_id"] for item in failures} == {"study0", "study1"}
    message = "\n".join(item["problem"] for item in failures)
    for field, original, choices in INVALID_DECLARATIONS:
        for expected in (field, original, *choices):
            assert expected in message
            assert expected in canonical_log


def test_invalid_row_does_not_block_reporting_of_independently_valid_rows(
    tmp_path, capfd,
):
    records = sheet(tmp_path, 2)
    records[0]["effect_type"] = "betaa"

    canonical_log, observations, _ = _invoke_public_cli(
        tmp_path, capfd, records,
    )

    failures = [item for item in observations if item["stage"] == "sample sheet" and item["status"] == "failed"]
    assert failures and {item["dataset_id"] for item in failures} == {"study0"}
    assert "betaa" in canonical_log
    assert any(item["dataset_id"] == "study1" and item["status"] == "passed" for item in observations)
