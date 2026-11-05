"""Replacement reports distinguish candidate selection from row-bound approval."""

from copy import deepcopy
import json
import os
import subprocess
import sys
import textwrap

import polars as pl
import pytest

from postgwas.modules.harmonisation import allele_frequency as af
from postgwas.modules.harmonisation.html_report import render_dataset_report
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.study_reporting import render_study_decisions
from test_harmonisation_replacement_af_integration import fixture, validate


def normal(value):
    return " ".join(value.split())


def reports(decisions):
    return (normal(render_study_decisions(decisions)),
            render_dataset_report({}, {"dataset": {"study_decisions": decisions}}))


@pytest.mark.parametrize("same", [False, True])
def test_passed_reports_show_final_approval_and_source_reuse_without_mutation(tmp_path, same):
    decisions = validate(fixture(tmp_path, same=same))
    before = deepcopy(decisions)
    for report in reports(decisions):
        assert "Replacement validation" in report
        assert "Passed; replacement EAF accepted for assignment" in report
        assert "assignment pending" not in report
        assert "Replacement correlation" in report
        assert "1; required &gt;= 0.95" in report or "1; required >= 0.95" in report
        assert "Replacement major frequency" in report
        assert "6 / 12 (50.00%); required &gt; 20.00%" in report or "6 / 12 (50.00%); required > 20.00%" in report
        assert "Original study AF preserved" in report
        assert "population-reference proxy" in report
        assert "Reused default-reference evidence" in report
        assert ("Yes" if same else "No") in report
        assert "full identities retained in structured evidence" in report.lower()
        assert str(tmp_path) in report
    identities = decisions["eaf_replacement_validation"]["source"]["identities"]
    assert set(identities) == {"1"}
    assert identities["1"]["path"] == str(tmp_path / ("default.tsv" if same else "replacement.tsv"))
    assert decisions == before


def test_provisional_selection_never_claims_acceptance_or_assignment_permission(tmp_path):
    decisions = fixture(tmp_path)[2]
    for report in reports(decisions):
        assert "Provisionally selected; replacement validation pending" in report
        assert "accepted for assignment" not in report
        assert "Assign selected reference EAF during chromosome processing" not in report


@pytest.mark.parametrize("override,reason", [
    ({"eaf.maf_reference_min_overlap": 13}, "insufficient overlap"),
    ({"eaf.maf_reference_min_overlap_fraction": 1.0}, "insufficient overlap fraction"),
    ({"eaf.study_reference_min_correlation": .99}, "correlation below minimum"),
    ({"eaf.maf_reference_max_mean_absolute_difference": .001}, "mean absolute difference above maximum"),
    ({"eaf.study_reference_min_major_fraction": 1.0}, "major fraction not above minimum"),
])
def test_numeric_failures_report_fixed_denominator_thresholds_and_all_reasons(tmp_path, override, reason):
    decisions = fixture(tmp_path)[2]
    frame = pl.DataFrame({"A": [.1, .2, .3, .4] * 3,
                          "B": [.8, .8, .7, .6] * 3})
    policies = default_policies().with_overrides({"eaf.maf_reference_min_overlap": 2, **override})
    decisions["eaf_replacement_validation"] = af.validate_replacement_frequency_pairs(
        frame, "A", "B", policies, study_variant_count=24,
        source=decisions["eaf_replacement"], excluded_counts={"unmatched_original_rows": 12})
    for report in reports(decisions):
        assert "Inconclusive; dataset stopped before chromosome analysis" in report
        assert reason in report
        assert "12 / 24 (50.00%)" in report
        assert "Replacement correlation" in report
        assert "Replacement mean absolute difference" in report
        assert "Replacement major frequency" in report
        assert "unmatched original rows" in report
        assert "accepted for assignment" not in report
        if not decisions["eaf_replacement_validation"]["shared_evidence_gate_passed"]:
            assert "not available; required" in report
        for failed in decisions["eaf_replacement_validation"]["failed_requirements"]:
            assert failed.replace("_", " ") in report


def test_late_identity_failure_overrides_earlier_numeric_pass(tmp_path):
    decisions = validate(fixture(tmp_path))
    evidence = decisions["eaf_replacement_validation"]
    evidence.update(status="invalid", failure_reason="replacement_inputs_changed")
    evidence.pop("contract")
    for report in reports(decisions):
        assert "Invalid; dataset stopped before chromosome analysis" in report
        assert "replacement inputs changed" in report
        assert "accepted for assignment" not in report
        assert "Replacement correlation" in report


def test_early_io_failure_keeps_replacement_unvalidated(tmp_path):
    decisions = fixture(tmp_path)[2]
    decisions["report_failure"] = "Permission denied while opening replacement source"
    for report in reports(decisions):
        assert "Not validated; dataset stopped before replacement approval" in report
        assert "accepted for assignment" not in report
    decisions.pop("report_failure")
    html = render_dataset_report({}, {"dataset": {"study_decisions": decisions},
        "failure": {"message": "Permission denied while opening replacement source"}})
    assert "Not validated; dataset stopped before replacement approval" in html
    assert "replacement validation pending" not in html


def test_numeric_status_without_adoption_contract_is_not_final_approval(tmp_path):
    decisions = validate(fixture(tmp_path))
    decisions["eaf_replacement_validation"].pop("contract")
    for report in reports(decisions):
        assert "Numeric gates passed; replacement approval not recorded" in report
        assert "accepted for assignment" not in report


@pytest.mark.parametrize("visible", [False, True])
def test_final_replacement_report_survives_visible_and_hidden_canonical_recording(tmp_path, visible):
    from postgwas.core.pipeline_logging import _PREFIX_WIDTH
    decisions = validate(fixture(tmp_path / "inputs"))
    (tmp_path / "decisions.json").write_text(json.dumps(decisions))
    script = textwrap.dedent("""
        import json
        from pathlib import Path
        import sys
        from postgwas.core.pipeline_logging import PipelineLogger
        from postgwas.core.screen_logging import ScreenSettings, record_screen
        from postgwas.modules.harmonisation.study_reporting import render_replacement_frequency_validation
        from postgwas.modules.harmonisation.service import _announce
        output = Path(sys.argv[1])
        decisions = json.loads((output / 'decisions.json').read_text())
        with PipelineLogger('study', 'dataset', str(output / 'logs'), screen_level='ERROR') as logger:
            with record_screen(ScreenSettings(sys.argv[2] == 'True', output / 'screen.log')):
                _announce(logger, render_replacement_frequency_validation(decisions))
    """)
    completed = subprocess.run([sys.executable, "-B", "-c", script, str(tmp_path), str(visible)],
        env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1", POLARS_MAX_THREADS="4"),
        capture_output=True, text=True, timeout=60, check=False)
    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    assert ("accepted for assignment" in completed.stdout) is visible
    for path in (tmp_path / "screen.log", tmp_path / "logs" / "study_dataset.log"):
        raw = path.read_text()
        if path.name == "study_dataset.log":
            raw = "\n".join(line[_PREFIX_WIDTH:] for line in raw.splitlines())
        report = normal(raw)
        assert "Passed; replacement EAF accepted for assignment" in report
        assert "\x1b" not in report


@pytest.mark.parametrize("mode", ["distinct", "inconclusive"])
def test_real_preparation_report_contains_only_terminal_replacement_outcome(tmp_path, monkeypatch, capsys, mode):
    from test_harmonisation_eaf_replacement_resources import test_real_dataset_preparation_selects_then_stages_replacement
    test_real_dataset_preparation_selects_then_stages_replacement(tmp_path, monkeypatch, mode)
    report = normal(capsys.readouterr().out)
    assert "Provisionally selected; replacement validation pending" not in report
    assert report.count("Study-wide decisions and evidence") == 1
    assert report.count("Selected build") == 1
    if mode == "distinct":
        assert "accepted for assignment" in report
        assert "INFO" in report.split("Study-wide decisions and evidence", 1)[1]
        log = (tmp_path / "logs" / "study_dataset.log").read_text()
        assert log.index("Split by chromosome") < log.index("Dataset-wide INFO score type")
        assert log.index("Dataset-wide INFO score type") < log.index("Study-wide decisions and evidence")
    else:
        assert "Inconclusive; dataset stopped before chromosome analysis" in report
        assert "accepted for assignment" not in report
