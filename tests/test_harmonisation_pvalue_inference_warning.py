"""Dataset-level p-value inference notices preserve scientific decisions."""

import math
from pathlib import Path
import subprocess
import sys
import textwrap

import polars as pl
import pytest

from postgwas.core.pipeline_logging import PipelineLogger, _PREFIX_WIDTH
from postgwas.modules.harmonisation.p_values import harmonise_p_values
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.study_properties import (
    StudyPropertyError,
    resolve_study_properties,
)


def _study_values(kind, size=1000):
    values = [(index + 0.5) / size for index in range(size)]
    return values if kind == "raw" else [-math.log10(value) for value in values]


@pytest.mark.parametrize("kind", ["raw", "neglog10"])
def test_automatic_pvalue_notice_is_structured_and_logged_once(
    kind, tmp_path, capsys,
):
    policies = default_policies()
    frame = pl.DataFrame({"ASSOCIATION_P": _study_values(kind)})
    original = frame.clone()
    mapping = {"pval_col": "ASSOCIATION_P"}
    logger = PipelineLogger(
        "study", "dataset", str(tmp_path), policies=policies,
        screen_level="ERROR",
    )
    try:
        decisions = resolve_study_properties(
            frame, mapping, policies=policies, logger=logger,
        )
    finally:
        logger.close()

    assert frame.equals(original)
    assert mapping == {"pval_col": "ASSOCIATION_P"}
    assert decisions["pvalue_type"] == kind
    assert decisions["pvalue_type_detected"] == kind
    assert decisions["pvalue_type_source"] == "detector"
    warning = decisions["pvalue_type_inference_warning"]
    assert "AUTOMATIC P-VALUE-TYPE INFERENCE" in warning
    assert "ASSOCIATION_P" in warning
    assert "p_value_type" in warning
    assert "pvalue.type" in warning

    screen = " ".join(capsys.readouterr().out.split())
    assert screen.count("P-value type automatically inferred") == 1
    for label in (
        "Input column", "Inferred type", "Evidence used", "Usable values",
        "Observed range", "Above threshold", "-log10 rule", "Median value",
        "Planned action", "Why this warning", "Recommended", "Run behaviour",
    ):
        assert label in screen
    assert "Continuing with the inferred type" in screen
    assert "does not determine whether tests are one- or two-sided" in screen
    assert "p_value_type = " + kind in screen
    if kind == "raw":
        assert "raw probability scale" in screen
        assert "not used to select raw scale" in screen
    else:
        assert "P = 10^(-input value)" in screen
        assert "inclusive" in screen
    log_text = Path(logger.log_path).read_text(encoding="utf-8")
    assert log_text.count("AUTOMATIC P-VALUE-TYPE INFERENCE") == 1
    # Canonical records wrap with a scope/time prefix on every continuation.
    log_content = " ".join(
        line[_PREFIX_WIDTH:].strip() for line in log_text.splitlines()
    )
    assert warning in log_content


def test_notice_uses_finite_evidence_and_resolved_thresholds(capsys):
    policies = default_policies().with_overrides({
        "pvalue.mlogp_detect_threshold": 2.0,
        "pvalue.mlogp_detect_proportion": 0.2,
        "pvalue.mlogp_detect_min_count": 2,
        "pvalue.mlogp_expected_median": 0.4,
        "pvalue.mlogp_median_tolerance": 0.1,
    })
    # Exactly 2.0 is not above the threshold; invalid tokens do not dilute 2/10.
    values = ["3", "4", "2"] + ["0.4"] * 7 + [
        None, "not-a-number", "NaN", "inf", "-inf",
    ]
    decisions = resolve_study_properties(
        pl.DataFrame({"P": values}), {"pval_col": "P"}, policies=policies,
    )

    assert decisions["pvalue_type"] == "neglog10"
    evidence = decisions["pvalue_evidence"]["detector_stats"]
    assert evidence["n_usable"] == 10
    assert evidence["n_above_threshold"] == 2
    assert evidence["above_threshold_fraction"] == pytest.approx(0.2)
    screen = " ".join(capsys.readouterr().out.split())
    assert "2 / 10 (20.0000%)" in screen
    assert "Excluded values" in screen
    assert "5 missing, non-numeric or non-finite values" in screen
    assert "20%" in screen
    assert "0.3 to 0.5 inclusive" in screen
    assert "1.05" not in screen


@pytest.mark.parametrize("median", [0.25, 0.75])
def test_notice_keeps_inclusive_median_and_count_fraction_boundaries(
    median, capsys,
):
    policies = default_policies().with_overrides({
        "pvalue.mlogp_detect_proportion": 0.1,
        "pvalue.mlogp_detect_min_count": 2,
        "pvalue.mlogp_expected_median": 0.5,
        "pvalue.mlogp_median_tolerance": 0.25,
    })
    decisions = resolve_study_properties(
        pl.DataFrame({"P": [2.0, 3.0] + [median] * 18}),
        {"pval_col": "P"}, policies=policies,
    )

    assert decisions["pvalue_type"] == "neglog10"
    evidence = decisions["pvalue_evidence"]["detector_stats"]
    assert evidence["count_supports_neglog10"] is True
    assert evidence["fraction_supports_neglog10"] is True
    screen = " ".join(capsys.readouterr().out.split())
    assert "2 / 20 (10.0000%)" in screen
    assert "0.25 to 0.75 inclusive" in screen


def test_raw_notice_does_not_claim_every_value_is_valid(capsys):
    decisions = resolve_study_properties(
        pl.DataFrame({"P": [0.5] * 10_000 + [999.0]}),
        {"pval_col": "P"}, policies=default_policies(),
    )

    assert decisions["pvalue_type"] == "raw"
    assert decisions["pvalue_evidence"]["detector_stats"]["max"] == 999.0
    screen = " ".join(capsys.readouterr().out.split())
    assert "1 / 10,001 (0.0100%)" in screen
    assert "999" in screen
    assert "chromosome QC still applies" in screen
    assert "not used to select raw scale" in screen
    assert "all values are valid" not in screen


@pytest.mark.parametrize(
    ("values", "message"),
    [
        ([], "no rows"),
        ([None, None], "no usable numeric value"),
        (["not-a-number", None], "no usable numeric value"),
        (["NaN", "inf", "-inf"], "No numeric P-values"),
        ([-0.2, 0.5], "negative usable value"),
        ([0.5] * 100 + [999.0], "single sentinel"),
        ([3.0] * 100, "does not agree"),
    ],
)
def test_failed_detection_never_announces_successful_inference(
    values, message, capsys,
):
    with pytest.raises(StudyPropertyError, match=message):
        resolve_study_properties(
            pl.DataFrame({"P": values}), {"pval_col": "P"},
            policies=default_policies(),
        )

    screen = capsys.readouterr().out
    assert "P-value type automatically inferred" not in screen
    assert "AUTOMATIC P-VALUE-TYPE INFERENCE" not in screen
    assert "Continuing with the inferred type" not in screen


@pytest.mark.parametrize("kind", ["raw", "neglog10"])
@pytest.mark.parametrize("source", ["policy", "sample_sheet"])
def test_explicit_pvalue_type_is_not_reported_as_applied_inference(
    kind, source, capsys,
):
    mapping = {"pval_col": "P"}
    policies = default_policies()
    if source == "sample_sheet":
        mapping["declared_pvalue_type"] = kind
    else:
        policies = policies.with_overrides({"pvalue.type": kind})
    decisions = resolve_study_properties(
        pl.DataFrame({"P": _study_values(kind)}), mapping, policies=policies,
    )

    assert decisions["pvalue_type"] == kind
    assert decisions["pvalue_type_source"] == source
    assert "pvalue_type_inference_warning" not in decisions
    assert "P-value type automatically inferred" not in capsys.readouterr().out


def test_declaration_mismatch_remains_separate_from_applied_inference(capsys):
    decisions = resolve_study_properties(
        pl.DataFrame({"P": _study_values("neglog10")}),
        {"pval_col": "P", "declared_pvalue_type": "raw"},
        policies=default_policies(),
    )

    assert decisions["pvalue_type"] == "raw"
    assert decisions["pvalue_type_detected"] == "neglog10"
    assert decisions["pvalue_type_matches_declaration"] is False
    assert "pvalue_type_inference_warning" not in decisions
    screen = capsys.readouterr().out
    assert "DECLARATION MISMATCH" in screen
    assert "P-value type automatically inferred" not in screen


def test_inconclusive_cross_check_keeps_explicit_declaration(capsys):
    decisions = resolve_study_properties(
        pl.DataFrame({"P": [3.0] * 100}),
        {"pval_col": "P", "declared_pvalue_type": "neglog10"},
        policies=default_policies(),
    )

    assert decisions["pvalue_type"] == "neglog10"
    assert decisions["pvalue_type_source"] == "sample_sheet"
    assert decisions["pvalue_type_detected"] is None
    assert "pvalue_type_inference_warning" not in decisions
    assert "P-value type automatically inferred" not in capsys.readouterr().out


def test_tolerated_negative_warning_is_preserved_alongside_inference(capsys):
    decisions = resolve_study_properties(
        pl.DataFrame({"P": [-0.2] + [0.5] * 999}),
        {"pval_col": "P"}, policies=default_policies(),
    )

    assert decisions["pvalue_type"] == "raw"
    assert "NEGATIVE P-VALUES" in decisions["pvalue_negative_warning"]
    assert "pvalue_type_inference_warning" in decisions
    screen = " ".join(capsys.readouterr().out.split())
    assert "NEGATIVE P-VALUES" in screen
    assert "Negative values" in screen
    assert "1 / 1,000" in screen
    assert screen.count("P-value type automatically inferred") == 1


@pytest.mark.parametrize("action", ["null", "fail"])
def test_tolerated_negative_notice_respects_configured_later_range_action(
    action, capsys,
):
    policies = default_policies().with_overrides({"pvalue.out_of_range": action})
    decisions = resolve_study_properties(
        pl.DataFrame({"P": [-0.2] + [0.5] * 999}),
        {"pval_col": "P"}, policies=policies,
    )

    assert decisions["pvalue_type"] == "raw"
    warning = decisions["pvalue_negative_warning"]
    assert "pvalue.out_of_range" in warning
    assert action in warning
    assert "will be rejected" not in warning
    assert "will be rejected" not in capsys.readouterr().out


def test_no_pvalue_column_does_not_emit_pvalue_inference(capsys):
    decisions = resolve_study_properties(
        pl.DataFrame({"OTHER": [1.0]}), {}, policies=default_policies(),
    )

    assert decisions["pvalue_type"] is None
    assert decisions["pvalue_type_source"] == "no_pvalue_column"
    assert "pvalue_type_inference_warning" not in decisions
    assert "P-value type automatically inferred" not in capsys.readouterr().out


def test_warning_reuses_the_existing_single_study_scan(monkeypatch, capsys):
    frame = pl.DataFrame({"P": [0.1, 0.5, 1.0]})
    original_select = pl.DataFrame.select
    scans = []

    def track_select(self, *args, **kwargs):
        if self is frame:
            scans.append(True)
        return original_select(self, *args, **kwargs)

    monkeypatch.setattr(pl.DataFrame, "select", track_select)
    decisions = resolve_study_properties(
        frame, {"pval_col": "P"}, policies=default_policies(),
    )

    assert decisions["pvalue_type"] == "raw"
    assert len(scans) == 1
    assert "P-value type automatically inferred" in capsys.readouterr().out


@pytest.mark.parametrize("verify", [False, True])
def test_chromosomes_apply_one_study_decision_without_repeating_notice(
    verify, capsys,
):
    policies = default_policies().with_overrides({
        "pvalue.verify_per_chromosome": verify,
    })
    frame = pl.DataFrame({"P": _study_values("neglog10")})
    decisions = resolve_study_properties(
        frame, {"pval_col": "P"}, policies=policies,
    )
    assert capsys.readouterr().out.count("P-value type automatically inferred") == 1

    for chromosome, subset in enumerate(frame.iter_slices(500), start=1):
        result, qc, mapping = harmonise_p_values(
            str(chromosome), subset, {"pval_col": "P"},
            policies=policies, decision=decisions["pvalue_type"],
        )
        expected = [10.0 ** -value for value in subset["P"].to_list()]
        assert result[mapping["pval_col"]].to_list() == pytest.approx(expected)
        assert qc["detected_scale"] == "neglog10"
        assert "pvalue_type_inference_warning" not in qc
    screen = capsys.readouterr().out
    assert "P-value type automatically inferred" not in screen


@pytest.mark.parametrize("show_screen", [False, True])
def test_notice_reaches_dataset_and_shared_transcripts_when_screen_is_hidden(
    show_screen, tmp_path,
):
    script = textwrap.dedent("""
        from pathlib import Path
        import sys
        import polars as pl
        from postgwas.core.pipeline_logging import PipelineLogger
        from postgwas.core.screen_logging import ScreenSettings, record_screen
        from postgwas.modules.harmonisation.cli import _DatasetScreenRouter
        from postgwas.modules.harmonisation.policies import default_policies
        from postgwas.modules.harmonisation.study_properties import resolve_study_properties

        output = Path(sys.argv[1])
        show_screen = sys.argv[2] == 'True'
        policies = default_policies()
        logger = PipelineLogger(
            'study', 'dataset', str(output), policies=policies, screen_level='ERROR',
        )
        try:
            with record_screen(ScreenSettings(show_screen, output / 'screen.log')):
                with _DatasetScreenRouter(
                    {'study': output / 'study_screen_report.txt'}, display=show_screen,
                ) as router:
                    router.select('study')
                    resolve_study_properties(
                        pl.DataFrame({'P': [0.1, 0.5, 1.0]}),
                        {'pval_col': 'P'}, policies=policies, logger=logger,
                    )
        finally:
            logger.close()
    """)
    completed = subprocess.run(
        [sys.executable, "-c", script, str(tmp_path), str(show_screen)],
        capture_output=True, text=True, timeout=30, check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    if show_screen:
        assert completed.stdout.count("P-value type automatically inferred") == 1
    else:
        assert completed.stdout == ""
    for path in (tmp_path / "screen.log", tmp_path / "study_screen_report.txt"):
        text = " ".join(path.read_text(encoding="utf-8").split())
        assert text.count("P-value type automatically inferred") == 1
        assert "Input column" in text
        assert "Evidence used" in text
        assert "Recommended" in text
        assert "Continuing with the inferred type" in text
        assert "\x1b" not in text
    log_text = (tmp_path / "study_dataset.log").read_text(encoding="utf-8")
    assert log_text.count("AUTOMATIC P-VALUE-TYPE INFERENCE") == 1
