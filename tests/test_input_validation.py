"""Transient validation reuse must preserve exact contracts and truthful reports."""

import json
import os

import pytest

from postgwas.core.input_validation import (
    InputValidationSession,
    capture_preflight_file_identities,
    current_validation_session,
    record_file_validation,
    require_unchanged_preflight_files,
    validate_once,
    validation_scope,
)
from postgwas.core.paths import require_nonempty_file


def _file(tmp_path, name="reference.txt"):
    source = tmp_path / name
    source.write_text("first\n", encoding="utf-8")
    return source


def test_direct_validation_does_not_cache_or_record(tmp_path):
    source = _file(tmp_path)
    calls = []
    with validation_scope("direct"):
        assert current_validation_session() is None
        record_file_validation(source, "Reference", metrics={"ignored": object()})
        assert require_nonempty_file(source, "Reference") == source
        for _ in range(2):
            validate_once([source], "check", lambda: calls.append(True))
    assert calls == [True, True]


def test_scope_records_truthful_availability_and_merges_consumers(tmp_path):
    source = _file(tmp_path)
    with InputValidationSession() as session:
        for module in ("first", "second", "first"):
            with session.scope(module):
                require_nonempty_file(source, "Reference")
        assert len(session.records) == 1
        record = session.records[0]
        assert record.path == str(source)
        assert record.consumers == ("first", "second")
        assert record.status == "passed"
        assert record.checks == ("regular file", "non-empty file")
        assert record.metrics == {}
    assert current_validation_session() is None


@pytest.mark.parametrize("kind", ["missing", "empty", "directory", "unspecified"])
def test_availability_failure_keeps_existing_error_and_reports(tmp_path, kind):
    source = tmp_path / kind
    if kind == "empty":
        source.touch()
    elif kind == "directory":
        source.mkdir()
    elif kind == "unspecified":
        source = None
    with InputValidationSession() as session, session.scope("consumer"):
        with pytest.raises(ValueError, match="Reference"):
            require_nonempty_file(source, "Reference")
    assert len(session.records) == 1
    assert session.records[0].status == "failed"
    assert session.records[0].consumers == ("consumer",)


def test_shared_identity_capture_records_exact_file_failure(tmp_path):
    source = _file(tmp_path)
    absent = tmp_path / "absent"
    with InputValidationSession() as session, session.scope("module"):
        with pytest.raises(RuntimeError, match="unavailable"):
            capture_preflight_file_identities([source, absent], label="LD bundle")
    assert [(record.path, record.status) for record in session.records] == [
        (str(source), "passed"), (str(absent), "failed"),
    ]
    assert all(record.consumers == ("module",) for record in session.records)


def test_exact_contract_cache_replays_validation_evidence_for_each_consumer(tmp_path):
    source = _file(tmp_path)
    calls = []

    def inspect():
        calls.append(True)
        record_file_validation(
            source, "Reference", checks=("all rows inspected",),
            metrics={"variants": 2, "identifier_type": "rsid"},
        )
        return "rsid"

    with InputValidationSession() as session:
        for module in ("gcta_cojo", "magma", "formatter"):
            with validation_scope(module):
                assert validate_once(
                    [source], {"validator": "BIM", "pattern": "rs"}, inspect,
                ) == "rsid"
    assert calls == [True]
    assert all(
        record.consumers == ("gcta_cojo", "magma", "formatter")
        for record in session.records
    )
    assert any(record.metrics.get("variants") == 2 for record in session.records)


def test_different_reference_or_scientific_contract_requires_own_validation(tmp_path):
    first = _file(tmp_path, "first")
    second = _file(tmp_path, "second")
    calls = []
    with InputValidationSession():
        for paths, contract in (
            ([first], {"validator": "table", "columns": ["beta"]}),
            ([second], {"validator": "table", "columns": ["beta"]}),
            ([first], {"validator": "table", "columns": ["se"]}),
        ):
            validate_once(paths, contract, lambda: calls.append(True))
    assert calls == [True, True, True]


def test_mapping_order_does_not_change_contract_identity(tmp_path):
    source = _file(tmp_path)
    calls = []
    with InputValidationSession():
        validate_once([source], {"validator": "x", "policy": "fail"},
                      lambda: calls.append(True))
        validate_once([source], {"policy": "fail", "validator": "x"},
                      lambda: calls.append(True))
    assert calls == [True]


def test_changed_file_cannot_reuse_or_refresh_stale_session_evidence(tmp_path):
    source = _file(tmp_path)
    calls = []
    with InputValidationSession() as session:
        validate_once([source], "check", lambda: calls.append(True))
        source.write_text("changed-size\n", encoding="utf-8")
        with pytest.raises(RuntimeError, match="changed after pipeline preflight"):
            validate_once([source], "check", lambda: calls.append(True))
    assert calls == [True]
    assert any(record.status == "failed" for record in session.records)


def test_file_changed_during_validation_is_not_cached(tmp_path):
    source = _file(tmp_path)

    def change_during_check():
        source.write_text("changed-size\n", encoding="utf-8")
        return "unsafe result"

    with InputValidationSession() as session:
        with pytest.raises(RuntimeError, match="changed after pipeline preflight"):
            validate_once([source], "check", change_during_check)
    assert any(record.status == "failed" for record in session.records)


def test_same_size_replacement_with_restored_mtime_does_not_reuse_evidence(tmp_path):
    source = _file(tmp_path)
    original = source.stat()
    replacement = tmp_path / "replacement"
    replacement.write_text("other\n", encoding="utf-8")
    os.utime(replacement, ns=(original.st_atime_ns, original.st_mtime_ns))
    with InputValidationSession():
        validate_once([source], "check", lambda: "first")
        replacement.replace(source)
        assert source.stat().st_size == original.st_size
        assert source.stat().st_mtime_ns == original.st_mtime_ns
        with pytest.raises(RuntimeError, match="changed after pipeline preflight"):
            validate_once([source], "check", lambda: "would be stale")


def test_failed_operation_preserves_exception_and_does_not_blame_each_file(tmp_path):
    first = _file(tmp_path, "study")
    second = _file(tmp_path, "reference")
    calls = []

    def incompatible():
        calls.append(True)
        raise ValueError("Study/reference genome builds differ")

    with InputValidationSession() as session, session.scope("analysis"):
        for _ in range(2):
            with pytest.raises(ValueError, match="genome builds differ"):
                validate_once([first, second], "build compatibility", incompatible)
    assert calls == [True, True]
    failures = [record for record in session.records if record.status == "failed"]
    assert len(failures) == 1
    assert failures[0].path is None
    assert failures[0].metrics["paths"] == [str(first), str(second)]


def test_validation_missing_file_reports_without_running_operation(tmp_path):
    source = tmp_path / "absent"
    calls = []
    with InputValidationSession() as session:
        with pytest.raises(RuntimeError, match="unavailable"):
            validate_once([source], "check", lambda: calls.append(True))
    assert calls == []
    assert session.records[0].path == str(source)
    assert session.records[0].status == "failed"


def test_nested_cached_validators_replay_all_observed_records(tmp_path):
    source = _file(tmp_path)
    calls = []

    def inner():
        calls.append("inner")
        record_file_validation(source, "Inner check", checks=("schema",))
        return 2

    def outer():
        calls.append("outer")
        count = validate_once([source], "inner", inner)
        record_file_validation(source, "Outer check", metrics={"rows": count})
        return count

    with InputValidationSession() as session:
        with session.scope("first"):
            assert validate_once([source], "outer", outer) == 2
        with session.scope("second"):
            assert validate_once([source], "outer", outer) == 2
    assert calls == ["outer", "inner"]
    assert all(record.consumers == ("first", "second") for record in session.records)


@pytest.mark.parametrize("warm_inner_cache", [False, True])
def test_parent_cache_refuses_changed_nested_resource(tmp_path, warm_inner_cache):
    parent = _file(tmp_path, "parent")
    child = _file(tmp_path, "child")
    calls = []

    def inner():
        calls.append("inner")
        return child.read_text(encoding="utf-8")

    def outer():
        calls.append("outer")
        return validate_once([child], "inner", inner)

    with InputValidationSession() as session:
        if warm_inner_cache:
            validate_once([child], "inner", inner)
        assert validate_once([parent], "outer", outer) == "first\n"
        child.write_text("changed child\n", encoding="utf-8")
        with pytest.raises(RuntimeError, match="changed after pipeline preflight"):
            validate_once([parent], "outer", outer)
    assert calls.count("inner") == calls.count("outer") == 1
    assert any(record.path == str(child) and record.status == "failed"
               for record in session.records)


def test_parent_operation_refuses_nested_resource_changed_after_child_read(tmp_path):
    parent = _file(tmp_path, "parent")
    child = _file(tmp_path, "child")

    def outer():
        result = validate_once([child], "inner", lambda: "old child result")
        child.write_text("changed child\n", encoding="utf-8")
        return result

    with InputValidationSession():
        with pytest.raises(RuntimeError, match="changed after pipeline preflight"):
            validate_once([parent], "outer", outer)


def test_module_scope_restored_after_exception_and_nested_session(tmp_path):
    source = _file(tmp_path)
    with InputValidationSession() as outer, outer.scope("parent"):
        with pytest.raises(ValueError):
            with outer.scope("child"):
                raise ValueError("test")
        record_file_validation(source, "After failure")
        with InputValidationSession() as inner:
            record_file_validation(source, "Other invocation")
        record_file_validation(source, "After nested session")
    assert [record.consumers for record in outer.records] == [
        ("parent",), ("parent",),
    ]
    assert inner.records[0].consumers == ()
    assert current_validation_session() is None


def test_pending_and_failed_checks_are_not_overwritten_by_availability_pass(tmp_path):
    source = _file(tmp_path)
    with InputValidationSession() as session, session.scope("consumer"):
        record_file_validation(None, "Generated association table", status="deferred")
        record_file_validation(source, "Compatibility", status="failed")
        require_nonempty_file(source, "Reference")
    assert [record.status for record in session.records] == [
        "deferred", "failed", "passed",
    ]


def test_report_metrics_are_json_safe_and_detached(tmp_path):
    source = _file(tmp_path)
    metrics = {"count": 2, "path": source, "columns": ["SNP", "BETA"]}
    with InputValidationSession() as session:
        record_file_validation(source, "Table", metrics=metrics)
    metrics["columns"].append("unrelated")
    report = session.to_dict()
    assert json.loads(json.dumps(report)) == report
    assert report["files"][0]["metrics"]["path"] == str(source)
    assert report["files"][0]["metrics"]["columns"] == ["SNP", "BETA"]
    session.records[0].metrics["columns"].append("cannot mutate evidence")
    report["files"][0]["metrics"]["count"] = 100
    assert session.records[0].metrics["columns"] == ["SNP", "BETA"]
    assert session.records[0].metrics["count"] == 2


@pytest.mark.parametrize("value", [float("nan"), float("inf"), object(), {1: "x"}])
def test_unsafe_or_nonfinite_report_values_are_rejected(value):
    with InputValidationSession(), pytest.raises((ValueError, TypeError)):
        record_file_validation(None, "Result", metrics={"value": value})


def test_unknown_status_and_inactive_or_empty_scope_rejected():
    session = InputValidationSession()
    with pytest.raises(RuntimeError, match="active session"):
        with session.scope("module"):
            pass
    with session:
        with pytest.raises(ValueError, match="module name"):
            with session.scope(" "):
                pass
        with pytest.raises(ValueError, match="Unknown input-validation status"):
            record_file_validation(None, "Result", status="valid enough")


def test_completed_session_cannot_carry_cache_into_another_invocation():
    session = InputValidationSession()
    with session:
        pass
    with pytest.raises(RuntimeError, match="fresh input-validation session"):
        with session:
            pass


def test_no_file_identity_means_no_cached_validation():
    calls = []
    with InputValidationSession():
        for _ in range(2):
            validate_once([], "configuration check", lambda: calls.append(True))
    assert calls == [True, True]


def test_moved_identity_api_keeps_paths_and_failure_behavior(tmp_path):
    source = _file(tmp_path)
    identity = capture_preflight_file_identities([source, source])
    assert len(identity) == 1
    assert require_unchanged_preflight_files(identity) == (source,)
    source.write_text("different\n", encoding="utf-8")
    with pytest.raises(ValueError, match="changed after pipeline preflight"):
        require_unchanged_preflight_files(identity, error_type=ValueError)
