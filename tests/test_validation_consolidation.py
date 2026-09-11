"""Whole-invocation validation facts, including incremental checks and audit history."""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import yaml
import pytest

from postgwas.config import load_configuration
from postgwas.core.input_validation import InputValidationSession, record_file_validation
from postgwas.core.validation_reporting import (
    FileValidationDisplay, flush_file_validation_display, run_with_file_validation,
    write_validation_audit, consolidate_validation_fields, file_validation_display,
    register_file_validation_bundle,
)
from postgwas.pipeline.cli import _print_pipeline_preflight_summary


def test_incremental_file_checks_do_not_repeat_filename_or_established_facts(tmp_path):
    source = tmp_path / "study.vcf.gz"
    source.write_text("fixture")

    def operation():
        record_file_validation(source, "Study", checks=("header",), metrics={"variants_from_index": 4})
        flush_file_validation_display()
        record_file_validation(source, "Study", checks=("header", "allele compatibility"),
                               metrics={"variants_from_index": 4, "sample": "STUDY"})
        flush_file_validation_display()
        record_file_validation(source, "Study", status="warning", message="Ancestry is a declaration")

    screen = StringIO()
    with redirect_stdout(screen):
        run_with_file_validation(operation, command="example", configuration=load_configuration(),
                                 output_directory=tmp_path / "output")
    text = screen.getvalue()
    assert text.count("study.vcf.gz") == 1
    assert text.count("Variants (from index)") == 1
    assert text.count("header") == 1
    assert "allele compatibility" in text and "STUDY" in text
    assert "Ancestry is a declaration" in text
    audit = yaml.safe_load((tmp_path / "output/run_metadata/example_input_validation.yaml").read_text())
    assert audit["report_version"] == 2
    assert len(audit["files"]) == 1
    assert audit["files"][0]["status"] == "warning"
    assert audit["files"][0]["checks"] == ["header", "allele compatibility"]
    assert audit["files"][0]["metrics"]["variants_from_index"] == 4


def test_combined_audit_preserves_conflicting_metrics_and_check_associations(tmp_path):
    source = tmp_path / "data"
    with InputValidationSession() as session:
        record_file_validation(source, "First contract", checks=("all rows",), metrics={"rows": 10})
        record_file_validation(source, "Subset contract", checks=("eligible rows",), metrics={"rows": 6})
    destination = tmp_path / "audit.yaml"
    write_validation_audit({"report_kind": "test", **session.to_dict()}, session.records, destination)
    records = yaml.safe_load(destination.read_text())["files"]
    assert len(records) == 1
    assert records[0]["metrics"] == {"rows": 10, "rows#2": 6}
    assert records[0]["observations"][1]["metric_keys"] == ["rows#2"]


def test_readiness_does_not_reprint_file_details(capsys):
    _print_pipeline_preflight_summary({"indexed": object()}, {}, ())
    text = capsys.readouterr().out
    assert "Pipeline readiness" in text and "Stage preflights" in text
    assert "Genome build" not in text and "Input GWAS-VCF" not in text


def test_new_output_gets_own_record_and_failed_check_is_not_hidden(tmp_path, capsys):
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        for filename, role in (("input.vcf.gz", "Input VCF"), ("output.vcf.gz", "Generated VCF")):
            record_file_validation(tmp_path / filename, role, checks=("header",))
            display.flush()
        record_file_validation(tmp_path / "output.vcf.gz", "Generated VCF", status="failed", message="Count changed")
        display.flush()
    text = capsys.readouterr().out
    assert text.count("input.vcf.gz") == text.count("output.vcf.gz") == 1
    assert "Count changed" in text and "FAILED" in text


def test_module_fields_join_file_card_and_leave_scientific_step_values_intact(tmp_path, capsys):
    from postgwas.core.pipeline_logging import PipelineLogger
    from postgwas.core.ui import StageProgress

    source = tmp_path / "study.vcf.gz"
    with InputValidationSession() as session, file_validation_display(session, load_configuration()):
        record_file_validation(source, "Input VCF", checks=("index",), metrics={"variants_from_index": 4})
        with PipelineLogger("STUDY", "run", str(tmp_path),
                            stage_progress=StageProgress("Example", enabled=True)) as logger:
            with logger.step(1, 1, "Validate input", "validate_input") as step:
                step.outcome("Input checked", fields=[("info", "Input file", str(source)),
                             ("count", "Total variants", 4), ("warning", "Population", "Declared only")],
                             input_vcf=str(source), variants_from_index=4)
            assert logger.summary()["steps"][0]["extra"]["outcome"]["input_vcf"] == str(source)
            text_log = Path(logger.log_path).read_text()
            assert str(source) not in text_log
            assert "variants_from_index" not in text_log
    text = capsys.readouterr().out
    assert text.count("study.vcf.gz") == 1
    assert text.count("Variants (from index)") == 1
    assert "Declared only" in text


def test_ambiguous_filename_and_cross_file_fields_are_not_guessed(tmp_path):
    fields = [("info", "Input file", "study.tsv"), ("count", "Genes", 4)]
    with InputValidationSession() as session, file_validation_display(session, load_configuration()):
        for folder in ("one", "two"):
            record_file_validation(tmp_path / folder / "study.tsv", "Table", metrics={"genes": 4})
        assert consolidate_validation_fields(fields) == fields
        cross_file = [("info", "First", str(tmp_path / "one/study.tsv")),
                      ("info", "Second", str(tmp_path / "two/study.tsv")),
                      ("count", "Shared genes", 2)]
        combined = consolidate_validation_fields(cross_file)
        assert combined[-1] == ("count", "Shared genes", 2)
        assert combined[0] == ("info", "First", "File 1")
        assert combined[1] == ("info", "Second", "File 2")


def test_direct_reuse_does_not_enable_pipeline_contracts_and_changed_files_revalidate(tmp_path):
    from postgwas.core.input_validation import current_validation_session, validate_once

    source = tmp_path / "table"
    source.write_text("first")
    calls = []

    def inspect():
        assert current_validation_session() is None
        calls.append(source.read_text())
        return calls[-1]

    with InputValidationSession(pipeline=False):
        assert validate_once((source,), "schema-one", inspect) == "first"
        assert validate_once((source,), "schema-one", inspect) == "first"
        validate_once((source,), "schema-two", inspect)
        source.write_text("new version")
        assert validate_once((source,), "schema-one", inspect) == "new version"
    assert calls == ["first", "first", "new version"]


def test_direct_nested_dependency_change_invalidates_parent_and_failures_are_not_cached(tmp_path):
    from postgwas.core.input_validation import validate_once

    manifest, data = tmp_path / "manifest", tmp_path / "data"
    manifest.write_text("data")
    data.write_text("valid")
    calls = []

    def read():
        calls.append(1)
        if data.read_text() == "invalid":
            raise ValueError("invalid child")
        return data.read_text()

    def parent():
        return validate_once((data,), "data", read)

    with InputValidationSession(pipeline=False):
        assert validate_once((manifest,), "parent", parent) == "valid"
        assert validate_once((manifest,), "parent", parent) == "valid"
        data.write_text("invalid")
        for _ in range(2):
            with pytest.raises(ValueError, match="invalid child"):
                validate_once((manifest,), "parent", parent)
    assert len(calls) == 3


def test_direct_missing_cached_dependency_keeps_original_validator_error(tmp_path):
    from postgwas.core.input_validation import validate_once

    manifest, child = tmp_path / "manifest", tmp_path / "child"
    manifest.write_text("child")
    child.write_text("data")

    def read():
        if not child.is_file():
            raise ValueError("Original child-file error")
        return child.read_text()

    with InputValidationSession(pipeline=False):
        validate_once((manifest,), "parent", lambda: validate_once((child,), "child", read))
        child.unlink()
        with pytest.raises(ValueError, match="Original child-file error"):
            validate_once((manifest,), "parent", lambda: validate_once((child,), "child", read))


def test_companion_updates_stay_with_the_parent_file(tmp_path, capsys):
    source = tmp_path / "study.vcf.gz"
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        record_file_validation(source, "Study", checks=("header",))
        record_file_validation(str(source) + ".csi", "Index", checks=("availability",))
        assert display.file_reference(str(source)) == "File 1"
        display.flush()
        record_file_validation(str(source) + ".csi", "Index", status="failed", message="Index is stale")
        display.flush()
    output = capsys.readouterr().out
    assert output.count("study.vcf.gz") == 1
    assert "Index is stale" in output and "findings for file 1" in output


def test_structured_validation_logs_combine_without_filtering_native_diagnostics(tmp_path):
    from postgwas.core.pipeline_logging import PipelineLogger

    source = tmp_path / "study.vcf.gz"
    with InputValidationSession() as session, file_validation_display(session, load_configuration()):
        record_file_validation(source, "Study", checks=("header",))
        with PipelineLogger("STUDY", "run", str(tmp_path)) as logger:
            for _ in range(2):
                logger.record("VALIDATE", "contract", input_vcf=str(source), status="PASSED", sample="S")
            logger.info("native diagnostic [literal]: " + str(source))
            log = Path(logger.log_path).read_text()
        assert log.count("contract") == 1
        assert log.count(str(source)) == 1  # Native diagnostic remains literal.
        assert "native diagnostic [literal]" in log
        assert any(record.metrics.get("logged_validation", {}).get("sample") == "S"
                   for record in session.records)


def test_nonfinite_analysis_details_keep_original_log_representation(tmp_path):
    fields = [("info", "Input", str(tmp_path / "file"))]
    details = {"diagnostic": float("nan")}
    with InputValidationSession() as session, file_validation_display(session, load_configuration()):
        record_file_validation(tmp_path / "file", "Table", checks=("header",))
        consolidate_validation_fields(fields, log_details=details)
    assert "diagnostic" in details


@pytest.mark.parametrize("renderer", ["filtering", "qc"])
def test_module_vcf_cards_consolidate_metadata_but_preserve_provenance(tmp_path, capsys, renderer):
    source = tmp_path / "study.vcf.gz"
    with InputValidationSession() as session, file_validation_display(session, load_configuration()):
        record_file_validation(source, "Study", checks=("header",), metrics={
            "declared_genome_build": "GRCh37", "contigs": ["1", "2"], "variants_from_index": 4,
        })
        flush_file_validation_display()
        if renderer == "filtering":
            from postgwas.modules.filtering.reporting import render_input_vcf_validation
            rendered = render_input_vcf_validation({
                "input_vcf": str(source), "status": "PASS", "genome_build": "GRCh37",
                "declared_contigs": 2, "total_variants": 4,
                "field_provenance": {"status": "AVAILABLE", "detail": "Study supplied AF"},
            }, label_width=42)
        else:
            from postgwas.modules.qc_summary.reporting import qc_input_validation_screen_lines
            rendered = qc_input_validation_screen_lines({
                "raw_vcf": str(source), "raw": {"num_records": 4},
                "vcf_header_validation": {"status": "passed", "genome_build": "GRCh37", "declared_contig_count": 2},
                "vcf_provenance": {"status": "available", "values": {"info_source": "Study supplied INFO"}},
            })
        assert not rendered
    output = capsys.readouterr().out
    assert output.count("study.vcf.gz") == 1
    assert output.count("GRCh37") == 1
    assert output.count("Variants (from index)") == 1
    assert "Study supplied" in output


def test_flames_runtime_probe_reused_with_native_messages_preserved(tmp_path, monkeypatch):
    import sys
    from postgwas.core.preflight import PreflightLogBuffer
    from postgwas.modules.flames import service

    config = load_configuration(cli_overrides={"resources.executables.python": sys.executable})
    calls = []

    def command(arguments, purpose, *, logger, **kwargs):
        calls.append(arguments)
        logger.info("native runtime diagnostic")
        return ""

    monkeypatch.setattr(service, "run_checked_command", command)
    with InputValidationSession():
        service._validate_runtime(config, config.modules.flames)
        logger = PreflightLogBuffer()
        service._validate_runtime(config, config.modules.flames, logger=logger)
        assert len(calls) == 1
        assert logger.events[0].arguments == ("native runtime diagnostic",)
        monkeypatch.setenv("PYTHONPATH", str(tmp_path))
        service._validate_runtime(config, config.modules.flames)
        assert len(calls) == 2


def test_manhattan_runtime_probe_reuses_exact_contract_only(tmp_path, monkeypatch):
    import sys
    from types import SimpleNamespace
    from postgwas.modules.manhattan import service

    source = tmp_path / "study.vcf.gz"
    source.write_text("fixture")
    indexed = SimpleNamespace(vcf=source, sample="STUDY", header=(
        '##FORMAT=<ID=LP,Number=1,Type=Float,Description="p">\n'
        '##FORMAT=<ID=AF,Number=1,Type=Float,Description="af">\n'
    ))
    monkeypatch.setattr(service, "validate_indexed_vcf", lambda *args, **kwargs: indexed)
    monkeypatch.setattr(service, "resolve_executable", lambda *args, **kwargs: sys.executable)
    calls = []

    def run(arguments, **kwargs):
        calls.append(arguments)
        return SimpleNamespace(stdout="runtime checked")

    monkeypatch.setattr(service.subprocess, "run", run)
    args = SimpleNamespace(vcf=source, dataset_id="STUDY")
    with InputValidationSession():
        service.resolve_plot_inputs(args)
        service.resolve_plot_inputs(args)
        assert len(calls) == 1
        monkeypatch.setenv("R_LIBS_USER", str(tmp_path))
        service.resolve_plot_inputs(args)
        assert len(calls) == 2


def test_unknown_display_alias_is_rejected_by_schema():
    from postgwas.core.errors import ConfigurationError
    with pytest.raises(ConfigurationError, match="outcome_metric_aliases"):
        load_configuration(cli_overrides={
            "logging.file_validation.outcome_metric_aliases": {"Wrong label": "unregistered_metric"},
        })


def test_screen_limit_never_hides_later_failure(tmp_path, capsys):
    configuration = load_configuration(cli_overrides={"logging.file_validation.max_screen_files": 1})
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, configuration)
        for name in ("one", "two"):
            record_file_validation(tmp_path / name, "Reference", checks=("header",))
        display.flush()
        record_file_validation(tmp_path / "two", "Reference", status="failed", message="Invalid header")
        display.flush()
    assert "Invalid header" in capsys.readouterr().out


def test_large_pops_feature_bundle_is_compact_on_screen_and_complete_in_audit(
    tmp_path, capsys,
):
    configuration = load_configuration()
    limit = configuration.logging.file_validation.max_screen_files
    assert limit == 20
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, configuration)
        prefix = tmp_path / "pops_features"
        record_file_validation(
            str(prefix) + ".rows.txt",
            "PoPS feature gene-order file",
            checks=("nonempty unique names",),
            metrics={"names": 18_383},
        )
        for chunk in range(116):
            record_file_validation(
                "%s.cols.%d.txt" % (prefix, chunk),
                "PoPS feature-column names",
                checks=("nonempty unique names",),
                metrics={"names": 500},
            )
            record_file_validation(
                "%s.mat.%d.npy" % (prefix, chunk),
                "Feature matrix",
                checks=("companion dimensions", "all values finite and numeric"),
                metrics={"rows": 18_383, "columns": 500, "dtype": "float32"},
            )
        bundle_paths = tuple(
            path
            for chunk in range(116)
            for path in (
                "%s.cols.%d.txt" % (prefix, chunk),
                "%s.mat.%d.npy" % (prefix, chunk),
            )
        )
        register_file_validation_bundle(
            bundle_paths,
            "PoPS feature matrix bundle",
            (
                ("count", "pops_feature_chunks", 116),
                ("success", "pops_feature_column_files", "116 / 116"),
                ("success", "pops_feature_matrix_files", "116 / 116"),
                ("count", "genes", 18_383),
                ("count", "features", 58_000),
                ("count", "pops_matrix_dtypes", ["float32"]),
                ("count", "pops_unique_bundle_files", 232),
                (
                    "success",
                    "pops_matrix_validation",
                    "unique names; dimensions matched; finite numeric values",
                ),
            ),
            covered_checks=(
                "nonempty unique names",
                "companion dimensions",
                "all values finite and numeric",
            ),
            covered_metric_keys=("names", "rows", "columns", "dtype"),
        )
        display.flush()
        destination = tmp_path / "input_validation.yaml"
        write_validation_audit(
            {"report_kind": "test", "status": "passed"},
            session.records,
            destination,
        )

    screen = capsys.readouterr().out
    assert screen.count("— CHECKS PASSED") == 2
    assert "Files 2–233 · PoPS feature-matrix resources — CHECKS PASSED" in screen
    assert "Matrix chunks" in screen and "116" in screen
    assert "Physical files checked" in screen and "232" in screen
    assert "Additional files" not in screen
    assert "pops_features.mat.115.npy" not in screen
    audit = yaml.safe_load(destination.read_text())
    assert len(audit["files"]) == 233
    assert any(
        row["path"].endswith("pops_features.mat.115.npy")
        for row in audit["files"]
    )


def test_direct_record_count_is_not_mislabelled_as_an_index_count(tmp_path, capsys):
    source = tmp_path / "study.vcf.gz"
    with InputValidationSession(pipeline=False) as session, file_validation_display(session, load_configuration()):
        record_file_validation(source, "Study", checks=("record count",), metrics={"total_variants": 4})
        consolidate_validation_fields([("info", "Input", str(source)), ("count", "Total variants", 4)],
                                      log_details={"total_variants": 4})
        flush_file_validation_display()
    output = capsys.readouterr().out
    assert "Variants (from index)" not in output
    assert output.count("Variants") == 1
    assert all("variants_from_index" not in record.metrics for record in session.records)


def test_incremental_long_list_change_and_same_basename_remain_distinguishable(tmp_path, capsys):
    configuration = load_configuration(cli_overrides={"logging.file_validation.max_list_items": 2})
    first, second = tmp_path / "one/table.tsv", tmp_path / "two/table.tsv"
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, configuration)
        record_file_validation(first, "Table", metrics={"columns": ["A", "B", "C"]})
        display.flush()
        capsys.readouterr()
        record_file_validation(first, "Table", metrics={"columns": ["A", "B", "D"]})
        record_file_validation(second, "Table", checks=("header",))
        display.flush()
    output = capsys.readouterr().out
    assert "Columns" in output
    assert str(second) in "".join(output.split())  # Shared path wrapping remains enabled.
    assert "findings for file 1" in output


@pytest.mark.parametrize("record_count", [4, 3])
def test_index_and_record_counts_combine_only_when_they_agree(tmp_path, capsys, record_count):
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        source = tmp_path / "study.vcf.gz"
        record_file_validation(source, "Index check", metrics={"variants_from_index": 4})
        display.flush()
        record_file_validation(source, "Record check", metrics={"total_variants": record_count})
        display.flush()
    output = capsys.readouterr().out
    assert output.count("Variants") == (1 if record_count == 4 else 2)
    assert output.count("Variants (from index)") == 1


def test_enrichment_worker_leaves_gene_count_to_service_without_changing_genes(tmp_path, monkeypatch):
    """Exercise the worker entry function without starting its isolated runtime."""
    import ast
    import sys
    from types import SimpleNamespace
    from postgwas.modules.enrichment.utils import load_gene_list

    source = Path(__file__).resolve().parents[1] / "src/postgwas/modules/enrichment/main.py"
    tree = ast.parse(source.read_text())
    main = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main")
    genes = tmp_path / "genes.tsv"
    genes.write_text("GENE\nTP53\nBRCA1\nTP53\n")
    args = SimpleNamespace(gene_input_file=genes, output_directory=tmp_path / "out", dataset_id="STUDY",
                           biogrid_key=None, david_email=None, dsigdb_gmt=None, reference_set="reference")
    calls, messages = [], []
    monkeypatch.setitem(sys.modules, "postgwas.modules.enrichment.service", SimpleNamespace(
        run_multisource_enrichment_pipeline=lambda **values: calls.append(values),
    ))
    namespace = {
        "Path": Path, "sys": sys, "os": SimpleNamespace(environ={}),
        "get_geneset_parser": lambda **kwargs: SimpleNamespace(parse_args=lambda: args),
        "_configured_enrichment_python": lambda: Path(sys.executable).resolve(),
        "_enrichment_runtime_environment": lambda target: {},
        "print_screen_message": lambda kind, message, **kwargs: messages.append((kind, message)),
    }
    exec(compile(ast.Module(body=[main], type_ignores=[]), str(source), "exec"), namespace)
    assert namespace["main"]() == 0
    assert len(calls) == 1 and calls[0]["gene_list"] == load_gene_list(genes) == ["BRCA1", "TP53"]
    assert not any(kind == "count" for kind, _ in messages)
