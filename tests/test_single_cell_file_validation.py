"""Shared file reads are reusable; generated-input compatibility is not skipped."""

import json
from pathlib import Path

import pytest

from postgwas.config import load_configuration
from postgwas.core import ldsc_validation, single_cell_validation
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.preflight import PipelinePreflightEvidence
from postgwas.core.validation_reporting import FileValidationDisplay
from postgwas.modules.single_cell.errors import SingleCellError
from postgwas.modules.single_cell.methods.scdrs import runner as scdrs
from postgwas.modules.single_cell.methods.ldsc_celltype import runner as ldsc
from postgwas.pipeline.resource_validation import record_pipeline_resource_validation
from test_single_cell import (
    _fake_ldsc, _fake_scdrs, _gene_sets, _h5ad, _identifier_map,
    _ldsc_inputs, _magma_gene_statistics, _scdrs_covariates,
)


def _count_calls(monkeypatch, module, name):
    original = getattr(module, name)
    calls = []

    def wrapped(*args, **kwargs):
        calls.append(True)
        return original(*args, **kwargs)

    monkeypatch.setattr(module, name, wrapped)
    return calls


def test_scdrs_static_inputs_read_once_across_pipeline_and_consumer_preflights(
    tmp_path, monkeypatch,
):
    method = load_configuration().modules.single_cell.scdrs.model_copy(deep=True)
    method.magma_gene_set.source = "magma"
    method.input.h5ad_file = _h5ad(tmp_path / "atlas.h5ad", gene_count=1000)
    method.input.gene_identifier_map_file = _identifier_map(tmp_path / "map.tsv")
    method.input.covariate_file = _scdrs_covariates(tmp_path / "cells.cov")
    executable = _fake_scdrs(tmp_path / "scdrs")
    reads = {
        name: _count_calls(monkeypatch, single_cell_validation, name)
        for name in ("_load_anndata", "_inspect_identifier_crosswalk", "_inspect_cell_covariates")
    }
    with InputValidationSession() as session:
        with session.scope("startup"):
            startup = scdrs.preflight_scdrs(
                method, executable, dataset_id="study", pipeline_pending_magma=True,
            )
        method.input.magma_gene_results_file = _magma_gene_statistics(tmp_path / "study.genes.out")
        with session.scope("consumer"):
            ready = scdrs.preflight_scdrs(method, executable, dataset_id="study")
    assert all(len(calls) == 1 for calls in reads.values())
    assert startup.h5ad_summary == ready.h5ad_summary
    assert startup.gene_universe == ready.gene_universe
    content = [record for record in session.records if "complete X-matrix value scan" in record.checks]
    assert len(content) == 1
    assert content[0].consumers == ("startup", "consumer")
    assert "gene_universe" not in content[0].metrics
    assert "cell_names" not in content[0].metrics
    json.dumps(session.to_dict())


def test_h5ad_changed_policy_rescans_and_changed_file_preserves_typed_error(
    tmp_path, monkeypatch,
):
    method = load_configuration().modules.single_cell.scdrs.model_copy(deep=True)
    source = _h5ad(tmp_path / "atlas.h5ad", fractional=True)
    method.validation.raw_count_integer_tolerance = 1
    reads = _count_calls(monkeypatch, single_cell_validation, "_load_anndata")
    with InputValidationSession():
        scdrs.validate_scdrs_h5ad(source, method)
        method.validation.raw_count_integer_tolerance = 0
        with pytest.raises(SingleCellError, match="non-integer"):
            scdrs.validate_scdrs_h5ad(source, method)
        _h5ad(source)
        with pytest.raises(SingleCellError, match="changed after pipeline preflight"):
            scdrs.validate_scdrs_h5ad(source, method)
    assert len(reads) == 2


@pytest.mark.parametrize("last_value,message", [
    (1.0, None),
    (float("nan"), "NaN or infinite"),
    (float("inf"), "NaN or infinite"),
    (-1.0, "negative expression"),
    (0.5, "non-integer"),
])
def test_backed_sparse_h5ad_scans_values_through_last_chunk(
    tmp_path, last_value, message,
):
    anndata = pytest.importorskip("anndata")
    sparse = pytest.importorskip("scipy.sparse")
    source = _h5ad(tmp_path / "sparse.h5ad")
    atlas = anndata.read_h5ad(source)
    atlas.X = sparse.csr_matrix(atlas.X)
    atlas.X.data[-1] = last_value
    atlas.write_h5ad(source)
    method = load_configuration().modules.single_cell.scdrs.model_copy(deep=True)
    method.validation.matrix_chunk_rows = 7
    with InputValidationSession():
        if message is not None:
            with pytest.raises(SingleCellError, match=message):
                scdrs.validate_scdrs_h5ad(source, method)
        else:
            summary = scdrs.validate_scdrs_h5ad(source, method)
            assert summary["cells_after_configured_filter"] == 60
            assert summary["genes_after_configured_filter"] == 300


def test_gene_set_parser_reuses_file_but_rechecks_new_gene_universe(tmp_path, monkeypatch):
    method = load_configuration().modules.single_cell.scdrs
    source = _gene_sets(tmp_path / "study.gs")
    reads = _count_calls(monkeypatch, single_cell_validation, "_inspect_gene_sets")
    with InputValidationSession():
        result = scdrs.validate_scdrs_gene_sets(
            source, method, gene_universe=tuple("GENE%03d" % index for index in range(300)),
        )
        with pytest.raises(SingleCellError, match="configured minimum"):
            scdrs.validate_scdrs_gene_sets(source, method, gene_universe=("GENE000",))
    assert result["traits"][0]["effective_genes"] == 51
    assert len(reads) == 1


def test_cell_covariate_parser_reuses_file_but_rechecks_new_cell_set(tmp_path, monkeypatch):
    method = load_configuration().modules.single_cell.scdrs
    source = _scdrs_covariates(tmp_path / "cells.cov")
    reads = _count_calls(monkeypatch, single_cell_validation, "_inspect_cell_covariates")
    cells = tuple("cell_%03d" % index for index in range(60))
    with InputValidationSession():
        scdrs.validate_scdrs_covariates(source, method, cell_names=cells)
        with pytest.raises(SingleCellError, match="exactly match"):
            scdrs.validate_scdrs_covariates(source, method, cell_names=cells[:-1])
    assert len(reads) == 1


@pytest.mark.parametrize("contents,message", [
    ("ENTREZID\tSYMBOL\n1\tA\n1\tB\n", "not one-to-one"),
    ("ENTREZID\tSYMBOL\n1\t\n", "empty source or target"),
    ("wrong\tcolumns\n1\tA\n", "missing columns"),
])
def test_crosswalk_invalid_content_fails_at_startup(tmp_path, contents, message):
    source = tmp_path / "mapping.tsv"
    source.write_text(contents, encoding="utf-8")
    method = load_configuration().modules.single_cell.scdrs.model_copy(deep=True)
    method.input.h5ad_file = _h5ad(tmp_path / "atlas.h5ad")
    method.input.gene_identifier_map_file = source
    method.magma_gene_set.source = "magma"
    executable = _fake_scdrs(tmp_path / "scdrs")
    with InputValidationSession(), pytest.raises(SingleCellError, match=message):
        scdrs.preflight_scdrs(
            method, executable, dataset_id="study", pipeline_pending_magma=True,
        )


def test_ldsc_manifest_and_reference_inventory_reused_without_claiming_content_checks(
    tmp_path, monkeypatch, capsys,
):
    inputs = _ldsc_inputs(tmp_path)
    method = load_configuration().modules.single_cell.ldsc_celltype.model_copy(deep=True)
    method.input.ldcts_file = Path(inputs["ldcts"])
    method.input.sumstats_file = Path(inputs["sumstats"])
    method.input.baseline_ld_prefixes = [inputs["baseline"]]
    method.input.weights_ld_prefix = inputs["weights"]
    executable = _fake_ldsc(tmp_path / "ldsc")
    reads = _count_calls(monkeypatch, ldsc_validation, "_inspect_ldcts")
    headers = _count_calls(monkeypatch, ldsc_validation, "_inspect_munged_header")
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        with session.scope("startup"):
            first = ldsc.preflight_ldsc_celltype(method, executable, executable)
        with session.scope("consumer"):
            second = ldsc.preflight_ldsc_celltype(method, executable, executable)
        record_pipeline_resource_validation(
            "single_cell",
            PipelinePreflightEvidence(
                "single_cell",
                {},
                {"methods": {"ldsc_celltype": first}},
            ),
        )
        display.flush()
    assert len(reads) == len(headers) == 1
    assert first.reference_inventory == second.reference_inventory
    records = [record for record in session.records if record.role == "LDSC reference companion"]
    assert len(records) == 154
    assert all(record.consumers == ("startup", "consumer") for record in records)
    assert all(record.checks == ("regular file", "non-empty file") for record in records)
    assert all("Availability only" in record.message for record in records)
    screen = capsys.readouterr().out
    assert (
        "LDSC cell-type reference resources — AVAILABLE — availability only"
        in screen
    )
    assert "Physical files checked" in screen and "154" in screen
    assert not any(Path(record.path).name in screen for record in records)


def test_ldsc_manifest_changed_contract_does_not_reuse_old_result(tmp_path):
    source = tmp_path / "brain.ldcts"
    source.write_text("Neuron\ttest.,control.\n", encoding="utf-8")
    with InputValidationSession():
        ldsc_validation.read_ldcts_manifest(source, prefix_separator=",", minimum_prefixes_per_cell_type=2)
        with pytest.raises(SingleCellError, match="at least 3"):
            ldsc_validation.read_ldcts_manifest(
                source, prefix_separator=",", minimum_prefixes_per_cell_type=3,
                error_type=SingleCellError,
            )


@pytest.mark.parametrize("name,contents,message", [
    ("bad.sumstats.gz", b"not gzip", "Cannot read"),
    ("missing.sumstats", b"SNP Z\nrs1 2\n", "missing required"),
    ("empty.sumstats", b"SNP Z N\n", "no variants"),
])
def test_ldsc_header_failures_preserve_consumer_error(tmp_path, name, contents, message):
    source = tmp_path / name
    source.write_bytes(contents)
    with InputValidationSession(), pytest.raises(SingleCellError, match=message):
        ldsc_validation.validate_ldsc_sumstats_header(
            source, required_columns=("SNP", "Z", "N"), error_type=SingleCellError,
        )


def test_ldsc_header_report_does_not_claim_full_value_validation(tmp_path):
    source = tmp_path / "study.sumstats"
    source.write_text("SNP Z N\nrs1 bad bad\n", encoding="utf-8")
    with InputValidationSession() as session:
        ldsc_validation.validate_ldsc_sumstats_header(source, required_columns=("SNP", "Z", "N"))
    report = next(record for record in session.records if record.role == "Munged LDSC summary-statistic header")
    assert report.status == "passed"
    assert "remaining records are not scanned" in report.message
