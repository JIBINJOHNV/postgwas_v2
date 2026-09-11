"""Resource reports preserve exact check scope and reuse only current evidence."""

from argparse import Namespace
import gzip
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from postgwas.config import load_configuration
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.interval_validation import validate_bed4_file
from postgwas.core.preflight import PipelinePreflightEvidence
from postgwas.core.validation_reporting import FileValidationDisplay
from postgwas.modules.ld_annotation.service import preflight_ld_annotation
from postgwas.modules.pops import service as pops_service
from postgwas.pipeline.registry import REGISTRY
from postgwas.pipeline.resource_validation import (
    RESOURCE_REPORTERS,
    record_pipeline_resource_validation,
)
from preflight_support import pipeline_input_vcf_evidence


def _bed(path, rows):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(rows)
    return path


def _pops_args(tmp_path):
    prefix = tmp_path / "features"
    genes = ["ENSG%03d" % index for index in range(10)]
    Path(str(prefix) + ".rows.txt").write_text("\n".join(genes) + "\n")
    Path(str(prefix) + ".cols.0.txt").write_text("feature_a\nfeature_b\n")
    np.save(str(prefix) + ".mat.0.npy", np.arange(20, dtype=float).reshape(10, 2))
    annotation = tmp_path / "genes.tsv"
    annotation.write_text(
        "ENSGID\tCHR\tTSS\n" + "".join("%s\t1\t%d\n" % (gene, i) for i, gene in enumerate(genes))
    )
    return Namespace(
        genome_build="GRCh37", feature_matrix_prefix=str(prefix),
        feature_matrix_chunks=1, pops_gene_location_file=str(annotation),
    )


def test_all_registered_pipeline_preflights_have_explicit_resource_reporters():
    expected = {
        name for name in REGISTRY.names(include_internal=True)
        if REGISTRY.get(name).pipeline_enabled and REGISTRY.get(name).runner is not None
    }
    assert set(RESOURCE_REPORTERS) == expected


def test_bed4_all_rows_and_legal_zero_length_interval(tmp_path):
    path = _bed(tmp_path / "blocks.bed.gz", "# comment\n1\t0\t10\tfirst\n1\t10\t10\tsecond\n")
    with InputValidationSession() as session, session.scope("annot_ldblock"):
        validated = validate_bed4_file(path)
    assert validated.block_count == 2
    assert validated.contigs == frozenset({"1"})
    assert any(record.metrics.get("rows") == 2 for record in session.records)


@pytest.mark.parametrize("row", (
    "1\tbad\t10\ta\n", "1\t1.5\t10\ta\n", "1\t-1\t10\ta\n",
    "1\t1_000\t2000\ta\n",
    "1\t11\t10\ta\n", "1\t0\t10\t\n", "1\t0\t10\n", "# only a comment\n",
))
def test_bed4_rejects_invalid_rows_and_reports_file(tmp_path, row):
    path = _bed(tmp_path / "invalid.bed.gz", row)
    with InputValidationSession() as session, session.scope("annot_ldblock"):
        with pytest.raises(ValueError):
            validate_bed4_file(path)
    assert any(record.path == str(path) and record.status == "failed" for record in session.records)


def test_bed4_rejects_truncated_compression(tmp_path):
    path = _bed(tmp_path / "truncated.bed.gz", "1\t0\t10\ta\n")
    path.write_bytes(path.read_bytes()[:-8])
    with pytest.raises(ValueError, match="Invalid gzip BED4"):
        validate_bed4_file(path)


def test_bed4_pipeline_reuses_scan_but_direct_revalidates(monkeypatch, tmp_path):
    path = _bed(tmp_path / "blocks.bed.gz", "1\t0\t10\ta\n")
    real_open = gzip.open
    calls = []

    def counted(*args, **kwargs):
        calls.append(args[0])
        return real_open(*args, **kwargs)

    monkeypatch.setattr(gzip, "open", counted)
    with InputValidationSession() as session:
        with session.scope("annot_ldblock"):
            validate_bed4_file(path)
        with session.scope("consumer"):
            validate_bed4_file(path)
    assert len(calls) == 1
    content = next(record for record in session.records if "rows" in record.metrics)
    assert content.consumers == ("annot_ldblock", "consumer")
    validate_bed4_file(path)
    validate_bed4_file(path)
    assert len(calls) == 3


def test_ld_annotation_startup_checks_exact_beds_and_naming(tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"already validated VCF")
    bed = _bed(tmp_path / "GRCh37_EUR_ldetect.bed.gz", "chr1\t0\t10\ta\n")
    args = Namespace(vcf=str(vcf), ld_region_dir=str(tmp_path), ld_block_populations=["EUR"])
    with InputValidationSession() as session, session.scope("annot_ldblock"):
        with pytest.raises(ValueError, match="Chromosome names are incompatible"):
            preflight_ld_annotation(args, preflight_evidence=pipeline_input_vcf_evidence())
    assert any(record.path == str(bed) and record.status == "failed" for record in session.records)


@pytest.mark.parametrize("kind", ("imputation", "heritability", "caldera"))
def test_availability_is_not_misreported_as_content_validation(tmp_path, kind):
    path = tmp_path / "opaque.resource"
    path.write_bytes(b"content not parsed by current module preflight")
    if kind == "imputation":
        resources = SimpleNamespace(reference_files=(path,), pred_ld_script=path)
    elif kind == "heritability":
        resources = SimpleNamespace(reference=SimpleNamespace(merge_alleles=path, required_files=(path,)))
    else:
        resources = (None, {key: path for key in ("coding_variants", "gene_locations", "model", "upstream_script", "adapter")})
    evidence = PipelinePreflightEvidence(kind, {}, resources)
    with InputValidationSession() as session, session.scope(kind):
        record_pipeline_resource_validation(kind, evidence)
    assert session.records
    assert all(record.checks == ("regular nonempty file",) for record in session.records)
    assert all("not validated" in record.message for record in session.records)


def test_pops_resource_scans_reused_with_truthful_bundle_report(
    monkeypatch, tmp_path, capsys,
):
    args = _pops_args(tmp_path)
    calls = []
    from postgwas.core import matrix_validation
    original = matrix_validation.matrix_values_are_finite

    def inspect(*values):
        calls.append(values)
        return original(*values)

    monkeypatch.setattr(matrix_validation, "matrix_values_are_finite", inspect)
    with InputValidationSession() as session, session.scope("pops"):
        display = FileValidationDisplay(session, load_configuration())
        evidence = pops_service.preflight_pops_pipeline(
            args, preflight_evidence=pipeline_input_vcf_evidence(),
        )
        record_pipeline_resource_validation("pops", evidence)
        _, features, annotation, _ = pops_service.validate_pops_configuration(
            args, pipeline=True,
        )
        assert features["row_count"] == annotation["gene_count"] == 10
        display.flush()
    assert len(calls) == 1
    matrix = next(record for record in session.records if record.role == "PoPS feature matrix")
    assert matrix.metrics["rows"] == 10
    assert "all values: finite numeric data" in matrix.checks
    screen = capsys.readouterr().out
    assert "PoPS feature-matrix resources — CHECKS PASSED" in screen
    assert "features.cols.0.txt" not in screen
    assert "features.mat.0.npy" not in screen
    json.dumps(session.to_dict())


def test_pops_changed_file_cannot_reuse_preflight(tmp_path):
    args = _pops_args(tmp_path)
    with InputValidationSession() as session, session.scope("pops"):
        pops_service.preflight_pops_pipeline(args, preflight_evidence=pipeline_input_vcf_evidence())
        matrix = Path(args.feature_matrix_prefix + ".mat.0.npy")
        matrix.write_bytes(b"changed matrix")
        with pytest.raises(RuntimeError, match="changed"):
            pops_service.validate_pops_configuration(args, pipeline=True)


def test_pops_changed_schema_cannot_reuse_preflight(tmp_path):
    args = _pops_args(tmp_path)
    with InputValidationSession() as session, session.scope("pops"):
        pops_service.preflight_pops_pipeline(args, preflight_evidence=pipeline_input_vcf_evidence())
        configuration = pops_service._resolved_configuration(args)
        configuration.modules.pops.input_schema.gene_annotation_tss_column = "absent_tss"
        with pytest.raises(pops_service.PopsError, match="absent_tss"):
            pops_service._validate_resolved_pops_configuration(args, configuration, pipeline=True)
