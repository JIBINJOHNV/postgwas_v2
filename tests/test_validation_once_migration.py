"""Central file checks must not repeat scans or bypass consumer requirements."""

from pathlib import Path

import numpy as np
import pytest

from postgwas.core.input_validation import InputValidationSession, validation_scope, validate_once
from postgwas.core import vcf
from postgwas.core import matrix_validation as matrices
from postgwas.core import gene_annotation_validation as annotations
from postgwas.core import gene_property_validation as properties


def test_vcf_consumers_share_queries_but_still_require_their_own_fields(tmp_path, monkeypatch):
    path = tmp_path / "study.vcf.gz"
    path.write_bytes(b"fixture")
    Path(str(path) + ".tbi").write_bytes(b"index")
    header = (
        "##fileformat=VCFv4.2\n##genome_build=GRCh37\n"
        "##contig=<ID=1>\n"
        '##FORMAT=<ID=ES,Number=1,Type=Float,Description="effect">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY\n"
    )
    calls = []

    def query(command, *args, **kwargs):
        operation = tuple(command[1:3])
        calls.append(operation)
        return {("view", "--header-only"): header,
                ("index", "-n"): "3\n", ("query", "-l"): "STUDY\n"}[operation]

    monkeypatch.setattr(vcf, "run_checked_command", query)
    contract = dict(genome_build_header="##genome_build={build}",
                    supported_genome_builds=["GRCh37"], required_fields=["FORMAT/ES"])
    with InputValidationSession() as session:
        for consumer in ("formatter", "qc_summary", "ld_clump"):
            with validation_scope(consumer):
                assert vcf.validate_indexed_vcf(path, "STUDY", "bcftools", **contract).variant_count == 3
        with pytest.raises(vcf.VcfQueryError, match="FORMAT/MISSING"):
            vcf.validate_indexed_vcf(path, "STUDY", "bcftools",
                                     **{**contract, "required_fields": ["FORMAT/MISSING"]})
    assert calls == [("view", "--header-only"), ("index", "-n"), ("query", "-l")]
    queries = [r for r in session.records if r.role == "Reading VCF header"]
    assert queries[0].consumers == ("formatter", "qc_summary", "ld_clump")


def test_vcf_queries_reject_changed_inputs_and_use_fresh_run_cache(tmp_path, monkeypatch):
    path = tmp_path / "study.vcf.gz"
    path.write_bytes(b"fixture")
    calls = []
    monkeypatch.setattr(vcf, "run_checked_command", lambda *a, **k: calls.append(True) or "header\n")
    with InputValidationSession():
        assert vcf.read_vcf_header(path, "bcftools") == "header\n"
        path.write_bytes(b"changed fixture")
        with pytest.raises(RuntimeError, match="changed after pipeline preflight"):
            vcf.read_vcf_header(path, "bcftools")
    with InputValidationSession():
        vcf.read_vcf_header(path, "bcftools")
    assert calls == [True, True]


def test_index_change_cannot_reuse_cached_count(tmp_path, monkeypatch):
    path = tmp_path / "study.vcf.gz"
    path.write_bytes(b"fixture")
    index = Path(str(path) + ".tbi")
    index.write_bytes(b"index")
    monkeypatch.setattr(vcf, "run_checked_command", lambda *a, **k: "3\n")
    with InputValidationSession():
        assert vcf.count_indexed_vcf_records(path, "bcftools") == 3
        index.write_bytes(b"changed index")
        with pytest.raises(RuntimeError, match="changed after pipeline preflight"):
            vcf.count_indexed_vcf_records(path, "bcftools")


def _matrix_bundle(tmp_path, values):
    rows = tmp_path / "rows.txt"
    columns = tmp_path / "columns.txt"
    matrix = tmp_path / "matrix.npy"
    rows.write_text("gene1\ngene2\n")
    columns.write_text("feature1\nfeature2\n")
    np.save(matrix, np.asarray(values))
    return rows, ((columns, matrix),)


def test_matrix_bundle_scanned_once_with_detached_evidence(tmp_path, monkeypatch):
    rows, chunks = _matrix_bundle(tmp_path, [[1., 2.], [3., 4.]])
    original = matrices.matrix_values_are_finite
    calls = []
    monkeypatch.setattr(matrices, "matrix_values_are_finite",
                        lambda value: calls.append(True) or original(value))
    with InputValidationSession():
        evidence = matrices.validate_feature_matrix_bundle(rows, chunks)
        evidence["chunks"][0]["rows"] = 999
        second = matrices.validate_feature_matrix_bundle(rows, chunks)
        assert second["chunks"][0]["rows"] == 2
        assert not second["rows"].flags.writeable
    assert calls == [True]


@pytest.mark.parametrize("values,match", [
    ([[1., np.nan], [3., 4.]], "non-finite"),
    ([[1., np.inf], [3., 4.]], "non-finite"),
    ([[1., 2.]], "shape"),
    ([["a", "b"], ["c", "d"]], "numeric"),
])
def test_invalid_feature_matrices_fail(tmp_path, values, match):
    rows, chunks = _matrix_bundle(tmp_path, values)
    with pytest.raises(ValueError, match=match):
        matrices.validate_feature_matrix_bundle(rows, chunks)


def test_duplicate_matrix_names_and_changed_companions_fail(tmp_path):
    rows, chunks = _matrix_bundle(tmp_path, [[1., 2.], [3., 4.]])
    with InputValidationSession():
        matrices.validate_feature_matrix_bundle(rows, chunks)
        chunks[0][0].write_text("feature1\nfeature1\n")
        with pytest.raises(ValueError, match="changed"):
            matrices.validate_feature_matrix_bundle(rows, chunks)
    with pytest.raises(ValueError, match="duplicate"):
        matrices.validate_feature_matrix_bundle(rows, chunks)


def test_missing_or_empty_matrix_bundle_fails(tmp_path):
    with pytest.raises(ValueError, match="companion pairs"):
        matrices.validate_feature_matrix_bundle(tmp_path / "rows", ())
    with pytest.raises(ValueError, match="missing"):
        matrices.read_unique_names(tmp_path / "missing", "Names")


def test_kernel_reuse_and_stricter_symmetry_contract(tmp_path, monkeypatch):
    path = tmp_path / "kernel.bin"
    np.asarray([[1., 0.1001], [0.1, 1.]], dtype=np.float32).tofile(path)
    kwargs = dict(dimension=2, dtype=np.float32, bytes_per_value=4,
                  chunk_rows=1, relative_tolerance=0., absolute_tolerance=0.001,
                  label="Test kernel")
    original = np.memmap
    calls = []
    monkeypatch.setattr(matrices.np, "memmap", lambda *a, **k: calls.append(True) or original(*a, **k))
    with InputValidationSession():
        assert matrices.validate_binary_kernel(path, **kwargs)["dimension"] == 2
        matrices.validate_binary_kernel(path, **kwargs)
        assert calls == [True]
        with pytest.raises(ValueError, match="symmetric"):
            matrices.validate_binary_kernel(path, **{**kwargs, "absolute_tolerance": 0.000001})
    assert len(calls) == 2


@pytest.mark.parametrize("values,dimension,match", [
    ([1.], 2, "bytes"),
    ([[1., np.inf], [np.inf, 1.]], 2, "non-finite"),
    ([[1., 1.], [0., 1.]], 2, "symmetric"),
])
def test_invalid_binary_kernels_fail(tmp_path, values, dimension, match):
    path = tmp_path / "kernel.bin"
    np.asarray(values, dtype=np.float32).tofile(path)
    with pytest.raises(ValueError, match=match):
        matrices.validate_binary_kernel(path, dimension=dimension, dtype=np.float32,
            bytes_per_value=4, chunk_rows=1, relative_tolerance=0., absolute_tolerance=0., label="Kernel")


def test_annotation_reuses_structure_not_configured_chromosome_policy(tmp_path, monkeypatch):
    from postgwas.config import load_configuration
    from postgwas.modules.pops.service import _validate_gene_annotation, PopsError

    path = tmp_path / "annotation.tsv"
    path.write_text("ENSGID\tCHR\tTSS\ng1\t1\t0\ng2\t2\t20\n")
    module = load_configuration().modules.pops
    module.gene_location_file = str(path)
    module.covariate_projection_chromosomes = ["1"]
    module.feature_selection_chromosomes = ["1"]
    module.training_chromosomes = ["2"]
    original = annotations.read_pandas_table
    calls = []
    monkeypatch.setattr(annotations, "read_pandas_table", lambda *a, **k: calls.append(True) or original(*a, **k))
    with InputValidationSession():
        first = _validate_gene_annotation(module)
        first["gene_names"]["g1"] = "changed by consumer"
        assert _validate_gene_annotation(module)["gene_names"]["g1"] == "g1"
        module.training_chromosomes = ["3"]
        with pytest.raises(PopsError, match="absent"):
            _validate_gene_annotation(module)
    assert calls == [True]


@pytest.mark.parametrize("rows,match", [
    ("g1\t1\t1\ng1\t2\t2\n", "duplicate"),
    ("g1\t1\tNA\n", "finite"),
    ("g1\t1\t-1\n", "non-negative"),
    ("g1\t1\tinf\n", "finite"),
])
def test_annotation_invalid_rows_fail(tmp_path, rows, match):
    path = tmp_path / "annotation.tsv"
    path.write_text("ID\tCHR\tTSS\n" + rows)
    with pytest.raises(ValueError, match=match):
        annotations.validate_gene_tss_annotation(
            path, delimiter="\t", id_column="ID", chromosome_column="CHR", tss_column="TSS",
            name_column="NAME", require_names=False, require_nonempty_identifiers=True, label="Annotation")


def test_gene_property_shared_scan_and_generated_universe_contract(tmp_path, monkeypatch):
    path = tmp_path / "properties.tsv"
    path.write_text("GENE Trait Average\ng1 0.1 0.2\ng2 0.3 0.4\ng3 NA 0.6\n")
    original = properties._inspect_magma_covariate_table
    calls = []
    monkeypatch.setattr(properties, "_inspect_magma_covariate_table", lambda *a, **k: calls.append(True) or original(*a, **k))
    policy = dict(minimum_genes=2, maximum_missing_fraction=0.5, missing_genes="drop")
    with InputValidationSession():
        first = properties.validate_magma_covariate_table(path, **policy)
        first["property_names"].clear()
        second = properties.validate_magma_covariate_table(path, **policy)
        assert second["property_names"] == ["Trait", "Average"]
        assert calls == [True]
        subset = properties.validate_magma_covariate_table(path, eligible_gene_ids={"g1", "g2"}, **policy)
        assert subset["overlapping_genes"] == 2
        assert subset["property_missingness"][0]["missing_fraction"] == 0
        with pytest.raises(ValueError, match="constant"):
            properties.validate_magma_covariate_table(path, eligible_gene_ids={"g1", "g3"}, **policy)
    assert len(calls) == 3


def test_native_gene_results_reuse_and_stricter_minimum(tmp_path, monkeypatch):
    path = tmp_path / "study.genes.raw"
    path.write_text("# VERSION = 110\ng1 1 0 10 1 1 1000\ng2 2 10 20 1 1 1000\n")
    original = properties._magma_gene_ids
    calls = []
    monkeypatch.setattr(properties, "_magma_gene_ids", lambda *a, **k: calls.append(True) or original(*a, **k))
    with InputValidationSession():
        _, genes = properties.validate_magma_gene_results(path, minimum_genes=1)
        genes.clear()
        assert len(properties.validate_magma_gene_results(path, minimum_genes=2)[1]) == 2
        with pytest.raises(ValueError, match="at least 3"):
            properties.validate_magma_gene_results(path, minimum_genes=3)
    assert calls == [True]


def test_identity_failure_uses_requested_error_type(tmp_path):
    class ConsumerError(ValueError):
        pass

    path = tmp_path / "resource"
    path.write_text("original")
    with InputValidationSession():
        validate_once((path,), "typed", lambda: 1, error_type=ConsumerError)
        path.write_text("changed")
        with pytest.raises(ConsumerError, match="changed"):
            validate_once((path,), "typed", lambda: 1, error_type=ConsumerError)


def test_explicit_vcf_evidence_cannot_bypass_session_file_identity(tmp_path, monkeypatch):
    import os

    path = tmp_path / "study.vcf.gz"
    path.write_bytes(b"fixture")
    Path(str(path) + ".tbi").write_bytes(b"index")
    header = "##fileformat=VCFv4.2\n##genome_build=GRCh37\n##contig=<ID=1>\n"
    monkeypatch.setattr(vcf, "run_checked_command", lambda command, *a, **k:
        header if command[1] == "view" else "1\n" if command[1] == "index" else "STUDY\n")
    contract = dict(genome_build_header="##genome_build={build}", supported_genome_builds=["GRCh37"])
    with InputValidationSession():
        cached = vcf.validate_indexed_vcf(path, "STUDY", "bcftools", **contract)
        metadata = path.stat()
        path.write_bytes(b"changed")  # Same byte count, with modification time restored.
        os.utime(path, ns=(metadata.st_atime_ns, metadata.st_mtime_ns))
        with pytest.raises(RuntimeError, match="changed"):
            vcf.validate_indexed_vcf(path, "STUDY", "bcftools", cached=cached, **contract)


@pytest.mark.parametrize("module_path,removed", [
    ("formatting/table.py", {"validate_harmonised_vcf_header", "validate_harmonised_indexed_vcf"}),
    ("filtering/sumstat_filter.py", {"_read_vcf_header", "_validate_vcf_contract"}),
    ("pops/service.py", {"_read_names", "_read_table", "_require_columns", "_matrix_values_are_finite",
                         "_inspect_feature_resources", "_inspect_gene_annotation"}),
    ("kpops/service.py", {"_read_table", "_require_columns"}),
    ("caldera/service.py", {"_read_table", "_require_columns"}),
    ("magmacovar/main.py", {"_magma_gene_ids", "validate_magma_gene_results", "validate_magma_covariate_table"}),
    ("ld_annotation/annot_ldblock.py", {"_contig_aliases", "_validate_contig_compatibility"}),
    ("gcta_gene/service.py", {"_inspect_fastbat_set_list"}),
    ("single_cell/methods/ldsc_celltype/runner.py", {"_read_ldcts", "_validate_references", "_reference_file_record"}),
])
def test_superseded_file_check_implementations_are_removed(module_path, removed):
    import ast

    source = Path(__file__).resolve().parents[1] / "src/postgwas/modules" / module_path
    tree = ast.parse(source.read_text(encoding="utf-8"))
    definitions = {node.name for node in ast.walk(tree)
                   if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))}
    assert definitions.isdisjoint(removed)


def test_binary_kernel_declared_byte_width_must_match_reader_dtype(tmp_path):
    path = tmp_path / "kernel.bin"
    np.ones((2, 2), dtype=np.float32).tofile(path)
    with pytest.raises(ValueError, match="dtype/byte-width"):
        matrices.validate_binary_kernel(path, dimension=2, dtype=np.float32,
            bytes_per_value=8, chunk_rows=1, relative_tolerance=0., absolute_tolerance=0., label="Kernel")


@pytest.mark.parametrize("family", ["vcf", "names", "features", "kernel", "annotation", "genes_raw", "properties"])
def test_shared_file_helpers_preserve_consumer_error_type_on_changed_input(tmp_path, monkeypatch, family):
    class ConsumerError(ValueError):
        pass

    path = tmp_path / "input"
    if family == "vcf":
        path.write_text("fixture")
        monkeypatch.setattr(vcf, "run_checked_command", lambda *a, **k: "header\n")
        validate = lambda: vcf.read_vcf_header(path, "bcftools", error_type=ConsumerError)
    elif family == "names":
        path.write_text("g1\ng2\n")
        validate = lambda: matrices.read_unique_names(path, "Names", error_type=ConsumerError)
    elif family == "features":
        rows, chunks = _matrix_bundle(tmp_path, [[1., 2.], [3., 4.]])
        path = chunks[0][0]
        validate = lambda: matrices.validate_feature_matrix_bundle(rows, chunks, error_type=ConsumerError)
    elif family == "kernel":
        np.eye(2, dtype=np.float32).tofile(path)
        validate = lambda: matrices.validate_binary_kernel(
            path, dimension=2, dtype=np.float32, bytes_per_value=4, chunk_rows=1,
            relative_tolerance=0., absolute_tolerance=0., label="Kernel", error_type=ConsumerError)
    elif family == "annotation":
        path.write_text("ID\tCHR\tTSS\ng1\t1\t1\n")
        validate = lambda: annotations.validate_gene_tss_annotation(
            path, delimiter="\t", id_column="ID", chromosome_column="CHR", tss_column="TSS",
            name_column="NAME", require_names=False, require_nonempty_identifiers=True,
            label="Annotation", error_type=ConsumerError)
    elif family == "genes_raw":
        path.write_text("# VERSION = 110\ng1 1 0 10 1 1 1000\n")
        validate = lambda: properties.validate_magma_gene_results(path, minimum_genes=1, error_type=ConsumerError)
    else:
        path.write_text("GENE Trait\ng1 0.1\ng2 0.2\n")
        validate = lambda: properties.validate_magma_covariate_table(
            path, minimum_genes=2, maximum_missing_fraction=0., missing_genes="drop", error_type=ConsumerError)
    with InputValidationSession():
        validate()
        path.write_bytes(path.read_bytes() + b"\n")
        with pytest.raises(ConsumerError, match="changed after pipeline preflight"):
            validate()
