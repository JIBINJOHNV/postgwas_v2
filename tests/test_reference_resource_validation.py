"""Shared static resource checks retain scope, syntax and run-local identities."""

from argparse import Namespace
import gzip
from pathlib import Path
from types import SimpleNamespace

from pydantic import BaseModel
import pytest

from postgwas.config import load_configuration, load_module_configuration
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.io.tables import read_pandas_table
from postgwas.core import reference_resources as readers
from postgwas.modules.caldera import service as caldera
from postgwas.modules.flames import service as flames
from postgwas.modules.ld_clumping import service as ld_clump
from postgwas.modules.mixer.results import validate_gsa_go_file
from preflight_support import pipeline_input_vcf_evidence


class Manifest(BaseModel):
    names: list[str]


def _file(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def test_inventory_reports_all_missing_files_without_claiming_content(tmp_path):
    valid = _file(tmp_path / "valid.bin", "opaque")
    empty = _file(tmp_path / "empty.bin", "")
    absent = [tmp_path / ("missing%d" % index) for index in range(7)]
    with InputValidationSession() as session, session.scope("resources"):
        with pytest.raises(ValueError) as failure:
            readers.require_file_inventory(
                (valid, empty, *absent), "Reference",
                missing_message="Reference files missing or empty",
            )
    assert all(str(path) in str(failure.value) for path in (empty, *absent))
    assert len([record for record in session.records if record.status == "failed"]) == 8
    assert all("header" not in " ".join(record.checks) for record in session.records)


@pytest.mark.parametrize("delimiter, engine, text", (
    (None, "python", "gene,score\ngene1,0.1\n"),
    ("\t", "c", "gene\tscore\ngene1\t0.1\n"),
    ("\t", None, "gene\tscore\ngene1\t0.1\n"),
))
def test_shared_pandas_reader_preserves_explicit_and_detected_delimiters(
    tmp_path, delimiter, engine, text,
):
    path = _file(tmp_path / "table.txt", text)
    result = read_pandas_table(path, delimiter, "Table", engine=engine, error_type=caldera.CalderaError)
    assert result.to_dict("records") == [{"gene": "gene1", "score": 0.1}]
    path.write_text(text.splitlines()[0] + "\n")
    with pytest.raises(caldera.CalderaError, match="contains no data rows"):
        read_pandas_table(path, delimiter, "Table", engine=engine, error_type=caldera.CalderaError)


def test_inventory_retains_order_and_rechecks_removed_file(tmp_path):
    first, second = (_file(tmp_path / name, "opaque") for name in ("b", "a"))
    with InputValidationSession():
        assert readers.require_file_inventory(
            (first, second), "Panel", missing_message="Missing panel",
        ) == (first, second)
        second.unlink()
        with pytest.raises(ValueError, match="Missing panel"):
            readers.require_file_inventory((first, second), "Panel", missing_message="Missing panel")


def test_yaml_schema_reader_runs_once_and_returns_detached_model(tmp_path, monkeypatch):
    path = _file(tmp_path / "manifest.yaml", "names: [gene1, gene2]\n")
    original = readers.yaml.safe_load
    calls = []
    monkeypatch.setattr(readers.yaml, "safe_load", lambda value: calls.append(value) or original(value))
    with InputValidationSession() as session:
        with session.scope("first"):
            result = readers.read_yaml_manifest(path, Manifest, "Manifest")
        result.names.append("not in file")
        with session.scope("second"):
            assert readers.read_yaml_manifest(path, Manifest, "Manifest").names == ["gene1", "gene2"]
    assert len(calls) == 1
    record = next(record for record in session.records if "YAML mapping" in record.checks)
    assert record.consumers == ("first", "second")
    assert "contents are not established" in record.message
    readers.read_yaml_manifest(path, Manifest, "Manifest")
    readers.read_yaml_manifest(path, Manifest, "Manifest")
    assert len(calls) == 3


@pytest.mark.parametrize("text, message", (
    ("- item\n", "must be a YAML mapping"),
    ("names: [\n", "Cannot read"),
    ("names: null\n", "Invalid Manifest"),
    ("# only comment\n", "Invalid Manifest"),
))
def test_yaml_manifest_failures_have_precise_file_evidence(tmp_path, text, message):
    path = _file(tmp_path / "invalid.yaml", text)
    with InputValidationSession() as session:
        with pytest.raises(ValueError, match=message):
            readers.read_yaml_manifest(path, Manifest, "Manifest")
    assert any(record.status == "failed" and record.path == str(path) for record in session.records)


def test_manifest_changed_file_and_distinct_schema_cannot_use_old_result(tmp_path):
    path = _file(tmp_path / "manifest.yaml", "names: [gene1]\n")

    class OtherManifest(BaseModel):
        count: int

    with InputValidationSession():
        readers.read_yaml_manifest(path, Manifest, "Manifest")
        with pytest.raises(ValueError, match="count"):
            readers.read_yaml_manifest(path, OtherManifest, "Manifest")
        path.write_text("names: [different]\n")
        with pytest.raises(ValueError, match="changed"):
            readers.read_yaml_manifest(path, Manifest, "Manifest")


def test_ld_clump_manifest_parser_is_shared_but_compatibility_still_applies(tmp_path, monkeypatch):
    module = load_configuration().modules.ld_clumping
    reference = module.reference
    document = {
        "format_version": reference.format_version,
        "genome_build": module.genome_build.value,
        "populations": [module.population.value],
        "orientation": reference.orientation,
        "file_pattern": reference.file_pattern,
        "reverse_file_pattern": reference.reverse_file_pattern,
        "variant_inventory_pattern": reference.variant_inventory_pattern,
        "columns": reference.columns,
        "variant_inventory_columns": reference.variant_inventory_columns,
        "window_kb": module.window_kb,
        "minimum_r2": min(module.clump_r2, module.lead_r2),
        "minimum_maf": module.minimum_reference_maf,
        "allele_order_preserved": True,
        "plink_version": "PLINK v1.90",
    }
    _file(tmp_path / reference.manifest_filename, readers.yaml.safe_dump(document))
    original = readers.yaml.safe_load
    calls = []
    monkeypatch.setattr(readers.yaml, "safe_load", lambda value: calls.append(value) or original(value))
    with InputValidationSession():
        ld_clump._validate_reference_contract(module, tmp_path)
        ld_clump._validate_reference_contract(module, tmp_path)
        changed = module.model_copy(update={"window_kb": module.window_kb + 1})
        with pytest.raises(ld_clump.LDClumpingError, match="scientifically incompatible"):
            ld_clump._validate_reference_contract(changed, tmp_path)
    assert len(calls) == 1


@pytest.mark.parametrize("compressed", (False, True))
def test_ontology_header_reused_with_exact_scope(tmp_path, monkeypatch, compressed):
    path = tmp_path / ("ontology.tsv.gz" if compressed else "ontology.tsv")
    # The existing contract intentionally does not inspect later values/widths.
    text = "GO\tGENE\nset1\tgene1\nlater\trow\thas extra fields\n"
    if compressed:
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    settings = SimpleNamespace(go_file_delimiter="\t", go_file_required_columns=["GO", "GENE"])
    original = readers.open_text
    calls = []
    monkeypatch.setattr(readers, "open_text", lambda value: calls.append(value) or original(value))
    with InputValidationSession() as session:
        validate_gsa_go_file(path, settings)
        validate_gsa_go_file(path, settings)
        assert len(calls) == 1
        settings.go_file_required_columns = ["GO", "ABSENT"]
        with pytest.raises(RuntimeError, match="ABSENT"):
            validate_gsa_go_file(path, settings)
    record = next(record for record in session.records if "columns" in record.metrics)
    assert record.checks == ("configured header columns", "at least one data record")
    assert "remaining records are not validated" in record.message


@pytest.mark.parametrize("text, message", (("GO\tGENE\n", "contains no records"), ("GO\nset1\n", "GENE")))
def test_ontology_header_rejects_missing_columns_or_no_records(tmp_path, text, message):
    path = _file(tmp_path / "ontology.tsv", text)
    with pytest.raises(RuntimeError, match=message):
        validate_gsa_go_file(path, SimpleNamespace(
            go_file_delimiter="\t", go_file_required_columns=["GO", "GENE"],
        ))


def test_line_name_reader_preserves_whole_lines_and_reuses_content(tmp_path, monkeypatch):
    path = _file(tmp_path / "names.txt", " feature one \n\n#literal_name\nfeature_two\n")
    original = Path.read_text
    calls = []

    def read(source, *args, **kwargs):
        if source == path:
            calls.append(source)
        return original(source, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", read)
    with InputValidationSession():
        first = flames._feature_names(path)
        first.append("not in file")
        assert flames._feature_names(path) == ["feature one", "#literal_name", "feature_two"]
    assert len(calls) == 1
    flames._feature_names(path)
    flames._feature_names(path)
    assert len(calls) == 3


@pytest.mark.parametrize("text", (" \n\t\n", "feature\n feature \n"))
def test_line_name_reader_rejects_empty_or_duplicate_names(tmp_path, text):
    with pytest.raises(flames.FlamesError, match="unique non-empty names"):
        flames._feature_names(_file(tmp_path / "names.txt", text))


def test_line_name_reader_rejects_changed_file_within_run(tmp_path):
    path = _file(tmp_path / "names.txt", "feature1\n")
    with InputValidationSession():
        flames._feature_names(path)
        path.write_text("feature2\n")
        with pytest.raises(flames.FlamesError, match="changed"):
            flames._feature_names(path)


def test_bundle_inventory_rechecks_new_missing_and_empty_members(tmp_path):
    first = _file(tmp_path / "directory" / "first", "opaque")
    contract = dict(label="Annotation resource", bundle="Bundle")
    with InputValidationSession():
        _, inventory = readers.validate_resource_directories(tmp_path, ["directory"], **contract)
        assert list(inventory.values()) == [first]
        second = _file(first.parent / "second", "also opaque")
        _, inventory = readers.validate_resource_directories(tmp_path, ["directory"], **contract)
        assert set(inventory.values()) == {first, second}
        second.write_text("")
        with pytest.raises(ValueError, match="does not exist or is empty"):
            readers.validate_resource_directories(tmp_path, ["directory"], **contract)
        first.unlink()
        second.unlink()
        with pytest.raises(ValueError, match="without any files"):
            readers.validate_resource_directories(tmp_path, ["directory"], **contract)
        first.parent.rmdir()
        with pytest.raises(ValueError, match="missing configured directories"):
            readers.validate_resource_directories(tmp_path, ["directory"], **contract)


@pytest.mark.parametrize("consumer", ("flames", "caldera", "bundle"))
def test_resource_containment_rejects_symlink_escape(tmp_path, consumer):
    root = tmp_path / "root"
    root.mkdir()
    outside = _file(tmp_path / "outside", "opaque")
    (root / "escaped").symlink_to(outside)
    with pytest.raises(ValueError if consumer == "bundle" else RuntimeError, match="outside|leaves"):
        if consumer == "flames":
            flames._module_resource(root, "escaped", "Resource")
        elif consumer == "caldera":
            caldera._repository_file(root, "escaped", "Resource")
        else:
            readers.validate_resource_directories(root, ["."], label="Annotation resource", bundle="Bundle")


def test_flames_preflight_no_longer_sets_configuration_only_namespace_cache(monkeypatch):
    monkeypatch.setattr(flames, "_validate_upstream_resources", lambda configuration: {"checked": True})
    args = Namespace(annotation_resource_directory="resources")
    result = flames.preflight_flames_pipeline(args, preflight_evidence=pipeline_input_vcf_evidence())
    assert not hasattr(args, "_flames_resource_preflight")
    assert result.resources[1] == {"checked": True}
    assert result.deferred_checks


def test_flames_resource_wrapper_rechecks_bundle_but_reuses_feature_parse(tmp_path, monkeypatch):
    module = load_module_configuration("flames")
    annotation = tmp_path / "annotations"
    for relative in module.upstream.required_annotation_directories:
        _file(annotation / relative / "resource.dat", "opaque\n")
    for pattern in module.upstream.required_annotation_file_patterns:
        _file(annotation / pattern.format(genome_build=module.genome_build.value.upper()), "opaque\n")
    model = tmp_path / "model"
    _file(model / module.upstream.model_file, "opaque model")
    names = _file(model / module.upstream.feature_file, "feature1\nfeature2\n")
    module = module.model_copy(update={"model_directory": str(model), "annotation_resource_directory": str(annotation)})
    config = SimpleNamespace(modules=SimpleNamespace(flames=module))
    probes = []
    monkeypatch.setattr(flames, "_validate_runtime", lambda *args, **kwargs: probes.append(True) or "python")
    with InputValidationSession():
        result = flames._validate_upstream_resources(config)
        new_file = _file(annotation / module.upstream.required_annotation_directories[0] / "new.dat", "new")
        repeated = flames._validate_upstream_resources(config)
        assert len(repeated["annotation_inventory"]) == len(result["annotation_inventory"]) + 1
        new_file.unlink()
        names.write_text("changed feature\n")
        with pytest.raises(flames.FlamesError, match="changed"):
            flames._validate_upstream_resources(config)
    assert len(probes) == 3
