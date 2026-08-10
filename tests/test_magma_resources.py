"""Regression tests for pinned MAGMA functional-mapping resource preparation."""

from copy import deepcopy
from pathlib import Path
import re
import zipfile

import pytest
import yaml

from postgwas.config import load_configuration, load_module_configuration
from postgwas.core.resource_preparation import sha256
from postgwas.resources import cli as resource_cli
from postgwas.resources.magma_functional_mapping import preparation as PREPARER

RESOURCE_PACKAGE = Path(PREPARER.__file__).parent
_generated_module_config = PREPARER._generated_module_config
_generated_run_config = PREPARER._generated_run_config
_read_gene_identifier_bridge = PREPARER._read_gene_identifier_bridge
_extract = PREPARER._extract


def _resource_configuration() -> dict:
    return yaml.safe_load(
        (RESOURCE_PACKAGE / "config.yaml").read_text(encoding="utf-8")
    )


def test_resource_archives_are_pinned_and_recursive_resources_are_selected():
    config = _resource_configuration()
    assert config["download_executable"] == "curl"
    for source in config["downloads"].values():
        assert source["url"].startswith("https://")
        assert re.fullmatch(r"[0-9a-f]{64}", source["sha256"])
        assert all(
            nested["archive_type"] in {"zip", "tar.gz"}
            for nested in source.get("nested", [])
        )

    assert "Codes/**/*" in config["downloads"]["h_magma_core"]["include"]
    assert "Input_Files/**/*" in config["downloads"]["h_magma_core"]["include"]
    assert (
        "Annotation_Files/**/*"
        in config["downloads"]["h_magma_protocol"]["include"]
    )
    assert (
        "Analysis/Rscripts/**/*"
        in config["downloads"]["chrom_magma"]["include"]
    )
    assert config["configuration"]["expected_mapping_definition_counts"] == {
        "positional": 1,
        "emagma": 47,
        "h_magma": 35,
        "n_magma": 4,
        "chrom_magma": 1,
    }
    protocol_group = config["configuration"]["annotation_groups"][2]
    assert len(protocol_group["context_aliases"]) == 28
    assert (
        protocol_group["context_aliases"]["X5628FC.transcript"]
        == "Dorsolateral prefrontal cortex"
    )


def test_archive_extraction_excludes_configured_packaging_artifacts(tmp_path):
    archive = tmp_path / "reference.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("reference/data.txt", "scientific data\n")
        handle.writestr("reference/.DS_Store", "metadata\n")
        handle.writestr("reference/__MACOSX/._data.txt", "metadata\n")

    destination = tmp_path / "extracted"
    _extract(
        archive,
        "zip",
        destination,
        _resource_configuration()["archive_member_exclude"],
    )

    assert (destination / "reference" / "data.txt").is_file()
    assert not (destination / "reference" / ".DS_Store").exists()
    assert not (destination / "reference" / "__MACOSX").exists()


def test_identifier_bridge_retains_only_unambiguous_symbol_links(tmp_path):
    exported = tmp_path / "ensembl_symbol.tsv"
    exported.write_text(
        "ensembl_gene_id\thgnc_symbol\n"
        "ENSG1.1\tGENE1\n"
        "ENSG2\tAMBIGUOUS\n"
        "ENSG3\tGENE3\n"
        "ENSG3\tOTHER\n",
        encoding="utf-8",
    )
    locations = tmp_path / "genes.loc"
    locations.write_text(
        "101 1 10 20 + GENE1\n"
        "201 1 30 40 + AMBIGUOUS\n"
        "202 1 50 60 + AMBIGUOUS\n"
        "303 1 70 80 + GENE3\n",
        encoding="utf-8",
    )
    specification = deepcopy(
        _resource_configuration()["identifier_harmonisation"]["emagma_networks"]
    )

    observed = _read_gene_identifier_bridge(exported, locations, specification)

    assert observed == [("ENSG1", "GENE1", "101")]


def test_generated_mapping_configuration_is_schema_validated(tmp_path):
    config = deepcopy(_resource_configuration()["configuration"])
    config["annotation_groups"] = [config["annotation_groups"][0]]
    config["expected_mapping_definition_counts"] = {
        "positional": 1,
        "emagma": 1,
        "n_magma": 1,
        "chrom_magma": 1,
    }
    config["n_magma"]["tissues"] = {
        "cortex": config["n_magma"]["tissues"]["cortex"]
    }
    annotation = tmp_path / "emagma" / "annotations" / "Brain_Cortex.genes.annot"
    annotation.parent.mkdir(parents=True)
    annotation.write_text("101 1:10:20 rs1\n", encoding="utf-8")
    network = tmp_path / "emagma" / "networks_entrez" / "Brain_Cortex.txt"
    network.parent.mkdir(parents=True)
    network.write_text("module 101\n", encoding="utf-8")
    for relative in [
        config["gene_location_file"],
        config["n_magma"]["gene_location_file"],
        *config["n_magma"]["tissues"]["cortex"]["components"],
        config["chrom_magma"]["regulatory_element_location_file"],
        config["chrom_magma"]["element_to_gene_file"],
    ]:
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("resource\n", encoding="utf-8")

    generated = _generated_module_config(tmp_path, tmp_path, config)
    generated_path = tmp_path / "magma.yaml"
    generated_path.write_text(
        yaml.safe_dump(generated, sort_keys=False), encoding="utf-8",
    )
    module = load_module_configuration("magma", generated_path)

    assert module.mapping.definitions["emagma_brain_cortex"].gene_id_type == "entrez"
    assert (
        module.mapping.definitions["emagma_brain_cortex"]
        .minimum_gene_id_overlap_fraction
        == pytest.approx(0.01)
    )
    assert module.mapping.definitions["n_magma_cortex"].gene_id_type == "symbol"
    assert (
        module.mapping.definitions["n_magma_cortex"]
        .annotation_window_upstream_kb
        == 0
    )
    assert (
        module.mapping.definitions["n_magma_cortex"]
        .annotation_window_downstream_kb
        == 0
    )
    assert (
        module.mapping.definitions["chrom_magma_ovarian_h3k27ac"]
        .result_statistic_type
        == "minimum_regulatory_element_p_value"
    )
    assert (
        module.mapping.definitions["chrom_magma_ovarian_h3k27ac"].gene_id_type
        == "mixed"
    )

    generated_run = _generated_run_config(tmp_path, tmp_path, config)
    run_path = tmp_path / "magma_pipeline.yaml"
    run_path.write_text(
        yaml.safe_dump(generated_run, sort_keys=False), encoding="utf-8",
    )
    pipeline = load_configuration(run_path)
    assert pipeline.modules.magma.mapping.primary == "positional"
    assert "emagma_brain_cortex" in pipeline.modules.magma.mapping.definitions


def test_metadata_refresh_converts_module_yaml_to_pipeline_yaml(
    tmp_path, monkeypatch,
):
    config = _resource_configuration()
    output = tmp_path / "functional_mapping"
    generated_path = output / config["generated_config_file"]
    generated_path.parent.mkdir(parents=True)
    generated_path.write_text("enabled: false\n", encoding="utf-8")
    resource = output / "reference.txt"
    resource.write_text("verified scientific resource\n", encoding="utf-8")
    manifest_path = output / config["manifest_file"]
    manifest_path.write_text(
        yaml.safe_dump({
            "schema_version": 1,
            "files": [
                {
                    "path": config["generated_config_file"],
                    "bytes": generated_path.stat().st_size,
                    "sha256": sha256(generated_path),
                },
                {
                    "path": "reference.txt",
                    "bytes": resource.stat().st_size,
                    "sha256": sha256(resource),
                },
            ],
        }, sort_keys=False),
        encoding="utf-8",
    )
    module_values = load_module_configuration("magma").model_dump(mode="json")
    monkeypatch.setattr(
        PREPARER,
        "_generated_run_config",
        lambda *_: {"modules": {"magma": module_values}},
    )

    PREPARER.refresh_metadata(RESOURCE_PACKAGE / "config.yaml", output)

    refreshed = load_configuration(generated_path)
    assert refreshed.modules.magma.mapping.primary == "positional"
    assert resource.read_text(encoding="utf-8") == "verified scientific resource\n"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    config_record = next(
        record for record in manifest["files"]
        if record["path"] == config["generated_config_file"]
    )
    assert config_record["sha256"] == sha256(generated_path)
    assert manifest["preparation_file_sha256"]["preparation.py"] == sha256(
        RESOURCE_PACKAGE / "preparation.py"
    )


def test_packaged_resource_paths_are_user_relative():
    config = _resource_configuration()

    assert config["default_output_directory"].startswith("~/")
    assert config["default_cache_directory"].startswith("~/")
    assert "/Users/" not in config["default_output_directory"]


def test_installed_resource_cli_dispatches_to_packaged_preparer(
    tmp_path, monkeypatch, capsys,
):
    output = tmp_path / "functional_mapping"
    generated = output / "configs" / "magma_functional_mapping.yaml"
    observed = {}

    def fake_run(action, **values):
        observed.update({"action": action, **values})
        return output, generated

    monkeypatch.setattr(resource_cli, "run", fake_run)

    status = resource_cli.main([
        "refresh", "magma", "--output-directory", str(output),
    ])

    assert status == 0
    assert observed["action"] == "refresh"
    assert observed["output_directory"] == output
    assert observed["cache_directory"] is None
    assert observed["config_path"] == RESOURCE_PACKAGE / "config.yaml"
    assert observed["configuration"]["schema_version"] == 1
    terminal = capsys.readouterr().out
    assert str(output) in terminal
    assert str(generated) in terminal
