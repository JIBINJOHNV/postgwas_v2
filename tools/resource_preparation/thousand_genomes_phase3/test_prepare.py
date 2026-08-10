from pathlib import Path
import importlib.util
import json
import shutil
import subprocess
import sys

import pytest
from pydantic import ValidationError


MODULE_PATH = Path(__file__).with_name("prepare.py")
SPEC = importlib.util.spec_from_file_location("prepare_1kg", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def base_config(tmp_path: Path) -> dict:
    return {
        "resource": {
            "name": "test",
            "genome_build": "GRCh37",
            "release": "test",
            "source_url": "https://example.org/release",
            "output_directory": tmp_path,
        },
        "panel": {
            "filename": "panel.tsv",
            "sample_column": "sample",
            "grouping_column": "super_pop",
            "groups": ["AFR", "EUR"],
        },
        "contigs": [{"label": "1", "filename": "chr1.vcf.gz"}],
        "output_sets": [
            {
                "name": "autosomes",
                "contigs": ["1"],
                "filename_template": "{group}.vcf.gz",
            }
        ],
        "processing": {
            "bcftools_executable": "bcftools",
            "curl_executable": "curl",
            "bcftools_threads": 1,
            "download_workers": 1,
            "download_retries": 1,
            "index_suffix": ".tbi",
            "output_index_type": "csi",
            "recalculate_info_tags": ["AC", "AN", "AF", "NS"],
            "remove_per_contig_outputs_after_success": True,
            "keep_downloads_after_success": True,
        },
        "directories": {
            "downloads": "downloads",
            "metadata": "metadata",
            "per_contig": "working/per_contig",
            "final": "final",
            "logs": "logs",
        },
        "naming": {
            "partial_suffix": ".partial",
            "group_sample_list_template": "{group}.samples.txt",
            "per_contig_vcf_template": "{group}.vcf.gz",
            "split_groups_filename": "groups.tsv",
            "split_prepared_vcf_template": "{group}.prepared.vcf.gz",
            "concat_input_list_filename": "inputs.txt",
            "log_filename": "preparation.jsonl",
            "manifest_filename": "resource_manifest.yaml",
        },
    }


def test_read_panel_uses_declared_super_population(tmp_path: Path) -> None:
    config = MODULE.WorkflowConfig.model_validate(base_config(tmp_path))
    panel = tmp_path / "panel.tsv"
    panel.write_text(
        "sample\tpop\tsuper_pop\tgender\n"
        "sample_a\tYRI\tAFR\tfemale\n"
        "sample_b\tGBR\tEUR\tmale\n",
        encoding="utf-8",
    )

    assert MODULE.read_panel(config, panel) == {
        "AFR": ["sample_a"],
        "EUR": ["sample_b"],
    }


def test_config_rejects_contig_in_multiple_output_sets(tmp_path: Path) -> None:
    raw = base_config(tmp_path)
    raw["output_sets"].append(
        {"name": "duplicate", "contigs": ["1"], "filename_template": "x.{group}.vcf.gz"}
    )

    with pytest.raises(ValidationError, match="occurs in output sets"):
        MODULE.WorkflowConfig.model_validate(raw)


def test_config_rejects_output_without_group_placeholder(tmp_path: Path) -> None:
    raw = base_config(tmp_path)
    raw["output_sets"][0]["filename_template"] = "fixed.vcf.gz"

    with pytest.raises(ValidationError, match="must contain"):
        MODULE.WorkflowConfig.model_validate(raw)


def test_panel_rejects_duplicate_samples(tmp_path: Path) -> None:
    config = MODULE.WorkflowConfig.model_validate(base_config(tmp_path))
    panel = tmp_path / "panel.tsv"
    panel.write_text(
        "sample\tsuper_pop\nshared\tAFR\nshared\tEUR\n",
        encoding="utf-8",
    )

    with pytest.raises(MODULE.PreparationError, match="duplicate sample"):
        MODULE.read_panel(config, panel)


def test_logger_records_unfinished_previous_run(tmp_path: Path) -> None:
    log_path = tmp_path / "preparation.jsonl"
    logger = MODULE.JsonLinesRunLogger(log_path)
    logger.write("run_started", stage="download")

    logger.record_unfinished_previous_run()

    records = [
        json.loads(line) for line in log_path.read_text(encoding="utf-8").splitlines()
    ]
    assert records[-1]["event"] == "prior_run_interruption_detected"


def test_source_validation_does_not_require_index_count_metadata(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = MODULE.WorkflowConfig.model_validate(base_config(tmp_path))
    paths = MODULE.layout(config)
    paths["downloads"].mkdir(parents=True)
    source = paths["downloads"] / "chr1.vcf.gz"
    source.write_bytes(b"source")
    Path(f"{source}.tbi").write_bytes(b"index")
    log_path = tmp_path / "source-validation.jsonl"
    logger = MODULE.JsonLinesRunLogger(log_path)
    monkeypatch.setattr(MODULE, "vcf_samples", lambda *_: ["sample_a"])
    monkeypatch.setattr(MODULE, "source_vcf_index_contigs", lambda *_: ["1"])
    monkeypatch.setattr(
        MODULE,
        "vcf_index_stats",
        lambda *_: (_ for _ in ()).throw(AssertionError("must not be called")),
    )

    MODULE.validate_sources(config, paths, "bcftools", logger)

    record = json.loads(log_path.read_text(encoding="utf-8").splitlines()[-1])
    assert record["records"] is None
    assert record["record_count_status"] == "unavailable_in_upstream_tbi"


@pytest.mark.skipif(
    shutil.which("bcftools") is None or shutil.which("bgzip") is None,
    reason="bcftools and bgzip are required for the integration test",
)
def test_split_and_single_contig_concat_preserve_population_genotypes(
    tmp_path: Path,
) -> None:
    raw = base_config(tmp_path)
    config = MODULE.WorkflowConfig.model_validate(raw)
    paths = MODULE.layout(config)
    paths["downloads"].mkdir(parents=True)
    source_text = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=1000>\n"
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n"
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Samples\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        "afr_1\tafr_2\teur_1\teur_2\n"
        "1\t10\t.\tA\tG\t.\tPASS\tAC=4;AN=8;AF=0.5;NS=4\tGT\t"
        "0/1\t1/1\t0/0\t0/1\n"
    )
    plain = tmp_path / "source.vcf"
    plain.write_text(source_text, encoding="utf-8")
    source = paths["downloads"] / "chr1.vcf.gz"
    with source.open("wb") as output:
        subprocess.run(["bgzip", "--stdout", str(plain)], check=True, stdout=output)
    subprocess.run(["bcftools", "index", "--tbi", str(source)], check=True)
    groups = {"AFR": ["afr_1", "afr_2"], "EUR": ["eur_1", "eur_2"]}
    logger = MODULE.JsonLinesRunLogger(tmp_path / "test.jsonl")

    metrics = MODULE.split_contig(
        config,
        config.contigs[0],
        groups,
        paths,
        "bcftools",
        logger,
    )
    final = MODULE.concat_outputs(
        config,
        config.output_sets[0],
        "AFR",
        paths,
        "bcftools",
        logger,
    )

    assert metrics["AFR"]["samples"] == 2
    assert metrics["AFR"]["contigs"] == ["1"]
    assert MODULE.vcf_samples("bcftools", final) == ["afr_1", "afr_2"]
    info = subprocess.run(
        ["bcftools", "query", "--format", "%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/NS\n", str(final)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout.strip()
    assert info == "3\t4\t0.75\t2"
