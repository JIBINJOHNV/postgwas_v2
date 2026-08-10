import bz2
import gzip
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


MODULE_PATH = Path(__file__).with_name("prepare.py")
SPEC = importlib.util.spec_from_file_location("prepare_dbsnp157", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def policy() -> object:
    return MODULE.ChromosomePolicyConfig(
        assembly_report_sequence_column="Sequence-Name",
        assembly_report_role_column="Sequence-Role",
        assembly_report_accession_column="RefSeq-Accn",
        assembly_report_length_column="Sequence-Length",
        included_sequence_role="assembled-molecule",
        expected_primary_labels=["1", "MT"],
    )


def test_parse_assembly_report_selects_declared_role(tmp_path: Path) -> None:
    report = tmp_path / "assembly_report.txt"
    report.write_text(
        "# Sequence-Name\tSequence-Role\tRefSeq-Accn\tSequence-Length\n"
        "1\tassembled-molecule\tNC_000001.10\t100\n"
        "PATCH\tfix-patch\tNW_000001.1\t10\n"
        "MT\tassembled-molecule\tNC_012920.1\t20\n",
        encoding="utf-8",
    )

    assert MODULE.parse_assembly_report(report, policy()) == [
        ("NC_000001.10", "1", 100),
        ("NC_012920.1", "MT", 20),
    ]


def test_merged_table_includes_historical_members(tmp_path: Path) -> None:
    source = tmp_path / "merged.json.bz2"
    record = {
        "refsnp_id": "332",
        "dbsnp1_merges": [
            {"merged_rsid": "33969899", "revision": "127", "merge_date": "date1"}
        ],
        "merged_snapshot_data": {
            "proxy_build_id": "133",
            "proxy_time": "date2",
            "merged_into": ["121909001"],
        },
    }
    with bz2.open(source, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")
    destination = tmp_path / "merged.tsv.gz"

    metrics = MODULE.build_merged_alias_table(source, destination, 6)

    with gzip.open(destination, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    assert metrics == {
        "path": str(destination),
        "source_records": 1,
        "output_rows": 2,
        "records_without_current_target": 0,
        "rows_without_current_target": 0,
    }
    assert lines[1].startswith("rs332\trs121909001\tmerged_record")
    assert lines[2].startswith("rs33969899\trs121909001\thistorical_merge_member")


def test_merged_table_retains_ids_without_a_current_target(tmp_path: Path) -> None:
    source = tmp_path / "merged.json.bz2"
    record = {
        "refsnp_id": "3044362",
        "dbsnp1_merges": [{"merged_rsid": "10608448", "revision": "137"}],
        "merged_snapshot_data": {
            "proxy_build_id": "152",
            "proxy_time": "date1",
            "merged_into": [],
        },
    }
    with bz2.open(source, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")
    destination = tmp_path / "merged.tsv.gz"

    metrics = MODULE.build_merged_alias_table(source, destination, 6)

    with gzip.open(destination, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    assert metrics["records_without_current_target"] == 1
    assert metrics["rows_without_current_target"] == 2
    assert lines[1].startswith("rs3044362\t\tmerged_without_current_target")
    assert lines[2].startswith(
        "rs10608448\t\thistorical_member_without_current_target"
    )


def test_withdrawn_table_propagates_status_to_historical_members(
    tmp_path: Path,
) -> None:
    source = tmp_path / "withdrawn.json.bz2"
    record = {
        "refsnp_id": "386",
        "last_update_build_id": "152",
        "dbsnp1_merges": [{"merged_rsid": "2387453"}],
        "withdrawn_snapshot_data": {"withdrawn_time": "2015-06-16T14:54Z"},
    }
    with bz2.open(source, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")
    destination = tmp_path / "withdrawn.tsv.gz"

    metrics = MODULE.build_withdrawn_status_table(source, destination, 6)

    with gzip.open(destination, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    assert metrics["output_rows"] == 2
    assert lines[1].startswith("rs386\trs386\twithdrawn_record")
    assert lines[2].startswith("rs2387453\trs386\thistorical_merge_member")


@pytest.mark.skipif(
    shutil.which("bcftools") is None or shutil.which("bgzip") is None,
    reason="bcftools and bgzip are required for the integration test",
)
def test_vcf_conversion_renames_primary_and_excludes_non_primary(
    tmp_path: Path,
) -> None:
    plain = tmp_path / "source.vcf"
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=NC_000001.10,length=100>\n"
        "##contig=<ID=NW_000001.1,length=10>\n"
        "##INFO=<ID=RS,Number=1,Type=Integer,Description=\"dbSNP ID\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "NC_000001.10\t1\trs1;rs2\tA\tG,T\t.\t.\tRS=1\n"
        "NW_000001.1\t2\trs2\tC\tT\t.\t.\tRS=2\n",
        encoding="utf-8",
    )
    source = tmp_path / "source.vcf.gz"
    with source.open("wb") as output:
        subprocess.run(["bgzip", "--stdout", str(plain)], check=True, stdout=output)
    subprocess.run(["bcftools", "index", "--tbi", str(source)], check=True)
    mapping = tmp_path / "map.tsv"
    mapping.write_text("NC_000001.10\t1\n", encoding="utf-8")
    targets = tmp_path / "targets.txt"
    targets.write_text("1\t1\t100\n", encoding="utf-8")
    destination = tmp_path / "standard.vcf.gz"

    metrics = MODULE.convert_vcf(
        source=source,
        destination=destination,
        chromosome_map=mapping,
        primary_targets=targets,
        expected_contigs=["1"],
        bcftools="bcftools",
        threads=1,
        index_type="csi",
    )

    assert metrics["records"] == 2
    record = subprocess.run(
        ["bcftools", "query", "--format", "%CHROM\t%ID\n", str(destination)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout.strip().splitlines()
    assert record == ["1\trs1;rs2", "1\trs1;rs2"]
    info = subprocess.run(
        ["bcftools", "query", "--format", "%INFO\n", str(destination)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout.strip()
    assert info.splitlines() == [".", "."]
    header = subprocess.run(
        ["bcftools", "view", "--header-only", str(destination)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout
    assert "##INFO=<" not in header
