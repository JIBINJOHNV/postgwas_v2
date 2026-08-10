"""Regression tests for exact-build, exact-chromosome resource preflight."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.resource_preflight import (
    ResourcePreflightError,
    recheck_preflighted_resource_map,
    validate_harmonisation_resource_maps,
)
from postgwas.modules.harmonisation.service import (
    build_resource_map,
    harmonise_chromosomes,
    preflight_harmonisation_resources,
)


MAPPING = {
    "chr": "CHROM",
    "pos": "POS",
    "a1": "ALT",
    "a2": "REF",
    "delimiter": "tab",
}

VCF_CONFIG = {
    "external_frequency_columns": [
        "CHROM", "POS", "REF", "ALT",
        "INFO/AFR", "INFO/EAS", "INFO/EUR", "INFO/SAS",
    ],
}


def _write(path: Path, text: str = "data\n") -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return str(path)


def _resource_maps(tmp_path: Path) -> dict[str, dict[str, str]]:
    annotation = _write(
        tmp_path / "genes.gff3",
        "##gff-version 3\n1\tsource\tgene\t1\t10\t.\t+\t.\tID=gene1\n",
    )
    chain = _write(
        tmp_path / "source_to_target.chain",
        "chain 1 1 100 + 0 100 1 100 + 0 100 1\n100\n\n"
        "chain 1 2 100 + 0 100 2 100 + 0 100 2\n100\n",
    )
    external_eaf = _write(
        tmp_path / "external_eaf.tsv",
        "CHROM\tPOS\tREF\tALT\tEAF\n1\t10\tA\tG\t0.2\n",
    )
    external_info = _write(
        tmp_path / "external_info.tsv",
        "CHROM\tPOS\tREF\tALT\tINFO\n1\t10\tA\tG\t0.9\n",
    )
    maps = {}
    for chromosome in ("1", "2"):
        comparison = Path(_write(tmp_path / ("comparison_chr%s.vcf.gz" % chromosome)))
        _write(Path(str(comparison) + ".tbi"))
        dbsnp = Path(_write(tmp_path / ("dbsnp_chr%s.vcf.gz" % chromosome)))
        _write(Path(str(dbsnp) + ".csi"))
        source_fasta = Path(_write(tmp_path / ("source_chr%s.fa" % chromosome)))
        _write(
            Path(str(source_fasta) + ".fai"),
            "%s\t100\t3\t100\t101\n" % chromosome,
        )
        target_fasta = Path(_write(tmp_path / ("target_chr%s.fa" % chromosome)))
        _write(
            Path(str(target_fasta) + ".fai"),
            "%s\t100\t3\t100\t101\n" % chromosome,
        )
        default_eaf = _write(
            tmp_path / ("default_eaf_chr%s.tsv" % chromosome),
            "CHROM\tPOS\tREF\tALT\tEUR\n%s\t10\tA\tG\t0.2\n"
            % chromosome,
        )
        maps[chromosome] = {
            "user_eaf_file": external_eaf,
            "default_eaf_file": default_eaf,
            "default_comparison_af_file": str(comparison),
            "user_eaf_column": "EAF",
            "default_comparison_af_column": "EUR",
            "user_info_file": external_info,
            "user_info_column": "INFO",
            "genome_fasta_file": str(source_fasta),
            "dbsnp_file": str(dbsnp),
            "annot_path": annotation,
            "target_fasta": str(target_fasta),
            "chain_file": chain,
        }
    return maps


def _mock_bcftools(arguments, _purpose, **_kwargs):
    path = Path(arguments[-1])
    if arguments[1:3] == ["index", "--stats"]:
        chromosome = path.name.split("chr", 1)[1].split(".", 1)[0]
        return "%s\t100\t1\n" % chromosome
    if arguments[1:3] == ["view", "--header-only"]:
        return "\n".join(
            "##INFO=<ID=%s,Number=1,Type=Float,Description=frequency>" % tag
            for tag in ("AFR", "EAS", "EUR", "SAS")
        )
    raise AssertionError("unexpected command: %r" % (arguments,))


def test_resource_preflight_checks_exact_resources_once(monkeypatch, tmp_path):
    maps = _resource_maps(tmp_path)
    monkeypatch.setattr(
        "postgwas.modules.harmonisation.resource_preflight.run_checked_command",
        _mock_bcftools,
    )

    summary = validate_harmonisation_resource_maps(
        maps,
        bcftools="bcftools",
        vcf_config=VCF_CONFIG,
        default_eaf_colmap=MAPPING,
        external_eaf_colmap=MAPPING,
        external_info_colmap=MAPPING,
        policies=default_policies(),
        require_default_eaf=True,
    )

    assert summary["status"] == "passed"
    assert summary["chromosomes"] == ["1", "2"]
    assert summary["unique_resource_files_checked"] == 14
    assert summary["indexed_vcfs_checked"] == 4
    assert summary["frequency_vcf_headers_checked"] == 2
    assert summary["fasta_indexes_checked"] == 4
    assert summary["table_headers_checked"] == 4
    assert summary["structured_text_files_checked"] == 2
    assert summary["resources_by_chromosome"]["1"]["user_info_column"] == "INFO"


def test_preflight_path_resolution_does_not_repeat_existence_checks(
    monkeypatch, tmp_path,
):
    monkeypatch.setattr(
        "postgwas.modules.harmonisation.service.validate_path",
        lambda **_kwargs: pytest.fail("preflight repeated build_resource_map file checks"),
    )
    resources = build_resource_map(
        chromosome="1",
        grch_version="GRCh37",
        target_build="GRCh38",
        resource_folder=str(tmp_path),
        user_eaf_file=str(tmp_path / "external_{chromosome}.tsv"),
        default_comparison_af_file="1000G",
        resource_layout={
            "default_eaf": "{build}/default_{chromosome}.tsv",
            "comparison_af": "{build}/comparison_{chromosome}.vcf.gz",
            "dbsnp": "{build}/{source}_{chromosome}.vcf.gz",
            "fasta": "{build}/{build}_{chromosome}.fa",
            "annotation": "{build}/genes.gff3.gz",
            "chain": "{source_build}_to_{target_build}.chain",
        },
        user_info_file=None,
        dbsnp="dbSNP157",
        user_eaf_column="EAF",
        default_comparison_af_column="EUR",
        user_info_column=None,
        require_default_eaf=True,
        validate_files=False,
    )

    assert resources["user_eaf_file"].endswith("external_1.tsv")
    assert resources["default_eaf_file"].endswith("GRCh37/default_1.tsv")


def test_dataset_preflight_reports_all_missing_resources_together(tmp_path):
    layout = {
        "default_eaf": "{build}/default_{chromosome}.tsv",
        "comparison_af": "{build}/comparison_{chromosome}.vcf.gz",
        "dbsnp": "{build}/{source}_{chromosome}.vcf.gz",
        "fasta": "{build}/{build}_{chromosome}.fa",
        "annotation": "{build}/genes.gff3.gz",
        "chain": "{source_build}_to_{target_build}.chain",
    }

    with pytest.raises(ResourcePreflightError) as caught:
        preflight_harmonisation_resources(
            ["1", "2"],
            grch_version="GRCh37",
            target_build="GRCh38",
            resource_folder=str(tmp_path),
            user_eaf_file=None,
            default_comparison_af_file="1000G",
            resource_layout=layout,
            user_info_file=None,
            dbsnp="dbSNP157",
            user_eaf_column=None,
            default_comparison_af_column="EUR",
            user_info_column=None,
            default_eaf_colmap=MAPPING,
            external_eaf_colmap=MAPPING,
            external_info_colmap=MAPPING,
            executables={"bcftools": "bcftools"},
            vcf_config=VCF_CONFIG,
            policies=default_policies(),
            require_default_eaf=True,
        )

    labels = {issue["resource"] for issue in caught.value.issues}
    assert "default allele-frequency table" in labels
    assert "comparison allele-frequency VCF" in labels
    assert "dbSNP VCF" in labels
    assert "source-build FASTA" in labels
    assert "target-build FASTA" in labels
    assert "gene annotation GFF" in labels
    chain_issue = next(
        issue for issue in caught.value.issues
        if issue["resource"] == "liftover chain"
    )
    assert chain_issue["chromosomes"] == ["1", "2"]


def test_resource_preflight_groups_shared_failures_by_chromosome(
    monkeypatch, tmp_path,
):
    maps = _resource_maps(tmp_path)
    Path(maps["1"]["genome_fasta_file"] + ".fai").unlink()
    Path(maps["1"]["user_info_file"]).write_text(
        "CHROM\tPOS\tREF\tALT\n1\t10\tA\tG\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "postgwas.modules.harmonisation.resource_preflight.run_checked_command",
        _mock_bcftools,
    )

    with pytest.raises(ResourcePreflightError) as caught:
        validate_harmonisation_resource_maps(
            maps,
            bcftools="bcftools",
            vcf_config=VCF_CONFIG,
            default_eaf_colmap=MAPPING,
            external_eaf_colmap=MAPPING,
            external_info_colmap=MAPPING,
            policies=default_policies(),
            require_default_eaf=True,
        )

    shared = [
        issue for issue in caught.value.issues
        if issue["resource"] == "external INFO table"
    ]
    assert len(shared) == 1
    assert shared[0]["chromosomes"] == ["1", "2"]
    assert "INFO" in shared[0]["problem"]
    assert "source-build FASTA" in str(caught.value)
    assert "valid .fai index" in str(caught.value)


def test_worker_recheck_detects_a_resource_removed_after_preflight(tmp_path):
    resources = _resource_maps(tmp_path)["1"]
    Path(resources["chain_file"]).unlink()

    with pytest.raises(ResourcePreflightError, match="liftover chain"):
        recheck_preflighted_resource_map(
            "1", resources, require_default_eaf=True,
        )


def test_chain_target_label_can_differ_when_target_fasta_matches(
    monkeypatch, tmp_path,
):
    resources = _resource_maps(tmp_path)["1"]
    Path(resources["chain_file"]).write_text(
        "chain 1 1 100 + 0 100 chr1 100 + 0 100 1\n100\n",
        encoding="utf-8",
    )
    Path(resources["target_fasta"] + ".fai").write_text(
        "chr1\t100\t3\t100\t101\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "postgwas.modules.harmonisation.resource_preflight.run_checked_command",
        _mock_bcftools,
    )

    summary = validate_harmonisation_resource_maps(
        {"1": resources},
        bcftools="bcftools",
        vcf_config=VCF_CONFIG,
        default_eaf_colmap=MAPPING,
        external_eaf_colmap=MAPPING,
        external_info_colmap=MAPPING,
        policies=default_policies(),
        require_default_eaf=True,
    )

    assert summary["status"] == "passed"


def test_frequency_vcf_tags_must_be_numeric_allele_values(
    monkeypatch, tmp_path,
):
    resources = _resource_maps(tmp_path)["1"]

    def invalid_header(arguments, purpose, **kwargs):
        if arguments[1:3] == ["view", "--header-only"]:
            return "\n".join(
                "##INFO=<ID=%s,Number=A,Type=%s,Description=frequency>"
                % (tag, "String" if tag == "EUR" else "Float")
                for tag in ("AFR", "EAS", "EUR", "SAS")
            )
        return _mock_bcftools(arguments, purpose, **kwargs)

    monkeypatch.setattr(
        "postgwas.modules.harmonisation.resource_preflight.run_checked_command",
        invalid_header,
    )

    with pytest.raises(ResourcePreflightError, match="EUR.*Type=String"):
        validate_harmonisation_resource_maps(
            {"1": resources},
            bcftools="bcftools",
            vcf_config=VCF_CONFIG,
            default_eaf_colmap=MAPPING,
            external_eaf_colmap=MAPPING,
            external_info_colmap=MAPPING,
            policies=default_policies(),
            require_default_eaf=True,
        )


def test_exact_resource_preflight_precedes_partition_and_worker_launch():
    source = inspect.getsource(harmonise_chromosomes)

    assert source.index("preflight_harmonisation_resources(") < source.index(
        "write_chromosome_partitions("
    )
    assert source.index("write_chromosome_partitions(") < source.index(
        "_run_one_round("
    )
