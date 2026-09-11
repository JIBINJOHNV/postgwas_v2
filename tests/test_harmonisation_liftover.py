"""Scientific and command contracts for allele-aware build liftover."""

import shutil
import subprocess

import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.policies import PolicyError, default_policies
from postgwas.modules.harmonisation.vcf_processing import (
    _parse_liftover_audit,
    _resolve_liftover_tag_arguments,
)


ANNOTATED_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Study ALT frequency\">\n"
    "##INFO=<ID=AFR,Number=A,Type=Float,Description=\"AFR ALT frequency\">\n"
    "##INFO=<ID=EAS,Number=A,Type=Float,Description=\"EAS ALT frequency\">\n"
    "##INFO=<ID=EUR,Number=A,Type=Float,Description=\"EUR ALT frequency\">\n"
    "##INFO=<ID=SAS,Number=A,Type=Float,Description=\"SAS ALT frequency\">\n"
    "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Study ALT frequency\">\n"
    "##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"ALT effect\">\n"
    "##FORMAT=<ID=EZ,Number=A,Type=Float,Description=\"ALT Z score\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tstudy\n"
)


def _vcf_config():
    return load_configuration().modules.harmonisation.vcf_processing.model_dump()


def test_liftover_fields_are_detected_from_canonical_roles_and_annotations():
    arguments, resolved = _resolve_liftover_tag_arguments(
        ANNOTATED_HEADER, _vcf_config(),
    )

    assert resolved == {
        "allele_frequency": [
            "INFO/AF", "FORMAT/AF", "INFO/AFR", "INFO/EAS", "INFO/EUR",
            "INFO/SAS",
        ],
        "signed_effect": ["FORMAT/ES", "FORMAT/EZ"],
    }
    assert arguments == [
        "--af-tags",
        "INFO/AF,FMT/AF,INFO/AFR,INFO/EAS,INFO/EUR,INFO/SAS",
        "--es-tags",
        "FMT/ES,FMT/EZ",
    ]
    assert all("AP1" not in value and "AP2" not in value and "ED" not in value
               for value in arguments)


@pytest.mark.parametrize(
    ("old", "new", "problem"),
    [
        (
            "##INFO=<ID=EUR,Number=A,Type=Float",
            "##INFO=<ID=EUR,Number=1,Type=Float",
            "INFO/EUR has Number=1",
        ),
        (
            "##FORMAT=<ID=EZ,Number=A,Type=Float",
            "##FORMAT=<ID=EZ,Number=A,Type=String",
            "FORMAT/EZ has Type=String",
        ),
    ],
)
def test_liftover_rejects_incompatible_allele_dependent_header_fields(
    old, new, problem,
):
    with pytest.raises(ValueError, match=problem):
        _resolve_liftover_tag_arguments(
            ANNOTATED_HEADER.replace(old, new), _vcf_config(),
        )


def test_liftover_audit_separates_swaps_rejects_splits_and_duplicates():
    audit = _parse_liftover_audit(
        "Lines   total/swapped/reference added/rejected:\t100/7/3/5\n"
        "Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:"
        "\t85/3/0/4/0/0/0\n"
        "Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:"
        "\t88/0/0/0/0/2/0\n"
    )

    assert audit["plugin"] == {
        "total": 100, "swapped": 7, "reference_added": 3, "rejected": 5,
    }
    assert audit["normalization"]["split"] == 3
    assert audit["deduplication"]["duplicate_removed"] == 2


def test_liftover_policy_surface_contains_only_exclude_and_keep():
    policies = default_policies()

    assert policies.get("vcf.liftover_swap") == "exclude"
    assert "vcf.liftover_update_tag_args" not in policies
    assert policies.with_overrides({"vcf.liftover_swap": "keep"}).get(
        "vcf.liftover_swap"
    ) == "keep"
    with pytest.raises(PolicyError, match="vcf.liftover_swap"):
        policies.with_overrides({"vcf.liftover_swap": "update_tags"})


def test_liftover_plugin_updates_configured_af_es_and_ez_fields(tmp_path):
    """Prove the configured roles against an installed liftover plugin."""
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("bcftools is required for the liftover integration contract")
    plugin_list = subprocess.run(
        [bcftools, "plugin", "-l"],
        check=False,
        capture_output=True,
        text=True,
    )
    if plugin_list.returncode or "liftover" not in plugin_list.stdout.splitlines():
        pytest.skip("the bcftools liftover plugin is required")

    source_fasta = tmp_path / "source.fa"
    target_fasta = tmp_path / "target.fa"
    chain = tmp_path / "source_to_target.chain"
    source_vcf = tmp_path / "source.vcf"
    target_vcf = tmp_path / "target.vcf"
    source_fasta.write_text(">1\nAAAAAAAAAAAAAAAAAAAA\n", encoding="utf-8")
    target_fasta.write_text(">1\nAAAAGAAAAAAAAAAAAAAA\n", encoding="utf-8")
    chain.write_text(
        "chain 1 1 20 + 0 20 1 20 + 0 20 1\n20\n\n",
        encoding="utf-8",
    )
    source_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=20>\n"
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Study ALT frequency\">\n"
        "##INFO=<ID=AFR,Number=A,Type=Float,Description=\"AFR ALT frequency\">\n"
        "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Study ALT frequency\">\n"
        "##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"ALT effect\">\n"
        "##FORMAT=<ID=EZ,Number=A,Type=Float,Description=\"ALT Z score\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tstudy\n"
        "1\t5\trs1\tA\tG\t.\tPASS\tAF=0.2;AFR=0.25\tAF:ES:EZ\t0.2:0.4:2\n",
        encoding="utf-8",
    )

    completed = subprocess.run(
        [
            bcftools, "+liftover", str(source_vcf), "--no-version",
            "--output-type", "v", "--output", str(target_vcf), "--",
            "--src-fasta-ref", str(source_fasta),
            "--fasta-ref", str(target_fasta),
            "--chain", str(chain),
            "--af-tags", "INFO/AF,FMT/AF,INFO/AFR",
            "--es-tags", "FMT/ES,FMT/EZ",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr

    query = subprocess.run(
        [
            bcftools, "query", "--format",
            "%REF\t%ALT\t%INFO/SWAP\t%INFO/AF\t%INFO/AFR"
            "\t[%AF\t%ES\t%EZ]\n",
            str(target_vcf),
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.rstrip("\n").split("\t")
    assert query[:3] == ["G", "A", "1"]
    assert [float(value) for value in query[3:]] == pytest.approx(
        [0.8, 0.75, 0.8, -0.4, -2.0]
    )
