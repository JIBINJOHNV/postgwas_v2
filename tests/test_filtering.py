"""Filtering configuration, scientific-contract, and execution regressions."""

from argparse import Namespace
from pathlib import Path
import shutil
import subprocess

import pytest

from postgwas.config import load_module_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.modules.filtering.cli import build_parser
from postgwas.modules.filtering.service import (
    resolve_filtering_configuration,
    run_sumstat_filter_direct,
)
from postgwas.modules.filtering.sumstat_filter import (
    _match_contig_naming,
    _validate_vcf_contract,
)


BCFTOOLS = shutil.which("bcftools")
TABIX = shutil.which("tabix")


def _write_filtering_vcf(path: Path) -> None:
    path.write_text(
        """##fileformat=VCFv4.2
##reference=GRCh37
##contig=<ID=6,length=171115067>
##INFO=<ID=AF,Number=A,Type=Float,Description="Study allele frequency">
##INFO=<ID=EUR,Number=A,Type=Float,Description="Reference allele frequency">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Study allele frequency">
##FORMAT=<ID=SI,Number=1,Type=Float,Description="Imputation quality">
##FORMAT=<ID=LP,Number=1,Type=Float,Description="Negative log10 P-value">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY
6\t100\trs1\tA\tC\t.\tPASS\tAF=0.10;EUR=0.12\tAF:SI:LP\t0.10:0.90:8
6\t200\trs2\tA\tT\t.\tPASS\tAF=0.50;EUR=0.48\tAF:SI:LP\t0.50:0.95:7
6\t300\trs4\tA\tAC\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP\t0.20:0.95:9
6\t24999999\trs_before_mhc\tG\tA\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP\t0.20:0.95:9
6\t25000000\trs_mhc_start\tG\tA\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP\t0.20:0.95:9
6\t26000000\trs_mhc_middle\tG\tA\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP\t0.20:0.95:9
6\t34000000\trs_mhc_end\tG\tA\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP\t0.20:0.95:9
6\t34000001\trs_after_mhc\tG\tA\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP\t0.20:0.95:9
""",
        encoding="utf-8",
    )


def test_filtering_cli_equals_syntax_overrides_yaml_without_default_leakage():
    args = build_parser().parse_args(
        ["--minimum-info=0.91", "--no-remove-mhc"]
    )

    configuration = resolve_filtering_configuration(args)

    assert configuration.modules.filtering.info_min == 0.91
    assert configuration.modules.filtering.remove_mhc is False
    assert configuration.modules.filtering.maf_min == 0.01


def test_pipeline_resolution_uses_filtering_yaml_values(tmp_path):
    config_path = tmp_path / "run.yaml"
    config_path.write_text(
        """modules:
  filtering:
    info_min: 0.83
    include_indels: true
    remove_palindromic: false
    remove_mhc: false
""",
        encoding="utf-8",
    )
    args = Namespace(
        run_config=config_path,
        vcf=None,
        dataset_id=None,
        output_directory=None,
    )

    module = resolve_filtering_configuration(args).modules.filtering

    assert module.info_min == 0.83
    assert module.include_indels is True
    assert module.remove_palindromic is False
    assert module.remove_mhc is False


@pytest.mark.parametrize(
    ("document", "message"),
    [
        (
            "inputs:\n  dataset_id: ../outside\n",
            "modules.filtering.inputs.dataset_id",
        ),
        (
            "reference_population_tag: 'EUR); touch marker #'\n",
            "valid VCF tag",
        ),
        (
            "output_layout:\n  filtered_vcf: '../{dataset_id}_{genome_build}.vcf.gz'\n",
            "stay inside",
        ),
    ],
)
def test_filtering_configuration_rejects_unsafe_shell_and_output_values(
    tmp_path, document, message,
):
    config_path = tmp_path / "filtering.yaml"
    config_path.write_text(document, encoding="utf-8")

    with pytest.raises(ConfigurationError, match=message):
        load_module_configuration("filtering", config_path)


def test_filtering_configuration_requires_matching_mhc_build(tmp_path):
    config_path = tmp_path / "filtering.yaml"
    config_path.write_text("genome_build: GRCh38\n", encoding="utf-8")

    with pytest.raises(ConfigurationError, match="mhc.genome_build must match"):
        load_module_configuration("filtering", config_path)


def test_vcf_contract_rejects_build_and_active_field_mismatches():
    header = """##fileformat=VCFv4.2
##reference=GRCh37
##contig=<ID=6>
##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY
"""

    with pytest.raises(ValueError, match="does not declare configured genome build"):
        _validate_vcf_contract(
            header=header,
            genome_build="GRCh38",
            genome_build_header_tokens={"GRCh38": ["GRCh38"]},
            required_fields=[],
        )
    with pytest.raises(ValueError, match="FORMAT/SI"):
        _validate_vcf_contract(
            header=header,
            genome_build="GRCh37",
            genome_build_header_tokens={"GRCh37": ["GRCh37"]},
            required_fields=["FORMAT/SI"],
        )

    misleading_header = header.replace(
        'Description="AF"', 'Description="AF derived on GRCh38"'
    )
    with pytest.raises(ValueError, match="does not declare configured genome build"):
        _validate_vcf_contract(
            header=misleading_header,
            genome_build="GRCh38",
            genome_build_header_tokens={"GRCh38": ["GRCh38"]},
            required_fields=[],
        )

    assembly_header = header.replace(
        "##reference=GRCh37\n##contig=<ID=6>",
        "##contig=<ID=6,assembly=GRCh37>",
    )
    assert _validate_vcf_contract(
        header=assembly_header,
        genome_build="GRCh37",
        genome_build_header_tokens={"GRCh37": ["GRCh37"]},
        required_fields=["FORMAT/AF"],
    ) == ["6"]


def test_mhc_filter_fails_when_configured_contig_is_absent():
    with pytest.raises(ValueError, match="MHC contig.*absent"):
        _match_contig_naming("6", ["1", "2"], lambda *_: None)


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX,
    reason="bcftools and tabix are required for the filtering integration test",
)
def test_filtering_runs_bcftools_and_applies_active_rules(tmp_path):
    vcf = tmp_path / "input.vcf"
    output = tmp_path / "output"
    _write_filtering_vcf(vcf)
    args = Namespace(vcf=vcf, dataset_id="study", output_directory=output)

    result = run_sumstat_filter_direct(args)

    filtered = Path(result["filtered_vcf"])
    assert filtered.name == "study_GRCh37_filtered.vcf.gz"
    assert filtered.is_file() and filtered.stat().st_size > 0
    assert Path(result["filtered_vcf_index"]).is_file()
    query = subprocess.run(
        [BCFTOOLS, "query", "-f", "%ID\\n", str(filtered)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert query.stdout.splitlines() == ["rs1", "rs_before_mhc", "rs_after_mhc"]
    assert Path(result["filter_reason_summary"]).is_file()
    log_text = Path(result["filter_log"]).read_text(encoding="utf-8")
    assert "genome_build=GRCh37" in log_text
    assert "STATUS: COMPLETED" in log_text


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX or not shutil.which("false"),
    reason="bcftools, tabix, and false are required for the failure integration test",
)
def test_filtering_failure_keeps_final_output_absent_and_writes_log(tmp_path):
    vcf = tmp_path / "input.vcf"
    output = tmp_path / "output"
    config_path = tmp_path / "run.yaml"
    _write_filtering_vcf(vcf)
    config_path.write_text(
        "config_version: 1\nresources:\n  executables:\n    tabix: %s\n"
        % shutil.which("false"),
        encoding="utf-8",
    )
    args = Namespace(
        run_config=config_path,
        vcf=vcf,
        dataset_id="failed_study",
        output_directory=output,
    )

    with pytest.raises(RuntimeError, match="filtering pipeline failed"):
        run_sumstat_filter_direct(args)

    final_vcf = output / "failed_study_GRCh37_filtered.vcf.gz"
    log_path = output / "logs" / (
        "failed_study_GRCh37_filter_gwas_vcf_bcftools.log"
    )
    assert not final_vcf.exists()
    assert not Path(str(final_vcf) + ".tbi").exists()
    assert log_path.is_file()
    assert "STATUS: FAILED" in log_path.read_text(encoding="utf-8")
    assert not list(output.glob(".*.vcf.gz*"))
