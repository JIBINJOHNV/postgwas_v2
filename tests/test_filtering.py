"""Filtering configuration, scientific-contract, and execution regressions."""

from argparse import Namespace
import csv
from pathlib import Path
import shutil
import subprocess
import sys
from unittest.mock import patch

import pytest

from postgwas.config import load_configuration, load_module_configuration
from postgwas.config.models.modules.qc_summary import QCSummaryConfig
from postgwas.core.errors import ConfigurationError
from postgwas.core.vcf import required_vcf_field_presence, validate_vcf_header_contract
from postgwas.modules.filtering.cli import build_parser
from postgwas.modules.filtering.service import (
    resolve_filtering_configuration,
    run_sumstat_filter_direct,
)
from postgwas.modules.filtering.sumstat_filter import (
    _filtering_field_provenance,
    _match_contig_naming,
    _summarize_filter_tags,
)
from postgwas.modules.qc_summary.assessment import run_vcf_qc_assessment
from postgwas.pipeline.runners import run_sumstat_filter_runner


BCFTOOLS = shutil.which("bcftools")
TABIX = shutil.which("tabix")


def _write_filtering_vcf(path: Path, build: str = "GRCh37") -> None:
    regions = {
        "GRCh37": (28477797, 33448354),
        "GRCh38": (28510120, 33480577),
    }
    mhc_start, mhc_end = regions[build]
    path.write_text(
        """##fileformat=VCFv4.2
##genome_build={build}
##postgwas_version="test-version"
##postgwas_dataset_id="study"
##postgwas_vcf_status="complete"
##postgwas_af_source="study column"
##postgwas_info_source="external user-provided reference"
##postgwas_population_af_fields="INFO/AF | INFO/EUR"
##postgwas_pvalue_harmonisation="validated negative-log10 P-value"
##contig=<ID=6,length=171115067>
##INFO=<ID=AF,Number=A,Type=Float,Description="Study allele frequency">
##INFO=<ID=EUR,Number=A,Type=Float,Description="Reference allele frequency">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Study allele frequency">
##FORMAT=<ID=SI,Number=1,Type=Float,Description="Imputation quality">
##FORMAT=<ID=LP,Number=1,Type=Float,Description="Negative log10 P-value">
##FORMAT=<ID=NEF,Number=1,Type=Float,Description="Effective sample size">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY
6\t100\trs1\tA\tC\t.\tPASS\tAF=0.10;EUR=0.12\tAF:SI:LP:NEF\t0.10:0.90:8:1000
6\t200\trs2\tA\tT\t.\tPASS\tAF=0.50;EUR=0.48\tAF:SI:LP:NEF\t0.50:0.95:7:1000
6\t300\trs4\tA\tAC\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP:NEF\t0.20:0.95:9:1000
6\t400\trs_mixed\tA\tT,AT\t.\tPASS\tAF=0.20,0.10;EUR=0.21,0.11\tAF:SI:LP:NEF\t0.20,0.10:0.95:9:1000
6\t{before}\trs_before_mhc\tG\tA\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP:NEF\t0.20:0.95:9:1000
6\t{mhc_start}\trs_mhc_start\tG\tA\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP:NEF\t0.20:0.95:9:1000
6\t30000000\trs_mhc_middle\tG\tA\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP:NEF\t0.20:0.95:9:1000
6\t{mhc_end}\trs_mhc_end\tG\tA\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP:NEF\t0.20:0.95:9:1000
6\t{after}\trs_after_mhc\tG\tA\t.\tPASS\tAF=0.20;EUR=0.21\tAF:SI:LP:NEF\t0.20:0.95:9:1000
""".format(
            build=build,
            before=mhc_start - 1,
            mhc_start=mhc_start,
            mhc_end=mhc_end,
            after=mhc_end + 1,
        ),
        encoding="utf-8",
    )


def test_filtering_cli_equals_syntax_overrides_yaml_without_default_leakage():
    args = build_parser().parse_args(
        [
            "--minimum-info=0.91",
            "--missing-af-action=keep",
            "--missing-pvalue-action=keep",
            "--remove-mhc",
            "--write-soft-filter-vcf",
        ]
    )

    configuration = resolve_filtering_configuration(args)

    assert configuration.modules.filtering.info_min == 0.91
    assert configuration.modules.filtering.missing_af_action == "keep"
    assert configuration.modules.filtering.missing_pvalue_action == "keep"
    assert configuration.modules.filtering.remove_mhc is True
    assert configuration.modules.filtering.maf_min == 0.01
    assert configuration.modules.filtering.write_soft_filter_vcf is True


def test_filtering_cli_uses_single_presence_flags():
    parser = build_parser()
    help_text = parser.format_help()
    omitted = parser.parse_args([])

    for option in (
        "include-indels",
        "remove-palindromic",
        "remove-mhc",
        "write-soft-filter-vcf",
    ):
        assert "--%s" % option in help_text
        assert "--%s BOOL" % option not in help_text
        assert "--no-%s" % option not in help_text
        assert not hasattr(omitted, option.replace("-", "_"))
    assert "Default: indels are not included" in help_text
    assert "Default: palindromic variants are not removed" in help_text
    assert "Default: MHC variants are not removed" in help_text

    args = parser.parse_args([
        "--include-indels",
        "--remove-palindromic",
        "--remove-mhc",
        "--write-soft-filter-vcf",
    ])
    assert args.include_indels is True
    assert args.remove_palindromic is True
    assert args.remove_mhc is True
    assert args.write_soft_filter_vcf is True

    with pytest.raises(SystemExit):
        parser.parse_args(["--remove-mhc", "false"])


def test_filtering_cli_resolves_an_input_build_mhc_region_override():
    args = build_parser().parse_args([
        "--remove-mhc",
        "--mhc-chrom", "chr6",
        "--mhc-start", "29000000",
        "--mhc-end", "33000000",
    ])

    module = resolve_filtering_configuration(args).modules.filtering

    assert module.remove_mhc is True
    assert module.mhc_region_override.model_dump() == {
        "chromosome": "chr6",
        "start": 29000000,
        "end": 33000000,
    }
    assert module.mhc_regions["GRCh37"].start == 28477797
    assert module.mhc_regions["GRCh38"].start == 28510120


def test_filtering_rejects_an_inactive_mhc_region_override():
    args = build_parser().parse_args(["--mhc-start", "29000000"])

    with pytest.raises(
        ConfigurationError,
        match="mhc_region_override requires remove_mhc: true or --remove-mhc",
    ):
        resolve_filtering_configuration(args)


def test_filtering_cli_uses_the_configured_bcftools_without_a_public_override():
    parser = build_parser()

    assert "--bcftools" not in parser.format_help()
    with pytest.raises(SystemExit):
        parser.parse_args(["--bcftools", "/custom/bcftools"])

    pipeline_help = subprocess.run(
        [
            sys.executable,
            "-m",
            "postgwas",
            "pipeline",
            "--modules",
            "sumstat_filter",
            "--help",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    assert "--bcftools" not in pipeline_help


def test_filtering_accepts_mach_rsq_maximum_through_two():
    args = build_parser().parse_args(["--maximum-info=2"])

    module = resolve_filtering_configuration(args).modules.filtering

    assert module.info_max == pytest.approx(2.0)


def test_qc_and_filtering_preserve_their_intentional_info_max_defaults():
    assert load_module_configuration("qc_summary").rules.info_max == 1.05
    assert load_module_configuration("filtering").info_max is None


def test_filtering_rejects_info_maximum_above_two():
    args = build_parser().parse_args(["--maximum-info=2.0001"])

    with pytest.raises(ConfigurationError, match="less than or equal to 2"):
        resolve_filtering_configuration(args)


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
    assert module.mhc_region_override.model_dump() == {
        "chromosome": None, "start": None, "end": None,
    }


def test_pipeline_keeps_hard_filter_vcf_as_the_downstream_input(
    tmp_path, monkeypatch,
):
    from postgwas.core.contracts import RunContext
    from postgwas.pipeline import runners

    args = Namespace(
        output_directory=str(tmp_path),
        vcf="input.vcf.gz",
        _step_num=2,
    )
    context = RunContext()
    outputs = {
        "filtered_vcf": "hard_filtered.vcf.gz",
        "soft_filtered_vcf": "soft_filtered.vcf.gz",
    }

    validated = []
    monkeypatch.setattr(
        runners,
        "_validate_current_pipeline_vcf",
        lambda current_args, _ctx: validated.append(current_args.vcf),
    )
    with patch(
        "postgwas.modules.filtering.service.run_sumstat_filter_direct",
        return_value=outputs,
    ):
        result = run_sumstat_filter_runner(args, context)

    assert result == outputs
    assert context["sumstat_filter"] == outputs
    assert args.vcf == "hard_filtered.vcf.gz"
    assert validated == ["hard_filtered.vcf.gz"]
    assert args.output_directory == str(tmp_path)


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
        (
            "filter_reason_ids:\n  missing_pvalue: PGWAS_MAF\n",
            "filter reason IDs must be unique",
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


def test_filtering_configuration_defines_build_specific_mhc_regions():
    module = load_module_configuration("filtering")

    assert module.mhc_regions["GRCh37"].model_dump() == {
        "chromosome": "6", "start": 28477797, "end": 33448354,
    }
    assert module.mhc_regions["GRCh38"].model_dump() == {
        "chromosome": "6", "start": 28510120, "end": 33480577,
    }
    assert module.write_soft_filter_vcf is False
    assert module.remove_palindromic is False
    assert module.remove_mhc is False
    assert module.filter_reason_ids.non_snp == "PGWAS_NON_SNP"
    assert len(set(module.filter_reason_ids.model_dump().values())) == 13
    assert module.output_layout.summary_csv.endswith("_filter_summary.csv")
    assert module.output_layout.html_report.endswith("_filter_report.html")


def test_vcf_contract_rejects_build_and_active_field_mismatches():
    header = """##fileformat=VCFv4.2
##genome_build=GRCh37
##contig=<ID=6>
##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY
"""

    with pytest.raises(ValueError, match="unsupported genome build metadata"):
        validate_vcf_header_contract(
            header=header.replace("GRCh37", "OtherBuild"),
            genome_build_header="##genome_build={build}",
            supported_genome_builds=["GRCh37", "GRCh38"],
            required_fields=[],
        )
    with pytest.raises(ValueError, match="FORMAT/SI"):
        validate_vcf_header_contract(
            header=header,
            genome_build_header="##genome_build={build}",
            supported_genome_builds=["GRCh37", "GRCh38"],
            required_fields=["FORMAT/SI"],
        )

    misleading_header = header.replace(
        "##genome_build=GRCh37\n",
        "##reference=GRCh37\n",
    )
    with pytest.raises(ValueError, match="exactly one ##genome_build=<build>"):
        validate_vcf_header_contract(
            header=misleading_header,
            genome_build_header="##genome_build={build}",
            supported_genome_builds=["GRCh37", "GRCh38"],
            required_fields=[],
        )

    duplicate_header = header.replace(
        "##genome_build=GRCh37",
        "##genome_build=GRCh37\n##genome_build=GRCh38",
    )
    with pytest.raises(ValueError, match="found 2"):
        validate_vcf_header_contract(
            header=duplicate_header,
            genome_build_header="##genome_build={build}",
            supported_genome_builds=["GRCh37", "GRCh38"],
            required_fields=[],
        )

    assert validate_vcf_header_contract(
        header=header,
        genome_build_header="##genome_build={build}",
        supported_genome_builds=["GRCh37", "GRCh38"],
        required_fields=["FORMAT/AF"],
    ) == ("GRCh37", ["6"])
    assert required_vcf_field_presence(
        header, ["FORMAT/AF", "FORMAT/SI"],
    ) == {"FORMAT/AF": True, "FORMAT/SI": False}


def test_filtering_reports_optional_postgwas_provenance_as_unavailable():
    provenance = load_configuration().modules.harmonisation.vcf_processing.provenance

    evidence, field_summaries = _filtering_field_provenance(
        header="##fileformat=VCFv4.2\n",
        provenance_headers=provenance.headers,
        provenance_missing_value=provenance.missing_value,
        pvalue_field=None,
        study_af_field="FORMAT/AF",
        study_info_af_field="INFO/AF",
        external_info_af_field="INFO/EUR",
        imputation_field="FORMAT/SI",
        af_uses=["MAF"],
    )

    assert evidence["status"] == "UNAVAILABLE"
    assert "not declared" in evidence["detail"]
    assert field_summaries == [
        {
            "kind": "genetic",
            "label": "Allele-frequency fields",
            "value": (
                "FORMAT/AF for MAF; INFO/AF versus INFO/EUR for external AF "
                "concordance"
            ),
        },
        {
            "kind": "analysis",
            "label": "Imputation quality score",
            "value": "FORMAT/SI · source unavailable",
        },
    ]


def test_mhc_filter_fails_when_configured_contig_is_absent():
    with pytest.raises(ValueError, match="MHC contig.*absent"):
        _match_contig_naming("6", ["1", "2"], lambda *_: None)


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX,
    reason="bcftools and tabix are required for the filtering integration test",
)
@pytest.mark.parametrize(
    ("build", "mhc_start", "mhc_end"),
    [
        ("GRCh37", 28477797, 33448354),
        ("GRCh38", 28510120, 33480577),
    ],
)
def test_filtering_runs_bcftools_and_applies_build_specific_mhc_rule(
    tmp_path, build, mhc_start, mhc_end, capsys,
):
    vcf = tmp_path / "input.vcf"
    output = tmp_path / "output"
    _write_filtering_vcf(vcf, build)
    args = Namespace(
        vcf=vcf,
        dataset_id="study",
        output_directory=output,
        threads=1,
        minimum_neglog10_p=8.0,
        remove_palindromic=True,
        remove_mhc=True,
    )

    result = run_sumstat_filter_direct(args)

    filtered = Path(result["filtered_vcf"])
    assert result["genome_build"] == build
    assert filtered.name == "study_%s_filtered.vcf.gz" % build
    assert filtered.is_file() and filtered.stat().st_size > 0
    assert Path(result["filtered_vcf_index"]).is_file()
    query = subprocess.run(
        [BCFTOOLS, "query", "-f", "%ID\\n", str(filtered)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert query.stdout.splitlines() == ["rs1", "rs_before_mhc", "rs_after_mhc"]
    assert Path(result["mhc_exclusion_bed"]).read_text(encoding="utf-8") == (
        "6\t%d\t%d\n" % (mhc_start - 1, mhc_end)
    )
    assert Path(result["filter_reason_summary"]).is_file()
    with Path(result["filter_reason_summary"]).open(encoding="utf-8") as handle:
        report = handle.read()
    assert "PGWAS_NON_SNP" in report
    assert "Indels and other non-SNP variants\tremove\t2\t2" in report
    summary_csv = Path(result["filter_summary_csv"])
    html_report = Path(result["filter_html_report"])
    assert summary_csv.name == "study_%s_filter_summary.csv" % build
    assert summary_csv.is_file() and summary_csv.stat().st_size > 0
    assert html_report.name == "study_%s_filter_report.html" % build
    assert html_report.is_file() and html_report.stat().st_size > 0
    with summary_csv.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["record_type"] == "overall"
    assert rows[0]["variants_before"] == "9"
    assert rows[0]["variants_after"] == "3"
    assert rows[0]["reconciliation_status"] == "PASS"
    assert rows[0]["requested_threads"] == "1"
    assert rows[0]["polars_thread_pool_size"] == "1"
    assert rows[0]["thread_budget_enforced"] == "True"
    html = html_report.read_text(encoding="utf-8")
    assert "Variant filtering report" in html
    assert "PGWAS_NON_SNP" in html
    assert "PGWAS_MHC" in html
    terminal = capsys.readouterr().out
    assert "Summary-statistics filtering progress" in terminal
    assert "Validate input VCF and count variants" in terminal
    assert "Count input VCF variants" not in terminal
    assert "Save and publish filtering reports" in terminal
    assert "All 5 stages completed" in terminal
    assert "Summary CSV" in terminal
    assert "Detailed HTML report" in terminal
    assert terminal.count("Input VCF validation") == 1
    assert "Input VCF" in terminal and "input.vcf" in terminal
    assert "Declared contigs" in terminal and "1" in terminal
    assert "Total variants" in terminal and "9" in terminal
    assert "Header contract" in terminal
    assert "PASSED · 5/5 required INFO/FORMAT fields declared" in terminal
    assert "Field provenance" in terminal and "AVAILABLE" in terminal
    assert "PostGWAS test-version; VCF status complete" in terminal
    assert "Allele-frequency fields" in terminal
    assert "FORMAT/AF for MAF and palindromic ambiguity" in terminal
    assert "INFO/AF versus INFO/EUR" in terminal
    assert "for external AF concordance" in terminal
    assert "Imputation quality score" in terminal
    assert "FORMAT/SI · external user-provided reference" in terminal
    assert "P-value evidence" in terminal
    assert "FORMAT/LP · validated negative-log10 P-value" in terminal
    assert "Filtering plan · applied in this order" in terminal
    assert "1. EXCLUDE · Missing FORMAT/LP" in terminal
    assert "3. EXCLUDE · Missing FORMAT/AF" in terminal
    assert "12. EXCLUDE · Variants in the MHC region" in terminal
    result_groups = (
        "Statistical significance · FORMAT/LP",
        "Allele frequency · FORMAT/AF",
        "Imputation quality · FORMAT/SI",
        "Study/reference frequency concordance · INFO/AF and INFO/EUR",
        "Variant type · TYPE",
        "Palindromic allele ambiguity · REF/ALT and FORMAT/AF",
        "Genomic region · CHROM/POS",
    )
    assert all(terminal.count(group) == 1 for group in result_groups)
    assert [terminal.index(group) for group in result_groups] == sorted(
        terminal.index(group) for group in result_groups
    )
    assert all(group in html for group in result_groups)
    assert terminal.index(
        "Completed 1/5 · Validate input VCF and count variants"
    ) < terminal.index("Input VCF validation")
    assert terminal.index("Input VCF validation") < terminal.index(
        "Filtering plan · applied in this order"
    )
    assert terminal.index(
        "Filtering plan · applied in this order"
    ) < terminal.index("Completed 2/5 · Audit active filter-rule failures")
    log_text = Path(result["filter_log"]).read_text(encoding="utf-8")
    assert "genome_build=%s" % build in log_text
    assert (
        "requested=1; bcftools_compression=1; polars_pool=1; "
        "polars_budget_enforced=True"
    ) in log_text
    assert "STATUS: COMPLETED" in log_text
    assert "Input VCF validation" in log_text
    assert "Header contract" in log_text
    assert "Filtering plan · applied in this order" in log_text
    assert result["thread_runtime"]["requested_threads"] == 1
    assert result["thread_runtime"]["polars_thread_pool_size"] == 1
    assert result["thread_runtime"]["thread_budget_enforced"] is True
    assert result["input_validation"]["status"] == "PASS"
    assert result["input_validation"]["declared_contigs"] == 1
    assert result["input_validation"]["total_variants"] == 9
    assert result["input_validation"]["field_provenance"]["status"] == (
        "AVAILABLE"
    )
    assert all(
        field["declared"]
        for field in result["input_validation"]["required_fields"]
    )
    assert "Input VCF validation" in html
    assert "Filtering plan and rule results" in html
    assert "Required field declarations" in html
    assert "A declared field may still be missing" not in html
    assert "does not mean every record contains a value" in html


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX,
    reason="bcftools and tabix are required for the filtering default test",
)
def test_palindromic_and_mhc_filters_are_opt_in(tmp_path):
    vcf = tmp_path / "input.vcf"
    output = tmp_path / "output"
    _write_filtering_vcf(vcf)

    result = run_sumstat_filter_direct(Namespace(
        vcf=vcf,
        dataset_id="study",
        output_directory=output,
        threads=1,
    ))

    query = subprocess.run(
        [BCFTOOLS, "query", "-f", "%ID\\n", result["filtered_vcf"]],
        check=True,
        capture_output=True,
        text=True,
    )
    assert query.stdout.splitlines() == [
        "rs1",
        "rs2",
        "rs_before_mhc",
        "rs_mhc_start",
        "rs_mhc_middle",
        "rs_mhc_end",
        "rs_after_mhc",
    ]
    assert result["mhc_exclusion_bed"] is None
    assert result["variant_qc_policy"]["remove_palindromic"] is False
    assert result["variant_qc_policy"]["remove_mhc"] is False
    reason_report = Path(result["filter_reason_summary"]).read_text(
        encoding="utf-8",
    )
    assert "PGWAS_PALINDROMIC" not in reason_report
    assert "PGWAS_MHC" not in reason_report


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX,
    reason="bcftools and tabix are required for the MHC override test",
)
def test_mhc_cli_region_overrides_the_inferred_build_default(tmp_path):
    vcf = tmp_path / "input.vcf"
    output = tmp_path / "output"
    _write_filtering_vcf(vcf)
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--dataset-id", "study",
        "--output-directory", str(output),
        "--remove-mhc",
        "--mhc-chrom", "6",
        "--mhc-start", "99",
        "--mhc-end", "101",
    ])

    result = run_sumstat_filter_direct(args)

    query = subprocess.run(
        [BCFTOOLS, "query", "-f", "%ID\\n", result["filtered_vcf"]],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "rs1" not in query.stdout.splitlines()
    assert "rs_mhc_middle" in query.stdout.splitlines()
    assert result["variant_qc_policy"]["mhc_chromosome"] == "6"
    assert result["variant_qc_policy"]["mhc_start"] == 99
    assert result["variant_qc_policy"]["mhc_end"] == 101
    assert Path(result["mhc_exclusion_bed"]).read_text(encoding="utf-8") == (
        "6\t98\t101\n"
    )
    log = Path(result["filter_log"]).read_text(encoding="utf-8")
    assert "MHC region selected: 6:99-101 (GRCh37; custom override)" in log


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX,
    reason="bcftools and tabix are required for the QC/filter parity test",
)
def test_qc_virtual_subset_matches_filtering_under_the_same_policy(tmp_path):
    vcf = tmp_path / "input.vcf"
    _write_filtering_vcf(vcf)
    full_configuration = load_configuration()
    qc_document = full_configuration.modules.qc_summary.model_dump()
    qc_document["rules"].update({
        "remove_palindromic": True,
        "remove_mhc": True,
    })
    qc_configuration = QCSummaryConfig.model_validate(qc_document)
    qc = run_vcf_qc_assessment(
        vcf_path=vcf,
        output_directory=tmp_path / "qc_output",
        dataset_id="study",
        external_af_name="EUR",
        configuration=qc_configuration,
        bcftools_bin=BCFTOOLS,
        genome_build_header=(
            full_configuration.modules.harmonisation.vcf_processing.genome_build_header
        ),
        supported_genome_builds=tuple(full_configuration.resources.genomes),
        threads=1,
    )
    filtered = run_sumstat_filter_direct(Namespace(
        vcf=vcf,
        dataset_id="study",
        output_directory=tmp_path / "filter_output",
        maximum_info=1.05,
        remove_palindromic=True,
        remove_mhc=True,
        threads=1,
    ))

    assert qc["variant_qc_policy"] == filtered["variant_qc_policy"]
    assert qc["qc_passed"]["num_records"] == filtered["variants_after"] == 3
    query = subprocess.run(
        [BCFTOOLS, "query", "-f", "%ID\\n", filtered["filtered_vcf"]],
        check=True,
        capture_output=True,
        text=True,
    )
    assert query.stdout.splitlines() == [
        "rs1", "rs_before_mhc", "rs_after_mhc",
    ]


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX,
    reason="bcftools and tabix are required for the soft-filter integration test",
)
def test_soft_filter_vcf_retains_every_record_and_uses_semantic_reason_ids(
    tmp_path,
):
    vcf = tmp_path / "input.vcf"
    output = tmp_path / "output"
    _write_filtering_vcf(vcf)
    args = Namespace(
        vcf=vcf,
        dataset_id="study",
        output_directory=output,
        remove_palindromic=True,
        remove_mhc=True,
        write_soft_filter_vcf=True,
    )

    result = run_sumstat_filter_direct(args)

    soft_vcf = Path(result["soft_filtered_vcf"])
    assert soft_vcf.name == "study_GRCh37_soft_filtered.vcf.gz"
    assert soft_vcf.is_file()
    assert Path(result["soft_filtered_vcf_index"]).is_file()
    query = subprocess.run(
        [BCFTOOLS, "query", "-f", "%ID\t%FILTER\n", str(soft_vcf)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert query.stdout.splitlines() == [
        "rs1\tPASS",
        "rs2\tPGWAS_PALINDROMIC",
        "rs4\tPGWAS_NON_SNP",
        "rs_mixed\tPGWAS_NON_SNP",
        "rs_before_mhc\tPASS",
        "rs_mhc_start\tPGWAS_MHC",
        "rs_mhc_middle\tPGWAS_MHC",
        "rs_mhc_end\tPGWAS_MHC",
        "rs_after_mhc\tPASS",
    ]
    header = subprocess.run(
        [BCFTOOLS, "view", "--header-only", str(soft_vcf)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    assert "##FILTER=<ID=PGWAS_NON_SNP" in header
    assert "##FILTER=<ID=PGWAS_MHC" in header
    assert result["variants_before"] == 9
    assert result["variants_after"] == 3
    assert result["filter_reason_reconciled"] is True


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX,
    reason="bcftools and tabix are required for the empty-policy test",
)
def test_empty_match_none_policy_is_reconciled_and_soft_tagged(tmp_path):
    vcf = tmp_path / "input.vcf"
    output = tmp_path / "output"
    config_path = tmp_path / "run.yaml"
    _write_filtering_vcf(vcf)
    config_path.write_text(
        """config_version: 1
modules:
  filtering:
    maf_min: null
    info_min: null
    info_max: null
    minimum_neglog10_p: null
    frequency_difference_max: null
    include_indels: true
    remove_palindromic: false
    remove_mhc: false
    empty_expression_action: match_none
    write_soft_filter_vcf: true
""",
        encoding="utf-8",
    )
    args = Namespace(
        run_config=config_path,
        vcf=vcf,
        dataset_id="study",
        output_directory=output,
    )

    result = run_sumstat_filter_direct(args)

    assert result["variants_before"] == 9
    assert result["variants_after"] == 0
    assert result["filter_reason_reconciled"] is True
    assert result["filter_summary"]["rules"][0]["display_group"] == (
        "Filtering policy"
    )
    query = subprocess.run(
        [
            BCFTOOLS,
            "query",
            "-f",
            "%FILTER\n",
            result["soft_filtered_vcf"],
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert query.stdout.splitlines() == ["PGWAS_EMPTY_POLICY"] * 9


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX,
    reason="bcftools and tabix are required for the reconciliation test",
)
def test_filtering_does_not_publish_when_reason_accounting_disagrees(tmp_path):
    vcf = tmp_path / "input.vcf"
    output = tmp_path / "output"
    _write_filtering_vcf(vcf)
    args = Namespace(vcf=vcf, dataset_id="study", output_directory=output)

    def incorrect_summary(*args, **kwargs):
        statistics = _summarize_filter_tags(*args, **kwargs)
        statistics["primary_removed_total"] += 1
        return statistics

    with patch(
        "postgwas.modules.filtering.sumstat_filter._summarize_filter_tags",
        side_effect=incorrect_summary,
    ), pytest.raises(RuntimeError, match="reconciliation failed"):
        run_sumstat_filter_direct(args)

    assert not (output / "study_GRCh37_filtered.vcf.gz").exists()
    assert not (
        output / "qc_summary" / "study_GRCh37_filter_reason_summary.tsv"
    ).exists()
    assert not (
        output / "qc_summary" / "study_GRCh37_filter_summary.csv"
    ).exists()
    assert not (
        output / "reports" / "study_GRCh37_filter_report.html"
    ).exists()
    log = (
        output / "logs" / "study_GRCh37_filter_gwas_vcf_bcftools.log"
    ).read_text(encoding="utf-8")
    assert "No filtering outputs were published" in log
    assert "STATUS: FAILED" in log


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX,
    reason="bcftools and tabix are required for the FILTER collision test",
)
def test_filtering_rejects_existing_configured_filter_ids(tmp_path):
    vcf = tmp_path / "input.vcf"
    output = tmp_path / "output"
    _write_filtering_vcf(vcf)
    vcf.write_text(
        vcf.read_text(encoding="utf-8").replace(
            "##genome_build=GRCh37\n",
            "##genome_build=GRCh37\n"
            "##FILTER=<ID=PGWAS_MAF,Description=\"Existing annotation\">\n",
        ),
        encoding="utf-8",
    )
    args = Namespace(vcf=vcf, dataset_id="study", output_directory=output)

    with pytest.raises(ValueError, match="already declares.*PGWAS_MAF"):
        run_sumstat_filter_direct(args)

    assert not list(output.glob("*_filtered.vcf.gz*"))


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX or not shutil.which("false"),
    reason="bcftools, tabix, and false are required for the failure integration test",
)
def test_filtering_failure_keeps_final_output_absent_and_writes_log(
    tmp_path, capsys,
):
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
    terminal = capsys.readouterr().out
    assert "Failed 3/5 · Create and index hard-filtered VCF" in terminal
    assert "All 5 stages completed" not in terminal


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX,
    reason="bcftools and tabix are required for the filtering validation test",
)
def test_filtering_rejects_missing_authoritative_build_header_and_writes_log(
    tmp_path,
):
    vcf = tmp_path / "input.vcf"
    output = tmp_path / "output"
    _write_filtering_vcf(vcf)
    vcf.write_text(
        vcf.read_text(encoding="utf-8").replace(
            "##genome_build=GRCh37", "##reference=GRCh37"
        ),
        encoding="utf-8",
    )
    args = Namespace(vcf=vcf, dataset_id="invalid_build", output_directory=output)

    with pytest.raises(ValueError, match="exactly one ##genome_build=<build>"):
        run_sumstat_filter_direct(args)

    log_path = output / "logs" / "invalid_build_filter_preflight.log"
    assert log_path.is_file()
    log_text = log_path.read_text(encoding="utf-8")
    assert "Input VCF preflight failed" in log_text
    assert "STATUS: FAILED" in log_text
    assert not list(output.glob("*_filtered.vcf.gz*"))


@pytest.mark.skipif(
    not BCFTOOLS or not TABIX,
    reason="bcftools and tabix are required for the field-validation test",
)
def test_filtering_reports_each_required_header_field_before_failure(
    tmp_path, capsys,
):
    vcf = tmp_path / "input.vcf"
    output = tmp_path / "output"
    _write_filtering_vcf(vcf)
    vcf.write_text(
        vcf.read_text(encoding="utf-8").replace(
            "##INFO=<ID=EUR,Number=A,Type=Float,"
            'Description="Reference allele frequency">\n',
            "",
        ),
        encoding="utf-8",
    )
    args = Namespace(vcf=vcf, dataset_id="missing_field", output_directory=output)

    with pytest.raises(ValueError, match="INFO/EUR"):
        run_sumstat_filter_direct(args)

    terminal = capsys.readouterr().out
    assert "Input VCF validation" in terminal
    assert "Header contract" in terminal
    assert "FAILED · 3/4 required INFO/FORMAT fields declared" in terminal
    assert "FORMAT/AF" in terminal and "present" in terminal
    assert "INFO/EUR" in terminal and "missing" in terminal
    assert "Failed 1/5 · Validate input VCF and count variants" in terminal
    assert "Count input VCF variants" not in terminal
    log = output / "logs" / "missing_field_filter_preflight.log"
    log_text = log.read_text(encoding="utf-8")
    assert "Header contract" in log_text
    assert "INFO/EUR" in log_text and "missing" in log_text
    assert "STATUS: FAILED" in log_text
    assert not list(output.glob("*_filtered.vcf.gz*"))
