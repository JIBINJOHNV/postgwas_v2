"""Scientific contract tests for the streaming merged-VCF QC assessment."""

import csv
from io import StringIO
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
from unittest.mock import patch

import polars as pl
import pytest
from rich.cells import cell_len
from rich.console import Console
import yaml
from pydantic import ValidationError

from postgwas.config import load_configuration
from postgwas.config.models.modules.qc_summary import QCSummaryConfig
from postgwas.core.errors import ConfigurationError
from postgwas.core.dataframes import collect_streaming
from postgwas.core.ui import StageProgress
from postgwas.core.ui.screen import style_screen_block, terminal_style
from postgwas.modules.harmonisation.qc_reporting import (
    harmonisation_qc_summary_lines,
    harmonisation_qc_takeaway_lines,
    summarise_strand_orientation,
)
from postgwas.modules.qc_summary.assessment import (
    VcfAssessmentError,
    assess_variant_table,
    extract_vcf_assessment_table,
    resolve_qc_rule_contract,
    run_vcf_qc_assessment,
    validate_vcf_assessment_header,
)
from postgwas.modules.qc_summary.cli import build_parser
from postgwas.modules.qc_summary.reporting import (
    QC_SUMMARY_COLUMNS,
    qc_summary_lines,
    qc_summary_screen_lines,
)
from postgwas.modules.qc_summary.service import (
    _emit_summary,
    resolve_qc_summary_configuration,
    run_qc_summary_direct,
)


def _assessment_table(path: Path) -> Path:
    path.write_text(
        "CHROM\tPOS\tREF\tALT\tINFO_AF\tINFO_EXTERNAL_AF\t"
        "FORMAT_AF\tFORMAT_SI\tFORMAT_LP\tFORMAT_NEF\n"
        "1\t100\tA\tG\t0.2\t0.21\t0.2\t0.9\t8\t100\n"
        "1\t101\tC\tT\t0.3\t.\t0.3\t0.9\t8\t100\n"
        "1\t102\tA\tC\t0.001\t0.001\t0.001\t0.9\t8\t100\n"
        "1\t103\tG\tA\t0.2\t0.2\t0.2\t0.5\t8\t100\n"
        "1\t104\tA\tAT\t0.2\t0.2\t0.2\t0.9\t8\t100\n"
        "1\t105\tA\tT\t0.5\t0.5\t0.5\t0.9\t8\t100\n"
        "6\t30000000\tA\tG\t0.2\t0.2\t0.2\t0.9\t8\t100\n"
        "1\t106\tA\tAT\t0.001\t0.5\t0.001\t0.5\t8\t1000\n",
        encoding="utf-8",
    )
    return path


def _sample_size_table(path: Path, values: list[str]) -> Path:
    header = (
        "CHROM\tPOS\tREF\tALT\tINFO_AF\tINFO_EXTERNAL_AF\t"
        "FORMAT_AF\tFORMAT_SI\tFORMAT_LP\tFORMAT_NEF\n"
    )
    rows = [
        "1\t%s\tA\tG\t0.2\t0.2\t0.2\t0.9\t8\t%s" % (100 + index, value)
        for index, value in enumerate(values)
    ]
    path.write_text(header + "\n".join(rows) + "\n", encoding="utf-8")
    return path


def _vcf_without_declarations(path: Path, *fields: str) -> Path:
    source = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    omitted = set(fields)
    lines = []
    for line in source.read_text(encoding="utf-8").splitlines():
        if any(
            line.startswith("##%s=<ID=%s," % tuple(field.split("/", 1)))
            for field in omitted
        ):
            continue
        lines.append(line)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class _RecordLogger:
    def __init__(self):
        self.records = []

    def record(self, *args, **kwargs):
        self.records.append((args, kwargs))


def _qc_config(**rule_overrides) -> QCSummaryConfig:
    module = load_configuration().modules.qc_summary
    if not rule_overrides:
        return module
    document = module.model_dump()
    document["rules"].update(rule_overrides)
    return QCSummaryConfig.model_validate(document)


def _provenance_headers() -> dict[str, str]:
    return dict(
        load_configuration()
        .modules.harmonisation.vcf_processing.provenance.headers
    )


def _vcf_header_contract() -> dict[str, object]:
    configuration = load_configuration()
    return {
        "genome_build_header": (
            configuration.modules.harmonisation.vcf_processing.genome_build_header
        ),
        "supported_genome_builds": tuple(configuration.resources.genomes),
    }


def _assess(path: Path, configuration: QCSummaryConfig | None = None):
    configuration = configuration or _qc_config(
        remove_palindromic=True,
        remove_mhc=True,
    )
    return assess_variant_table(
        path,
        configuration,
        "GRCh37",
        configuration.reference_af_column,
        configuration.vcf_fields.model_dump(),
        delimiter=configuration.table.delimiter,
        null_values=configuration.table.null_values,
    )


def test_assessment_applies_every_rule_to_raw_data_and_combines_failures(tmp_path):
    assessment = _assess(
        _assessment_table(tmp_path / "raw_vcf.tsv"),
    )

    assert assessment["raw"]["num_records"] == 8
    assert assessment["raw"]["num_snps"] == 6
    assert assessment["raw"]["num_non_snps"] == 2
    assert assessment["raw"]["ts_tv_ratio"] == 2.0
    assert assessment["raw"]["external_af_missing"] == 1
    assert assessment["raw"]["af_comparable"] == 7
    assert assessment["raw"]["af_difference_above_cutoff"] == 1
    assert assessment["raw"]["effective_sample_size_available"] == 8
    assert assessment["raw"]["effective_sample_size_mean"] == pytest.approx(212.5)
    assert assessment["raw"]["effective_sample_size_above_outlier_threshold"] == 0
    assert assessment["raw"][
        "effective_sample_size_reference_quantile_value"
    ] == pytest.approx(370.0)
    assert assessment["raw"][
        "effective_sample_size_minimum_threshold"
    ] == pytest.approx(370.0 * (2.0 / 3.0))
    assert assessment["raw"][
        "effective_sample_size_below_minimum_threshold"
    ] == 7
    assert assessment["raw"][
        "effective_sample_size_below_minimum_threshold_fraction"
    ] == pytest.approx(7.0 / 8.0)

    rules = {rule["key"]: rule for rule in assessment["rules"]}
    assert rules["allele_frequency"]["failed_raw"] == 2
    assert rules["imputation_quality"]["failed_raw"] == 2
    assert rules["external_af_concordance"]["failed_raw"] == 2
    assert rules["variant_type"]["failed_raw"] == 2
    assert rules["palindromic_variants"]["failed_raw"] == 1
    assert rules["mhc_region"]["failed_raw"] == 1
    assert rules["allele_frequency"]["unique_only_raw"] == 1
    assert rules["allele_frequency"]["overlap_raw"] == 1
    assert rules["imputation_quality"]["unique_only_raw"] == 1
    assert rules["imputation_quality"]["overlap_raw"] == 1
    assert rules["palindromic_variants"]["unique_only_raw"] == 1
    assert rules["palindromic_variants"]["overlap_raw"] == 0
    assert rules["imputation_quality"]["category"] == (
        "Score quality and completeness"
    )
    info_details = {
        detail["key"]: detail
        for detail in rules["imputation_quality"]["details"]
    }
    assert info_details["missing"]["matched_raw"] == 0
    assert info_details["below_minimum"]["matched_raw"] == 2
    assert info_details["above_maximum"]["matched_raw"] == 0
    assert info_details["below_minimum"]["decision"] == (
        "exclude_from_virtual_subset"
    )

    assert assessment["qc_passed"]["num_records"] == 1
    assert assessment["qc_passed"]["effective_sample_size_available"] == 1
    assert assessment["qc_passed"]["effective_sample_size_mean"] == 100.0
    assert assessment["qc_passed"]["effective_sample_size_standard_deviation"] is None
    assert assessment["qc_passed"][
        "effective_sample_size_minimum_threshold"
    ] == assessment["raw"]["effective_sample_size_minimum_threshold"]
    assert assessment["qc_passed"][
        "effective_sample_size_below_minimum_threshold"
    ] == 1
    assert assessment["excluded_total"] == 7
    assert assessment["unique_rule_only_total"] == 6
    assert assessment["active_rule_count"] == 6
    assert assessment["accounting_balanced"] is True
    assert assessment["overlap_variants"] == 1
    assert "overlapping_rule_match_total" not in assessment
    assert "rule_match_total" not in assessment
    assert "extra_rule_matches" not in assessment
    assert assessment["inactive_rules"] == [{
        "key": "significance",
        "label": "Association significance",
        "category": "Analysis scope",
        "reason": "No minimum −log10(P) threshold is configured.",
    }]
    assert "no filtered VCF is created" in assessment["definition"]


def test_packaged_qc_boolean_policies_are_inactive_until_opted_in(tmp_path):
    configuration = _qc_config()
    assessment = _assess(
        _assessment_table(tmp_path / "default_policy.tsv"),
        configuration,
    )

    assert configuration.rules.include_indels is False
    assert configuration.rules.remove_palindromic is False
    assert configuration.rules.remove_mhc is False
    assert assessment["active_rule_count"] == 4
    assert {rule["key"] for rule in assessment["inactive_rules"]} == {
        "significance",
        "palindromic_variants",
        "mhc_region",
    }
    assert assessment["qc_passed"]["num_records"] == 3


def test_directional_score_details_separate_low_and_high_values(tmp_path):
    table = _assessment_table(tmp_path / "directional.tsv")
    contents = table.read_text(encoding="utf-8").replace(
        "0.2\t0.9\t8\t100\n",
        "0.2\t1.2\t8\t100\n",
        1,
    )
    table.write_text(contents, encoding="utf-8")

    assessment = _assess(table)
    rule = next(
        rule for rule in assessment["rules"]
        if rule["key"] == "imputation_quality"
    )
    details = {detail["key"]: detail for detail in rule["details"]}

    assert details["below_minimum"]["matched_raw"] == 2
    assert details["above_maximum"]["matched_raw"] == 1
    assert rule["failed_raw"] == 3


def test_sample_size_outliers_are_reported_before_and_after_combined_qc(tmp_path):
    configuration = _qc_config(sample_size_outlier_standard_deviations=1.5)
    assessment = _assess(
        _assessment_table(tmp_path / "raw_vcf.tsv"),
        configuration,
    )

    assert assessment["sample_size_outlier_standard_deviations"] == 1.5
    assert assessment["raw"]["effective_sample_size_outlier_threshold"] == pytest.approx(
        212.5 + 1.5 * 318.1980515339464
    )
    assert assessment["raw"]["effective_sample_size_above_outlier_threshold"] == 1
    assert assessment["qc_passed"]["effective_sample_size_above_outlier_threshold"] == 0


def test_low_neff_uses_one_raw_reference_for_both_assessment_stages(tmp_path):
    assessment = _assess(
        _sample_size_table(
            tmp_path / "low_neff.tsv",
            ["150000"] * 95 + ["50000"] * 5,
        )
    )

    assert assessment["sample_size_reference_quantile"] == 0.9
    assert assessment["sample_size_minimum_fraction_of_reference"] == pytest.approx(
        2.0 / 3.0
    )
    for stage in ("raw", "qc_passed"):
        metrics = assessment[stage]
        assert metrics["effective_sample_size_reference_quantile_value"] == 150000
        assert metrics["effective_sample_size_minimum_threshold"] == pytest.approx(
            100000
        )
        assert metrics["effective_sample_size_below_minimum_threshold"] == 5
        assert metrics[
            "effective_sample_size_below_minimum_threshold_fraction"
        ] == pytest.approx(0.05)
        assert metrics["effective_sample_size_above_outlier_threshold"] == 0

    assert assessment["qc_passed"]["num_records"] == 100
    assert assessment["excluded_total"] == 0


def test_unusable_sample_sizes_are_reported_but_not_used_in_distribution(tmp_path):
    table = _assessment_table(tmp_path / "raw_vcf.tsv")
    contents = table.read_text(encoding="utf-8")
    contents = contents.replace("\t100\n", "\t.\n", 1).replace("\t100\n", "\t0\n", 1)
    table.write_text(contents, encoding="utf-8")
    assessment = _assess(table)

    assert assessment["raw"]["effective_sample_size_available"] == 6
    assert assessment["raw"]["effective_sample_size_missing_or_invalid"] == 2


@pytest.mark.parametrize(
    ("values", "available", "mean"),
    [
        ([".", "0", "-5", "nan", "inf"], 0, None),
        ([".", "0", "100"], 1, 100.0),
    ],
)
def test_two_pass_sample_size_statistics_preserve_sparse_edge_cases(
    tmp_path,
    values,
    available,
    mean,
):
    assessment = _assess(
        _sample_size_table(tmp_path / "sample_sizes.tsv", values)
    )

    for stage in ("raw", "qc_passed"):
        metrics = assessment[stage]
        assert metrics["effective_sample_size_available"] == available
        assert metrics["effective_sample_size_mean"] == mean
        assert metrics["effective_sample_size_standard_deviation"] is None
        assert metrics["effective_sample_size_outlier_threshold"] is None
        assert metrics["effective_sample_size_above_outlier_threshold"] == 0
        assert metrics["effective_sample_size_reference_quantile_value"] == mean
        assert metrics["effective_sample_size_minimum_threshold"] == (
            None if mean is None else pytest.approx(mean * (2.0 / 3.0))
        )
        assert metrics["effective_sample_size_below_minimum_threshold"] == 0
        assert metrics[
            "effective_sample_size_below_minimum_threshold_fraction"
        ] == (None if available == 0 else 0.0)


def test_two_pass_sample_sd_remains_stable_for_nearly_constant_values(tmp_path):
    assessment = _assess(
        _sample_size_table(
            tmp_path / "nearly_constant_sample_sizes.tsv",
            ["150000", "150001", "149999"],
        )
    )

    for stage in ("raw", "qc_passed"):
        metrics = assessment[stage]
        assert metrics["effective_sample_size_mean"] == pytest.approx(150000.0)
        assert metrics["effective_sample_size_standard_deviation"] == pytest.approx(
            1.0
        )


def test_assessment_uses_two_streaming_collections_without_window_broadcasts(
    tmp_path,
):
    with patch(
        "postgwas.modules.qc_summary.assessment.collect_streaming",
        wraps=collect_streaming,
    ) as collect, patch.object(
        pl.Expr,
        "over",
        side_effect=AssertionError("QC aggregation must not use window broadcasts"),
    ):
        assessment = _assess(_assessment_table(tmp_path / "raw_vcf.tsv"))

    assert collect.call_count == 2
    assert assessment["aggregation"]["aggregation_strategy"] == (
        "two_pass_streaming"
    )
    assert assessment["aggregation"]["temporary_table_scans"] == 2
    assert assessment["aggregation"]["streaming_collection_api"] in {
        "engine=streaming",
        "streaming=True",
    }
    assert assessment["aggregation"]["polars_version"] == pl.__version__


def test_bcftools_extracts_all_qc_fields_in_one_data_query(tmp_path):
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("bcftools is not installed")
    vcf = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    table = tmp_path / "assessment.tsv"
    configuration = _qc_config()

    extract_vcf_assessment_table(
        vcf_path=vcf,
        table_path=table,
        dataset_id="study",
        external_af_name="EUR",
        vcf_fields=configuration.vcf_fields.model_dump(),
        table_delimiter=configuration.table.delimiter,
        io_buffer_bytes=configuration.table.io_buffer_bytes,
        bcftools_bin=bcftools,
        **_vcf_header_contract(),
    )

    rows = table.read_text(encoding="utf-8").splitlines()
    assert len(rows) == 3
    assert rows[0].split("\t") == [
        "CHROM", "POS", "REF", "ALT", "INFO_AF",
        "INFO_EXTERNAL_AF", "FORMAT_AF", "FORMAT_SI", "FORMAT_LP", "FORMAT_NEF",
    ]
    assert rows[1].split("\t") == [
        "1", "100", "A", "G", "0.2", "0.21", "0.2", "0.9", "8", "1000",
    ]
    assert rows[2].split("\t")[5] == "."
    assessment = _assess(table)
    assert assessment["raw"]["external_af_missing"] == 1


def test_header_preflight_reports_all_undeclared_fields_before_query(tmp_path):
    vcf = _vcf_without_declarations(
        tmp_path / "missing_fields.vcf",
        "INFO/EUR",
        "FORMAT/SI",
    )
    configuration = _qc_config()
    logger = _RecordLogger()

    with patch(
        "postgwas.modules.qc_summary.assessment.read_vcf_header",
        return_value=vcf.read_text(encoding="utf-8"),
    ), patch(
        "postgwas.modules.qc_summary.assessment.extract_vcf_table",
    ) as extract, pytest.raises(VcfAssessmentError) as error:
        extract_vcf_assessment_table(
            vcf_path=vcf,
            table_path=tmp_path / "assessment.tsv",
            dataset_id="study",
            external_af_name="EUR",
            vcf_fields=configuration.vcf_fields.model_dump(),
            table_delimiter=configuration.table.delimiter,
            io_buffer_bytes=configuration.table.io_buffer_bytes,
            bcftools_bin=shutil.which("bcftools") or "bcftools",
            **_vcf_header_contract(),
            logger=logger,
        )

    extract.assert_not_called()
    message = str(error.value)
    assert "INFO/EUR" in message
    assert "FORMAT/SI" in message
    validation_records = [
        kwargs for args, kwargs in logger.records
        if args == ("VALIDATION", "vcf_qc_header_contract")
    ]
    assert len(validation_records) == 1
    assert validation_records[0]["status"] == "FAILED"
    assert validation_records[0]["missing_fields"] == "INFO/EUR,FORMAT/SI"


def test_scientific_provenance_is_reported_but_not_required_for_contract(
    tmp_path,
):
    source = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    vcf = tmp_path / "without_provenance.vcf"
    vcf.write_text(
        "\n".join(
            line for line in source.read_text(encoding="utf-8").splitlines()
            if not line.startswith("##postgwas_")
        ) + "\n",
        encoding="utf-8",
    )
    configuration = _qc_config()

    with patch(
        "postgwas.modules.qc_summary.assessment.read_vcf_header",
        return_value=vcf.read_text(encoding="utf-8"),
    ):
        validation = validate_vcf_assessment_header(
            vcf_path=vcf,
            external_af_name="EUR",
            vcf_fields=configuration.vcf_fields.model_dump(),
            bcftools_bin="configured-bcftools",
            **_vcf_header_contract(),
            provenance_headers=_provenance_headers(),
        )

    assert validation.as_report()["status"] == "passed"
    assert validation.as_report()["genome_build"] == "GRCh37"
    assert validation.as_report()["genome_build_source"] == "vcf_header"
    assert validation.provenance_report()["status"] == "not_declared"
    assert validation.provenance_report()["values"] == {}


def test_header_preflight_rejects_output_build_provenance_mismatch(tmp_path):
    source = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    vcf = tmp_path / "inconsistent_build.vcf"
    vcf.write_text(
        source.read_text(encoding="utf-8").replace(
            '##postgwas_output_genome_build="GRCh37"',
            '##postgwas_output_genome_build="GRCh38"',
            1,
        ),
        encoding="utf-8",
    )
    configuration = _qc_config()

    with pytest.raises(
        VcfAssessmentError,
        match="assessed GWAS-VCF build GRCh37 conflicts.*GRCh38",
    ):
        validate_vcf_assessment_header(
            vcf_path=vcf,
            external_af_name="EUR",
            vcf_fields=configuration.vcf_fields.model_dump(),
            bcftools_bin="configured-bcftools",
            header=vcf.read_text(encoding="utf-8"),
            **_vcf_header_contract(),
            provenance_headers=_provenance_headers(),
        )


def test_qc_vcf_query_fields_are_configuration_driven(tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"placeholder")
    table = tmp_path / "assessment.tsv"
    configuration = _qc_config()
    fields = configuration.vcf_fields.model_dump()
    fields.update({
        "study_info_af": "%INFO/STUDY_AF",
        "external_info_af": "%INFO/REF_{external_af}",
        "study_format_af": "[%EAF]",
        "imputation_format": "[%INFO_SCORE]",
        "log_pvalue_format": "[%MLOG10P]",
        "effective_sample_size_format": "[%EFFECTIVE_N]",
    })
    header = (
        "##genome_build=GRCh37\n"
        "##INFO=<ID=STUDY_AF,Number=A,Type=Float>\n"
        "##INFO=<ID=REF_EUR,Number=A,Type=Float>\n"
        "##FORMAT=<ID=EAF,Number=A,Type=Float>\n"
        "##FORMAT=<ID=INFO_SCORE,Number=A,Type=Float>\n"
        "##FORMAT=<ID=MLOG10P,Number=A,Type=Float>\n"
        "##FORMAT=<ID=EFFECTIVE_N,Number=1,Type=Float>\n"
    )
    with patch(
        "postgwas.modules.qc_summary.assessment.extract_vcf_table",
        return_value=str(table),
    ) as query, patch(
        "postgwas.modules.qc_summary.assessment.read_vcf_header",
        return_value=header,
    ):
        extract_vcf_assessment_table(
            vcf_path=vcf,
            table_path=table,
            dataset_id="study",
            external_af_name="EUR",
            vcf_fields=fields,
            table_delimiter=configuration.table.delimiter,
            io_buffer_bytes=configuration.table.io_buffer_bytes,
            bcftools_bin="configured-bcftools",
            **_vcf_header_contract(),
        )
    columns = query.call_args.args[3]
    assert columns["INFO_AF"] == "%INFO/STUDY_AF"
    assert columns["INFO_EXTERNAL_AF"] == "%INFO/REF_EUR"
    assert columns["FORMAT_AF"] == "[%EAF]"
    assert columns["FORMAT_NEF"] == "[%EFFECTIVE_N]"
    assert "allow_undefined_tags" not in query.call_args.kwargs


def test_sample_size_outlier_threshold_is_validated_and_defaults_to_five():
    assert _qc_config().rules.sample_size_outlier_standard_deviations == 5.0
    assert _qc_config().rules.sample_size_reference_quantile == 0.9
    assert _qc_config().rules.sample_size_minimum_fraction_of_reference == pytest.approx(
        2.0 / 3.0
    )
    with pytest.raises(ValidationError, match="sample_size_outlier_standard_deviations"):
        _qc_config(sample_size_outlier_standard_deviations=0)
    with pytest.raises(ValidationError, match="sample_size_reference_quantile"):
        _qc_config(sample_size_reference_quantile=0)
    with pytest.raises(
        ValidationError,
        match="sample_size_minimum_fraction_of_reference",
    ):
        _qc_config(sample_size_minimum_fraction_of_reference=1.01)


def test_configuration_rejects_an_invalid_info_interval_before_analysis():
    with pytest.raises(ValidationError, match="info_max"):
        _qc_config(info_min=0.9, info_max=0.8)


def test_qc_info_max_supports_mach_rsq_but_rejects_values_above_two():
    assert _qc_config(info_max=2.0).rules.info_max == 2.0
    with pytest.raises(ValidationError, match="less than or equal to 2"):
        _qc_config(info_max=2.0001)


@pytest.mark.parametrize(
    ("field", "invalid", "message"),
    [
        ("summary_csv", "qc/{dataset_id}_{build}.tsv", "summary_csv"),
        ("html_report", "reports/{dataset_id}_{build}.txt", "html_report"),
    ],
)
def test_qc_report_output_extensions_are_schema_validated(
    field, invalid, message,
):
    document = _qc_config().model_dump()
    document["output_layout"][field] = invalid
    with pytest.raises(ValidationError, match=message):
        QCSummaryConfig.model_validate(document)


def test_mhc_regions_are_build_specific_and_use_grc_coordinates():
    configuration = _qc_config()
    grch37 = configuration.mhc_region("GRCh37")
    grch38 = configuration.mhc_region("GRCh38")

    assert (grch37.chromosome, grch37.start, grch37.end) == (
        "6", 28477797, 33448354,
    )
    assert (grch38.chromosome, grch38.start, grch38.end) == (
        "6", 28510120, 33480577,
    )


def test_enabled_mhc_assessment_requires_every_configured_genome_build():
    with pytest.raises(ConfigurationError, match="missing: GRCh38"):
        load_configuration(cli_overrides={
            "modules.qc_summary.rules.remove_mhc": True,
            "modules.qc_summary.rules.mhc_regions": {
                "GRCh37": {
                    "chromosome": "6",
                    "start": 28477797,
                    "end": 33448354,
                },
            },
        })


def test_run_writes_reports_but_never_writes_a_filtered_vcf(tmp_path):
    raw_vcf = tmp_path / "study_<unsafe&>_GRCh37_merged.vcf.gz"
    raw_vcf.write_bytes(b"placeholder")
    document = _qc_config(
        remove_palindromic=True,
        remove_mhc=True,
    ).model_dump()
    document["output_layout"].update({
        "metric_report": "reports/{dataset_id}/{build}/metrics.tsv",
        "rule_report": "reports/{dataset_id}/{build}/rules.tsv",
        "assessment_json": "reports/{dataset_id}/{build}/assessment.json",
        "summary_csv": "reports/{dataset_id}/{build}/summary.csv",
        "html_report": "reports/{dataset_id}/{build}/report.html",
        "temporary_table_prefix": "reports/{dataset_id}/{build}/.working_",
    })
    configuration = QCSummaryConfig.model_validate(document)
    logger = _RecordLogger()

    def fake_extract(**kwargs):
        _assessment_table(Path(kwargs["table_path"]))
        return str(kwargs["table_path"])

    with patch(
        "postgwas.modules.qc_summary.assessment.extract_vcf_assessment_table",
        side_effect=fake_extract,
    ), patch(
        "postgwas.modules.qc_summary.assessment.read_vcf_header",
        return_value=(
            Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
        ).read_text(encoding="utf-8"),
    ):
        result = run_vcf_qc_assessment(
            vcf_path=raw_vcf,
            output_directory=tmp_path,
            dataset_id="study",
            external_af_name="EUR",
            configuration=configuration,
            bcftools_bin="configured-bcftools",
            **_vcf_header_contract(),
            threads=1,
            provenance_headers=_provenance_headers(),
            logger=logger,
        )

    assert Path(result["reports"]["summary"]).is_file()
    assert Path(result["reports"]["rules"]).is_file()
    assert Path(result["reports"]["json"]).is_file()
    assert Path(result["reports"]["csv"]).is_file()
    assert Path(result["reports"]["html"]).is_file()
    assert Path(result["reports"]["summary"]) == (
        tmp_path / "reports" / "study" / "GRCh37" / "metrics.tsv"
    )
    assert not list(tmp_path.rglob("*filtered*.vcf*"))
    assert not list(tmp_path.rglob(".working_*.tsv"))
    saved = json.loads(Path(result["reports"]["json"]).read_text(encoding="utf-8"))
    assert saved["qc_passed"]["num_records"] == 1
    assert "extra_rule_matches" not in saved
    assert "rule_match_total" not in saved
    assert "overlapping_rule_match_total" not in saved
    assert saved["vcf_header_validation"]["status"] == "passed"
    assert saved["vcf_header_validation"]["required_field_count"] == 6
    assert saved["vcf_provenance"]["status"] == "available"
    assert saved["vcf_provenance"]["values"]["info_source"] == (
        "external user-provided reference"
    )
    assert not {
        "input_genome_build", "output_genome_build", "liftover",
    } & saved["vcf_provenance"]["values"].keys()
    assert not {
        "input_genome_build", "output_genome_build", "liftover",
    } & set(saved["vcf_provenance"]["expected_fields"])
    assert saved["aggregation"]["aggregation_strategy"] == "two_pass_streaming"
    assert saved["aggregation"]["requested_threads"] == 1
    assert saved["aggregation"]["polars_thread_pool_size"] == 1
    assert saved["aggregation"]["thread_budget_enforced"] is True
    aggregation_records = [
        kwargs for args, kwargs in logger.records
        if args == ("PARAM", "qc_summary_aggregation")
    ]
    assert aggregation_records == [saved["aggregation"]]
    result_records = [
        kwargs for args, kwargs in logger.records
        if args == ("RESULT", "vcf_qc_assessment")
    ]
    assert len(result_records) == 1
    assert result_records[0]["sample_size_reference_quantile"] == 0.9
    assert result_records[0]["raw_sample_size_below_minimum_threshold"] == 7
    assert result_records[0][
        "qc_passed_sample_size_below_minimum_threshold"
    ] == 1
    summary = Path(result["reports"]["summary"]).read_text(encoding="utf-8")
    assert "vcf_header_validation\tstatus\tpassed" in summary
    assert "vcf_header_validation\tgenome_build\tGRCh37" in summary
    assert "vcf_header_validation\tgenome_build_source\tvcf_header" in summary
    assert "vcf_header_validation\trequired_field_count\t6" in summary
    assert "raw\teffective_sample_size_mean\t212.5" in summary
    assert "qc_passed\teffective_sample_size_mean\t100.0" in summary
    assert saved["raw"][
        "effective_sample_size_reference_quantile_value"
    ] == pytest.approx(370.0)
    assert "raw\teffective_sample_size_reference_quantile_value\t" in summary
    assert "raw\teffective_sample_size_below_minimum_threshold\t7" in summary
    assert (
        "assessment\tsample_size_reference_quantile\t0.9" in summary
    )
    assert "aggregation\taggregation_strategy\ttwo_pass_streaming" in summary
    assert "aggregation\ttemporary_table_scans\t2" in summary
    assert "extra_rule_matches" not in summary
    assert "rule_match_total" not in summary
    assert "overlapping_rule_match_total" not in summary
    assert "input_genome_build" not in summary
    assert "output_genome_build" not in summary
    assert "liftover" not in summary
    with Path(result["reports"]["csv"]).open(
        encoding="utf-8", newline="",
    ) as handle:
        reader = csv.DictReader(handle)
        csv_rows = list(reader)
    assert tuple(reader.fieldnames or ()) == QC_SUMMARY_COLUMNS
    assert csv_rows[0]["record_type"] == "overall"
    assert csv_rows[0]["raw_variants"] == "8"
    assert csv_rows[0]["qc_passed_variants"] == "1"
    assert csv_rows[0]["excluded_variants"] == "7"
    assert csv_rows[0]["variants_unique_to_one_rule"] == "6"
    assert csv_rows[0]["variants_matching_multiple_rules"] == "1"
    assert csv_rows[0]["accounting_balanced"] == "True"
    assert csv_rows[0]["requested_threads"] == "1"
    assert any(
        row["record_type"] == "rule"
        and row["key"] == "imputation_quality"
        and row["variants_matching_raw"] == "2"
        and row["variants_unique_to_rule"] == "1"
        and row["variants_overlapping_other_rules"] == "1"
        for row in csv_rows
    )
    assert any(
        row["record_type"] == "vcf_provenance"
        and row["key"] == "info_interpretation"
        and "not a study-measured" in row["provenance_value"]
        for row in csv_rows
    )
    assert not {
        row["key"]
        for row in csv_rows
        if row["record_type"] == "vcf_provenance"
    } & {"input_genome_build", "output_genome_build", "liftover"}
    html = Path(result["reports"]["html"]).read_text(encoding="utf-8")
    assert html.startswith("<!doctype html>")
    assert "GWAS-VCF quality-control report" in html
    assert "Report-only assessment:" in html
    assert "What these results mean" in html
    assert "Input versus virtual QC-passed metrics" in html
    assert "Scientific field provenance" in html
    assert "Validated input summary-statistics GWAS-VCF" in html
    assert "Input VCF" in html
    assert "Genome build" in html
    assert "Declared contigs" in html
    assert "Total variants" in html
    assert html.index("Validated input summary-statistics GWAS-VCF") < html.index(
        "What these results mean"
    )
    assert "Original summary-statistics build" not in html
    assert "Harmonised GWAS-VCF build" not in html
    assert "Genome-build conversion" not in html
    assert "Applied: GRCh38 to GRCh37" not in html
    assert "Build source" not in html
    assert "Input genome build" not in html
    assert "Output genome build" not in html
    assert "They The build-history" not in html
    assert "Imputation quality score" in html
    assert "Configured quality score" not in html
    assert "External reference proxy; not a study-measured" in html
    assert "Unique to this rule" in html
    assert "Also matched another rule" in html
    assert "FORMAT/SI below 0.7" in html
    assert "FORMAT/SI above 1.05" in html
    assert "Inactive QC rules" in html
    assert "No minimum −log10(P) threshold is configured." in html
    assert "8 input = 1 virtual QC-passed + 7 excluded" in html
    assert "Rule failures beyond the first" not in html
    assert "Additional overlapping matches" not in html
    assert "All rule matches" not in html
    assert "study_&lt;unsafe&amp;&gt;_GRCh37_merged.vcf.gz" in html
    assert "study_<unsafe&>_GRCh37_merged.vcf.gz" not in html
    assert saved["reports"] == result["reports"]


def test_qc_assessment_reports_four_truthful_data_and_report_stages(tmp_path):
    raw_vcf = tmp_path / "study_GRCh37_merged.vcf.gz"
    raw_vcf.write_bytes(b"placeholder")
    configuration = _qc_config()
    logger = _RecordLogger()
    stream = StringIO()
    progress = StageProgress(
        "GWAS-VCF QC assessment progress",
        enabled=True,
        outcome_label_width=42,
        console=Console(file=stream, force_terminal=False, width=140),
    )

    def fake_extract(**kwargs):
        _assessment_table(Path(kwargs["table_path"]))
        return str(kwargs["table_path"])

    with patch(
        "postgwas.modules.qc_summary.assessment.extract_vcf_assessment_table",
        side_effect=fake_extract,
    ), patch(
        "postgwas.modules.qc_summary.assessment.read_vcf_header",
        return_value=(
            Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
        ).read_text(encoding="utf-8"),
    ):
        result = run_vcf_qc_assessment(
            vcf_path=raw_vcf,
            output_directory=tmp_path / "output",
            dataset_id="study",
            external_af_name="EUR",
            configuration=configuration,
            bcftools_bin="configured-bcftools",
            **_vcf_header_contract(),
            threads=1,
            logger=logger,
            stage_progress=progress,
        )

    rendered = stream.getvalue()
    assert "GWAS-VCF QC assessment progress" in rendered
    assert (
        "Completed 1/4 · Convert required GWAS-VCF fields to a temporary TSV"
        in rendered
    )
    assert (
        "Completed 2/4 · Evaluate the configured QC rules and calculate metrics"
        in rendered
    )
    assert "Completed 3/4 · Generate and validate the QC reports" in rendered
    assert "Completed 4/4 · Publish the validated QC report set" in rendered
    assert "All 4 stages completed" in rendered
    stage_statuses = [
        values for marker, values in logger.records
        if marker == ("STATUS", "qc_summary_stage")
    ]
    assert [values["status"] for values in stage_statuses] == [
        "COMPLETED", "COMPLETED", "COMPLETED", "COMPLETED",
    ]
    assert result["accounting_balanced"] is True


def test_resolved_rule_contract_contains_only_validated_policy_metadata():
    contract = resolve_qc_rule_contract(
        _qc_config(), "GRCh37", "EUR",
    )

    assert contract["active_rule_count"] == len(contract["rules"])
    assert [rule["key"] for rule in contract["rules"]] == [
        "allele_frequency",
        "imputation_quality",
        "external_af_concordance",
        "variant_type",
    ]
    assert contract["field_labels"]["external_info_af"] == "INFO/EUR"
    assert all("failure" not in rule for rule in contract["rules"])
    assert [rule["display_group"] for rule in contract["rules"]] == [
        "Allele frequency · FORMAT/AF",
        "Imputation quality · FORMAT/SI",
        "Study/reference frequency concordance · INFO/AF and INFO/EUR",
        "Variant type · REF/ALT",
    ]
    assert [
        (item["action"], item["label"])
        for item in contract["decision_plan"]
    ] == [
        ("EXCLUDE", "Missing FORMAT/AF"),
        ("EXCLUDE", "Minor allele frequency outside 0.01 ≤ AF ≤ 0.99"),
        ("EXCLUDE", "Missing FORMAT/SI"),
        ("EXCLUDE", "Imputation quality below FORMAT/SI 0.7"),
        ("EXCLUDE", "Imputation quality above FORMAT/SI 1.05"),
        ("EXCLUDE", "Missing INFO/AF"),
        ("EXCLUDE", "Missing INFO/EUR"),
        ("EXCLUDE", "External AF absolute difference |AF − EUR| > 0.2"),
        ("EXCLUDE", "Indels and other non-SNP variants"),
    ]
    assert [rule["key"] for rule in contract["inactive_rules"]] == [
        "significance", "palindromic_variants", "mhc_region",
    ]


def test_resolved_rule_plan_marks_retained_missing_values_as_keep():
    contract = resolve_qc_rule_contract(
        _qc_config(
            missing_pvalue_action="keep",
            missing_af_action="keep",
            missing_info_action="keep",
        ),
        "GRCh37",
        "EUR",
    )

    missing_actions = {
        item["label"]: item["action"]
        for item in contract["decision_plan"]
        if item["label"].startswith("Missing ")
    }
    assert missing_actions == {
        "Missing FORMAT/AF": "KEEP",
        "Missing FORMAT/SI": "KEEP",
        "Missing INFO/AF": "KEEP",
        "Missing INFO/EUR": "KEEP",
    }


def test_qc_report_rendering_failure_publishes_no_partial_report_set(tmp_path):
    raw_vcf = tmp_path / "study_GRCh37_merged.vcf.gz"
    raw_vcf.write_bytes(b"placeholder")
    configuration = _qc_config()
    stream = StringIO()
    progress = StageProgress(
        "GWAS-VCF QC assessment progress",
        enabled=True,
        outcome_label_width=42,
        console=Console(file=stream, force_terminal=False, width=140),
    )

    def fake_extract(**kwargs):
        _assessment_table(Path(kwargs["table_path"]))
        return str(kwargs["table_path"])

    with patch(
        "postgwas.modules.qc_summary.assessment.extract_vcf_assessment_table",
        side_effect=fake_extract,
    ), patch(
        "postgwas.modules.qc_summary.assessment.read_vcf_header",
        return_value=(
            Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
        ).read_text(encoding="utf-8"),
    ), patch(
        "postgwas.modules.qc_summary.assessment.render_qc_summary_html",
        side_effect=RuntimeError("simulated HTML failure"),
    ), pytest.raises(RuntimeError, match="simulated HTML failure"):
        run_vcf_qc_assessment(
            vcf_path=raw_vcf,
            output_directory=tmp_path / "output",
            dataset_id="study",
            external_af_name="EUR",
            configuration=configuration,
            bcftools_bin="configured-bcftools",
            **_vcf_header_contract(),
            threads=1,
            stage_progress=progress,
        )

    rendered = stream.getvalue()
    assert "Failed 3/4 · Generate and validate the QC reports" in rendered
    assert "Completed 4/4" not in rendered
    assert "All 4 stages completed" not in rendered
    assert not list((tmp_path / "output").rglob("*.tsv"))
    assert not list((tmp_path / "output").rglob("*.json"))
    assert not list((tmp_path / "output").rglob("*.csv"))
    assert not list((tmp_path / "output").rglob("*.html"))


def test_direct_qc_uses_format_af_without_genotype_derived_metrics(tmp_path):
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("bcftools is not installed")
    vcf = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--dataset-id", "study",
        "--output-directory", str(tmp_path),
        "--threads", "1",
        "--hide-screen",
    ])

    assessment = run_qc_summary_direct(args)
    assert assessment["raw"]["num_records"] == 2
    assert assessment["raw"]["format_af_missing"] == 0
    assert assessment["vcf_header_validation"]["status"] == "passed"
    assert assessment["vcf_header_validation"]["genome_build"] == "GRCh37"
    assert assessment["vcf_header_validation"]["genome_build_source"] == (
        "vcf_header"
    )
    assert assessment["vcf_provenance"]["status"] == "available"
    assert assessment["vcf_provenance"]["values"]["info_source"] == (
        "external user-provided reference"
    )
    assert assessment["aggregation"]["aggregation_strategy"] == "two_pass_streaming"
    assert assessment["aggregation"]["requested_threads"] == 1
    assert assessment["aggregation"]["polars_thread_pool_size"] == 1
    assert assessment["aggregation"]["thread_budget_enforced"] is True
    assert "singletons" not in assessment["raw"]
    report = Path(assessment["reports"]["summary"]).read_text(encoding="utf-8")
    assert "singleton" not in report.lower()
    assert "af=0" not in report.lower()
    log = (
        tmp_path / "qc_summary" / "study_GRCh37_qc.log"
    ).read_text(encoding="utf-8")
    assert "qc_summary_runtime" in log
    assert "genome_build=GRCh37" in log
    assert "genome_build_source=vcf_header" in log
    assert "aggregation_strategy=two_pass_streaming" in log
    assert "temporary_table_scans=2" in log
    assert "polars_version=" in log
    assert "requested_threads=1" in log
    assert "polars_thread_pool_size=1" in log
    assert "thread_budget_enforced=true" in log
    assert "sample_size_reference_quantile=0.9" in log
    assert "raw_sample_size_below_minimum_threshold=0" in log
    assert "qc_passed_sample_size_below_minimum_threshold=0" in log
    assert "qc_summary_stage number=1 total=6" in log
    assert "qc_summary_stage number=6 total=6" in log
    assert "vcf_scientific_provenance status=available" in log
    assert "Input VCF validation" in log
    assert "PASSED · 6/6 required INFO/FORMAT fields declared" in log
    assert "QC assessment plan · conditions evaluated independently" in log
    assert "1. EXCLUDE · Missing FORMAT/AF" in log
    assert "GWAS-VCF quality-control summary" in log
    assert "QC results by rule" in log
    assert "Exact virtual-subset accounting" in log
    assert "Saved reports" in log
    assert "Final virtual QC outcome" in log
    assert "Detailed QC audit" in log
    assert "1. Validated input summary-statistics GWAS-VCF" in log
    assert "Input VCF" in log
    assert "Genome build" in log
    assert "Declared contigs" in log
    assert "Total variants" in log
    assert "2. Profile input GWAS-VCF" in log
    assert "Original summary-statistics build" not in log
    assert "Harmonised GWAS-VCF build" not in log
    assert "Genome-build conversion" not in log
    assert "input_genome_build=" not in log
    assert "output_genome_build=" not in log
    assert "liftover=" not in log
    assert "qc_summary_rule number=2 key=imputation_quality" in log
    assert "variants_unique_to_rule=" in log
    assert "qc_summary_inactive_rule key=significance" in log
    assert "qc_summary_decision_plan number=1 action=EXCLUDE" in log
    assert "qc_summary_decision_plan number=5 action=EXCLUDE" in log
    assert "Imputation quality above FORMAT/SI 1.05" in log
    assert "4. Final virtual QC assessment" in log
    assert "the input GWAS-VCF is unchanged; no QC-filtered VCF is created" in log
    assert "Detailed HTML report" in log
    assert log.count("STATUS   qc_summary_stage") == 6
    assert log.count("status=COMPLETED") >= 6
    completion = yaml.safe_load(
        (
            tmp_path / "qc_summary" / "study_GRCh37_qc_complete.yaml"
        ).read_text(encoding="utf-8")
    )
    assert set(completion["outputs"]) == {
        "summary", "rules", "json", "csv", "html",
    }
    assert completion["metrics"]["aggregation_strategy"] == "two_pass_streaming"
    assert completion["metrics"]["temporary_table_scans"] == 2
    assert completion["metrics"]["polars_version"] == pl.__version__
    assert completion["metrics"]["requested_threads"] == 1
    assert completion["metrics"]["polars_thread_pool_size"] == 1
    assert completion["metrics"]["thread_budget_enforced"] is True
    assert completion["metrics"]["bcftools_version"]
    assert completion["metrics"]["scientific_provenance_status"] == "available"
    assert completion["metrics"]["variants_unique_to_one_rule"] == 1
    assert completion["metrics"]["variants_matching_multiple_rules"] == 0
    assert completion["metrics"]["sample_size_reference_quantile"] == 0.9
    assert completion["metrics"]["raw_low_neff_variants"] == 0
    assert completion["metrics"]["raw_low_neff_fraction"] == 0.0
    assert completion["metrics"]["qc_passed_low_neff_variants"] == 0
    assert completion["metrics"]["qc_passed_low_neff_fraction"] == 0.0
    assert not Path(str(vcf) + ".stats").exists()

    with patch(
        "postgwas.modules.qc_summary.service.run_qc_assessment",
        side_effect=AssertionError("a validated resume must not rescan the VCF"),
    ):
        resumed = run_qc_summary_direct(args)
    assert resumed["raw"]["num_records"] == 2
    resumed_log = (
        tmp_path / "qc_summary" / "study_GRCh37_qc.log"
    ).read_text(encoding="utf-8")
    assert "qc_summary_stage number=1 total=4" in resumed_log
    assert "qc_summary_stage number=4 total=4" in resumed_log
    assert "Validate the completed QC report set" in resumed_log
    assert "Load and verify the completed QC assessment" in resumed_log
    assert "reason=provenance_validated_complete_outputs" in resumed_log


def test_direct_qc_interleaves_six_real_stages_with_stage_summaries(
    tmp_path, capsys,
):
    if shutil.which("bcftools") is None:
        pytest.skip("bcftools is not installed")
    vcf = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--dataset-id", "study",
        "--output-directory", str(tmp_path),
        "--threads", "1",
        "--overwrite",
    ])

    run_qc_summary_direct(args)

    output = " ".join(capsys.readouterr().out.split())
    milestones = (
        "Completed 1/6 · Validate the input summary-statistics GWAS-VCF",
        "Input VCF validation",
        "Completed 2/6 · Resolve and validate the configured QC rules",
        "QC assessment plan · conditions evaluated independently",
        "Completed 3/6 · Convert required GWAS-VCF fields to a temporary TSV",
        "Completed 4/6 · Evaluate the configured QC rules and calculate metrics",
        "Completed 5/6 · Generate and validate the QC reports",
        "Completed 6/6 · Publish the validated QC report set",
        "GWAS-VCF quality-control summary",
        "QC results by rule",
        "Exact virtual-subset accounting",
        "Saved reports",
        "Final virtual QC outcome",
    )
    positions = [output.index(milestone) for milestone in milestones]
    assert positions == sorted(positions)
    assert "Total variants" not in output[
        output.index("Input VCF validation"):output.index("QC assessment plan")
    ]
    assert "All 6 stages completed" in output
    assert output.count("Input VCF validation") == 1
    assert output.count("QC assessment plan") == 1
    assert output.count("QC results by rule") == 1
    assert "1. EXCLUDE · Missing FORMAT/AF" in output
    assert "5. EXCLUDE · Imputation quality above FORMAT/SI 1.05" in output
    assert "9. EXCLUDE · Indels and other non-SNP variants" in output
    assert "REMOVE ·" not in output


def test_direct_qc_missing_declaration_logs_failure_without_outputs(tmp_path):
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("bcftools is not installed")
    vcf = _vcf_without_declarations(tmp_path / "missing_eur.vcf", "INFO/EUR")
    output = tmp_path / "output"
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--dataset-id", "study",
        "--output-directory", str(output),
        "--hide-screen",
    ])

    with pytest.raises(VcfAssessmentError, match="INFO/EUR"):
        run_qc_summary_direct(args)

    qc_directory = output / "qc_summary"
    log = (qc_directory / "study_GRCh37_qc.log").read_text(encoding="utf-8")
    assert "vcf_qc_header_contract status=FAILED" in log
    assert "missing_fields=INFO/EUR" in log
    assert "qc_summary_run status=FAILED" in log
    assert "status=COMPLETED" not in log
    assert not (qc_directory / "study_GRCh37_vcf_qc_metrics.tsv").exists()
    assert not (qc_directory / "study_GRCh37_vcf_qc_rule_results.tsv").exists()
    assert not (qc_directory / "study_GRCh37_qc_assessment.json").exists()
    assert not (qc_directory / "study_GRCh37_qc_summary.csv").exists()
    assert not (output / "reports" / "study_GRCh37_qc_report.html").exists()
    assert not (qc_directory / "study_GRCh37_qc_complete.yaml").exists()


@pytest.mark.parametrize(
    ("replacement", "message"),
    [
        ("", "exactly one ##genome_build=<build>"),
        (
            "##genome_build=GRCh37\n##genome_build=GRCh38",
            "found 2",
        ),
        ("##genome_build=GRCh36", "unsupported genome build"),
    ],
)
def test_direct_qc_rejects_unusable_build_metadata_before_reports(
    tmp_path,
    replacement,
    message,
):
    if shutil.which("bcftools") is None:
        pytest.skip("bcftools is not installed")
    source = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    vcf = tmp_path / "invalid_build.vcf"
    vcf.write_text(
        source.read_text(encoding="utf-8").replace(
            "##genome_build=GRCh37",
            replacement,
            1,
        ),
        encoding="utf-8",
    )
    output = tmp_path / "output"
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--dataset-id", "study",
        "--output-directory", str(output),
        "--hide-screen",
    ])

    with pytest.raises(VcfAssessmentError, match=message):
        run_qc_summary_direct(args)

    preflight = output / "qc_summary" / "study_qc_preflight.log"
    assert preflight.is_file()
    assert "cannot infer a supported genome build" in preflight.read_text(
        encoding="utf-8"
    )
    assert "status=FAILED" in preflight.read_text(encoding="utf-8")
    assert not list((output / "qc_summary").glob("study_GRCh*_qc_assessment.json"))
    assert not (output / "reports").exists()


def test_direct_qc_uses_the_grch38_build_declared_by_the_vcf(tmp_path):
    if shutil.which("bcftools") is None:
        pytest.skip("bcftools is not installed")
    source = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    vcf = tmp_path / "study.vcf"
    vcf.write_text(
        source.read_text(encoding="utf-8").replace("GRCh37", "GRCh38"),
        encoding="utf-8",
    )
    output = tmp_path / "output"
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--dataset-id", "study",
        "--output-directory", str(output),
        "--threads", "1",
        "--hide-screen",
    ])

    assessment = run_qc_summary_direct(args)

    assert assessment["genome_build"] == "GRCh38"
    assert assessment["vcf_header_validation"]["genome_build_source"] == (
        "vcf_header"
    )
    assert Path(assessment["reports"]["json"]).name == (
        "study_GRCh38_qc_assessment.json"
    )
    assert (output / "qc_summary" / "study_GRCh38_qc.log").is_file()


def test_header_preflight_runs_before_resume_and_preserves_completed_outputs(tmp_path):
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("bcftools is not installed")
    source = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    vcf = tmp_path / "study.vcf"
    vcf.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    output = tmp_path / "output"
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--dataset-id", "study",
        "--output-directory", str(output),
        "--hide-screen",
    ])
    completed = run_qc_summary_direct(args)
    assessment_path = Path(completed["reports"]["json"])
    completed_text = assessment_path.read_text(encoding="utf-8")

    _vcf_without_declarations(vcf, "INFO/EUR")
    with patch(
        "postgwas.modules.qc_summary.service.resolve_completion_resume",
        side_effect=AssertionError("resume must not be considered before preflight"),
    ), pytest.raises(VcfAssessmentError, match="INFO/EUR"):
        run_qc_summary_direct(args)

    assert assessment_path.read_text(encoding="utf-8") == completed_text


def test_direct_run_config_is_canonical_and_cli_only_overrides_explicit_values(
    tmp_path,
):
    vcf = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    run_config = tmp_path / "qc.yaml"
    run_config.write_text(
        "reference_af_column: AFR\n"
        "rules:\n"
        "  maximum_af_difference: 0.15\n",
        encoding="utf-8",
    )
    args = build_parser().parse_args([
        "--run-config", str(run_config),
        "--vcf", str(vcf),
        "--dataset-id", "study",
        "--output-directory", str(tmp_path / "out"),
    ])

    configuration = resolve_qc_summary_configuration(args)

    assert configuration.modules.qc_summary.reference_af_column == "AFR"
    assert configuration.modules.qc_summary.rules.maximum_af_difference == 0.15
    assert configuration.modules.qc_summary.inputs.vcf == vcf


def test_qc_cli_overrides_the_shared_policy_and_low_neff_settings(tmp_path):
    run_config = tmp_path / "qc.yaml"
    run_config.write_text(
        "rules:\n"
        "  include_indels: false\n"
        "  remove_palindromic: false\n"
        "  remove_mhc: false\n",
        encoding="utf-8",
    )
    args = build_parser().parse_args([
        "--run-config", str(run_config),
        "--minimum-neglog10-p", "7.30103",
        "--minimum-maf", "0.02",
        "--reference-af-column", "AFR",
        "--maximum-af-difference", "0.1",
        "--minimum-info", "0.8",
        "--maximum-info", "2",
        "--missing-pvalue-action", "keep",
        "--missing-af-action", "keep",
        "--missing-info-action", "keep",
        "--include-indels",
        "--remove-palindromic",
        "--palindromic-af-lower", "0.35",
        "--palindromic-af-upper", "0.65",
        "--remove-mhc",
        "--sample-size-reference-quantile", "0.8",
        "--sample-size-minimum-fraction", "0.5",
    ])

    module = resolve_qc_summary_configuration(args).modules.qc_summary

    assert module.reference_af_column == "AFR"
    assert module.rules.minimum_neglog10_p == pytest.approx(7.30103)
    assert module.rules.maf_min == pytest.approx(0.02)
    assert module.rules.maximum_af_difference == pytest.approx(0.1)
    assert module.rules.info_min == pytest.approx(0.8)
    assert module.rules.info_max == pytest.approx(2.0)
    assert module.rules.missing_pvalue_action == "keep"
    assert module.rules.missing_af_action == "keep"
    assert module.rules.missing_info_action == "keep"
    assert module.rules.include_indels is True
    assert module.rules.remove_palindromic is True
    assert module.rules.palindromic_lower == pytest.approx(0.35)
    assert module.rules.palindromic_upper == pytest.approx(0.65)
    assert module.rules.remove_mhc is True
    assert module.rules.sample_size_reference_quantile == pytest.approx(0.8)
    assert (
        module.rules.sample_size_minimum_fraction_of_reference
        == pytest.approx(0.5)
    )


def test_omitted_qc_presence_flags_preserve_enabled_yaml_values(tmp_path):
    run_config = tmp_path / "qc.yaml"
    run_config.write_text(
        "rules:\n"
        "  include_indels: true\n"
        "  remove_palindromic: true\n"
        "  remove_mhc: true\n",
        encoding="utf-8",
    )

    module = resolve_qc_summary_configuration(
        build_parser().parse_args(["--run-config", str(run_config)])
    ).modules.qc_summary

    assert module.rules.include_indels is True
    assert module.rules.remove_palindromic is True
    assert module.rules.remove_mhc is True


def test_qc_cli_omits_unset_policy_overrides():
    args = build_parser().parse_args([])

    for destination in (
        "minimum_neglog10_p",
        "minimum_maf",
        "missing_pvalue_action",
        "missing_af_action",
        "minimum_info",
        "maximum_info",
        "missing_info_action",
        "include_indels",
        "remove_palindromic",
        "remove_mhc",
        "sample_size_reference_quantile",
        "sample_size_minimum_fraction",
    ):
        assert not hasattr(args, destination)


def test_qc_cli_uses_positive_presence_flags_and_infers_runtime_resources():
    parser = build_parser()
    help_text = parser.format_help()

    for option in (
        "--include-indels",
        "--remove-palindromic",
        "--remove-mhc",
    ):
        assert option in help_text
        assert "--no-" + option.removeprefix("--") not in help_text
    assert "--genome-build" not in help_text
    assert "--bcftools" not in help_text
    assert "target_build" not in QCSummaryConfig.model_fields

    with pytest.raises(SystemExit):
        parser.parse_args(["--genome-build", "GRCh37"])
    with pytest.raises(SystemExit):
        parser.parse_args(["--bcftools", "/custom/bcftools"])
    with pytest.raises(SystemExit):
        parser.parse_args(["--remove-mhc", "false"])


def test_qc_cli_rejects_info_maximum_above_two():
    args = build_parser().parse_args(["--maximum-info", "2.0001"])

    with pytest.raises(ConfigurationError, match="less than or equal to 2"):
        resolve_qc_summary_configuration(args)


def test_screen_report_has_raw_rule_and_combined_final_sections(tmp_path):
    assessment = _assess(_assessment_table(tmp_path / "raw.tsv"))
    assessment["raw_vcf"] = str(tmp_path / "study_GRCh37_merged.vcf.gz")
    assessment["vcf_header_validation"] = {
        "status": "passed",
        "genome_build": "GRCh37",
        "genome_build_source": "vcf_header",
        "declared_contig_count": 23,
        "required_field_count": 6,
        "required_fields": [
            "INFO/AF", "INFO/EUR", "FORMAT/AF", "FORMAT/SI", "FORMAT/LP",
            "FORMAT/NEF",
        ],
        "declared_info_field_count": 2,
        "declared_format_field_count": 4,
    }
    assessment["reports"] = {
        "summary": str(tmp_path / "study_qc_assessment.tsv"),
        "rules": str(tmp_path / "study_qc_filter_rules.tsv"),
        "json": str(tmp_path / "study_qc_assessment.json"),
        "csv": str(tmp_path / "study_qc_summary.csv"),
        "html": str(tmp_path / "study_qc_report.html"),
    }
    lines = harmonisation_qc_summary_lines(
        {
            "total_variant_infile": 10,
            "total_variant_read": 10,
            "total_variant_removed_missing_values": 2,
            "total_variant_remaining_for_harmonisation": 10,
            "total_variant_removed_palindromic_ambiguous": 1,
            "total_variant_removed_reference_unmatched": 1,
            "total_variant_removed_reference_ambiguous": 0,
            "total_variant_removed_chromosome_harmonisation": 2,
            "total_variant_passed_chromosome_harmonisation": 8,
            "total_variant_in_vcf_input": 8,
        },
        assessment,
    )
    output = "\n".join(lines)
    compact_output = " ".join(output.split())

    assert "1. Before VCF creation" in output
    assert "2. Validated raw merged GWAS-VCF" in output
    assert "Input VCF" in output
    assert "VCF header field contract" in output
    assert "Genome build" in output
    assert "Declared contigs" in output
    assert "Total variants" in output
    assert "GRCh37 (inferred from the VCF header)" not in compact_output
    assert "PASSED · 6 required INFO/FORMAT fields declared" in compact_output
    assert "INFO/AF, INFO/EUR, FORMAT/AF" in compact_output
    assert "3. Profile raw merged GWAS-VCF" in output
    assert "2 temporary-table scans using" in compact_output
    assert "Polars %s" % pl.__version__ in compact_output
    assert "Effective sample-size distribution" in output
    assert "Upper-tail diagnostic count" in output
    upper_tail_line = next(
        line for line in lines if "Upper-tail diagnostic count" in line
    )
    assert "🔬" in upper_tail_line
    assert "❗" not in upper_tail_line
    assert "Low-Neff reference (raw q=0.900)" in output
    assert "Below low-Neff threshold" in output
    low_neff_line = next(
        line for line in lines if "Below low-Neff threshold" in line
    )
    assert "❗" in low_neff_line
    assert "4. QC conditions assessed on the input GWAS-VCF" in output
    assert "Rule 1 · Study allele-frequency range" in output
    assert "Imputation quality score" in output
    assert "5. Final virtual QC assessment" in output
    assert "Frequency evidence" in output
    assert "Unique to this rule" in output
    assert "Also matched another rule" in output
    assert "Inactive QC conditions" in output
    assert "No minimum −log10(P) threshold is configured." in output
    assert "FORMAT/NEF upper-tail diagnostic count" in output
    assert "FORMAT/NEF values below low-Neff threshold" in output
    assert "no QC-filtered VCF is created" in compact_output
    assert "tested against every active condition" in compact_output
    assert "all 6 active conditions are applied together" in compact_output
    assert "additional overlapping matches" not in compact_output
    assert "8 input = 1 virtual QC-passed + 7 excluded" in compact_output
    assert "Passed initial input validation" in output
    assert "Missing read-stage mandatory values" in output
    assert "Palindromic ambiguous removed" in output
    assert "Reference unmatched removed" in output
    assert "Reference ambiguous removed" in output
    assert "All chromosome-stage removals" in output
    assert "Passed chromosome harmonisation" in output
    assert "10 initial = 8 sent to GWAS-to-VCF + 2 removed" in compact_output
    assert "Summary CSV" in output and "study_qc_summary.csv" in output
    assert "Detailed HTML report" in output and "study_qc_report.html" in output

    # Every field at the same hierarchy depth uses one colon column. Remove
    # the warning symbol's variation selector before measuring terminal text.
    aligned_fields = []
    for indentation in (6, 8, 12):
        fields = [
            line.replace("\ufe0f", "")
            for line in lines
            if len(line) - len(line.lstrip(" ")) == indentation and " : " in line
        ]
        assert fields
        assert len({line.index(" : ") for line in fields}) == 1
        aligned_fields.extend(fields)
    assert len({line.index(" : ") for line in aligned_fields}) == 1

    assert next(line for line in lines if "Frequency and imputation fields" in line).startswith(
        " " * 10
    )
    assert next(line for line in lines if "Rule 1 ·" in line).startswith(" " * 10)
    rule_detail = next(
        line for line in lines
        if "Missing FORMAT/AF" in line and "of input" in line
    )
    assert len(rule_detail) - len(rule_detail.lstrip(" ")) == 12
    assert next(line for line in lines if "Final virtual subset metrics" in line).startswith(
        " " * 10
    )


def test_standalone_screen_report_uses_the_shared_qc_sections(tmp_path):
    assessment = _assess(_assessment_table(tmp_path / "raw.tsv"))
    assessment["raw_vcf"] = str(tmp_path / "study_GRCh37_merged.vcf.gz")
    assessment["reports"] = {}

    lines = qc_summary_lines(assessment, label_width=42)
    output = "\n".join(lines)

    assert "Before VCF creation" not in output
    assert "1. Validated input summary-statistics GWAS-VCF" in output
    assert "Input VCF" in output
    assert "Total variants" in output
    assert "2. Profile input GWAS-VCF" in output
    assert "3. QC conditions assessed on the input GWAS-VCF" in output
    assert "4. Final virtual QC assessment" in output
    aligned = [
        line.replace("\ufe0f", "") for line in lines if " : " in line
    ]
    assert len(aligned) > 10
    assert len({cell_len(line.split(" : ", 1)[0]) for line in aligned}) == 1


def test_standalone_screen_summary_is_concise_and_uses_shared_cards(tmp_path):
    assessment = _assess(_assessment_table(tmp_path / "raw.tsv"))
    assessment["raw_vcf"] = str(tmp_path / "study_GRCh37_merged.vcf.gz")
    assessment["dataset_id"] = "study"
    assessment["genome_build"] = "GRCh37"
    assessment["vcf_header_validation"] = {
        "status": "passed",
        "genome_build": "GRCh37",
        "declared_contig_count": 23,
        "required_field_count": 6,
        "required_fields": [
            "INFO/AF", "INFO/EUR", "FORMAT/AF", "FORMAT/SI", "FORMAT/LP",
            "FORMAT/NEF",
        ],
        "declared_info_field_count": 2,
        "declared_format_field_count": 4,
    }
    assessment["vcf_provenance"] = {
        "status": "available",
        "expected_fields": ["info_source", "info_interpretation"],
        "missing_fields": [],
        "values": {
            "info_source": "external user-provided reference",
            "info_interpretation": (
                "External reference proxy; not a study-measured score"
            ),
        },
    }
    assessment["reports"] = {
        "summary": str(tmp_path / "study_metrics.tsv"),
        "rules": str(tmp_path / "study_rules.tsv"),
        "json": str(tmp_path / "study_assessment.json"),
        "csv": str(tmp_path / "study_summary.csv"),
        "html": str(tmp_path / "study_report.html"),
    }

    lines = qc_summary_screen_lines(assessment, label_width=42)
    output = " ".join("\n".join(lines).split())

    assert "GWAS-VCF quality-control summary" in output
    assert "Input VCF validation" in output
    assert "study_GRCh37_merged.vcf.gz" in output
    assert "Genome build" in output
    assert "GRCh37" in output
    assert "Declared contigs" in output
    assert "23" in output
    assert "Total variants" in output
    assert "8" in output
    assert "PASSED · 6/6 required INFO/FORMAT fields declared" in output
    assert "Field provenance" in output
    assert "FORMAT/AF for MAF; INFO/AF versus INFO/EUR" in output
    assert "Imputation quality score" in output
    assert "QC assessment plan · conditions evaluated independently" in output
    assert "EXCLUDE affects only the reported virtual QC-passed subset" in output
    applied_plan = (
        "1. EXCLUDE · Missing FORMAT/AF",
        "2. EXCLUDE · Minor allele frequency outside 0.01 ≤ AF ≤ 0.99",
        "3. EXCLUDE · Missing FORMAT/SI",
        "4. EXCLUDE · Imputation quality below FORMAT/SI 0.7",
        "5. EXCLUDE · Imputation quality above FORMAT/SI 1.05",
        "6. EXCLUDE · Missing INFO/AF",
        "7. EXCLUDE · Missing INFO/EUR",
        "8. EXCLUDE · External AF absolute difference |AF − EUR| > 0.2",
        "9. EXCLUDE · Indels and other non-SNP variants",
        (
            "10. EXCLUDE · Palindromic SNPs with ambiguous frequency "
            "(0.4 ≤ AF ≤ 0.6)"
        ),
        (
            "11. EXCLUDE · Variants in the MHC region "
            "(6:28,477,797-33,448,354)"
        ),
    )
    positions = [output.index(item) for item in applied_plan]
    assert positions == sorted(positions)
    assert "REMOVE ·" not in output
    assert "QC results by rule" in output
    result_groups = (
        "Allele frequency · FORMAT/AF",
        "Imputation quality · FORMAT/SI",
        "Study/reference frequency concordance · INFO/AF and INFO/EUR",
        "Variant type · REF/ALT",
        "Palindromic allele ambiguity · REF/ALT and FORMAT/AF",
        "Genomic region · CHROM/POS",
    )
    assert all(output.count(group) == 1 for group in result_groups)
    assert [output.index(group) for group in result_groups] == sorted(
        output.index(group) for group in result_groups
    )
    assert "Study allele-frequency range" in output
    assert "2 variants failed this rule (25.00% of input)" in output
    assert "Reference-frequency concordance" in output
    assert "1 variant failed this rule" in output
    assert "Trigger counts: Missing FORMAT/AF (0)" in output
    assert "FORMAT/SI above 1.05 (0)" in output
    assert "Missing INFO/AF (0)" in output
    assert "Missing INFO/EUR (1)" in output
    assert "Absolute frequency difference above 0.2 (1)" in output
    assert "Inactive QC rules" in output
    assert "Association significance" in output
    assert "Exact virtual-subset accounting" in output
    assert "Variants failing two or more QC rules : 1" in output
    assert "Rule failures beyond the first" not in output
    assert "1 passed + 7 excluded = 8 input" in output
    assert "Final virtual QC outcome" in output
    assert "Virtual QC-passed subset : 1" in output
    assert "Variants retained : 12.50%" in output
    assert "8 usable · 0 missing/invalid · 7 below 246.67" in output
    assert "unchanged · no filtered VCF created" in output
    assert "Saved reports" in output
    assert "study_summary.csv" in output
    assert "study_report.html" in output
    assert "study_metrics.tsv" in output
    assert "study_rules.tsv" in output
    assert "study_assessment.json" in output

    top_level_headings = (
        "Input VCF validation",
        "QC assessment plan",
        "GWAS-VCF quality-control summary",
    )
    for heading in top_level_headings:
        line = next(
            line for line in lines if heading in line and " : " not in line
        )
        assert len(line) - len(line.lstrip(" ")) == 4
    for heading in (
        "QC results by rule",
        "Exact virtual-subset accounting",
        "Saved reports",
        "Final virtual QC outcome",
    ):
        line = next(line for line in lines if heading in line)
        assert len(line) - len(line.lstrip(" ")) == 6
    detail_lines = [
        line.replace("\ufe0f", "")
        for line in lines
        if " : " in line
    ]
    assert detail_lines
    assert len({cell_len(line.split(" : ", 1)[0]) for line in detail_lines}) == 1
    label_line = next(line for line in lines if "Effective sample size" in line)
    assert " : " in label_line


def test_qc_screen_uses_the_global_theme_and_neutral_values(tmp_path):
    assessment = _assess(_assessment_table(tmp_path / "raw.tsv"))
    assessment["raw_vcf"] = str(tmp_path / "study_GRCh37_merged.vcf.gz")
    assessment["reports"] = {}
    block = "\n".join(qc_summary_screen_lines(assessment))

    styled = style_screen_block(block)
    styled_fragments = [
        (styled.plain[span.start:span.end], str(span.style))
        for span in styled.spans
    ]

    assert styled.plain == block
    for heading in ("GWAS-VCF quality-control summary", "QC results by rule", "Saved reports"):
        assert any(
            heading in fragment and style == terminal_style("analysis")
            for fragment, style in styled_fragments
        )
    for neutral_label in (
        "Dataset",
        "Genome build",
        "Allele frequency · FORMAT/AF",
        "Variants failing two or more QC rules",
    ):
        assert not any(
            neutral_label in fragment and style != "default"
            for fragment, style in styled_fragments
        )
    assert any(
        "Study allele-frequency range" in fragment
        and style == terminal_style("warning")
        for fragment, style in styled_fragments
    )
    assert any(
        "Variants retained" in fragment and style == "bold green"
        for fragment, style in styled_fragments
    )
    value = str(assessment["raw"]["num_records"])
    assert not any(
        value == fragment.strip()
        for fragment, _style in styled_fragments
    )


def test_concise_screen_rule_list_follows_the_resolved_policy(tmp_path):
    configuration = _qc_config(
        minimum_neglog10_p=7.30103,
        include_indels=True,
        remove_palindromic=False,
        remove_mhc=False,
    )
    assessment = _assess(
        _assessment_table(tmp_path / "raw.tsv"), configuration,
    )

    output = " ".join("\n".join(qc_summary_screen_lines(assessment)).split())

    assert assessment["active_rule_count"] == 4
    assert "QC assessment plan" in output
    assert "P-value evidence below FORMAT/LP 7.30103" in output
    assert "QC results by rule" in output
    assert "Association significance" in output
    assert "Variant-type scope" not in output
    assert "SNP-only scope" in output
    assert "Frequency-ambiguous palindromic variants" in output
    assert "MHC-region scope" in output


def test_direct_summary_splits_concise_screen_from_detailed_log(
    tmp_path, capsys,
):
    assessment = _assess(_assessment_table(tmp_path / "raw.tsv"))
    assessment["raw_vcf"] = str(tmp_path / "study_GRCh37_merged.vcf.gz")
    assessment["reports"] = {}
    messages = []

    class SummaryLogger:
        def log(self, *args, **kwargs):
            messages.append((args, kwargs))

    configuration = SimpleNamespace(logging=SimpleNamespace(
        terminal_label_width=42,
        show_screen=True,
    ))

    _emit_summary(assessment, configuration, SummaryLogger())

    screen = capsys.readouterr().out
    assert len(messages) == 1
    assert messages[0][0][0] == "RESULT"
    logged = messages[0][0][1]
    assert messages[0][1] == {"wrap": False, "screen": False}
    assert "\x1b[" not in logged
    for expected in (
        "GWAS-VCF quality-control summary",
        "Input VCF validation",
        "QC assessment plan",
        "QC results by rule",
        "Exact virtual-subset accounting",
        "Saved reports",
        "Final virtual QC outcome",
    ):
        assert expected in screen
    for unexpected in (
        "1. Validated input summary-statistics GWAS-VCF",
        "2. Profile input GWAS-VCF",
        "3. QC conditions assessed on the input GWAS-VCF",
        "4. Final virtual QC assessment",
    ):
        assert unexpected not in screen
    for expected in (
        "GWAS-VCF quality-control summary",
        "QC results by rule",
        "Exact virtual-subset accounting",
        "Saved reports",
        "Final virtual QC outcome",
        "Detailed QC audit",
        "1. Validated input summary-statistics GWAS-VCF",
        "2. Profile input GWAS-VCF",
        "3. QC conditions assessed on the input GWAS-VCF",
        "4. Final virtual QC assessment",
    ):
        assert expected in logged
    for screen_only in (
        "Input VCF validation", "QC assessment plan",
    ):
        assert screen_only not in logged
    for removed_metric in (
        "Rule failures beyond the first",
        "Rule failure matches",
        "Additional overlapping matches",
    ):
        assert removed_metric not in screen
        assert removed_metric not in logged


def test_hidden_screen_summary_is_still_emitted_for_the_durable_recorder(
    tmp_path, capsys, monkeypatch,
):
    assessment = _assess(_assessment_table(tmp_path / "raw.tsv"))
    assessment["raw_vcf"] = str(tmp_path / "study_GRCh37_merged.vcf.gz")

    class SummaryLogger:
        def log(self, *args, **kwargs):
            return None

    configuration = SimpleNamespace(logging=SimpleNamespace(
        terminal_label_width=42,
        show_screen=False,
    ))
    monkeypatch.setattr(
        "postgwas.modules.qc_summary.service.screen_recording_active",
        lambda: True,
    )

    _emit_summary(assessment, configuration, SummaryLogger())

    screen_stream = capsys.readouterr().out
    assert "GWAS-VCF quality-control summary" in screen_stream
    assert "QC assessment plan" in screen_stream
    assert "QC results by rule" in screen_stream


def test_final_takeaways_reuse_existing_counts_and_show_precise_af_percentages(tmp_path):
    assessment = _assess(_assessment_table(tmp_path / "raw.tsv"))
    lines = harmonisation_qc_takeaway_lines(
        {
            "total_variant_infile": 10,
            "total_variant_read": 10,
            "total_variant_in_vcf_input": 8,
            "total_variant_with_invalid_beta_se": 0,
        },
        assessment,
        genome_build_info={
            "inferred_build": "GRCh37",
            "matches": {"GRCh37": 99},
            "testable_variants": 100,
            "input_variants": 100,
            "reference_marker_counts": {"GRCh37": 120},
            "matched_reference_marker_counts": {"GRCh37": 90},
        },
        study_decisions={
            "effect_type": "odds_ratio",
            "effect_type_detected": "odds_ratio",
            "effect_type_source": "detector",
            "pvalue_type": "raw",
            "pvalue_type_detected": "raw",
            "pvalue_type_source": "detector",
            "frequency_type": "effect_allele_frequency",
            "strand": "forward",
            "strand_consensus": {"dominant_fraction": 1.0},
        },
        population_frequency={
            "closest_population": "EUR",
            "selected_population_check": {
                "normalized_population": "EUR",
                "status": "match",
            },
            "external_file_checks": [],
        },
        dataset_summary={
            "completed": ["1"],
            "failed": [],
            "chromosomes": {"1": {"status": "ok"}},
            "chromosome_summaries": {
                "1": {
                    "stage_qc": {
                        "eaf_qc": {
                            "strand_orientation": {
                                "status": "success",
                                "initial_variants": 10,
                                "final_variants": 8,
                                "reference_file": "/reference/panel_chr1.tsv.gz",
                                "reference_population_column": "EUR",
                                "study_strand_consensus": "forward",
                                "actions": {
                                    "forward": 7,
                                    "forward_swapped": 1,
                                },
                                "reference_unmatched": 1,
                                "palindromic_ambiguous": 1,
                                "reference_ambiguous": 0,
                            }
                        }
                    }
                }
            },
            "reconciliation": {
                "balanced": True,
                "complete": True,
                "unprocessed_failed_chromosome_rows": 0,
            },
            "rejected_variants_rows": 2,
            "rejected_variants_file": str(tmp_path / "rejected.tsv.gz"),
        },
        merge_summary={
            "merge_status": "OK",
            "required_merge_failures": [],
        },
        primary_outputs={
            "GRCh37": "input.vcf.gz",
            "GRCh38": "target.vcf.gz",
            "gwas2vcf": "raw.vcf.gz",
        },
        final_status="OK",
        reference_source="1000G",
        reference_population="EUR",
    )
    output = " ".join("\n".join(lines).split())

    assert "Final QC takeaways" in output
    assert "10 / 10 variants" in output
    assert "8 variants used → 8 variants in final VCF" in output
    assert "6 SNPs · 2 indels / other variants" in output
    assert "1 / 6 SNPs passed all 6 active rules (16.67%)" in output
    assert "5 SNPs failed ≥1 active rule (83.33%)" in output
    assert "6 / 7 concordant (85.7143%)" in output
    assert "1 mismatched (14.2857%)" in output
    assert (
        "0 missing study AF · 1 missing reference AF · "
        "|AF difference| ≤ 0.2"
    ) in output
    assert "2 variants removed with reasons recorded" in output
    assert "GRCh37 · automatically detected" in output
    assert "99 / 100 input variants matched (99.00% input-wide)" in output
    assert "90 / 120 relevant reference markers matched (75.00% coverage)" in output
    assert "1000G · EUR · 1 chromosome reference file" in output
    assert "8 total · forward 7 · forward-swapped 1" in output
    assert "10 evaluated = 8 retained + 2 removed" in output
    assert "odds_ratio · automatically detected" in output
    assert "raw · automatically detected" in output
    assert "7 merged-VCF variants failed ≥1 of 6 active rules" in output
    assert (
        "Every parsed variant reached a terminal bucket and all chromosomes "
        "completed"
    ) in output
    assert "1 / 1 completed" in output
    assert "Final merged VCF" in output
    assert "Indexed and validated" in output
    assert "3 primary outputs completed · dataset status OK" in output
    assert "1 low Neff below" in output

    card_headings = [
        "Variant flow",
        "Reference alignment",
        "Allele frequency",
        "Statistical quality",
        "Audit trail",
        "VCF integrity",
    ]
    for heading in card_headings:
        line = next(line for line in lines if heading in line)
        assert len(line) - len(line.lstrip(" ")) == 8
    assert lines.count("") == len(card_headings)
    detail_lines = [
        line.replace("\ufe0f", "")
        for line in lines
        if len(line) - len(line.lstrip(" ")) == 12 and " : " in line
    ]
    assert detail_lines
    assert len({line.index(" : ") for line in detail_lines}) == 1
    assert next(line for line in lines if "Input / read" in line).startswith(" " * 12)


def test_upper_neff_diagnostic_does_not_set_statistical_warning(tmp_path):
    assessment = _assess(
        _sample_size_table(tmp_path / "upper_neff.tsv", ["100"] * 10)
    )
    assessment["qc_passed"][
        "effective_sample_size_above_outlier_threshold"
    ] = 1
    lines = harmonisation_qc_takeaway_lines(
        {"total_variant_with_invalid_beta_se": 0},
        assessment,
        genome_build_info={},
        study_decisions={"effect_type": "beta", "pvalue_type": "raw"},
        population_frequency={},
        dataset_summary={},
        merge_summary={},
        primary_outputs={},
        final_status="OK",
    )

    statistical_quality = next(
        line for line in lines if "Statistical quality" in line
    )
    assert "✅" in statistical_quality
    assert "1 informational upper-tail diagnostic" in " ".join(lines)


def test_strand_summary_combines_completed_chromosomes_without_rescanning_rows():
    def chromosome(initial, final, actions, unmatched, palindromic, ambiguous, name):
        return {
            "stage_qc": {
                "eaf_qc": {
                    "strand_orientation": {
                        "status": "success",
                        "initial_variants": initial,
                        "final_variants": final,
                        "reference_file": "/reference/%s.tsv.gz" % name,
                        "reference_population_column": "EUR",
                        "study_strand_consensus": "forward",
                        "actions": actions,
                        "reference_unmatched": unmatched,
                        "palindromic_ambiguous": palindromic,
                        "reference_ambiguous": ambiguous,
                    }
                }
            }
        }

    summary = summarise_strand_orientation({
        "completed": ["1", "2"],
        "failed": [],
        "chromosomes": {"1": {"status": "ok"}, "2": {"status": "ok"}},
        "chromosome_summaries": {
            "1": chromosome(
                10, 8,
                {"forward": 5, "forward_swapped": 2, "reverse_complement": 1},
                1, 0, 1, "panel_chr1",
            ),
            "2": chromosome(
                20, 17,
                {
                    "forward": 15,
                    "forward_swapped": 1,
                    "reverse_complement_swapped": 1,
                },
                1, 1, 1, "panel_chr2",
            ),
        },
    })

    assert summary == {
        "status": "success",
        "chromosomes_expected": 2,
        "chromosomes_completed": 2,
        "chromosomes_summarized": 2,
        "complete_chromosome_coverage": True,
        "metadata_complete": True,
        "reference_population": "EUR",
        "reference_file_count": 2,
        "reference_files": [
            "/reference/panel_chr1.tsv.gz",
            "/reference/panel_chr2.tsv.gz",
        ],
        "study_consensus": "forward",
        "variants_evaluated": 30,
        "variants_retained": 25,
        "variants_matched": 25,
        "reference_unmatched_detected": 2,
        "reference_unmatched_retained": 0,
        "forward": 20,
        "forward_swapped": 3,
        "reverse_complement": 1,
        "reverse_complement_swapped": 1,
        "reference_unmatched": 2,
        "palindromic_frequency_conflict": 0,
        "palindromic_orientation_unavailable": 0,
        "palindromic_frequency_discordant": 0,
        "palindromic_ambiguous": 1,
        "reference_ambiguous": 2,
        "removed_total": 5,
        "accounting_balanced": True,
    }


def test_strand_summary_accounts_for_opt_in_reference_unmatched_retention():
    orientation = {
        "status": "success",
        "initial_variants": 3,
        "final_variants": 2,
        "reference_file": "/reference/panel_chr1.tsv.gz",
        "reference_population_column": "EUR",
        "study_strand_consensus": "mixed",
        "actions": {
            "forward": 1,
            "reference_unmatched_retained": 1,
        },
        "reference_unmatched_detected": 2,
        "reference_unmatched_retained": 1,
        "reference_unmatched": 0,
        "palindromic_orientation_unavailable": 1,
        "palindromic_ambiguous": 0,
        "reference_ambiguous": 0,
    }

    summary = summarise_strand_orientation({
        "completed": ["1"],
        "failed": [],
        "chromosomes": {"1": {"status": "ok"}},
        "chromosome_summaries": {
            "1": {"stage_qc": {"eaf_qc": {
                "strand_orientation": orientation,
            }}},
        },
    })

    assert summary["status"] == "success"
    assert summary["variants_evaluated"] == 3
    assert summary["variants_retained"] == 2
    assert summary["variants_matched"] == 1
    assert summary["reference_unmatched_detected"] == 2
    assert summary["reference_unmatched_retained"] == 1
    assert summary["reference_unmatched"] == 0
    assert summary["palindromic_orientation_unavailable"] == 1
    assert summary["removed_total"] == 1
    assert summary["accounting_balanced"] is True


def test_strand_summary_rejects_legacy_disabled_and_malformed_evidence():
    legacy_disabled = {
        "stage_qc": {
            "eaf_qc": {
                "strand_orientation": {
                    "status": "disabled",
                    "initial_variants": 4,
                    "final_variants": 4,
                    "actions": {"disabled": 4},
                }
            }
        }
    }
    result = summarise_strand_orientation({
        "completed": ["1"],
        "failed": [],
        "chromosomes": {"1": {"status": "ok"}},
        "chromosome_summaries": {"1": legacy_disabled},
    })
    assert result["status"] == "inconsistent"
    assert result["removed_total"] is None
    assert result["accounting_balanced"] is False

    malformed = {
        "stage_qc": {
            "eaf_qc": {
                "strand_orientation": {
                    "status": "success",
                    "initial_variants": 4.5,
                    "final_variants": 4,
                    "actions": {"forward": 4},
                    "reference_unmatched": 0,
                    "palindromic_ambiguous": 0,
                    "reference_ambiguous": 0,
                }
            }
        }
    }
    result = summarise_strand_orientation({
        "completed": ["1"],
        "failed": ["2"],
        "chromosomes": {
            "1": {"status": "ok"},
            "2": {"status": "failed"},
        },
        "chromosome_summaries": {"1": malformed},
    })
    assert result["status"] == "partial"
    assert result["complete_chromosome_coverage"] is False
    assert result["variants_evaluated"] is None
    assert result["accounting_balanced"] is False
