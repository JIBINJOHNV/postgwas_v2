"""Population-frequency QC on the unfiltered concatenated GWAS-VCF."""

from pathlib import Path
from unittest.mock import patch

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.population_frequency import (
    assess_population_frequency_table,
    check_external_population_filenames,
    check_selected_population,
    run_population_frequency_qc,
)
from postgwas.modules.harmonisation.service import _population_frequency_block


def _settings():
    values = (
        load_configuration()
        .modules.harmonisation.population_frequency_qc.model_dump()
    )
    values.update({
        "minimum_comparable_variants": 4,
        "minimum_correlation": -1.0,
        "minimum_correlation_gap": 0.001,
    })
    return values


def _frequency_table(path: Path) -> Path:
    path.write_text(
        "STUDY_AF\tAFR\tEAS\tEUR\tSAS\n"
        "0.05\t0.10\t0.70\t0.05\t0.10\n"
        "0.10\t0.40\t0.10\t0.11\t0.20\n"
        "0.20\t0.15\t0.60\t0.19\t0.15\n"
        "0.35\t0.70\t0.20\t0.36\t0.40\n"
        "0.50\t0.30\t0.80\t0.49\t0.45\n"
        "0.65\t0.80\t0.30\t0.66\t0.70\n"
        "0.80\t0.55\t0.90\t0.79\t0.75\n"
        "0.90\t0.95\t0.40\t0.91\t.\n",
        encoding="utf-8",
    )
    return path


def test_similarity_uses_common_complete_rows_and_reports_missingness(tmp_path):
    result = assess_population_frequency_table(
        _frequency_table(tmp_path / "frequencies.tsv"),
        _settings(),
        delimiter="\t",
        null_values=["", "."],
    )

    assert result["total_records"] == 8
    assert result["comparable_variants"] == 7
    assert result["fields"]["SAS"]["missing"] == 1
    assert result["fields"]["SAS"]["missing_fraction"] == 0.125
    assert result["fields"]["EUR"]["missing"] == 0
    assert all(
        metrics["comparable_variants"] == 7
        for metrics in result["populations"].values()
    )
    assert result["closest_population"] == "EUR"
    assert result["status"] == "complete"


def test_only_the_best_population_must_pass_the_minimum_correlation(tmp_path):
    settings = _settings()
    settings["minimum_correlation"] = 0.99

    result = assess_population_frequency_table(
        _frequency_table(tmp_path / "frequencies.tsv"),
        settings,
        delimiter="\t",
        null_values=["", "."],
    )

    assert result["closest_population"] == "EUR"


def test_external_eaf_and_info_filenames_receive_the_same_mismatch_warning():
    checks, warnings = check_external_population_filenames(
        "EUR",
        {
            "external EAF": "/references/1000G_AFR_frequency.tsv.gz",
            "external INFO": "/references/imputation_EAS_info.tsv.gz",
        },
        ["AFR", "EAS", "EUR", "SAS"],
    )

    assert [check["status"] for check in checks] == ["mismatch", "mismatch"]
    assert len(warnings) == 2
    assert checks[0]["normalized_filename"] == "1000G_AFR_FREQUENCY_TSV_GZ"
    assert "external EAF filename 1000G_AFR_frequency.tsv.gz indicates AFR" in warnings[0]
    assert "external INFO filename imputation_EAS_info.tsv.gz indicates EAS" in warnings[1]


def test_filename_population_requires_an_exact_token():
    checks, warnings = check_external_population_filenames(
        "EUR",
        {"external EAF": "/references/AFRICA_frequency.tsv.gz"},
        ["AFR", "EAS", "EUR", "SAS"],
    )

    assert checks[0]["status"] == "no_population_token"
    assert warnings == []


def test_selected_and_external_population_matches_are_rendered():
    selected, warnings = check_selected_population(
        "EUR", "eur", ["AFR", "EAS", "EUR", "SAS"],
    )
    checks, filename_warnings = check_external_population_filenames(
        "EUR",
        {"external EAF": "/references/1000G-eur-frequency.tsv.gz"},
        ["AFR", "EAS", "EUR", "SAS"],
    )
    block = _population_frequency_block({
        "status": "complete",
        "closest_population": "EUR",
        "selected_population_check": selected,
        "external_file_checks": checks,
    })

    assert warnings == []
    assert filename_warnings == []
    assert "EUR matches closest population EUR" in block
    assert "indicates EUR and matches closest population EUR" in block

    mismatch, mismatch_warnings = check_selected_population(
        "EUR", "AFR", ["AFR", "EAS", "EUR", "SAS"],
    )
    assert mismatch["status"] == "mismatch"
    assert "selected comparison-AF column is AFR" in mismatch_warnings[0]


def test_vcf_extraction_contains_only_the_five_configured_info_fields(tmp_path):
    raw_vcf = tmp_path / "study_GRCh37_merged.vcf.gz"
    raw_vcf.write_bytes(b"placeholder")
    config = load_configuration().modules.harmonisation

    def fake_extract(_vcf, table, _dataset, _columns, _bcftools, **_kwargs):
        _frequency_table(Path(table))
        return str(table)

    with patch(
        "postgwas.modules.harmonisation.population_frequency.extract_vcf_table",
        side_effect=fake_extract,
    ) as extract:
        result = run_population_frequency_qc(
            vcf_path=raw_vcf,
            output_directory=tmp_path,
            dataset_id="study",
            genome_build="GRCh37",
            settings=_settings(),
            output_layout=dict(config.output_layout.root),
            table_delimiter=config.vcf_processing.table_delimiter,
            table_null_values=config.vcf_processing.table_null_values,
            temporary_table_suffix=config.vcf_processing.temporary_table_suffix,
            io_buffer_bytes=config.vcf_processing.io_buffer_bytes,
            bcftools_bin="configured-bcftools",
            selected_population="EUR",
            external_files={},
        )

    assert extract.call_args.args[3] == {
        "STUDY_AF": "%INFO/AF",
        "AFR": "%INFO/AFR",
        "EAS": "%INFO/EAS",
        "EUR": "%INFO/EUR",
        "SAS": "%INFO/SAS",
    }
    assert Path(result["report"]).is_file()
    assert result["raw_merged_vcf"] == str(raw_vcf.resolve())
    assert result["selected_population_check"]["status"] == "match"
