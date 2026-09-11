"""Regression tests for v2 harmonisation sample-sheet generation."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

from postgwas.config import load_module_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.modules.harmonisation.input_validation import validate_header
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.sample_sheet import (
    HarmonisationSampleSheetRow,
    load_harmonisation_sample_sheet,
)
from postgwas.modules.harmonisation.sample_sheet_generator import (
    _generation_summary,
    generate_sample_sheet,
    main,
)


def _summary_statistics(
    path: Path,
    *,
    frequency: str | None = "EAF",
    effect: str = "OR",
    p_value: str = "P",
    sample_columns: tuple[str, ...] = ("N_CONTROL", "N_CASE"),
) -> Path:
    columns = ["CHR", "BP", "SNP", "A1", "A2"]
    values = ["1", "100", "rs1", "A", "G"]
    if frequency is not None:
        columns.append(frequency)
        values.append("0.2")
    columns.extend((effect, "SE", p_value, *sample_columns, "INFO"))
    values.extend((
        "1.1" if effect == "OR" else "0.1",
        "0.03",
        "8" if p_value == "MLOG10P" else "1e-8",
        *("1000" for _ in sample_columns),
        "0.99",
    ))
    path.write_text(
        "\t".join(columns) + "\n" + "\t".join(values) + "\n",
        encoding="utf-8",
    )
    return path


def _read_row(path: Path) -> tuple[list[str], dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or ()), next(reader)


def test_generator_writes_exact_v2_schema_and_validates_output(tmp_path):
    _summary_statistics(tmp_path / "study.meta")
    output = tmp_path / "manifest.csv"

    result = generate_sample_sheet(tmp_path, output)
    fields, generated = _read_row(output)
    loaded = load_harmonisation_sample_sheet(output)[0]

    assert fields == list(HarmonisationSampleSheetRow.model_fields)
    assert result.dataset_count == 1
    assert result.candidate_count == 1
    assert result.rejected_files == ()
    assert result.rejection_report_file.read_text(encoding="utf-8") == (
        "input_file\treason\n"
    )
    assert result.requires_completion is False
    assert loaded.dataset_id == "study"
    assert loaded.effect_column == "OR"
    assert loaded.trait_type == "auto"
    assert loaded.effect_type == "auto"
    assert loaded.p_value_type == "auto"
    assert loaded.delimiter == "tab"
    assert generated["input_file"] == str((tmp_path / "study.meta").resolve())
    assert "sumstat_file" not in generated
    assert "resource_folder" not in generated
    assert "output_folder" not in generated
    assert "liftover" not in generated


def test_generator_removes_complete_vcf_gzip_suffix_from_dataset_id(tmp_path):
    source = tmp_path / "study.vcf.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(
            "CHR\tBP\tSNP\tA1\tA2\tEAF\tOR\tSE\tP\tN_CONTROL\tN_CASE\tINFO\n"
            "1\t100\trs1\tA\tG\t0.2\t1.1\t0.03\t1e-8\t1000\t1000\t0.99\n"
        )
    output = tmp_path / "manifest.csv"

    generate_sample_sheet(tmp_path, output)
    _, generated = _read_row(output)

    assert generated["dataset_id"] == "study"


def test_generator_keeps_trait_effect_and_p_value_types_auto(tmp_path):
    _summary_statistics(
        tmp_path / "quantitative.tsv",
        effect="BETA",
        p_value="MLOG10P",
        sample_columns=("N",),
    )
    output = tmp_path / "manifest.csv"

    generate_sample_sheet(tmp_path, output)
    loaded = load_harmonisation_sample_sheet(output)[0]

    assert loaded.effect_column == "BETA"
    assert loaded.trait_type == "auto"
    assert loaded.effect_type == "auto"
    assert loaded.p_value_column == "MLOG10P"
    assert loaded.p_value_type == "auto"
    assert loaded.control_count_column == "N"
    assert loaded.case_count_column is None
    assert loaded.trait_type == "auto"


@pytest.mark.parametrize("chromosome_header", ("#chrom", "#chr"))
def test_generator_maps_alt_aligned_sequencing_headers(
    tmp_path,
    chromosome_header,
):
    source = tmp_path / "summary_stats.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(
            "\t".join((
                chromosome_header, "pos", "rsid", "ref", "alt", "neg_log_pvalue",
                "beta", "stderr_beta", "alt_allele_freq",
            ))
            + "\n"
            + "\t".join((
                "1", "13668", "rs2691328", "G", "A", "1.22268",
                "-0.560559", "0.297911", "0.00599738",
            ))
            + "\n"
        )
    output = tmp_path / "manifest.csv"

    result = generate_sample_sheet(tmp_path, output)
    fields, generated = _read_row(output)

    assert fields == list(HarmonisationSampleSheetRow.model_fields)
    assert result.requires_completion is True
    assert generated["chromosome_column"] == chromosome_header
    assert generated["effect_allele_column"] == "alt"
    assert generated["other_allele_column"] == "ref"
    assert generated["effect_allele_frequency_column"] == "alt_allele_freq"
    assert generated["standard_error_column"] == "stderr_beta"
    assert generated["p_value_column"] == "neg_log_pvalue"
    assert generated["effect_type"] == "auto"
    assert generated["p_value_type"] == "auto"
    assert generated["control_count_column"] == "NA"
    assert generated["case_count_column"] == "NA"
    assert generated["control_count"] == "NA"
    assert generated["case_count"] == "NA"
    assert result.warnings == ()
    assert len(result.datasets) == 1
    assert result.datasets[0].missing_sample_size is True
    assert result.datasets[0].missing_info is True
    summary = _generation_summary(result)
    assert "Need attention    : 1" in summary
    assert "summary_stats (summary_stats.gz): sample size, imputation INFO" in summary
    assert "How to complete sample size" in summary
    assert "How to provide imputation INFO" in summary
    valid, problems = validate_header(
        {
            "delimiter": "tab",
            "chr_col": generated["chromosome_column"],
            "pos_col": generated["position_column"],
            "snp_id_col": generated["variant_id_column"],
            "ea_col": generated["effect_allele_column"],
            "oa_col": generated["other_allele_column"],
            "eaf_col": generated["effect_allele_frequency_column"],
            "beta_or_col": generated["effect_column"],
            "se_col": generated["standard_error_column"],
            "pval_col": generated["p_value_column"],
        },
        str(source),
        policies=default_policies().with_overrides({"input.delimiter": "tab"}),
    )
    assert valid, [str(problem) for problem in problems]
    with pytest.raises(ConfigurationError, match="provide control_count_column"):
        load_harmonisation_sample_sheet(output)


def test_cli_groups_attention_for_ten_input_files(tmp_path, capsys):
    for index in range(9):
        _summary_statistics(tmp_path / ("ready_%02d.tsv" % index))
    _summary_statistics(
        tmp_path / "missing_sample_size.tsv",
        sample_columns=(),
    )
    output = tmp_path / "manifest.csv"

    assert main([
        "--input-directory", str(tmp_path),
        "--output", str(output),
    ]) == 0
    terminal = capsys.readouterr().out
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 10
    assert "Input files       : 10" in terminal
    assert "Rows written      : 10" in terminal
    assert "Ready             : 9" in terminal
    assert "Need attention    : 1" in terminal
    assert "Trait type        : auto" in terminal
    assert "Effect type       : auto" in terminal
    assert "P-value type      : auto" in terminal
    assert "missing_sample_size (missing_sample_size.tsv): sample size" in terminal
    assert terminal.count("How to complete sample size") == 1


def test_generator_can_restore_strict_missing_sample_size_policy(tmp_path):
    _summary_statistics(
        tmp_path / "missing_n.tsv",
        sample_columns=(),
    )
    output = tmp_path / "manifest.csv"
    config = tmp_path / "strict.yaml"
    config.write_text(
        "sample_sheet_generator:\n  on_missing_sample_size: fail\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="no recognized control"):
        generate_sample_sheet(tmp_path, output, config_file=config)

    assert not output.exists()


def test_generator_rejects_alt_frequency_when_effect_is_not_alt(tmp_path):
    _summary_statistics(
        tmp_path / "misaligned.tsv",
        frequency="ALT_ALLELE_FREQ",
        effect="BETA",
        sample_columns=("N",),
    )
    output = tmp_path / "manifest.csv"

    with pytest.raises(
        ConfigurationError,
        match="Do not use ALT frequency as EAF unless the reported effect is aligned to ALT",
    ):
        generate_sample_sheet(tmp_path, output)

    assert not output.exists()


def test_generator_marks_maf_for_required_reference_confirmation(tmp_path):
    _summary_statistics(tmp_path / "maf.tsv", frequency="MAF")
    output = tmp_path / "manifest.csv"

    result = generate_sample_sheet(tmp_path, output)

    assert load_harmonisation_sample_sheet(output)[0].effect_allele_frequency_column == "MAF"
    assert len(result.warnings) == 1
    assert "MAF-like, not confirmed EAF" in result.warnings[0]


def test_generator_writes_draft_when_study_frequency_is_missing(tmp_path):
    _summary_statistics(tmp_path / "missing_frequency.tsv", frequency=None)
    output = tmp_path / "manifest.csv"

    result = generate_sample_sheet(tmp_path, output)
    _, generated = _read_row(output)
    summary = _generation_summary(result)

    assert result.requires_completion is True
    assert result.datasets[0].missing_eaf is True
    assert generated["effect_allele_frequency_column"] == "NA"
    assert generated["external_eaf_file"] == "NA"
    assert generated["external_eaf_column"] == "NA"
    assert "missing_frequency (missing_frequency.tsv): allele frequency" in summary
    assert "How to provide allele frequency" in summary
    assert "do not invent a fixed frequency" in summary
    with pytest.raises(ConfigurationError, match="has no allele-frequency source"):
        load_harmonisation_sample_sheet(output)


def test_missing_frequency_does_not_block_other_dataset_rows(tmp_path):
    _summary_statistics(tmp_path / "complete.tsv")
    _summary_statistics(tmp_path / "missing_frequency.tsv", frequency=None)
    output = tmp_path / "manifest.csv"

    result = generate_sample_sheet(tmp_path, output)
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert result.dataset_count == 2
    assert len(rows) == 2
    assert [row["dataset_id"] for row in rows] == ["complete", "missing_frequency"]
    assert rows[0]["effect_allele_frequency_column"] == "EAF"
    assert rows[1]["effect_allele_frequency_column"] == "NA"
    assert [dataset.missing_eaf for dataset in result.datasets] == [False, True]


def test_generator_does_not_relabel_reference_frequency_as_study_eaf(tmp_path):
    _summary_statistics(
        tmp_path / "reference.tsv",
        frequency="FREQ_EUROPEAN_1000GENOMES",
    )
    output = tmp_path / "manifest.csv"

    result = generate_sample_sheet(tmp_path, output)
    _, generated = _read_row(output)

    assert output.is_file()
    assert result.requires_completion is True
    assert result.datasets[0].missing_eaf is True
    assert generated["effect_allele_frequency_column"] == "NA"
    assert len(result.warnings) == 1
    assert "was not used as study EAF" in result.warnings[0]


def test_generator_does_not_treat_total_n_as_controls_when_cases_exist(tmp_path):
    _summary_statistics(
        tmp_path / "incomplete_binary.tsv",
        sample_columns=("N", "N_CASE"),
    )
    output = tmp_path / "manifest.csv"

    with pytest.raises(ConfigurationError, match="cannot be assumed to contain controls"):
        generate_sample_sheet(tmp_path, output)

    assert not output.exists()


def test_default_policy_writes_valid_rows_and_reports_unmapped_files(tmp_path):
    _summary_statistics(tmp_path / "valid.tsv")
    (tmp_path / "invalid.tsv").write_text("A\tB\n1\t2\n", encoding="utf-8")
    output = tmp_path / "manifest.csv"

    result = generate_sample_sheet(tmp_path, output)
    _, generated = _read_row(output)
    with result.rejection_report_file.open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        rejections = list(csv.DictReader(handle, delimiter="\t"))

    assert load_module_configuration(
        "harmonisation"
    ).sample_sheet_generator.failure_policy == "write_valid"
    assert result.candidate_count == 2
    assert result.dataset_count == 1
    assert generated["dataset_id"] == "valid"
    assert len(result.rejected_files) == 1
    assert rejections == [{
        "input_file": str((tmp_path / "invalid.tsv").resolve()),
        "reason": result.rejected_files[0].reason,
    }]
    summary = _generation_summary(result)
    assert "Rejected files    : 1" in summary
    assert "invalid.tsv:" in summary

    rerun = generate_sample_sheet(tmp_path, output)

    assert rerun.candidate_count == 2
    assert rerun.dataset_count == 1


def test_fail_all_policy_remains_available_explicitly(tmp_path):
    _summary_statistics(tmp_path / "valid.tsv")
    (tmp_path / "invalid.tsv").write_text("A\tB\n1\t2\n", encoding="utf-8")
    output = tmp_path / "manifest.csv"
    config = tmp_path / "strict.yaml"
    config.write_text(
        "sample_sheet_generator:\n  failure_policy: fail_all\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="no output was written"):
        generate_sample_sheet(tmp_path, output, config_file=config)

    assert not output.exists()
    assert not (tmp_path / "manifest.csv.rejected_files.tsv").exists()


def test_default_policy_fails_when_no_candidate_can_be_mapped(tmp_path):
    (tmp_path / "invalid.tsv").write_text("A\tB\n1\t2\n", encoding="utf-8")
    output = tmp_path / "manifest.csv"

    with pytest.raises(ConfigurationError, match="none of the 1 candidate"):
        generate_sample_sheet(tmp_path, output)

    assert not output.exists()
    assert not (tmp_path / "manifest.csv.rejected_files.tsv").exists()
