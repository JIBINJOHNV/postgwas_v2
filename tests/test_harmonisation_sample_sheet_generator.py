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
    assert "Missing sample N  : 1" in summary
    assert "Missing INFO      : 1" in summary
    assert "summary_stats.gz" not in summary
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
    assert "Missing sample N  : 1" in terminal
    assert "missing_sample_size.tsv" not in terminal


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

    with pytest.raises(ConfigurationError, match="no sample sheet"):
        generate_sample_sheet(tmp_path, output, config_file=config)

    assert not output.exists()


def test_generator_retains_alt_frequency_with_alignment_review(tmp_path):
    _summary_statistics(tmp_path / "misaligned.tsv", frequency="ALT_ALLELE_FREQ")
    result = generate_sample_sheet(tmp_path, tmp_path / "manifest.csv")
    _, row = _read_row(result.output_file)
    assert row["effect_allele_frequency_column"] == "ALT_ALLELE_FREQ"
    assert result.datasets[0].needs_attention
    assert "confirm allele alignment" in result.warnings[0]


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
    assert "Missing frequency : 1" in summary
    assert "missing_frequency.tsv" not in summary
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


@pytest.mark.parametrize("frequency", ["EAF_HRC", "FREQ_EUROPEAN_1000GENOMES"])
def test_generator_retains_reference_frequency_with_provenance_review(tmp_path, frequency):
    _summary_statistics(tmp_path / "reference.tsv", frequency=frequency)
    result = generate_sample_sheet(tmp_path, tmp_path / "manifest.csv")
    _, row = _read_row(result.output_file)
    assert row["effect_allele_frequency_column"] == frequency
    assert not result.datasets[0].missing_eaf
    assert result.datasets[0].needs_attention
    assert "reference-panel frequency candidate retained" in result.warnings[0]


def test_generator_retains_cases_without_relabelling_total_n_as_controls(tmp_path):
    _summary_statistics(tmp_path / "incomplete_binary.tsv", sample_columns=("N", "N_CASE"))
    result = generate_sample_sheet(tmp_path, tmp_path / "manifest.csv")
    _, row = _read_row(result.output_file)
    assert row["case_count_column"] == "N_CASE"
    assert row["control_count_column"] == "NA"
    assert "cannot be assumed to contain controls" in result.warnings[0]
    with pytest.raises(ConfigurationError, match="provide control_count_column"):
        load_harmonisation_sample_sheet(result.output_file)


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
    assert "invalid.tsv:" not in summary

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

    with pytest.raises(ConfigurationError, match="no sample sheet was written"):
        generate_sample_sheet(tmp_path, output, config_file=config)

    assert not output.exists()
    assert (tmp_path / "manifest.csv.rejected_files.tsv").exists()


def test_default_policy_fails_when_no_candidate_can_be_mapped(tmp_path):
    (tmp_path / "invalid.tsv").write_text("A\tB\n1\t2\n", encoding="utf-8")
    output = tmp_path / "manifest.csv"

    with pytest.raises(ConfigurationError, match="None of the 1 candidate"):
        generate_sample_sheet(tmp_path, output)

    assert not output.exists()
    assert (tmp_path / "manifest.csv.rejected_files.tsv").exists()


@pytest.mark.parametrize("name", ["study?invalid", "study with spaces", "_study"])
@pytest.mark.parametrize("frequency", ["EAF", None])
def test_invalid_generated_ids_follow_per_file_rejection_policy(tmp_path, name, frequency):
    source = _summary_statistics(tmp_path / f"{name}.tsv", frequency=frequency)
    original = source.read_bytes()
    _summary_statistics(tmp_path / "valid.tsv")
    output = tmp_path / "manifest.csv"

    result = generate_sample_sheet(tmp_path, output)

    assert result.candidate_count == 2
    assert result.dataset_count == 1
    assert load_harmonisation_sample_sheet(output)[0].dataset_id == "valid"
    assert len(result.rejected_files) == 1
    assert result.rejected_files[0].input_file == source.resolve()
    assert "dataset_id" in result.rejected_files[0].reason
    assert "use only letters" in result.rejected_files[0].reason
    assert "pydantic" not in result.rejected_files[0].reason
    assert source.read_bytes() == original
    assert str(source.resolve()) in result.rejection_report_file.read_text()
    assert source.name not in _generation_summary(result)


@pytest.mark.parametrize("policy", ["write_valid", "fail_all"])
def test_invalid_generated_ids_produce_actionable_cli_failure(tmp_path, policy, capsys):
    source = _summary_statistics(tmp_path / "study?invalid.tsv", frequency=None)
    if policy == "fail_all":
        _summary_statistics(tmp_path / "valid.tsv")
    config = tmp_path / "strict.yaml"
    config.write_text(f"sample_sheet_generator:\n  failure_policy: {policy}\n")
    output = tmp_path / "manifest.csv"

    with pytest.raises(SystemExit) as failure:
        main(["--input-directory", str(tmp_path), "--output", str(output),
              "--run-config", str(config)])

    assert failure.value.code == 2
    error = capsys.readouterr().err
    assert source.name not in error
    report = (tmp_path / "manifest.csv.rejected_files.tsv").read_text()
    assert source.name in report
    assert "dataset_id" in report
    assert "Traceback" not in error
    assert "pydantic" not in error
    assert not output.exists()
    assert (tmp_path / "manifest.csv.rejected_files.tsv").exists()


@pytest.mark.parametrize("frequency", ["EAF", None])
def test_plus_in_dataset_id_is_preserved_in_generated_and_loaded_rows(tmp_path, frequency):
    from postgwas.modules.harmonisation.sample_sheet import to_harmonisation_input

    name = "EVER_SMOKER_GWAS_MA_UKB+TAG"
    source = _summary_statistics(tmp_path / f"{name}.tsv", frequency=frequency)
    original = source.read_bytes()
    output = tmp_path / "manifest.csv"
    result = generate_sample_sheet(tmp_path, output)
    _, generated = _read_row(output)
    assert result.dataset_count == 1
    assert result.rejected_files == ()
    assert generated["dataset_id"] == name
    assert source.read_bytes() == original
    row = HarmonisationSampleSheetRow.validate_generated_draft(generated)
    if frequency:
        row = load_harmonisation_sample_sheet(output)[0]
    inputs = to_harmonisation_input(
        row, resource_directory=tmp_path, output_directory=tmp_path / "results",
        output_layout={"dataset_directory": "{dataset_id}/harmonisation"},
    )
    assert inputs["gwas_outputname"] == name
    assert Path(inputs["output_folder"]) == tmp_path / "results" / name / "harmonisation"


@pytest.mark.parametrize("count,retained", [(2, False), (3, True), (6, True)])
def test_partial_mapping_screening_boundary_and_analysis_validation(tmp_path, count, retained):
    headers = ["CHR", "BP", "A1", "A2", "BETA", "P"][:count]
    (tmp_path / "partial.tsv").write_text("\t".join(headers) + "\n" + "\t".join(["1"] * count) + "\n")
    output = tmp_path / "manifest.csv"
    if not retained:
        with pytest.raises(ConfigurationError, match="draft screening"):
            generate_sample_sheet(tmp_path, output)
        assert not output.exists()
        return
    result = generate_sample_sheet(tmp_path, output)
    assert result.dataset_count == 1
    assert result.requires_completion
    with pytest.raises(ConfigurationError):
        load_harmonisation_sample_sheet(output)
    if count == 3:
        _, row = _read_row(output)
        assert row["other_allele_column"] == "NA"
        assert "other_allele_column" in result.datasets[0].missing_mappings


def test_alias_priority_is_configurable_and_reports_alternatives(tmp_path):
    source = _summary_statistics(tmp_path / "study.tsv")
    source.write_text(source.read_text().replace("OR\t", "BETA\tOR\t").replace("1.1\t", "0.1\t1.1\t"))
    result = generate_sample_sheet(tmp_path, tmp_path / "manifest.csv")
    _, row = _read_row(result.output_file)
    assert row["effect_column"] == "BETA"
    assert any("alternatives: OR" in note for note in result.warnings)
    source.write_text(source.read_text().replace("BETA\t", "B\tBETA\t").replace("0.1\t", "0.2\t0.1\t", 1))
    config = tmp_path / "priority.yaml"
    config.write_text("sample_sheet_generator:\n  column_aliases:\n    beta_effect_column: [BETA, B]\n")
    result = generate_sample_sheet(tmp_path, result.output_file, config_file=config)
    _, row = _read_row(result.output_file)
    assert row["effect_column"] == "BETA"
    assert any("alternatives: B, OR" in note for note in result.warnings)


def test_review_table_covers_ready_incomplete_frequency_and_rejected_files(tmp_path):
    _summary_statistics(tmp_path / "ready.tsv")
    missing = _summary_statistics(tmp_path / "missing.tsv", sample_columns=(), frequency=None)
    missing.write_text(missing.read_text().replace("\tINFO", "").replace("\t0.99", ""))
    _summary_statistics(tmp_path / "maf.tsv", frequency="MAF")
    _summary_statistics(tmp_path / "reference.tsv", frequency="EAF_HRC")
    (tmp_path / "unmapped.tsv").write_text("X\tY\n1\t2\n")
    result = generate_sample_sheet(tmp_path, tmp_path / "manifest.csv")
    with result.review_report_file.open() as handle:
        rows = {r["dataset_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    assert len(rows) == result.candidate_count == 5
    assert rows["ready"]["status"] == "ready"
    assert rows["missing"]["status"] == "needs_attention"
    for field in ("missing_info", "missing_sample_size", "missing_frequency"):
        assert rows["missing"][field] == "True"
    assert "MAF-like, not confirmed EAF" in rows["maf"]["review_notes"]
    assert "reference-panel" in rows["reference"]["review_notes"]
    assert rows["unmapped"]["status"] == "rejected"
    assert "screening mappings" in rows["unmapped"]["review_notes"]
    assert generate_sample_sheet(tmp_path, result.output_file).candidate_count == 5


@pytest.mark.parametrize("settings", [
    "minimum_mapped_columns: 0", "minimum_mapped_columns: 7",
    "screening_columns: [not_a_field]", "screening_columns: [effect_column, effect_column]",
    "review_report_suffix: .rejected_files.tsv",
])
def test_invalid_draft_screening_configuration_fails(tmp_path, settings):
    _summary_statistics(tmp_path / "study.tsv")
    config = tmp_path / "invalid.yaml"
    config.write_text("sample_sheet_generator:\n  " + settings + "\n")
    with pytest.raises((ConfigurationError, ValueError)):
        generate_sample_sheet(tmp_path, tmp_path / "manifest.csv", config_file=config)
    assert not (tmp_path / "manifest.csv").exists()


def test_console_summary_has_counts_without_per_file_details(tmp_path, capsys):
    _summary_statistics(tmp_path / "ready.tsv")
    _summary_statistics(tmp_path / "frequency_review.tsv", frequency="MAF", sample_columns=())
    (tmp_path / "bad.tsv").write_text("X\tY\n1\t2\n")
    output = tmp_path / "manifest.csv"
    assert main(["--input-directory", str(tmp_path), "--output", str(output)]) == 0
    terminal = capsys.readouterr().out
    for line in ("Input files       : 3", "Rows written      : 2", "Rejected files    : 1",
                 "Missing sample N  : 1", "Case/control cols : 1", "Frequency mapped  : 2",
                 "Coordinates mapped: 2", "Mapping review    : 1"):
        assert line in terminal
    assert len(terminal.splitlines()) <= 25
    for detail in ("ready.tsv", "frequency_review.tsv", "bad.tsv", "MAF-like"):
        assert detail not in terminal
    assert "MAF-like" in Path(str(output) + ".mapping_review.tsv").read_text()


def test_public_short_command_matches_direct_module(tmp_path):
    import os
    import subprocess
    import sys
    source = tmp_path / 'inputs'
    source.mkdir()
    _summary_statistics(source / 'UKB+TAG.tsv', frequency='MAF')
    outputs = [tmp_path / 'short.csv', tmp_path / 'direct.csv']
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    for command, output in zip(
        [['postgwas', 'sample_sheet_generator'],
         ['postgwas.modules.harmonisation.sample_sheet_generator']], outputs,
    ):
        result = subprocess.run(
            [sys.executable, '-m', *command, '--input-directory', str(source), '--output', str(output)],
            cwd=tmp_path, env=env, capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert 'Rows written      : 1' in result.stdout
        assert 'UKB+TAG.tsv' not in result.stdout
    assert outputs[0].read_bytes() == outputs[1].read_bytes()
    assert not (tmp_path / 'postgwas_output').exists()


def test_short_command_reports_all_missing_inputs(tmp_path):
    import os
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, '-m', 'postgwas', 'sample_sheet_generator'], cwd=tmp_path,
        env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'), capture_output=True, text=True,
    )
    assert result.returncode == 2
    assert 'Required argument not provided: --input-directory.' in result.stderr
    assert 'Required argument not provided: --output.' in result.stderr
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize('frequency', [None, 'EAF', 'MAF', 'EAF_HRC'])
@pytest.mark.parametrize('has_info', [False, True])
def test_external_defaults_fill_only_absent_mappings(tmp_path, frequency, has_info):
    source = _summary_statistics(tmp_path / 'study.tsv', frequency=frequency)
    if not has_info:
        source.write_text(source.read_text().replace('\tINFO', '').replace('\t0.99', ''))
    reference = tmp_path / 'reference.tsv.gz'
    with gzip.open(reference, 'wt') as handle:
        handle.write('CHR\tBP\tA1\tA2\tAF\tINFO\n1\t100\tA\tG\t0.2\t0.9\n')
    output = tmp_path / 'studies.csv'
    result = generate_sample_sheet(
        tmp_path, output, external_eaf_file=reference, external_eaf_column='AF',
        external_info_file=reference, external_info_column='INFO',
    )
    _, row = _read_row(output)
    assert result.candidate_count == 1  # Reference is not treated as a GWAS.
    assert row['effect_allele_frequency_column'] == (frequency or 'NA')
    assert row['external_eaf_file'] == ('NA' if frequency else str(reference))
    assert row['external_info_file'] == ('NA' if has_info else str(reference))
    status = result.datasets[0]
    assert not status.missing_eaf and not status.missing_info
    assert status.eaf_source == ('internal' if frequency else 'external')
    assert status.info_source == ('internal' if has_info else 'external')
    with result.review_report_file.open() as handle:
        review = next(csv.DictReader(handle, delimiter='\t'))
    assert review['eaf_source'] == status.eaf_source
    assert review['info_source'] == status.info_source
    assert 'Missing INFO      : 0' in _generation_summary(result)
    assert len(load_harmonisation_sample_sheet(output)) == 1
    again = generate_sample_sheet(
        tmp_path, output, external_eaf_file=reference, external_eaf_column='AF',
        external_info_file=reference, external_info_column='INFO',
    )
    assert again.candidate_count == 1


def test_external_templates_yaml_and_cli_override(tmp_path, monkeypatch):
    data = tmp_path / 'data'
    data.mkdir()
    _summary_statistics(data / 'study.tsv', frequency=None)
    for build in ('GRCh37', 'GRCh38'):
        folder = tmp_path / build
        folder.mkdir()
        (folder / 'chr1.tsv').write_text('CHR\tBP\tA1\tA2\tAF\n1\t100\tA\tG\t.2\n')
    config = tmp_path / 'run.yaml'
    config.write_text('sample_sheet_generator:\n  external_eaf_file: "{build}/chr{chromosome}.tsv"\n  external_eaf_column: OLD\n')
    monkeypatch.chdir(data)
    output = data / 'studies.csv'
    generate_sample_sheet(data, output, config_file=config, external_eaf_column='AF')
    _, row = _read_row(output)
    assert row['external_eaf_file'] == str(tmp_path / '{build}/chr{chromosome}.tsv')
    assert row['external_eaf_column'] == 'AF'
    from postgwas.modules.harmonisation.resource_paths import resolve_resource_file
    for build in ('GRCh37', 'GRCh38'):
        assert resolve_resource_file(row['external_eaf_file'], build, '1') == str(tmp_path / build / 'chr1.tsv')


def test_external_pairs_report_all_missing_options_before_outputs(tmp_path, capsys):
    with pytest.raises(SystemExit) as caught:
        main(['--input-directory', str(tmp_path), '--output', str(tmp_path / 'studies.csv'),
              '--external-eaf-file', 'ref.tsv', '--external-info-column', 'INFO'])
    assert caught.value.code == 2
    error = capsys.readouterr().err
    assert 'Required argument not provided: --external-eaf-column' in error
    assert 'Required argument not provided: --external-info-file' in error
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize('spec', ['missing.tsv', '{unsupported}/ref.tsv', '{build}/missing.tsv'])
def test_external_invalid_paths_fail_before_output(tmp_path, spec):
    _summary_statistics(tmp_path / 'study.tsv')
    output = tmp_path / 'studies.csv'
    with pytest.raises(ConfigurationError, match='Invalid argument --external-eaf-file'):
        generate_sample_sheet(tmp_path, output, external_eaf_file=tmp_path / spec,
                              external_eaf_column='AF')
    assert not output.exists()
    assert not list(tmp_path.glob('*.mapping_review.tsv'))


def test_reference_cannot_be_overwritten_by_output(tmp_path):
    _summary_statistics(tmp_path / 'study.tsv', frequency=None)
    reference = tmp_path / 'reference.csv'
    reference.write_text('CHR,BP,A1,A2,AF\n1,100,A,G,.2\n')
    before = reference.read_bytes()
    with pytest.raises(ConfigurationError, match='must not overwrite'):
        generate_sample_sheet(tmp_path, reference, external_eaf_file=reference, external_eaf_column='AF')
    assert reference.read_bytes() == before


def test_external_short_cli_with_relative_path(tmp_path):
    import os
    import subprocess
    import sys
    _summary_statistics(tmp_path / 'study.tsv', frequency=None)
    reference = tmp_path / 'reference.tsv'
    reference.write_text('CHR\tBP\tA1\tA2\tAF\n1\t100\tA\tG\t.2\n')
    output = tmp_path / 'studies.csv'
    result = subprocess.run(
        [sys.executable, '-m', 'postgwas', 'sample_sheet_generator',
         '--input-directory', '.', '--output', str(output),
         '--external-eaf-file', 'reference.tsv', '--external-eaf-column', 'AF'],
        cwd=tmp_path, capture_output=True, text=True,
        env={**os.environ, 'PYTHONPATH': str(Path(__file__).resolve().parents[1] / 'src'),
             'PYTHONDONTWRITEBYTECODE': '1'},
    )
    assert result.returncode == 0, result.stderr
    assert 'External EAF      : 1' in result.stdout
    assert 'Missing frequency : 0' in result.stdout
    _, row = _read_row(output)
    assert row['external_eaf_file'] == str(reference)
    assert row['external_info_file'] == 'NA'


@pytest.mark.parametrize('field', ['external_eaf_file', 'external_eaf_column', 'external_info_file', 'external_info_column'])
def test_blank_external_values_are_invalid(tmp_path, field):
    with pytest.raises(ConfigurationError, match='non-empty or null'):
        generate_sample_sheet(tmp_path, tmp_path / 'out.csv', **{field: '   '})


@pytest.mark.parametrize('suffix', ['.tbl', '.tbl.gz', '.tbl.bgz', '.tbl.bz2', '.tbl.xz', '.tbl.zip', '.gz', '.zip', '.TBL.GZ'])
def test_table_and_compressed_suffixes(tmp_path, suffix):
    import bz2
    import lzma
    import zipfile
    plain = _summary_statistics(tmp_path / 'fixture.tsv')
    data = plain.read_bytes()
    plain.unlink()
    source = tmp_path / ('study' + suffix)
    if suffix.lower().endswith('.zip'):
        with zipfile.ZipFile(source, 'w') as archive:
            archive.writestr('study.tbl', data)
            archive.writestr('__MACOSX/._study.tbl', 'metadata')
    else:
        opener = (gzip.open if suffix.lower().endswith(('.gz', '.bgz')) else
                  bz2.open if suffix.endswith('.bz2') else
                  lzma.open if suffix.endswith('.xz') else open)
        with opener(source, 'wb') as handle:
            handle.write(data)
    output = tmp_path / 'studies.csv'
    result = generate_sample_sheet(tmp_path, output)
    _, row = _read_row(output)
    assert result.dataset_count == 1 and not result.rejected_files
    assert row['dataset_id'] == 'study'
    assert row['input_file'] == str(source)
    assert row['effect_allele_frequency_column'] == 'EAF'
    assert row['p_value_column'] == 'P'


def test_multifile_zip_is_reported_not_silently_selected(tmp_path):
    import zipfile
    _summary_statistics(tmp_path / 'valid.tbl')
    with zipfile.ZipFile(tmp_path / 'multiple.zip', 'w') as archive:
        archive.writestr('first.tbl', 'CHR\tBP\tA1\n1\t100\tA\n')
        archive.writestr('second.tbl', 'CHR\tBP\tA1\n1\t200\tA\n')
    result = generate_sample_sheet(tmp_path, tmp_path / 'studies.csv')
    assert result.dataset_count == 1
    assert len(result.rejected_files) == 1
    assert 'exactly one data file' in result.rejected_files[0].reason


@pytest.mark.parametrize('weight_header', ['Weight', 'WEIGHT', 'weight'])
@pytest.mark.parametrize('sample_columns,expected,cases', [
    ((), 'fallback', 'NA'),
    (('N',), 'N', 'NA'),
    (('NEFF',), 'NEFF', 'NA'),
    (('N_CONTROL',), 'N_CONTROL', 'NA'),
    (('N_CASE',), 'fallback', 'N_CASE'),
    (('N_CONTROL', 'N_CASE'), 'N_CONTROL', 'N_CASE'),
    (('N', 'N_CASE'), 'fallback', 'N_CASE'),
])
def test_weight_fallback_mapping_and_priority(tmp_path, weight_header, sample_columns, expected, cases):
    _summary_statistics(tmp_path / 'study.tbl', sample_columns=(*sample_columns, weight_header))
    output = tmp_path / 'studies.csv'
    result = generate_sample_sheet(tmp_path, output)
    _, row = _read_row(output)
    assert row['control_count_column'] == (weight_header if expected == 'fallback' else expected)
    assert row['case_count_column'] == cases
    assert not result.datasets[0].missing_sample_size
    notes = result.review_report_file.read_text()
    assert ('as a fallback' in notes) == (expected == 'fallback')
    if expected == 'fallback':
        assert result.datasets[0].needs_attention
        assert 'Missing sample N  : 0' in _generation_summary(result)
    assert load_harmonisation_sample_sheet(output)[0].control_count_column == row['control_count_column']


def test_weight_satisfies_fail_on_missing_and_alias_is_configurable(tmp_path):
    _summary_statistics(tmp_path / 'study.tbl', sample_columns=('CUSTOM_WEIGHT',))
    config = tmp_path / 'run.yaml'
    config.write_text('sample_sheet_generator:\n  on_missing_sample_size: fail\n  column_aliases:\n    fallback_sample_size_column: [CUSTOM_WEIGHT]\n')
    output = tmp_path / 'studies.csv'
    generate_sample_sheet(tmp_path, output, config_file=config)
    assert _read_row(output)[1]['control_count_column'] == 'CUSTOM_WEIGHT'


_MIXED_TABLE = (
    'CHR \t BP \t A1 \t A2 \t EAF \t BETA \t SE \t P \t N \t INFO\n'
    '1\t100\tA\tG\t0.2\t0.1\t0.03\t0.01 1000 0.99\n'
    '1\t200\tC\tT\t0.3\t0.2\t0.04\t0.02   1000 0.98\n'
)


@pytest.mark.parametrize('suffix', ['.txt', '.txt.gz'])
@pytest.mark.parametrize('conversion_enabled', [True, False])
def test_space_padded_tabs_use_validated_preparation(
    tmp_path, suffix, conversion_enabled,
):
    from postgwas.modules.harmonisation.summary_statistics_io import load_summary_statistics_table

    records = [
        ['CHR', 'BP', 'SNP', 'A1', 'A2', 'EAF', 'BETA', 'SE', 'P', 'N', 'INFO'],
        ['1', '79033', 'rs2462495', 'A', 'G', '0.0012', '0.0617', '0.1017', '0.5439', '216611', '0.99'],
        ['1', '79034', 'rs2', 'C', 'T', '0.3', '-0.1127', '0.104', '1e-400', '190387', '0.98'],
    ]
    text = ''.join(' \t '.join(record) + '\n' for record in records)
    source = tmp_path / ('study' + suffix)
    opener = gzip.open if suffix.endswith('.gz') else open
    with opener(source, 'wt', encoding='utf-8') as handle:
        handle.write(text)
    original = source.read_bytes()
    config = tmp_path / 'config.yaml'
    config.write_text(
        'sample_sheet_generator:\n  convert_mixed_whitespace: '
        + str(conversion_enabled).lower() + '\n'
    )
    output = tmp_path / 'studies.csv'
    if not conversion_enabled:
        with pytest.raises(ConfigurationError, match='None of the 1 candidate'):
            generate_sample_sheet(tmp_path, output, config_file=config)
        with Path(str(output) + '.rejected_files.tsv').open() as handle:
            rejection = next(csv.DictReader(handle, delimiter='\t'))
        assert 'header' in rejection['reason'].lower()
        assert 'whitespace' in rejection['reason'].lower()
        assert source.read_bytes() == original
        assert not output.exists()
        assert not (tmp_path / 'prepared_inputs').exists()
        return
    result = generate_sample_sheet(tmp_path, output, config_file=config)
    _, row = _read_row(result.output_file)
    prepared = tmp_path / 'prepared_inputs' / (source.name + '.converted.txt')

    assert result.dataset_count == 1 and not result.rejected_files
    assert row['input_file'] == str(prepared)
    assert row['delimiter'] == 'space'
    assert row['effect_column'] == 'BETA'
    assert row['standard_error_column'] == 'SE'
    assert row['p_value_column'] == 'P'
    assert not result.datasets[0].missing_mappings
    assert result.datasets[0].prepared_input_file == prepared
    assert result.datasets[0].converted_rows == 2
    assert 'Conversions       : 1 attempted; 1 accepted; 0 rejected' in _generation_summary(result)
    assert prepared.read_text() == ''.join(' '.join(record) + '\n' for record in records)
    assert source.read_bytes() == original
    with result.review_report_file.open() as handle:
        review = next(csv.DictReader(handle, delimiter='\t'))
    assert review['conversion_status'] == 'converted'

    mapping = {
        'chr_col': row['chromosome_column'], 'pos_col': row['position_column'],
        'snp_id_col': row['variant_id_column'], 'ea_col': row['effect_allele_column'],
        'oa_col': row['other_allele_column'], 'eaf_col': row['effect_allele_frequency_column'],
        'beta_or_col': row['effect_column'], 'se_col': row['standard_error_column'],
        'pval_col': row['p_value_column'], 'ncontrol_col': row['control_count_column'],
    }
    frame, _ = load_summary_statistics_table(
        str(prepared), mapping,
        policies=default_policies().with_overrides({'input.delimiter': row['delimiter']}),
    )
    assert frame.height == 2
    assert frame[row['effect_column']].to_list() == ['0.0617', '-0.1127']
    assert frame[row['standard_error_column']].to_list() == ['0.1017', '0.104']
    assert frame[row['p_value_column']].to_list() == ['0.5439', '1e-400']


def test_padded_mixed_separators_recover_after_invalid_space_candidate(tmp_path):
    from postgwas.modules.harmonisation.summary_statistics_io import load_summary_statistics_table

    headers = ['CHR', 'BP', 'A1', 'A2', 'EAF', 'BETA', 'SE', 'P', 'N', 'INFO']
    values = [
        ['1', '79033', 'A', 'G', '0.0012', '0.0617', '0.1017', '0.5439', '216611', '0.99'],
        ['1', '79034', 'C', 'T', '0.3', '-0.1127', '0.104', '1e-400', '190387', '0.98'],
    ]
    # Tab counts differ (10/9/9), but literal spaces produce consistent 19-field
    # records with tab-only artificial header names. This is recoverable padding,
    # unlike genuinely empty or duplicate source header fields.
    source = tmp_path / 'study.tbl'
    source.write_text(
        ' \t '.join(headers) + '\n'
        + ''.join(' \t '.join(row).replace(' \t ', '  ', 1) + '\n' for row in values)
    )
    original = source.read_bytes()

    result = generate_sample_sheet(tmp_path, tmp_path / 'studies.csv')
    _, row = _read_row(result.output_file)
    prepared = tmp_path / 'prepared_inputs' / (source.name + '.converted.txt')

    assert result.dataset_count == 1 and not result.rejected_files
    assert row['input_file'] == str(prepared)
    assert row['delimiter'] == 'space'
    assert result.datasets[0].converted_rows == 2
    assert source.read_bytes() == original
    assert prepared.read_text() == ' '.join(headers) + '\n' + ''.join(' '.join(value) + '\n' for value in values)
    mapping = {
        'chr_col': row['chromosome_column'], 'pos_col': row['position_column'],
        'ea_col': row['effect_allele_column'], 'oa_col': row['other_allele_column'],
        'eaf_col': row['effect_allele_frequency_column'], 'beta_or_col': row['effect_column'],
        'se_col': row['standard_error_column'], 'pval_col': row['p_value_column'],
        'ncontrol_col': row['control_count_column'], 'imp_info_col': row['imputation_info_column'],
    }
    frame, _ = load_summary_statistics_table(
        str(prepared), mapping,
        policies=default_policies().with_overrides({'input.delimiter': row['delimiter']}),
    )
    assert frame.select(headers).rows() == [tuple(value) for value in values]


@pytest.mark.parametrize('defect', [
    'empty_header', 'duplicate_header', 'placental_width',
    'comma_empty_header', 'semicolon_empty_header',
])
def test_header_candidate_validation_keeps_malformed_sources_rejected(tmp_path, defect):
    _summary_statistics(tmp_path / 'valid.tsv')
    source = tmp_path / 'malformed.tsv'
    if defect == 'placental_width':
        headers = [
            'MarkerName', 'Allele1', 'Allele2', 'Freq1', 'FreqSE', 'MinFreq', 'MaxFreq',
            'Effect', 'StdErr', 'P-value', 'Direction', 'HetISq', 'HetChiSq', 'HetDf',
            'HetPVal', 'TOTALSAMPLESIZEchr', 'pos',
        ]
        values = [
            'rs10399793', 'T', 'C', '0.4104', '0.0174', '0.3622', '0.4507',
            '-0.0143', '0.0308', '0.6437', '??-+++', '9.9', '5.552', '5',
            '0.3523', '5093', '1', '49298',
        ]
    elif defect in ('comma_empty_header', 'semicolon_empty_header'):
        headers = ['CHR', 'BP', '', 'A1', 'A2', 'EAF', 'BETA', 'SE', 'P', 'N', 'INFO']
        values = ['1', '100', 'unlabelled', 'A', 'G', '0.2', '0.1', '0.03', '0.01', '1000', '0.99']
    else:
        headers = ['CHR', 'BP', 'A1', 'A2', 'EAF', 'BETA', 'SE', 'P', 'N', 'INFO']
        values = ['1', '100', 'A', 'G', '0.2', '0.1', '0.03', '0.01', '1000', '0.99']
        headers.insert(6, '' if defect == 'empty_header' else 'BETA-')
        values.insert(6, '0.2')
    separator = ', ' if defect == 'comma_empty_header' else '; ' if defect == 'semicolon_empty_header' else '\t'
    source.write_text(separator.join(headers) + '\n' + separator.join(values) + '\n')
    original = source.read_bytes()

    result = generate_sample_sheet(tmp_path, tmp_path / 'studies.csv')

    assert result.dataset_count == 1
    assert len(result.rejected_files) == 1
    assert result.rejected_files[0].input_file == source
    assert source.read_bytes() == original
    assert load_harmonisation_sample_sheet(result.output_file)[0].dataset_id == 'valid'
    assert not list(tmp_path.glob('prepared_inputs/*.converted.txt'))


def test_native_tab_quoted_text_and_empty_cells_are_preserved(tmp_path):
    from postgwas.modules.harmonisation.summary_statistics_io import load_summary_statistics_table

    source = _summary_statistics(tmp_path / 'study.tsv', effect='BETA', sample_columns=('N',))
    header, data = source.read_text().splitlines()
    source.write_text(
        header + '\t"Study note"\n'
        + data.replace('\t0.03\t', '\t\t') + '\t"brain age"\n'
    )
    original = source.read_bytes()
    result = generate_sample_sheet(tmp_path, tmp_path / 'studies.csv')
    _, row = _read_row(result.output_file)

    assert row['delimiter'] == 'tab'
    assert result.datasets[0].prepared_input_file is None
    frame, _ = load_summary_statistics_table(
        str(source), {'beta_or_col': row['effect_column'], 'se_col': row['standard_error_column'],
                      'pval_col': row['p_value_column']},
        policies=default_policies().with_overrides({'input.delimiter': row['delimiter']}),
    )
    assert frame['Study note'].to_list() == ['brain age']
    assert frame['BETA'].to_list() == ['0.1']
    assert frame['SE'].to_list() == [None]
    assert frame['P'].to_list() == ['1e-8']
    assert source.read_bytes() == original
    assert not (tmp_path / 'prepared_inputs').exists()


@pytest.mark.parametrize('suffix', ['.tbl', '.csv', '.txt', '.gz', '.tbl.gz', '.zip', '.tsv.bz2', '.txt.xz'])
def test_mixed_whitespace_conversion_and_downstream_reader(tmp_path, suffix):
    import bz2
    import lzma
    import zipfile
    from postgwas.modules.harmonisation.summary_statistics_io import load_summary_statistics_table
    source = tmp_path / ('study' + suffix)
    if suffix.endswith('.zip'):
        with zipfile.ZipFile(source, 'w') as archive:
            archive.writestr('study.txt', _MIXED_TABLE)
    else:
        opener = gzip.open if suffix.endswith('.gz') else bz2.open if suffix.endswith('.bz2') else lzma.open if suffix.endswith('.xz') else open
        with opener(source, 'wt') as handle:
            handle.write(_MIXED_TABLE)
    original = source.read_bytes()
    output = tmp_path / 'studies.csv'
    result = generate_sample_sheet(tmp_path, output)
    _, row = _read_row(output)
    converted = tmp_path / 'prepared_inputs' / (source.name + '.converted.txt')
    assert source.read_bytes() == original
    assert row['input_file'] == str(converted)
    assert row['dataset_id'] == 'study'
    assert row['delimiter'] == 'space'
    assert result.datasets[0].converted_rows == 2
    assert result.datasets[0].needs_attention
    assert 'Conversions       : 1 attempted; 1 accepted; 0 rejected' in _generation_summary(result)
    assert converted.read_text() == __import__('re').sub(' +', ' ', _MIXED_TABLE.replace('\t', ' '))
    mapping = {'delimiter': row['delimiter'], 'chr_col': row['chromosome_column'],
               'pos_col': row['position_column'], 'ea_col': row['effect_allele_column'],
               'oa_col': row['other_allele_column'], 'beta_or_col': row['effect_column'],
               'pval_col': row['p_value_column'], 'ncontrol_col': row['control_count_column']}
    frame, _mapping = load_summary_statistics_table(str(converted), mapping)
    assert frame.height == 2
    assert frame['BETA'].to_list() == ['0.1', '0.2']
    assert frame['N'].to_list() == ['1000', '1000']
    mtime = converted.stat().st_mtime_ns
    again = generate_sample_sheet(tmp_path, output)
    assert again.candidate_count == 1
    assert converted.stat().st_mtime_ns == mtime
    with result.review_report_file.open() as handle:
        report = next(csv.DictReader(handle, delimiter='\t'))
    assert report['input_file'] == str(source)
    assert report['prepared_input_file'] == str(converted)
    assert report['conversion_status'] == 'converted'
    assert report['converted_rows'] == '2'
    assert not list(converted.parent.glob('.preparing-*'))


@pytest.mark.parametrize('bad_line', [
    '1\t100\tA\tG\t0.2\t0.1\t0.03\t0.01 1000\n',
    '1\t100\tA\tG\t\t0.1\t0.03\t0.01 1000 0.99\n',
    '1\t100\tA\tG\t0.2\t0.1\t0.03\t0.01 "1000" 0.99\n',
    ' 1\t100\tA\tG\t0.2\t0.1\t0.03\t0.01 1000 0.99\n',
])
def test_conversion_rejects_late_ragged_or_ambiguous_rows(tmp_path, bad_line):
    _summary_statistics(tmp_path / 'valid.tsv')
    source = tmp_path / 'broken.tbl'
    # Ensure the late bad row is checked even beyond delimiter sampling.
    source.write_text(_MIXED_TABLE + _MIXED_TABLE.splitlines(True)[1] * 100 + bad_line)
    result = generate_sample_sheet(tmp_path, tmp_path / 'studies.csv')
    assert result.dataset_count == 1
    assert result.rejected_files[0].conversion_attempted
    assert 'Line ' in result.rejected_files[0].reason
    assert 'Conversions       : 1 attempted; 0 accepted; 1 rejected' in _generation_summary(result)
    assert not list((tmp_path / 'prepared_inputs').iterdir())


def test_conversion_disabled_and_no_conversion_for_mapping_failures(tmp_path):
    _summary_statistics(tmp_path / 'valid.tsv')
    (tmp_path / 'mixed.tbl').write_text(_MIXED_TABLE)
    (tmp_path / 'unrelated.tsv').write_text('X\tY\n1\t2\n')
    config = tmp_path / 'config.yaml'
    config.write_text('sample_sheet_generator:\n  convert_mixed_whitespace: false\n')
    result = generate_sample_sheet(tmp_path, tmp_path / 'studies.csv', config_file=config)
    assert len(result.rejected_files) == 2
    assert not any(item.conversion_attempted for item in result.rejected_files)
    assert not (tmp_path / 'prepared_inputs').exists()


def test_conversion_refuses_modified_output(tmp_path):
    (tmp_path / 'mixed.tbl').write_text(_MIXED_TABLE)
    output = tmp_path / 'studies.csv'
    generate_sample_sheet(tmp_path, output)
    converted = tmp_path / 'prepared_inputs/mixed.tbl.converted.txt'
    converted.write_text('user modified content\n')
    with pytest.raises(ConfigurationError, match='None of the 1'):
        generate_sample_sheet(tmp_path, output)
    assert converted.read_text() == 'user modified content\n'
    assert 'different content' in (tmp_path / 'studies.csv.rejected_files.tsv').read_text()


def test_configured_conversion_names_and_symlink_refusal(tmp_path):
    (tmp_path / 'mixed.tbl').write_text(_MIXED_TABLE)
    config = tmp_path / 'config.yaml'
    config.write_text('sample_sheet_generator:\n  prepared_input_directory: prepared\n  converted_file_suffix: .fixed.txt\n')
    output = tmp_path / 'studies.csv'
    generate_sample_sheet(tmp_path, output, config_file=config)
    assert Path(_read_row(output)[1]['input_file']).name == 'mixed.tbl.fixed.txt'
    converted = tmp_path / 'prepared/mixed.tbl.fixed.txt'
    converted.unlink()
    converted.symlink_to(tmp_path / 'mixed.tbl')
    with pytest.raises(ConfigurationError):
        generate_sample_sheet(tmp_path, output, config_file=config)
    assert converted.is_symlink()


def test_whitespace_only_rows_are_not_silently_changed(tmp_path):
    (tmp_path / 'mixed.tbl').write_text(_MIXED_TABLE + '  \t  \n')
    with pytest.raises(ConfigurationError):
        generate_sample_sheet(tmp_path, tmp_path / 'studies.csv')
    assert 'whitespace-only row' in (tmp_path / 'studies.csv.rejected_files.tsv').read_text()
    assert not list((tmp_path / 'prepared_inputs').iterdir())


@pytest.mark.parametrize('key,value', [('prepared_input_directory','../escape'), ('prepared_input_directory','/tmp/out'), ('converted_file_suffix',''), ('converted_file_suffix','/escape.txt')])
def test_preparation_paths_cannot_escape_output(tmp_path, key, value):
    import yaml
    config = tmp_path / 'config.yaml'
    config.write_text(yaml.safe_dump({'sample_sheet_generator': {key: value}}))
    with pytest.raises(ConfigurationError):
        generate_sample_sheet(tmp_path, tmp_path / 'studies.csv', config_file=config)
    assert not (tmp_path / 'studies.csv').exists()


def test_prepared_files_are_excluded_and_mapping_rejection_leaves_no_copy(tmp_path):
    _summary_statistics(tmp_path / 'valid.tsv')
    _summary_statistics(tmp_path / 'old.tbl.converted.txt')
    (tmp_path / 'unrelated.tbl').write_text('X \t Y\n1 2\n')
    result = generate_sample_sheet(tmp_path, tmp_path / 'studies.csv')
    assert result.candidate_count == 2
    assert result.dataset_count == 1
    assert result.rejected_files[0].conversion_attempted
    assert not list((tmp_path / 'prepared_inputs').iterdir())
