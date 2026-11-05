"""A combined-coordinate source may also supply a study identifier safely."""

import csv
import gzip
from pathlib import Path

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.modules.harmonisation.coordinates import harmonise_coordinates_and_alleles
from postgwas.modules.harmonisation.input_validation import validate_header
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.shared.column_mapping import (
    validate_scientific_column_mappings,
)
from postgwas.modules.harmonisation.summary_statistics_io import read_summary_statistics
from postgwas.modules.harmonisation.variant_identifiers import harmonise_variant_identifiers


def _mapping(source_column="SNP", **updates):
    return {
        "gwas_outputname": "shared_coordinate_id", "delimiter": "tab",
        "chr_pos_col": source_column, "snp_id_col": source_column,
        "ea_col": "A1", "oa_col": "A2", "beta_or_col": "Beta",
        "se_col": "SE", "imp_z_col": "Z", "pval_col": "P", **updates,
    }


def _write_study(tmp_path, identifiers, *, source_column="SNP", extra_columns=None):
    values = {
        source_column: identifiers, "A1": ["A", "G"], "A2": ["G", "C"],
        "Z": ["0.841", "-1.172"], "P": ["0.4003", "0.2411"],
        "N": ["56994", "56994"], "SE": ["0.00593735", "0.00709625"],
        "Beta": ["0.00499331", "-0.0083168"], **(extra_columns or {}),
    }
    source = tmp_path / "study.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(values)
        writer.writerows(zip(*values.values()))
    return source, pl.DataFrame(values)


def _read(source, tmp_path, mapping):
    policies = default_policies()
    valid, problems = validate_header(mapping, str(source), policies=policies)
    assert valid, problems
    return read_summary_statistics(
        str(source), str(tmp_path / "read"), mapping,
        dict(load_configuration().modules.harmonisation.output_layout.root),
        "\t", policies=policies, input_line_count=2,
    )


@pytest.mark.parametrize(
    ("source_column", "identifiers", "chromosomes", "positions", "working_ids"),
    [
        ("SNP", ["10:10000018", "10:100000625"], ["10", "10"],
         [10000018, 100000625], ["10_10000018", "10_100000625"]),
        ("marker", ["1:100:A:G", "1:100:G:C"], ["1", "1"],
         [100, 100], ["1_100_A_G", "1_100_G_C"]),
    ],
)
def test_shared_source_reaches_reader_and_identifier_steps_without_changing_source(
    tmp_path, source_column, identifiers, chromosomes, positions, working_ids,
):
    source, expected = _write_study(
        tmp_path, identifiers, source_column=source_column,
        extra_columns={"snp_col": ["occupied", "also occupied"],
                       "snp_col_": ["preserve one", "preserve two"]},
    )
    original_bytes = source.read_bytes()
    read = _read(source, tmp_path, _mapping(source_column))
    frame, snapshot, _, rows_read, mapping = read[:5]

    assert rows_read == frame.height == 2
    assert all(count == 0 for count in read[5:10])
    assert snapshot.select(expected.columns).equals(expected)
    assert [name for name in snapshot.columns if name in expected.columns] == expected.columns
    assert frame[mapping["chr_col"]].to_list() == chromosomes
    assert frame[mapping["pos_col"]].to_list() == positions
    destination = mapping["snp_id_col"]
    assert destination not in expected.columns
    assert frame[destination].to_list() == identifiers
    validate_scientific_column_mappings(mapping, stage="shared source resolved")

    result, resolved = harmonise_variant_identifiers("all", frame, mapping)
    assert result[destination].to_list() == working_ids
    assert result[source_column].to_list() == identifiers
    for column in expected.columns:
        if column in {"Z", "P", "SE", "Beta"}:
            assert result[column].cast(pl.Float64).to_list() == (
                expected[column].cast(pl.Float64).to_list()
            )
        else:
            assert result[column].to_list() == expected[column].to_list()
    assert snapshot.select(expected.columns).equals(expected)
    assert source.read_bytes() == original_bytes

    repeated, repeated_mapping = harmonise_variant_identifiers("all", result, dict(resolved))
    assert repeated.equals(result)
    assert repeated_mapping == resolved


@pytest.mark.parametrize(
    "selected",
    [{}, {"chr_col": "CHR"}, {"pos_col": "BP"}, {"chr_col": "CHR", "pos_col": "BP"}],
)
def test_shared_source_keeps_explicit_coordinates_and_is_resolved_once(tmp_path, selected):
    frame = pl.DataFrame({"SNP": ["1:100", "1:101"], "CHR": ["1", "1"], "BP": [100, 101]})
    original = frame.clone()
    logger = PipelineLogger(
        "shared_coordinate_id", "all", str(tmp_path),
        policies=default_policies(), screen_level="ERROR",
    )
    try:
        result, mapping, _ = harmonise_coordinates_and_alleles(
            "all", frame, {"chr_pos_col": "SNP", "snp_id_col": "SNP", **selected},
            logger=logger, return_qc_info=True,
        )
    finally:
        logger.close()

    assert frame.equals(original)
    assert result.select(frame.columns).equals(frame)
    for role, name in selected.items():
        assert mapping[role] == name
    assert result[mapping["chr_col"]].to_list() == ["1", "1"]
    assert result[mapping["pos_col"]].to_list() == [100, 101]
    destination = mapping["snp_id_col"]
    assert destination not in frame.columns
    assert result[destination].to_list() == ["1:100", "1:101"]
    log = Path(logger.log_path).read_text(encoding="utf-8")
    assert "source_column=SNP" in log
    assert "identifier_column=%s" % destination in log
    assert "copied_rows=2" in log

    repeated, repeated_mapping = harmonise_coordinates_and_alleles("all", result, dict(mapping))
    assert repeated.equals(result)
    assert repeated_mapping == mapping


def test_shared_identifier_copy_avoids_existing_and_reserved_column_names():
    frame = pl.DataFrame({
        "SNP": ["1:100"], "SNP_": ["occupied"], "SNP__": ["also occupied"],
    })
    mapping = {"chr_pos_col": "SNP", "snp_id_col": "SNP", "se_col": "SNP___"}
    result, resolved = harmonise_coordinates_and_alleles("all", frame, mapping)
    destination = resolved["snp_id_col"]
    assert destination not in {*frame.columns, "SNP___"}
    assert result.select(frame.columns).equals(frame)
    assert result[destination].to_list() == ["1:100"]


@pytest.mark.parametrize("selected", [{"chr_col": "CHR"}, {"pos_col": "BP"}])
def test_shared_source_does_not_hide_conflicting_explicit_coordinates(selected):
    frame = pl.DataFrame({"SNP": ["1:100"], "CHR": ["2"], "BP": [101]})
    original = frame.clone()
    with pytest.raises(ValueError, match="conflict"):
        harmonise_coordinates_and_alleles(
            "all", frame, {"chr_pos_col": "SNP", "snp_id_col": "SNP", **selected},
        )
    assert frame.equals(original)


@pytest.mark.parametrize("invalid", [None, "1:bad"])
def test_shared_source_keeps_existing_invalid_coordinate_rejections(tmp_path, invalid):
    source, expected = _write_study(tmp_path, ["1:100", invalid])
    read = _read(source, tmp_path, _mapping())
    frame, snapshot, _, rows_read, mapping = read[:5]
    assert rows_read == 2
    assert frame.height == 1
    assert frame["SNP"].to_list() == ["1:100"]
    assert frame[mapping["snp_id_col"]].to_list() == ["1:100"]
    assert snapshot["SNP"].to_list() == expected["SNP"].to_list()
    assert sum(read[5:10]) == 1


def test_shared_source_with_unparseable_coordinates_still_fails(tmp_path):
    source, _ = _write_study(tmp_path, ["rs1", "rs2"])
    with pytest.raises(RuntimeError, match="Could not split"):
        _read(source, tmp_path, _mapping())


@pytest.mark.parametrize(
    "updates",
    [{"beta_or_col": "SNP"}, {"pos_col": "SNP"}, {"se_col": "Beta"}],
)
def test_shared_source_exception_does_not_allow_other_role_collisions(tmp_path, updates):
    source, _ = _write_study(tmp_path, ["1:100", "1:101"])
    mapping = _mapping(**updates)
    valid, problems = validate_header(mapping, str(source), policies=default_policies())
    assert not valid
    assert any(problem.category == "ambiguous_mapping" for problem in problems)
    policies = default_policies().with_overrides({"validation.on_ambiguous_column_mapping": "warn"})
    with pytest.raises(ValueError, match="incompatible scientific column mappings"):
        read_summary_statistics(
            str(source), str(tmp_path / "read"), mapping,
            dict(load_configuration().modules.harmonisation.output_layout.root),
            "\t", policies=policies, input_line_count=2,
        )
