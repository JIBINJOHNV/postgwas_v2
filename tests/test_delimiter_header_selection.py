"""Header-aware detection is opt-in; general table contracts stay unchanged."""

import pytest

from postgwas.core.io.delimiters import detect_delimiter


@pytest.mark.parametrize("header_check, expected", [(None, " "), (lambda row: all(value.strip() for value in row), "\t")])
def test_header_check_excludes_artificial_empty_columns(tmp_path, header_check, expected):
    source = tmp_path / "padded.tsv"
    source.write_text("CHR \t BP \t BETA\n1 \t 100 \t 0.1\n")
    original = source.read_bytes()

    result = detect_delimiter(
        source, candidates=["tab", "space"], minimum_columns=2,
        maximum_columns=10, sample_lines=2, header_check=header_check,
    )

    assert result.value == expected
    assert source.read_bytes() == original


def test_header_check_can_reject_all_consistent_candidates(tmp_path):
    source = tmp_path / "empty_header.tsv"
    source.write_text("CHR\t\tBETA\n1\t100\t0.1\n")

    result = detect_delimiter(
        source, candidates=["tab", "space"], minimum_columns=2,
        maximum_columns=10, sample_lines=2, fail_policy="unknown",
        header_check=lambda row: all(value.strip() for value in row),
    )

    assert result.kind == "unknown"


def test_headerless_detection_keeps_numeric_first_row(tmp_path):
    source = tmp_path / "headerless.txt"
    source.write_text("1|100|0.1\n2|200|0.2\n")

    result = detect_delimiter(
        source, candidates=["tab", "pipe"], minimum_columns=2,
        maximum_columns=10, sample_lines=2,
    )

    assert result.value == "|"
    assert result.column_count == 3
