"""Metadata is excluded before parsing without changing scientific table data."""

import bz2
import gzip
import lzma
from pathlib import Path
import re
from types import SimpleNamespace
import zipfile

import polars as pl
import pytest

from postgwas.core.io.delimiters import open_binary
from postgwas.modules.harmonisation import summary_statistics_io
from postgwas.modules.harmonisation.concordance.analysis import (
    INPUT_SOURCE_ROW_COLUMN,
    PARTITION_COLUMN,
)
from postgwas.modules.harmonisation.concordance.service import _stage_input
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.summary_statistics_io import (
    load_summary_statistics_table,
)


MAPPING = {
    "chr_col": "#CHROM", "pos_col": "POS", "snp_id_col": "ID",
    "beta_or_col": "BETA", "se_col": "SE", "pval_col": "PVAL",
}
HEADER = "#CHROM\tPOS\tID\tBETA\tSE\tPVAL\tNOTE\n"
TABLE = HEADER + "2\t20\trs_second\t-0.15\t0.07\t1e-400\tfirst\n" + (
    "1\t10\trs_first\t0.2\t0.1\t0.0455\tsecond\n"
)
EXPECTED_ROWS = [
    ("2", "20", "rs_second", "-0.15", "0.07", "1e-400", "first"),
    ("1", "10", "rs_first", "0.2", "0.1", "0.0455", "second"),
]


def _source(tmp_path, suffix, content):
    path = tmp_path / ("study" + suffix)
    payload = content.encode("utf-8")
    if suffix in (".gz", ".bgz"):
        path.write_bytes(gzip.compress(payload))
    elif suffix == ".bz2":
        path.write_bytes(bz2.compress(payload))
    elif suffix == ".xz":
        path.write_bytes(lzma.compress(payload))
    elif suffix == ".zip":
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("study.tsv", payload)
    else:
        path.write_bytes(payload)
    return path


def _read(path, **overrides):
    policies = default_policies().with_overrides({
        "input.delimiter": "tab", **overrides,
    })
    return load_summary_statistics_table(str(path), MAPPING.copy(), policies=policies)


@pytest.mark.parametrize("suffix", [".tsv", ".gz", ".bgz", ".bz2", ".xz", ".zip"])
@pytest.mark.parametrize("line_count,bad_index", [(1, 0), (4, 2), (79, 71)])
def test_malformed_metadata_is_skipped_at_variable_positions_and_lengths(
    tmp_path, suffix, line_count, bad_index,
):
    # Changing a fixed line offset or repairing only a particular metadata ID
    # cannot satisfy all cases. The malformed content is deliberately opaque.
    metadata = ["##ordinary metadata\n"] * line_count
    metadata[bad_index] = '##arbitrary=<Description="example">tail">\n'
    path = _source(tmp_path, suffix, "".join(metadata) + TABLE)
    original = path.read_bytes()

    frame, mapping = _read(path)

    assert frame.columns == ["#CHROM", "POS", "ID", "BETA", "SE", "PVAL", "NOTE"]
    assert frame.rows() == EXPECTED_ROWS
    assert mapping == MAPPING
    assert path.read_bytes() == original
    assert list(tmp_path.iterdir()) == [path]


@pytest.mark.parametrize("suffix", [".tsv", ".gz", ".bz2", ".xz", ".zip"])
@pytest.mark.parametrize("metadata", ["", "##valid metadata\n"])
def test_valid_tables_work_with_or_without_a_metadata_preamble(tmp_path, suffix, metadata):
    path = _source(tmp_path, suffix, metadata + TABLE)
    original = path.read_bytes()

    frame, _ = _read(path)

    assert frame.rows() == EXPECTED_ROWS
    assert path.read_bytes() == original


@pytest.mark.parametrize("suffix", [".tsv", ".gz"])
def test_crlf_preamble_and_quoted_metadata_like_data_remain_distinct(tmp_path, suffix):
    # Filtering every physical line beginning ## would corrupt this quoted
    # field. Quoting also carries an embedded tab and escaped double quote.
    content = '\r\n##broken="description\r\n\r\n' + (
        HEADER.replace("\n", "\r\n")
        + '1\t10\trs1\t0.2\t0.1\t0.05\t"first line\r\n'
        + '##this is data\r\nlast\tline ""quoted"""\r\n'
    )
    path = _source(tmp_path, suffix, content)
    original = path.read_bytes()

    frame, _ = _read(path)

    assert frame.rows() == [(
        "1", "10", "rs1", "0.2", "0.1", "0.05",
        'first line\r\n##this is data\r\nlast\tline "quoted"',
    )]
    assert path.read_bytes() == original


def test_disabled_metadata_policy_does_not_silently_remove_the_preamble(tmp_path):
    path = _source(tmp_path, ".tsv", "##metadata\n" + TABLE)
    original = path.read_bytes()

    with pytest.raises((pl.exceptions.PolarsError, RuntimeError)):
        _read(path, **{"input.strip_double_hash_lines": False})

    assert path.read_bytes() == original


@pytest.mark.parametrize("content", ["", "##metadata\n", "##metadata\n" + HEADER])
def test_empty_or_header_only_input_cannot_become_a_successful_table(tmp_path, content):
    path = _source(tmp_path, ".gz", content)

    with pytest.raises((pl.exceptions.NoDataError, RuntimeError), match="[Ee]mpty"):
        _read(path)


def test_malformed_data_rows_still_fail_after_metadata_is_excluded(tmp_path):
    path = _source(tmp_path, ".gz", '##metadata="unclosed\n' + TABLE + (
        "1\t30\trs3\t0.3\t0.1\t0.01\tannotation\textra\n"
    ))
    original = path.read_bytes()

    with pytest.raises(pl.exceptions.ComputeError, match="more fields"):
        _read(path)

    assert path.read_bytes() == original


def test_metadata_exclusion_reports_the_count_observed_in_this_file(tmp_path, capsys):
    path = _source(tmp_path, ".tsv", "##opaque metadata\n" * 7 + TABLE)

    frame, _ = _read(path)

    assert frame.rows() == EXPECTED_ROWS
    transcript = capsys.readouterr().out
    assert re.search(r"(?:\b7\b[^\n]*metadata|metadata[^\n]*\b7\b)", transcript, re.I)


@pytest.mark.parametrize("suffix", [".tsv", ".gz"])
def test_memory_retry_rewinds_a_consumed_table_stream(tmp_path, monkeypatch, suffix):
    metadata = ["##opaque metadata\n"] * 11
    metadata[3] = '##another_field="unterminated description\n'
    path = _source(tmp_path, suffix, "".join(metadata) + TABLE)
    original = path.read_bytes()
    native_reader = pl.read_csv
    consumed = []

    def fail_after_consuming_first_source(source, **options):
        if not consumed:
            # Old production code passes a path, so keep its expected failure
            # at the parser boundary rather than manufacturing a test error.
            if isinstance(source, (str, Path)):
                with open_binary(source) as handle:
                    consumed.append(handle.read())
            else:
                consumed.append(source.read())
            raise MemoryError("allocation failed after reading the stream")
        return native_reader(source, **options)

    monkeypatch.setattr(summary_statistics_io.pl, "read_csv", fail_after_consuming_first_source)

    frame, _ = _read(path)

    assert frame.rows() == EXPECTED_ROWS
    assert consumed == [TABLE.encode("utf-8")]
    assert path.read_bytes() == original


@pytest.mark.parametrize("suffix", [".tsv", ".gz"])
@pytest.mark.parametrize("metadata_lines", [2, 39])
def test_concordance_staging_preserves_tokens_and_source_rows_after_metadata(
    tmp_path, suffix, metadata_lines,
):
    metadata = ["##ordinary metadata\n"] * metadata_lines
    metadata[-1] = '##any_identifier=<Description="text">broken">\n'
    path = _source(tmp_path, suffix, "".join(metadata) + TABLE)
    original = path.read_bytes()
    row = SimpleNamespace(
        input_file=path, delimiter="tab", chromosome_column="#CHROM",
        position_column="POS", chromosome_position_column=None,
        effect_allele_column=None, other_allele_column=None,
        effect_allele_frequency_column=None, effect_column="BETA",
        standard_error_column="SE", z_score_column=None, p_value_column="PVAL",
    )
    destination = tmp_path / "staged.parquet"

    count = _stage_input(
        row, destination, default_policies(),
        SimpleNamespace(info=lambda *_args, **_kwargs: None),
        batch_rows=1, compression="zstd",
    )

    staged = pl.read_parquet(destination)
    assert count == 2
    assert staged.select([
        INPUT_SOURCE_ROW_COLUMN, "#CHROM", "POS", "BETA", "SE", "PVAL", PARTITION_COLUMN,
    ]).rows() == [
        (1, "2", "20", "-0.15", "0.07", "1e-400", "2"),
        (2, "1", "10", "0.2", "0.1", "0.0455", "1"),
    ]
    assert path.read_bytes() == original
