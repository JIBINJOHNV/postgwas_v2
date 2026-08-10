import gzip
import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest


MODULE_PATH = Path(__file__).with_name("prepare_info_scores.py")
SPEC = importlib.util.spec_from_file_location("prepare_finngen_info", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_prepares_original_mean_and_median_tracks(tmp_path: Path) -> None:
    source = tmp_path / "input.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(
            "#variant\tchr\tpos\tref\talt\tINFO\tINFO_batch1\tINFO_batch2\t"
            "INFO_batch3\n"
            "1:10:A:G\t1\t10\tA\tG\t0.8\t0.6\t0.8\t1.0\n"
            "1:20:C:T\t1\t20\tC\tT\t0.4\t0.2\tNA\t0.6\n"
            "23:30:G:A\t23\t30\tG\tA\t0.7\t0.5\t0.7\t0.9\n"
        )
    output = tmp_path / "output"

    counts = MODULE.prepare_info_scores(source, output)

    assert counts["1"] == 2
    assert counts["X"] == 1
    expected = {
        ("original", "1"): [("1", "10", "A", "G", 0.8), ("1", "20", "C", "T", 0.4)],
        ("mean", "1"): [("1", "10", "A", "G", 0.8), ("1", "20", "C", "T", 0.4)],
        ("median", "1"): [("1", "10", "A", "G", 0.8), ("1", "20", "C", "T", 0.4)],
        ("original", "X"): [("X", "30", "G", "A", 0.7)],
        ("mean", "X"): [("X", "30", "G", "A", 0.7)],
        ("median", "X"): [("X", "30", "G", "A", 0.7)],
    }
    for (statistic, chrom), rows in expected.items():
        path = MODULE.output_paths(output, chrom)[statistic]
        observed = subprocess.run(
            [MODULE.TABIX_EXECUTABLE, str(path), chrom],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        ).stdout.splitlines()
        assert len(observed) == len(rows)
        for line, expected_row in zip(observed, rows, strict=True):
            fields = line.split("\t")
            assert tuple(fields[:4]) == expected_row[:4]
            assert float(fields[4]) == pytest.approx(expected_row[4])
