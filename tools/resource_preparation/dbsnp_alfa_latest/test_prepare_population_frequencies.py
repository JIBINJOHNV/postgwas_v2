"""Focused tests for ALFA population-frequency conversion."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest


SCRIPT = Path(__file__).with_name("prepare_population_frequencies.py")
SPEC = spec_from_file_location("prepare_population_frequencies", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_population_af_pools_counts_and_distinguishes_zero_from_missing() -> None:
    assert MODULE.population_af(["100:1", "300:9"], ["AN", "AC"], 0) == "0.025"
    assert MODULE.population_af(["100:0"], ["AN", "AC"], 0) == "0"
    assert MODULE.population_af(["0:0", ".:."], ["AN", "AC"], 0) == "."


def test_convert_record_splits_multiallelic_alt_frequencies() -> None:
    line = "NC_000022.10\t16050075\trs1\tA\tG,T\t.\t.\t.\tAN:AC\t100:1,2\t300:9,12\n"
    indexes = {"AMR": (9, 10)}
    records = MODULE.convert_record(line, indexes, {"NC_000022.10": "22"})
    assert records == [
        "22\t16050075\trs1\tA\tG\t.\t.\tAMR=0.025\n",
        "22\t16050075\trs1\tA\tT\t.\t.\tAMR=0.035\n",
    ]


def test_population_af_rejects_impossible_counts() -> None:
    with pytest.raises(RuntimeError, match="invalid ALFA allele counts"):
        MODULE.population_af(["10:11"], ["AN", "AC"], 0)
