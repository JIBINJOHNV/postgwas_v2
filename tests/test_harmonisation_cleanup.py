"""Regression tests for retry-safe harmonisation cleanup."""

from pathlib import Path

import polars as pl

from postgwas.config import load_configuration
from postgwas.core.paths import configured_output_path
from postgwas.modules.harmonisation.cleanup import (
    remove_partial_chromosome_outputs,
)


def test_retry_cleanup_preserves_inputs_and_removes_partial_outputs(tmp_path):
    layout = dict(load_configuration().modules.harmonisation.output_layout.root)
    values = {
        "dataset_id": "study",
        "chromosome": "1",
        "build": "GRCh37",
        "target_build": "GRCh38",
    }
    chromosome_table = configured_output_path(
        tmp_path, layout["chromosome_table"], **values,
    )
    source_snapshot = configured_output_path(
        tmp_path, layout["chromosome_source_snapshot"], **values,
    )
    rejected = configured_output_path(
        tmp_path, layout["chromosome_reject"], **values,
    )
    adapter_input = configured_output_path(
        tmp_path, layout["adapter_input"], **values,
    )
    partial_vcf = configured_output_path(
        tmp_path, layout["chromosome_raw_vcf"], **values,
    )
    partial_index = Path(str(partial_vcf) + ".tbi")

    partial_files = (
        chromosome_table, source_snapshot, rejected, adapter_input, partial_vcf,
    )
    for path in partial_files:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"partial\n")
    partial_index.write_bytes(b"index\n")

    removed = remove_partial_chromosome_outputs(
        tmp_path, "study", "1", layout,
    )

    assert removed == 4
    assert chromosome_table.is_file()
    retained = pl.read_csv(
        chromosome_table, separator="\t", has_header=False,
    )
    assert retained.height == 1
    assert source_snapshot.is_file()
    assert not rejected.exists()
    assert not adapter_input.exists()
    assert not partial_vcf.exists()
    assert not partial_index.exists()
