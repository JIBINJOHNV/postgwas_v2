"""Regression tests for one-pass user EAF/INFO whole-file staging."""

import gzip
from pathlib import Path

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation import external_reference_staging as staging_module
from postgwas.modules.harmonisation.external_reference_staging import (
    ExternalReferenceStagingError,
    stage_shared_external_reference_files,
)
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.shared import (
    variant_columns as variant_columns_module,
)
from postgwas.modules.harmonisation.shared.variant_columns import (
    canonicalize_variant_frame,
    read_reference_variant_table,
)


MAPPING = {
    "chr": "CHROM",
    "pos": "POS",
    "a1": "ALT",
    "a2": "REF",
    "delimiter": "tab",
}


def _settings_and_layout():
    harmonisation = load_configuration().modules.harmonisation
    return (
        harmonisation.external_reference_staging.model_dump(),
        dict(harmonisation.output_layout.root),
    )


def _maps(chromosomes, *, eaf=None, info=None):
    return {
        chromosome: {
            "user_eaf_file": None if eaf is None else str(eaf),
            "user_info_file": None if info is None else str(info),
        }
        for chromosome in chromosomes
    }


def test_whole_user_eaf_is_projected_and_partitioned_once(tmp_path):
    source = tmp_path / "whole_eaf.tsv"
    source.write_text(
        "CHROM\tPOS\tREF\tALT\tEAF\tUNUSED\n"
        "chr1\t10\tA\tG\t0.1\tone\n"
        "02\t20\tC\tT\t0.2\ttwo\n"
        "23\t30\tG\tA\t0.3\tthree\n"
        "3\t40\tT\tC\t0.4\tfour\n",
        encoding="utf-8",
    )
    chromosomes = ("1", "2", "X")
    settings, layout = _settings_and_layout()
    settings["batch_rows"] = 1

    staged, summary = stage_shared_external_reference_files(
        _maps(chromosomes, eaf=source),
        chromosomes=chromosomes,
        dataset_id="study",
        output_directory=tmp_path / "output",
        output_layout=layout,
        user_eaf_specification=str(source),
        user_eaf_column="EAF",
        external_eaf_mapping=MAPPING,
        user_info_specification=None,
        user_info_column=None,
        external_info_mapping=MAPPING,
        settings=settings,
        policies=default_policies(),
    )

    assert summary["status"] == "staged"
    assert summary["source_files_read"] == 1
    resource = summary["resources"][0]
    assert resource["rows_read"] == 4
    assert resource["batch_rows"] == 1
    assert resource["rows_by_chromosome"] == {"1": 1, "2": 1, "X": 1}
    assert resource["rows_not_in_observed_chromosomes"] == 1
    assert resource["selected_columns"] == ["CHROM", "POS", "ALT", "REF", "EAF"]

    for chromosome, original_label in (("1", "chr1"), ("2", "02"), ("X", "23")):
        partition = Path(staged[chromosome]["user_eaf_file"])
        assert partition != source
        assert partition.is_file()
        frame = pl.read_parquet(partition)
        assert frame.columns == ["CHROM", "POS", "ALT", "REF", "EAF"]
        assert frame.get_column("CHROM").to_list() == [original_label]
        assert "UNUSED" not in frame.columns


def test_staged_partition_uses_existing_reference_validation_path(tmp_path):
    source = tmp_path / "whole_eaf.tsv"
    source.write_text(
        "CHROM\tPOS\tREF\tALT\tEAF\n"
        "chr01\t10\tA\tG\t0.25\n"
        "2\t20\tC\tT\t0.40\n",
        encoding="utf-8",
    )
    chromosomes = ("1", "2")
    settings, layout = _settings_and_layout()
    policies = default_policies()
    staged, _summary = stage_shared_external_reference_files(
        _maps(chromosomes, eaf=source),
        chromosomes=chromosomes,
        dataset_id="study",
        output_directory=tmp_path / "output",
        output_layout=layout,
        user_eaf_specification=str(source),
        user_eaf_column="EAF",
        external_eaf_mapping=MAPPING,
        user_info_specification=None,
        user_info_column=None,
        external_info_mapping=MAPPING,
        settings=settings,
        policies=policies,
    )

    frame, detected = read_reference_variant_table(
        staged["1"]["user_eaf_file"],
        MAPPING,
        policies,
        value_columns=["EAF"],
        description="external frequency file",
    )
    normalized, stats = canonicalize_variant_frame(
        frame,
        {"chr": "CHROM", "pos": "POS", "ea": "ALT", "oa": "REF"},
        policies,
        label="external frequency file",
    )

    assert detected.method == "dataset-level projected Parquet partition"
    assert normalized.row(0, named=True) == {
        "CHROM": "1", "POS": 10, "ALT": "G", "REF": "A", "EAF": "0.25",
    }
    assert stats["positions_unreadable"] == 0


def test_direct_parquet_reference_projects_required_columns_at_read(
    monkeypatch, tmp_path,
):
    source = tmp_path / "wide_reference.parquet"
    pl.DataFrame({
        "CHROM": ["1"],
        "POS": [10],
        "REF": ["A"],
        "ALT": ["G"],
        "EUR": [0.25],
        "AFR": [0.35],
        "EAS": [0.15],
        "UNUSED": ["payload"],
    }).write_parquet(source)
    observed = {}
    original = variant_columns_module.read_parquet_table

    def tracked_read_parquet_table(*args, **kwargs):
        observed["columns"] = kwargs.get("columns")
        return original(*args, **kwargs)

    monkeypatch.setattr(
        variant_columns_module,
        "read_parquet_table",
        tracked_read_parquet_table,
    )
    frame, detected = read_reference_variant_table(
        source,
        MAPPING,
        default_policies(),
        value_columns=["EUR"],
        description="population frequency reference",
    )

    assert observed["columns"] == ["CHROM", "POS", "ALT", "REF", "EUR"]
    assert frame.columns == ["CHROM", "POS", "ALT", "REF", "EUR"]
    assert detected.method == "dataset-level projected Parquet partition"


def test_direct_compressed_reference_projects_required_columns_at_read(
    monkeypatch, tmp_path,
):
    source = tmp_path / "wide_reference.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(
            "CHROM\tPOS\tREF\tALT\tEUR\tAFR\tEAS\tUNUSED\n"
            "1\t10\tA\tG\t0.25\t0.35\t0.15\tpayload\n"
        )
    observed = {}
    original = variant_columns_module.read_delimited_table

    def tracked_read_delimited_table(*args, **kwargs):
        observed["columns"] = kwargs.get("columns")
        return original(*args, **kwargs)

    monkeypatch.setattr(
        variant_columns_module,
        "read_delimited_table",
        tracked_read_delimited_table,
    )
    frame, detected = read_reference_variant_table(
        source,
        MAPPING,
        default_policies(),
        value_columns=["EUR"],
        description="population frequency reference",
    )

    assert observed["columns"] == ["CHROM", "POS", "ALT", "REF", "EUR"]
    assert frame.columns == ["CHROM", "POS", "ALT", "REF", "EUR"]
    assert detected.value == "\t"


def test_compressed_whole_user_file_is_streamed_once(tmp_path):
    source = tmp_path / "whole_eaf.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(
            "CHROM\tPOS\tREF\tALT\tEAF\tUNUSED\n"
            "1\t10\tA\tG\t0.1\tone\n"
            "2\t20\tC\tT\t0.2\ttwo\n"
        )
    chromosomes = ("1", "2")
    settings, layout = _settings_and_layout()

    staged, summary = stage_shared_external_reference_files(
        _maps(chromosomes, eaf=source),
        chromosomes=chromosomes,
        dataset_id="study",
        output_directory=tmp_path / "output",
        output_layout=layout,
        user_eaf_specification=str(source),
        user_eaf_column="EAF",
        external_eaf_mapping=MAPPING,
        user_info_specification=None,
        user_info_column=None,
        external_info_mapping=MAPPING,
        settings=settings,
        policies=default_policies(),
    )

    assert summary["source_files_read"] == 1
    assert summary["resources"][0]["rows_read"] == 2
    assert pl.read_parquet(staged["2"]["user_eaf_file"]).row(0, named=True) == {
        "CHROM": "2", "POS": "20", "ALT": "T", "REF": "C", "EAF": "0.2",
    }


def test_eaf_and_info_in_same_whole_file_are_read_as_one_projection(
    monkeypatch, tmp_path,
):
    source = tmp_path / "whole_annotations.tsv"
    source.write_text(
        "CHROM\tPOS\tREF\tALT\tEAF\tINFO\tUNUSED\n"
        "1\t10\tA\tG\t0.1\t0.9\tone\n"
        "2\t20\tC\tT\t0.2\t0.8\ttwo\n",
        encoding="utf-8",
    )
    chromosomes = ("1", "2")
    settings, layout = _settings_and_layout()
    source_scans = []
    original_batches = staging_module._reference_batches

    def tracked_batches(*args, **kwargs):
        source_scans.append(Path(args[0]))
        yield from original_batches(*args, **kwargs)

    monkeypatch.setattr(staging_module, "_reference_batches", tracked_batches)

    staged, summary = stage_shared_external_reference_files(
        _maps(chromosomes, eaf=source, info=source),
        chromosomes=chromosomes,
        dataset_id="study",
        output_directory=tmp_path / "output",
        output_layout=layout,
        user_eaf_specification=str(source),
        user_eaf_column="EAF",
        external_eaf_mapping=MAPPING,
        user_info_specification=str(source),
        user_info_column="INFO",
        external_info_mapping={**MAPPING, "delimiter": "auto"},
        settings=settings,
        policies=default_policies(),
    )

    assert summary["source_files_read"] == 1
    assert source_scans == [source]
    assert summary["resources"][0]["roles"] == [
        "user_eaf_file", "user_info_file",
    ]
    assert summary["resources"][0]["selected_columns"] == [
        "CHROM", "POS", "ALT", "REF", "EAF", "INFO",
    ]
    for chromosome in chromosomes:
        assert staged[chromosome]["user_eaf_file"] == staged[chromosome]["user_info_file"]
        assert pl.read_parquet(staged[chromosome]["user_eaf_file"]).columns == [
            "CHROM", "POS", "ALT", "REF", "EAF", "INFO",
        ]


@pytest.mark.parametrize(
    ("chromosomes", "specification"),
    [
        (("1", "2"), "/references/eaf_chr{chromosome}.tsv"),
        (("1",), "/references/whole_eaf.tsv"),
    ],
)
def test_template_and_single_chromosome_inputs_are_not_staged(
    tmp_path, chromosomes, specification,
):
    maps = {
        chromosome: {
            "user_eaf_file": specification.format(chromosome=chromosome),
            "user_info_file": None,
        }
        for chromosome in chromosomes
    }
    settings, layout = _settings_and_layout()

    staged, summary = stage_shared_external_reference_files(
        maps,
        chromosomes=chromosomes,
        dataset_id="study",
        output_directory=tmp_path / "output",
        output_layout=layout,
        user_eaf_specification=specification,
        user_eaf_column="EAF",
        external_eaf_mapping=MAPPING,
        user_info_specification=None,
        user_info_column=None,
        external_info_mapping=MAPPING,
        settings=settings,
        policies=default_policies(),
    )

    assert staged == maps
    assert summary["status"] == "not_needed"
    assert summary["source_files_read"] == 0
    assert not (tmp_path / "output").exists()


def test_staging_failure_is_actionable_and_does_not_publish_partitions(tmp_path):
    source = tmp_path / "whole_eaf.tsv"
    source.write_text(
        "CHROM\tPOS\tREF\tALT\n1\t10\tA\tG\n2\t20\tC\tT\n",
        encoding="utf-8",
    )
    chromosomes = ("1", "2")
    settings, layout = _settings_and_layout()

    with pytest.raises(ExternalReferenceStagingError, match="EAF"):
        stage_shared_external_reference_files(
            _maps(chromosomes, eaf=source),
            chromosomes=chromosomes,
            dataset_id="study",
            output_directory=tmp_path / "output",
            output_layout=layout,
            user_eaf_specification=str(source),
            user_eaf_column="EAF",
            external_eaf_mapping=MAPPING,
            user_info_specification=None,
            user_info_column=None,
            external_info_mapping=MAPPING,
            settings=settings,
            policies=default_policies(),
        )

    assert not list((tmp_path / "output").rglob("*.parquet"))
    assert not list((tmp_path / "output").rglob("*.tmp"))
