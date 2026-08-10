from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from postgwas.config import load_run_configuration_for_module
from postgwas.core.errors import ConfigurationError
from postgwas.modules.fine_mapping.engines.susie.summary_preparation import (
    prepare_locus_summary_statistics,
)
from postgwas.modules.fine_mapping.output_layout import resolve_output_paths


def _settings():
    config = load_run_configuration_for_module("fine_mapping")
    return config.modules.fine_mapping.summary_statistics_preparation.model_dump()


def _directories(output: Path):
    module = load_run_configuration_for_module("fine_mapping").modules.fine_mapping
    paths = resolve_output_paths(
        output, module.output_layout.model_dump(), "study"
    )
    return {
        "chromosome_cache_directory": paths["chromosome_cache_directory"],
        "locus_input_directory": paths["locus_summary_statistics_directory"],
        "manifest_directory": paths["summary_statistics_directory"],
        "quality_control_directory": paths["quality_control_directory"],
    }


def _write_sumstats(path: Path):
    pd.DataFrame(
        {
            "CHR": [1, 1, 1, 1, 2],
            "BP": [90, 100, 150, 210, 100],
            "SNP": ["rs90", "rs100", "rs150", "rs210", "rs2"],
            "REF": ["A"] * 5,
            "ALT": ["G"] * 5,
            "EZ": [1.0, 2.0, 3.0, 1.0, 2.0],
            "NEF": [1000] * 5,
            "LP": [2.0, 8.0, 9.0, 3.0, 8.0],
        }
    ).to_csv(path, sep="\t", index=False)


def _write_loci(path: Path):
    pd.DataFrame(
        {
            "CHR": [1, 1, 2],
            "START": [90, 150, 100],
            "END": [150, 210, 100],
            "GenomicLocus": ["first", "second", "third"],
        }
    ).to_csv(path, sep="\t", index=False)


def test_primary_preparation_reads_source_once_and_preserves_each_boundary(tmp_path):
    source = tmp_path / "sumstats.tsv"
    loci = tmp_path / "loci.tsv"
    _write_sumstats(source)
    _write_loci(loci)
    real_read_csv = pd.read_csv
    source_reads = 0

    def tracked_read_csv(path, *args, **kwargs):
        nonlocal source_reads
        if Path(path).resolve() == source.resolve():
            source_reads += 1
        return real_read_csv(path, *args, **kwargs)

    with patch(
        "postgwas.modules.fine_mapping.engines.susie.summary_preparation.pd.read_csv",
        side_effect=tracked_read_csv,
    ):
        result = prepare_locus_summary_statistics(
            source_file=source,
            locus_files=[loci],
            **_directories(tmp_path / "out"),
            settings=_settings(),
        )

    assert source_reads == 1
    manifest = pd.read_csv(result["locus_manifest"], sep="\t")
    assert manifest[["START", "END"]].values.tolist() == [
        [90, 150], [150, 210], [100, 100]
    ]
    first = pd.read_csv(manifest.loc[0, "Filename"], sep="\t")
    second = pd.read_csv(manifest.loc[1, "Filename"], sep="\t")
    assert first["BP"].tolist() == [90, 100, 150]
    assert second["BP"].tolist() == [150, 210]
    assert 150 in first["BP"].tolist() and 150 in second["BP"].tolist()
    assert result["n_variant_memberships"] == 6
    summary = pd.read_csv(result["preparation_summary"], sep="\t")
    assert summary.loc[0, "source_rows_scanned"] == 5
    assert summary.loc[0, "chromosome_cache_variants"] == 5
    assert summary.loc[0, "source_rows_outside_analysed_chromosomes"] == 0


def test_joint_preparation_reuses_primary_chromosome_cache(tmp_path):
    source = tmp_path / "sumstats.tsv"
    primary_loci = tmp_path / "primary.tsv"
    joint_loci = tmp_path / "joint.tsv"
    _write_sumstats(source)
    _write_loci(primary_loci)
    primary = prepare_locus_summary_statistics(
        source_file=source,
        locus_files=[primary_loci],
        **_directories(tmp_path / "primary"),
        settings=_settings(),
    )
    pd.DataFrame(
        {
            "CHR": [1], "START": [90], "END": [210],
            "GenomicLocus": ["joint"],
        }
    ).to_csv(joint_loci, sep="\t", index=False)

    with patch(
        "postgwas.modules.fine_mapping.engines.susie.summary_preparation._stream_chromosome_cache",
        side_effect=AssertionError("source must not be reread for the joint round"),
    ):
        joint = prepare_locus_summary_statistics(
            source_file=source,
            locus_files=[joint_loci],
            **_directories(tmp_path / "joint"),
            settings=_settings(),
            existing_chromosome_manifest=primary["chromosome_manifest"],
        )

    manifest = pd.read_csv(joint["locus_manifest"], sep="\t")
    prepared = pd.read_csv(manifest.loc[0, "Filename"], sep="\t")
    assert joint["cache_reused"] is True
    assert prepared["BP"].tolist() == [90, 100, 150, 210]


def test_preparation_writes_header_only_file_for_locus_without_variants(tmp_path):
    source = tmp_path / "sumstats.tsv"
    loci = tmp_path / "empty_locus.tsv"
    _write_sumstats(source)
    pd.DataFrame(
        {"CHR": [1], "START": [500], "END": [600], "GenomicLocus": ["empty"]}
    ).to_csv(loci, sep="\t", index=False)
    result = prepare_locus_summary_statistics(
        source_file=source,
        locus_files=[loci],
        **_directories(tmp_path / "out"),
        settings=_settings(),
    )
    manifest = pd.read_csv(result["locus_manifest"], sep="\t")
    prepared = pd.read_csv(manifest.loc[0, "Filename"], sep="\t")
    assert manifest.loc[0, "n_variants"] == 0
    assert prepared.empty
    assert list(prepared.columns) == [
        "CHR", "BP", "SNP", "REF", "ALT", "EZ", "NEF", "LP"
    ]


def test_preparation_rejects_invalid_rows_on_analyzed_chromosome(tmp_path):
    source = tmp_path / "sumstats.tsv"
    loci = tmp_path / "loci.tsv"
    _write_sumstats(source)
    table = pd.read_csv(source, sep="\t")
    table.loc[0, "NEF"] = 0
    table.to_csv(source, sep="\t", index=False)
    pd.DataFrame(
        {"CHR": [1], "START": [1], "END": [500], "GenomicLocus": ["bad"]}
    ).to_csv(loci, sep="\t", index=False)
    with pytest.raises(ValueError, match="invalid row"):
        prepare_locus_summary_statistics(
            source_file=source,
            locus_files=[loci],
            **_directories(tmp_path / "out"),
            settings=_settings(),
        )


def test_preparation_configuration_rejects_unsafe_filename_template(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "summary_statistics_preparation:\n"
        "  locus_filename_template: ../locus_{index}.tsv\n"
    )
    with pytest.raises(ConfigurationError, match="must produce basenames"):
        load_run_configuration_for_module("fine_mapping", config_file)
