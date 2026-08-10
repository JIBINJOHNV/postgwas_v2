"""Regression tests for post-primary connected-locus joint reruns."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytest

from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.modules.fine_mapping.overlap_resolution import (
    _read_primary_result_loci,
    build_connected_overlap_plan,
    resolve_overlapping_results,
)
from postgwas.modules.fine_mapping.output_layout import resolve_output_paths
from postgwas.modules.fine_mapping.service import _resolve_fine_mapping_arguments


def _overlap_args(settings, output_directory, **overrides):
    module = load_configuration().modules.fine_mapping
    values = {
        "overlap_resolution": settings,
        "fine_mapping_output_layout": module.output_layout.model_dump(),
        "dataset_id": "study",
        "finemap_method": "susie",
        "lp_threshold": 7.3,
        "locus_file": "unused.tsv",
        "locus_type": "range",
        "window_kb": 1500,
        "output_directory": str(output_directory),
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def _write_handoff(base: Path, rows: list[tuple[str, str]], label: str) -> dict:
    flames = base / "flames_input"
    flames.mkdir(parents=True)
    index_rows = []
    manifest_rows = []
    locus_rows = []
    for sequence, (locus, coordinates) in enumerate(rows, start=1):
        chromosome, interval = coordinates.split(":")
        start, end = interval.split("-")
        filename = flames / f"{label}_{sequence}.txt"
        filename.write_text("index\tcred1\tprob1\n1\t1:1:A_G\t0.95\n", encoding="utf-8")
        index_rows.append({
            "Filename": str(filename.resolve()),
            "GenomicLocus": locus,
            "Annotfiles": str((flames / f"annotated_{label}_{sequence}.txt").resolve()),
        })
        manifest_rows.append({
            "Filename": filename.name,
            "GenomicLocus": locus,
            "engine": label,
            "credible_set": sequence,
        })
        locus_rows.append({
            "GenomicLocus": locus,
            "chr": chromosome.replace("chr", ""),
            "start": int(start),
            "end": int(end),
        })
    index = flames / "indexfile.txt"
    manifest = flames / f"{label}_manifest.tsv"
    genomic_loci = flames / "genomic_loci.tsv"
    status = base / f"{label}_status.tsv"
    pd.DataFrame(index_rows).to_csv(index, sep="\t", index=False)
    pd.DataFrame(manifest_rows).to_csv(manifest, sep="\t", index=False)
    pd.DataFrame(locus_rows).to_csv(genomic_loci, sep="\t", index=False)
    pd.DataFrame([
        {
            "genomic_locus": locus,
            "warning_reason": f"{label}_warning_{sequence}",
            "failure_reason": "",
        }
        for sequence, (locus, _) in enumerate(rows, start=1)
    ]).to_csv(status, sep="\t", index=False)
    return {
        "status": "success",
        "output_dir": str(base),
        "flames_input": str(flames),
        "flames_index": str(index),
        "credible_set_manifests": [str(manifest)],
        "locus_status": str(status),
        "genomic_loci": str(genomic_loci),
    }


def test_connected_overlap_plan_is_transitive_and_uses_closed_boundaries():
    loci = pd.DataFrame({
        "GenomicLocus": ["A", "B", "C", "D"],
        "chr": ["1", "1", "1", "2"],
        "start": [100, 200, 300, 150],
        "end": [200, 300, 400, 250],
    })

    plan, joints = build_connected_overlap_plan(
        loci,
        maximum_joint_region_kb=None,
        overlap_group_prefix="group_",
        provenance_delimiter=";",
    )

    chromosome_one = plan[plan["chromosome"].eq("1")]
    assert chromosome_one["overlap_group_id"].nunique() == 1
    assert chromosome_one["source_primary_genomic_loci"].unique().tolist() == [
        "A;B;C"
    ]
    assert joints.loc[0, ["START", "END"]].tolist() == [100, 400]
    assert plan.loc[plan["primary_genomic_locus"].eq("D"), "resolution_status"].item() == (
        "primary_retained"
    )


def test_joint_region_limit_excludes_group_before_rerun():
    loci = pd.DataFrame({
        "GenomicLocus": ["A", "B"],
        "chr": ["1", "1"],
        "start": [1, 900],
        "end": [1000, 2000],
    })

    plan, joints = build_connected_overlap_plan(
        loci,
        maximum_joint_region_kb=1,
        overlap_group_prefix="group_",
        provenance_delimiter=";",
    )

    assert not joints["eligible_for_rerun"].item()
    assert set(plan["resolution_status"]) == {"excluded_before_joint_rerun"}


def test_converged_primary_locus_without_a_credible_set_still_enters_overlap_plan(
    tmp_path,
):
    result = _write_handoff(
        tmp_path,
        [("with_set", "chr1:100-200")],
        "primary",
    )
    pd.DataFrame([
        {
            "genomic_locus": "with_set",
            "locus_chr": "1",
            "locus_start": 100,
            "locus_end": 200,
            "converged": True,
        },
        {
            "genomic_locus": "converged_without_set",
            "locus_chr": "1",
            "locus_start": 150,
            "locus_end": 250,
            "converged": True,
        },
    ]).to_csv(result["locus_status"], sep="\t", index=False)

    loci = _read_primary_result_loci(result)
    plan, _ = build_connected_overlap_plan(
        loci,
        maximum_joint_region_kb=None,
        overlap_group_prefix="group_",
        provenance_delimiter=";",
    )

    assert set(plan["primary_genomic_locus"]) == {
        "with_set", "converged_without_set",
    }
    assert plan["overlap_group_id"].ne("").all()


def test_all_no_set_primary_loci_can_produce_an_authoritative_joint_result(tmp_path):
    primary = _write_handoff(
        tmp_path,
        [("no_set_A", "chr1:100-200"), ("no_set_B", "chr1:150-250")],
        "primary",
    )
    pd.DataFrame(columns=["Filename", "GenomicLocus", "Annotfiles"]).to_csv(
        primary["flames_index"], sep="\t", index=False
    )
    Path(primary["genomic_loci"]).unlink()
    pd.DataFrame([
        {
            "genomic_locus": locus,
            "locus_chr": "1",
            "locus_start": start,
            "locus_end": end,
            "converged": True,
        }
        for locus, start, end in (
            ("no_set_A", 100, 200),
            ("no_set_B", 150, 250),
        )
    ]).to_csv(primary["locus_status"], sep="\t", index=False)
    settings = (
        load_configuration()
        .modules.fine_mapping.overlap_resolution.model_dump()
    )
    args = _overlap_args(settings, tmp_path)

    def fake_joint_runner(joint_args):
        locus = pd.read_csv(joint_args.locus_file, sep="\t").iloc[0]
        return _write_handoff(
            Path(joint_args.output_directory),
            [(locus["GenomicLocus"], "chr1:100-250")],
            "joint",
        )

    result = resolve_overlapping_results(args, primary, fake_joint_runner)
    combined = pd.read_csv(result["final_combined_credible_sets"], sep="\t")

    assert combined["analysis_round"].tolist() == ["joint"]
    assert combined["source_primary_genomic_loci"].item() == "no_set_A;no_set_B"
    assert pd.isna(combined["source_primary_credible_set_files"].item())


def test_non_overlapping_no_set_result_remains_a_valid_empty_outcome(tmp_path):
    primary = _write_handoff(
        tmp_path,
        [("no_set", "chr1:100-200")],
        "primary",
    )
    pd.DataFrame(columns=["Filename", "GenomicLocus", "Annotfiles"]).to_csv(
        primary["flames_index"], sep="\t", index=False
    )
    pd.DataFrame([{
        "genomic_locus": "no_set",
        "locus_chr": "1",
        "locus_start": 100,
        "locus_end": 200,
        "converged": True,
    }]).to_csv(primary["locus_status"], sep="\t", index=False)
    settings = (
        load_configuration()
        .modules.fine_mapping.overlap_resolution.model_dump()
    )
    args = _overlap_args(settings, tmp_path)

    result = resolve_overlapping_results(
        args,
        primary,
        lambda _: (_ for _ in ()).throw(AssertionError("rerun was not required")),
    )

    assert result["status"] == "completed_no_credible_sets"
    assert result["flames_input"] is None
    assert result["n_final_credible_sets"] == 0
    assert Path(result["final_combined_credible_sets"]).is_file()
    output_paths = resolve_output_paths(
        tmp_path,
        args.fine_mapping_output_layout,
        args.dataset_id,
    )
    assert not output_paths["downstream_flames_directory"].exists()
    assert not output_paths["prepared_loci_directory"].exists()


def test_overlap_output_directories_cannot_escape_the_run_directory(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "output_layout:\n  downstream_flames_directory: ..\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="must be"):
        _resolve_fine_mapping_arguments(
            argparse.Namespace(finemap_method="susie", run_config=str(config_file))
        )


def test_final_directory_cannot_replace_primary_handoff(tmp_path):
    primary = _write_handoff(
        tmp_path,
        [("primary", "chr1:100-200")],
        "primary",
    )
    settings = (
        load_configuration()
        .modules.fine_mapping.overlap_resolution.model_dump()
    )
    module = load_configuration().modules.fine_mapping
    layout = module.output_layout.model_dump()
    layout["downstream_flames_directory"] = "flames_input"
    layout["downstream_flames_annotations_directory"] = (
        "flames_input/annotations"
    )
    args = _overlap_args(settings, tmp_path)
    args.fine_mapping_output_layout = layout

    with pytest.raises(ValueError, match="cannot replace primary"):
        resolve_overlapping_results(args, primary, lambda _: None)

    assert Path(primary["flames_index"]).is_file()


def test_failed_joint_rerun_retains_exact_variant_limit_reason(tmp_path):
    primary = _write_handoff(
        tmp_path,
        [("primary_A", "chr1:100-200"), ("primary_B", "chr1:150-250")],
        "primary",
    )
    settings = (
        load_configuration()
        .modules.fine_mapping.overlap_resolution.model_dump()
    )
    args = _overlap_args(settings, tmp_path)

    def failed_joint_runner(joint_args):
        joint_locus = pd.read_csv(joint_args.locus_file, sep="\t").iloc[0]
        joint_output = Path(joint_args.output_directory)
        joint_output.mkdir(parents=True)
        joint_paths = resolve_output_paths(
            joint_output,
            joint_args.fine_mapping_output_layout,
            joint_args.dataset_id,
        )
        joint_paths["quality_control_directory"].mkdir(parents=True)
        pd.DataFrame([{
            "genomic_locus": joint_locus["GenomicLocus"],
            "status": "resource_limit_exceeded",
            "warning_reason": "",
            "failure_reason": "maximum_variants_per_locus_exceeded",
        }]).to_csv(
            joint_paths["susie_qc_file"],
            sep="\t",
            index=False,
        )
        raise RuntimeError("joint engine reported no eligible loci")

    with pytest.raises(RuntimeError, match="no authoritative credible sets"):
        resolve_overlapping_results(args, primary, failed_joint_runner)

    output_paths = resolve_output_paths(
        tmp_path,
        args.fine_mapping_output_layout,
        args.dataset_id,
    )
    plan = pd.read_csv(
        output_paths["quality_control_directory"] / settings["plan_filename"],
        sep="\t",
    )
    assert set(plan["resolution_status"]) == {"joint_rerun_failed_excluded"}
    assert set(plan["joint_failure_reason"]) == {
        "maximum_variants_per_locus_exceeded"
    }
    assert set(plan["resolution_reason"]) == {
        "maximum_variants_per_locus_exceeded"
    }


def test_final_combined_table_retains_primary_and_joint_provenance(tmp_path):
    primary = _write_handoff(
        tmp_path,
        [
            ("primary_A", "chr1:100-200"),
            ("primary_B", "chr1:150-250"),
            ("primary_C", "chr1:500-600"),
        ],
        "primary",
    )
    primary_cache_manifest = tmp_path / "chromosome_cache_manifest.tsv"
    primary["summary_statistics_chromosome_manifest"] = str(
        primary_cache_manifest
    )
    settings = (
        load_configuration()
        .modules.fine_mapping.overlap_resolution.model_dump()
    )
    args = _overlap_args(settings, tmp_path)
    calls = []

    def fake_joint_runner(joint_args):
        calls.append(joint_args)
        joint_loci = pd.read_csv(joint_args.locus_file, sep="\t")
        assert joint_args.window_kb == 0
        assert joint_args.locus_type == "range"
        assert joint_args.summary_statistics_chromosome_manifest == str(
            primary_cache_manifest
        )
        assert joint_loci[["START", "END"]].values.tolist() == [[100, 250]]
        return _write_handoff(
            Path(joint_args.output_directory),
            [(joint_loci.loc[0, "GenomicLocus"], "chr1:100-250")],
            "joint",
        )

    result = resolve_overlapping_results(args, primary, fake_joint_runner)

    assert len(calls) == 1
    combined = pd.read_csv(result["final_combined_credible_sets"], sep="\t")
    assert set(combined["analysis_round"]) == {"primary", "joint"}
    primary_row = combined[combined["analysis_round"].eq("primary")].iloc[0]
    joint_row = combined[combined["analysis_round"].eq("joint")].iloc[0]
    assert primary_row["source_primary_genomic_loci"] == "primary_C"
    assert joint_row["source_primary_genomic_loci"] == "primary_A;primary_B"
    assert "primary_1.txt" in joint_row["source_primary_credible_set_files"]
    assert "primary_2.txt" in joint_row["source_primary_credible_set_files"]
    assert "primary_3.txt" not in joint_row["source_primary_credible_set_files"]
    assert joint_row["warning_reason"] == "joint_warning_1"
    assert (
        joint_row["source_primary_warning_reasons"]
        == "primary_A=primary_warning_1;primary_B=primary_warning_2"
    )
    assert Path(result["overlap_resolution_configuration"]).is_file()
    assert Path(result["overlap_resolution_summary"]).is_file()

    final_index = pd.read_csv(
        Path(result["flames_input"]) / settings["index_filename"], sep="\t"
    )
    assert len(final_index) == 2
    assert set(final_index["GenomicLocus"]) == {"primary_C", "chr1:100-250"}
    assert not any("primary_1.txt" in value for value in final_index["Filename"])
    assert not any("primary_2.txt" in value for value in final_index["Filename"])

    plan = pd.read_csv(result["overlap_resolution"], sep="\t")
    replaced = plan[plan["resolution_status"].eq("replaced_by_joint_rerun")]
    assert set(replaced["primary_genomic_locus"]) == {"primary_A", "primary_B"}
