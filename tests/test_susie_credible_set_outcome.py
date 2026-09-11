"""Regression tests for SuSiE credible-set outcome reporting."""

from __future__ import annotations

import argparse

import pandas as pd
import pytest

from postgwas.modules.fine_mapping.engines.susie.adapter import (
    SUSIE_STATUS_NO_CREDIBLE_SETS,
    SUSIE_STATUS_SUCCESS,
    _summarize_credible_set_handoff,
    _summarize_recovery_audit,
    _write_flames_index,
)
from postgwas.pipeline.runners import run_flames_runner


def test_no_worker_credible_sets_remove_stale_index_and_disable_handoff(tmp_path):
    worker = tmp_path / "worker_1"
    (worker / "flames_input").mkdir(parents=True)
    final = tmp_path / "final"
    final_flames = final / "flames_input"
    final_flames.mkdir(parents=True)
    stale_index = final_flames / "indexfile.txt"
    stale_index.write_text(
        "Filename\tGenomicLocus\tAnnotfiles\nold.txt\told\told.annot\n",
        encoding="utf-8",
    )

    index_file = _write_flames_index(
        [worker], final_flames, final_flames / "annotations", "indexfile.txt",
        "annotated_",
    )
    summary = _summarize_credible_set_handoff(
        index_file,
        {"chr1:100-200", "chr2:300-400"},
    )

    assert index_file is None
    assert not stale_index.exists()
    assert summary["status"] == SUSIE_STATUS_NO_CREDIBLE_SETS
    assert summary["flames_input"] is None
    assert summary["n_credible_sets"] == 0
    assert summary["n_loci_with_credible_sets"] == 0
    assert summary["n_converged_without_credible_sets"] == 2


def test_handoff_counts_sets_and_converged_loci_without_sets(tmp_path):
    worker = tmp_path / "worker_1"
    worker_flames = worker / "flames_input"
    worker_flames.mkdir(parents=True)
    credible_set = worker_flames / "study_CS_L1.txt"
    credible_set.write_text(
        "index cred1 prob1\n1 1:110:A_G 1.0\n",
        encoding="utf-8",
    )
    pd.DataFrame(
        [{"Filename": credible_set.name, "GenomicLocus": "chr1:100-200"}]
    ).to_csv(worker_flames / "indexfile_rows.tsv", sep="\t", index=False)
    final = tmp_path / "final"
    final_flames = final / "flames_input"
    final_flames.mkdir(parents=True)
    final_credible_set = final_flames / credible_set.name
    final_credible_set.write_text(
        credible_set.read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    index_file = _write_flames_index(
        [worker], final_flames, final_flames / "annotations", "indexfile.txt",
        "annotated_",
    )
    summary = _summarize_credible_set_handoff(
        index_file,
        {"chr1:100-200", "chr2:300-400"},
    )

    assert summary["status"] == SUSIE_STATUS_SUCCESS
    assert summary["flames_input"] == str(final_flames)
    assert summary["n_credible_sets"] == 1
    assert summary["n_loci_with_credible_sets"] == 1
    assert summary["n_converged_without_credible_sets"] == 1


def test_handoff_rejects_an_index_with_a_missing_credible_set(tmp_path):
    index_file = tmp_path / "indexfile.txt"
    pd.DataFrame(
        [
            {
                "Filename": "missing_CS_L1.txt",
                "GenomicLocus": "chr1:100-200",
                "Annotfiles": "annotated.txt",
            }
        ]
    ).to_csv(index_file, sep="\t", index=False)

    with pytest.raises(FileNotFoundError, match="missing or empty"):
        _summarize_credible_set_handoff(index_file, {"chr1:100-200"})


def test_pipeline_flames_runner_rejects_missing_credible_set_handoff(tmp_path):
    args = argparse.Namespace(output_directory=str(tmp_path))
    context = {
        "finemap": {
            "status": SUSIE_STATUS_NO_CREDIBLE_SETS,
            "flames_input": None,
        }
    }

    with pytest.raises(
        ValueError,
        match="requires at least one retained fine-mapping credible set",
    ):
        run_flames_runner(args, context)


def test_pipeline_flames_runner_uses_native_not_corrected_magmacovar_result(
    tmp_path, monkeypatch,
):
    args = argparse.Namespace(output_directory=str(tmp_path))
    context = {
        "finemap": {"status": SUSIE_STATUS_SUCCESS, "flames_input": "credible"},
        "magma": {
            "primary_mapping": "positional",
            "mapping_analyses": {
                "positional": {
                    "result_statistic_type": "calibrated_gene_p_value",
                },
            },
            "magma_genes_out": "study.genes.out",
        },
        "magma_covar": {
            "raw_results": "study.gsa.out",
            "corrected_results": "study_magmacovar_corrected.tsv",
        },
        "pops_output": "study.pops.tsv",
    }
    observed = {}

    def fake_flames(runtime_args):
        observed["magmacovar"] = runtime_args.magma_covariate_results_file
        return "completed"

    monkeypatch.setattr(
        "postgwas.modules.flames.service.run_flames_direct", fake_flames,
    )

    result = run_flames_runner(args, context)

    assert result == "completed"
    assert observed["magmacovar"] == "study.gsa.out"
    assert context["flames"] == "completed"
    assert args.output_directory == str(tmp_path)


def test_recovery_audit_summary_counts_scientific_outcomes(tmp_path):
    audit_file = tmp_path / "recovery.tsv"
    pd.DataFrame(
        [
            {
                "genomic_locus": "chr1:1-2",
                "stage": "ld_validation",
                "status": "failed",
                "reason": "ld_non_psd",
                "repairable": True,
            },
            {
                "genomic_locus": "chr1:1-2",
                "stage": "recovery_fit",
                "status": "success",
                "reason": "susie_converged",
                "repairable": False,
            },
            {
                "genomic_locus": "chr2:1-2",
                "stage": "ld_validation",
                "status": "failed",
                "reason": "ld_z_mismatch",
                "repairable": False,
            },
            {
                "genomic_locus": "chr3:1-2",
                "stage": "recovery_fit",
                "status": "failed",
                "reason": "susie_timeout",
                "repairable": False,
            },
        ]
    ).to_csv(audit_file, sep="\t", index=False)

    summary = _summarize_recovery_audit(audit_file)

    assert summary["n_recovery_successes"] == 1
    assert summary["n_recovery_failures"] == 1
    assert summary["n_fatal_ld_failures"] == 1
