"""Regression coverage for fine-mapping input and resource preflight."""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from postgwas.core.checkpointing import checkpoint_namespace
from postgwas.core.errors import MissingRequiredArgumentsError
from postgwas.modules.fine_mapping.preflight import (
    FineMappingPreflightError,
    FineMappingResourcePreflight,
    run_fine_mapping_preflight,
    run_fine_mapping_resource_preflight,
)
from postgwas.modules.fine_mapping.engines.susie.adapter import split_locus_file
from postgwas.modules.fine_mapping.service import (
    _require_fine_mapping_inputs,
    _resolve_fine_mapping_arguments,
    preflight_fine_mapping_pipeline,
    run_fine_mapping,
)
from postgwas.pipeline.registry import REGISTRY
from preflight_support import pipeline_input_vcf_evidence


def _executable(path: Path) -> Path:
    path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    path.chmod(0o755)
    return path


def _plink_reference(tmp_path: Path, *, bed_size_delta: int = 0) -> Path:
    prefix = tmp_path / "reference"
    prefix.with_suffix(".bim").write_text(
        "1 rs1 0 100 A G\n1 rs2 0 200 C T\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".fam").write_text(
        "F1 I1 0 0 0 -9\nF2 I2 0 0 0 -9\n",
        encoding="utf-8",
    )
    # Two samples require one genotype byte per variant in SNP-major BED.
    prefix.with_suffix(".bed").write_bytes(
        bytes((0x6C, 0x1B, 0x01)) + b"\x00" * (2 + bed_size_delta)
    )
    return prefix


def _susie_inputs(tmp_path: Path) -> tuple[Path, Path]:
    summary = tmp_path / "study_susie.tsv"
    summary.write_text(
        "SNP\tCHR\tBP\tREF\tALT\tEZ\tNEF\tLP\n"
        "rs1\t1\t100\tG\tA\t3.0\t1000\t8.0\n"
        "rs2\t1\t200\tT\tC\t2.0\t1000\t7.5\n",
        encoding="utf-8",
    )
    loci = tmp_path / "loci.tsv"
    loci.write_text(
        "CHR\tSTART\tEND\tLP\tGenomicLocus\n"
        "1\t90\t210\t8.0\tchr1:90-210\n",
        encoding="utf-8",
    )
    return summary, loci


def _resolved_susie_args(
    tmp_path: Path,
    *,
    summary: Path | None = None,
    reference: Path | None = None,
):
    default_summary, loci = _susie_inputs(tmp_path)
    return _resolve_fine_mapping_arguments(
        argparse.Namespace(
            finemap_method="susie",
            susie_input_file=summary or default_summary,
            locus_file=loci,
            finemap_ld_reference=reference or _plink_reference(tmp_path),
            plink=_executable(tmp_path / "plink"),
        )
    )


def test_susie_preflight_validates_inputs_reference_and_runtime(
    tmp_path, monkeypatch,
):
    args = _resolved_susie_args(tmp_path)
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._tool_version",
        lambda *_args, **_kwargs: "PLINK v1.90",
    )
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight.resolve_susie_r_runtime",
        lambda _rscript, _timeout: {
            "rscript": str(_executable(tmp_path / "Rscript")),
            "version": "R version 4.4.0; susieR available",
            "library_paths": [str(tmp_path / "R-library")],
            "environment": {"R_LIBS_USER": str(tmp_path / "R-library")},
        },
    )

    preflight = run_fine_mapping_preflight(args)

    assert preflight.inputs.summary_statistics_rows == 2
    assert preflight.inputs.input_loci == 1
    assert preflight.inputs.eligible_loci == 1
    assert preflight.reference.variants == 2
    assert preflight.reference.samples == 2
    assert preflight.reference.bed_bytes == 5
    assert preflight.reference.coordinate_concordant_variant_matches == 2
    assert preflight.reference.loci_with_reference_variants == 1
    assert [tool.name for tool in preflight.tools] == ["PLINK", "R / susieR"]
    assert "File-structure and value validation" in {
        label for _, label, _ in preflight.input_screen_fields()
    }
    assert "Reference and runtime validation" in {
        label for _, label, _ in preflight.resource_screen_fields()
    }
    assert {row["category"] for row in preflight.audit_rows()} == {
        "input", "reference", "runtime",
    }


def test_pipeline_preflight_validates_resources_before_generated_inputs(
    tmp_path, monkeypatch,
):
    reference = _plink_reference(tmp_path)
    plink = _executable(tmp_path / "plink")
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._tool_version",
        lambda *_args, **_kwargs: "PLINK v1.90",
    )
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight.resolve_susie_r_runtime",
        lambda _rscript, _timeout: {
            "rscript": str(_executable(tmp_path / "Rscript")),
            "version": "R version 4.4.0; susieR available",
            "library_paths": [str(tmp_path / "R-library")],
            "environment": {"R_LIBS_USER": str(tmp_path / "R-library")},
        },
    )

    args = argparse.Namespace(
        finemap_method="susie",
        finemap_ld_reference=reference,
        plink=plink,
    )
    evidence = preflight_fine_mapping_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )

    assert REGISTRY.get("finemap").preflight.endswith(
        ":preflight_fine_mapping_pipeline"
    )
    assert evidence.module == "finemap"
    assert isinstance(evidence.resources, FineMappingResourcePreflight)
    assert evidence.resources.reference_prefix == str(reference.resolve())
    assert "formatter summary statistics" in evidence.deferred_checks[0]
    assert not hasattr(args, "susie_r_environment")
    assert not hasattr(args, "_susie_r_environment")


def test_full_preflight_reuses_tools_and_rejects_changed_reference(
    tmp_path, monkeypatch,
):
    args = _resolved_susie_args(tmp_path)
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._tool_version",
        lambda *_args, **_kwargs: "PLINK v1.90",
    )
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight.resolve_susie_r_runtime",
        lambda _rscript, _timeout: {
            "rscript": str(_executable(tmp_path / "Rscript")),
            "version": "R version 4.4.0; susieR available",
            "library_paths": [str(tmp_path / "R-library")],
            "environment": {"R_LIBS_USER": str(tmp_path / "R-library")},
        },
    )
    resources = run_fine_mapping_resource_preflight(args)
    Path(str(args.finemap_ld_reference) + ".fam").write_text(
        "F1 I1 0 0 0 -9\nF2 I2 0 0 0 -9\nF3 I3 0 0 0 -9\n",
        encoding="utf-8",
    )

    with pytest.raises(FineMappingPreflightError, match="changed after pipeline"):
        run_fine_mapping_preflight(args, resource_preflight=resources)


def test_full_preflight_restores_cached_runtime_without_version_subprocesses(
    tmp_path, monkeypatch,
):
    args = _resolved_susie_args(tmp_path)
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._tool_version",
        lambda *_args, **_kwargs: "PLINK v1.90",
    )
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight.resolve_susie_r_runtime",
        lambda _rscript, _timeout: {
            "rscript": str(_executable(tmp_path / "Rscript")),
            "version": "R version 4.4.0; susieR available",
            "library_paths": [str(tmp_path / "R-library")],
            "environment": {"R_LIBS_USER": str(tmp_path / "R-library")},
        },
    )
    resources = run_fine_mapping_resource_preflight(args)
    args.plink = "unresolved-plink"
    args.rscript = "unresolved-rscript"
    args._susie_r_environment = {}
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._tool_version",
        lambda *_args, **_kwargs: pytest.fail("tool version check was repeated"),
    )
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight.resolve_susie_r_runtime",
        lambda *_args, **_kwargs: pytest.fail("R runtime check was repeated"),
    )

    preflight = run_fine_mapping_preflight(
        args,
        resource_preflight=resources,
    )

    assert args.plink == resources.tools[0].path
    assert args.rscript == resources.r_runtime["rscript"]
    assert args._susie_r_environment == resources.r_runtime["environment"]
    assert not hasattr(args, "susie_r_environment")
    assert "_susie_r_environment" not in checkpoint_namespace(args)
    assert preflight.tools == resources.tools


def test_preflight_rejects_invalid_summary_values_before_runtime(
    tmp_path, monkeypatch,
):
    args = _resolved_susie_args(tmp_path)
    summary = Path(args.susie_input_file)
    summary.write_text(
        "SNP\tCHR\tBP\tREF\tALT\tEZ\tNEF\tLP\n"
        "rs1\t1\t100\tG\tA\t3.0\t0\t8.0\n",
        encoding="utf-8",
    )
    runtime_called = False

    def unexpected_runtime(*_args, **_kwargs):
        nonlocal runtime_called
        runtime_called = True

    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight.resolve_susie_r_runtime",
        unexpected_runtime,
    )

    with pytest.raises(FineMappingPreflightError, match="invalid row"):
        run_fine_mapping_preflight(args)

    assert runtime_called is False


def test_preflight_rejects_bed_bim_fam_dimension_mismatch(tmp_path, monkeypatch):
    reference = _plink_reference(tmp_path, bed_size_delta=1)
    args = _resolved_susie_args(tmp_path, reference=reference)
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._tool_version",
        lambda *_args, **_kwargs: "PLINK v1.90",
    )

    with pytest.raises(FineMappingPreflightError, match="BED size is inconsistent"):
        run_fine_mapping_preflight(args)


def test_preflight_records_loci_excluded_by_the_lp_threshold(
    tmp_path, monkeypatch,
):
    args = _resolved_susie_args(tmp_path)
    args.window_kb = 0
    Path(args.locus_file).write_text(
        "CHR\tSTART\tEND\tLP\tGenomicLocus\n"
        "1\t90\t110\t8.0\tchr1:90-110\n"
        "1\t190\t210\t6.0\tchr1:190-210\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._tool_version",
        lambda *_args, **_kwargs: "PLINK v1.90",
    )
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight.resolve_susie_r_runtime",
        lambda _rscript, _timeout: {
            "rscript": str(_executable(tmp_path / "Rscript")),
            "version": "R version 4.4.0; susieR available",
            "library_paths": [str(tmp_path / "R-library")],
            "environment": {"R_LIBS_USER": str(tmp_path / "R-library")},
        },
    )

    preflight = run_fine_mapping_preflight(args)

    assert preflight.inputs.input_loci == 2
    assert preflight.inputs.eligible_loci == 1
    assert preflight.inputs.loci_below_lp_threshold == 1
    assert preflight.inputs.locus_summary_variants == 1


def test_preflight_reports_partial_locus_reference_coverage(
    tmp_path, monkeypatch,
):
    args = _resolved_susie_args(tmp_path)
    summary = Path(args.susie_input_file)
    loci = Path(args.locus_file)
    summary.write_text(
        "SNP\tCHR\tBP\tREF\tALT\tEZ\tNEF\tLP\n"
        "rs1\t1\t100\tG\tA\t3.0\t1000\t8.0\n"
        "rs3\t1\t1000\tT\tC\t2.0\t1000\t7.5\n",
        encoding="utf-8",
    )
    loci.write_text(
        "CHR\tSTART\tEND\tLP\tGenomicLocus\n"
        "1\t90\t110\t8.0\tchr1:90-110\n"
        "1\t990\t1010\t7.5\tchr1:990-1010\n",
        encoding="utf-8",
    )
    args.window_kb = 0
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._tool_version",
        lambda *_args, **_kwargs: "PLINK v1.90",
    )
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight.resolve_susie_r_runtime",
        lambda _rscript, _timeout: {
            "rscript": str(_executable(tmp_path / "Rscript")),
            "version": "R version 4.4.0; susieR available",
            "library_paths": [str(tmp_path / "R-library")],
            "environment": {"R_LIBS_USER": str(tmp_path / "R-library")},
        },
    )

    preflight = run_fine_mapping_preflight(args)

    assert preflight.reference.locus_variant_id_matches == 1
    assert preflight.reference.coordinate_concordant_variant_matches == 1
    assert preflight.reference.loci_with_reference_variants == 1
    assert preflight.reference.loci_requested == 2
    locus_field = next(
        field
        for field in preflight.resource_screen_fields()
        if field[1] == "Eligible loci with reference variants"
    )
    assert locus_field == ("warning", "Eligible loci with reference variants", "1/2")
    overlap_audit = next(
        row
        for row in preflight.audit_rows()
        if row["check"] == "summary_reference_overlap"
    )
    assert overlap_audit["status"] == "warning"


def test_susie_execution_uses_normalized_chromosomes_for_mhc_exclusion(tmp_path):
    loci = tmp_path / "loci.tsv"
    loci.write_text(
        "CHR\tSTART\tEND\tLP\n"
        "chr6\t26000000\t27000000\t8.0\n"
        "chr1\t90\t110\t8.0\n",
        encoding="utf-8",
    )

    chunks = split_locus_file(
        str(loci),
        tmp_path / "chunks",
        2,
        finemap_skip_mhc=True,
        mhc_chr=6,
        mhc_start=25000000,
        mhc_end=35000000,
    )

    assert len(chunks) == 1
    retained = chunks[0].read_text(encoding="utf-8")
    assert "chr1" in retained
    assert "chr6" not in retained


def test_finemap_toolchain_comes_from_resolved_resource_yaml(tmp_path, monkeypatch):
    reference = _plink_reference(tmp_path)
    summary = tmp_path / "study_finemap.tsv"
    summary.write_text(
        "rsid\tchromosome\tposition\tallele1\tallele2\tmaf\tbeta\tse\tNEF\n"
        "rs1\t1\t100\tA\tG\t0.1\t0.2\t0.1\t1000\n",
        encoding="utf-8",
    )
    loci = tmp_path / "loci.tsv"
    loci.write_text(
        "CHR\tSTART\tEND\tLP\n1\t90\t110\t8.0\n",
        encoding="utf-8",
    )
    executables = {
        name: _executable(tmp_path / name)
        for name in ("plink2", "bgenix", "ldstore", "finemap")
    }
    config = tmp_path / "run.yaml"
    config.write_text(
        "config_version: 1\n"
        "resources:\n"
        "  executables:\n"
        + "".join(
            "    %s: %s\n" % (name, path)
            for name, path in executables.items()
        ),
        encoding="utf-8",
    )
    args = _resolve_fine_mapping_arguments(
        argparse.Namespace(
            finemap_method="finemap",
            run_config=str(config),
            finemap_in_files=summary,
            locus_file=loci,
            finemap_ld_reference=reference,
        )
    )
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._tool_version",
        lambda executable, *_args, **_kwargs: (
            "PLINK v2.00"
            if Path(executable).name == "plink2"
            else "%s version" % Path(executable).name
        ),
    )

    preflight = run_fine_mapping_preflight(args)

    assert {tool.path for tool in preflight.tools} == {
        str(path.resolve()) for path in executables.values()
    }
    assert args.bgenix == str(executables["bgenix"].resolve())
    assert args.ldstore == str(executables["ldstore"].resolve())
    assert args.finemap_executable == str(executables["finemap"].resolve())


def test_required_fine_mapping_inputs_are_reported_together():
    args = _resolve_fine_mapping_arguments(
        argparse.Namespace(finemap_method="susie")
    )

    with pytest.raises(MissingRequiredArgumentsError) as caught:
        _require_fine_mapping_inputs(args)

    message = str(caught.value)
    assert "--locus-file" in message
    assert "--susie-input-file" in message
    assert "--finemap-ld-reference" in message
    assert "modules.fine_mapping.input.locus_file" in message


def test_service_displays_preflight_before_engine_and_publishes_qc(
    tmp_path, monkeypatch, capsys,
):
    summary, loci = _susie_inputs(tmp_path)
    reference = _plink_reference(tmp_path)
    plink = _executable(tmp_path / "plink")
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._tool_version",
        lambda *_args, **_kwargs: "PLINK v1.90",
    )
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight.resolve_susie_r_runtime",
        lambda _rscript, _timeout: {
            "rscript": str(_executable(tmp_path / "Rscript")),
            "version": "R version 4.4.0; susieR available",
            "library_paths": [str(tmp_path / "R-library")],
            "environment": {"R_LIBS_USER": str(tmp_path / "R-library")},
        },
    )

    def fake_engine(args, screen=None):
        assert args._fine_mapping_preflight.inputs.summary_statistics_rows == 2
        for stage in range(3, 8):
            screen.complete(stage, [("success", "Synthetic stage", "passed")])
        return {
            "status": "success",
            "output_dir": str(Path(args.output_directory).resolve()),
            "flames_input": str(Path(args.output_directory) / "flames"),
            "locus_status": str(Path(args.output_directory) / "locus.tsv"),
            "n_attempted": 1,
            "n_successful": 1,
            "n_failed": 0,
            "n_warnings": 0,
            "warning_reason_counts": {},
            "failure_reason_counts": {},
            "n_credible_sets": 1,
        }

    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.engines.susie.adapter.run_parallel_susie",
        fake_engine,
    )
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.service.resolve_overlapping_results",
        lambda _args, primary, _engine: {
            **primary,
            "status": "success",
            "n_overlap_groups": 0,
            "n_joint_rerun_successful": 0,
            "n_final_credible_sets": 1,
            "final_combined_credible_sets": str(
                Path(primary["output_dir"]) / "final.tsv"
            ),
            "overlap_resolution_summary": str(
                Path(primary["output_dir"]) / "overlap.tsv"
            ),
        },
    )
    report_arguments = {}

    def fake_report(_args, _preflight, _result, paths, **kwargs):
        report_arguments.update(kwargs)
        return paths["html_report_file"]

    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.service.write_fine_mapping_html_report",
        fake_report,
    )
    upstream = {"ld_clump": {"html_report": str(tmp_path / "clumping.html")}}
    output = tmp_path / "output"
    result = run_fine_mapping(
        argparse.Namespace(
            finemap_method="susie",
            susie_input_file=summary,
            locus_file=loci,
            finemap_ld_reference=reference,
            plink=plink,
            dataset_id="STUDY",
            output_directory=output,
        ),
        upstream_results=upstream,
    )

    report = Path(result["preflight_validation"])
    assert report_arguments["upstream_results"] == upstream
    assert report.is_file()
    assert result["html_report"].endswith("STUDY_fine_mapping_report.html")
    assert "summary_statistics" in report.read_text(encoding="utf-8")
    terminal = capsys.readouterr().out
    assert "Validate summary statistics and locus definitions" in terminal
    assert "Summary-statistic variants" in terminal
    assert "Locus variant IDs found in BIM" in terminal
    assert "Reference and runtime validation" in terminal
    assert "Scientific HTML report" in terminal


def test_service_persists_failed_preflight_without_starting_the_engine(
    tmp_path, monkeypatch,
):
    summary, loci = _susie_inputs(tmp_path)
    summary.write_text(
        "SNP\tCHR\tBP\tREF\tALT\tEZ\tNEF\tLP\n"
        "rs1\t1\t100\tG\tA\t3.0\t0\t8.0\n",
        encoding="utf-8",
    )
    engine_called = False

    def unexpected_engine(*_args, **_kwargs):
        nonlocal engine_called
        engine_called = True

    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.engines.susie.adapter.run_parallel_susie",
        unexpected_engine,
    )
    output = tmp_path / "output"

    with pytest.raises(FineMappingPreflightError, match="invalid row"):
        run_fine_mapping(
            argparse.Namespace(
                finemap_method="susie",
                susie_input_file=summary,
                locus_file=loci,
                finemap_ld_reference=_plink_reference(tmp_path),
                plink=_executable(tmp_path / "plink"),
                dataset_id="STUDY",
                output_directory=output,
            )
        )

    report = output / "quality_control" / "input_and_resource_validation.tsv"
    log = output / "run_metadata" / "pipeline_summary.log"
    assert engine_called is False
    assert "\tfailed\t" in report.read_text(encoding="utf-8")
    assert "Fine-mapping preflight failed" in log.read_text(encoding="utf-8")
