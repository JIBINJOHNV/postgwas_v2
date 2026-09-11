"""Regression tests for validated PRED-LD execution boundaries."""

from argparse import Namespace
from pathlib import Path
import sys

import pytest

from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.validation_reporting import FileValidationDisplay
from postgwas.modules.imputation import service
from postgwas.modules.imputation.engines.pred_ld.pred_ld_runner import (
    build_interleaved_chr_order,
    process_pred_ld_results_all_parallel,
    run_pred_ld_parallel,
    validate_pred_ld_reference,
)
from postgwas.modules.imputation.service import (
    PredLDPipelineResources,
    preflight_imputation_pipeline,
    run_sumstat_imputation_direct,
)
from preflight_support import pipeline_input_vcf_evidence


def _write_pred_ld_reference(root: Path) -> tuple[Path, tuple[Path, ...]]:
    configuration = load_configuration()
    pred_ld = configuration.modules.imputation.engines.pred_ld
    reference = root / "pred_ld_reference"
    resource_directory = reference / pred_ld.reference_subdirectory_template.format(
        mode=pred_ld.mode,
        population="EUR",
    )
    resource_directory.mkdir(parents=True)
    paths = []
    for chromosome in configuration.modules.formatting.chromosomes:
        for kind in pred_ld.reference_file_kinds:
            path = resource_directory / pred_ld.reference_file_template.format(
                population="EUR",
                chromosome=chromosome,
                kind=kind,
            )
            path.write_bytes(b"reference\n")
            paths.append(path.resolve())
    return reference, tuple(paths)


def _pipeline_args(tmp_path: Path, reference: Path, **values) -> Namespace:
    supplied = {
        "imputation_ld_reference": str(reference),
        "resource_directory": str(tmp_path / "postgwas_resources"),
        "genome_build": "GRCh37",
        "population": "EUR",
        "dataset_id": "STUDY",
        "output_directory": str(tmp_path / "output"),
        "threads": 1,
    }
    supplied.update(values)
    return Namespace(**supplied)


def test_validate_pred_ld_reference_enumerates_exact_supported_files(
    tmp_path, capsys,
):
    reference, expected = _write_pred_ld_reference(tmp_path)
    pred_ld = load_configuration().modules.imputation.engines.pred_ld

    configuration = load_configuration()
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, configuration)
        validated_root, files = validate_pred_ld_reference(
            reference,
            population="EUR",
            chromosomes=configuration.modules.formatting.chromosomes,
            mode="TOP_LD",
            subdirectory_template=pred_ld.reference_subdirectory_template,
            file_template=pred_ld.reference_file_template,
            file_kinds=pred_ld.reference_file_kinds,
        )
        display.flush()

    assert validated_root == reference.resolve()
    assert files == expected
    screen = capsys.readouterr().out
    assert (
        "PRED-LD chromosome reference resources — AVAILABLE — availability only"
        in screen
    )
    assert "Chromosome reference files" in screen
    assert all(path.name not in screen for path in expected)


def test_pipeline_preflight_validates_pred_ld_and_harmonisation_resources(
    tmp_path, monkeypatch,
):
    reference, expected = _write_pred_ld_reference(tmp_path)
    resource_root = tmp_path / "postgwas_resources"
    resource_root.mkdir()
    args = _pipeline_args(tmp_path, reference)
    monkeypatch.setattr(
        service,
        "require_binaries",
        lambda *_args, **_kwargs: {"python": sys.executable},
    )

    evidence = preflight_imputation_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )

    resources = evidence.resources
    assert isinstance(resources, PredLDPipelineResources)
    assert resources.reference_files == expected
    assert resources.resource_root == resource_root.resolve()
    assert resources.configuration.modules.imputation.genome_build.value == "GRCh37"
    identity_paths = {identity.path for identity in resources.file_identities}
    assert set(expected).issubset(identity_paths)
    assert resources.pred_ld_script in identity_paths
    assert Path(sys.executable).resolve() in identity_paths


def test_pipeline_imputation_rejects_reference_changed_after_preflight(
    tmp_path, monkeypatch,
):
    reference, expected = _write_pred_ld_reference(tmp_path)
    (tmp_path / "postgwas_resources").mkdir()
    formatted = tmp_path / "formatted"
    formatted.mkdir()
    args = _pipeline_args(
        tmp_path,
        reference,
        pred_ld_input_directory=str(formatted),
    )
    monkeypatch.setattr(
        service,
        "require_binaries",
        lambda *_args, **_kwargs: {"python": sys.executable},
    )
    evidence = preflight_imputation_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )
    expected[0].write_bytes(b"reference changed after validation\n")

    with pytest.raises(RuntimeError, match="changed after pipeline preflight"):
        run_sumstat_imputation_direct(
            args,
            pipeline_resources=evidence.resources,
        )

    assert not Path(args.output_directory).exists()


def test_pipeline_imputation_build_mismatch_fails_before_resource_checks(
    tmp_path, monkeypatch,
):
    reference, _ = _write_pred_ld_reference(tmp_path)
    (tmp_path / "postgwas_resources").mkdir()
    args = _pipeline_args(tmp_path, reference, genome_build="GRCh38")
    monkeypatch.setattr(
        service,
        "require_binaries",
        lambda *_args, **_kwargs: pytest.fail(
            "software validation ran after a scientific build mismatch"
        ),
    )

    with pytest.raises(ValueError, match="declares genome build GRCh37"):
        preflight_imputation_pipeline(
            args,
            preflight_evidence=pipeline_input_vcf_evidence(),
        )


def test_pred_ld_runner_rejects_any_missing_configured_chromosome(tmp_path):
    formatted = tmp_path / "formatted"
    reference = tmp_path / "reference"
    output = tmp_path / "output"
    formatted.mkdir()
    reference.mkdir()

    pred_ld = load_configuration().modules.imputation.engines.pred_ld
    with pytest.raises(
        RuntimeError,
        match="did not complete every configured chromosome",
    ):
        run_pred_ld_parallel(
            predld_input_dir=str(formatted),
            output_folder=str(output),
            output_prefix="STUDY",
            pred_ld_ref=str(reference),
            chromosomes=("1",),
            r2threshold=0.8,
            maf=0.001,
            population="EUR",
            ref="TOP_LD",
            threads=1,
            memory_gb=load_configuration().execution.memory_gb,
            pred_ld_script=Path(__file__),
            python_executable=sys.executable,
            memory_gb_per_worker=pred_ld.memory_gb_per_worker,
            free_memory_threshold_gb=pred_ld.free_memory_threshold_gb,
            free_memory_threshold_fraction=(
                pred_ld.free_memory_threshold_fraction
            ),
            memory_poll_seconds=pred_ld.memory_poll_seconds,
            worker_poll_seconds=pred_ld.worker_poll_seconds,
            large_chromosomes=pred_ld.large_chromosomes,
            preferred_chromosome_order=pred_ld.preferred_chromosome_order,
        )

    combined_log = output / "STUDY_predld_combined.log"
    assert combined_log.is_file()
    assert "missing or empty for chr1" in combined_log.read_text(encoding="utf-8")
    assert not (output / ".predld_work").exists()


def test_configured_chromosome_schedule_never_drops_unlisted_input():
    assert build_interleaved_chr_order(
        ("chr2", "chrY", "1"), ("1", "2"),
    ) == ["1", "2", "Y"]


def test_imputation_budget_failure_precedes_reference_work(tmp_path, monkeypatch):
    (tmp_path / "postgwas_resources").mkdir()
    monkeypatch.setattr(
        service, "validate_pred_ld_reference",
        lambda *_args, **_kwargs: pytest.fail(
            "infeasible worker must not start reference work"
        ),
    )
    with pytest.raises(ValueError, match="cannot fit one worker"):
        preflight_imputation_pipeline(
            _pipeline_args(tmp_path, tmp_path / "unused", memory_gb=16),
            preflight_evidence=pipeline_input_vcf_evidence(),
        )


def test_pred_ld_runner_budget_failure_creates_no_outputs(tmp_path):
    pred_ld = load_configuration().modules.imputation.engines.pred_ld
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="cannot fit one worker"):
        run_pred_ld_parallel(
            predld_input_dir=str(tmp_path / "unused"),
            output_folder=str(output), output_prefix="STUDY",
            pred_ld_ref=str(tmp_path / "unused-reference"), chromosomes=("1",),
            r2threshold=pred_ld.minimum_r2, maf=pred_ld.minimum_maf,
            population="EUR", ref=pred_ld.mode, threads=3, memory_gb=16,
            pred_ld_script=Path(__file__), python_executable=sys.executable,
            memory_gb_per_worker=pred_ld.memory_gb_per_worker,
            free_memory_threshold_gb=pred_ld.free_memory_threshold_gb,
            free_memory_threshold_fraction=pred_ld.free_memory_threshold_fraction,
            memory_poll_seconds=pred_ld.memory_poll_seconds,
            worker_poll_seconds=pred_ld.worker_poll_seconds,
            large_chromosomes=pred_ld.large_chromosomes,
            preferred_chromosome_order=pred_ld.preferred_chromosome_order,
        )
    assert not output.exists()


def test_pred_ld_reference_mode_is_schema_restricted_to_supported_handoff():
    with pytest.raises(ConfigurationError, match="TOP_LD"):
        load_configuration(cli_overrides={
            "modules.imputation.engines.pred_ld.mode": "SNP_LD",
        })


def test_pred_ld_scheduler_controls_are_schema_validated():
    with pytest.raises(ConfigurationError, match="greater than 0"):
        load_configuration(cli_overrides={
            "modules.imputation.engines.pred_ld.memory_gb_per_worker": 0,
        })
    with pytest.raises(ConfigurationError, match="unique chromosome labels"):
        load_configuration(cli_overrides={
            "modules.imputation.engines.pred_ld.preferred_chromosome_order": [
                "1", "chr1",
            ],
        })
    with pytest.raises(ConfigurationError, match="exactly these placeholders"):
        load_configuration(cli_overrides={
            "modules.imputation.engines.pred_ld.reference_file_template": (
                "{population}_{chromosome}.csv.gz"
            ),
        })


def _write_postprocessing_inputs(root: Path) -> None:
    root.mkdir()
    rows = [
        ("rs1", 1, 1.0), ("rs1", 0, 1.0),
        ("rs2", 1, 2.0), ("rs2", 0, 2.0),
        ("rs3", 1, 10.0), ("rs3", 0, 3.0),
    ]
    result = root / "imputation_results_chr1.txt"
    result.write_text(
        "chr\tsnp\tA1\tA2\tpos\tbeta\tSE\tz\timputed\tR2\tNC\tSS\tAF\tLP\tSI\n"
        + "".join(
            "1\t%s\tA\tG\t%d\t%s\t0.1\t%s\t%d\t0.9\t500\t1000\t0.2\t2\t0.95\n"
            % (snp, index + 1, beta, beta, imputed)
            for index, (snp, imputed, beta) in enumerate(rows)
        ),
        encoding="utf-8",
    )
    information = root / "LD_info_TOP_LD_chr1.txt"
    information.write_text(
        "pos1\tpos2\tR2\tDprime\tALT_AF1\tALT_AF2\t+/-corr\t"
        "rsID1\trsID2\tREF1\tALT1\tREF2\tALT2\n"
        "1\t2\t0.9\t1\t0.2\t0.2\t+\trs1\trs1\tG\tA\tG\tA\n"
        "2\t3\t0.9\t1\t0.2\t0.2\t+\trs2\trs2\tG\tA\tG\tA\n"
        "3\t4\t0.9\t1\t0.2\t0.2\t+\trs3\trs3\tG\tA\tG\tA\n",
        encoding="utf-8",
    )


def test_postprocessing_uses_the_configured_correlation_method(tmp_path):
    observed = {}
    for method in ("pearson", "spearman"):
        root = tmp_path / method
        _write_postprocessing_inputs(root)
        _, correlations, sample_sheet = process_pred_ld_results_all_parallel(
            folder_path=str(root),
            output_path=str(root),
            output_prefix="STUDY",
            harmonised_dataset_id="STUDY_imputed",
            corr_method=method,
            threads=1,
        )
        observed[method] = correlations.row(0, named=True)["beta_corr"]
        assert Path(sample_sheet).is_file()
        assert (root / "STUDY_PREDLD_allchr.tsv.gz").is_file()
        assert (root / "STUDY_imputation_results.txt.gz").is_file()

    assert observed["pearson"] < 1.0
    assert observed["spearman"] == pytest.approx(1.0)


def test_postprocessing_records_chromosome_failure_before_stopping(tmp_path):
    root = tmp_path / "invalid"
    _write_postprocessing_inputs(root)
    (root / "LD_info_TOP_LD_chr1.txt").write_text(
        "invalid\ncontent\n",
        encoding="utf-8",
    )
    combined_log = root / "STUDY_predld_combined.log"
    combined_log.write_text("PRED-LD execution completed\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="post-processing failed for chromosome"):
        process_pred_ld_results_all_parallel(
            folder_path=str(root),
            output_path=str(root),
            output_prefix="STUDY",
            harmonised_dataset_id="STUDY_imputed",
            corr_method="pearson",
            threads=1,
        )

    log_text = combined_log.read_text(encoding="utf-8")
    assert "PRED-LD post-processing failed" in log_text
    assert "chr1:" in log_text
