"""Detailed reports reuse real preparation results and preserve early failures."""

from pathlib import Path
import json
import time
from types import SimpleNamespace

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.modules.harmonisation import service
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.study_properties import (
    StudyPropertyError, resolve_study_properties,
)


def test_early_property_failure_retains_completed_decision_evidence():
    frame = pl.DataFrame({"BETA": [-.1, .1, .2, -.2], "P": [-.1, .2, .3, .4]})
    with pytest.raises(StudyPropertyError) as caught:
        resolve_study_properties(frame, {
            "beta_or_col": "BETA", "pval_col": "P", "declared_effect_type": "beta",
        })
    error = caught.value
    assert "negative" in str(error).lower()
    assert error.study_decisions["effect_type"] == "beta"
    assert error.study_decisions["effect_evidence"]["non_null"] == 4
    assert error.study_decisions["pvalue_evidence"]["negative"] == 1
    assert frame["P"].to_list() == [-.1, .2, .3, .4]


def _prepare(tmp_path, monkeypatch, *, fail_frequency=False, failure_stage=None,
             dataset_id="study", automatic_inference=False):
    tmp_path.mkdir(parents=True, exist_ok=True)
    module = load_configuration().modules.harmonisation
    source = tmp_path / "study.tsv"
    pl.DataFrame({
        "CHR": ["1"] * 4, "POS": [10, 20, 30, 40], "EA": ["G"] * 4,
        "OA": ["A"] * 4, "BETA": [-.1, .1, .2, -.2], "SE": [.1] * 4,
        "P": [.2, .3, .4, .5], "AF": [.1, .2, .7, .8],
        "INFO": [100.] * 4 if failure_stage == "info" else [.9] * 4,
    }).write_csv(source, separator="\t")
    mapping = {
        "gwas_outputname": dataset_id, "sumstat_file": str(source),
        "output_folder": str(tmp_path / "out"), "chr_col": "CHR", "pos_col": "POS",
        "ea_col": "EA", "oa_col": "OA", "beta_or_col": "BETA", "se_col": "SE",
        "pval_col": "P", "eaf_col": "AF", "declared_effect_type": "beta",
        "declared_pvalue_type": "raw", "ncontrol": "1000",
        "imp_info_col": "INFO", "info_source": "internal",
    }
    if automatic_inference:
        mapping.pop("declared_effect_type")
        mapping.pop("declared_pvalue_type")
    resource_maps = {"1": {"default_eaf_file": str(tmp_path / "unused.tsv"),
        "default_eaf_reference_column": "EUR", "user_eaf_file": None,
        "user_info_file": None}}
    summary = dict.fromkeys(("unique_resource_files_checked", "indexed_vcfs_checked",
        "frequency_vcf_headers_checked", "fasta_indexes_checked", "table_headers_checked",
        "structured_text_files_checked"), 0)
    summary.update(chromosomes=["1"], file_identities=[], resources_by_chromosome=resource_maps)
    monkeypatch.setattr(service, "infer_genome_build", lambda *_a, **_k: {
        "inferred_build": "GRCh37", "input_variants": 4, "testable_variants": 4,
        "untestable_variants": 0,
        "matches": {"GRCh37": 4, "GRCh38": 0},
        "per_chromosome": {"1": {"assessment": "supports_selected_build",
            "testable_variants": 4, "matches": {"GRCh37": 4, "GRCh38": 0},
            "exclusive_snp_loci": {"GRCh37": 4, "GRCh38": 0}}},
    })
    monkeypatch.setattr(service, "resolve_strand_consensus", lambda *_a, **_k: {
        "strand": "forward", "strand_mode": "auto", "forward": 4, "reverse": 0,
        "ambiguous": 0, "informative": 4, "dominant_fraction": 1.,
        "minimum_informative_variants": 1, "consensus_threshold": .99,
        "reason": "known-forward reference fixture",
    })
    monkeypatch.setattr(service, "preflight_harmonisation_resources", lambda *_a, **_k: (resource_maps, summary))
    if failure_stage == "build":
        monkeypatch.setattr(service, "infer_genome_build", lambda *_a, **_k: {
            "inferred_build": "Ambiguous", "input_variants": 4, "testable_variants": 0,
            "untestable_variants": 4, "ambiguous_reason": "no informative build matches",
        })
    if failure_stage == "resources":
        def resources_fail(*_args, **_kwargs):
            raise service.ResourcePreflightError([{
                "resource": "FASTA", "chromosomes": ["1"],
                "path": str(tmp_path / "missing.fa"), "problem": "missing fixture reference",
            }])
        monkeypatch.setattr(service, "preflight_harmonisation_resources", resources_fail)
    if failure_stage == "replacement":
        def replacement_fail(*_args, **_kwargs):
            raise service.PipelineError("replacement reference fixture failed")
        monkeypatch.setattr(service, "_activate_eaf_replacement_resources", replacement_fail)
    if failure_stage == "staging":
        def staging_fail(*_args, **_kwargs):
            raise service.ExternalReferenceStagingError("staging fixture failure")
        monkeypatch.setattr(service, "stage_shared_external_reference_files", staging_fail)
    if failure_stage == "partition":
        def partition_fail(*_args, **_kwargs):
            raise OSError("partition fixture failure")
        monkeypatch.setattr(service, "write_chromosome_partitions", partition_fail)
    if failure_stage == "chromosome_set":
        from postgwas.modules.harmonisation.chromosome_partition import ChromosomePartitions
        writer = service.write_chromosome_partitions
        def mismatched_partitions(*args, **kwargs):
            parts = writer(*args, **kwargs)
            return ChromosomePartitions({"2": parts.paths["1"]}, {"2": parts.rows["1"]})
        monkeypatch.setattr(service, "write_chromosome_partitions", mismatched_partitions)
    if fail_frequency:
        from postgwas.modules.harmonisation.allele_frequency import StudyAlleleFrequencyError
        def fail(frame, columns, decisions, *_args, **_kwargs):
            raise StudyAlleleFrequencyError("reference fixture failure", study_decisions={
                **decisions, "eaf_reference_validation": {
                    "status": "failed", "comparable_variants": 3,
                    "minimum_overlap": 1000, "decision_reason": "insufficient_overlap",
                },
            })
        monkeypatch.setattr(service, "validate_study_allele_frequency", fail)
    with PipelineLogger(dataset_id, "dataset", str(tmp_path / "logs")) as logger:
        return service._prepare_dataset_for_chromosome_processing(
            sumstat_file=str(source), sample_column_dict=mapping,
            output_dir_path=tmp_path / "out", resource_folder=str(tmp_path),
            default_eaf_reference_source="fixture", default_comparison_af_file="fixture",
            dbsnp="fixture", resource_layout={}, output_layout=dict(module.output_layout.root),
            gwas2vcf_input=module.gwas2vcf_input.model_dump(), executables={},
            vcf_config=module.vcf_processing.model_dump(), user_eaf_file=None,
            user_eaf_column=None, user_info_file=None, user_info_column=None,
            default_eaf_reference_column="EUR", default_comparison_af_column="EUR",
            build_reference_files={}, build_check_colmap={}, default_eaf_colmap={},
            external_eaf_colmap={}, external_info_colmap={},
            external_reference_staging=module.external_reference_staging.model_dump(),
            policies=default_policies(), logger=logger, input_line_count=4,
        )


def test_resolved_report_and_manifest_evidence_follow_real_af_decision(tmp_path, monkeypatch, capsys):
    prepared = _prepare(tmp_path, monkeypatch)
    screen = capsys.readouterr().out
    assert screen.count("Study-wide decisions and evidence") == 1
    assert "Not required" in screen
    assert "Forward" in screen or "forward" in screen
    for key in ("effect_evidence", "pvalue_evidence", "beta_z_sign_consistency"):
        assert prepared.dataset_study_decisions[key] == prepared.study_decisions[key]
    assert prepared.dataset_study_decisions["effect_type_matches_declaration"] is True
    assert prepared.dataset_study_decisions["pvalue_type_matches_declaration"] is True
    assert prepared.dataset_study_decisions["se_scale_matches_declaration"] is None
    assert prepared.dataset_study_decisions["declaration_mismatches"] == []
    assert prepared.dataset_study_decisions["strand_consensus"]["minimum_informative_variants"] == 1
    partition = pl.read_parquet(prepared.chromosome_files["1"])
    assert partition["AF"].to_list() == [.1, .2, .7, .8]
    assert partition["BETA"].to_list() == [-.1, .1, .2, -.2]
    log = next((tmp_path / "logs").glob("*dataset.log")).read_text()
    assert log.index("Study-wide decisions and evidence") < log.index("Study AF accepted as declared EAF")
    assert log.index("Study AF accepted as declared EAF") < log.index("Declared study EAF retained")
    assert log.index("Selected build") < log.index("Prepare chromosome files")
    assert log.index("Imputation-quality score") < log.index("typed Parquet chromosome work")
    assert screen.count("Selected build") == 1
    assert "Chromosome 1 build evidence" not in screen
    assert "Imputation-quality score" in screen.split("Study-wide decisions and evidence", 1)[1]
    assert prepared.study_decisions["info_score_type"] == "standard_info"
    assert prepared.genome_build_info["per_chromosome"]["1"]["matches"]["GRCh37"] == 4


def test_af_failure_renders_known_evidence_before_propagating(tmp_path, monkeypatch, capsys):
    from postgwas.modules.harmonisation.allele_frequency import StudyAlleleFrequencyError
    with pytest.raises(StudyAlleleFrequencyError, match="reference fixture failure"):
        _prepare(tmp_path, monkeypatch, fail_frequency=True)
    screen = capsys.readouterr().out
    assert screen.count("Study-wide decisions and evidence") == 1
    assert "insufficient_overlap" in screen or "insufficient overlap" in screen.lower()
    assert "reference fixture failure" in screen
    assert not list((tmp_path / "out").glob("*chr*.parquet"))


def test_failed_manifest_keeps_existing_measured_study_evidence(tmp_path):
    from postgwas.modules.harmonisation.allele_frequency import StudyAlleleFrequencyError
    error = StudyAlleleFrequencyError("insufficient overlap", study_decisions={
        "effect_type": "beta", "eaf_reference_validation": {
            "status": "failed", "comparable_variants": 3,
        },
    })
    error.genome_build_info = {"inferred_build": "GRCh37", "matches": {"GRCh37": 4}}
    path = tmp_path / "manifest.json"
    assert service._finalise_failed_run_manifest(
        path, {"status": "RUNNING"}, error, time.time(),
        SimpleNamespace(error=lambda *_args: None),
    )
    manifest = json.loads(path.read_text())
    assert manifest["status"] == "FAILED"
    assert manifest["failure"]["message"] == "insufficient overlap"
    assert manifest["dataset"]["study_decisions"] == error.study_decisions
    assert manifest["dataset"]["genome_build"] == error.genome_build_info


@pytest.mark.parametrize("failure_stage", ["build", "resources", "replacement"])
def test_preflight_failure_renders_and_preserves_completed_build_evidence(
    tmp_path, monkeypatch, capsys, failure_stage,
):
    with pytest.raises(Exception) as caught:
        _prepare(tmp_path, monkeypatch, failure_stage=failure_stage)
    screen = capsys.readouterr().out
    assert screen.count("Study-wide decisions and evidence") == 1
    assert caught.value.genome_build_info["input_variants"] == 4
    assert caught.value.genome_build_info["inferred_build"] == (
        "Ambiguous" if failure_stage == "build" else "GRCh37"
    )
    assert "Dataset stopped" in screen
    assert not list((tmp_path / "out").glob("*chr*.parquet"))


def test_report_reuses_existing_file_numbers_without_source_reads(tmp_path, monkeypatch):
    from postgwas.core.input_validation import InputValidationSession, record_file_validation
    from postgwas.core.validation_reporting import FileValidationDisplay
    paths = [tmp_path / "one" / "reference.tsv", tmp_path / "two" / "reference.tsv"]
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        for path in paths:
            record_file_validation(path, "reference", checks=("schema",), status="passed")
        assert display.file_reference(str(paths[0])) == "File 1"
        assert display.file_reference(str(paths[1])) == "File 2"
        monkeypatch.setattr(Path, "open", lambda *_a, **_k: pytest.fail("report read a reference"))
        text = service._study_decisions_block({
            "eaf_column": "AF", "eaf_is_maf": True,
            "eaf_reference_validation": {"status": "passed", "comparable_variants": 2},
        }, genome_build_info={"inferred_build": "GRCh37"},
            build_reference_files={"GRCh37": str(paths[0])},
            resource_maps={"1": {"default_eaf_file": str(paths[1]), "default_eaf_reference_column": "EUR"}},
        )
    assert "File 1" in text and "File 2" in text
    assert str(paths[0]) not in text and str(paths[1]) not in text


@pytest.mark.parametrize("failure_stage,reason", [
    ("info", "INFO validation stopped"),
    ("staging", "staging fixture failure"),
    ("partition", "partition fixture failure"),
    ("chromosome_set", "different chromosome set"),
])
def test_late_preparation_failure_has_one_terminal_report_with_known_evidence(
        tmp_path, monkeypatch, capsys, failure_stage, reason):
    with pytest.raises(Exception, match=reason) as caught:
        _prepare(tmp_path, monkeypatch, failure_stage=failure_stage)
    screen = capsys.readouterr().out
    assert screen.count("Study-wide decisions and evidence") == 1
    assert screen.count("Selected build") == 1
    assert "Chromosome 1 build evidence" not in screen
    report = " ".join(screen.split()).split("Study-wide decisions and evidence", 1)[1]
    assert "Dataset stopped" in report and reason in report
    assert caught.value.study_decisions["effect_type"] == "beta"
    assert caught.value.study_decisions["pvalue_type"] == "raw"
    assert caught.value.genome_build_info["matches"]["GRCh37"] == 4
    assert "Final analysis completed" not in report
    if failure_stage != "info":
        assert caught.value.study_decisions["info_score_type"] == "standard_info"
    else:
        assert "info_score_type" not in caught.value.study_decisions


@pytest.mark.parametrize("failure", ["external_af", "replacement"])
def test_real_late_af_failure_does_not_replay_pending_study_report(tmp_path, monkeypatch, capsys, failure):
    if failure == "external_af":
        from test_harmonisation_external_af_preflight import (
            test_real_service_missing_native_af_enters_global_preparation_and_fails_before_workers,
        )
        test_real_service_missing_native_af_enters_global_preparation_and_fails_before_workers(
            tmp_path, monkeypatch, "inconclusive")
    else:
        from test_harmonisation_unified_frequency_policy import (
            test_service_activates_and_stages_external_maf_replacement_before_workers,
        )
        test_service_activates_and_stages_external_maf_replacement_before_workers(tmp_path, monkeypatch, True)
    screen = capsys.readouterr().out
    assert screen.count("Study-wide decisions and evidence") == 1
    assert screen.count("Selected build") == 1
    report = " ".join(screen.split()).split("Study-wide decisions and evidence", 1)[1]
    assert "Dataset stopped" in report
    assert "Effect type" in report and "P-value type" in report
    assert "accepted for assignment" not in report
    if failure == "replacement":
        assert "replacement frequency maf like" in report
        assert "Replacement distribution" in report
    else:
        assert "both rules failed" in report


def test_automatic_inference_notices_belong_to_single_orchestrated_report(tmp_path, monkeypatch, capsys):
    prepared = _prepare(tmp_path, monkeypatch, automatic_inference=True)
    screen = capsys.readouterr().out
    assert screen.count("Study-wide decisions and evidence") == 1
    assert "Effect type automatically inferred" not in screen
    assert "P-value type automatically inferred" not in screen
    report = " ".join(screen.split()).split("Study-wide decisions and evidence", 1)[1]
    assert "warning" in report.lower()
    assert "effect_type" in report and "pvalue.type" in report
    assert prepared.study_decisions["effect_type_source"] == "detector"
    assert prepared.study_decisions["pvalue_type_source"] == "detector"
    assert "AUTOMATIC EFFECT-TYPE INFERENCE" in prepared.study_decisions["effect_type_inference_warning"]
    assert "AUTOMATIC P-VALUE-TYPE INFERENCE" in prepared.study_decisions["pvalue_type_inference_warning"]


def test_direct_property_notices_remain_default_and_deferral_keeps_scientific_decisions(capsys):
    frame = pl.DataFrame({"BETA": [-.1, .1, .2, -.2], "P": [.2, .3, .4, .5]})
    columns = {"beta_or_col": "BETA", "pval_col": "P"}
    direct = resolve_study_properties(frame, columns)
    notice = capsys.readouterr().out
    assert notice.count("Effect type automatically inferred") == 1
    assert notice.count("P-value type automatically inferred") == 1
    deferred = resolve_study_properties(frame, columns, report_inference_warnings=False)
    assert capsys.readouterr().out == ""
    assert deferred == direct
    assert frame["BETA"].to_list() == [-.1, .1, .2, -.2]
    assert frame["P"].to_list() == [.2, .3, .4, .5]


def test_fatal_declaration_mismatch_is_immediate_even_when_routine_notices_deferred(capsys):
    frame = pl.DataFrame({"BETA": [-.1, .1, .2, -.2]})
    with pytest.raises(StudyPropertyError, match="DECLARATION MISMATCH") as caught:
        resolve_study_properties(frame, {"beta_or_col": "BETA", "declared_effect_type": "odds_ratio"},
            policies=default_policies().with_overrides({"validation.declaration_mismatch_action": "fail"}),
            report_inference_warnings=False)
    assert "DECLARATION MISMATCH" in capsys.readouterr().out
    assert caught.value.study_decisions["effect_type"] == "odds_ratio"
    assert caught.value.study_decisions["effect_type_matches_declaration"] is False
