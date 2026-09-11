from argparse import Namespace
from io import StringIO
import logging

import pytest
from rich.cells import cell_len
from rich.console import Console

from postgwas.modules.fine_mapping.logging_utils import detailed_file_logging
from postgwas.modules.fine_mapping.presentation import FineMappingScreen
from postgwas.pipeline import runners


@pytest.mark.parametrize("width", [80, 120])
def test_final_summary_paths_keep_value_alignment_at_console_width(tmp_path, width):
    stream = StringIO()
    screen = FineMappingScreen(
        "FINEMAP", show_progress=False, label_width=42,
        console=Console(file=stream, force_terminal=False, color_system=None, width=width),
    )
    root = tmp_path / "long_parent_directory" / "study_with_a_long_reproducible_name"
    paths = {
        "Scientific HTML report": ("html_report", "results/study_fine_mapping_report.html"),
        "Final combined result": ("final_combined_credible_sets", "results/combined_results/final_combined_credible_sets.tsv"),
        "FLAMES input": ("flames_input", "downstream_inputs/flames"),
        "Input and resource QC": ("preflight_validation", "quality_control/input_and_resource_validation.tsv"),
        "Locus QC": ("locus_status", "quality_control/finemap_locus_status.tsv"),
        "QC summary": ("overlap_resolution_summary", "quality_control/overlap_resolution_summary.tsv"),
    }
    screen.print_final_summary(
        {"status": "success", "output_dir": root, "n_attempted": 8, "n_successful": 8,
         "n_final_credible_sets": 14, "n_warnings": 8,
         "warning_reason_counts": {"nef_relative_range_exceeds_threshold;ld_z_mismatch_warning": 8},
         **{key: root / relative for key, relative in paths.values()}},
        root / "run_metadata/pipeline_summary.log",
    )
    lines = stream.getvalue().splitlines()
    assert all(cell_len(line) <= width for line in lines)
    expected_paths = {label: relative for label, (_, relative) in paths.items()}
    expected_paths.update({"Output directory": str(root), "Full detailed log": "run_metadata/pipeline_summary.log"})
    for label, expected in expected_paths.items():
        index = next(index for index, line in enumerate(lines) if label in line)
        prefix, first = lines[index].split(" : ", 1)
        assert first
        indentation = " " * cell_len(prefix + " : ")
        fragments = [first]
        for continuation in lines[index + 1:]:
            if not continuation.startswith(indentation):
                break
            fragments.append(continuation[len(indentation):])
        assert "".join(fragments) == expected


def test_fine_mapping_screen_is_concise_and_decision_focused(tmp_path):
    stream = StringIO()
    screen = FineMappingScreen(
        "SuSiE-RSS",
        show_progress=True,
        label_width=42,
        console=Console(
            file=stream,
            force_terminal=True,
            color_system=None,
            width=120,
        ),
    )

    screen.start()
    for stage in range(1, 9):
        screen.complete(stage, [("count", "Completed items", stage)])
    screen.print_final_summary(
        {
            "status": "success",
            "n_attempted": 4,
            "n_successful": 3,
            "n_failed": 1,
            "n_warnings": 4,
            "warning_reason_counts": {
                "nef_relative_range_exceeds_threshold;ld_z_mismatch_warning": 4
            },
            "failure_reason_counts": {"finemap_timeout": 1},
            "n_overlap_groups": 0,
            "n_final_credible_sets": 4,
            "output_dir": tmp_path,
            "html_report": tmp_path / "report.html",
            "final_combined_credible_sets": tmp_path / "final.tsv",
            "flames_input": tmp_path / "flames",
            "preflight_validation": tmp_path / "preflight.tsv",
            "locus_status": tmp_path / "locus_qc.tsv",
            "overlap_resolution_summary": tmp_path / "overlap_qc.tsv",
        },
        tmp_path / "pipeline_summary.log",
    )

    text = stream.getvalue()
    normalized_text = " ".join(text.split())
    for stage in range(1, 9):
        assert f"Completed {stage}/8" in text
    assert "All 8 stages completed" in text
    assert "Fine-mapping completed" in text
    assert "Primary loci successful" in text
    assert "Final credible sets" in text
    assert (
        "nef relative range exceeds threshold; ld z mismatch warning"
        in normalized_text
    )
    assert "finemap timeout (1 locus)" in normalized_text
    assert "final.tsv" in text
    assert "pipeline_summary.log" in text
    assert "preflight.tsv" in text
    assert "report.html" in text
    assert "Full detailed log" in text
    assert "[STAGE]" not in text
    assert "[PROGRESS]" not in text


def test_detailed_fine_mapping_events_are_file_only(tmp_path, capsys):
    log_file = tmp_path / "pipeline_summary.log"
    logger = logging.getLogger("postgwas.modules.fine_mapping.test")

    with detailed_file_logging(
        "postgwas.modules.fine_mapping",
        log_file,
        "INFO",
    ):
        logger.info(
            "[STAGE] stage=worker_metadata status=copied source=/intermediate"
        )
        logger.info(
            "[PROGRESS] scope=pipeline stage=complete status=success"
        )

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""
    audit = log_file.read_text(encoding="utf-8")
    assert "stage=worker_metadata" in audit
    assert "scope=pipeline stage=complete" in audit


def test_pipeline_runners_do_not_print_raw_fine_mapping_objects(
    monkeypatch, tmp_path, capsys
):
    ld_result = {
        "ld_clump_standard": {
            "ldpruned_sig_file": str(tmp_path / "loci.tsv")
        }
    }
    formatter = {
        "susie": {"susie_input": str(tmp_path / "susie.tsv.gz")}
    }
    args = Namespace(
        output_directory=str(tmp_path),
        _step_num=4,
        finemap_method="susie",
    )
    context = {"formatter": formatter}

    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.run_ld_clump_direct",
        lambda _args, **_kwargs: ld_result,
    )
    assert runners.run_ld_clump_runner(args, context) == ld_result

    result = {
        "status": "success",
        "output_dir": str(tmp_path / "04_finemap"),
    }
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.service.run_fine_mapping",
        lambda _args, **_kwargs: result,
    )
    assert runners.run_finemap_runner(args, context) == result

    terminal = capsys.readouterr().out
    assert terminal == ""
    assert "ld clumping completed" not in terminal
    assert "Finemap analysis started" not in terminal
    assert str(ld_result) not in terminal
    assert str(result) not in terminal
    assert str(tmp_path / "loci.tsv") not in terminal


def test_ld_clumping_runner_passes_current_vcf_evidence_to_the_service(
    monkeypatch, tmp_path,
):
    from postgwas.core.contracts import RunContext

    indexed = object()
    captured = {}
    args = Namespace(
        output_directory=str(tmp_path),
        vcf="study.vcf.gz",
        _step_num=1,
    )
    context = RunContext()
    monkeypatch.setattr(
        runners,
        "_validate_current_pipeline_vcf",
        lambda *_args: {"indexed": indexed},
    )

    def run(_args, **kwargs):
        captured.update(kwargs)
        return {"ld_clump_standard": {"ldpruned_sig_file": "loci.tsv"}}

    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.run_ld_clump_direct",
        run,
    )

    runners.run_ld_clump_runner(args, context)

    assert captured["cached_vcf"] is indexed
    assert captured["pipeline"] is True


@pytest.mark.parametrize("engine,input_key", [("susie", "susie_input"), ("finemap", "finemap_input")])
def test_finemap_runner_forwards_prior_reports_for_both_engines(tmp_path, monkeypatch, engine, input_key):
    from postgwas.core.contracts import RunContext

    context = RunContext({
        "ld_clump": {"ld_clump_standard": {"ldpruned_sig_file": "loci.tsv"},
                     "html_report": "clumping.html"},
        "formatter": {engine: {input_key: "summary.tsv", "html_report": "formatting.html"}},
        "finemap": {"html_report": "previous_finemap.html"},
    })
    captured = {}

    def run(_args, **kwargs):
        captured.update(kwargs)
        return {"status": "success"}

    monkeypatch.setattr("postgwas.modules.fine_mapping.service.run_fine_mapping", run)
    args = Namespace(output_directory=str(tmp_path), finemap_method=engine, _step_num=3)
    runners.run_finemap_runner(args, context)
    assert set(captured["upstream_results"]) == {"ld_clump", "formatter"}
    assert captured["upstream_results"]["ld_clump"]["html_report"] == "clumping.html"
    assert args.output_directory == str(tmp_path)
