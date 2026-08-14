from argparse import Namespace
from io import StringIO
import logging

from rich.console import Console

from postgwas.modules.fine_mapping.logging_utils import detailed_file_logging
from postgwas.modules.fine_mapping.presentation import FineMappingScreen
from postgwas.pipeline import runners


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
            "final_combined_credible_sets": tmp_path / "final.tsv",
            "flames_input": tmp_path / "flames",
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
        lambda _args: result,
    )
    assert runners.run_finemap_runner(args, context) == result

    terminal = capsys.readouterr().out
    assert terminal == ""
    assert "ld clumping completed" not in terminal
    assert "Finemap analysis started" not in terminal
    assert str(ld_result) not in terminal
    assert str(result) not in terminal
    assert str(tmp_path / "loci.tsv") not in terminal
