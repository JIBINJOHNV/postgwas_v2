from pathlib import Path

import polars as pl

from postgwas.modules.ld_clumping import ld_prune_standard


def test_standard_clumping_prints_summary_and_logs_variant_details(
    monkeypatch, tmp_path, capsys
):
    formatted = tmp_path / "STUDY_formatted.tsv"
    pl.DataFrame(
        {
            "chrcol": ["1"],
            "poscol": [101],
            "neacol": ["A"],
            "eacol": ["G"],
            "rsIDcol": ["rs1"],
            "pcol": [1e-9],
            "becol": [0.2],
            "secol": [0.1],
            "eafcol": [0.3],
        }
    ).write_csv(formatted, separator="\t")

    monkeypatch.setattr(
        ld_prune_standard,
        "vcf_to_standard_ldclump",
        lambda *args, **kwargs: str(formatted),
    )
    monkeypatch.setattr(
        ld_prune_standard,
        "process_chromosome",
        lambda *args, **kwargs: {
            "chrom": "1",
            "logs": [
                "[INFO] selected_variant=1_101_A_G",
                "[WARNING] action=self-only | expected_id=1_101_A_G",
            ],
            "progress": {
                "status": "no_lead_variants",
                "significant": 1,
                "independent": 1,
                "lead": 0,
                "loci": 0,
            },
        },
    )
    monkeypatch.setattr(
        "postgwas.core.execution.runtime.auto_detect_workers", lambda: 1
    )
    monkeypatch.setattr(
        "postgwas.core.execution.runtime.safe_thread_count", lambda *args, **kwargs: 1
    )
    (tmp_path / "EUR_chr1.ld.gz").write_bytes(b"ld")
    (tmp_path / "EUR_chr1.ld.gz.tbi").write_bytes(b"index")

    result = ld_prune_standard.ld_clump_standard(
        vcf_path=tmp_path / "study.vcf.gz",
        ld_folder=tmp_path,
        dataset_id="STUDY",
        output_directory=tmp_path,
        threads=1,
    )

    assert result["status"] == "no_loci"
    screen = capsys.readouterr().out
    assert "chr1" in screen
    assert "1 significant · 1 independent · 0 lead · 0 risk loci" in screen
    assert "Significant but no risk locus" in screen
    assert "1 warning" in screen
    assert "Detailed log" in screen
    assert str(tmp_path) not in screen
    assert "selected_variant=1_101_A_G" not in screen
    assert "===== CHROMOSOME" not in screen
    detailed_log = Path(
        ld_prune_standard._standard_conversion_paths(tmp_path, "STUDY")[1]
    )
    log_text = detailed_log.read_text(encoding="utf-8")
    assert "Chromosome 1 clumping diagnostics" in log_text
    assert "selected_variant=1_101_A_G" in log_text
    assert "action=self-only" in log_text


def test_standard_clumping_groups_chromosomes_without_significant_variants(
    monkeypatch, tmp_path, capsys
):
    formatted = tmp_path / "STUDY_formatted.tsv"
    pl.DataFrame(
        {
            "chrcol": ["1", "2", "3", "5"],
            "poscol": [101, 201, 301, 501],
            "neacol": ["A"] * 4,
            "eacol": ["G"] * 4,
            "rsIDcol": ["rs1", "rs2", "rs3", "rs5"],
            "pcol": [1e-9, 0.5, 0.6, 0.7],
            "becol": [0.2] * 4,
            "secol": [0.1] * 4,
            "eafcol": [0.3] * 4,
        }
    ).write_csv(formatted, separator="\t")

    monkeypatch.setattr(
        ld_prune_standard,
        "vcf_to_standard_ldclump",
        lambda *args, **kwargs: str(formatted),
    )

    def process_chromosome(chrom, *args, **kwargs):
        has_signal = str(chrom) == "1"
        return {
            "chrom": str(chrom),
            "logs": [f"[INFO] chromosome={chrom} diagnostics"],
            "progress": {
                "status": (
                    "no_lead_variants" if has_signal else "no_significant_variants"
                ),
                "significant": int(has_signal),
                "independent": int(has_signal),
                "lead": 0,
                "loci": 0,
            },
        }

    monkeypatch.setattr(
        ld_prune_standard,
        "process_chromosome",
        process_chromosome,
    )
    monkeypatch.setattr(
        "postgwas.core.execution.runtime.auto_detect_workers", lambda: 1
    )
    monkeypatch.setattr(
        "postgwas.core.execution.runtime.safe_thread_count", lambda *args, **kwargs: 1
    )
    (tmp_path / "EUR_chr1.ld.gz").write_bytes(b"ld")
    (tmp_path / "EUR_chr1.ld.gz.tbi").write_bytes(b"index")

    result = ld_prune_standard.ld_clump_standard(
        vcf_path=tmp_path / "study.vcf.gz",
        ld_folder=tmp_path,
        dataset_id="STUDY",
        output_directory=tmp_path,
        threads=1,
    )

    assert result["status"] == "no_loci"
    screen = capsys.readouterr().out
    assert screen.count("chr2") == 1
    assert screen.count("chr3") == 0
    assert screen.count("chr5") == 1
    assert "No significant variants" in screen
    assert "chr2–3, chr5 (3)" in screen
    assert "[2/4]" not in screen
    assert "[3/4]" not in screen
    assert "[4/4]" not in screen
    detailed_log = Path(
        ld_prune_standard._standard_conversion_paths(tmp_path, "STUDY")[1]
    )
    log_text = detailed_log.read_text(encoding="utf-8")
    assert "Chromosome 2 clumping diagnostics" in log_text
    assert "Chromosome 3 clumping diagnostics" in log_text
    assert "Chromosome 5 clumping diagnostics" in log_text


def test_compact_chromosome_labels_orders_ranges_and_sex_chromosomes():
    assert ld_prune_standard._compact_chromosome_labels(
        ["chr22", "2", "3", "1", "X", "5", "CHR2"]
    ) == "chr1–3, chr5, chr22, chrX"


def test_chromosome_failure_is_preserved_in_detailed_log(tmp_path):
    detailed_log = tmp_path / "clumping.log"
    error = ld_prune_standard.PipelineStageError(
        "04 LD reference validation",
        "process_chromosome",
        "LD reference file was not found",
        chromosome="2",
    )

    ld_prune_standard._append_chromosome_diagnostics(
        detailed_log,
        "2",
        ["[INFO] chromosome processing started"],
        error=error,
    )

    log_text = detailed_log.read_text(encoding="utf-8")
    assert "Chromosome 2 clumping diagnostics" in log_text
    assert "chromosome processing started" in log_text
    assert "[ERROR] [STAGE: 04 LD reference validation]" in log_text
    assert "LD reference file was not found" in log_text
