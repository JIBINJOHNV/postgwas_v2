from pathlib import Path

import polars as pl
import pytest
from rich.cells import cell_len

from postgwas.modules.ld_clumping import ld_prune_standard


def _stub_conversion(monkeypatch, formatted):
    """Stand in for the VCF projection, which returns the table and the frame."""
    frame = pl.read_csv(
        formatted,
        separator="\t",
        schema_overrides={"chrcol": pl.String, "rsIDcol": pl.String},
    )
    monkeypatch.setattr(
        ld_prune_standard,
        "vcf_to_standard_ldclump",
        lambda *args, **kwargs: (str(formatted), frame),
    )


def _stub_reference_files(directory, chromosome):
    for name in (
        f"EUR_chr{chromosome}.ld.gz",
        f"EUR_chr{chromosome}.reverse.ld.gz",
        f"EUR_chr{chromosome}.variants.tsv.gz",
    ):
        path = directory / name
        path.write_bytes(b"reference")
        Path(str(path) + ".tbi").write_bytes(b"index")


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

    _stub_conversion(monkeypatch, formatted)
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
                "input_significant": 1,
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
    _stub_reference_files(tmp_path, "1")

    result = ld_prune_standard.ld_clump_standard(
        vcf_path=tmp_path / "study.vcf.gz",
        ld_folder=tmp_path,
        dataset_id="STUDY",
        output_directory=tmp_path,
        threads=1,
    )

    assert result["status"] == "no_loci"
    screen = " ".join(capsys.readouterr().out.split())
    assert "chr1" in screen
    assert "chr1 GWS variants" in screen
    assert "1 input · 1 retained · 0 excluded" in screen
    assert "Clumping results" in screen
    assert "independent 1 · lead 0 · loci 0" in screen
    assert "Significant but no risk locus" in screen
    assert "Other warnings" in screen
    assert "1 · see detailed log" in screen
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


def test_chromosome_progress_is_compact_aligned_and_visually_separated(capsys):
    ld_prune_standard._print_chromosome_progress(
        "1",
        {
            "status": "completed",
            "significant": 1825,
            "input_significant": 1835,
            "independent": 70,
            "lead": 32,
            "loci": 24,
        },
        10,
        {"missing_from_reference": 8, "allele_mismatch": 2},
    )

    rendered = capsys.readouterr().out
    fields = [line for line in rendered.splitlines() if " : " in line]
    scheduled = ld_prune_standard.screen_field(
        "count", "Chromosomes scheduled", 22, indent=6, label_width=24,
    )

    assert len(fields) == 3
    assert len({cell_len(line.split(" : ", 1)[0]) for line in fields}) == 1
    assert cell_len(fields[0].split(" : ", 1)[0]) == cell_len(
        scheduled.split(" : ", 1)[0]
    )
    assert "chr1 GWS variants" in rendered
    assert "1,835 input · 1,825 retained · 10 excluded" in rendered
    assert "Clumping results" in rendered
    assert "independent 70 · lead 32 · loci 24" in rendered
    assert "Excluded GWS variants" in rendered
    assert "not found 8 · allele mismatch 2" in rendered
    lines = rendered.splitlines()
    parent = next(line for line in lines if "chr1 GWS variants" in line)
    clumping = next(line for line in lines if "Clumping results" in line)
    exclusions = next(line for line in lines if "Excluded GWS variants" in line)
    assert len(clumping) - len(clumping.lstrip()) == (
        len(parent) - len(parent.lstrip()) + 2
    )
    assert len(exclusions) - len(exclusions.lstrip()) == (
        len(parent) - len(parent.lstrip()) + 2
    )
    assert "exact chromosome/position" not in rendered
    assert rendered.endswith("\n\n")


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

    _stub_conversion(monkeypatch, formatted)

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
    _stub_reference_files(tmp_path, "1")

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


def test_missing_significant_chromosome_reference_warns_and_skips_only_that_chromosome(
    monkeypatch, tmp_path, capsys,
):
    formatted = tmp_path / "STUDY_formatted.tsv"
    pl.DataFrame({
        "chrcol": ["1", "X", "2"],
        "poscol": [101, 201, 301],
        "neacol": ["A", "C", "G"],
        "eacol": ["G", "T", "A"],
        "rsIDcol": ["rs1", "rsX", "rs2"],
        "pcol": [1e-9, 2e-9, 0.5],
        "becol": [0.2, -0.3, 0.1],
        "secol": [0.1, 0.2, 0.1],
        "eafcol": [0.3, 0.4, 0.2],
    }).write_csv(formatted, separator="\t")
    configuration = (
        ld_prune_standard.load_configuration()
        .modules.ld_clumping.model_copy(deep=True)
    )
    configuration.reporting.result_table_columns = {
        **configuration.reporting.result_table_columns,
        "genomic_risk_loci": ["Genomic_locus", "CHR"],
        "lead_snps": ["Genomic_locus", "lead_SNP_id"],
        "independent_significant_snps": [
            "Genomic_locus", "ind_sig_SNP_id",
        ],
    }
    observed_chromosomes = []
    canonical_records = []

    class CaptureLogger:
        def record(self, marker, subject, **values):
            canonical_records.append((marker, subject, values))

    _stub_conversion(monkeypatch, formatted)

    def process_chromosome(chrom, *_args, **_kwargs):
        observed_chromosomes.append(str(chrom))
        significant = int(str(chrom) == "1")
        result = {
            "chrom": str(chrom),
            "logs": [],
            "progress": {
                "status": (
                    "completed"
                    if significant
                    else "no_significant_variants"
                ),
                "significant": significant,
                "independent": significant,
                "lead": significant,
                "loci": significant,
                "reference_ld_rows_before_maf": 10 * significant,
                "reference_ld_rows_retained": 8 * significant,
                "low_maf_partner_rows_excluded": 2 * significant,
                "low_maf_partner_variants_excluded": significant,
            },
        }
        if significant:
            result.update({
                "summ": pl.DataFrame({
                    "Genomic_locus": [1],
                    "CHR": ["1"],
                    "START": [90],
                    "END": [120],
                }),
                "hier": pl.DataFrame({"Genomic_locus": [1], "CHR": ["1"]}),
                "is_c": pl.DataFrame({
                    "Genomic_locus": [1],
                    "ind_sig_SNP_id": ["1_101_A_G"],
                }),
                "l_un": pl.DataFrame({
                    "Genomic_locus": [1],
                    "lead_SNP_id": ["1_101_A_G"],
                }),
            })
        return result

    monkeypatch.setattr(ld_prune_standard, "process_chromosome", process_chromosome)
    _stub_reference_files(tmp_path, "1")

    result = ld_prune_standard.ld_clump_standard(
        vcf_path=tmp_path / "study.vcf.gz",
        ld_folder=tmp_path,
        dataset_id="STUDY",
        output_directory=tmp_path,
        threads=1,
        memory_gb=2,
        include_report_tables=True,
        configuration=configuration,
        logger=CaptureLogger(),
    )

    assert result["status"] == "partial_reference"
    assert result["reference_coverage_status"] == "partial"
    assert result["input_significant_variants"] == 2
    assert result["significant_variants"] == 1
    assert result["skipped_significant_variants"] == 1
    assert result["skipped_chromosomes"] == ["X"]
    assert len(result["missing_reference_resources"]) == 6
    assert result["reference_ld_rows_before_maf"] == 10
    assert result["reference_ld_rows_retained"] == 8
    assert result["low_maf_partner_rows_excluded"] == 2
    assert result["low_maf_partner_variants_excluded"] == 1
    assert result[
        "reference_exclusions_within_reported_locus_boundaries"
    ] == 0
    assert result[
        "reference_exclusions_outside_reported_locus_boundaries"
    ] == 1
    assert list(result["_report_tables"]) == [
        "genomic_risk_loci",
        "lead_snps",
        "independent_significant_snps",
    ]
    assert result["_report_tables"]["genomic_risk_loci"]["columns"] == [
        "Genomic_locus", "CHR",
    ]
    assert result["_report_tables"]["genomic_risk_loci"]["rows"] == [
        (1, "1"),
    ]
    assert result["_reference_exclusion_details"][0][
        "reported_locus_boundary_status"
    ] == "outside_all_reported_locus_boundaries"
    assert Path(result["ldpruned_sig_file"]).is_file()
    assert sorted(observed_chromosomes) == ["1", "2"]
    screen = " ".join(capsys.readouterr().out.split())
    assert "chrX skipped" in screen
    assert "1 index variant · 6 LD resources missing/empty · no chrX loci" in screen
    assert "results are not genome-wide" not in screen
    detailed_log = Path(
        ld_prune_standard._standard_conversion_paths(
            tmp_path, "STUDY", configuration,
        )[1]
    ).read_text(encoding="utf-8")
    assert "significant_variants_skipped=1" in detailed_log
    assert "consequence=results are not genome-wide" in detailed_log
    assert "Excluded-index relationship to reported loci" in detailed_log
    assert "canonical_id=X_201_C_T" in detailed_log
    assert "reported_locus_boundary_status=outside_all_reported_locus_boundaries" in (
        detailed_log
    )
    audit = pl.read_csv(
        tmp_path / "STUDY_EUR_LD_Reference_Exclusions.tsv",
        separator="\t",
    )
    assert audit.select(
        "canonical_id",
        "reported_locus_boundary_status",
        "overlapping_genomic_loci",
    ).row(0) == (
        "X_201_C_T",
        "outside_all_reported_locus_boundaries",
        "",
    )
    assert any(
        marker == "WARNING"
        and subject == "standard_ld_reference_chromosome_skipped"
        and values["chromosome"] == "X"
        and values["significant_variants_skipped"] == 1
        for marker, subject, values in canonical_records
    )
    assert any(
        marker == "OUTPUT"
        and subject == "standard_reference_exclusions"
        and values[
            "reference_exclusions_outside_reported_locus_boundaries"
        ] == 1
        for marker, subject, values in canonical_records
    )
    assert any(
        marker == "OUTPUT"
        and subject == "standard_clumping"
        and values["low_maf_partner_rows_excluded"] == 2
        and "_report_tables" not in values
        and "_reference_exclusion_details" not in values
        for marker, subject, values in canonical_records
    )


def test_missing_significant_chromosome_reference_can_remain_strict(
    monkeypatch, tmp_path,
):
    formatted = tmp_path / "STUDY_formatted.tsv"
    pl.DataFrame({
        "chrcol": ["X"],
        "poscol": [201],
        "neacol": ["C"],
        "eacol": ["T"],
        "rsIDcol": ["rsX"],
        "pcol": [2e-9],
        "becol": [-0.3],
        "secol": [0.2],
        "eafcol": [0.4],
    }).write_csv(formatted, separator="\t")
    configuration = (
        ld_prune_standard.load_configuration()
        .modules.ld_clumping.model_copy(deep=True)
    )
    configuration.missing_chromosome_action = "error"
    _stub_conversion(monkeypatch, formatted)

    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match=r"Missing or empty LD resources.*EUR_chrX\.ld\.gz",
    ):
        ld_prune_standard.ld_clump_standard(
            vcf_path=tmp_path / "study.vcf.gz",
            ld_folder=tmp_path,
            dataset_id="STUDY",
            output_directory=tmp_path,
            threads=1,
            memory_gb=2,
            configuration=configuration,
        )


def test_all_significant_variants_skipped_fails_with_exclusion_audit(
    monkeypatch, tmp_path,
):
    formatted = tmp_path / "STUDY_formatted.tsv"
    pl.DataFrame({
        "chrcol": ["X"],
        "poscol": [201],
        "neacol": ["C"],
        "eacol": ["T"],
        "rsIDcol": ["rsX"],
        "pcol": [2e-9],
        "becol": [-0.3],
        "secol": [0.2],
        "eafcol": [0.4],
    }).write_csv(formatted, separator="\t")
    configuration = (
        ld_prune_standard.load_configuration()
        .modules.ld_clumping.model_copy(deep=True)
    )
    _stub_conversion(monkeypatch, formatted)

    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match="All significant variants were excluded",
    ):
        ld_prune_standard.ld_clump_standard(
            vcf_path=tmp_path / "study.vcf.gz",
            ld_folder=tmp_path,
            dataset_id="STUDY",
            output_directory=tmp_path,
            threads=1,
            memory_gb=2,
            configuration=configuration,
        )

    audit = pl.read_csv(
        tmp_path / "STUDY_EUR_LD_Reference_Exclusions.tsv",
        separator="\t",
    )
    assert audit.select("canonical_id", "reason", "action").rows() == [
        ("X_201_C_T", "missing_chromosome_reference", "warning_skip")
    ]
    assert audit.select(
        "reported_locus_boundary_status",
        "overlapping_genomic_loci",
        "locus_membership_interpretation",
    ).row(0) == (
        "outside_all_reported_locus_boundaries",
        "",
        "outside all reported locus boundaries",
    )


def test_non_autosomal_chromosomes_survive_the_projection(
    monkeypatch, tmp_path,
):
    formatted = tmp_path / "STUDY_formatted.tsv"
    row_count = 101
    pl.DataFrame({
        "chrcol": ["1"] * row_count + ["X"],
        "poscol": list(range(1, row_count + 2)),
        "neacol": ["A"] * (row_count + 1),
        "eacol": ["G"] * (row_count + 1),
        "rsIDcol": ["rs%d" % index for index in range(row_count)] + ["rsX"],
        "pcol": [0.5] * (row_count + 1),
        "becol": [0.2] * (row_count + 1),
        "secol": [0.1] * (row_count + 1),
        "eafcol": [0.3] * (row_count + 1),
    }).write_csv(formatted, separator="\t")
    configuration = (
        ld_prune_standard.load_configuration()
        .modules.ld_clumping.model_copy(deep=True)
    )
    observed_chromosomes = []

    _stub_conversion(monkeypatch, formatted)

    def process_chromosome(chrom, *_args, **_kwargs):
        observed_chromosomes.append(str(chrom))
        return {
            "chrom": str(chrom),
            "logs": [],
            "progress": {
                "status": "no_significant_variants",
                "significant": 0,
                "independent": 0,
                "lead": 0,
                "loci": 0,
            },
        }

    monkeypatch.setattr(ld_prune_standard, "process_chromosome", process_chromosome)

    result = ld_prune_standard.ld_clump_standard(
        vcf_path=tmp_path / "study.vcf.gz",
        ld_folder=tmp_path,
        dataset_id="STUDY",
        output_directory=tmp_path,
        threads=1,
        memory_gb=2,
        configuration=configuration,
    )

    assert result["status"] == "no_loci"
    assert sorted(observed_chromosomes) == ["1", "X"]
