"""Tests for dataset-level user-facing harmonisation output."""

import inspect
from io import StringIO

import polars as pl
from rich.console import Console

from postgwas.config import load_configuration
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.modules.harmonisation.cli import (
    _dataset_start_block,
    _engine_defaults,
    _preparation_block,
    _run_validated_rows,
)
from postgwas.modules.harmonisation.service import (
    DATASET_STEP_TOTAL,
    POST_MERGE_STEP_TOTAL,
    _HarmonisationProgress,
    _input_validation_summary_block,
    _ready_variant_type_counts,
    _study_decisions_block,
    _validated_column_mapping_block,
)


def test_run_preparation_context_precedes_dataset_processing():
    text = _preparation_block(
        sample_sheet="/data/studies.csv",
        dataset_count=7,
        resource_directory="/references",
        output_directory="/results",
    )

    assert text.splitlines() == [
        "    ▶️  Preparing PostGWAS harmonisation",
        "        🔹  Sample sheet         : /data/studies.csv",
        "        🧮  Datasets selected    : 7",
        "        🔹  Resource directory   : /references",
        "        🔹  Output directory     : /results",
        "        🔹  Current stage        : Validating configuration and required resources",
    ]
    source = inspect.getsource(_run_validated_rows)
    assert source.index("_preparation_block(") < source.index(
        "for dataset_index, row in enumerate(rows, start=1):"
    )


def test_direct_mode_preparation_does_not_invent_a_sample_sheet_path():
    text = _preparation_block(
        sample_sheet=None,
        dataset_count=1,
        resource_directory="/references",
        output_directory="/results",
    )

    assert "Sample sheet         : Direct API input — no sample-sheet path" in text


def test_each_dataset_start_shows_position_identity_and_input_file():
    text = _dataset_start_block(
        dataset_index=3,
        dataset_count=7,
        dataset_id="ADHD_female",
        input_file="/data/ADHD_female.sumstats.tsv.gz",
    )

    assert text.splitlines() == [
        "    ▶️  Starting dataset 3 of 7 — ADHD_female",
        "        🔹  Summary statistics   : /data/ADHD_female.sumstats.tsv.gz",
    ]
    source = inspect.getsource(_run_validated_rows)
    assert "for dataset_index, row in enumerate(rows, start=1):" in source
    assert source.index("_dataset_start_block(") < source.index("_write_run_metadata(")


def test_harmonisation_progress_uses_canonical_setting_and_validated_boundaries(
    tmp_path,
):
    configuration = load_configuration()
    engine = _engine_defaults(configuration)
    assert engine["terminal_progress"] == {
        "enabled": configuration.logging.show_progress,
        "outcome_label_width": configuration.logging.terminal_label_width,
    }

    stream = StringIO()
    progress = _HarmonisationProgress(
        enabled=True,
        outcome_label_width=configuration.logging.terminal_label_width,
        console=Console(
            file=stream,
            force_terminal=False,
            color_system=None,
            width=120,
        ),
    )
    logger = PipelineLogger(
        "study",
        "dataset",
        str(tmp_path),
        stage_progress=progress.dataset_stages,
    )

    for number in range(1, DATASET_STEP_TOTAL + 1):
        with logger.step(
            number,
            DATASET_STEP_TOTAL,
            "Dataset stage %d" % number,
            "test.dataset_stage",
        ):
            pass

    progress.start_chromosomes(("1", "2"))
    progress.record_chromosome_result("1", "ok", 1)
    progress.record_chromosome_result("2", "failed", 1)
    progress.record_chromosome_result("2", "ok", 2)
    progress.finish_chromosomes()

    logger.set_stage_progress(progress.start_post_merge())
    for number in range(1, POST_MERGE_STEP_TOTAL + 1):
        with logger.step(
            number,
            POST_MERGE_STEP_TOTAL,
            "Post-merge stage %d" % number,
            "test.post_merge_stage",
        ):
            pass
    logger.close()
    progress.close()

    text = stream.getvalue()
    assert "Harmonisation dataset stages" in text
    assert "Completed 8/8 · Dataset stage 8" in text
    assert "All 8 stages completed" in text
    assert "Harmonisation chromosome progress" in text
    assert "Progress 1/2 · Chromosome 1 completed on attempt 1 · 50%" in text
    assert "Current 1/2 · Chromosome 2 failed on attempt 1" in text
    assert (
        "Completed 2/2 · All chromosome outputs and required counts validated"
        in text
    )
    assert "Harmonisation post-merge stages" in text
    assert "Completed 5/5 · Post-merge stage 5" in text
    assert "All 5 stages completed" in text


def test_harmonisation_chromosome_progress_fails_below_completion():
    stream = StringIO()
    progress = _HarmonisationProgress(
        enabled=True,
        outcome_label_width=24,
        console=Console(
            file=stream,
            force_terminal=False,
            color_system=None,
            width=120,
        ),
    )

    progress.start_chromosomes(("1", "2", "X"))
    progress.record_chromosome_result("1", "ok", 1)
    progress.record_chromosome_result("2", "failed", 1)
    progress.finish_chromosomes(("2", "X"))
    progress.close()

    text = stream.getvalue()
    assert "Progress 1/3 · Chromosome 1 completed on attempt 1 · 33%" in text
    assert "Failed 1/3 · 2 chromosomes incomplete: 2, X" in text
    assert "Completed 3/3" not in text


def test_study_decisions_are_shown_as_aligned_user_facing_values():
    text = _study_decisions_block({
        "effect_type": "odds_ratio",
        "se_scale": "log_odds",
        "pvalue_type": "raw",
        "eaf_is_maf": False,
        "strand": "forward",
    })

    assert text.splitlines() == [
        "    🧠  Study-wide decisions",
        "        🔬  Effect type      : odds_ratio",
        "        🔬  SE scale         : log_odds",
        "        🔬  P-value type     : raw",
        "        🧬  Frequency type   : effect allele frequency",
        "        🧬  Strand consensus : forward",
    ]


def test_input_validation_reports_ready_snp_and_indel_other_counts():
    frame = pl.DataFrame({
        "EA": ["A", "A", "AC", "g"],
        "OA": ["G", "AT", "GT", "c"],
    })
    counts = _ready_variant_type_counts(
        frame, {"ea_col": "EA", "oa_col": "OA"},
    )

    assert counts == {"snps": 2, "indels_or_other_variants": 2}
    text = _input_validation_summary_block(
        input_variants=6,
        variants_read=6,
        missing_required=1,
        invalid_coordinates=0,
        unsupported_chromosomes=0,
        non_standard_alleles=0,
        duplicate_variants=1,
        ready_variants=4,
        ready_snps=counts["snps"],
        ready_indels_or_other=counts["indels_or_other_variants"],
    )

    assert "Ready for harmonisation        : 4" in text
    assert "Missing read-stage mandatory values" in " ".join(text.split())
    assert "Unsupported chromosomes" in text
    assert "Ready SNPs                     : 2" in text
    assert "Ready indels / other variants  : 2" in text


def test_maf_like_study_decision_is_displayed_as_provisional():
    text = _study_decisions_block({
        "effect_type": "beta",
        "pvalue_type": "raw",
        "eaf_is_maf": True,
        "strand": "forward",
    })

    assert "MAF-like; reference confirmation required" in text
    assert "SE scale         : not applicable to beta" in text


def test_validated_column_mapping_uses_public_sample_sheet_fields():
    columns = {
        "gwas_outputname": "ADHD2022", "sumstat_file": "/data/adhd.tsv.gz",
        "trait_type": "case_control", "delimiter": "whitespace",
        "chr_col": "CHR", "pos_col": "BP", "chr_pos_col": "NA",
        "snp_id_col": "SNP", "ea_col": "A1", "oa_col": "A2",
        "eaf_col": "FRQ_U", "beta_or_col": "OR", "se_col": "SE",
        "imp_z_col": "NA", "pval_col": "P", "ncontrol_col": "Nco",
        "ncase_col": "Nca", "imp_info_col": "INFO",
        "info_source_detail": "internal column INFO", "fixed_info": None,
        "declared_effect_type": "odds_ratio",
        "declared_pvalue_type": "raw",
        "ncontrol": "NA", "ncase": "NA",
        "infofile": "NA", "infocolumn": "NA",
        "eaffile": "NA", "eafcolumn": "NA",
    }
    decisions = {
        "effect_type": "odds_ratio", "pvalue_type": "raw", "eaf_is_maf": False,
    }

    text = _validated_column_mapping_block(columns, decisions)

    assert "    🔬  Validated sample-sheet values" in text
    assert "        🔹  Dataset" in text
    assert "        🧬  Variant columns" in text
    assert "        🔬  Association columns" in text
    assert "        🧮  Sample size" in text
    assert "        🔬  Quality and external sources" in text
    for field in (
        "dataset_id", "input_file", "trait_type", "chromosome_column",
        "position_column", "chromosome_position_column", "variant_id_column",
        "effect_allele_column", "other_allele_column",
        "effect_allele_frequency_column", "effect_column", "effect_type",
        "standard_error_column", "z_score_column", "p_value_column",
        "p_value_type", "control_count_column", "case_count_column",
        "control_count", "case_count", "imputation_info_column",
        "resolved_info_source", "fixed_info",
        "external_info_file", "external_info_column", "external_eaf_file",
        "external_eaf_column", "delimiter",
    ):
        assert field in text
    assert "dataset_id                     : ADHD2022" in text
    assert "effect_column                  : OR" in text
    assert "effect_type                    : odds_ratio" in text
    assert "z_score_column                 : Not provided — will not be used" in text
    assert "external_info_file             : Not provided — will not be used" in text
    assert "          🔹  dataset_id" in text
    assert "│" not in text
