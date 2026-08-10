"""Scientific regression tests for study consensus and chromosome orientation."""

import math
from pathlib import Path

import polars as pl
import pytest

from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.modules.harmonisation.allele_frequency import harmonise_allele_frequency
from postgwas.modules.harmonisation.allele_frequency import confirm_eaf_with_reference
from postgwas.modules.harmonisation.effect_type import harmonise_effect_estimates
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.rejects import SOURCE_INPUT_ROW_COLUMN
from postgwas.modules.harmonisation.strand import (
    StrandOrientationError,
    harmonise_strand_orientation,
    resolve_strand_consensus,
)


MAPPING = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}
COLUMNS = {
    "chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA",
    "eaf_col": "EAF", "beta_or_col": "EFFECT", "imp_z_col": "Z",
}


def _reference(tmp_path):
    path = tmp_path / "GRCh37_1000G_freq_chr1.tsv"
    path.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.30\n"
        "1\t200\tC\tT\t0.28\n"
        "1\t300\tA\tG\t0.24\n"
        "1\t400\tT\tC\t0.20\n"
        "1\t500\tC\tG\t0.10\n"
        "1\t600\tA\tT\t0.48\n"
        "1\t700\tG\tT\t0.35\n",
        encoding="utf-8",
    )
    return str(path)


def _study():
    return pl.DataFrame({
        "CHR": ["1"] * 7,
        "POS": [100, 200, 300, 400, 500, 600, 700],
        "EA": ["A", "C", "C", "A", "C", "A", "A"],
        "OA": ["G", "T", "T", "G", "G", "T", "G"],
        "EAF": [0.31, 0.72, 0.25, 0.81, 0.11, 0.49, 0.40],
        "EFFECT": [0.20, 0.40, -0.15, 0.30, 0.12, 0.10, -0.20],
        "Z": [2.0, 4.0, -2.0, 3.0, 1.2, 1.0, -2.0],
    })


def test_raw_reference_single_join_resolves_four_actions_and_rejects_unsafe_rows(tmp_path):
    study = _study().with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    out, qc, _ = harmonise_strand_orientation(
        "1", study, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "reverse", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.get_column("POS").to_list() == [100, 200, 300, 400, 500]
    assert out.get_column("strand_action").to_list() == [
        "forward", "forward_swapped", "reverse_complement",
        "reverse_complement_swapped", "reverse_complement",
    ]
    assert out.get_column("EA").to_list() == ["A", "T", "G", "C", "G"]
    assert out.get_column("OA").to_list() == ["G", "C", "A", "T", "C"]
    assert out.get_column("EAF").to_list() == pytest.approx([0.31, 0.28, 0.25, 0.19, 0.11])
    assert out.get_column("EFFECT").to_list() == pytest.approx([0.20, -0.40, -0.15, -0.30, 0.12])
    assert out.get_column("Z").to_list() == pytest.approx([2.0, -4.0, -2.0, -3.0, 1.2])
    assert out.get_column(SOURCE_INPUT_ROW_COLUMN).to_list() == [1, 2, 3, 4, 5]
    assert qc["palindromic_ambiguous"] == 1
    assert qc["reference_unmatched"] == 1


def test_join_normalizes_integer_study_chromosome_to_reference_type(tmp_path):
    frame = _study().filter(pl.col("POS") == 100).with_columns(
        pl.col("CHR").cast(pl.Int64)
    )
    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "forward", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.height == 1
    assert out.item(0, "CHR") == "1"
    assert qc["actions"] == {"forward": 1}


def test_chromosome_log_identifies_exact_strand_reference_and_mapping(tmp_path):
    reference = _reference(tmp_path)
    logger = PipelineLogger(
        "study", "1", str(tmp_path), policies=default_policies(),
        level="INFO", screen_level="ERROR",
    )
    with logger.step(3, 16, "Allele frequency", "test", rows_in=1) as ctx:
        output, _, _ = harmonise_strand_orientation(
            "1", _study().filter(pl.col("POS") == 100), dict(COLUMNS),
            reference, "EUR", MAPPING,
            study_decision={
                "strand": "forward", "eaf_is_maf": False, "effect_type": "beta",
            },
            policies=default_policies(), ctx=ctx,
        )
        ctx.set_rows(output.height)
    logger.close()

    text = Path(logger.log_path).read_text(encoding="utf-8")
    assert "INPUT    Chromosome strand reference" in text
    assert "file=%s" % reference in text
    assert "population_frequency_column=EUR" in text
    assert "reference_allele_column=REF" in text
    assert "alternate_allele_column=ALT" in text
    assert "OBSERVED Applied study strand consensus value=forward" in text


def test_odds_ratio_is_reciprocated_but_maf_is_not_inverted(tmp_path):
    frame = _study().filter(pl.col("POS") == 200).with_columns([
        pl.lit(2.0).alias("EFFECT"),
        pl.lit(0.20).alias("EAF"),
    ])
    out, _, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "forward", "eaf_is_maf": True, "effect_type": "odds_ratio"},
        policies=default_policies(),
    )

    assert out.item(0, "EFFECT") == pytest.approx(0.5)
    assert out.item(0, "EAF") == pytest.approx(0.20)
    assert out.item(0, "Z") == pytest.approx(-4.0)


def test_raw_or_scale_se_is_converted_before_swapped_allele_orientation(tmp_path):
    frame = _study().filter(pl.col("POS") == 200).with_columns([
        pl.lit(2.0).alias("EFFECT"),
        pl.lit(0.30).alias("SE"),
    ])
    columns = dict(COLUMNS, se_col="SE")
    policies = default_policies()

    frame, effect_qc, columns = harmonise_effect_estimates(
        "1",
        frame,
        columns,
        decision={"effect_type": "odds_ratio", "se_scale": "as_given"},
        policies=policies,
    )
    out, _, _ = harmonise_strand_orientation(
        "1",
        frame,
        columns,
        _reference(tmp_path),
        "EUR",
        MAPPING,
        study_decision={
            "strand": "forward",
            "eaf_is_maf": False,
            "effect_type": "beta",
            "input_effect_type": "odds_ratio",
            "se_scale": "as_given",
        },
        policies=policies,
    )

    assert effect_qc["se_rescaled_by_or"] is True
    assert out.item(0, "beta") == pytest.approx(-math.log(2.0))
    assert out.item(0, "SE") == pytest.approx(0.15)


def test_strand_rejects_unstandardized_raw_or_scale_se(tmp_path):
    frame = _study().filter(pl.col("POS") == 200).with_columns(
        pl.lit(0.30).alias("SE")
    )
    with pytest.raises(StrandOrientationError, match="harmonise_effect_estimates"):
        harmonise_strand_orientation(
            "1",
            frame,
            dict(COLUMNS, se_col="SE"),
            _reference(tmp_path),
            "EUR",
            MAPPING,
            study_decision={
                "strand": "forward",
                "eaf_is_maf": False,
                "effect_type": "odds_ratio",
                "se_scale": "as_given",
            },
            policies=default_policies(),
        )


def test_consensus_uses_non_palindromic_orientation_evidence():
    policies = default_policies().with_overrides({
        "strand": {"min_informative_variants": 10, "consensus_threshold": 0.99},
    })
    decision = resolve_strand_consensus({
        "strand_evidence": {"GRCh37": {"forward": 1, "reverse": 999, "ambiguous": 4}},
    }, "GRCh37", policies)

    assert decision["strand"] == "reverse"
    assert decision["informative"] == 1000
    assert decision["dominant_fraction"] == pytest.approx(0.999)


def test_fail_policy_stops_on_reference_unmatched(tmp_path):
    policies = default_policies().with_overrides({
        "strand": {"unmatched_action": "fail"},
    })
    with pytest.raises(StrandOrientationError, match="reference-unmatched"):
        harmonise_strand_orientation(
            "1", _study().filter(pl.col("POS") == 700), dict(COLUMNS),
            _reference(tmp_path), "EUR", MAPPING,
            study_decision={"strand": "forward", "eaf_is_maf": False, "effect_type": "beta"},
            policies=policies,
        )


def test_absent_coordinate_is_rejected_as_reference_unmatched(tmp_path):
    frame = _study().filter(pl.col("POS") == 100).with_columns(
        pl.lit(999).alias("POS")
    )
    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "forward", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.is_empty()
    assert qc["reference_unmatched"] == 1


def test_disabled_strand_policy_skips_consensus_and_chromosome_orientation(tmp_path):
    policies = default_policies().with_overrides({"strand": {"enabled": False}})
    decision = resolve_strand_consensus({}, "GRCh37", policies)
    out, qc, _ = harmonise_strand_orientation(
        "1", _study(), dict(COLUMNS), "does-not-exist.tsv", "EUR", MAPPING,
        policies=policies,
    )

    assert decision["strand"] == "disabled"
    assert out.equals(_study())
    assert qc["status"] == "disabled"


def test_maf_confirmation_reuses_aligned_reference_without_reading_file():
    frame = pl.DataFrame({
        "CHR": ["1", "1"], "POS": [1, 2], "EA": ["A", "G"], "OA": ["G", "A"],
        "MAF": [0.20, 0.10], "strand_reference_af": [0.80, 0.90],
    })
    policies = default_policies().with_overrides({
        "eaf": {"maf_reference_min_overlap": 2},
    })
    decision, stats = confirm_eaf_with_reference(
        frame, "MAF", "does-not-exist.tsv", "EUR", MAPPING,
        {"chr": "CHR", "pos": "POS", "ea": "EA", "oa": "OA"},
        policies=policies,
        aligned_reference_col="strand_reference_af",
    )

    assert decision == "maf"
    assert stats["comparable_variants"] == 2


def test_final_eaf_is_compared_with_aligned_reference_af(tmp_path):
    frame = _study().filter(pl.col("POS") == 100).with_columns(
        pl.lit(0.90).alias("EAF")
    )
    columns = dict(COLUMNS, gwas_outputname="study", output_folder=str(tmp_path))
    layout = {
        "frequency_qc_directory": "eaf_qc",
        "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
        "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
    }

    out, qc, _ = harmonise_allele_frequency(
        "1", frame, columns, layout, "\t",
        external_eaf_colmap=MAPPING,
        default_eaf_file=_reference(tmp_path),
        default_eaf_column="EUR",
        default_eaf_colmap=MAPPING,
        study_decision={"strand": "forward", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.height == 1
    assert out.item(0, "strand_af_difference") == pytest.approx(0.60)
    assert qc["strand_af_comparable"] == 1
    assert qc["strand_af_discordant"] == 1


def test_maf_like_column_proven_to_be_eaf_is_then_inverted_for_swapped_row(tmp_path):
    frame = _study().filter(pl.col("POS") == 200).with_columns(
        pl.lit(0.72).alias("EAF")
    )
    # Reference ALT frequency is 0.28, so the original REF effect allele has
    # frequency 0.72. After the proof that this is EAF, standardized ALT EAF is 0.28.
    columns = dict(COLUMNS, gwas_outputname="study", output_folder=str(tmp_path))
    layout = {
        "frequency_qc_directory": "eaf_qc",
        "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
        "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
    }
    policies = default_policies().with_overrides({
        "eaf": {"maf_reference_min_overlap": 1},
    })

    out, qc, _ = harmonise_allele_frequency(
        "1", frame, columns, layout, "\t",
        external_eaf_colmap=MAPPING,
        default_eaf_file=_reference(tmp_path),
        default_eaf_column="EUR",
        default_eaf_colmap=MAPPING,
        study_decision={"strand": "forward", "eaf_is_maf": True, "effect_type": "beta"},
        policies=policies,
    )

    assert qc["maf_reference_decision"] == "eaf"
    assert out.item(0, "EAF") == pytest.approx(0.28)
    assert out.item(0, "strand_af_difference") == pytest.approx(0.0)
