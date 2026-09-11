"""Scientific regression tests for study consensus and chromosome orientation."""

import math
from pathlib import Path

import polars as pl
import pytest

from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.modules.harmonisation import strand as strand_module
from postgwas.modules.harmonisation.allele_frequency import (
    AlleleFrequencyError,
    harmonise_allele_frequency,
)
from postgwas.modules.harmonisation.allele_frequency import confirm_eaf_with_reference
from postgwas.modules.harmonisation.effect_type import harmonise_effect_estimates
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    SOURCE_INPUT_ROW_COLUMN,
)
from postgwas.modules.harmonisation.strand import (
    REFERENCE_AF_COLUMN,
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


def _step_four_layout():
    return {
        "frequency_qc_directory": "eaf_qc",
        "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
        "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
        "post_orientation_duplicates": (
            "qc/{dataset_id}_duplicates_chr{chromosome}.tsv"
        ),
    }


def _reject_collector(frame, tmp_path):
    return RejectCollector(
        source_snapshot=frame,
        logger=None,
        out_path=str(tmp_path / "rejected.tsv"),
        delimiter="\t",
        compress=False,
    )


def _frequency_alignment_reference(tmp_path, first=0.10, second=0.20):
    path = tmp_path / "frequency_alignment_reference.tsv"
    path.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t%s\n"
        "1\t200\tC\tT\t%s\n"
        "1\t300\tA\tT\t0.10\n" % (first, second),
        encoding="utf-8",
    )
    return str(path)


def _frequency_alignment_study(first, second, include_palindrome=False):
    frame = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [100, 200],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
        "EAF": [first, second],
        "EFFECT": [0.20, -0.10],
        "Z": [2.0, -1.0],
    })
    if include_palindrome:
        frame = pl.concat([
            frame,
            pl.DataFrame({
                "CHR": ["1"], "POS": [300], "EA": ["A"], "OA": ["T"],
                "EAF": [0.90], "EFFECT": [0.30], "Z": [3.0],
            }),
        ])
    return frame


def _frequency_alignment_policies():
    return default_policies().with_overrides({
        "eaf.non_effect_frequency_min_overlap": 2,
    })


def test_correct_effect_allele_frequency_passes_pre_orientation_check(tmp_path):
    out, qc, _ = harmonise_strand_orientation(
        "1",
        _frequency_alignment_study(0.10, 0.80),
        dict(COLUMNS),
        _frequency_alignment_reference(tmp_path),
        "EUR",
        MAPPING,
        study_decision={
            "strand": "forward", "eaf_is_maf": False, "effect_type": "beta",
        },
        policies=_frequency_alignment_policies(),
    )

    evidence = qc["study_frequency_allele_alignment"]
    assert evidence["decision"] == "effect_allele_frequency"
    assert evidence["comparable_variants"] == 2
    assert evidence["declared_eaf_mean_absolute_error"] == pytest.approx(0.0)
    assert evidence["inverted_eaf_mean_absolute_error"] == pytest.approx(0.7)
    assert out.get_column("EAF").to_list() == pytest.approx([0.10, 0.20])


def test_non_effect_allele_frequency_fails_before_any_frequency_transform(tmp_path):
    with pytest.raises(
        StrandOrientationError,
        match=(
            "appears to contain the frequency of the non-effect allele.*"
            "PostGWAS did not modify the data"
        ),
    ):
        harmonise_strand_orientation(
            "1",
            _frequency_alignment_study(0.90, 0.20, include_palindrome=True),
            dict(COLUMNS),
            _frequency_alignment_reference(tmp_path),
            "EUR",
            MAPPING,
            study_decision={
                "strand": "mixed", "eaf_is_maf": False, "effect_type": "beta",
            },
            policies=_frequency_alignment_policies(),
        )


def test_near_half_frequency_evidence_is_inconclusive_and_not_changed(tmp_path):
    frame = _frequency_alignment_study(0.48, 0.52).with_columns([
        pl.Series("EA", ["A", "T"]),
        pl.Series("OA", ["G", "C"]),
    ])
    out, qc, _ = harmonise_strand_orientation(
        "1",
        frame,
        dict(COLUMNS),
        _frequency_alignment_reference(tmp_path, first=0.48, second=0.52),
        "EUR",
        MAPPING,
        study_decision={
            "strand": "forward", "eaf_is_maf": False, "effect_type": "beta",
        },
        policies=_frequency_alignment_policies(),
    )

    evidence = qc["study_frequency_allele_alignment"]
    assert evidence["decision"] == "inconclusive"
    assert evidence["decision_reason"] == (
        "evidence_did_not_support_one_interpretation"
    )
    assert out.get_column("EAF").to_list() == pytest.approx([0.48, 0.52])


def test_population_frequency_difference_does_not_false_diagnose_inversion(
    tmp_path,
):
    out, qc, _ = harmonise_strand_orientation(
        "1",
        _frequency_alignment_study(0.30, 0.65),
        dict(COLUMNS),
        _frequency_alignment_reference(tmp_path),
        "EUR",
        MAPPING,
        study_decision={
            "strand": "forward", "eaf_is_maf": False, "effect_type": "beta",
        },
        policies=_frequency_alignment_policies(),
    )

    evidence = qc["study_frequency_allele_alignment"]
    assert evidence["decision"] == "inconclusive"
    assert evidence["declared_eaf_correlation"] == pytest.approx(1.0)
    assert evidence["declared_eaf_mean_absolute_error"] == pytest.approx(0.175)
    assert out.get_column("EAF").to_list() == pytest.approx([0.30, 0.35])


def test_frequency_inversion_diagnosis_respects_yaml_error_margin(tmp_path):
    policies = _frequency_alignment_policies().with_overrides({
        "eaf.non_effect_frequency_error_margin": 0.80,
    })
    out, qc, _ = harmonise_strand_orientation(
        "1",
        _frequency_alignment_study(0.90, 0.20),
        dict(COLUMNS),
        _frequency_alignment_reference(tmp_path),
        "EUR",
        MAPPING,
        study_decision={
            "strand": "forward", "eaf_is_maf": False, "effect_type": "beta",
        },
        policies=policies,
    )

    evidence = qc["study_frequency_allele_alignment"]
    assert evidence["decision"] == "inconclusive"
    assert evidence["minimum_error_margin"] == pytest.approx(0.80)
    assert out.get_column("EAF").to_list() == pytest.approx([0.90, 0.80])


def test_raw_reference_single_join_resolves_four_actions_and_rejects_unsafe_rows(tmp_path):
    study = _study().with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    out, qc, _ = harmonise_strand_orientation(
        "1", study, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "reverse", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.get_column("POS").to_list() == [100, 200, 300, 400, 500, 600]
    assert out.get_column("strand_action").to_list() == [
        "forward", "forward_swapped", "reverse_complement",
        "reverse_complement_swapped", "reverse_complement",
        "reverse_complement",
    ]
    assert out.get_column("EA").to_list() == ["A", "T", "G", "C", "G", "T"]
    assert out.get_column("OA").to_list() == ["G", "C", "A", "T", "C", "A"]
    assert out.get_column("EAF").to_list() == pytest.approx(
        [0.31, 0.28, 0.25, 0.19, 0.11, 0.49]
    )
    assert out.get_column("EFFECT").to_list() == pytest.approx(
        [0.20, -0.40, -0.15, -0.30, 0.12, 0.10]
    )
    assert out.get_column("Z").to_list() == pytest.approx(
        [2.0, -4.0, -2.0, -3.0, 1.2, 1.0]
    )
    assert out.get_column(SOURCE_INPUT_ROW_COLUMN).to_list() == [1, 2, 3, 4, 5, 6]
    assert qc["palindromic_ambiguous"] == 0
    assert qc["reference_unmatched"] == 1
    assert qc["palindromic_orientation_basis"] == (
        "study_wide_non_palindromic_consensus_or_internal_eaf"
    )


def test_step_four_collapses_consistent_reverse_complement_duplicates(tmp_path):
    frame = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [100, 100],
        "EA": ["A", "T"],
        "OA": ["G", "C"],
        "EAF": [0.30, 0.30],
        "EFFECT": [0.20, 0.20],
        "Z": [2.0, 2.0],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = _reject_collector(frame, tmp_path)
    columns = dict(
        COLUMNS,
        gwas_outputname="study",
        output_folder=str(tmp_path),
    )

    out, qc, _ = harmonise_allele_frequency(
        "1",
        frame,
        columns,
        _step_four_layout(),
        "\t",
        external_eaf_colmap=MAPPING,
        default_eaf_file=_reference(tmp_path),
        default_eaf_column="EUR",
        default_eaf_colmap=MAPPING,
        study_decision={
            "strand": "mixed",
            "eaf_is_maf": False,
            "effect_type": "beta",
        },
        policies=default_policies(),
        rejects=collector,
    )

    assert out.height == 1
    assert out.select(["EA", "OA", "EFFECT", "Z", "EAF"]).row(0) == (
        "A", "G", 0.20, 2.0, 0.30,
    )
    assert qc["post_orientation_duplicates"]["consistent_groups"] == 1
    assert qc["post_orientation_duplicates"]["rows_removed"] == 1
    assert collector.counts() == {"duplicate_variant": 1}
    rejected = collector.frame()
    assert rejected.select([
        "EA", "OA", "reject_step", "reject_reason",
    ]).row(0) == (
        "T", "C", "04 strand_orientation", "duplicate_variant",
    )

    report = pl.read_csv(
        tmp_path / "qc" / "study_duplicates_chr1.tsv",
        separator="\t",
    )
    assert report["strand_action"].to_list() == [
        "forward", "reverse_complement",
    ]
    assert report["duplicate_action"].to_list() == ["kept", "removed"]


def test_step_four_waits_for_deferred_effect_frequency_alignment(tmp_path):
    frame = pl.DataFrame({
        "CHR": ["1", "1", "1"],
        "POS": [100, 100, 200],
        "EA": ["A", "C", "C"],
        "OA": ["G", "T", "T"],
        "EAF": [0.30, 0.70, 0.72],
        "EFFECT": [0.20, -0.20, 0.40],
        "Z": [2.0, -2.0, 4.0],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = _reject_collector(frame, tmp_path)
    columns = dict(
        COLUMNS,
        gwas_outputname="study",
        output_folder=str(tmp_path),
    )
    policies = default_policies().with_overrides({
        "eaf": {"maf_reference_min_overlap": 2},
    })

    out, qc, _ = harmonise_allele_frequency(
        "1",
        frame,
        columns,
        _step_four_layout(),
        "\t",
        external_eaf_colmap=MAPPING,
        default_eaf_file=_reference(tmp_path),
        default_eaf_column="EUR",
        default_eaf_colmap=MAPPING,
        study_decision={
            "strand": "mixed",
            "eaf_is_maf": True,
            "effect_type": "beta",
        },
        policies=policies,
        rejects=collector,
    )

    assert qc["maf_reference_decision"] == "eaf"
    assert qc["post_orientation_duplicates"]["consistent_groups"] == 1
    assert out.height == 2
    assert out.select(["EA", "OA", "EAF", "EFFECT", "Z"]).row(0) == (
        "A", "G", 0.30, 0.20, 2.0,
    )
    assert collector.counts() == {"duplicate_variant": 1}
    report = pl.read_csv(
        tmp_path / "qc" / "study_duplicates_chr1.tsv",
        separator="\t",
    )
    assert report["strand_action"].to_list() == [
        "forward", "reverse_complement_swapped",
    ]
    assert report["EAF"].to_list() == pytest.approx([0.30, 0.30])
    assert report.item(0, "EAF") != report.item(1, "EAF")
    assert report["duplicate_action"].to_list() == ["kept", "removed"]


def test_step_four_rejects_conflicting_reverse_complement_duplicates(tmp_path):
    frame = pl.DataFrame({
        "CHR": ["1", "1", "1"],
        "POS": [100, 100, 200],
        "EA": ["A", "T", "T"],
        "OA": ["G", "C", "C"],
        "EAF": [0.30, 0.30, 0.28],
        "EFFECT": [0.20, 0.25, 0.40],
        "Z": [2.0, 2.5, 4.0],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = _reject_collector(frame, tmp_path)
    columns = dict(
        COLUMNS,
        gwas_outputname="study",
        output_folder=str(tmp_path),
    )

    out, qc, _ = harmonise_allele_frequency(
        "1",
        frame,
        columns,
        _step_four_layout(),
        "\t",
        external_eaf_colmap=MAPPING,
        default_eaf_file=_reference(tmp_path),
        default_eaf_column="EUR",
        default_eaf_colmap=MAPPING,
        study_decision={
            "strand": "mixed",
            "eaf_is_maf": False,
            "effect_type": "beta",
        },
        policies=default_policies(),
        rejects=collector,
    )

    assert out.get_column("POS").to_list() == [200]
    assert qc["post_orientation_duplicates"]["conflicting_groups"] == 1
    assert qc["post_orientation_duplicates"]["conflicting_rows_removed"] == 2
    assert collector.counts() == {"conflicting_duplicate": 2}
    report = pl.read_csv(
        tmp_path / "qc" / "study_duplicates_chr1.tsv",
        separator="\t",
    )
    assert report["duplicate_class"].unique().to_list() == ["conflicting"]
    assert report["duplicate_conflicting_fields"].unique().to_list() == [
        "beta,zscore"
    ]
    assert report["duplicate_action"].unique().to_list() == ["remove_all"]


def test_step_four_duplicate_conflict_can_fail_after_writing_evidence(tmp_path):
    frame = pl.DataFrame({
        "CHR": ["1", "1", "1"],
        "POS": [100, 100, 200],
        "EA": ["A", "T", "T"],
        "OA": ["G", "C", "C"],
        "EAF": [0.30, 0.30, 0.28],
        "EFFECT": [0.20, 0.25, 0.40],
        "Z": [2.0, 2.5, 4.0],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = _reject_collector(frame, tmp_path)
    columns = dict(
        COLUMNS,
        gwas_outputname="study",
        output_folder=str(tmp_path),
    )
    policies = default_policies().with_overrides({
        "duplicates.conflicting_action": "fail_dataset",
    })

    with pytest.raises(
        AlleleFrequencyError,
        match="unsafe reference-aligned duplicate group",
    ):
        harmonise_allele_frequency(
            "1",
            frame,
            columns,
            _step_four_layout(),
            "\t",
            external_eaf_colmap=MAPPING,
            default_eaf_file=_reference(tmp_path),
            default_eaf_column="EUR",
            default_eaf_colmap=MAPPING,
            study_decision={
                "strand": "mixed",
                "eaf_is_maf": False,
                "effect_type": "beta",
            },
            policies=policies,
            rejects=collector,
        )

    assert collector.counts() == {"conflicting_duplicate": 2}
    report = pl.read_csv(
        tmp_path / "qc" / "study_duplicates_chr1.tsv",
        separator="\t",
    )
    assert report.height == 2
    assert report["duplicate_action"].unique().to_list() == ["fail_dataset"]


def test_step_four_keeps_distinct_alternates_at_one_coordinate(tmp_path):
    reference = tmp_path / "multiallelic_reference.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.30\n"
        "1\t100\tG\tT\t0.20\n",
        encoding="utf-8",
    )
    frame = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [100, 100],
        "EA": ["A", "T"],
        "OA": ["G", "G"],
        "EAF": [0.30, 0.20],
        "EFFECT": [0.20, -0.10],
        "Z": [2.0, -1.0],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = _reject_collector(frame, tmp_path)
    columns = dict(
        COLUMNS,
        gwas_outputname="study",
        output_folder=str(tmp_path),
    )

    out, qc, _ = harmonise_allele_frequency(
        "1",
        frame,
        columns,
        _step_four_layout(),
        "\t",
        external_eaf_colmap=MAPPING,
        default_eaf_file=str(reference),
        default_eaf_column="EUR",
        default_eaf_colmap=MAPPING,
        study_decision={
            "strand": "forward",
            "eaf_is_maf": False,
            "effect_type": "beta",
        },
        policies=default_policies(),
        rejects=collector,
    )

    assert out.select(["EA", "OA"]).rows() == [("A", "G"), ("T", "G")]
    assert qc["post_orientation_duplicates"]["duplicate_groups"] == 0
    assert collector.total() == 0
    report = pl.read_csv(
        tmp_path / "qc" / "study_duplicates_chr1.tsv",
        separator="\t",
    )
    assert report.is_empty()


def test_palindromic_variant_requires_strong_study_consensus(tmp_path):
    frame = _study().filter(pl.col("POS") == 600)
    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "mixed", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.is_empty()
    assert qc["palindromic_ambiguous"] == 1
    assert qc["palindromic_orientation_unavailable"] == 0


def test_palindromic_variant_without_consensus_or_internal_eaf_is_unavailable(
    tmp_path,
):
    frame = _study().filter(pl.col("POS") == 500).drop("EAF")
    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "mixed", "eaf_is_maf": None, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.is_empty()
    assert qc["palindromic_orientation_unavailable"] == 1
    assert qc["palindromic_frequency_discordant"] == 0
    assert qc["palindromic_ambiguous"] == 0


def test_internal_effect_allele_frequency_resolves_palindrome_without_consensus(tmp_path):
    frame = _study().filter(pl.col("POS") == 500)
    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "mixed", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.height == 1
    assert out.item(0, "strand_action") == "reverse_complement"
    assert out.item(0, "EA") == "G"
    assert out.item(0, "OA") == "C"
    assert out.item(0, "EAF") == pytest.approx(0.11)
    assert out.item(0, "EFFECT") == pytest.approx(0.12)
    assert qc["palindromic_resolution_counts"] == {"internal_eaf": 1}
    assert qc["palindromic_frequency_state_counts"] == {"resolved": 1}


def test_internal_frequency_conflicting_with_strong_consensus_is_rejected(tmp_path):
    frame = _study().filter(pl.col("POS") == 500)
    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "forward", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.is_empty()
    assert qc["palindromic_frequency_conflict"] == 1
    assert qc["palindromic_ambiguous"] == 0


def test_internal_frequency_confirming_strong_consensus_is_recorded(tmp_path):
    frame = _study().filter(pl.col("POS") == 500)
    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "reverse", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.height == 1
    assert out.item(0, "strand_action") == "reverse_complement"
    assert qc["palindromic_resolution_counts"] == {
        "consensus_and_internal_eaf": 1
    }


@pytest.mark.parametrize(
    ("overrides", "expected_state"),
    [
        (
            {
                "strand.palindromic_af_ambiguity_lower": 0.10,
                "strand.palindromic_af_ambiguity_upper": 0.90,
            },
            "ambiguous_frequency_band",
        ),
        (
            {"strand.palindromic_af_max_difference": 0.005},
            "maximum_difference_exceeded",
        ),
        (
            {"strand.palindromic_af_min_error_margin": 0.80},
            "insufficient_error_margin",
        ),
    ],
)
def test_palindromic_frequency_rules_use_yaml_overrides(
    tmp_path, overrides, expected_state,
):
    frame = _study().filter(pl.col("POS") == 500)
    policies = default_policies().with_overrides(overrides)

    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "mixed", "eaf_is_maf": False, "effect_type": "beta"},
        policies=policies,
    )

    assert out.is_empty()
    assert qc["palindromic_frequency_state_counts"] == {expected_state: 1}
    if expected_state == "maximum_difference_exceeded":
        assert qc["palindromic_frequency_discordant"] == 1
        assert qc["palindromic_ambiguous"] == 0
    else:
        assert qc["palindromic_frequency_discordant"] == 0
        assert qc["palindromic_ambiguous"] == 1


def test_palindromic_frequency_boundary_is_ambiguous_without_consensus(tmp_path):
    frame = _study().filter(pl.col("POS") == 500).with_columns(
        pl.lit(0.40).alias("EAF")
    )
    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "mixed", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.is_empty()
    assert qc["palindromic_frequency_state_counts"] == {
        "ambiguous_frequency_band": 1
    }
    assert qc["palindromic_ambiguous"] == 1


def test_strong_consensus_orients_near_half_palindromic_frequency(tmp_path):
    frame = _study().filter(pl.col("POS") == 600)
    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={"strand": "reverse", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )

    assert out.height == 1
    assert out.item(0, "strand_action") == "reverse_complement"
    assert out.item(0, "EAF") == pytest.approx(0.49)
    assert qc["palindromic_ambiguous"] == 0


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


@pytest.mark.parametrize(
    ("study_label", "reference_label"),
    [("1", "01"), ("X", "23")],
)
def test_strand_reference_uses_shared_chromosome_policies(
    tmp_path, study_label, reference_label,
):
    reference = tmp_path / ("strand_%s.tsv" % reference_label)
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n%s\t100\tG\tA\t0.30\n"
        % reference_label,
        encoding="utf-8",
    )
    frame = _study().filter(pl.col("POS") == 100).with_columns(
        pl.lit(study_label).alias("CHR")
    )

    out, qc, _ = harmonise_strand_orientation(
        study_label,
        frame,
        dict(COLUMNS),
        str(reference),
        "EUR",
        MAPPING,
        study_decision={
            "strand": "forward", "eaf_is_maf": False, "effect_type": "beta",
        },
        policies=default_policies(),
    )

    assert out.height == 1
    assert out.item(0, "CHR") == study_label
    assert qc["actions"] == {"forward": 1}


def test_strand_reference_rejects_position_cast_that_loses_every_value(tmp_path):
    reference = tmp_path / "invalid_positions.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n1\tnot-a-position\tG\tA\t0.30\n",
        encoding="utf-8",
    )

    with pytest.raises(
        StrandOrientationError,
        match="Every strand reference position value became empty",
    ):
        harmonise_strand_orientation(
            "1", _study().filter(pl.col("POS") == 100), dict(COLUMNS),
            str(reference), "EUR", MAPPING,
            study_decision={
                "strand": "forward", "eaf_is_maf": False,
                "effect_type": "beta",
            },
            policies=default_policies(),
        )


def test_strand_study_rejects_position_cast_that_loses_every_value(tmp_path):
    frame = _study().filter(pl.col("POS") == 100).with_columns(
        pl.lit("not-a-position").alias("POS")
    )

    with pytest.raises(
        StrandOrientationError,
        match="Every study before strand-reference matching position value became empty",
    ):
        harmonise_strand_orientation(
            "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
            study_decision={
                "strand": "forward", "eaf_is_maf": False,
                "effect_type": "beta",
            },
            policies=default_policies(),
        )


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


def test_strand_refuses_raw_odds_ratio_before_reference_or_allele_changes(
    tmp_path,
    monkeypatch,
):
    frame = _study().filter(pl.col("POS") == 200).with_columns([
        pl.lit(2.0).alias("EFFECT"),
        pl.lit(0.20).alias("EAF"),
    ])
    reference = _reference(tmp_path)
    reference_read = False

    original = strand_module._load_reference

    def observe_reference_read(*args, **kwargs):
        nonlocal reference_read
        reference_read = True
        return original(*args, **kwargs)

    monkeypatch.setattr(strand_module, "_load_reference", observe_reference_read)

    with pytest.raises(
        StrandOrientationError,
        match="execution-order error.*harmonise_effect_estimates",
    ):
        harmonise_strand_orientation(
            "1", frame, dict(COLUMNS), reference, "EUR", MAPPING,
            study_decision={
                "strand": "forward",
                "eaf_is_maf": True,
                "effect_type": "odds_ratio",
            },
            policies=default_policies(),
        )

    assert reference_read is False
    assert frame.item(0, "EA") == "C"
    assert frame.item(0, "OA") == "T"
    assert frame.item(0, "EFFECT") == pytest.approx(2.0)
    assert frame.item(0, "EAF") == pytest.approx(0.20)


def test_internal_strand_function_is_not_exported():
    assert "harmonise_strand_orientation" not in strand_module.__all__


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


def test_fail_policy_allows_a_fully_reference_matched_chromosome(tmp_path):
    policies = default_policies().with_overrides({
        "strand": {"unmatched_action": "fail"},
    })

    out, qc, _ = harmonise_strand_orientation(
        "1", _study().filter(pl.col("POS") == 100), dict(COLUMNS),
        _reference(tmp_path), "EUR", MAPPING,
        study_decision={
            "strand": "forward",
            "eaf_is_maf": False,
            "effect_type": "beta",
        },
        policies=policies,
    )

    assert out.height == 1
    assert qc["reference_unmatched_action"] == "fail"
    assert qc["reference_unmatched_detected"] == 0
    assert qc["reference_unmatched_retained"] == 0


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
    assert qc["reference_unmatched_action"] == "reject"
    assert qc["reference_unmatched_detected"] == 1
    assert qc["reference_unmatched"] == 1
    assert qc["reference_unmatched_retained"] == 0


def test_retain_policy_preserves_forward_unmatched_variant_for_fasta_gate(
    tmp_path,
):
    frame = _study().filter(pl.col("POS") == 100).with_columns(
        pl.lit(999).alias("POS")
    )
    policies = default_policies().with_overrides({
        "strand": {"unmatched_action": "retain"},
    })

    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={
            "strand": "forward",
            "eaf_is_maf": False,
            "effect_type": "beta",
        },
        policies=policies,
    )

    assert out.height == 1
    assert out.item(0, "EA") == "A"
    assert out.item(0, "OA") == "G"
    assert out.item(0, "EFFECT") == pytest.approx(0.20)
    assert out.item(0, "Z") == pytest.approx(2.0)
    assert out.item(0, "EAF") == pytest.approx(0.31)
    assert out.item(0, REFERENCE_AF_COLUMN) is None
    assert out.item(0, "strand_action") == "reference_unmatched_retained"
    assert qc["reference_unmatched_action"] == "retain"
    assert qc["reference_unmatched_detected"] == 1
    assert qc["reference_unmatched"] == 0
    assert qc["reference_unmatched_retained"] == 1


def test_retain_policy_applies_strong_reverse_consensus_without_swapping(
    tmp_path,
):
    frame = _study().filter(pl.col("POS") == 100).with_columns(
        pl.lit(999).alias("POS")
    )
    policies = default_policies().with_overrides({
        "strand": {"unmatched_action": "retain"},
    })

    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={
            "strand": "reverse",
            "eaf_is_maf": False,
            "effect_type": "beta",
        },
        policies=policies,
    )

    assert out.height == 1
    assert out.item(0, "EA") == "T"
    assert out.item(0, "OA") == "C"
    assert out.item(0, "EFFECT") == pytest.approx(0.20)
    assert out.item(0, "Z") == pytest.approx(2.0)
    assert out.item(0, "EAF") == pytest.approx(0.31)
    assert out.item(0, REFERENCE_AF_COLUMN) is None
    assert out.item(0, "strand_action") == (
        "reference_unmatched_retained_reverse_complement"
    )
    assert qc["reference_unmatched_detected"] == 1
    assert qc["reference_unmatched"] == 0
    assert qc["reference_unmatched_retained"] == 1


def test_retain_policy_does_not_guess_unmatched_palindromic_orientation(
    tmp_path,
):
    frame = _study().filter(pl.col("POS") == 600).with_columns(
        pl.lit(999).alias("POS")
    )
    policies = default_policies().with_overrides({
        "strand": {"unmatched_action": "retain"},
    })

    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={
            "strand": "mixed",
            "eaf_is_maf": False,
            "effect_type": "beta",
        },
        policies=policies,
    )

    assert out.is_empty()
    assert qc["reference_unmatched_detected"] == 1
    assert qc["reference_unmatched"] == 0
    assert qc["reference_unmatched_retained"] == 0
    assert qc["palindromic_orientation_unavailable"] == 1


def test_indels_use_direct_matching_but_never_reverse_complement(tmp_path):
    reference = tmp_path / "indel.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n1\t800\tA\tAG\t0.20\n",
        encoding="utf-8",
    )
    base = {
        "CHR": ["1"], "POS": [800], "EAF": [0.20],
        "EFFECT": [0.10], "Z": [1.0],
    }
    direct = pl.DataFrame(dict(base, EA=["AG"], OA=["A"]))
    out, _, _ = harmonise_strand_orientation(
        "1", direct, dict(COLUMNS), str(reference), "EUR", MAPPING,
        study_decision={"strand": "reverse", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )
    assert out.height == 1
    assert out.item(0, "strand_action") == "forward"

    reverse_complement_only = pl.DataFrame(dict(base, EA=["CT"], OA=["T"]))
    out, qc, _ = harmonise_strand_orientation(
        "1", reverse_complement_only, dict(COLUMNS), str(reference), "EUR", MAPPING,
        study_decision={"strand": "reverse", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )
    assert out.is_empty()
    assert qc["reference_unmatched"] == 1


def test_reference_aligned_mode_still_validates_every_row_against_reference(tmp_path):
    policies = default_policies().with_overrides({
        "strand": {"mode": "reference_aligned"}
    })
    decision = resolve_strand_consensus({}, "GRCh37", policies)
    out, qc, _ = harmonise_strand_orientation(
        "1", _study(), dict(COLUMNS), _reference(tmp_path), "EUR", MAPPING,
        study_decision={
            "strand": "reference_aligned",
            "eaf_is_maf": False,
            "effect_type": "beta",
        },
        policies=policies,
    )

    assert decision["strand"] == "reference_aligned"
    assert out.get_column("POS").to_list() == [100, 200, 500, 600]
    assert out.get_column("strand_action").to_list() == [
        "forward", "forward_swapped", "forward_swapped", "forward_swapped",
    ]
    assert qc["status"] == "success"
    assert qc["strand_mode"] == "reference_aligned"
    assert qc["reference_unmatched"] == 3
    assert qc["palindromic_orientation_basis"] == (
        "reference_aligned_declaration"
    )


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
    logger = PipelineLogger("study", "1", str(tmp_path), screen_level="ERROR")
    frame = _study().filter(pl.col("POS") == 100).with_columns(
        pl.lit(0.90).alias("EAF")
    )
    columns = dict(COLUMNS, gwas_outputname="study", output_folder=str(tmp_path))
    layout = {
        "frequency_qc_directory": "eaf_qc",
        "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
        "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
        "post_orientation_duplicates": "qc/{dataset_id}_duplicates_chr{chromosome}.tsv",
    }
    out, qc, _ = harmonise_allele_frequency(
        "1", frame, columns, layout, "\t",
        external_eaf_colmap=MAPPING,
        default_eaf_file=_reference(tmp_path),
        default_eaf_column="EUR",
        default_eaf_colmap=MAPPING,
        study_decision={"strand": "forward", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
        logger=logger,
    )

    actions = logger.summary()["qc_actions"]
    logger.close()
    warning = next(action for action in actions if action["check"] == "non-palindromic study/reference frequency concordance")
    assert warning["matched"] == 1
    assert warning["changed"] == 0
    assert warning["removed"] == 0
    assert warning["outcome"] == "retained_without_frequency_change"
    assert out.height == 1
    assert "zmaf" not in out.columns
    assert out.item(0, "strand_af_difference") == pytest.approx(0.60)
    assert qc["strand_af_comparable"] == 1
    assert qc["strand_af_discordant"] == 1
    assert qc["strand_non_palindromic_af_discordant"] == 1
    assert qc["strand_palindromic_af_discordant"] == 0


def test_palindromic_af_discordance_is_rejected_after_consensus_orientation(tmp_path):
    frame = _study().filter(pl.col("POS") == 600).with_columns(
        pl.lit(0.10).alias("EAF")
    )
    columns = dict(COLUMNS, gwas_outputname="study", output_folder=str(tmp_path))
    layout = {
        "frequency_qc_directory": "eaf_qc",
        "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
        "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
        "post_orientation_duplicates": "qc/{dataset_id}_duplicates_chr{chromosome}.tsv",
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

    assert out.is_empty()
    assert qc["strand_palindromic_af_comparable"] == 1
    assert qc["strand_palindromic_af_discordant"] == 1
    assert qc["strand_palindromic_af_discordance_action"] == "reject"


def test_external_palindromic_eaf_is_aligned_after_consensus_without_beta_flip(tmp_path):
    external = tmp_path / "external_eaf.tsv"
    external.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n1\t500\tC\tG\t0.10\n",
        encoding="utf-8",
    )
    frame = _study().filter(pl.col("POS") == 500).drop("EAF")
    columns = dict(
        COLUMNS,
        eaf_col=None,
        eafcolumn="EUR",
        gwas_outputname="study",
        output_folder=str(tmp_path),
    )
    layout = {
        "frequency_qc_directory": "eaf_qc",
        "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
        "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
        "post_orientation_duplicates": "qc/{dataset_id}_duplicates_chr{chromosome}.tsv",
    }
    policies = default_policies().with_overrides({
        # This one-row fixture exercises palindromic alignment, not a
        # chromosome-scale MAF screen. Disable suspicion at the strict upper
        # boundary so the test isolates its intended scientific behavior.
        "eaf.maf_decision_cutoff": 1.0,
    })

    out, qc, _ = harmonise_allele_frequency(
        "1", frame, columns, layout, "\t",
        eaffile=str(external),
        external_eaf_colmap=MAPPING,
        default_eaf_file=_reference(tmp_path),
        default_eaf_column="EUR",
        default_eaf_colmap=MAPPING,
        study_decision={"strand": "reverse", "eaf_is_maf": None, "effect_type": "beta"},
        policies=policies,
    )

    assert out.height == 1
    assert out.item(0, "strand_action") == "reverse_complement"
    assert out.item(0, "EFFECT") == pytest.approx(0.12)
    assert out.item(0, "Z") == pytest.approx(1.2)
    assert out.item(0, "EUR") == pytest.approx(0.10)
    assert qc["external_merge_palindromic_orientation_resolved_rows"] == 1


def test_strand_reference_duplicate_policy_keeps_exact_and_discards_conflict(tmp_path):
    exact = tmp_path / "exact.tsv"
    exact.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.30\n"
        "1\t100\tG\tA\t0.30\n",
        encoding="utf-8",
    )
    frame = _study().filter(pl.col("POS") == 100)
    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), str(exact), "EUR", MAPPING,
        study_decision={"strand": "forward", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )
    assert out.height == 1
    assert qc["reference_exact_duplicate_groups"] == 1
    assert qc["reference_non_identical_duplicate_groups"] == 0

    conflict = tmp_path / "conflict.tsv"
    conflict.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.30\n"
        "1\t100\tG\tA\t0.45\n",
        encoding="utf-8",
    )
    out, qc, _ = harmonise_strand_orientation(
        "1", frame, dict(COLUMNS), str(conflict), "EUR", MAPPING,
        study_decision={"strand": "forward", "eaf_is_maf": False, "effect_type": "beta"},
        policies=default_policies(),
    )
    assert out.is_empty()
    assert qc["reference_non_identical_duplicate_groups"] == 1
    assert qc["reference_duplicate_groups_discarded"] == 1
    assert qc["reference_unmatched"] == 1


def test_maf_like_column_proven_to_be_eaf_is_then_inverted_for_swapped_row(tmp_path):
    frame = _study().filter(pl.col("POS").is_in([100, 200]))
    # Reference ALT frequency is 0.28, so the original REF effect allele has
    # frequency 0.72. After the proof that this is EAF, standardized ALT EAF is 0.28.
    columns = dict(COLUMNS, gwas_outputname="study", output_folder=str(tmp_path))
    layout = {
        "frequency_qc_directory": "eaf_qc",
        "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
        "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
        "post_orientation_duplicates": "qc/{dataset_id}_duplicates_chr{chromosome}.tsv",
    }
    policies = default_policies().with_overrides({
        "eaf": {"maf_reference_min_overlap": 2},
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
    swapped = out.filter(pl.col("POS") == 200)
    assert swapped.item(0, "EAF") == pytest.approx(0.28)
    assert swapped.item(0, "strand_af_difference") == pytest.approx(0.0)
