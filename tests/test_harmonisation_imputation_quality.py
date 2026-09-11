"""Dataset-wide INFO/MaCH Rsq detection and chromosome range handling."""

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.imputation_quality import (
    ImputationQualityError,
    harmonise_imputation_quality,
    resolve_dataset_info_score_type,
    resolve_info_score_type_from_batches,
)
from postgwas.modules.harmonisation.policies import PolicyError, default_policies
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    SOURCE_INPUT_ROW_COLUMN,
)
from postgwas.modules.harmonisation.summary_statistics_io import (
    normalise_imputation_quality_column,
)


class _CaptureLogger:
    def __init__(self):
        self.infos = []
        self.warnings = []

    def info(self, message, indent=0):
        self.infos.append((message, indent))

    def warn(self, message, indent=0):
        self.warnings.append((message, indent))


def _study(values, *, resolved_type=None):
    mapping = {
        "chr_col": "CHR",
        "pos_col": "POS",
        "ea_col": "EA",
        "oa_col": "OA",
        "imp_info_col": "INFO",
        "info_source": "internal",
    }
    if resolved_type is not None:
        mapping["resolved_info_score_type"] = resolved_type
        mapping["info_score_type_source"] = "automatic_full_dataset"
    return (
        pl.DataFrame({
            "CHR": ["1"] * len(values),
            "POS": list(range(1, len(values) + 1)),
            "EA": ["A"] * len(values),
            "OA": ["G"] * len(values),
            "INFO": values,
        }),
        mapping,
    )


def test_info_score_type_defaults_and_relationships_are_schema_validated():
    policies = default_policies()

    assert policies.get("info.multi_value_delimiter") == ","
    assert policies.get("info.multi_value_aggregation") == "median"
    assert policies.get("info.multi_value_output_column") == (
        "__postgwas_multi_value_info"
    )
    assert policies.get("info.multi_value_invalid_token_action") == "ignore"
    assert policies.get("info.score_type") == "auto"
    assert policies.get("info.out_of_range") == "reject"
    assert policies.get("info.clip_tolerance") == pytest.approx(1.05)
    assert policies.get("info.mach_rsq_max") == pytest.approx(2.0)
    assert policies.get("info.auto_mach_rsq_fraction") == pytest.approx(0.001)
    assert policies.get("info.maximum_invalid_fraction") == pytest.approx(0.001)
    with pytest.raises(PolicyError, match="info.score_type"):
        policies.with_overrides({"info.score_type": "mach_r2"})
    with pytest.raises(PolicyError, match="info.clip_tolerance"):
        policies.with_overrides({"info.clip_tolerance": 0.9})
    with pytest.raises(PolicyError, match="info.clip_tolerance"):
        policies.with_overrides({
            "info.clip_tolerance": 2.0,
            "info.mach_rsq_max": 2.0,
        })
    with pytest.raises(PolicyError, match="info.multi_value_aggregation"):
        policies.with_overrides({"info.multi_value_aggregation": "weighted"})
    with pytest.raises(PolicyError, match="info.multi_value_delimiter"):
        policies.with_overrides({"info.multi_value_delimiter": "||"})
    with pytest.raises(PolicyError, match="info.multi_value_output_column"):
        policies.with_overrides({"info.multi_value_output_column": "INFO_MEAN"})
    with pytest.raises(
        PolicyError, match="info.multi_value_invalid_token_action",
    ):
        policies.with_overrides({
            "info.multi_value_invalid_token_action": "warn",
        })


def test_multi_value_info_uses_configured_unweighted_aggregation_and_logger():
    frame = pl.DataFrame({
        "INFO": ["0.9,0.7,0.1", "0.8", "NA,0.6,bad", None],
    })
    logger = _CaptureLogger()
    mapping = {"imp_info_col": "INFO", "info_source": "internal"}
    policies = default_policies().with_overrides({
        "info.multi_value_aggregation": "mean",
    })

    result, resolved_mapping = normalise_imputation_quality_column(
        frame,
        mapping,
        policies=policies,
        logger=logger,
    )

    output_column = "__postgwas_multi_value_info"
    assert result["INFO"].to_list() == frame["INFO"].to_list()
    assert result[output_column].to_list()[:3] == pytest.approx([
        (0.9 + 0.7 + 0.1) / 3.0,
        0.8,
        0.6,
    ])
    assert result[output_column].to_list()[3] is None
    assert resolved_mapping["imp_info_col"] == output_column
    assert resolved_mapping["info_multi_value_source_column"] == "INFO"
    assert resolved_mapping["info_multi_value_aggregation"] == "mean"
    assert any(
        "before duplicate validation" in message
        and "unweighted" in message
        for message, _indent in logger.infos
    )
    assert any(
        "invalid_token_action='ignore'" in message
        for message, _indent in logger.warnings
    )
    assert any(
        "1 missing scalar" in message
        for message, _indent in logger.warnings
    )


def test_multi_value_info_fail_policies_stop_before_scientific_use():
    frame = pl.DataFrame({"INFO": ["0.9,0.7", "0.8,bad"]})
    mapping = {"imp_info_col": "INFO", "info_source": "internal"}

    with pytest.raises(ValueError, match="multi_value_aggregation is 'fail'"):
        normalise_imputation_quality_column(
            frame,
            dict(mapping),
            policies=default_policies().with_overrides({
                "info.multi_value_aggregation": "fail",
            }),
        )

    with pytest.raises(
        ValueError, match="multi_value_invalid_token_action is 'fail'",
    ):
        normalise_imputation_quality_column(
            frame,
            dict(mapping),
            policies=default_policies().with_overrides({
                "info.multi_value_invalid_token_action": "fail",
            }),
        )


def test_multi_value_info_refuses_to_overwrite_configured_working_column():
    frame = pl.DataFrame({
        "INFO": ["0.9,0.7"],
        "__postgwas_multi_value_info": [0.1],
    })

    with pytest.raises(ValueError, match="already exists"):
        normalise_imputation_quality_column(
            frame,
            {"imp_info_col": "INFO", "info_source": "internal"},
            policies=default_policies(),
        )


def test_one_high_value_cannot_flip_a_large_standard_info_source():
    decision = resolve_info_score_type_from_batches(
        [pl.DataFrame({"INFO": [0.8] * 9_999 + [1.5]})],
        "INFO",
        policies=default_policies(),
    )

    assert decision["resolved_type"] == "standard_info"
    assert decision["mach_range_values"] == 1
    assert decision["mach_range_fraction"] == pytest.approx(0.0001)


def test_auto_selects_mach_rsq_at_configured_agreement_fraction():
    decision = resolve_info_score_type_from_batches(
        [pl.DataFrame({"INFO": [0.8] * 999 + [1.5]})],
        "INFO",
        policies=default_policies(),
    )

    assert decision["resolved_type"] == "mach_rsq"
    assert decision["mach_range_fraction"] == pytest.approx(0.001)


def test_values_above_mach_max_stop_only_above_configured_fraction():
    allowed = resolve_info_score_type_from_batches(
        [pl.DataFrame({"INFO": [0.8] * 999 + [87.0]})],
        "INFO",
        policies=default_policies(),
    )

    assert allowed["resolved_type"] == "standard_info"
    assert allowed["above_mach_maximum"] == 1
    assert allowed["invalid_fraction"] == pytest.approx(0.001)

    with pytest.raises(
        ImputationQualityError,
        match="stopped before chromosome processing.*external_info_column",
    ):
        resolve_info_score_type_from_batches(
            [pl.DataFrame({"INFO": [0.8] * 998 + [87.0, 99.0]})],
            "INFO",
            policies=default_policies(),
        )


def test_explicit_type_overrides_detection_but_not_gross_invalid_guard():
    policies = default_policies().with_overrides({
        "info.score_type": "standard_info",
    })
    decision = resolve_info_score_type_from_batches(
        [pl.DataFrame({"INFO": [0.8, 1.5, 1.6]})],
        "INFO",
        policies=policies,
    )

    assert decision["resolved_type"] == "standard_info"
    assert decision["decision_source"] == "explicit_configuration"
    with pytest.raises(ImputationQualityError, match="info.mach_rsq_max"):
        resolve_info_score_type_from_batches(
            [pl.DataFrame({"INFO": [87.0, 99.0]})],
            "INFO",
            policies=policies,
        )


def test_chromosome_worker_requires_dataset_decision_and_preserves_mach_rsq():
    frame, mapping = _study([0.9, 1.02, 1.5])
    with pytest.raises(ImputationQualityError, match="complete-dataset level"):
        harmonise_imputation_quality(
            "1", frame, mapping, policies=default_policies(),
        )

    standard_frame, standard_mapping = _study(
        [0.9, 1.02, 1.5], resolved_type="standard_info",
    )
    standard, standard_qc, _ = harmonise_imputation_quality(
        "1", standard_frame, standard_mapping, policies=default_policies(),
    )
    assert standard["INFO"].to_list() == pytest.approx([0.9, 1.0])
    assert standard_qc["rejected_out_of_range"] == 1

    mach_frame, mach_mapping = _study(
        [0.9, 1.02, 1.5], resolved_type="mach_rsq",
    )
    mach, mach_qc, _ = harmonise_imputation_quality(
        "1", mach_frame, mach_mapping, policies=default_policies(),
    )
    assert mach["INFO"].to_list() == pytest.approx([0.9, 1.02, 1.5])
    assert mach_qc["rescaled_within_tolerance"] == 0
    assert mach_qc["rejected_out_of_range"] == 0


@pytest.mark.parametrize(
    (
        "resolved_type",
        "values",
        "retained_positions",
        "out_of_range_positions",
    ),
    [
        ("standard_info", [None, -0.1, 0.8, 1.2], [3], {2, 4}),
        ("mach_rsq", [None, -0.1, 0.8, 1.2, 2.1], [3, 4], {2, 5}),
    ],
)
def test_nullified_out_of_range_info_retains_its_rejection_reason(
    tmp_path,
    resolved_type,
    values,
    retained_positions,
    out_of_range_positions,
):
    frame, mapping = _study(values, resolved_type=resolved_type)
    frame = frame.with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    reject_path = tmp_path / ("%s_rejected.tsv" % resolved_type)
    collector = RejectCollector(
        source_snapshot=frame,
        logger=None,
        out_path=reject_path,
        delimiter="\t",
        compress=False,
    )
    policies = default_policies().with_overrides({
        "info": {
            "score_type": resolved_type,
            "out_of_range": "null",
            "on_missing": "reject",
        }
    })

    retained, qc, _ = harmonise_imputation_quality(
        "1",
        frame,
        mapping,
        policies=policies,
        rejects=collector,
    )
    collector.flush()
    rejected = pl.read_csv(reject_path, separator="\t")
    reasons = dict(zip(rejected["POS"], rejected["reject_reason"]))
    details = dict(zip(rejected["POS"], rejected["reject_detail"]))

    assert retained["POS"].to_list() == retained_positions
    assert reasons[1] == "info_missing"
    assert {
        position for position, reason in reasons.items()
        if reason == "info_out_of_range"
    } == out_of_range_positions
    assert all(
        "info.out_of_range='null'" in details[position]
        and "info.on_missing='reject'" in details[position]
        for position in out_of_range_positions
    )
    assert qc["missing_info"] == 1
    assert qc["out_of_range"] == len(out_of_range_positions)
    assert qc["nullified_out_of_range"] == len(out_of_range_positions)
    assert qc["missing_info_after"] == 1 + len(out_of_range_positions)
    assert qc["rejected_out_of_range"] == len(out_of_range_positions)
    assert qc["rejected_missing"] == 1


def test_missing_info_failure_reports_original_and_nullified_origins():
    frame, mapping = _study(
        [None, -0.1, 0.8, 1.2], resolved_type="standard_info",
    )
    policies = default_policies().with_overrides({
        "info": {
            "score_type": "standard_info",
            "out_of_range": "null",
            "on_missing": "fail",
        }
    })

    with pytest.raises(
        ImputationQualityError,
        match=(
            "3 variants have no imputation quality after range handling: "
            "1 were already missing or non-finite and 2 were set to missing"
        ),
    ):
        harmonise_imputation_quality(
            "1", frame, mapping, policies=policies,
        )


def test_external_detection_scans_relevant_chromosomes_as_one_dataset(tmp_path):
    first = tmp_path / "info_chr1.tsv"
    second = tmp_path / "info_chr2.tsv"
    first.write_text(
        "CHROM\tINFO\n1\t0.8\n1\t0.9\n3\t99\n",
        encoding="utf-8",
    )
    second.write_text(
        "CHROM\tINFO\n2\t0.7\n2\t1.5\n3\t99\n",
        encoding="utf-8",
    )
    settings = (
        load_configuration().modules.harmonisation
        .external_reference_staging.model_dump()
    )
    policies = default_policies().with_overrides({
        "info.auto_mach_rsq_fraction": 0.25,
    })

    decision = resolve_dataset_info_score_type(
        df=None,
        sample_column_dict={"info_source": "external"},
        resource_maps={
            "1": {"user_info_file": str(first)},
            "2": {"user_info_file": str(second)},
        },
        chromosomes=["1", "2"],
        external_info_column="INFO",
        external_info_colmap={"chr": "CHROM", "delimiter": "tab"},
        external_reference_staging=settings,
        policies=policies,
    )

    assert decision["finite_values"] == 4
    assert decision["mach_range_values"] == 1
    assert decision["above_mach_maximum"] == 0
    assert decision["resolved_type"] == "mach_rsq"
