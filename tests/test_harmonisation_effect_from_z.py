"""Scientific and numerical-safety tests for effects derived from Z."""

import math

import polars as pl
import pytest

from postgwas.modules.harmonisation import effect_from_z
from postgwas.modules.harmonisation.effect_from_z import (
    BetaZSignMismatchError,
    XChromosomeZOnlyReconstructionError,
    derive_effect_and_standard_error_from_z,
)
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    SOURCE_INPUT_ROW_COLUMN,
)


class _CaptureLogger:
    def __init__(self):
        self.infos = []
        self.warnings = []

    def settings(self, _keys):
        pass

    def info(self, message):
        self.infos.append(message)

    def warn(self, message):
        self.warnings.append(message)

    def qc(self, *_args, **_kwargs):
        pass


def test_partial_supplied_se_uses_z_only_for_null_cells():
    result, qc, mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({
            "BETA": [0.1, 0.2, 0.3],
            "SE": [0.02, None, None],
            "Z": [5.0, 4.0, 0.0],
        }),
        {
            "beta_or_col": "BETA",
            "beta_col": "BETA",
            "se_col": "SE",
            "imp_z_col": "Z",
        },
        policies=default_policies(),
    )

    assert result["SE"].to_list() == [0.02, 0.05, None]
    assert result[effect_from_z.SE_UNAVAILABLE_FROM_Z_COLUMN].to_list() == [
        False, False, True,
    ]
    assert mapping["se_col"] == "SE"
    assert qc["se_missing_on_entry"] == 2
    assert qc["se_recovered_from_z"] == 1
    assert qc["se_unavailable_from_z"] == 1


def test_missing_se_is_positive_for_consistent_positive_and_negative_pairs():
    result, qc, mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({
            "BETA": [0.2, -0.2],
            "Z": [2.0, -2.0],
        }),
        {
            "beta_or_col": "BETA",
            "beta_col": "BETA",
            "se_col": None,
            "imp_z_col": "Z",
        },
        policies=default_policies(),
    )

    assert result["SE"].to_list() == pytest.approx([0.1, 0.1])
    assert mapping["se_col"] == "SE"
    assert qc["se_recovered_from_z"] == 2
    assert qc["beta_z_sign_mismatches"] == 0
    assert qc["removed_beta_z_sign_mismatches"] == 0


def test_default_rejects_beta_z_sign_mismatches_with_provenance(tmp_path):
    source = pl.DataFrame({
        "POS": [1, 2, 3, 4],
        "BETA": [0.2, -0.2, 0.3, 0.0],
        "Z": [2.0, -2.0, -3.0, 4.0],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = RejectCollector(
        source_snapshot=source,
        logger=None,
        out_path=tmp_path / "rejected.tsv",
        delimiter="\t",
        compress=False,
    )

    result, qc, _mapping = derive_effect_and_standard_error_from_z(
        "1",
        source,
        {
            "beta_or_col": "BETA",
            "beta_col": "BETA",
            "se_col": None,
            "imp_z_col": "Z",
        },
        policies=default_policies(),
        rejects=collector,
    )
    rejected = collector.frame()

    assert result["POS"].to_list() == [1, 2]
    assert result["SE"].to_list() == pytest.approx([0.1, 0.1])
    assert rejected["POS"].to_list() == ["3", "4"]
    assert rejected["reject_reason"].to_list() == [
        "beta_z_sign_discordant",
        "beta_z_sign_discordant",
    ]
    assert all(
        "effect_from_z.beta_z_sign_mismatch='reject'" in detail
        for detail in rejected["reject_detail"]
    )
    assert qc["beta_z_sign_mismatch_action"] == "reject"
    assert qc["beta_z_sign_mismatches"] == 2
    assert qc["removed_beta_z_sign_mismatches"] == 2
    assert qc["variants_removed_total"] == 2


def test_sign_policy_applies_only_to_missing_cells_in_supplied_se_column():
    result, qc, _mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({
            "POS": [1, 2, 3],
            "BETA": [0.2, 0.3, -0.4],
            "SE": [0.02, None, None],
            "Z": [-2.0, -3.0, -4.0],
        }),
        {
            "beta_or_col": "BETA",
            "beta_col": "BETA",
            "se_col": "SE",
            "imp_z_col": "Z",
        },
        policies=default_policies(),
    )

    assert result["POS"].to_list() == [1, 3]
    assert result["SE"].to_list() == pytest.approx([0.02, 0.1])
    assert qc["se_missing_on_entry"] == 2
    assert qc["se_recovered_from_z"] == 1
    assert qc["beta_z_sign_mismatches"] == 1
    assert qc["removed_beta_z_sign_mismatches"] == 1


def test_beta_z_sign_mismatch_fail_policy_stops_before_se_calculation():
    with pytest.raises(
        BetaZSignMismatchError,
        match=(
            "Chromosome 7: 2 variant.*No positive SE can satisfy "
            "Z = BETA / SE.*beta_z_sign_mismatch is 'fail'"
        ),
    ):
        derive_effect_and_standard_error_from_z(
            "7",
            pl.DataFrame({
                "BETA": [0.2, 0.0],
                "Z": [-2.0, 4.0],
            }),
            {
                "beta_or_col": "BETA",
                "beta_col": "BETA",
                "se_col": None,
                "imp_z_col": "Z",
            },
            policies=default_policies().with_overrides({
                "effect_from_z.beta_z_sign_mismatch": "fail",
            }),
        )


def test_z_only_effect_is_labelled_and_warned_as_standardized_approximation():
    logger = _CaptureLogger()
    policies = default_policies().with_overrides({
        "sample_size.trait_type": "binary",
    })

    result, qc, mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({"Z": [2.0], "EAF": [0.25], "Neff": [1000]}),
        {
            "imp_z_col": "Z",
            "eaf_col": "EAF",
            "beta_or_col": None,
            "se_col": None,
        },
        logger=logger,
        policies=policies,
    )

    expected_se = 1.0 / math.sqrt(2.0 * 0.25 * 0.75 * 1000.0)
    assert result["SE"][0] == pytest.approx(expected_se)
    assert result["BETA"][0] == pytest.approx(expected_se * 2.0)
    assert mapping["beta_col"] == "BETA"
    assert qc["effect_estimate_method"] == "metal_large_n"
    assert qc["effect_estimate_scale"] == "standardized_effect_estimate"
    assert "PMID 20616382" in qc["effect_estimate_citation"]
    assert "https://pubmed.ncbi.nlm.nih.gov/30038396/" in (
        qc["effect_estimate_pubmed_urls"]
    )
    assert qc["effect_estimate_trait_type"] == "binary"
    assert qc["effect_estimate_eaf_assumption"] == "same_analysed_samples_as_z"
    assert qc["effect_estimate_binary_interpretation"] == (
        "approximate_standardized_not_log_odds_or_liability_scale"
    )

    information = " ".join(logger.infos)
    warnings = " ".join(logger.warnings)
    assert "effect_from_z.method='metal_large_n'" in information
    assert "PMID 20616382" in information
    assert "https://pubmed.ncbi.nlm.nih.gov/20616382/" in information
    assert "same analysed samples" in warnings
    assert "external reference panel" in warnings
    assert "logistic-regression log odds ratio" in warnings
    assert "liability-scale effect" in warnings
    assert "29763751" not in (effect_from_z.__doc__ or "") + information


def test_z_only_reconstruction_stops_chromosome_x_by_default():
    with pytest.raises(
        XChromosomeZOnlyReconstructionError,
        match=(
            "cannot scientifically reconstruct chromosome-X BETA and SE.*"
            "No chromosome-X BETA or SE was calculated.*"
            "x_chromosome_z_only_action.*allow_autosomal_assumption"
        ),
    ):
        derive_effect_and_standard_error_from_z(
            "X",
            pl.DataFrame({"Z": [2.0], "EAF": [0.25], "Neff": [1000.0]}),
            {
                "imp_z_col": "Z",
                "eaf_col": "EAF",
                "beta_or_col": None,
                "se_col": None,
            },
            policies=default_policies(),
        )


def test_explicit_x_override_reproduces_autosomal_approximation_with_provenance():
    logger = _CaptureLogger()
    result, qc, _mapping = derive_effect_and_standard_error_from_z(
        "X",
        pl.DataFrame({"Z": [2.0], "EAF": [0.25], "Neff": [1000.0]}),
        {
            "imp_z_col": "Z",
            "eaf_col": "EAF",
            "beta_or_col": None,
            "se_col": None,
        },
        logger=logger,
        policies=default_policies().with_overrides({
            "effect_from_z.x_chromosome_z_only_action": (
                "allow_autosomal_assumption"
            ),
        }),
    )

    expected_se = 1.0 / math.sqrt(2.0 * 0.25 * 0.75 * 1000.0)
    assert result["SE"][0] == pytest.approx(expected_se)
    assert result["BETA"][0] == pytest.approx(2.0 * expected_se)
    assert qc["x_chromosome_z_only_policy"] == {
        "applies": True,
        "chromosome": "X",
        "action": "allow_autosomal_assumption",
        "assumption": "diploid_autosomal_genotype_variance_2p1mp",
    }
    assert "SCIENTIFIC OVERRIDE" in " ".join(logger.warnings)
    assert "not an X-specific reconstruction" in " ".join(logger.warnings)


def test_chromosome_x_exact_beta_from_supplied_se_and_z_is_still_allowed():
    result, qc, _mapping = derive_effect_and_standard_error_from_z(
        "X",
        pl.DataFrame({"Z": [2.0], "SE": [0.1]}),
        {
            "imp_z_col": "Z",
            "beta_or_col": None,
            "se_col": "SE",
        },
        policies=default_policies(),
    )

    assert result["BETA"][0] == pytest.approx(0.2)
    assert qc["status"] == "BETA derived from Z and SE"


def test_zhu_2016_method_uses_finite_sample_adjustment_and_pubmed_provenance():
    policies = default_policies().with_overrides({
        "effect_from_z.method": "zhu_2016",
    })

    result, qc, _mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({"Z": [2.0], "EAF": [0.25], "Neff": [1000.0]}),
        {
            "imp_z_col": "Z",
            "eaf_col": "EAF",
            "beta_or_col": None,
            "se_col": None,
        },
        policies=policies,
    )

    expected_se = 1.0 / math.sqrt(2.0 * 0.25 * 0.75 * (1000.0 + 4.0))
    assert result["SE"][0] == pytest.approx(expected_se)
    assert result["BETA"][0] == pytest.approx(2.0 * expected_se)
    assert qc["effect_estimate_method"] == "zhu_2016"
    assert qc["effect_estimate_citation"] == "Zhu et al. 2016 (PMID 27019110)"
    assert qc["effect_estimate_pubmed_urls"] == [
        "https://pubmed.ncbi.nlm.nih.gov/27019110/"
    ]


def test_configured_phenotype_standard_deviation_scales_beta_and_se_together():
    baseline, _baseline_qc, _ = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({"Z": [3.0], "EAF": [0.4], "Neff": [2000.0]}),
        {"imp_z_col": "Z", "eaf_col": "EAF"},
        policies=default_policies(),
    )
    scaled, qc, _ = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({"Z": [3.0], "EAF": [0.4], "Neff": [2000.0]}),
        {"imp_z_col": "Z", "eaf_col": "EAF"},
        policies=default_policies().with_overrides({
            "effect_from_z.phenotype_standard_deviation": 2.5,
            "sample_size.trait_type": "binary",
        }),
    )

    assert scaled["BETA"][0] == pytest.approx(2.5 * baseline["BETA"][0])
    assert scaled["SE"][0] == pytest.approx(2.5 * baseline["SE"][0])
    assert scaled["BETA"][0] / scaled["SE"][0] == pytest.approx(3.0)
    assert qc["effect_estimate_scale"] == "phenotype_scale_effect_estimate"
    assert qc["effect_estimate_scale_source"] == (
        "configured_phenotype_standard_deviation"
    )
    assert qc["effect_estimate_phenotype_standard_deviation"] == 2.5
    assert qc["effect_estimate_binary_interpretation"] == (
        "approximate_phenotype_scaled_not_log_odds_or_liability_scale"
    )


def test_reconstruction_fails_instead_of_exporting_numeric_overflow():
    policies = default_policies().with_overrides({
        "effect_from_z.phenotype_standard_deviation": 1e308,
    })

    with pytest.raises(
        ValueError,
        match="produced a non-finite BETA or non-positive/non-finite SE",
    ):
        derive_effect_and_standard_error_from_z(
            "1",
            pl.DataFrame({"Z": [1e308], "EAF": [0.25], "Neff": [1000.0]}),
            {"imp_z_col": "Z", "eaf_col": "EAF"},
            policies=policies,
        )


def test_z_only_reconstruction_rejects_endpoint_eaf_even_when_general_policy_keeps_it():
    policies = default_policies().with_overrides({"eaf.degenerate": "keep"})

    result, qc, _mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({
            "Z": [1.0, -1.0, 2.0],
            "EAF": [0.0, 1.0, 0.25],
            "Neff": [1000.0, 1000.0, 1000.0],
        }),
        {
            "imp_z_col": "Z",
            "eaf_col": "EAF",
            "beta_or_col": None,
            "se_col": None,
        },
        policies=policies,
    )

    assert result["EAF"].to_list() == [0.25]
    assert qc["removed_eaf_invalid"] == 2
    assert math.isfinite(result.item(0, "BETA"))
    assert math.isfinite(result.item(0, "SE"))


def test_z_only_reconstruction_rejects_invalid_z_even_when_general_policy_keeps_it():
    policies = default_policies().with_overrides({
        "validation.z_invalid": "keep",
    })

    result, qc, _mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({
            "Z": [float("nan"), 1.0],
            "EAF": [0.25, 0.25],
            "Neff": [1000.0, 1000.0],
        }),
        {
            "imp_z_col": "Z",
            "eaf_col": "EAF",
            "beta_or_col": None,
            "se_col": None,
        },
        policies=policies,
    )

    assert result["Z"].to_list() == [1.0]
    assert qc["removed_z_null"] == 1
    assert qc["removed_invalid_denominator"] == 0


def test_z_only_reconstruction_applies_configured_effective_variance_floor():
    result, qc, _mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({
            "Z": [1.0, 1.0],
            "EAF": [1e-8, 0.01],
            "Neff": [100_000.0, 100_000.0],
        }),
        {
            "imp_z_col": "Z",
            "eaf_col": "EAF",
            "beta_or_col": None,
            "se_col": None,
        },
        policies=default_policies(),
    )

    assert result["EAF"].to_list() == [0.01]
    assert qc["minimum_effective_variance"] == pytest.approx(1.0)
    assert qc["low_effective_variance_detected"] == 1
    assert qc["removed_low_effective_variance"] == 1
    assert qc["removed_invalid_denominator"] == 0


def test_effective_variance_floor_can_be_explicitly_disabled_without_disabling_denominator_guard():
    policies = default_policies().with_overrides({
        "effect_from_z.method": "zhu_2016",
        "effect_from_z.minimum_effective_variance": 0.0,
    })

    result, qc, _mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({
            "Z": [1.0, 1e308],
            "EAF": [1e-8, 0.25],
            "Neff": [100_000.0, 1000.0],
        }),
        {
            "imp_z_col": "Z",
            "eaf_col": "EAF",
            "beta_or_col": None,
            "se_col": None,
        },
        policies=policies,
    )

    assert result["EAF"].to_list() == [1e-8]
    assert qc["removed_low_effective_variance"] == 0
    assert qc["removed_invalid_denominator"] == 1
    assert result.item(0, "BETA") == pytest.approx(22.360568)
    assert math.isfinite(result.item(0, "BETA"))
    assert math.isfinite(result.item(0, "SE"))


def test_low_effective_variance_can_fail_with_actionable_message():
    policies = default_policies().with_overrides({
        "effect_from_z.low_effective_variance_action": "fail",
    })

    with pytest.raises(
        ValueError,
        match=r"2\*EAF\*\(1-EAF\)\*Neff below the configured minimum 1.0",
    ):
        derive_effect_and_standard_error_from_z(
            "1",
            pl.DataFrame({"Z": [1.0], "EAF": [1e-8], "Neff": [100_000.0]}),
            {
                "imp_z_col": "Z",
                "eaf_col": "EAF",
                "beta_or_col": None,
                "se_col": None,
            },
            policies=policies,
        )
