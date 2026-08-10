"""Scientific provenance tests for standardized effects derived from Z."""

import polars as pl

from postgwas.modules.harmonisation import effect_from_z
from postgwas.modules.harmonisation.effect_from_z import (
    derive_effect_and_standard_error_from_z,
)
from postgwas.modules.harmonisation.policies import default_policies


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

    assert result["BETA"][0] == result["SE"][0] * 2.0
    assert mapping["beta_col"] == "BETA"
    assert qc["effect_estimate_scale"] == "standardized_effect_estimate"
    assert qc["effect_estimate_citation"] == "Zhu et al. 2016; PMID 27019110"
    assert qc["effect_estimate_trait_type"] == "binary"
    assert qc["effect_estimate_eaf_assumption"] == "same_analysed_samples_as_z"
    assert qc["effect_estimate_binary_interpretation"] == (
        "approximate_standardized_not_log_odds_or_liability_scale"
    )

    information = " ".join(logger.infos)
    warnings = " ".join(logger.warnings)
    assert "standardized BETA estimate" in information
    assert "PMID 27019110" in information
    assert "same analysed samples" in warnings
    assert "external reference panel" in warnings
    assert "logistic-regression log odds ratio" in warnings
    assert "liability-scale effect" in warnings
    assert "29763751" not in (effect_from_z.__doc__ or "") + information
