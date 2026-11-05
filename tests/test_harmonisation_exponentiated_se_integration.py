"""Encoded SE decisions reach reporting, provenance and independent concordance."""

import math

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.concordance.analysis import prepare_input_table
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.study_properties import resolve_study_properties
from postgwas.modules.harmonisation.study_properties import StudyPropertyError
from postgwas.modules.harmonisation.study_reporting import _se_fields
from postgwas.modules.harmonisation.summary_statistics_io import read_summary_statistics
from postgwas.modules.harmonisation.vcf_provenance import _se_provenance

from test_harmonisation_concordance import _row, _settings


def _source():
    # Full-precision source examples; the expected SE values are independently
    # hand-checked inverses, rather than reconstructed by the function tested.
    return pl.DataFrame({
        "CHR": ["1", "1"], "BP": [758443, 794707],
        "EA": ["C", "C"], "OA": ["G", "T"],
        "BETA": [0.893986328923368, 0.855499389674865],
        "SE": [1.13037330491772, 1.18700177512822],
        "P": [0.360476858647774, 0.362613225398603],
        "Z": [None, None], "EAF": [0.135199240986717, 0.0607210626185958],
    })


def _prepare(frame):
    return prepare_input_table(
        frame, _row(effect_type="odds_ratio"), effect_type="odds_ratio",
        se_scale="exponentiated_log_odds", p_value_type="raw", eaf_is_maf=False,
        settings=_settings(), policies=default_policies(),
    )


def test_new_se_scale_is_schema_validated_and_study_decision_is_recorded():
    config = load_configuration(cli_overrides={
        "modules.harmonisation.policies.effect.se_scale": "exponentiated_log_odds",
    })
    policies = default_policies().with_overrides({
        "effect.type": "odds_ratio", "effect.se_scale": "exponentiated_log_odds",
        "pvalue.type": "raw",
    })
    decisions = resolve_study_properties(
        _source(), {"beta_or_col": "BETA", "se_col": "SE", "pval_col": "P"},
        policies=policies,
    )
    assert config.modules.harmonisation.policies["effect"]["se_scale"] == "exponentiated_log_odds"
    assert default_policies().get("effect.se_scale") == "auto"
    assert decisions["se_scale"] == "exponentiated_log_odds"
    assert decisions["se_scale_source"] == "policy"


def test_concordance_decodes_source_se_and_recovers_matching_z_without_mutating_input():
    source = _source()
    original = source.clone()
    result = _prepare(source)
    assert source.equals(original)
    assert result["input_se"].to_list() == pytest.approx([
        0.1225479365022153, 0.1714306111009031,
    ])
    assert result["input_beta"][0] == pytest.approx(-0.11206479595757401)
    reconstructed_p = math.erfc(abs(result["input_z"][0]) / math.sqrt(2))
    assert reconstructed_p == pytest.approx(0.360476858647774, abs=1e-12)
    assert result["input_se_source"].to_list() == ["supplied_canonical_scale"] * 2
    assert result["valid_input_variant"].to_list() == [True, True]


@pytest.mark.parametrize("invalid", [1.0, 0.9, 0.0, -1.0, float("inf"), float("nan")])
def test_concordance_rejects_invalid_encoded_values_instead_of_recovering_them(invalid):
    source = _source().with_columns(pl.Series("SE", [invalid, 1.18700177512822]))
    with pytest.raises(ValueError, match="exponentiated|encoded|greater than 1"):
        _prepare(source)


def test_concordance_keeps_existing_recovery_for_missing_encoded_se():
    source = _source().with_columns(pl.Series("SE", [None, 1.18700177512822]))
    result = _prepare(source)
    assert result["input_se"][0] == pytest.approx(0.1225479365022153, abs=1e-12)
    assert result["input_se_source"][0] == "calculated_abs_beta_div_z"


def test_vcf_provenance_describes_encoding_and_preserves_source_column_name():
    values = _se_provenance(
        {"beta_or_col": "OR", "se_col": "standard_error"},
        {"effect_type": "odds_ratio", "se_scale": "exponentiated_log_odds"},
        reconstruction={}, pvalue_tail=2, derive_partial_missing_se=False,
        missing="not available",
    )
    assert values["se_input_column"] == "standard_error"
    assert values["se_formula"] == "SE_logOR = ln(stored_SE)"
    assert "exponentiated" in values["se_input_scale"]


def test_study_report_shows_third_candidate_and_decoding_action():
    fields = _se_fields({
        "effect_type": "odds_ratio", "se_column": "standard_error",
        "se_scale": "exponentiated_log_odds", "se_scale_source": "policy",
        "se_scale_evidence": {
            "informative_variants": 120, "log_odds_matches": 0,
            "raw_odds_ratio_matches": 0, "exponentiated_log_odds_matches": 120,
        },
    }, default_policies())
    by_label = {label: value for _, label, value in fields}
    assert "120 / 120" in by_label["Exponentiated-SE agreement"]
    assert "ln" in by_label["Planned action"]


@pytest.mark.parametrize("scale", ["auto", "exponentiated_log_odds"])
@pytest.mark.parametrize("first_se", ["bad", None])
def test_reader_distinguishes_unparseable_and_missing_encoded_se(tmp_path, scale, first_se):
    ratios = [0.5, 2.0] * 60
    z = [math.log(value) / 0.15 for value in ratios]
    source = pl.DataFrame({
        "CHR": ["1"] * 120, "BP": list(range(1, 121)),
        "EA": ["A"] * 120, "OA": ["G"] * 120, "OR": ratios,
        "SE": [first_se] + [str(math.exp(0.15))] * 119,
        "Z": z, "P": [math.erfc(abs(value) / math.sqrt(2)) for value in z],
    })
    path = tmp_path / "encoded.tsv"
    source.write_csv(path, separator="\t")
    policies = default_policies().with_overrides({
        "effect.type": "odds_ratio", "effect.se_scale": scale, "pvalue.type": "raw",
    })
    read = read_summary_statistics(
        str(path), str(tmp_path / "read"), {
            "gwas_outputname": "encoded", "chr_col": "CHR", "pos_col": "BP",
            "ea_col": "EA", "oa_col": "OA", "beta_or_col": "OR",
            "se_col": "SE", "imp_z_col": "Z", "pval_col": "P",
        }, dict(load_configuration().modules.harmonisation.output_layout.root),
        "\t", policies=policies, input_line_count=120,
    )
    frame, snapshot, _, _, mapping = read[:5]
    assert snapshot["SE"][0] == first_se
    assert frame["SE"][0] is None
    if first_se is None:
        decisions = resolve_study_properties(frame, mapping, policies=policies)
        assert decisions["se_scale"] == "exponentiated_log_odds"
    else:
        with pytest.raises(StudyPropertyError, match="SE.*unparseable|unparseable.*SE"):
            resolve_study_properties(frame, mapping, policies=policies)
