"""Contracts for valid null effects and returned harmonisation artifacts."""

from pathlib import Path

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.z_score import (
    derive_z_score_from_effect_and_standard_error,
)
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.service import (
    PipelineError,
    _resolved_harmonisation_outputs,
)


def test_zero_beta_is_kept_and_produces_zero_z_by_default(tmp_path):
    frame = pl.DataFrame({"BETA": [0.0, 0.1], "SE": [0.2, 0.2]})
    columns = {
        "beta_col": "BETA",
        "se_col": "SE",
        "imp_z_col": "NA",
        "gwas_outputname": "study",
        "output_folder": str(tmp_path),
    }

    result, qc, resolved = derive_z_score_from_effect_and_standard_error(
        "1", frame, columns, policies=default_policies()
    )

    assert default_policies().get("validation.beta_zero") == "keep"
    assert result.height == 2
    assert result[resolved["imp_z_col"]].to_list() == [0.0, 0.5]
    assert qc["variants_with_zero_beta"] == 1
    assert qc["beta_zero_action"] == "keep"
    assert qc["beta_zero_removed_flag"] is False


def test_configured_se_division_floor_rejects_only_numerically_unsafe_rows(tmp_path):
    frame = pl.DataFrame(
        {"BETA": [0.2, 0.2, 0.2], "SE": [0.1, 1.0e-13, None]}
    )
    columns = {
        "beta_col": "BETA",
        "se_col": "SE",
        "imp_z_col": None,
        "gwas_outputname": "study",
        "output_folder": str(tmp_path),
    }

    result, qc, resolved = derive_z_score_from_effect_and_standard_error(
        "1", frame, columns, policies=default_policies()
    )

    assert result.height == 1
    assert result[resolved["imp_z_col"]].to_list() == [2.0]
    assert qc["variants_removed_due_to_missing_se"] == 1
    assert qc["variants_removed_due_to_invalid_se"] == 1
    assert qc["variants_removed_invalid_beta_se"] == 2
    assert qc["se_division_floor"] == pytest.approx(1.0e-12)


def test_primary_output_contract_returns_exact_existing_files(tmp_path):
    sample = "study"
    harmonisation = load_configuration().modules.harmonisation
    expected = {
        "GRCh37": tmp_path / "study_GRCh37_merged.vcf.gz",
        "GRCh38": tmp_path / "study_GRCh38_merged.vcf.gz",
        "gwas2vcf": tmp_path / "study_gwas2vcf_GRCh37_merged.vcf.gz",
    }
    for path in expected.values():
        path.write_bytes(b"vcf\n")

    outputs = _resolved_harmonisation_outputs(
        tmp_path,
        sample,
        "GRCh37",
        harmonisation.output_layout.root,
        harmonisation.vcf_processing.model_dump(mode="python"),
    )
    assert outputs == {name: str(path) for name, path in expected.items()}


def test_primary_output_contract_rejects_missing_promised_file(tmp_path):
    harmonisation = load_configuration().modules.harmonisation
    with pytest.raises(PipelineError, match="cannot return all promised VCF outputs"):
        _resolved_harmonisation_outputs(
            tmp_path,
            "study",
            "GRCh37",
            harmonisation.output_layout.root,
            harmonisation.vcf_processing.model_dump(mode="python"),
        )
