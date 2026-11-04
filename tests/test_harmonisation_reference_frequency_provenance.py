"""Reference AF headers describe complete-site attachment without export guards."""

import json

import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.vcf_processing import _render_provenance_headers
from postgwas.modules.harmonisation.vcf_provenance import build_harmonisation_vcf_provenance


def _provenance(tmp_path, source):
    replacement = tmp_path / "reference.tsv"
    sample = {
        "gwas_outputname": "study", "sumstat_file": str(tmp_path / "study.tsv"),
        "output_root": str(tmp_path / "output"), "beta_or_col": "BETA",
        "se_col": "SE", "imp_z_col": "Z", "pval_col": "P",
        "eaf_col": None if source == "external" else "EAF",
        "eaffile": str(replacement) if source == "external" else None,
        "eafcolumn": "EUR" if source == "external" else None,
        "info_source": "fixed_cli", "fixed_info": .99,
    }
    decisions = {"effect_type": "beta", "pvalue_type": "raw",
                 "frequency_type": "effect_allele_frequency"}
    if source == "replacement":
        sample["eaf_col"] = "MAF"
        decisions.update(
            eaf_is_maf=True, eaf_input_interpretation="maf",
            eaf_selected_source="replacement_reference",
            eaf_reference_validation={"status": "replacement_selected"},
            eaf_replacement={"file": str(replacement), "column": "EUR",
                             "resolved_files": {"1": str(replacement)}},
        )
    config = load_configuration().modules.harmonisation.vcf_processing.model_dump()
    values = build_harmonisation_vcf_provenance(
        sample=sample, study_decisions=decisions,
        resource_preflight={"require_default_eaf": False,
                            "resources_by_chromosome": {"1": {"user_eaf_file": str(replacement)}}},
        completed_chromosomes=["1"], policies=default_policies(),
        resource_directory=str(tmp_path), output_directory=str(tmp_path / "output"),
        run_manifest=str(tmp_path / "manifest.json"), vcf_config=config,
    )
    return values, config


@pytest.mark.parametrize("source", ["external", "replacement"])
def test_reference_frequency_formula_is_recorded_without_export_guard_claims(tmp_path, source):
    values, config = _provenance(tmp_path, source)

    formula = values["af_swapped_formula"]
    assert "sum(all distinct ALT AF)" in formula
    assert "blocked" not in formula.lower()
    assert "input AF" not in formula
    if source == "replacement":
        assert "original study MAF" in formula
        assert values["af_input_column"] == "MAF"
    headers = _render_provenance_headers(
        values, vcf_config=config, input_build="GRCh37", output_build="GRCh37",
        merge_group="grch37",
    )
    name = config["provenance"]["headers"]["af_swapped_formula"]
    matches = [line for line in headers if line.startswith("##%s=" % name)]
    assert len(matches) == 1
    assert json.loads(matches[0].split("=", 1)[1]) == formula


def test_internal_study_frequency_header_formula_is_unchanged(tmp_path):
    values, _ = _provenance(tmp_path, "study")

    assert values["af_swapped_formula"] == "AF = 1 - input AF for an allele-swapped match"
