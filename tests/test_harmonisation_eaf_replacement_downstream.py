"""Original MAF provenance and selected reference EAF stay distinct."""

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation import study_properties
from postgwas.modules.harmonisation.concordance import service
from postgwas.modules.harmonisation.concordance.analysis import compare_input_to_vcf
from postgwas.modules.harmonisation.html_report import render_dataset_report
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.report_steps import summarise_chromosome_evidence
from postgwas.modules.harmonisation.vcf_provenance import build_harmonisation_vcf_provenance


def decision(tmp_path):
    return {
        "eaf_is_maf": True, "eaf_is_maf_initial": True,
        "eaf_is_maf_source": "dataset_reference_concordance",
        "eaf_input_interpretation": "maf",
        "eaf_selected_source": "replacement_reference",
        "eaf_replacement": {"file": str(tmp_path / "replacement_{chromosome}.tsv"),
                            "column": "PANEL_AF",
                            "resolved_files": {"1": str(tmp_path / "replacement_1.tsv")}},
        "eaf_reference_validation": {"status": "replacement_selected"},
    }


def summaries():
    return {"1": {"status": "ok", "stage_qc": {"eaf_qc": {
        "eaf_provenance": "reference_imputed", "maf_reference_decision": "maf"}}}}


def row():
    return SimpleNamespace(dataset_id="study", input_file=Path("input.tsv"),
        delimiter="tab", chromosome_column="CHR", position_column="POS",
        chromosome_position_column=None, effect_allele_column="EA",
        other_allele_column="OA", effect_allele_frequency_column="MAF",
        effect_column="BETA", standard_error_column="SE", z_score_column="Z",
        p_value_column="P", effect_type="beta", p_value_type="raw",
        external_eaf_file=None, external_eaf_column=None)


def test_selected_replacement_finalises_output_eaf_and_preserves_original_maf(tmp_path):
    original = decision(tmp_path)
    result = study_properties.finalise_eaf_decision_from_chromosomes(original, summaries(), ["1"])
    assert result["eaf_is_maf"] is False
    assert result["frequency_type"] == "effect_allele_frequency"
    assert result["eaf_input_interpretation"] == "maf"
    assert result["eaf_is_maf_initial"] is True
    assert result["eaf_selected_source"] == "replacement_reference"
    assert result["eaf_provenance"] == "reference_imputed"
    assert original["eaf_is_maf"] is True


def test_selected_replacement_without_completed_outputs_remains_unresolved(tmp_path):
    result = study_properties.finalise_eaf_decision_from_chromosomes(decision(tmp_path), {}, [])
    assert result["eaf_is_maf"] is None
    assert result["eaf_input_interpretation"] == "maf"


def test_selected_replacement_requires_complete_descriptor(tmp_path):
    selected = decision(tmp_path)
    selected["eaf_replacement"].pop("column")
    with pytest.raises(study_properties.StudyPropertyError, match="replacement.*column"):
        study_properties.finalise_eaf_decision_from_chromosomes(selected, summaries(), ["1"])


@pytest.mark.parametrize("invalid_decision", [
    {"eaf_input_interpretation":"eaf"},
    {"eaf_reference_validation":{"status":"failed"}},
])
def test_selected_replacement_requires_recorded_maf_support(tmp_path, invalid_decision):
    selected = decision(tmp_path)
    selected.update(invalid_decision)
    with pytest.raises(study_properties.StudyPropertyError, match="replacement.*MAF"):
        study_properties.finalise_eaf_decision_from_chromosomes(selected, summaries(), ["1"])


def test_vcf_provenance_records_original_maf_and_original_replacement_paths(tmp_path):
    selected = decision(tmp_path)
    selected.update(effect_type="beta", pvalue_type="raw", frequency_type="effect_allele_frequency")
    values = build_harmonisation_vcf_provenance(
        sample={"gwas_outputname":"study", "sumstat_file":str(tmp_path/"input.tsv"),
                "output_root":str(tmp_path), "eaf_col":"MAF", "eaffile":None,
                "beta_or_col":"BETA", "pval_col":"P", "info_source":"fixed_cli", "fixed_info":.9},
        study_decisions=selected,
        resource_preflight={"require_default_eaf":False, "resources_by_chromosome":{"1":{
            "user_eaf_file":str(tmp_path/"staged.parquet")}}},
        completed_chromosomes=["1"], policies=default_policies(),
        resource_directory=str(tmp_path), output_directory=str(tmp_path),
        run_manifest=str(tmp_path/"manifest.json"),
        vcf_config=load_configuration().modules.harmonisation.vcf_processing.model_dump())
    assert values["af_input_column"] == "MAF"
    assert values["af_input_type"] == "minor_allele_frequency"
    assert "reference proxy" in values["af_source"]
    assert values["af_reference_template"] == selected["eaf_replacement"]["file"]
    assert values["af_reference_files"] == selected["eaf_replacement"]["resolved_files"]["1"]
    assert "PANEL_AF" in values["af_harmonisation"]
    assert "input AF" not in values["af_swapped_formula"]


def test_concordance_selected_row_preserves_original_mapping_and_uses_reference(tmp_path):
    original = row()
    selected = service._selected_eaf_row(original, decision(tmp_path))
    assert selected is not original
    assert selected.effect_allele_frequency_column == "MAF"
    assert selected.external_eaf_file == Path(decision(tmp_path)["eaf_replacement"]["file"])
    assert selected.external_eaf_column == "PANEL_AF"
    assert original.external_eaf_file is None
    config = load_configuration().modules.harmonisation
    source = pl.DataFrame({"CHR":["1"], "POS":[10], "EA":["A"], "OA":["G"],
                           "MAF":[.2], "BETA":[.2], "SE":[.1], "Z":[2.], "P":[.1]})
    reference = pl.DataFrame({"CHROM":["1"], "POS":[10], "REF":["G"], "ALT":["A"], "PANEL_AF":[.8]})
    vcf = pl.DataFrame({"CHROM":["1"], "POS":[10], "ID":["v1"], "REF":["G"], "ALT":["A"],
                        "ES":[.2], "SE":[.1], "EZ":[2.], "AF":[.8], "LP":[1.]})
    for observed, expected_pass in ((.8, True), (.2, False)):
        result = compare_input_to_vcf(source, vcf.with_columns(pl.lit(observed).alias("AF")), selected,
            effect_type="beta", p_value_type="raw", eaf_is_maf=False, input_eaf_is_maf=True,
            settings=config.concordance_validation, policies=default_policies(),
            external_eaf_frame=reference, external_eaf_mapping=config.external_eaf_mapping)
        assert (result.metric_summary["allele_frequency"]["concordant"] == 1) is expected_pass
        assert result.matched.item(0, "input_af_original") == .2
        assert result.matched.item(0, "input_frequency_is_maf") is True
    assert source["MAF"].to_list() == [.2]


@pytest.mark.parametrize("observed,concordant", [(.8, 1), (.2, 0)])
def test_selected_replacement_concordance_checks_reversed_reference_alleles(tmp_path, observed, concordant):
    config = load_configuration().modules.harmonisation
    selected = service._selected_eaf_row(row(), decision(tmp_path))
    source = pl.DataFrame({"CHR":["1"], "POS":[10], "EA":["A"], "OA":["G"],
                           "MAF":[.2], "BETA":[.2], "SE":[.1], "Z":[2.], "P":[.1]})
    reference = pl.DataFrame({"CHROM":["1"], "POS":[10], "REF":["A"], "ALT":["G"], "PANEL_AF":[.2]})
    vcf = pl.DataFrame({"CHROM":["1"], "POS":[10], "ID":["v1"], "REF":["G"], "ALT":["A"],
                        "ES":[.2], "SE":[.1], "EZ":[2.], "AF":[observed], "LP":[1.]})
    result = compare_input_to_vcf(source, vcf, selected,
        effect_type="beta", p_value_type="raw", eaf_is_maf=False, input_eaf_is_maf=True,
        settings=config.concordance_validation, policies=default_policies(),
        external_eaf_frame=reference, external_eaf_mapping=config.external_eaf_mapping)
    metric = result.metric_summary["allele_frequency"]
    assert metric["checked"] == 1
    assert metric["concordant"] == concordant
    assert result.matched.item(0, "expected_allele_frequency") == pytest.approx(.8)
    assert result.matched.item(0, "input_af_original") == .2
    assert result.matched.item(0, "input_af") == .2
    assert result.matched.item(0, "reference_expected_alt_af") == pytest.approx(.8)
    assert result.status == ("PASS" if concordant else "FAIL")


@pytest.mark.parametrize("observed,expected_status", [(0.7, "PASS"), (0.8, "FAIL")])
def test_selected_replacement_concordance_checks_complete_multiallelic_sum(tmp_path, observed, expected_status):
    config = load_configuration().modules.harmonisation
    selected = service._selected_eaf_row(row(), decision(tmp_path))
    source = pl.DataFrame({"CHR":["1"], "POS":[10], "EA":["A"], "OA":["G"],
                           "MAF":[.2], "BETA":[.2], "SE":[.1], "Z":[2.], "P":[.1]})
    reference = pl.DataFrame({"CHROM":["1","1"], "POS":[10,10], "REF":["A","A"],
                             "ALT":["G","C"], "PANEL_AF":[.2,.1]})
    vcf = pl.DataFrame({"CHROM":["1"], "POS":[10], "ID":["v1"], "REF":["G"], "ALT":["A"],
                        "ES":[.2], "SE":[.1], "EZ":[2.], "AF":[observed], "LP":[1.]})
    result = compare_input_to_vcf(source, vcf, selected,
        effect_type="beta", p_value_type="raw", eaf_is_maf=False, input_eaf_is_maf=True,
        settings=config.concordance_validation, policies=default_policies(),
        external_eaf_frame=reference, external_eaf_mapping=config.external_eaf_mapping)
    assert result.status == expected_status
    assert result.matched.item(0, "reference_expected_alt_af") == pytest.approx(0.7)
    assert result.matched.item(0, "input_af_original") == 0.2


@pytest.mark.parametrize("original_is_maf", [None, False])
def test_ordinary_external_concordance_retains_direct_reference_contract(tmp_path, original_is_maf):
    config = load_configuration().modules.harmonisation
    ordinary = row()
    ordinary.effect_allele_frequency_column = None
    ordinary.external_eaf_file = tmp_path / "external.tsv"
    ordinary.external_eaf_column = "PANEL_AF"
    source = pl.DataFrame({"CHR":["1"], "POS":[10], "EA":["A"], "OA":["G"],
                           "BETA":[.2], "SE":[.1], "Z":[2.], "P":[.1]})
    reference = pl.DataFrame({"CHROM":["1"], "POS":[10], "REF":["A"], "ALT":["G"], "PANEL_AF":[.2]})
    vcf = pl.DataFrame({"CHROM":["1"], "POS":[10], "ID":["v1"], "REF":["G"], "ALT":["A"],
                        "ES":[.2], "SE":[.1], "EZ":[2.], "AF":[.8], "LP":[1.]})
    result = compare_input_to_vcf(source, vcf, ordinary,
        effect_type="beta", p_value_type="raw", eaf_is_maf=False, input_eaf_is_maf=original_is_maf,
        settings=config.concordance_validation, policies=default_policies(),
        external_eaf_frame=reference, external_eaf_mapping=config.external_eaf_mapping)
    assert result.metric_summary["allele_frequency"]["checked"] == 0


def test_selected_replacement_direct_multiallelic_alt_preserves_original_maf(tmp_path):
    config = load_configuration().modules.harmonisation
    selected = service._selected_eaf_row(row(), decision(tmp_path))
    source = pl.DataFrame({"CHR":["1"], "POS":[10], "EA":["G"], "OA":["A"],
                           "MAF":[.1], "BETA":[.2], "SE":[.1], "Z":[2.], "P":[.1]})
    reference = pl.DataFrame({"CHROM":["1","1"], "POS":[10,10], "REF":["G","G"],
                             "ALT":["A","C"], "PANEL_AF":[.8,.1]})
    vcf = pl.DataFrame({"CHROM":["1"], "POS":[10], "ID":["v1"], "REF":["G"], "ALT":["A"],
                        "ES":[-.2], "SE":[.1], "EZ":[-2.], "AF":[.8], "LP":[1.]})
    result = compare_input_to_vcf(source, vcf, selected,
        effect_type="beta", p_value_type="raw", eaf_is_maf=False, input_eaf_is_maf=True,
        settings=config.concordance_validation, policies=default_policies(),
        external_eaf_frame=reference, external_eaf_mapping=config.external_eaf_mapping)
    assert result.status == "PASS"
    assert result.metric_summary["allele_frequency"]["checked"] == 1
    assert result.metric_summary["allele_frequency"]["concordant"] == 1
    assert result.matched.item(0, "input_af") == .1
    assert result.matched.item(0, "input_af_original") == .1
    assert result.matched.item(0, "reference_expected_alt_af") == .8
    assert result.matched.item(0, "expected_allele_frequency") == .8


def test_replacement_position_only_diagnostics_do_not_compare_original_maf(tmp_path):
    config = load_configuration().modules.harmonisation
    selected = service._selected_eaf_row(row(), decision(tmp_path))
    source = pl.DataFrame({"CHR":["1"], "POS":[10], "EA":["A"], "OA":["G"],
                          "MAF":[.2], "BETA":[.2], "SE":[.1], "Z":[2.], "P":[.1]})
    reference = pl.DataFrame({"CHROM":["1"], "POS":[10], "REF":["G"], "ALT":["A"], "PANEL_AF":[.8]})
    vcf = pl.DataFrame({"CHROM":["1"], "POS":[10], "ID":["v1"], "REF":["T"], "ALT":["G"],
                       "ES":[.2], "SE":[.1], "EZ":[2.], "AF":[.8], "LP":[1.]})
    result = compare_input_to_vcf(source, vcf, selected, effect_type="beta", p_value_type="raw",
        eaf_is_maf=False, input_eaf_is_maf=True, settings=config.concordance_validation,
        policies=default_policies(), external_eaf_frame=reference,
        external_eaf_mapping=config.external_eaf_mapping)
    assert result.position_matches.height == 1
    assert result.position_matches.item(0, "allele_frequency_checked") is False
    assert result.position_matches.item(0, "expected_allele_frequency") is None
    assert result.position_matches.item(0, "effect_checked") is True


@pytest.mark.parametrize("template", [False, True])
def test_selected_replacement_stages_with_original_internal_column(tmp_path, template):
    selected = decision(tmp_path)
    if not template:
        selected["eaf_replacement"]["file"] = str(tmp_path / "replacement_1.tsv")
    path = tmp_path / "replacement_1.tsv"
    path.write_text("CHROM\tPOS\tREF\tALT\tPANEL_AF\n1\t10\tG\tA\t0.8\n")
    source_row = service._selected_eaf_row(row(), selected)
    messages = []
    log = SimpleNamespace(info=lambda message, *args, **kwargs: messages.append(message))
    config = load_configuration().modules.harmonisation
    staged = service._stage_external_eaf(source_row, config.external_eaf_mapping,
        tmp_path / "staged.parquet", default_policies(), log, "GRCh37",
        batch_rows=1, compression="zstd")
    if template:
        assert staged is None
        assert "per concordance chromosome" in " ".join(messages)
    else:
        assert pl.read_parquet(staged)["PANEL_AF"].cast(pl.Float64).to_list() == [.8]
    assert source_row.effect_allele_frequency_column == "MAF"


def test_report_labels_original_maf_and_selected_reference_proxy(tmp_path):
    selected = decision(tmp_path)
    selected.update(frequency_type="effect_allele_frequency", eaf_provenance="reference_imputed")
    report = render_dataset_report({"dataset_id":"study", "status":"OK"},
        {"dataset":{"study_decisions":selected}, "config":{"eaf_col":"MAF"}})
    assert "Original input frequency" in report
    assert "reference proxy" in report
    assert "PANEL_AF" in report
    assert str(tmp_path / "replacement_{chromosome}.tsv") in report


def test_reports_preserve_reference_evidence_denominator_and_actual_replacement_count(tmp_path):
    selected = decision(tmp_path)
    selected["eaf_evidence"] = {"usable":190, "unusable":10}
    selected["eaf_reference_validation"].update(
        reference_effect_allele_minor_fraction=.8, comparable_variants=100)
    manifest = {"dataset":{"study_decisions":selected, "completed":["1"],
        "chromosome_summaries":{"1":{"stage_qc":{"eaf_qc":{
            "eaf_replacement_rows":3, "eaf_replacement_source_column":"PANEL_AF",
            "eaf_replacement_source_file":str(tmp_path/"replacement_1.tsv")}}}}},
        "config":{"eaf_col":"MAF"}}
    evidence = summarise_chromosome_evidence(manifest)["frequency"]
    assert evidence["counts"]["eaf_replacement_rows"]["value"] == 3
    assert evidence["counts"]["eaf_replacement_rows"]["complete"] is True
    report = render_dataset_report({"dataset_id":"study", "status":"OK"}, manifest)
    assert "Reference EAF values retained" in report
    assert "Comparable variants" in report and ">100<" in report
    assert "Detailed evidence" in report
    assert "Original input frequencies excluded from evidence" in report
    assert "Missing selected output AF" in report
