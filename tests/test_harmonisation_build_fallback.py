"""Real reference joins exercise the configurable dominant-build coverage exception."""

import json
from copy import deepcopy
from types import SimpleNamespace

import polars as pl
import pytest
import yaml

from postgwas.config import load_configuration
from postgwas.config.exporter import render_module_configuration
from postgwas.config.loader import load_run_configuration_for_module
from postgwas.core.completion import configuration_digest
from postgwas.core.direct_execution import _checkpoint_configuration
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.modules.harmonisation.genome_build import GenomeBuildError, infer_genome_build
from postgwas.modules.harmonisation.policies import PolicyError, default_policies
from postgwas.modules.harmonisation.study_reporting import render_study_decisions


def row(position, chromosome="1", ea="A", oa="G"):
    """Convert a zero-based fixture index to a valid one-based coordinate."""
    return (chromosome, position + 1, ea, oa)


def infer(tmp_path, study, first, second, *, overrides=None, enabled=True, logger=None, extra=None):
    paths = {}
    for name, rows in {"GRCh37": first, "GRCh38": second, **(extra or {})}.items():
        path = tmp_path / (name + ".csv")
        path.write_text("C,P,R,A\n" + "".join("%s,%s,%s,%s\n" % item for item in rows))
        paths[name] = str(path)
    settings = {"build.min_match_count": 1}
    if enabled is not None:
        settings.update({"build.low_coverage_fallback_enabled": enabled,
                         "build.low_coverage_min_unique_loci": 1})
    settings.update(overrides or {})
    frame = pl.DataFrame(study, schema={"chr": pl.String, "pos": pl.Int64,
                                      "ea": pl.String, "oa": pl.String}, orient="row")
    result = infer_genome_build(frame, paths,
        {"chr": "C", "pos": "P", "a1": "R", "a2": "A", "delimiter": "comma"},
        {"chr_col": "chr", "pos_col": "pos", "ea_col": "ea", "oa_col": "oa"},
        policies=default_policies().with_overrides(settings), logger=logger)
    assert frame.rows() == study  # Inference must not filter or mutate the study.
    return result


def padded(rows, count=300):
    return rows + [row(100000 + index) for index in range(count)]


@pytest.mark.parametrize("winner", ["GRCh37", "GRCh38"])
def test_low_coverage_accepts_either_build_and_can_be_explicitly_disabled(tmp_path, winner):
    study = [row(i) for i in range(20)]
    refs = [padded(study), [row(99999)]]
    if winner == "GRCh38":
        refs.reverse()
    disabled = infer(tmp_path, study, *refs, enabled=False)
    assert disabled["inferred_build"] == "Ambiguous"
    assert disabled["low_coverage_fallback"]["status"] == "disabled"
    accepted = infer(tmp_path, study, *refs)
    assert accepted["inferred_build"] == winner
    assert accepted["acceptance_path"] == "low_coverage_fallback"
    assert accepted["normal_failed_gates"] == ["build.min_reference_match_fraction"]
    assert accepted["distinct_loci"]["exclusive_snp_loci"][winner] == 20
    assert accepted["strand_evidence"][winner]["forward"] == 20


@pytest.mark.parametrize("winning,competing,accepted,gate", [
    (95, 5, True, None), (94, 5, False, "build.low_coverage_min_match_fraction"),
    (95, 6, False, "build.low_coverage_max_competitor_fraction"),
])
def test_same_testable_denominator_and_exact_fraction_boundaries(tmp_path, winning, competing, accepted, gate):
    study = [row(i) for i in range(100)]
    first = padded(study[:winning], 2000)
    # All coordinates are testable, including incompatible allele records.
    second = study[-competing:] + [row(i, ea="A", oa="C") for i in range(100)]
    result = infer(tmp_path, study, first, second)
    assert result["testable_variants"] == 100
    assert result["matches"] == {"GRCh37": winning, "GRCh38": competing}
    assert result["low_coverage_fallback"]["accepted"] is accepted
    if gate:
        assert gate in result["low_coverage_fallback"]["failed_gates"]
        assert gate in result["ambiguous_reason"]


@pytest.mark.parametrize("count,accepted", [(9999, False), (10000, True)])
def test_actual_default_unique_winner_only_snp_floor(tmp_path, count, accepted):
    study = [row(i) for i in range(count)]
    result = infer(tmp_path, study, padded(study, 100001), [row(999999)],
                   enabled=None)
    assert result["low_coverage_fallback"]["status"] == ("accepted" if accepted else "rejected")
    assert result["low_coverage_fallback"]["accepted"] is accepted
    assert result["distinct_loci"]["exclusive_snp_loci"]["GRCh37"] == count


@pytest.mark.parametrize("deduplicate", [True, False])
def test_duplicates_cannot_inflate_unique_evidence_but_row_fractions_are_preserved(tmp_path, deduplicate):
    study = [row(1)] * 20 + [row(2)]
    result = infer(tmp_path, study, padded([row(1)] * 4 + [row(2)]), [row(9999)],
        overrides={"build.low_coverage_min_unique_loci": 3,
                   "build.deduplicate_reference": deduplicate})
    assert result["matches"]["GRCh37"] == 21
    assert result["distinct_loci"]["exclusive_snp_loci"]["GRCh37"] == 2
    assert result["low_coverage_fallback"]["accepted"] is False


def test_both_compatible_and_opposing_alleles_at_same_locus_are_not_exclusive(tmp_path):
    study = [row(i) for i in range(100)] + [row(0, ea="A", oa="C")]
    result = infer(tmp_path, study, padded(study[:100], 2000),
                   [row(0, ea="A", oa="C"), row(1)],
                   overrides={"build.low_coverage_min_unique_loci": 99})
    assert result["matches"] == {"GRCh37": 100, "GRCh38": 2}
    assert result["distinct_loci"]["compatible_snp_loci"]["GRCh37"] == 100
    assert result["distinct_loci"]["exclusive_snp_loci"]["GRCh37"] == 98
    assert result["distinct_loci"]["exclusive_snp_loci"]["GRCh38"] == 0
    assert not result["low_coverage_fallback"]["accepted"]


def test_non_snp_identical_and_unknown_alleles_do_not_supply_unique_snp_floor(tmp_path):
    study = [row(1), row(2, ea="AA", oa="A"), row(3, ea="A", oa="A"), row(4, ea="N", oa="G")]
    result = infer(tmp_path, study, padded(study), [row(9999)],
                   overrides={"build.low_coverage_min_unique_loci": 2})
    assert result["matches"]["GRCh37"] == 4
    assert result["distinct_loci"]["exclusive_snp_loci"]["GRCh37"] == 1
    assert not result["low_coverage_fallback"]["accepted"]


@pytest.mark.parametrize("count,total,veto", [(99, 100, False), (100, 105, True), (100, 106, False)])
def test_chromosome_veto_requires_both_count_and_fraction_boundaries(tmp_path, count, total, veto):
    main = [row(i) for i in range(3000)]
    minor = [row(i, "2") for i in range(total)]
    first = padded(main + [row(99999, "2")], 40000)
    second = [row(99999)] + minor[:count] + [row(i, "2", ea="A", oa="C") for i in range(total)]
    result = infer(tmp_path, main + minor + [row(0, "3")], first, second)
    assert result["low_coverage_fallback"]["accepted"] is (not veto)
    assert result["per_chromosome"]["2"]["competitor_veto"] is veto
    assert result["per_chromosome"]["3"]["assessment"] == "unassessed"
    if veto:
        assert result["low_coverage_fallback"]["veto_chromosomes"] == ["2"]


def test_exact_chromosome_veto_fraction_and_normal_acceptance_scope(tmp_path):
    main = [row(i) for i in range(4000)]
    minor = [row(i, "2") for i in range(2000, 2100)]
    # 95 / 100 is exact, with a configurable 95-locus floor for this case.
    second = [row(99999)] + minor[:95] + [row(i, "2", ea="A", oa="C") for i in range(2000, 2100)]
    result = infer(tmp_path, main + minor, main + [row(99999, "2")], second,
                   overrides={"build.low_coverage_chromosome_min_competitor_loci": 95})
    assert result["acceptance_path"] == "normal"
    assert result["per_chromosome"]["2"]["competitor_veto"] is True
    assert result["low_coverage_fallback"]["status"] == "not_needed"


@pytest.mark.parametrize("overrides,gate", [
    ({"build.min_match_count": 101}, "build.min_match_count"),
    ({"build.min_match_fraction": 1.0}, "build.min_match_fraction"),
    ({"build.confidence_ratio": 1.0}, "build.confidence_ratio"),
])
def test_only_reference_coverage_may_be_waived(tmp_path, overrides, gate):
    study = [row(i) for i in range(100)]
    result = infer(tmp_path, study, padded(study[:95], 2000), study[95:], overrides=overrides)
    assert result["inferred_build"] == "Ambiguous"
    assert gate in result["normal_failed_gates"]
    assert result["low_coverage_fallback"]["status"] == "not_applicable"
    assert "sole failed" in result["ambiguous_reason"]


@pytest.mark.parametrize("study", [[], [row(123)], [row(1, ea="A", oa="C")]])
def test_empty_untestable_or_allele_incompatible_evidence_never_accepts(tmp_path, study):
    result = infer(tmp_path, study, [row(1)], [row(2)])
    assert result["inferred_build"] == "Ambiguous"
    assert not result["low_coverage_fallback"]["accepted"]


def test_declared_build_bypasses_fallback_and_two_reference_constraint(tmp_path):
    result = infer(tmp_path, [row(1)], padded([row(1)]), [row(2)],
                   overrides={"build.mode": "GRCh37"}, extra={"Third": [row(3)]})
    assert result["acceptance_path"] == "declared"
    assert result["forced"] is True
    assert result["low_coverage_fallback"]["status"] == "not_applicable"


def test_enabled_automatic_fallback_rejects_more_than_two_references(tmp_path):
    with pytest.raises(GenomeBuildError, match="exactly two.*low_coverage_fallback_enabled"):
        infer(tmp_path, [row(1)], [row(1)], [row(2)], enabled=None,
              extra={"Third": [row(3)]})
    assert infer(tmp_path, [row(1)], [row(1)], [row(2)], enabled=False,
                 extra={"Third": [row(3)]})["inferred_build"] == "GRCh37"


@pytest.mark.parametrize("key,invalid", [
    ("low_coverage_fallback_enabled", "sometimes"),
    ("low_coverage_min_match_fraction", -0.1),
    ("low_coverage_max_competitor_fraction", 1.01),
    ("low_coverage_min_unique_loci", 0),
    ("low_coverage_chromosome_min_competitor_loci", 1.5),
    ("low_coverage_chromosome_min_competitor_fraction", float("nan")),
])
def test_invalid_fallback_policies_fail_configuration_validation(key, invalid):
    with pytest.raises(PolicyError, match=key):
        default_policies().with_overrides({"build." + key: invalid})


def test_configuration_export_run_config_and_checkpoint_identity(tmp_path):
    exported = yaml.safe_load(render_module_configuration("harmonisation", style="values"))
    assert exported["policies"]["build"]["low_coverage_fallback_enabled"] is True
    assert exported["policies"]["build"]["low_coverage_min_unique_loci"] == 10000
    default = load_configuration()
    changed = load_configuration(cli_overrides={
        "modules.harmonisation.policies.build.low_coverage_fallback_enabled": False})
    assert configuration_digest(_checkpoint_configuration(default, "harmonisation", "harmonisation", ())) != configuration_digest(
        _checkpoint_configuration(changed, "harmonisation", "harmonisation", ()))


def test_saved_evidence_canonical_log_and_report_preserve_strand_and_acceptance(tmp_path):
    study = [row(1, ea="T", oa="C")]  # reverse complement of A/G
    policies = default_policies()
    logger = PipelineLogger("study", "dataset", str(tmp_path), policies=policies,
                            level="INFO", screen_level="ERROR")
    result = infer(tmp_path, study, padded([row(1)]), [row(9999)], logger=logger)
    logger.close()
    assert result["strand_evidence"]["GRCh37"]["reverse"] == 1
    assert result["strand_evidence"]["GRCh37"]["forward"] == 0
    assert json.loads(json.dumps(result))["low_coverage_fallback"]["accepted"] is True
    text = " ".join(render_study_decisions({}, genome_build_info=result).split())
    assert "low coverage fallback" in text.lower()
    assert "exclusive SNP loci" in text
    assert "not a posterior probability" in text
    log = open(logger.log_path).read()
    assert "low_coverage_fallback" in log
    assert "exclusive SNP loci" in log
    assert "Genome-build chromosome" not in log
    stage = next(item for item in logger.summary()["steps"] if "genome_build" in item.get("extra", {}))
    assert stage["extra"]["genome_build"]["per_chromosome"] == result["per_chromosome"]
    assert result["per_chromosome"]["1"]["testable_variants"] == 1


def test_missing_alleles_remain_in_coordinate_testable_denominator(tmp_path):
    study = [row(i) for i in range(20)] + [row(20, ea=None)]
    result = infer(tmp_path, study, padded([row(i) for i in range(21)]), [row(9999)])
    assert result["testable_variants"] == 21
    assert result["matches"]["GRCh37"] == 20
    assert result["low_coverage_fallback"]["accepted"] is True
    assert result["distinct_loci"]["exclusive_snp_loci"]["GRCh37"] == 20


def test_compatible_non_snp_competitor_evidence_excludes_winner_snp_locus(tmp_path):
    study = [row(1), row(1, ea="AA", oa="A")]
    result = infer(tmp_path, study, padded([row(1)]), [row(1, ea="AA", oa="A")])
    assert result["distinct_loci"]["compatible_snp_loci"]["GRCh37"] == 1
    assert result["distinct_loci"]["exclusive_snp_loci"]["GRCh37"] == 0


@pytest.mark.parametrize("key,value", [
    ("low_coverage_fallback_enabled", False),
    ("low_coverage_min_match_fraction", 0.96),
    ("low_coverage_max_competitor_fraction", 0.04),
    ("low_coverage_min_unique_loci", 10001),
    ("low_coverage_chromosome_min_competitor_loci", 101),
    ("low_coverage_chromosome_min_competitor_fraction", 0.96),
])
def test_each_fallback_policy_flows_from_run_yaml_into_scientific_checkpoint(tmp_path, key, value):
    path = tmp_path / "run.yaml"
    path.write_text(yaml.safe_dump({"modules": {"harmonisation": {
        "policies": {"build": {key: value}}}}}))
    changed = load_run_configuration_for_module("harmonisation", path)
    assert changed.modules.harmonisation.policies["build"][key] == value
    default = load_configuration()
    assert configuration_digest(_checkpoint_configuration(default, "harmonisation", "harmonisation", ())) != configuration_digest(
        _checkpoint_configuration(changed, "harmonisation", "harmonisation", ()))


def test_accepted_and_rejected_build_evidence_is_visible_in_html(tmp_path):
    from postgwas.modules.harmonisation.html_report import render_dataset_report
    study = [row(i) for i in range(20)]
    for floor, status in [(20, "accepted"), (21, "rejected")]:
        result = infer(tmp_path, study, padded(study), [row(9999)],
                       overrides={"build.low_coverage_min_unique_loci": floor})
        text = render_dataset_report({"dataset_id": "study", "status": "FAILED",
            "genome_build": result["inferred_build"]}, {"status": "FAILED", "dataset": {
                "genome_build": result}})
        assert "Genome-build acceptance path" in text
        assert "Low coverage fallback" in text
        assert status in text
        assert "Exclusive snp loci" in text
        if status == "rejected":
            assert "build.low_coverage_min_unique_loci" in text


@pytest.mark.parametrize("second", [[row(9999, "2")], [row(9999)]])
def test_fallback_requires_comparable_reference_scope_on_testable_chromosomes(tmp_path, second):
    study = [row(i) for i in range(100)] + [row(1, "2")]
    result = infer(tmp_path, study, padded(study, 2000), second)
    assert result["inferred_build"] == "Ambiguous"
    assert "comparable_reference_scope" in result["low_coverage_fallback"]["failed_gates"]
    assert "comparable_reference_scope" in result["ambiguous_reason"]


def test_nonpositive_and_missing_positions_cannot_inflate_distinct_snp_floor(tmp_path):
    study = [("1", -1, "A", "G"), ("1", 0, "A", "G"), ("1", None, "A", "G"), row(1)]
    result = infer(tmp_path, study, padded(study), [row(9999)],
                   overrides={"build.low_coverage_min_unique_loci": 2})
    assert result["distinct_loci"]["exclusive_snp_loci"]["GRCh37"] == 1
    assert not result["low_coverage_fallback"]["accepted"]


@pytest.mark.parametrize("path,source", [("normal", "detector"), ("low_coverage_fallback", "low_coverage_fallback")])
def test_csv_summary_exposes_fallback_method_without_changing_normal_source(tmp_path, path, source):
    from postgwas.modules.harmonisation import cli
    record = cli._run_summary_record(
        SimpleNamespace(dataset_id="study", input_file="study.tsv"),
        status="FAILED", manifest_path=tmp_path / "manifest.json",
        manifest={"dataset": {"genome_build": {
            "inferred_build": "GRCh37", "acceptance_path": path}}},
        strand_reference_panel="panel", strand_reference_population="EUR",
        comparison_af_panel="panel", comparison_af_population="EUR")
    assert record["genome_build_source"] == source
    assert record["genome_build_detected"] == "GRCh37"


@pytest.mark.parametrize("count,veto", [(94, False), (95, True)])
def test_exact_chromosome_veto_fraction_applies_to_low_coverage_fallback(tmp_path, count, veto):
    main = [row(i) for i in range(2000)]
    minor = [row(i, "2") for i in range(100)]
    result = infer(tmp_path, main + minor,
        padded(main + [row(99999, "2")], 30000),
        [row(99999)] + minor[:count] + [row(i, "2", ea="A", oa="C") for i in range(100)],
        overrides={"build.low_coverage_chromosome_min_competitor_loci": 94})
    assert result["per_chromosome"]["2"]["competitor_veto"] is veto
    assert result["low_coverage_fallback"]["accepted"] is (not veto)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_shared_numeric_validation_preserves_nonfinite_execution_rejection(value):
    with pytest.raises(PolicyError, match="execution.memory_budget_gb"):
        default_policies().with_overrides({"execution.memory_budget_gb": value})


@pytest.mark.parametrize("report_kind", ["dataset", "combined"])
@pytest.mark.parametrize("evidence,expected_note", [
    ({"acceptance_path": "low_coverage_fallback"}, True),
    ({"low_coverage_fallback": {"accepted": True}}, True),
    ({"acceptance_path": "normal", "low_coverage_fallback": {"enabled": True, "accepted": False}}, False),
    ({"acceptance_path": "ambiguous", "low_coverage_fallback": {"status": "rejected", "accepted": False}}, False),
    ({}, False),
])
def test_exported_html_retains_recorded_fallback_interpretation_limits(report_kind, evidence, expected_note):
    from postgwas.modules.harmonisation.html_report import render_dataset_report, render_run_report
    record = {"dataset_id": "study", "genome_build": "GRCh37"}
    manifest = {"dataset": {"genome_build": evidence}}
    original = deepcopy(manifest)
    text = (render_dataset_report(record, manifest) if report_kind == "dataset"
            else render_run_report([record], {"study": manifest}))
    assert ("not a posterior probability" in text) is expected_note
    assert ("small mixed components can pass" in text) is expected_note
    assert manifest == original
