"""Original-EA GWAS-row screening and one strictly screened AF replacement."""
from pathlib import Path
import shutil

import polars as pl
import pytest

from postgwas.modules.harmonisation import allele_frequency as af, strand
from postgwas.modules.harmonisation.rejects import SOURCE_INPUT_ROW_COLUMN
from postgwas.modules.harmonisation.shared import allele_join, statistics
from postgwas.modules.harmonisation.shared.runtime import NullStepContext
from postgwas.modules.harmonisation.study_properties import _decide_eaf_is_maf
from test_harmonisation_external_af_preflight import (
    MAPPING, adopt as adopt_external, fixture as external_fixture,
    validate as validate_external,
)
from test_harmonisation_replacement_af_integration import (
    adopt as adopt_replacement, fixture as replacement_fixture,
    validate as validate_replacement,
)


def _forbid_worker_reanalysis(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("approved chromosome worker repeated frequency screening/classification/mapping")
    for module, names in (
        (af, ("read_reference_variant_table", "allele_oriented_left_join",
              "merge_external_allele_frequencies", "_compare_frequency_pairs")),
        (strand, ("_load_reference",)),
        (statistics, ("frequency_maf_screen",)),
        (allele_join, ("frequency_maf_screen",)),
    ):
        for name in names:
            monkeypatch.setattr(module, name, forbidden)
    if hasattr(af, "frequency_maf_screen"):
        monkeypatch.setattr(af, "frequency_maf_screen", forbidden)


def _screen(evidence, *, rows, usable, low, maf_like):
    screen = evidence["effective_frequency_screen"]
    assert screen["rows"] == rows
    assert screen["usable"] == usable
    assert screen["unusable"] == rows - usable
    assert screen["at_or_below_0_5"] == low
    assert screen["low_fraction_of_usable"] == pytest.approx(low / usable)
    assert screen["maf_like"] is maf_like
    assert screen["scope"] == "merged_gwas_rows"
    assert screen["allele_direction"] == "original_gwas_effect_allele"
    return screen


@pytest.mark.parametrize("external_swap,study_swap", [
    (False, False), (False, True), (True, False), (True, True),
])
def test_effective_external_screen_uses_original_gwas_ea_once_and_workers_reuse(
        tmp_path, monkeypatch, external_swap, study_swap):
    raw = [.1 + i % 40 / 100 for i in range(1200)]
    parts = external_fixture(tmp_path, external_swap=external_swap, study_swap=study_swap,
        values=raw, default_values=[1-x for x in raw] if external_swap else raw)
    original_is_minor = external_swap == study_swap
    if not original_is_minor:
        monkeypatch.setattr(af, "_compare_frequency_pairs",
            lambda *a, **k: pytest.fail("non-MAF-like merged GWAS AF entered semantic comparator"))
    selected = validate_external(parts)
    _screen(selected["eaf_external_validation"], rows=1200, usable=1200,
        low=1200 if original_is_minor else 0, maf_like=original_is_minor)
    assert selected["eaf_selected_source"] == "external_eaf"
    assert parts[1]["eaf_col"] is None and selected["eaf_column"] is None
    _forbid_worker_reanalysis(monkeypatch)
    for chromosome, offset in (("1", 0), ("2", 1)):
        result, qc, columns = adopt_external(parts, selected, chromosome)
        expected = raw[offset::2]
        assert result[columns["eaf_col"]].to_list() == pytest.approx(
            [1-x for x in expected] if external_swap else expected)
        assert result["BETA"].to_list() == ([-.2] * 600 if study_swap else [.2] * 600)
        assert result.height == 600
        assert qc["eaf_external_validation"]["status"] == "passed"


def test_external_reference_only_rows_cannot_turn_merged_non_maf_into_maf(tmp_path):
    # 24,000 panel-only minor rows would make the raw panel 24,000/25,200 > .95 minor.
    # The GWAS receives 1,200 major values, so its actual screen has zero minor values.
    values = [.6 + i % 30 / 100 for i in range(1200)]
    def append_panel_only(source, default, external):
        with external.open("a") as handle:
            handle.writelines("1\t%d\tA\tC\t0.1\n" % (100000+i) for i in range(24000))
        return source
    parts = external_fixture(tmp_path, values=values, default_values=values, modify=append_panel_only)
    selected = validate_external(parts)
    evidence = selected["eaf_external_validation"]
    _screen(evidence, rows=1200, usable=1200, low=0, maf_like=False)
    for chromosome, offset in (("1", 0), ("2", 1)):
        result, _, columns = adopt_external(parts, selected, chromosome)
        assert result[columns["eaf_col"]].to_list() == pytest.approx(values[offset::2])


def test_external_screen_counts_all_merged_rows_before_default_diagnostics(tmp_path):
    # The 1,000 diagnostic rows are minor, but 200 unmatched-default major rows
    # belong to the merged GWAS screen: 1,000/1,200, not 1,000/1,000.
    values = [.1 + i % 30 / 100 for i in range(1000)] + [.7] * 200
    def remove_default_only(source, default, external):
        pl.read_csv(default, separator="\t").filter(pl.col("POS") < 1100).write_csv(default, separator="\t")
        return source
    parts = external_fixture(tmp_path, values=values, default_values=values,
        modify=remove_default_only, policies={"strand.unmatched_action": "retain"})
    selected = validate_external(parts)
    _screen(selected["eaf_external_validation"], rows=1200, usable=1200, low=1000, maf_like=False)
    result, _, columns = adopt_external(parts, selected, "1")
    assert result.height == 600
    assert result[columns["eaf_col"]].to_list() == pytest.approx(values[::2])


@pytest.mark.parametrize("role", ["native", "external"])
def test_resolved_multiallelic_default_diagnostics_keep_raw_a_and_all_pairs(tmp_path, role):
    # Each site has a second, zero-frequency ALT. The matched effect AF remains
    # .1-.49; explicit other=1-f and min(f,1-f) both retain that diagnostic value.
    raw = [.1 + i % 40 / 100 for i in range(1200)]
    def add_alternates(source, default, external):
        with default.open("a") as handle:
            handle.writelines("%s\t%d\tA\tG\t0\n" % (1+i % 2, 100+i) for i in range(1200))
        return source
    parts = external_fixture(tmp_path, values=raw, default_values=raw, modify=add_alternates)
    if role == "external":
        selected = validate_external(parts)
        evidence = selected["eaf_external_validation"]
        result, _, columns = adopt_external(parts, selected, "1")
        assert result[columns["eaf_col"]].to_list() == pytest.approx(raw[::2])
    else:
        source = parts[0].with_columns(pl.Series("INPUT_AF", raw))
        columns = dict(parts[1], eaf_col="INPUT_AF")
        decisions = {"strand": "forward", "effect_type": "beta"}
        _decide_eaf_is_maf(source, columns, parts[4], NullStepContext(), decisions, {})
        selected, prepared = af.validate_study_allele_frequency(source, columns, decisions,
            parts[3], MAPPING, policies=parts[4])
        evidence = selected["eaf_reference_validation"]
        assert prepared.select(source.columns).equals(source)
        assert prepared["INPUT_AF"].to_list() == raw
    assert evidence["status"] == "passed"
    assert evidence["comparable_variants"] == evidence["study_variant_count"] == 1200
    assert evidence["mean_eaf_error"] == pytest.approx(0.)
    assert evidence["rule_1_passed"] and not evidence["rule_2_passed"]


def test_replacement_maf_like_after_merge_stops_before_consistency_comparator(tmp_path, monkeypatch):
    # Replacement values are .1-.4 in the original GWAS EA direction. Half need
    # an allele swap to reach that direction; screen the merged values, not ALT.
    original = [.8, .7, .6, .1, .2, .3] * 2
    b = [min(x, 1-x) for x in original]
    source_alt = [1-x if i % 4 in (1, 3) else x for i, x in enumerate(b)]
    parts = replacement_fixture(tmp_path, replacement_values=source_alt)
    before = {p: Path(p).read_bytes() for p in parts[5].paths.values()}
    monkeypatch.setattr(af, "validate_replacement_frequency_pairs",
        lambda *a, **k: pytest.fail("MAF-like replacement reached consistency comparison"))
    monkeypatch.setattr(af, "_compare_frequency_pairs",
        lambda *a, **k: pytest.fail("MAF-like replacement reran default semantic comparison"))
    with pytest.raises(af.StudyAlleleFrequencyError, match="(?i)replacement.*MAF") as caught:
        validate_replacement(parts)
    evidence = caught.value.study_decisions["eaf_replacement_validation"]
    _screen(evidence, rows=12, usable=12, low=12, maf_like=True)
    assert evidence["status"] != "passed" and "contract" not in evidence
    assert {p: Path(p).read_bytes() for p in before} == before
    assert not Path(parts[1]["output_folder"]).exists()


def test_replacement_non_maf_screen_keeps_consistency_gate_and_original_a(tmp_path, monkeypatch):
    parts = replacement_fixture(tmp_path)
    selected = validate_replacement(parts)
    evidence = selected["eaf_replacement_validation"]
    _screen(evidence, rows=12, usable=12, low=6, maf_like=False)
    assert evidence["correlation"] == pytest.approx(1.)
    assert evidence["mean_absolute_difference"] == pytest.approx(0.)
    assert evidence["major_count"] == 6
    _forbid_worker_reanalysis(monkeypatch)
    result, qc, columns = adopt_replacement(parts, selected)
    assert result["MAF"].equals(parts[0]["MAF"])
    expected = [1-x if i % 4 in (1, 3) else x for i, x in enumerate(parts[6])]
    assert result[columns["eaf_col"]].to_list() == pytest.approx(expected)
    assert qc["eaf_replacement_validation"]["status"] == "passed"


def _external_maf_fixture(tmp_path, *, replacement_file=True, replacement_column=True):
    default_ea = [.8, .7, .6, .1, .2, .3] * 200
    raw_a = [min(x, 1-x) for x in default_ea]
    parts = list(external_fixture(tmp_path, values=raw_a, default_values=default_ea))
    replacement = tmp_path / "replacement.tsv"
    pl.read_csv(tmp_path / "default.tsv", separator="\t").write_csv(replacement, separator="\t")
    selected, work = af.validate_study_allele_frequency(parts[0], parts[1],
        {"eaf_column": None, "eaf_is_maf": None, "eaf_is_maf_source": "no_study_frequency_column",
         "effect_type": "beta", "strand": "forward"}, parts[3], MAPPING, policies=parts[4],
        eaf_replacement_file=str(replacement) if replacement_file else None,
        eaf_replacement_column="AF" if replacement_column else None)
    parts[2] = selected
    for chromosome, path in parts[5].paths.items():
        work.filter(pl.col("CHR") == chromosome).write_parquet(path)
    return parts, replacement, raw_a, default_ea


@pytest.mark.parametrize("file_given,column_given", [(False, False), (True, False), (False, True)])
def test_external_maf_transition_requires_both_public_replacement_arguments(
        tmp_path, file_given, column_given):
    missing = []
    if not file_given:
        missing.append("--eaf-replacement-file")
    if not column_given:
        missing.append("--eaf-replacement-column")
    with pytest.raises(af.StudyAlleleFrequencyError) as caught:
        parts, _, _, _ = _external_maf_fixture(tmp_path,
            replacement_file=file_given, replacement_column=column_given)
        validate_external(parts)
    message = str(caught.value)
    for option in missing:
        assert option in message
        assert "modules.harmonisation." + option[2:].replace("-", "_") in message


def test_external_maf_to_replacement_preserves_absent_native_and_frozen_external_a(tmp_path, monkeypatch):
    parts, replacement, raw_a, default_ea = _external_maf_fixture(tmp_path)
    selected = validate_external(parts)
    assert selected["eaf_selected_source"] == "replacement_reference"
    assert selected["eaf_reference_validation"]["status"] == "replacement_selected"
    assert selected["eaf_input_interpretation"] == "maf"
    assert selected["eaf_column"] is None and parts[1]["eaf_col"] is None
    assert "AF" not in parts[0].columns
    initial_column = selected["eaf_initial_frequency"]["column"]
    for chromosome, offset in (("1", 0), ("2", 1)):
        original = pl.read_parquet(parts[5].paths[chromosome])
        assert original[initial_column].to_list() == pytest.approx(raw_a[offset::2])
    # Resource activation/staging has its own service integration tests. Rebind
    # the declared replacement here before invoking its real global validator.
    maps = {c: dict(v, user_eaf_file=str(replacement)) for c, v in parts[3].items()}
    selected["eaf_replacement"]["resolved_files"] = {c: str(replacement) for c in maps}
    selected["eaf_replacement"]["declared_files"] = {c: str(replacement) for c in maps}
    approved = af.validate_replacement_allele_frequency(parts[5], parts[1], selected,
        maps, MAPPING, policies=parts[4])
    evidence = approved["eaf_replacement_validation"]
    _screen(evidence, rows=1200, usable=1200, low=600, maf_like=False)
    assert evidence["comparable_variants"] == 1200
    assert evidence["correlation"] == pytest.approx(1.)
    assert evidence["mean_absolute_difference"] == pytest.approx(0.)
    for chromosome, offset in (("1", 0), ("2", 1)):
        prepared = pl.read_parquet(parts[5].paths[chromosome])
        assert prepared[initial_column].to_list() == pytest.approx(raw_a[offset::2])
    _forbid_worker_reanalysis(monkeypatch)
    adopted = (parts[0], parts[1], parts[2], maps, parts[4], parts[5])
    for chromosome, offset in (("1", 0), ("2", 1)):
        result, qc, columns = adopt_external(adopted, approved, chromosome)
        assert result[columns["eaf_col"]].to_list() == pytest.approx(default_ea[offset::2])
        assert result.height == 600 and qc["eaf_selected_source"] == "replacement_reference"
        # A is private evidence, never fabricated as a native public AF column.
        assert approved["eaf_column"] is None and parts[1]["eaf_col"] is None
        assert qc["eaf_input_interpretation"] == "maf"
    assert parts[0].columns == [SOURCE_INPUT_ROW_COLUMN, "CHR", "POS", "EA", "OA", "BETA", "Z"]
    assert raw_a == [min(x, 1-x) for x in default_ea]


@pytest.mark.parametrize("replacement_is_maf", [False, True])
def test_service_activates_and_stages_external_maf_replacement_before_workers(
        tmp_path, monkeypatch, replacement_is_maf, _staged_b_mutation=False):
    from postgwas.config import load_configuration
    from postgwas.core.pipeline_logging import PipelineLogger
    from postgwas.modules.harmonisation import service
    from postgwas.modules.harmonisation.resource_preflight import validate_harmonisation_resource_maps

    default_ea = [.8, .7, .6, .1, .2, .3] * 200
    raw_a = [min(x, 1-x) for x in default_ea]
    parts = external_fixture(tmp_path / "inputs", values=raw_a, default_values=default_ea)
    source, columns, _, maps, policies, _ = parts
    replacement = tmp_path / "replacement.tsv"
    candidate_b = raw_a if replacement_is_maf else default_ea
    pl.read_csv(tmp_path / "inputs" / "default.tsv", separator="\t").with_columns(
        pl.Series("AF", candidate_b)).write_csv(replacement, separator="\t")
    replacement_snapshot = replacement.read_bytes()
    study = tmp_path / "study.tsv"
    source.drop(SOURCE_INPUT_ROW_COLUMN).with_columns(pl.lit(.1).alias("SE"),
        pl.lit(.0455002639).alias("P"), pl.lit(.9).alias("INFO")).write_csv(study, separator="\t")
    module = load_configuration().modules.harmonisation
    monkeypatch.setattr(service, "infer_genome_build", lambda *a, **k: {
        "inferred_build": "GRCh37", "input_variants": 1200,
        "testable_variants": 1200, "untestable_variants": 0})
    monkeypatch.setattr(service, "resolve_strand_consensus", lambda *a, **k: {
        "strand": "forward", "forward": 1200, "reverse": 0, "ambiguous": 0,
        "dominant_fraction": 1., "reason": "known-forward fixture"})

    def preflight(*args, **kwargs):
        for resources in maps.values():
            resources["user_eaf_column"] = "AF"
        checked = validate_harmonisation_resource_maps(maps, grch_version="GRCh37", bcftools="unused",
            vcf_config=module.vcf_processing.model_dump(), default_eaf_colmap=MAPPING,
            external_eaf_colmap=MAPPING, external_info_colmap={}, policies=policies,
            require_default_eaf=True, resource_keys=("default_eaf_file",),
            validate_frequency_semantics=kwargs.get("validate_frequency_semantics", True))
        return maps, checked

    monkeypatch.setattr(service, "preflight_harmonisation_resources", preflight)
    monkeypatch.setattr(service, "_run_one_round",
        lambda *a, **k: pytest.fail("frequency preparation failure launched workers or retried"))
    stage = service.stage_shared_external_reference_files
    staged_sources = []

    def observed_stage(resource_maps, **kwargs):
        staged_sources.append({v.get("user_eaf_file") for v in resource_maps.values()})
        staged_maps, summary = stage(resource_maps, **kwargs)
        if _staged_b_mutation and kwargs["user_eaf_specification"] == str(replacement):
            # The producer summary remains unchanged: its recorded identities
            # must detect modification after it returns, before measuring B.
            for path in {v["user_eaf_file"] for v in staged_maps.values()}:
                staged = pl.read_parquet(path)
                staged.with_columns((1. - pl.col("AF").cast(pl.Float64)).alias("AF")).write_parquet(path)
        return staged_maps, summary

    monkeypatch.setattr(service, "stage_shared_external_reference_files", observed_stage)
    columns.update(sumstat_file=str(study), eaffile=maps["1"]["user_eaf_file"],
        se_col="SE", pval_col="P", declared_effect_type="beta", declared_pvalue_type="raw",
        ncontrol="1000", imp_info_col="INFO", info_source="internal", output_folder=str(tmp_path / "out"))
    bcftools = shutil.which("bcftools")
    assert bcftools is not None, "Service integration requires the installed PostGWAS executable environment."
    with PipelineLogger("study", "dataset", str(tmp_path / "logs")) as logger:
        arguments = dict(sumstat_file=str(study), sample_column_dict=columns,
            output_dir=str(tmp_path / "out"), resource_folder=str(tmp_path),
            default_eaf_reference_source="fixture", default_comparison_af_file="fixture", dbsnp="fixture",
            resource_layout={}, output_layout=dict(module.output_layout.root),
            gwas2vcf_input=module.gwas2vcf_input.model_dump(), executables={"bcftools": bcftools},
            vcf_config=module.vcf_processing.model_dump(), user_eaf_file=maps["1"]["user_eaf_file"],
            user_eaf_column="AF", user_info_file=None, user_info_column=None,
            default_eaf_reference_column="AF", default_comparison_af_column="AF", build_reference_files={},
            build_check_colmap={}, default_eaf_colmap=MAPPING, external_eaf_colmap=MAPPING,
            external_info_colmap={}, external_reference_staging=module.external_reference_staging.model_dump(),
            policies=policies, logger=logger, input_line_count=source.height,
            max_workers=1, gwas2vcf_main_script_path="unused", eaf_replacement_file=str(replacement),
            eaf_replacement_column="AF")
        if _staged_b_mutation:
            arguments["output_dir_path"] = Path(arguments.pop("output_dir"))
            arguments.pop("max_workers")
            arguments.pop("gwas2vcf_main_script_path")
            with pytest.raises(af.StudyAlleleFrequencyError) as caught:
                service._prepare_dataset_for_chromosome_processing(**arguments)
            evidence = caught.value.study_decisions["eaf_replacement_validation"]
            assert evidence["status"] == "invalid" and "contract" not in evidence
            assert evidence.get("failure_reason", evidence.get("decision_reason")) in (
                "replacement_inputs_changed", "replacement_staged_inputs_changed")
            assert any(str(replacement) in names for names in staged_sources)
            assert replacement.read_bytes() == replacement_snapshot
            assert not list((tmp_path / "out").rglob("*.vcf*"))
            assert logger.summary()["steps"][-1]["status"] == "failed"
            return
        if replacement_is_maf:
            with pytest.raises(af.StudyAlleleFrequencyError, match="(?i)replacement.*MAF") as caught:
                service.harmonise_chromosomes(**arguments)
            decisions = caught.value.study_decisions
            evidence = decisions["eaf_replacement_validation"]
            _screen(evidence, rows=1200, usable=1200, low=1200, maf_like=True)
            assert evidence["status"] != "passed" and "contract" not in evidence
            assert decisions["eaf_column"] is None
            assert not list((tmp_path / "out").rglob("*.vcf*"))
            assert logger.summary()["steps"][-1]["status"] == "failed"
            return
        arguments["output_dir_path"] = Path(arguments.pop("output_dir"))
        arguments.pop("max_workers")
        arguments.pop("gwas2vcf_main_script_path")
        prepared = service._prepare_dataset_for_chromosome_processing(**arguments)

    selected = prepared.study_decisions
    assert selected["eaf_selected_source"] == "replacement_reference"
    assert selected["eaf_replacement_validation"]["status"] == "passed"
    assert selected["eaf_column"] is None and prepared.sample_columns["eaf_col"] is None
    for key, value in selected.items():
        if key.startswith("eaf_"):
            assert key in prepared.dataset_study_decisions
            assert prepared.dataset_study_decisions[key] == value, "stale dataset decision: " + key
    assert selected["eaf_replacement"]["resolved_files"] == {c: str(replacement) for c in maps}
    assert any(str(replacement) in names for names in staged_sources)
    _forbid_worker_reanalysis(monkeypatch)
    for chromosome, offset in (("1", 0), ("2", 1)):
        resources = prepared.resource_maps[chromosome]
        assert resources["user_eaf_file"] != str(replacement)
        assert pl.read_parquet(resources["user_eaf_file"])["AF"].cast(pl.Float64).to_list() == pytest.approx(
            candidate_b[offset::2])
        frame = pl.read_parquet(prepared.chromosome_files[chromosome])
        initial_column = selected["eaf_initial_frequency"]["column"]
        assert frame[initial_column].to_list() == pytest.approx(raw_a[offset::2])
        result, qc, output_columns = af.harmonise_allele_frequency(chromosome, frame,
            dict(prepared.sample_columns), dict(module.output_layout.root), "\t",
            eaffile=resources["user_eaf_file"], external_eaf_colmap=MAPPING,
            default_eaf_file=resources["default_eaf_file"], default_eaf_column="AF",
            default_eaf_colmap=MAPPING, policies=policies, study_decision=selected)
        assert result.height == 600
        assert result[output_columns["eaf_col"]].to_list() == pytest.approx(default_ea[offset::2])
        assert initial_column not in result.columns
        assert qc["eaf_selected_source"] == "replacement_reference"


def test_service_rejects_replacement_partition_modified_after_staging_proof(tmp_path, monkeypatch):
    test_service_activates_and_stages_external_maf_replacement_before_workers(
        tmp_path, monkeypatch, False, _staged_b_mutation=True)


def test_external_a_private_frequency_mutation_stops_before_replacement_measurement(tmp_path, monkeypatch):
    parts, replacement, _, _ = _external_maf_fixture(tmp_path)
    selected = validate_external(parts)
    initial_column = selected["eaf_initial_frequency"]["column"]
    for path in parts[5].paths.values():
        frame = pl.read_parquet(path)
        frame.with_columns((pl.col(initial_column) + .001).alias(initial_column)).write_parquet(path)
    changed = {p: Path(p).read_bytes() for p in parts[5].paths.values()}
    maps = {c: dict(v, user_eaf_file=str(replacement)) for c, v in parts[3].items()}
    selected["eaf_replacement"]["resolved_files"] = {c: str(replacement) for c in maps}
    selected["eaf_replacement"]["declared_files"] = {c: str(replacement) for c in maps}
    monkeypatch.setattr(af, "validate_replacement_frequency_pairs",
        lambda *a, **k: pytest.fail("mutated frozen A reached replacement measurement"))
    with pytest.raises(af.StudyAlleleFrequencyError) as caught:
        af.validate_replacement_allele_frequency(parts[5], parts[1], selected,
            maps, MAPPING, policies=parts[4])
    evidence = caught.value.study_decisions.get("eaf_replacement_validation") or {}
    assert evidence.get("status") != "passed" and "contract" not in evidence
    assert {p: Path(p).read_bytes() for p in changed} == changed


def test_unrelated_frequency_named_columns_do_not_change_external_selection(tmp_path):
    raw = [.1 + i % 40 / 100 for i in range(1200)]
    def extra_columns(source, default, external):
        return source.with_columns(pl.lit(.91).alias("effect_frequency"),
            pl.lit(.92).alias("default_frequency"))
    parts = external_fixture(tmp_path, values=raw, default_values=raw, modify=extra_columns)
    selected = validate_external(parts)
    evidence = selected["eaf_external_validation"]
    _screen(evidence, rows=1200, usable=1200, low=1200, maf_like=True)
    assert evidence["comparable_variants"] == 1200
    assert evidence["mean_eaf_error"] == pytest.approx(0.)
    for chromosome, offset in (("1", 0), ("2", 1)):
        result, _, columns = adopt_external(parts, selected, chromosome)
        assert result[columns["eaf_col"]].to_list() == pytest.approx(raw[offset::2])
        assert result["effect_frequency"].to_list() == [.91] * 600
        assert result["default_frequency"].to_list() == [.92] * 600


@pytest.mark.parametrize("frequency_column", ["effect_frequency", "default_frequency"])
def test_native_a_frequency_name_cannot_shadow_merged_replacement_b(tmp_path, frequency_column):
    parts = replacement_fixture(tmp_path, frequency_column=frequency_column)
    selected = validate_replacement(parts)
    evidence = selected["eaf_replacement_validation"]
    _screen(evidence, rows=12, usable=12, low=6, maf_like=False)
    assert evidence["comparable_variants"] == 12
    assert evidence["major_count"] == 6
    assert evidence["correlation"] == pytest.approx(1.)
    assert evidence["mean_absolute_difference"] == pytest.approx(0.)
    result, _, columns = adopt_replacement(parts, selected)
    assert result[frequency_column].equals(parts[0][frequency_column])
    expected = [1-x if i % 4 in (1, 3) else x for i, x in enumerate(parts[6])]
    assert result[columns["eaf_col"]].to_list() == pytest.approx(expected)


@pytest.mark.parametrize("key", ["chr", "pos", "ea", "oa"])
@pytest.mark.parametrize("name", ["effect_frequency", "default_frequency"])
def test_mapped_variant_column_names_cannot_shadow_frequency_diagnostic_values(tmp_path, key, name):
    raw = [.1 + i % 40 / 100 for i in range(1200)]
    parts = list(external_fixture(tmp_path, values=raw, default_values=raw))
    original_name = parts[1][key + "_col"]
    parts[0] = parts[0].rename({original_name: name})
    parts[1] = dict(parts[1], **{key + "_col": name})
    for path in parts[5].paths.values():
        pl.read_parquet(path).rename({original_name: name}).write_parquet(path)
    selected = validate_external(parts)
    evidence = selected["eaf_external_validation"]
    _screen(evidence, rows=1200, usable=1200, low=1200, maf_like=True)
    assert evidence["comparable_variants"] == 1200
    assert evidence["mean_eaf_error"] == pytest.approx(0.)
    for chromosome, offset in (("1", 0), ("2", 1)):
        result, _, columns = adopt_external(parts, selected, chromosome)
        assert result[columns["eaf_col"]].to_list() == pytest.approx(raw[offset::2])
        original = parts[0].filter(pl.col(parts[1]["chr_col"]) == chromosome)
        assert result[name].equals(original[name])
