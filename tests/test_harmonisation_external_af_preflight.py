"""Initial external AF is approved globally, without inventing native study AF."""
from copy import deepcopy
from pathlib import Path

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation import allele_frequency as af
from postgwas.modules.harmonisation.chromosome_partition import ChromosomePartitions
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.rejects import SOURCE_INPUT_ROW_COLUMN

MAPPING = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}


def fixture(tmp_path, *, count=1200, matched=None, external_swap=False, study_swap=False,
            values=None, default_values=None, same=False, policies=None, modify=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    matched = count if matched is None else matched
    chrom = [str(1 + i % 2) for i in range(count)]
    raw = [.1 + (i % 40) / 100 for i in range(count)] if values is None else values
    default_af = raw if default_values is None else default_values
    source = pl.DataFrame({"CHR": chrom, "POS": list(range(100, 100 + count)),
        "EA": ["A" if study_swap else "C"] * count, "OA": ["C" if study_swap else "A"] * count,
        "BETA": [.2] * count, "Z": [2.] * count}).with_row_index(SOURCE_INPUT_ROW_COLUMN)
    columns = {"chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA", "eaf_col": None,
        "eafcolumn": "AF", "beta_or_col": "BETA", "imp_z_col": "Z", "gwas_outputname": "study",
        "output_folder": str(tmp_path / "out")}
    ref = pl.DataFrame({"CHROM": chrom[:matched], "POS": list(range(100, 100 + matched)),
        "REF": ["A"] * matched, "ALT": ["C"] * matched, "AF": default_af[:matched]})
    default = tmp_path / "default.tsv"
    ref.write_csv(default, separator="\t")
    external = default if same else tmp_path / "external.tsv"
    if not same:
        ref.with_columns(pl.Series("AF", raw[:matched]),
            pl.lit("C" if external_swap else "A").alias("REF"),
            pl.lit("A" if external_swap else "C").alias("ALT")).write_csv(external, separator="\t")
    if modify:
        source = modify(source, default, external)
    maps = {c: {"default_eaf_file": str(default), "default_eaf_reference_column": "AF",
        "user_eaf_file": str(external)} for c in ("1", "2")}
    pol = default_policies().with_overrides(policies or {})
    decisions, work = af.validate_study_allele_frequency(source, columns,
        {"eaf_column": None, "eaf_is_maf": None, "eaf_is_maf_source": "no_study_frequency_column",
         "effect_type": "beta", "strand": "forward"}, maps, MAPPING, policies=pol)
    paths, rows = {}, {}
    for c in maps:
        path = tmp_path / ("work_" + c + ".parquet")
        subset = work.filter(pl.col("CHR") == c)
        subset.write_parquet(path)
        paths[c], rows[c] = str(path), subset.height
    return source, columns, decisions, maps, pol, ChromosomePartitions(paths, rows)


def validate(parts):
    _, columns, decisions, maps, pol, partitions = parts
    return af.validate_external_allele_frequency(partitions, columns, decisions, maps, MAPPING,
        MAPPING, original_resource_maps=maps, policies=pol)


def adopt(parts, decisions, chromosome="1", frame=None):
    _, columns, _, maps, pol, partitions = parts
    return af.harmonise_allele_frequency(chromosome,
        pl.read_parquet(partitions.paths[chromosome]) if frame is None else frame,
        dict(columns), dict(load_configuration().modules.harmonisation.output_layout.root), "\t",
        eaffile=maps[chromosome]["user_eaf_file"], external_eaf_colmap=MAPPING,
        default_eaf_file=maps[chromosome]["default_eaf_file"], default_eaf_column="AF",
        default_eaf_colmap=MAPPING, policies=pol, study_decision=decisions)


@pytest.mark.parametrize("external_swap,study_swap", [(False, False), (False, True), (True, False), (True, True)])
def test_pooled_external_screen_uses_original_gwas_effect_allele_and_never_reclassifies(tmp_path, monkeypatch, external_swap, study_swap):
    raw = [.1 + i % 40 / 100 for i in range(1200)]
    parts = fixture(tmp_path, external_swap=external_swap, study_swap=study_swap, values=raw,
        default_values=[1-x for x in raw] if external_swap else raw)
    decisions = validate(parts)
    evidence = decisions["eaf_external_validation"]
    assert evidence["status"] == "passed"
    original_ea_minor = external_swap == study_swap
    screen = evidence["effective_frequency_screen"]
    assert screen["rows"] == screen["usable"] == evidence["study_variant_count"] == 1200
    assert screen["at_or_below_0_5"] == (1200 if original_ea_minor else 0)
    assert screen["maf_like"] is original_ea_minor
    assert screen["scope"] == "merged_gwas_rows"
    assert screen["allele_direction"] == "original_gwas_effect_allele"
    if original_ea_minor:
        assert evidence["comparable_variants"] == 1200
        assert evidence["rule_1_passed"] and not evidence["rule_2_passed"]
        assert evidence["direct_eaf_correlation"] == pytest.approx(1.)
        assert evidence["mean_eaf_error"] == pytest.approx(0.)
    else:
        assert evidence["decision_reason"] == "external_frequency_not_predominantly_minor"
        assert "comparable_variants" not in evidence and "rule_1_passed" not in evidence
    assert evidence["raw_frequency_screen"]["usable"] == 1200
    assert decisions["eaf_column"] is None and decisions["eaf_is_maf"] is None
    assert "eaf_replacement" not in decisions
    for name in ("read_reference_variant_table", "allele_oriented_left_join", "_compare_frequency_pairs"):
        monkeypatch.setattr(af, name, lambda *a, **k: pytest.fail("approved worker repeated source classification/matching"))
    for c in ("1", "2"):
        result, qc, columns = adopt(parts, decisions, c)
        expected = parts[0].filter(pl.col("CHR") == c)
        values = raw[0 if c == "1" else 1::2]
        assert result[columns["eaf_col"]].to_list() == pytest.approx([1-x for x in values] if external_swap else values)
        assert result["BETA"].to_list() == ([-.2] * 600 if study_swap else [.2] * 600)
        assert result[SOURCE_INPUT_ROW_COLUMN].equals(expected[SOURCE_INPUT_ROW_COLUMN])
        assert qc["eaf_external_validation"]["status"] == "passed"


def test_fixed_denominator_coverage_gate_stops_before_metrics_and_keeps_shards(tmp_path, monkeypatch):
    parts = fixture(tmp_path, count=2402, matched=1200)
    before = {p: Path(p).read_bytes() for p in parts[5].paths.values()}
    monkeypatch.setattr(af, "_measure_frequency_pairs", lambda *a, **k: pytest.fail("coverage failure ran metrics"))
    with pytest.raises(af.StudyAlleleFrequencyError, match="insufficient_overlap_fraction") as caught:
        validate(parts)
    evidence = caught.value.study_decisions["eaf_external_validation"]
    assert evidence["study_variant_count"] == 2402 and evidence["comparable_variants"] == 1200
    assert evidence["observed_overlap_fraction"] == 1200 / 2402
    assert not evidence["shared_evidence_gate_passed"] and "contract" not in evidence
    assert {p: Path(p).read_bytes() for p in before} == before


def test_fixed_denominator_exact_half_passes(tmp_path):
    parts = fixture(tmp_path, count=2400, matched=1200)
    assert validate(parts)["eaf_external_validation"]["observed_overlap_fraction"] == .5


def test_missing_external_global_approval_cannot_fall_back(tmp_path, monkeypatch):
    parts = fixture(tmp_path)
    monkeypatch.setattr(af, "harmonise_strand_orientation", lambda *a, **k: pytest.fail("missing approval reached strand"))
    with pytest.raises(af.AlleleFrequencyError, match="external.*validation"):
        adopt(parts, parts[2])


@pytest.mark.parametrize("kind,reason", [("maf", "maf_interpretation_supported"), ("constant", "both_rules_failed"),
    ("both", "both_rules_passed"), ("missing", "no_usable_external_frequencies"), ("invalid", "no_usable_external_frequencies")])
def test_external_maf_inconclusive_and_unusable_sources_stop_without_replacement(tmp_path, kind, reason):
    x = [.1, .2, .3, .4, .49, .48] * 200
    r = [.1, .2, .3, .4, .51, .52] * 200
    policies = {}
    if kind == "maf":
        r = [.8, .7, .6, .1, .2, .3] * 200
        x = [min(v, 1-v) for v in r]
    elif kind == "constant":
        x = [.25] * 1200
    elif kind == "both":
        x = r
        policies = {"eaf.study_reference_min_minor_fraction": .6, "eaf.maf_decision_cutoff": .6}
    elif kind in ("missing", "invalid"):
        x = [None if kind == "missing" else -1.] * 1200
    parts = fixture(tmp_path, values=x, default_values=r, policies=policies)
    original = {p: Path(p).read_bytes() for p in parts[5].paths.values()}
    message = "(?s)MAF.*eaf-replacement-file.*eaf-replacement-column" if kind == "maf" else reason
    with pytest.raises(af.StudyAlleleFrequencyError, match=message) as caught:
        validate(parts)
    decisions = caught.value.study_decisions
    assert decisions["eaf_external_validation"]["decision_reason"] == reason
    if kind == "maf":
        assert decisions["eaf_external_validation"]["failure_reason"] == "replacement_arguments_missing"
    assert "contract" not in decisions["eaf_external_validation"] and "eaf_replacement" not in decisions
    assert {p: Path(p).read_bytes() for p in original} == original


@pytest.mark.parametrize("alias,population", [(False, False), (True, False), (True, True)])
def test_same_physical_file_never_supplies_independent_semantic_evidence(tmp_path, alias, population):
    parts = fixture(tmp_path, same=True)
    path = Path(parts[3]["1"]["user_eaf_file"])
    if population:
        pl.read_csv(path, separator="\t").with_columns(pl.col("AF").alias("OTHER")).write_csv(path, separator="\t")
        parts[1]["eafcolumn"] = "OTHER"
    if alias:
        link = tmp_path / "alias.tsv"
        link.symlink_to(path)
        for resources in parts[3].values():
            resources["user_eaf_file"] = str(link)
    with pytest.raises(af.StudyAlleleFrequencyError, match="independent_reference_required"):
        validate(parts)


def test_descriptive_panel_counts_and_decisive_gwas_screen_ignore_unobserved_chromosomes(tmp_path):
    parts = fixture(tmp_path)
    path = Path(parts[3]["1"]["user_eaf_file"])
    reference = pl.read_csv(path, separator="\t")
    extra = reference.with_columns(pl.lit(3).alias("CHROM"), pl.lit(.9).alias("AF"))
    pl.concat([reference, reference, extra, extra], how="vertical_relaxed").write_csv(path, separator="\t")
    result = validate(parts)
    screen = result["eaf_external_validation"]["raw_frequency_screen"]
    assert screen["rows"] == screen["usable"] == screen["at_or_below_0_5"] == 1200
    assert screen["maf_like"] is True
    effective = result["eaf_external_validation"]["effective_frequency_screen"]
    assert effective["rows"] == effective["usable"] == effective["at_or_below_0_5"] == 1200
    assert effective["scope"] == "merged_gwas_rows" and effective["maf_like"] is True


@pytest.mark.parametrize("change", ["source", "partition", "mapping", "column", "descriptor", "cache", "approval"])
def test_changed_external_binding_stops_before_strand_without_source_retry(tmp_path, monkeypatch, change):
    parts = fixture(tmp_path)
    decisions = validate(parts)
    frame = pl.read_parquet(parts[5].paths["1"])
    if change == "source":
        path = Path(parts[3]["1"]["user_eaf_file"])
        path.write_text(path.read_text() + "\n")
    elif change == "partition":
        frame.with_columns(pl.lit(.3).alias("BETA")).write_parquet(parts[5].paths["1"])
    elif change == "mapping":
        decisions["eaf_external_source"]["column"] = "DIFFERENT"
    elif change == "column":
        parts[1]["eafcolumn"] = "DIFFERENT"
    elif change == "descriptor":
        decisions["eaf_external_validation"]["contract"]["chromosomes"]["1"]["assignment_cache"]["assignment"]["value"] = "changed"
    elif change == "cache":
        cache = decisions["eaf_external_validation"]["contract"]["chromosomes"]["1"]["assignment_cache"]["column"]
        frame = frame.with_columns(pl.lit(None).cast(frame.schema[cache]).alias(cache))
    else:
        decisions.pop("eaf_external_validation")
    monkeypatch.setattr(af, "harmonise_strand_orientation", lambda *a, **k: pytest.fail("stale approval reached strand"))
    monkeypatch.setattr(af, "read_reference_variant_table", lambda *a, **k: pytest.fail("stale approval reread source"))
    with pytest.raises(af.AlleleFrequencyError):
        adopt(parts, decisions, frame=frame)


@pytest.mark.parametrize("when", ["after_metrics", "after_stream", "partial"])
def test_external_publication_never_approves_mutated_or_partly_written_shards(tmp_path, monkeypatch, when):
    parts = fixture(tmp_path)
    path = Path(parts[5].paths["1"])
    changed = pl.read_parquet(path).with_columns(pl.lit(.3).alias("BETA"))
    if when == "after_metrics":
        compare = af._compare_frequency_pairs
        def mutation(*a, **k):
            result = compare(*a, **k)
            changed.write_parquet(path)
            return result
        monkeypatch.setattr(af, "_compare_frequency_pairs", mutation)
    elif when == "after_stream":
        sink = pl.LazyFrame.sink_parquet
        def mutation(*a, **k):
            result = sink(*a, **k)
            changed.write_parquet(path)
            return result
        monkeypatch.setattr(pl.LazyFrame, "sink_parquet", mutation)
    else:
        writer = af._write_parquet_atomically
        calls = []
        def mutation(*a, **k):
            calls.append(1)
            if len(calls) == 2:
                raise OSError("simulated publication failure")
            return writer(*a, **k)
        monkeypatch.setattr(af, "_write_parquet_atomically", mutation)
    with pytest.raises(af.StudyAlleleFrequencyError) as caught:
        validate(parts)
    evidence = caught.value.study_decisions["eaf_external_validation"]
    assert evidence["status"] == "invalid"
    assert "contract" not in evidence
    assert "replacement" not in str(caught.value).lower()
    assert evidence["failure_reason"] == ("external_inputs_changed" if when == "after_metrics" else "external_cache_publication_failed")


@pytest.mark.parametrize("mode", ["inconclusive", "passed", "coverage", "same_passed", "same_maf", "mutate_content", "mutate_alias", "mutate_default_alias", "mutate_staged", "mutate_staged_alias", "mutate_during_staging", "spearman_budget"])
def test_real_service_missing_native_af_enters_global_preparation_and_fails_before_workers(tmp_path, monkeypatch, mode):
    from postgwas.core.pipeline_logging import PipelineLogger
    from postgwas.modules.harmonisation import service
    parts = fixture(tmp_path / "inputs", count=2402 if mode == "coverage" else 1200,
        matched=1200, same=mode.startswith("same") or mode == "mutate_staged_alias",
        policies={"eaf.maf_reference_correlation_method": "spearman", "execution.memory_budget_gb": .001} if mode == "spearman_budget" else None,
        values=([.1 + i % 20 / 100 if i % 2 == 0 else .6 + i % 20 / 100 for i in range(1200)]
            if mode in ("same_passed", "mutate_staged_alias") else [.25] * 1200 if mode == "inconclusive" else None))
    source, columns, _, maps, pol, _ = parts
    if mode == "same_passed":
        path = Path(maps["1"]["user_eaf_file"])
        with path.open("a") as handle:
            for position in range(10000, 60000):
                handle.write("3\t%s\tA\tC\t.1\n" % position)
    if mode in ("mutate_alias", "mutate_default_alias"):
        role = "default_eaf_file" if mode == "mutate_default_alias" else "user_eaf_file"
        original = Path(maps["1"][role])
        alias = tmp_path / "source_alias.tsv"
        alias.symlink_to(original)
        for resources in maps.values():
            resources[role] = str(alias)
    study = tmp_path / "study.tsv"
    source.drop(SOURCE_INPUT_ROW_COLUMN).with_columns(pl.lit(.1).alias("SE"),
        pl.lit(.0455002639).alias("P"), pl.lit(.9).alias("INFO")).write_csv(study, separator="\t")
    module = load_configuration().modules.harmonisation
    summary = dict.fromkeys(("unique_resource_files_checked", "indexed_vcfs_checked",
        "frequency_vcf_headers_checked", "fasta_indexes_checked", "table_headers_checked", "structured_text_files_checked"), 0)
    summary.update(chromosomes=["1", "2"], file_identities=[], resources_by_chromosome=deepcopy(maps))
    monkeypatch.setattr(service, "infer_genome_build", lambda *a, **k: {"inferred_build": "GRCh37",
        "input_variants": 1200, "testable_variants": 1200, "untestable_variants": 0})
    monkeypatch.setattr(service, "resolve_strand_consensus", lambda *a, **k: {"strand": "forward",
        "forward": 1200, "reverse": 0, "ambiguous": 0, "dominant_fraction": 1., "reason": "known-forward fixture"})
    from postgwas.modules.harmonisation.resource_preflight import validate_harmonisation_resource_maps
    def preflight(*a, **k):
        for resources in maps.values():
            resources["user_eaf_column"] = "AF"
        checked = validate_harmonisation_resource_maps(maps, grch_version="GRCh37", bcftools="unused",
            vcf_config=module.vcf_processing.model_dump(), default_eaf_colmap=MAPPING,
            external_eaf_colmap=MAPPING, external_info_colmap={}, policies=pol, require_default_eaf=True,
            resource_keys=("default_eaf_file",), validate_frequency_semantics=k.get("validate_frequency_semantics", True))
        return maps, checked
    monkeypatch.setattr(service, "preflight_harmonisation_resources", preflight)
    if mode == "mutate_during_staging":
        from postgwas.modules.harmonisation import external_reference_staging as staging
        batches = staging.projected_reference_batches
        def change_during_read(*a, **k):
            for batch in batches(*a, **k):
                yield batch
            path = Path(maps["1"]["user_eaf_file"])
            path.write_text(path.read_text() + "\n")
        monkeypatch.setattr(staging, "projected_reference_batches", change_during_read)
    if mode.startswith("mutate"):
        stage = service.stage_shared_external_reference_files
        def mutate(*a, **k):
            result = stage(*a, **k)
            path = Path(maps["1"]["default_eaf_file" if mode == "mutate_default_alias" else "user_eaf_file"])
            if mode in ("mutate_staged", "mutate_staged_alias"):
                staged = result[0]["1"]["user_eaf_file"]
                if mode == "mutate_staged_alias":
                    Path(staged).unlink()
                    Path(staged).symlink_to(path)
                else:
                    pl.read_parquet(staged).with_columns((pl.col("AF").cast(pl.Float64) + .005).alias("AF")).write_parquet(staged)
            elif mode in ("mutate_alias", "mutate_default_alias"):
                replacement = tmp_path / "different_original.tsv"
                replacement.write_bytes(path.read_bytes())
                path.unlink()
                path.symlink_to(replacement)
            else:
                path.write_text(path.read_text() + "\n")
            return result
        monkeypatch.setattr(service, "stage_shared_external_reference_files", mutate)
    monkeypatch.setattr(service, "_run_one_round", lambda *a, **k: pytest.fail("global failure launched workers/retry"))
    columns.update(sumstat_file=str(study), eaffile=maps["1"]["user_eaf_file"],
        se_col="SE", pval_col="P", declared_effect_type="beta", declared_pvalue_type="raw",
        ncontrol="1000", imp_info_col="INFO", info_source="internal", output_folder=str(tmp_path / "out"))
    with PipelineLogger("study", "dataset", str(tmp_path / "logs")) as logger:
        arguments = dict(sumstat_file=str(study), sample_column_dict=columns,
                output_dir=str(tmp_path / "out"), resource_folder=str(tmp_path),
                default_eaf_reference_source="fixture", default_comparison_af_file="fixture", dbsnp="fixture",
                resource_layout={}, output_layout=dict(module.output_layout.root), gwas2vcf_input=module.gwas2vcf_input.model_dump(),
                executables={"bcftools": "unused"}, vcf_config=module.vcf_processing.model_dump(),
                user_eaf_file=maps["1"]["user_eaf_file"], user_eaf_column="AF", user_info_file=None, user_info_column=None,
                default_eaf_reference_column="AF", default_comparison_af_column="AF", build_reference_files={},
                build_check_colmap={}, default_eaf_colmap=MAPPING, external_eaf_colmap=MAPPING, external_info_colmap={},
                external_reference_staging=module.external_reference_staging.model_dump(), policies=pol,
                logger=logger, input_line_count=source.height, max_workers=1, gwas2vcf_main_script_path="unused")
        if mode in ("passed", "same_passed"):
            arguments["output_dir_path"] = Path(arguments.pop("output_dir"))
            arguments.pop("max_workers")
            arguments.pop("gwas2vcf_main_script_path")
            result = service._prepare_dataset_for_chromosome_processing(**arguments)
            assert result.study_decisions["eaf_external_validation"]["status"] == "passed"
            assert result.dataset_study_decisions["eaf_external_validation"] == result.study_decisions["eaf_external_validation"]
            assert result.study_decisions["eaf_external_source"]["resolved_files"] == {
                c: maps[c]["user_eaf_file"] for c in maps}
            assert result.resource_preflight["external_reference_staging"]["source_files_read"] == 1
            assert all(result.resource_maps[c]["user_eaf_file"] != maps[c]["user_eaf_file"] for c in maps)
            return
        with pytest.raises(af.StudyAlleleFrequencyError, match="Spearman.*memory budget" if mode == "spearman_budget" else "external AF validation") as caught:
            service.harmonise_chromosomes(**arguments)
        if mode == "coverage":
            evidence = caught.value.study_decisions["eaf_external_validation"]
            assert evidence["study_variant_count"] == 2402 and evidence["comparable_variants"] == 1200
            assert evidence["decision_reason"] == "insufficient_overlap_fraction"
        assert caught.value.study_decisions["eaf_column"] is None
        evidence = caught.value.study_decisions["eaf_external_validation"]
        assert evidence["status"] == ("invalid" if mode.startswith("mutate") or mode == "spearman_budget" else "inconclusive")
        if mode == "spearman_budget":
            assert evidence["failure_reason"] == "spearman_memory_budget_exceeded"
            assert evidence["study_variant_count"] == evidence["comparable_variants"] == 1200
            assert evidence["raw_frequency_screen"]["usable"] == 1200 and "contract" not in evidence
            assert evidence["whole_study_rules_evaluated"] is False
        if mode.startswith("mutate"):
            assert evidence["failure_reason"] == ("external_staging_failed" if mode == "mutate_during_staging" else "external_inputs_changed") and "contract" not in evidence
        elif mode == "same_maf":
            assert evidence["decision_reason"] == "independent_reference_required"
        assert logger.summary()["steps"][-1]["status"] == "failed"
        assert "Whole-study external AF validation" in Path(logger.log_path).read_text()


def test_partition_metadata_cannot_shrink_original_fixed_denominator(tmp_path):
    parts = list(fixture(tmp_path, count=2400, matched=1200))
    parts[5] = ChromosomePartitions({"1": parts[5].paths["1"]}, {"1": 1200})
    parts[3].pop("2")
    with pytest.raises(af.StudyAlleleFrequencyError, match="reconcile"):
        validate(parts)


@pytest.mark.parametrize("case", ["ordinary", "direct_only", "subset", "source_name", "suffix", "missing", "invalid", "degenerate", "palindrome", "multiallelic", "unresolved", "unmatched", "default_multiallelic"])
def test_approved_assignment_matches_fresh_declared_alt_output_schema_and_all_row_qc(tmp_path, monkeypatch, case):
    def modify(source, default, external):
        if case in ("source_name", "suffix"):
            name = "__postgwas_reference_value__" if case == "source_name" else "AF"
            pl.read_csv(external, separator="\t").rename({"AF": name}).write_csv(external, separator="\t")
            return source.with_columns(pl.lit(.42).alias(name + "_" * i) for i in range(3))
        if case in ("missing", "invalid", "degenerate"):
            value = {"missing": None, "invalid": -1., "degenerate": 0.}[case]
            pl.read_csv(external, separator="\t").with_columns(
                pl.when(pl.col("POS") == 100).then(pl.lit(value)).otherwise(pl.col("AF")).alias("AF")
            ).write_csv(external, separator="\t")
        elif case == "palindrome":
            for path in (default, external):
                pl.read_csv(path, separator="\t").with_columns(pl.when(pl.col("POS") == 100)
                    .then(pl.lit("T")).otherwise(pl.col("ALT")).alias("ALT")).write_csv(path, separator="\t")
            source = source.with_columns(pl.when(pl.col("POS") == 100)
                .then(pl.lit("T")).otherwise(pl.col("EA")).alias("EA"))
        elif case in ("multiallelic", "unresolved", "default_multiallelic"):
            with (default if case == "default_multiallelic" else external).open("a") as handle:
                handle.write("1\t100\tA\tG\t%s\n" % ("NA" if case == "unresolved" else ".1"))
        elif case == "unmatched":
            pl.read_csv(default, separator="\t").filter(pl.col("POS") != 100).write_csv(default, separator="\t")
        return source
    parts = fixture(tmp_path, modify=modify, policies={"strand.unmatched_action": "retain",
        **({"eaf.match_orientations": ["direct"]} if case == "direct_only" else {})})
    if case in ("source_name", "suffix"):
        parts[1]["eafcolumn"] = "__postgwas_reference_value__" if case == "source_name" else "AF"
    decisions = validate(parts)
    evidence = decisions["eaf_external_validation"]
    assert evidence["comparable_variants"] == (1199 if case in (
        "missing", "invalid", "palindrome", "unmatched") else 1200)
    frame = pl.read_parquet(parts[5].paths["1"])
    if case == "subset":
        frame = frame.slice(1, 30).reverse()
    result, qc, columns = adopt(parts, decisions, frame=frame)
    def fresh(frame, cache, descriptor, keys, policies, **kwargs):
        return af.merge_external_allele_frequencies(frame, descriptor["statistics"]["source_file"],
            descriptor["statistics"]["source_column"], MAPPING, keys,
            chromosome=kwargs["chromosome"], policies=policies, logger=kwargs["logger"], ctx=kwargs["ctx"],
            rejects=kwargs["rejects"], study_columns_canonical=True, declared_alt_frequency=True)
    monkeypatch.setattr(af, "_consume_replacement_assignments", fresh)
    expected, expected_qc, expected_columns = adopt(parts, decisions, frame=frame)
    assert result.schema == expected.schema and result.equals(expected)
    assert columns == expected_columns
    ignored = {"external_merge_external_rows_loaded", "external_merge_reused_prepared_replacement_assignments", "external_merge_reused_prepared_external_assignments"}
    assert {k: v for k, v in qc.items() if k not in ignored} == {
        k: v for k, v in expected_qc.items() if k not in ignored}
    assert "MAF-like" not in str(qc.get("external_maf_reference_decision", ""))


@pytest.mark.parametrize("method", ["pearson", "spearman"])
def test_global_batch_metrics_equal_dataframe_and_obey_rank_budget(method, monkeypatch):
    from postgwas.modules.harmonisation.shared.runtime import NullStepContext
    frame = pl.DataFrame({"x": [.1, .2, .3, .4, .45, .48] * 200,
                          "r": [.1, .2, .3, .4, .45, .48] * 200})
    pol = default_policies().with_overrides({"eaf.maf_reference_correlation_method": method})
    expected = af._compare_frequency_pairs(frame, "x", "r", pol, NullStepContext(),
        accept_direct_eaf=True, study_variant_count=1200)
    batches = lambda: (frame.slice(i, 137) for i in range(0, 1200, 137))
    observed = af._compare_frequency_pairs(batches, "x", "r", pol, NullStepContext(),
        accept_direct_eaf=True, study_variant_count=1200, pair_count=1200)
    assert observed[0] == expected[0] == "eaf"
    for key in ("direct_eaf_correlation", "mean_eaf_error", "mean_maf_error", "reference_effect_allele_minor_fraction"):
        assert observed[1][key] == pytest.approx(expected[1][key])
    if method == "spearman":
        monkeypatch.setattr(af, "_measure_frequency_pairs", lambda *a, **k: pytest.fail("budget failure measured ranks"))
        with pytest.raises(af.StudyAlleleFrequencyError, match="memory budget"):
            af._compare_frequency_pairs(batches, "x", "r", pol.with_overrides({"execution.memory_budget_gb": .001}),
                NullStepContext(), accept_direct_eaf=True, study_variant_count=1200, pair_count=1200)


@pytest.mark.parametrize("passed", [False, True])
def test_external_terminal_and_html_report_recorded_global_evidence_without_reanalysis(tmp_path, monkeypatch, passed):
    from postgwas.modules.harmonisation.study_reporting import render_study_decisions, render_external_frequency_validation
    from postgwas.modules.harmonisation.html_report import render_dataset_report
    parts = fixture(tmp_path, values=None if passed else [.25] * 1200)
    if passed:
        decisions = validate(parts)
    else:
        with pytest.raises(af.StudyAlleleFrequencyError) as caught:
            validate(parts)
        decisions = caught.value.study_decisions
    monkeypatch.setattr(af, "_measure_frequency_pairs", lambda *a, **k: pytest.fail("report remeasured evidence"))
    reports = [render_study_decisions(decisions), render_external_frequency_validation(decisions),
        render_dataset_report({}, {"dataset": {"study_decisions": decisions}})]
    for report in reports:
        assert "External validation" in report or "External EAF approved" in report or "Inconclusive; dataset stopped" in report
        assert "1,200 / 1,200 (100.00%)" in report
        assert "Population-reference proxy" in report
        assert ("external EAF accepted for assignment" in report or "External EAF approved" in report) is passed
        assert "Replacement validation" not in report


def test_resolved_multiallelic_rows_keep_valid_private_diagnostic_frequencies(tmp_path, monkeypatch):
    def modify(source, default, external):
        with external.open("a") as handle:
            handle.write("1\t100\tA\tG\t.1\n")
        return source
    parts = fixture(tmp_path, modify=modify)
    writer = pl.DataFrame.write_parquet
    observed = []
    def inspect(frame, *a, **k):
        if frame.columns == ["study_frequency", "replacement_frequency", "eligible"]:
            observed.append(frame.clone())
            assert frame["eligible"].all()
            assert frame["study_frequency"].null_count() == 0
            assert frame["replacement_frequency"].null_count() == 0
        return writer(frame, *a, **k)
    monkeypatch.setattr(pl.DataFrame, "write_parquet", inspect)
    assert validate(parts)["eaf_external_validation"]["comparable_variants"] == 1200
    paired = pl.concat(observed)
    assert paired.height == 1200
    expected = [.1 + i % 40 / 100 for i in range(1200)]
    assert sorted(paired["study_frequency"].to_list()) == pytest.approx(sorted(expected))
    assert paired.select((pl.col("study_frequency") - pl.col("replacement_frequency")).abs().max()).item() == pytest.approx(0.)


def test_absolute_count_equality_pools_two_individually_insufficient_shards(tmp_path):
    parts = fixture(tmp_path, count=1000)
    assert list(parts[5].rows.values()) == [500, 500]
    evidence = validate(parts)["eaf_external_validation"]
    assert evidence["comparable_variants"] == evidence["minimum_overlap"] == 1000
    assert evidence["status"] == "passed"


@pytest.mark.parametrize("same", [False, True])
def test_non_maf_like_external_retains_declared_contract_without_comparison(tmp_path, monkeypatch, same):
    parts = fixture(tmp_path, same=same, values=[.1, .2, .3, .7, .8, .9] * 200)
    monkeypatch.setattr(af, "_compare_frequency_pairs", lambda *a, **k: pytest.fail("non-MAF-like source was reclassified"))
    decisions = validate(parts)
    evidence = decisions["eaf_external_validation"]
    assert evidence["status"] == "passed" and evidence["raw_frequency_screen"]["maf_like"] is False
    assert "rule_1_passed" not in evidence
    output, qc, columns = adopt(parts, decisions)
    assert output.height == 600 and qc["eaf_provenance"] == "reference_imputed"
    # Compare the existing direct-library declared source path, including the
    # identical-source retained column and effect/row/output contracts.
    original = {"effect_type": "beta", "strand": "forward", "eaf_is_maf": None}
    expected, _, expected_columns = adopt(parts, original, frame=parts[0].filter(pl.col("CHR") == "1"))
    assert output.schema == expected.schema and output.equals(expected)
    assert columns == expected_columns


@pytest.mark.parametrize("staged", [False, True])
def test_identical_declared_source_reuses_default_preparation_without_external_parse_or_join(tmp_path, monkeypatch, staged):
    parts = fixture(tmp_path, same=True, values=[.1, .2, .3, .7, .8, .9] * 200)
    originals = deepcopy(parts[3])
    source_bindings = af.capture_external_frequency_sources(originals)
    staging_summary = None
    if staged:
        from postgwas.modules.harmonisation.external_reference_staging import stage_shared_external_reference_files
        module = load_configuration().modules.harmonisation
        staged_maps, staging_summary = stage_shared_external_reference_files(parts[3], chromosomes=["1", "2"],
            dataset_id="study", output_directory=tmp_path / "staged", output_layout=dict(module.output_layout.root),
            user_eaf_specification=originals["1"]["user_eaf_file"], user_eaf_column="AF", external_eaf_mapping=MAPPING,
            user_info_specification=None, user_info_column=None, external_info_mapping={},
            settings=module.external_reference_staging.model_dump(), policies=parts[4])
        parts[3].update(staged_maps)
    monkeypatch.setattr(af, "read_reference_variant_table", lambda *a, **k: pytest.fail("identical source parsed twice"))
    monkeypatch.setattr(af, "allele_oriented_left_join", lambda *a, **k: pytest.fail("identical source rematched"))
    decisions = af.validate_external_allele_frequency(parts[5], parts[1], parts[2], parts[3], MAPPING,
        MAPPING, original_resource_maps=originals, policies=parts[4],
        source_bindings=source_bindings, staging_summary=staging_summary)
    assert decisions["eaf_external_validation"]["raw_frequency_screen"]["usable"] == 1200
    for c in parts[3]:
        result, qc, columns = adopt(parts, decisions, c)
        assert result.height == 600 and qc["external_merge_reused_aligned_strand_reference"] is True
        assert columns["eaf_col"] == ("AF" if staged else af.REFERENCE_AF_COLUMN)
        assert "assignment_cache" not in decisions["eaf_external_validation"]["contract"]["chromosomes"][c]


def test_old_directional_margin_failure_becomes_pooled_global_eaf_acceptance(tmp_path, monkeypatch):
    # Per 1000 rows: raw 96% minor, default counted-allele 94% minor;
    # 20 near-half and 40 major pairs each differ by .02. Direct MAE=.0012,
    # raw-AF-versus-reference-MAF MAE=.0112, so the .01 error gap fails the legacy .02 margin.
    block_x = [.1 + i % 31 / 100 for i in range(940)] + [.49] * 20 + [.63] * 40
    block_r = block_x[:940] + [.51] * 20 + [.65] * 40
    x, r = block_x * 4, block_r * 4
    parts = fixture(tmp_path, count=4000, values=x, default_values=r)
    std = {key: parts[1][key + "_col"] for key in ("chr", "pos", "ea", "oa")}
    legacy_frame = parts[0].with_columns(pl.Series("external", x), pl.Series("default", r))
    for c in ("1", "2"):
        decision, stats = af.confirm_eaf_with_reference(legacy_frame.filter(pl.col("CHR") == c),
            "external", parts[3][c]["default_eaf_file"], "AF", MAPPING, std,
            policies=parts[4], aligned_reference_col="default")
        assert decision == "inconclusive" and stats["decision_reason"] == "error_difference_below_required_margin"
    decisions = validate(parts)
    evidence = decisions["eaf_external_validation"]
    assert evidence["status"] == "passed" and evidence["rule_1_passed"] and not evidence["rule_2_passed"]
    assert evidence["study_variant_count"] == evidence["comparable_variants"] == 4000
    assert evidence["mean_eaf_error"] == pytest.approx(.0012)
    assert evidence["mean_maf_error"] == pytest.approx(.0112)
    assert evidence["reference_effect_allele_minor_fraction"] == pytest.approx(.94)
    monkeypatch.setattr(af, "_compare_frequency_pairs", lambda *a, **k: pytest.fail("approved worker used legacy classifier"))
    monkeypatch.setattr(af, "read_reference_variant_table", lambda *a, **k: pytest.fail("approved worker reread source"))
    for c in ("1", "2"):
        result, qc, columns = adopt(parts, decisions, c)
        assert result.height == 2000
        assert result[columns["eaf_col"]].to_list() == pytest.approx(x[0 if c == "1" else 1::2])
        assert result["BETA"].to_list() == [.2] * 2000


@pytest.mark.parametrize("visible", [False, True])
def test_initial_external_approval_survives_canonical_visible_and_hidden_screen_recording(tmp_path, visible):
    import json
    import os
    import subprocess
    import sys
    import textwrap
    decisions = validate(fixture(tmp_path / "inputs"))
    (tmp_path / "decisions.json").write_text(json.dumps(decisions))
    script = textwrap.dedent('''
        import json
        from pathlib import Path
        import sys
        from postgwas.core.pipeline_logging import PipelineLogger
        from postgwas.core.screen_logging import ScreenSettings, record_screen
        from postgwas.modules.harmonisation.study_reporting import render_external_frequency_validation
        from postgwas.modules.harmonisation.service import _announce
        output = Path(sys.argv[1])
        decisions = json.loads((output / "decisions.json").read_text())
        with PipelineLogger("study", "dataset", str(output / "logs"), screen_level="ERROR") as logger:
            with record_screen(ScreenSettings(sys.argv[2] == "True", output / "screen.log")):
                _announce(logger, render_external_frequency_validation(decisions))
    ''')
    completed = subprocess.run([sys.executable, "-B", "-c", script, str(tmp_path), str(visible)],
        env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1", POLARS_MAX_THREADS="4"),
        capture_output=True, text=True, timeout=60, check=False)
    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    assert ("External EAF approved" in completed.stdout) is visible
    for path in (tmp_path / "screen.log", tmp_path / "logs" / "study_dataset.log"):
        from postgwas.core.pipeline_logging import _PREFIX_WIDTH
        text = path.read_text()
        if path.name == "study_dataset.log":
            text = "\n".join(line[_PREFIX_WIDTH:] for line in text.splitlines())
        normalized = " ".join(text.split())
        assert "External EAF approved" in normalized
        assert "Whole-study external AF validation" in normalized and "\x1b" not in text


def test_spearman_resource_failure_preserves_measured_external_evidence_and_reports(tmp_path):
    from postgwas.modules.harmonisation.study_reporting import external_frequency_fields, render_external_frequency_validation
    from postgwas.modules.harmonisation.html_report import render_dataset_report
    parts = fixture(tmp_path, policies={"eaf.maf_reference_correlation_method": "spearman",
        "execution.memory_budget_gb": .001})
    before = {p: Path(p).read_bytes() for p in parts[5].paths.values()}
    with pytest.raises(af.StudyAlleleFrequencyError, match="Spearman.*memory budget") as caught:
        validate(parts)
    decisions = caught.value.study_decisions
    evidence = decisions["eaf_external_validation"]
    assert evidence["status"] == "invalid" and "contract" not in evidence
    assert evidence["failure_reason"] == "spearman_memory_budget_exceeded"
    assert evidence["decision_reason"] == "spearman_memory_budget_exceeded"
    assert evidence["study_variant_count"] == evidence["comparable_variants"] == 1200
    assert evidence["observed_overlap_fraction"] == 1.
    assert evidence["minimum_overlap"] == 1000 and evidence["minimum_overlap_fraction"] == .5
    assert evidence["raw_frequency_screen"]["usable"] == evidence["raw_frequency_screen"]["at_or_below_0_5"] == 1200
    assert evidence["source"]["column"] == "AF" and set(evidence["source"]["identities"]) == {"1", "2"}
    assert evidence["source"]["mapping"] == MAPPING
    assert evidence["maf_correlation_method"] == "spearman"
    assert evidence["shared_evidence_gate_passed"] is True
    assert evidence["whole_study_rules_evaluated"] is False
    assert evidence["estimated_rank_memory_bytes"] > 0
    assert evidence["available_rank_memory_bytes"] < evidence["estimated_rank_memory_bytes"]
    assert evidence["direct_eaf_correlation"] is None and evidence["mean_eaf_error"] is None
    assert "memory budget" in evidence["failure_message"]
    fields = {label: value for _, label, value in external_frequency_fields(decisions)}
    assert fields["Study rows entering AF preflight (fixed N)"] == "1,200"
    assert fields["Raw external AF at or below 0.5"] == "1,200 / 1,200 (100.00%)"
    assert fields["External diagnostic coverage"] == "1,200 / 1,200 (100.00%); required >= 50.00%"
    assert fields["External shared evidence gate"] == "Passed"
    assert fields["External rule 1: Retain declared external EAF"] == "Not evaluated"
    for report in (render_external_frequency_validation(decisions),
        render_dataset_report({}, {"dataset": {"study_decisions": decisions}})):
        assert "Invalid; dataset stopped before chromosome analysis" in report
        normalized = " ".join(report.split())
        assert "estimated rank memory" in normalized.lower() and "available rank memory" in normalized.lower()
        assert "external EAF accepted for assignment" not in report
    assert {p: Path(p).read_bytes() for p in before} == before


@pytest.mark.parametrize("moment", ["before_reference", "during_measurement"])
def test_typed_initial_preparation_failure_finalizes_known_evidence(tmp_path, monkeypatch, moment):
    parts = fixture(tmp_path)
    def fail(*a, **k):
        raise af.StudyAlleleFrequencyError("synthetic typed preparation failure")
    monkeypatch.setattr(af, "prepared_work_fingerprint" if moment == "before_reference" else "_measure_frequency_pairs", fail)
    with pytest.raises(af.StudyAlleleFrequencyError, match="synthetic typed preparation failure") as caught:
        validate(parts)
    evidence = caught.value.study_decisions["eaf_external_validation"]
    assert evidence["status"] == "invalid" and "contract" not in evidence
    assert evidence["failure_reason"] == "external_preparation_failed"
    assert evidence["study_variant_count"] == 1200
    assert evidence["source"]["column"] == "AF" and set(evidence["source"]["identities"]) == {"1", "2"}
    assert evidence["failure_message"] == "synthetic typed preparation failure"
    if moment == "during_measurement":
        assert evidence["raw_frequency_screen"]["usable"] == 1200
        assert evidence["comparable_variants"] == 1200 and evidence["observed_overlap_fraction"] == 1.
        assert evidence["direct_eaf_correlation"] is None
        assert evidence["shared_evidence_gate_passed"] is True and evidence["whole_study_rules_evaluated"] is False


@pytest.mark.parametrize("same", [False, True])
def test_previous_chromosome_packed_proof_is_released_before_next_preparation(tmp_path, monkeypatch, same):
    import gc
    import weakref
    parts = fixture(tmp_path, same=same, values=[.1, .2, .3, .7, .8, .9] * 200)
    prepare = af.prepare_study_reference_matches
    references = []
    def inspect(*a, **k):
        gc.collect()
        assert all(reference() is None for reference in references), "previous chromosome packed proof remains live"
        result = prepare(*a, **k)
        references.append(weakref.ref(result[1]))
        return result
    monkeypatch.setattr(af, "prepare_study_reference_matches", inspect)
    evidence = validate(parts)["eaf_external_validation"]
    assert evidence["status"] == "passed" and len(references) == 2
    gc.collect()
    assert all(reference() is None for reference in references)
    for c, path in parts[5].paths.items():
        assert pl.read_parquet(path).select(parts[0].columns).equals(parts[0].filter(pl.col("CHR") == c))


def test_typed_second_shard_failure_retains_partial_raw_counts_without_claiming_global_screen(tmp_path, monkeypatch):
    from postgwas.modules.harmonisation.study_reporting import external_frequency_fields
    parts = fixture(tmp_path)
    prepare = af.prepare_study_reference_matches
    calls = []
    def fail_second(*a, **k):
        calls.append(1)
        if len(calls) == 2:
            raise af.StudyAlleleFrequencyError("second chromosome reference preparation failed")
        return prepare(*a, **k)
    monkeypatch.setattr(af, "prepare_study_reference_matches", fail_second)
    with pytest.raises(af.StudyAlleleFrequencyError, match="second chromosome") as caught:
        validate(parts)
    evidence = caught.value.study_decisions["eaf_external_validation"]
    assert evidence["status"] == "invalid" and "contract" not in evidence
    assert evidence["study_variant_count"] == 1200
    assert "raw_frequency_screen" not in evidence and "comparable_variants" not in evidence
    partial = evidence["raw_frequency_screens_by_chromosome"]
    assert set(partial) == {"1"}
    assert partial["1"]["usable"] == partial["1"]["at_or_below_0_5"] == 600
    assert evidence["excluded_counts"]["original_orientation_ineligible_rows"] == 0
    report = {label: value for _, label, value in external_frequency_fields(caught.value.study_decisions)}
    assert "do not decide the GWAS-row screen" in report["Partial source-panel description"]
    assert "1 chromosome(s) described" in report["Partial source-panel description"]
    assert not any(label.startswith("Chromosome ") for label in report)
    assert "effective_frequency_screen" not in evidence
    effective = evidence["effective_frequency_screens_by_chromosome"]
    assert set(effective) == {"1"}
    assert effective["1"]["usable"] == effective["1"]["at_or_below_0_5"] == 600
