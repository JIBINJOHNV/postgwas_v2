"""The selected AF source must pass a row-bound whole-study gate before adoption."""
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


def test_prepared_work_fingerprint_preserves_each_callers_error_type(tmp_path):
    from postgwas.core.input_validation import InputValidationSession
    from postgwas.modules.harmonisation import service
    from postgwas.modules.harmonisation.shared import runtime
    path = tmp_path / "work.parquet"
    with InputValidationSession(pipeline=False):
        for error_type in (service.PipelineError, af.StudyAlleleFrequencyError):
            with pytest.raises(error_type, match="work.parquet"):
                runtime.prepared_work_fingerprint(path, error_type=error_type)
        path.write_bytes(b"validated work content")
        first = runtime.prepared_work_fingerprint(path, error_type=service.PipelineError)
        assert runtime.prepared_work_fingerprint(path, error_type=af.StudyAlleleFrequencyError) == first
        path.unlink()
        with pytest.raises(af.StudyAlleleFrequencyError, match="work.parquet"):
            runtime.prepared_work_fingerprint(path, error_type=af.StudyAlleleFrequencyError)


def fixture(tmp_path, *, same=False, chromosomes=("1",), replacement_values=None, modify=None, policy_overrides=None, frequency_column="MAF"):
    tmp_path.mkdir(parents=True, exist_ok=True)
    original_b = [.8, .7, .6, .1, .2, .3] * 2
    alleles = [("C", "A"), ("A", "C"), ("G", "T"), ("T", "G")] * 3
    chrom = [chromosomes[i % len(chromosomes)] for i in range(12)]
    reference_b = [1 - b if i % 4 in (1, 3) else b for i, b in enumerate(original_b)]
    source = pl.DataFrame({"CHR": chrom, "POS": list(range(100, 112)),
        "EA": [x[0] for x in alleles], "OA": [x[1] for x in alleles],
        "MAF": [min(b, 1-b) for b in original_b], "BETA": [.2] * 12, "Z": [2.] * 12,
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN).rename({"MAF": frequency_column})
    columns = {"chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA",
        "eaf_col": frequency_column, "beta_or_col": "BETA", "imp_z_col": "Z",
        "gwas_outputname": "study", "output_folder": str(tmp_path / "out")}
    reference = pl.DataFrame({"CHROM": chrom, "POS": list(range(100,112)),
        "REF": ["A"] * 12, "ALT": ["C"] * 12, "AF": reference_b})
    default = tmp_path / "default.tsv"
    reference.write_csv(default, separator="\t")
    replacement = default if same else tmp_path / "replacement.tsv"
    if not same:
        reference.with_columns(pl.Series("AF", reference_b if replacement_values is None else replacement_values)).write_csv(replacement, separator="\t")
    if modify:
        source = modify(source, default, replacement)
    policies = default_policies().with_overrides({"eaf.maf_reference_min_overlap": 2, **(policy_overrides or {})})
    maps = {c: {"default_eaf_file": str(default), "default_eaf_reference_column": "AF",
        "user_eaf_file": None if same else str(replacement)} for c in chromosomes}
    decisions, work = af.validate_study_allele_frequency(source, columns,
        {"eaf_column": frequency_column, "eaf_is_maf": True, "effect_type": "beta", "strand": "forward"},
        maps, MAPPING, policies=policies, eaf_replacement_file=str(replacement), eaf_replacement_column="AF")
    decisions["eaf_replacement"]["resolved_files"] = {c: str(replacement) for c in chromosomes}
    paths, rows = {}, {}
    for c in chromosomes:
        path = tmp_path / ("work_" + c + ".parquet")
        subset = work.filter(pl.col("CHR") == c)
        subset.write_parquet(path)
        paths[c], rows[c] = str(path), subset.height
    return source, columns, decisions, maps, policies, ChromosomePartitions(paths, rows), original_b


def validate(parts):
    _, columns, decisions, maps, policies, partitions, _ = parts
    return af.validate_replacement_allele_frequency(partitions, columns, decisions, maps, MAPPING, policies=policies)


def adopt(parts, decisions, frame=None, *, mapping=MAPPING, policies=None):
    _, columns, _, maps, default_pol, partitions, _ = parts
    c = next(iter(partitions.paths))
    return af.harmonise_allele_frequency(c, pl.read_parquet(partitions.paths[c]) if frame is None else frame,
        dict(columns), dict(load_configuration().modules.harmonisation.output_layout.root), "\t",
        eaffile=maps[c]["user_eaf_file"], external_eaf_colmap=mapping,
        default_eaf_file=maps[c]["default_eaf_file"], default_eaf_column="AF", default_eaf_colmap=MAPPING,
        policies=default_pol if policies is None else policies, study_decision=decisions)


def fresh_assignments(frame, cache, descriptor, keys, policies, **kwargs):
    return af.merge_external_allele_frequencies(frame, descriptor["statistics"]["source_file"],
        descriptor["statistics"]["source_column"], MAPPING, keys,
        chromosome=kwargs["chromosome"], policies=policies, logger=kwargs["logger"], ctx=kwargs["ctx"],
        rejects=kwargs["rejects"], study_columns_canonical=True, declared_alt_frequency=True)


def test_nonselected_replacement_gate_never_reads_partitions_or_candidate(monkeypatch):
    monkeypatch.setattr(pl, "read_parquet", lambda *a, **k: pytest.fail("candidate was opened"))
    original = {"eaf_selected_source": "study_supplied"}
    assert af.validate_replacement_allele_frequency(None, {}, original, {}, {}) == original


@pytest.mark.parametrize("same", [False, True])
@pytest.mark.parametrize("chromosomes", [("1",), ("1", "2")])
def test_gate_pools_original_ea_orientation_without_mutation_or_repeated_merge(tmp_path, monkeypatch, same, chromosomes):
    parts = fixture(tmp_path, same=same, chromosomes=chromosomes)
    source, columns, original, maps, policies, partitions, original_b = parts
    before = {p: pl.read_parquet(p) for p in partitions.paths.values()}
    calls = []
    merge = af.merge_external_allele_frequencies
    def observed(*args, **kwargs):
        calls.append(kwargs)
        return merge(*args, **kwargs)
    monkeypatch.setattr(af, "merge_external_allele_frequencies", observed)
    result = validate(parts)
    evidence = result["eaf_replacement_validation"]
    assert evidence["status"] == "passed"
    assert evidence["study_variant_count"] == evidence["comparable_variants"] == 12
    assert evidence["major_count"] == 6
    assert evidence["source"]["same_source_evidence_reused"] is same
    assert "evidence_is_independent_of_default" not in evidence["source"]
    assert evidence["correlation"] == pytest.approx(1.)
    assert evidence["mean_absolute_difference"] == pytest.approx(0.)
    assert result["eaf_reference_validation"] == original["eaf_reference_validation"]
    assert "eaf_replacement_validation" not in original
    assert len(calls) == (0 if same else len(chromosomes))
    assert all(call["min_match_fraction"] == 0.0 and call["declared_alt_frequency"] for call in calls)
    for path, original_frame in before.items():
        published = pl.read_parquet(path)
        assert published.select(original_frame.columns).equals(original_frame)
        assert published.width == original_frame.width + (0 if same else 1)
    assert set(tmp_path.iterdir()) == {Path(p) for p in partitions.paths.values()} | {tmp_path / "default.tsv"} | (set() if same else {tmp_path / "replacement.tsv"})


def test_failed_gate_carries_original_and_new_evidence(tmp_path):
    parts = fixture(tmp_path, replacement_values=[.25] * 12)
    with pytest.raises(af.StudyAlleleFrequencyError, match="replacement.*inconclusive") as caught:
        validate(parts)
    decisions = caught.value.study_decisions
    assert decisions["eaf_reference_validation"] == parts[2]["eaf_reference_validation"]
    assert decisions["eaf_replacement_validation"]["status"] == "inconclusive"
    assert decisions["eaf_replacement_validation"]["failed_requirements"]


def test_provisional_or_status_only_approval_cannot_adopt(tmp_path):
    parts = fixture(tmp_path)
    for approval in (None, {"status": "passed"}):
        decisions = deepcopy(parts[2])
        if approval:
            decisions["eaf_replacement_validation"] = approval
        with pytest.raises(af.AlleleFrequencyError, match="replacement.*validation"):
            adopt(parts, decisions)


def test_passed_adoption_accepts_validated_subset_and_preserves_original_maf(tmp_path):
    parts = fixture(tmp_path)
    decisions = validate(parts)
    frame = pl.read_parquet(parts[5].paths["1"]).slice(2, 7)
    result, qc, columns = adopt(parts, decisions, frame)
    assert result["MAF"].to_list() == frame["MAF"].to_list()
    assert result[columns["eaf_col"]].to_list() == pytest.approx([
        1-b if i % 4 in (1,3) else b for i,b in enumerate(parts[6])][2:9])
    assert qc["eaf_replacement_validation"]["status"] == "passed"


@pytest.mark.parametrize("change", ["frequency", "row_id", "alleles", "policy", "mapping", "source", "partition"])
def test_adoption_rejects_stale_or_changed_validation_contract(tmp_path, change):
    parts = fixture(tmp_path)
    decisions = validate(parts)
    frame = pl.read_parquet(parts[5].paths["1"])
    kwargs = {}
    if change == "frequency":
        frame = frame.with_columns(pl.lit(.123).alias("MAF"))
    elif change == "row_id":
        frame = frame.with_columns((pl.col(SOURCE_INPUT_ROW_COLUMN) + 100).alias(SOURCE_INPUT_ROW_COLUMN))
    elif change == "alleles":
        frame = frame.with_columns(pl.lit("T").alias("EA"))
    elif change == "policy":
        kwargs["policies"] = parts[4].with_overrides({"eaf.study_reference_min_correlation": .99})
    elif change == "mapping":
        kwargs["mapping"] = dict(MAPPING, a1="REF", a2="ALT")
    elif change == "source":
        path = tmp_path / "replacement.tsv"
        path.write_text(path.read_text().replace("0.8", "0.81"))
    else:
        frame.with_columns(pl.lit(.123).alias("MAF")).write_parquet(parts[5].paths["1"])
    with pytest.raises(af.AlleleFrequencyError, match="replacement.*validation|replacement.*changed"):
        adopt(parts, decisions, frame, **kwargs)


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), -1.])
def test_unchanged_invalid_original_values_remain_provenance(tmp_path, value):
    def modify(source, default, replacement):
        return source.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit(value, dtype=pl.Float64)).otherwise(pl.col("MAF")).alias("MAF"))
    parts = fixture(tmp_path, modify=modify)
    decisions = validate(parts)
    assert decisions["eaf_replacement_validation"]["comparable_variants"] == 11
    frame = pl.read_parquet(parts[5].paths["1"]).reverse()
    result, _, columns = adopt(parts, decisions, frame)
    assert result["MAF"].equals(frame["MAF"])
    assert result.height == 12 and result[columns["eaf_col"]].null_count() == 0


@pytest.mark.parametrize("extra", [.1, None, float("nan"), .9])
@pytest.mark.parametrize("swapped", [False, True])
def test_multiallelic_diagnostics_exclude_rows_without_changing_attachment(tmp_path, extra, swapped):
    def modify(source, default, replacement):
        ref = pl.read_csv(replacement, separator="\t", infer_schema_length=0)
        if swapped:
            ref = ref.with_columns(pl.when(pl.col("POS") == "100").then(pl.lit("C")).otherwise(pl.col("REF")).alias("REF"),
                pl.when(pl.col("POS") == "100").then(pl.lit("A")).otherwise(pl.col("ALT")).alias("ALT"),
                pl.when(pl.col("POS") == "100").then(pl.lit("0.2")).otherwise(pl.col("AF")).alias("AF"))
        ref.write_csv(replacement, separator="\t")
        with replacement.open("a") as handle:
            handle.write("1\t100\t%s\tG\t%s\n" % ("C" if swapped else "A", "NA" if extra is None else extra))
        return source
    parts = fixture(tmp_path, modify=modify)
    decisions = validate(parts)
    evidence = decisions["eaf_replacement_validation"]
    assert evidence["study_variant_count"] == 12 and evidence["comparable_variants"] == 11
    assert evidence["excluded_counts"]["replacement_non_biallelic_or_unresolved_rows"] == 1
    result, _, columns = adopt(parts, decisions)
    expected = (None if extra is None or extra != extra or extra == .9 else .7) if swapped else .8
    if expected is None:
        assert result[columns["eaf_col"]][0] is None
    else:
        assert result[columns["eaf_col"]][0] == pytest.approx(expected)
    assert result.height == 12 and result["MAF"].equals(parts[0]["MAF"])


def test_global_half_coverage_does_not_override_worker_eighty_percent_gate(tmp_path):
    def modify(source, default, replacement):
        pl.read_csv(replacement, separator="\t").head(6).write_csv(replacement, separator="\t")
        return source
    parts = fixture(tmp_path, modify=modify)
    decisions = validate(parts)
    assert decisions["eaf_replacement_validation"]["observed_overlap_fraction"] == .5
    with pytest.raises(af.AlleleFrequencyError, match="external_min_match_fraction"):
        adopt(parts, decisions)


def test_palindromic_and_unmatched_original_orientation_remain_in_fixed_n(tmp_path):
    def modify(source, default, replacement):
        # Two rows have no diagnostic orientation; other rows still support R2.
        source = source.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit("T")).otherwise(pl.col("EA")).alias("EA"))
        reference = pl.read_csv(default, separator="\t").filter(pl.col("POS") != 101)
        reference.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit("T")).otherwise(pl.col("ALT")).alias("ALT")).write_csv(default, separator="\t")
        return source
    parts = fixture(tmp_path, modify=modify)
    decisions = validate(parts)
    evidence = decisions["eaf_replacement_validation"]
    assert evidence["study_variant_count"] == 12 and evidence["comparable_variants"] == 10
    assert evidence["excluded_counts"]["original_orientation_ineligible_rows"] == 2
    assert pl.read_parquet(parts[5].paths["1"]).height == 12


def test_indel_diagnostic_keys_share_default_reference_canonicalization(tmp_path):
    def modify(source, default, replacement):
        for path in (default, replacement):
            ref = pl.read_csv(path, separator="\t", infer_schema_length=0)
            ref.with_columns(pl.when(pl.col("POS") == "100").then(pl.lit("aac")).otherwise(pl.col("ALT")).alias("ALT"),
                pl.when(pl.col("POS") == "100").then(pl.lit("aa")).otherwise(pl.col("REF")).alias("REF")).write_csv(path, separator="\t")
        return source.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit(" AAC ")).otherwise(pl.col("EA")).alias("EA"),
            pl.when(pl.col("POS") == 100).then(pl.lit(" aa ")).otherwise(pl.col("OA")).alias("OA"))
    parts = fixture(tmp_path, modify=modify)
    result = validate(parts)
    assert result["eaf_replacement_validation"]["comparable_variants"] == 12
    assert parts[0]["EA"][0] == " AAC "


def test_gate_refuses_source_changed_during_numeric_replay(tmp_path, monkeypatch):
    parts = fixture(tmp_path)
    original = af.validate_replacement_frequency_pairs
    def mutate(*args, **kwargs):
        result = original(*args, **kwargs)
        path = tmp_path / "replacement.tsv"
        path.write_text(path.read_text().replace("0.8", "0.81"))
        return result
    monkeypatch.setattr(af, "validate_replacement_frequency_pairs", mutate)
    with pytest.raises(af.StudyAlleleFrequencyError, match="changed") as caught:
        validate(parts)
    evidence = caught.value.study_decisions["eaf_replacement_validation"]
    assert evidence["status"] == "invalid" and "contract" not in evidence
    assert evidence["comparable_variants"] == 12


def test_distinct_partition_directories_can_share_a_basename(tmp_path):
    parts = list(fixture(tmp_path, chromosomes=("1", "2")))
    paths = {}
    for chromosome, original in parts[5].paths.items():
        destination = tmp_path / chromosome / "work.parquet"
        destination.parent.mkdir()
        Path(original).rename(destination)
        paths[chromosome] = str(destination)
    parts[5] = ChromosomePartitions(paths, parts[5].rows)
    assert validate(parts)["eaf_replacement_validation"]["status"] == "passed"
    assert list((tmp_path / "1").iterdir()) == [Path(paths["1"])]


def test_inconclusive_evidence_cannot_be_approved_by_changing_status(tmp_path):
    parts = fixture(tmp_path, replacement_values=[.25] * 12)
    with pytest.raises(af.StudyAlleleFrequencyError) as caught:
        validate(parts)
    decisions = caught.value.study_decisions
    decisions["eaf_replacement_validation"]["status"] = "passed"
    with pytest.raises(af.AlleleFrequencyError, match="replacement validation"):
        adopt(parts, decisions)


@pytest.mark.parametrize("alias_kind", ["symlink", "hardlink"])
def test_default_source_alias_reuses_prepared_af_without_source_read(tmp_path, monkeypatch, alias_kind):
    alias = tmp_path / "alias.tsv"
    def modify(source, default, replacement):
        getattr(alias, "symlink_to" if alias_kind == "symlink" else "hardlink_to")(default)
        return source
    parts = fixture(tmp_path, same=True, modify=modify)
    parts[2]["eaf_replacement"].update(file=str(alias), resolved_files={"1": str(alias)})
    monkeypatch.setattr(af, "read_reference_variant_table", lambda *a, **k: pytest.fail("identical source reread"))
    decisions = validate(parts)
    assert decisions["eaf_replacement_validation"]["source"]["reused_default_chromosomes"] == ["1"]
    assert adopt(parts, decisions)[0].height == 12


@pytest.mark.parametrize("route", ["direct", "activated", "staged"])
@pytest.mark.parametrize("when", ["after_measurement", "after_approval", "during_strand"])
def test_retargeted_alias_never_adopts_unvalidated_frequencies(tmp_path, monkeypatch, route, when):
    from postgwas.modules.harmonisation import service
    parts = list(fixture(tmp_path, same=route != "staged"))
    alias = tmp_path / "alias_GRCh37_1.tsv"
    original = tmp_path / ("replacement.tsv" if route == "staged" else "default.tsv")
    alias.symlink_to(original)
    alternate = tmp_path / "unvalidated.tsv"
    pl.read_csv(original, separator="\t").with_columns(pl.lit(.25).alias("AF")).write_csv(alternate, separator="\t")
    if route == "staged":
        parts[3]["1"]["user_eaf_file"] = str(alias)
    else:
        parts[2]["eaf_replacement"].update(file=str(alias), resolved_files={"1": str(alias)})
    if route == "activated":
        parts[2]["eaf_replacement"]["file"] = str(tmp_path / "alias_{build}_{chromosome}.tsv")
        parts[3], _ = service._activate_eaf_replacement_resources(
            parts[3], parts[2], grch_version="GRCh37", default_eaf_colmap=MAPPING,
            external_eaf_colmap=MAPPING, executables={"bcftools": "must-not-run"},
            vcf_config=load_configuration().modules.harmonisation.vcf_processing.model_dump(), policies=parts[4])
    def retarget():
        alias.unlink()
        alias.symlink_to(alternate)
    if when == "after_measurement":
        measure = af.validate_replacement_frequency_pairs
        def measured(*args, **kwargs):
            result = measure(*args, **kwargs)
            retarget()
            return result
        monkeypatch.setattr(af, "validate_replacement_frequency_pairs", measured)
        with pytest.raises(af.StudyAlleleFrequencyError, match="changed") as caught:
            validate(parts)
        evidence = caught.value.study_decisions["eaf_replacement_validation"]
        assert evidence["status"] == "invalid" and "contract" not in evidence
        assert evidence["comparable_variants"] == 12 and evidence["correlation"] == pytest.approx(1.)
        return
    decisions = validate(parts)
    if when == "after_approval":
        retarget()
        monkeypatch.setattr(af, "harmonise_strand_orientation", lambda *a, **k: pytest.fail("strand ran after alias retarget"))
        with pytest.raises(af.AlleleFrequencyError, match="changed"):
            adopt(parts, decisions)
        return
    orient = af.harmonise_strand_orientation
    def oriented(*args, **kwargs):
        result = orient(*args, **kwargs)
        retarget()
        return result
    monkeypatch.setattr(af, "harmonise_strand_orientation", oriented)
    # Once verified, attachment must use the pinned resolved source, never the
    # new target of a mutable alias passed between orientation and attachment.
    result, _, selected = adopt(parts, decisions)
    expected = pl.read_csv(original, separator="\t")["AF"].to_list()
    assert result[selected["eaf_col"]].to_list() == pytest.approx(expected)


def test_changed_frequency_column_same_file_runs_shared_merge(tmp_path, monkeypatch):
    def modify(source, default, replacement):
        pl.read_csv(default, separator="\t").with_columns(pl.col("AF").alias("AF2")).write_csv(default, separator="\t")
        return source
    parts = fixture(tmp_path, same=True, modify=modify)
    parts[2]["eaf_replacement"]["column"] = "AF2"
    calls = []
    original = af.merge_external_allele_frequencies
    def merge(*args, **kwargs):
        calls.append(args[2])
        return original(*args, **kwargs)
    monkeypatch.setattr(af, "merge_external_allele_frequencies", merge)
    assert validate(parts)["eaf_replacement_validation"]["status"] == "passed"
    assert calls == ["AF2"]


def test_changed_allele_mapping_same_file_runs_shared_merge(tmp_path, monkeypatch):
    def modify(source, default, replacement):
        pl.read_csv(default, separator="\t").with_columns(pl.col("ALT").alias("ALT2")).write_csv(default, separator="\t")
        return source
    parts = fixture(tmp_path, same=True, modify=modify)
    calls = []
    original = af.merge_external_allele_frequencies
    def merge(*args, **kwargs):
        calls.append(args[3])
        return original(*args, **kwargs)
    mapping = dict(MAPPING, a1="ALT2")
    monkeypatch.setattr(af, "merge_external_allele_frequencies", merge)
    decisions = af.validate_replacement_allele_frequency(parts[5], parts[1], parts[2], parts[3], mapping, policies=parts[4])
    assert decisions["eaf_replacement_validation"]["status"] == "passed"
    assert calls == [mapping]


def test_unreadable_candidate_preserves_io_failure_and_original_evidence(tmp_path, monkeypatch):
    parts = fixture(tmp_path)
    def unreadable(*args, **kwargs):
        raise PermissionError("selected replacement cannot be read")
    monkeypatch.setattr(af, "read_reference_variant_table", unreadable)
    with pytest.raises(PermissionError, match="selected replacement cannot be read") as caught:
        validate(parts)
    assert caught.value.study_decisions["eaf_reference_validation"] == parts[2]["eaf_reference_validation"]
    assert "eaf_replacement_validation" not in caught.value.study_decisions


def test_chromosome_without_diagnostic_candidates_does_not_shrink_fixed_n(tmp_path):
    def modify(source, default, replacement):
        for path in (default, replacement):
            pl.read_csv(path, separator="\t", infer_schema_length=0).with_columns(
                pl.when(pl.col("CHROM") == "2").then(pl.lit("T")).otherwise(pl.col("ALT")).alias("ALT")
            ).write_csv(path, separator="\t")
        return source.with_columns(
            pl.when(pl.col("CHR") == "2").then(pl.lit("T")).otherwise(pl.col("EA")).alias("EA"),
            pl.when(pl.col("CHR") == "2").then(pl.lit("A")).otherwise(pl.col("OA")).alias("OA"))
    parts = fixture(tmp_path, chromosomes=("1", "2"), modify=modify)
    evidence = validate(parts)["eaf_replacement_validation"]
    assert evidence["comparable_variants"] == 6 and evidence["study_variant_count"] == 12
    assert evidence["observed_overlap_fraction"] == .5


def test_chromosome_template_matches_one_global_whole_file_decision(tmp_path):
    parts = fixture(tmp_path, chromosomes=("1", "2"))
    whole = validate(parts)["eaf_replacement_validation"]
    reference = pl.read_csv(tmp_path / "replacement.tsv", separator="\t", infer_schema_length=0)
    for chromosome in parts[5].paths:
        path = tmp_path / ("replacement_" + chromosome + ".tsv")
        reference.filter(pl.col("CHROM") == chromosome).write_csv(path, separator="\t")
        parts[2]["eaf_replacement"]["resolved_files"][chromosome] = str(path)
        parts[3][chromosome]["user_eaf_file"] = str(path)
    parts[2]["eaf_replacement"]["file"] = str(tmp_path / "replacement_{chromosome}.tsv")
    templated = validate(parts)["eaf_replacement_validation"]
    for key in ("status", "comparable_variants", "major_count", "observed_overlap_fraction", "correlation", "mean_absolute_difference"):
        assert templated[key] == whole[key]


@pytest.mark.parametrize("problem", ["missing", "malformed", "column", "mapping", "wrong_build"])
def test_invalid_selected_sources_never_supply_approval(tmp_path, problem):
    parts = fixture(tmp_path)
    path = tmp_path / "replacement.tsv"
    if problem == "missing":
        path.unlink()
    elif problem == "malformed":
        path.write_text("not a valid table\n")
    elif problem == "column":
        parts[2]["eaf_replacement"]["column"] = "ABSENT_AF"
    elif problem == "mapping":
        with pytest.raises(af.AlleleFrequencyError, match="mapping.*missing|missing.*mapping"):
            af.validate_replacement_allele_frequency(parts[5], parts[1], parts[2], parts[3], {}, policies=parts[4])
        return
    else:
        pl.read_csv(path, separator="\t").with_columns((pl.col("POS") + 100000).alias("POS")).write_csv(path, separator="\t")
    with pytest.raises(af.AlleleFrequencyError) as caught:
        validate(parts)
    if problem == "wrong_build":
        assert caught.value.study_decisions["eaf_replacement_validation"]["comparable_variants"] == 0
    else:
        assert "eaf_replacement_validation" not in caught.value.study_decisions


@pytest.mark.parametrize("same", [False, True])
def test_preparation_and_adoption_read_and_match_distinct_source_once(tmp_path, monkeypatch, same):
    parts = fixture(tmp_path, same=same)
    calls = {"read": 0, "join": 0}
    for name, key in (("read_reference_variant_table", "read"), ("allele_oriented_left_join", "join")):
        original = getattr(af, name)
        def observed(*args, _original=original, _key=key, **kwargs):
            calls[_key] += 1
            return _original(*args, **kwargs)
        monkeypatch.setattr(af, name, observed)
    decisions = validate(parts)
    expected = 0 if same else 1
    assert calls == {"read": expected, "join": expected}
    result, qc, columns = adopt(parts, decisions)
    assert calls == {"read": expected, "join": expected}
    assert result.height == 12
    assert result[columns["eaf_col"]].null_count() == 0


def test_streaming_atomic_partition_writer_checks_rows_and_preserves_old_file(tmp_path):
    from postgwas.modules.harmonisation.chromosome_partition import _write_parquet_atomically
    path = tmp_path / "work.parquet"
    original = pl.DataFrame({"row": [1, 2], "value": ["a", "b"]})
    original.write_parquet(path)
    before = path.read_bytes()
    nested = pl.DataFrame({"row": [1, 2], "cache": [[{"value": .8}], [{"value": .2}]]})
    joined = pl.scan_parquet(path).join(nested.lazy(), on="row", how="left", coalesce=True)
    with pytest.raises(ValueError, match="row-count"):
        _write_parquet_atomically(joined, path, compression="zstd", expected_rows=3)
    assert path.read_bytes() == before
    _write_parquet_atomically(joined, path, compression="zstd", expected_rows=2)
    assert pl.read_parquet(path).select(original.columns).equals(original)
    assert pl.read_parquet(path)["cache"].to_list() == nested["cache"].to_list()
    assert list(tmp_path.iterdir()) == [path]


@pytest.mark.parametrize("change", ["missing", "null", "truncated", "value", "key", "disk_value"])
def test_cached_assignments_are_row_bound_and_cannot_fall_back(tmp_path, monkeypatch, change):
    parts = fixture(tmp_path)
    decisions = validate(parts)
    binding = decisions["eaf_replacement_validation"]["contract"]["chromosomes"]["1"]
    cache_column = binding["assignment_cache"]["column"]
    frame = pl.read_parquet(parts[5].paths["1"])
    if change == "missing":
        frame = frame.drop(cache_column)
    elif change == "null":
        frame = frame.with_columns(pl.lit(None, dtype=frame.schema[cache_column]).alias(cache_column))
    elif change == "truncated":
        frame = frame.with_columns(pl.col(cache_column).list.head(0))
    else:
        values = frame[cache_column].to_list()
        field = "EA" if change == "key" else binding["assignment_cache"]["assignment"]["value"]
        values[0][0][field] = "N" if change == "key" else .123
        frame = frame.with_columns(pl.Series(cache_column, values, dtype=frame.schema[cache_column]))
        if change == "disk_value":
            frame.write_parquet(parts[5].paths["1"])
    monkeypatch.setattr(af, "read_reference_variant_table", lambda *a, **k: pytest.fail("cache failure reread source"))
    with pytest.raises(af.AlleleFrequencyError, match="replacement.*validation|replacement.*cache"):
        adopt(parts, decisions, frame)


@pytest.mark.parametrize("expected", [None, True, -1, 2.0, "2"])
def test_streaming_atomic_writer_requires_integer_known_count(tmp_path, expected):
    from postgwas.modules.harmonisation.chromosome_partition import _write_parquet_atomically
    with pytest.raises(ValueError, match="expected row count"):
        _write_parquet_atomically(pl.DataFrame({"a": [1, 2]}).lazy(), tmp_path / "work.parquet",
            compression="zstd", expected_rows=expected)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("failure", ["sink", "second_partition", "source_after_write"])
def test_publication_failure_never_issues_approval_or_leaves_temporary_files(tmp_path, monkeypatch, failure):
    parts = fixture(tmp_path, chromosomes=("1", "2"))
    before = {path: Path(path).read_bytes() for path in parts[5].paths.values()}
    if failure == "sink":
        def failed_sink(*args, **kwargs):
            raise RuntimeError("streaming publication failed")
        monkeypatch.setattr(pl.LazyFrame, "sink_parquet", failed_sink)
    else:
        writer = af._write_parquet_atomically
        calls = []
        def failed_write(*args, **kwargs):
            calls.append(args)
            if failure == "second_partition" and len(calls) == 2:
                raise RuntimeError("second partition publication failed")
            writer(*args, **kwargs)
            if failure == "source_after_write":
                source = tmp_path / "replacement.tsv"
                source.write_text(source.read_text().replace("0.8", "0.81"))
        monkeypatch.setattr(af, "_write_parquet_atomically", failed_write)
    with pytest.raises((RuntimeError, af.StudyAlleleFrequencyError)) as caught:
        validate(parts)
    evidence = caught.value.study_decisions["eaf_replacement_validation"]
    assert evidence["status"] == "invalid" and "contract" not in evidence
    assert evidence["failure_reason"] == "replacement_cache_publication_failed"
    assert set(tmp_path.iterdir()) == {Path(p) for p in before} | {tmp_path / "default.tsv", tmp_path / "replacement.tsv"}
    if failure == "sink":
        assert {path: Path(path).read_bytes() for path in before} == before
    elif failure == "second_partition":
        assert Path(parts[5].paths["2"]).read_bytes() == before[parts[5].paths["2"]]


@pytest.mark.parametrize("subset", [False, True])
@pytest.mark.parametrize("case", ["palindrome", "unmatched", "unmatched_reverse", "multiallelic", "unresolved", "unmatched_site", "discarded_site", "collision", "missing_allele"])
def test_complete_cache_matches_fresh_assignment_values_and_scientific_qc(tmp_path, monkeypatch, subset, case):
    def modify(source, default, replacement):
        if case == "palindrome":
            source = source.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit("T")).otherwise(pl.col("EA")).alias("EA"))
            for path in (default, replacement):
                ref = pl.read_csv(path, separator="\t")
                ref.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit("T")).otherwise(pl.col("ALT")).alias("ALT")).write_csv(path, separator="\t")
        elif case in ("unmatched", "unmatched_reverse"):
            pl.read_csv(default, separator="\t").filter(pl.col("POS") != 100).write_csv(default, separator="\t")
            if case == "unmatched_reverse":
                source = source.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit("G")).otherwise(pl.col("EA")).alias("EA"),
                    pl.when(pl.col("POS") == 100).then(pl.lit("T")).otherwise(pl.col("OA")).alias("OA"))
        elif case in ("multiallelic", "unresolved", "unmatched_site", "discarded_site"):
            reference = pl.read_csv(replacement, separator="\t")
            if case == "unmatched_site":
                reference = reference.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit("G")).otherwise(pl.col("ALT")).alias("ALT"))
            else:
                reference = reference.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit("C")).otherwise(pl.col("REF")).alias("REF"),
                    pl.when(pl.col("POS") == 100).then(pl.lit("A")).otherwise(pl.col("ALT")).alias("ALT"),
                    pl.when(pl.col("POS") == 100).then(pl.lit(.2)).otherwise(pl.col("AF")).alias("AF"))
            reference.write_csv(replacement, separator="\t")
            with replacement.open("a") as handle:
                handle.write("1\t100\t%s\t%s\t%s\n" % (
                    "A" if case == "unmatched_site" else "C",
                    "T" if case == "unmatched_site" else "A" if case == "discarded_site" else "G",
                    "NA" if case == "unresolved" else ".1"))
        elif case == "missing_allele":
            source = source.with_columns(pl.when(pl.col("POS") == 100).then(None).otherwise(pl.col("EA")).alias("EA"))
        elif case == "collision":
            source = source.with_columns(pl.lit("original cache-like column").alias("__postgwas_replacement_assignments"),
                pl.lit("original raw-value-like column").alias("__postgwas_reference_value__"),
                pl.lit(.42).alias("AF"))
        return source
    parts = fixture(tmp_path, modify=modify, policy_overrides={"strand.unmatched_action": "retain"})
    if case == "unmatched_reverse":
        parts[2]["strand"] = "reverse"
    # Explicit retention exercises both possible unmatched output keys.
    decisions = validate(parts)
    frame = pl.read_parquet(parts[5].paths["1"])
    if subset:
        frame = frame.slice(1, 7).reverse()
    consume = af._consume_replacement_assignments
    monkeypatch.setattr(af, "_consume_replacement_assignments", fresh_assignments)
    expected, expected_qc, expected_columns = adopt(parts, decisions, frame)
    monkeypatch.setattr(af, "_consume_replacement_assignments", consume)
    monkeypatch.setattr(af, "read_reference_variant_table", lambda *a, **k: pytest.fail("worker reread replacement"))
    observed, observed_qc, observed_columns = adopt(parts, decisions, frame)
    assert observed.equals(expected)
    assert observed_columns == expected_columns
    for key in set(observed_qc) | set(expected_qc):
        if key in ("external_merge_external_rows_loaded", "external_merge_reused_prepared_replacement_assignments"):
            continue
        assert observed_qc.get(key) == expected_qc.get(key), key
    assert observed.select(parts[0].columns).drop("EA", "OA", "BETA", "Z").equals(
        frame.select(parts[0].columns).drop("EA", "OA", "BETA", "Z"))


@pytest.mark.parametrize("when", ["after_measurement", "during_write"])
def test_partition_alias_cannot_redirect_cache_publication(tmp_path, monkeypatch, when):
    parts = list(fixture(tmp_path))
    physical = Path(parts[5].paths["1"])
    alias = tmp_path / "alias.parquet"
    alias.symlink_to(physical)
    alternate = tmp_path / "unvalidated.parquet"
    alternate.write_bytes(physical.read_bytes())
    untouched = alternate.read_bytes()
    parts[5] = ChromosomePartitions({"1": str(alias)}, parts[5].rows)
    def retarget():
        alias.unlink()
        alias.symlink_to(alternate)
    if when == "after_measurement":
        measure = af.validate_replacement_frequency_pairs
        def measured(*args, **kwargs):
            result = measure(*args, **kwargs)
            retarget()
            return result
        monkeypatch.setattr(af, "validate_replacement_frequency_pairs", measured)
    else:
        writer = af._write_parquet_atomically
        def written(frame, destination, **kwargs):
            assert destination == physical
            retarget()
            writer(frame, destination, **kwargs)
        monkeypatch.setattr(af, "_write_parquet_atomically", written)
    with pytest.raises(af.StudyAlleleFrequencyError, match="partition changed") as caught:
        validate(parts)
    assert "contract" not in caught.value.study_decisions["eaf_replacement_validation"]
    assert alternate.read_bytes() == untouched


def test_assignment_descriptor_mutation_fails_before_orientation(tmp_path, monkeypatch):
    parts = fixture(tmp_path)
    decisions = validate(parts)
    descriptor = decisions["eaf_replacement_validation"]["contract"]["chromosomes"]["1"]["assignment_cache"]
    descriptor["assignment"]["value"] = "different_field"
    monkeypatch.setattr(af, "harmonise_strand_orientation", lambda *a, **k: pytest.fail("corrupted descriptor reached strand"))
    with pytest.raises(af.AlleleFrequencyError, match="cache descriptor changed"):
        adopt(parts, decisions)


@pytest.mark.parametrize("change", [
    "none", "value", "null_list", "inner_null", "boundary", "row_order", "entry_order",
    "nan_to_null", "nan_same", "original_af", "original_input", "dtype", "nested_dtype", "numeric_kind", "field_name",
    "column_order", "row_count", "empty", "all_null_lists", "all_empty_lists",
])
def test_prepared_frequency_row_comparison_matches_exact_native_semantics(change):
    dtype = pl.List(pl.Struct([pl.Field("key", pl.Int64), pl.Field("af", pl.Float64),
        pl.Field("allele", pl.String), pl.Field("eligible", pl.Boolean)]))
    record = {"key": 1, "af": .2, "allele": "A", "eligible": True}
    other = dict(record, key=2, af=float("nan"))
    packed = [[record, None], [], None, [other, record]]
    lhs = pl.DataFrame({"source": [0, 1, 2, 3], "original_af": [.2, None, .3, float("nan")],
        "original_input": [{"ea": "A", "oa": "C"}] * 4,
        "reference": pl.Series("reference", packed, dtype=dtype),
        "replacement": pl.Series("replacement", packed, dtype=dtype)})
    rhs = lhs.clone()
    changed = deepcopy(packed)
    if change == "value":
        changed[0][0]["af"] = .25
    elif change == "null_list":
        changed[1] = None
    elif change == "inner_null":
        # Polars 0.20 treats these parent-null representations as equal.
        changed[0][1] = {name: None for name in record}
    elif change == "boundary":
        changed[0], changed[1] = [record], [None]
    elif change == "row_order":
        changed[0], changed[1] = changed[1], changed[0]
    elif change == "entry_order":
        changed[3] = list(reversed(changed[3]))
    elif change == "nan_to_null":
        changed[3][0]["af"] = None
    elif change == "nan_same":
        changed[3][0]["af"] = float("nan")
    rhs = rhs.with_columns(pl.Series("replacement", changed, dtype=dtype))
    if change == "original_af":
        rhs = rhs.with_columns(pl.lit(.25).alias("original_af"))
    elif change == "original_input":
        rhs = rhs.with_columns(pl.struct(pl.lit("T").alias("ea"), pl.lit("C").alias("oa")).alias("original_input"))
    elif change == "dtype":
        rhs = rhs.with_columns(pl.col("source").cast(pl.Int32))
    elif change == "nested_dtype":
        rhs = rhs.with_columns(pl.col("replacement").cast(pl.List(pl.Struct([
            pl.Field("key", pl.Int32), pl.Field("af", pl.Float64),
            pl.Field("allele", pl.String), pl.Field("eligible", pl.Boolean)]))))
    elif change == "numeric_kind":
        rhs = rhs.with_columns(pl.col("source").cast(pl.Float64))
    elif change == "field_name":
        rhs = rhs.rename({"original_af": "different_af"})
    elif change == "column_order":
        rhs = rhs.select(list(reversed(rhs.columns)))
    elif change == "row_count":
        rhs = rhs.head(3)
    elif change == "empty":
        lhs, rhs = lhs.head(0), rhs.head(0)
    elif change in ("all_null_lists", "all_empty_lists"):
        values = [None if change == "all_null_lists" else []] * 4
        lhs = lhs.with_columns(pl.Series("replacement", values, dtype=dtype))
        rhs = rhs.with_columns(pl.Series("replacement", values, dtype=dtype))
    expected = lhs.equals(rhs, null_equal=True)
    assert af._prepared_frequency_rows_equal(lhs, rhs, ["reference", "replacement"]) is expected


def test_matching_source_order_guard_skips_join_and_sort_without_omitting_values(tmp_path, monkeypatch):
    parts = fixture(tmp_path)
    decisions = validate(parts)
    frame = pl.read_parquet(parts[5].paths["1"])
    monkeypatch.setattr(pl.DataFrame, "join", lambda *a, **k: pytest.fail("ordered guard joined"))
    monkeypatch.setattr(pl.DataFrame, "sort", lambda *a, **k: pytest.fail("ordered guard sorted"))
    assert af._require_replacement_validation("1", frame, parts[1], decisions,
        parts[3]["1"]["user_eaf_file"], MAPPING, parts[4])
    with pytest.raises(af.AlleleFrequencyError, match="rows changed"):
        af._require_replacement_validation("1", frame.with_columns(pl.lit(.123).alias("MAF")), parts[1], decisions,
            parts[3]["1"]["user_eaf_file"], MAPPING, parts[4])


@pytest.mark.parametrize("frequency_column", [
    "__postgwas_reference_value__", "__postgwas_allele_match_orientation__",
    "__postgwas_reference_frequency", "__postgwas_reference_reason",
    "__postgwas_reference_multiallelic", "__postgwas_multiallelic_resolved_site",
    "__postgwas_unresolved_site", "__postgwas_external_eaf_assignment_reason",
    "__postgwas_external_eaf_biallelic", "__postgwas_external_eaf_row", "AF",
])
@pytest.mark.parametrize("valid", [False, True])
def test_mapped_original_frequency_private_names_preserve_gate_and_adoption(tmp_path, monkeypatch, frequency_column, valid):
    values = None if valid else [.25] * 12
    ordinary = fixture(tmp_path / "ordinary", replacement_values=values)
    parts = fixture(tmp_path / "collision", replacement_values=values, frequency_column=frequency_column)
    original = pl.read_parquet(parts[5].paths["1"])
    if valid:
        expected_decisions = validate(ordinary)
        expected, _, expected_columns = adopt(ordinary, expected_decisions)
    else:
        with pytest.raises(af.StudyAlleleFrequencyError, match="inconclusive") as caught:
            validate(ordinary)
        expected_decisions = caught.value.study_decisions
    calls = {"read": 0, "join": 0}
    for name, key in (("read_reference_variant_table", "read"), ("allele_oriented_left_join", "join")):
        implementation = getattr(af, name)
        def observed(*args, _implementation=implementation, _key=key, **kwargs):
            calls[_key] += 1
            return _implementation(*args, **kwargs)
        monkeypatch.setattr(af, name, observed)
    if valid:
        decisions = validate(parts)
    else:
        with pytest.raises(af.StudyAlleleFrequencyError, match="inconclusive") as caught:
            validate(parts)
        decisions = caught.value.study_decisions
    evidence = decisions["eaf_replacement_validation"]
    baseline = expected_decisions["eaf_replacement_validation"]
    for field in ("status", "study_variant_count", "comparable_variants", "major_count",
                  "correlation", "mean_absolute_difference", "failed_requirements"):
        assert evidence.get(field) == baseline.get(field), field
    assert decisions["eaf_reference_validation"] == parts[2]["eaf_reference_validation"]
    published = pl.read_parquet(parts[5].paths["1"])
    assert published.select(original.columns).equals(original)
    if not valid:
        assert "contract" not in evidence
        assert published.equals(original)
        assert calls == {"read": 1, "join": 1}
        return
    result, qc, columns = adopt(parts, decisions)
    assert calls == {"read": 1, "join": 1}
    assert result[frequency_column].equals(parts[0][frequency_column])
    # Only the explicitly mapped original/destination labels differ from MAF.
    expected_names = [frequency_column if name == "MAF" else
        name + "_" if frequency_column == "AF" and name in ("AF", "AF_") else name
        for name in expected.columns]
    assert result.columns == expected_names
    normalized = result.select(pl.col(actual).alias(baseline)
        for actual, baseline in zip(result.columns, expected.columns))
    assert normalized.schema == expected.schema
    assert normalized.equals(expected)
    monkeypatch.setattr(af, "_consume_replacement_assignments", fresh_assignments)
    fresh_result, fresh_qc, fresh_columns = adopt(parts, decisions)
    assert result.schema == fresh_result.schema and result.equals(fresh_result)
    assert columns == fresh_columns
    for key in set(qc) | set(fresh_qc):
        if key not in ("external_merge_external_rows_loaded", "external_merge_reused_prepared_replacement_assignments"):
            assert qc.get(key) == fresh_qc.get(key), key


@pytest.mark.parametrize("when", ["before_writer", "before_stream", "after_stream"])
def test_partition_scientific_mutation_during_publication_never_approves(tmp_path, monkeypatch, when):
    parts = fixture(tmp_path)
    path = Path(parts[5].paths["1"])
    original = pl.read_parquet(path)
    changed = original.with_columns(pl.lit(.9).alias("BETA"))
    if when == "before_writer":
        writer = af._write_parquet_atomically
        def mutate_before(frame, destination, **kwargs):
            changed.write_parquet(path)
            return writer(frame, destination, **kwargs)
        monkeypatch.setattr(af, "_write_parquet_atomically", mutate_before)
    else:
        sink = pl.LazyFrame.sink_parquet
        def mutate_after(frame, *args, **kwargs):
            if when == "before_stream":
                changed.write_parquet(path)
            result = sink(frame, *args, **kwargs)
            if when == "after_stream":
                changed.write_parquet(path)
            return result
        monkeypatch.setattr(pl.LazyFrame, "sink_parquet", mutate_after)
    with pytest.raises(af.StudyAlleleFrequencyError, match="partition changed") as caught:
        validate(parts)
    evidence = caught.value.study_decisions["eaf_replacement_validation"]
    assert evidence["status"] == "invalid" and "contract" not in evidence
    assert evidence["failure_reason"] == "replacement_cache_publication_failed"
    # Refuse the cache replacement, leaving the independently mutated file alone.
    assert pl.read_parquet(path).equals(changed)
    assert set(tmp_path.iterdir()) == {path, tmp_path / "default.tsv", tmp_path / "replacement.tsv"}


@pytest.mark.parametrize("frequency_column,replacement_column,occupied", [
    ("MAF", "__postgwas_cached_order", 0),
    ("MAF", "__postgwas_cached_assignment", 0),
    ("MAF", "__postgwas_reference_value__", 0),
    ("MAF", "__postgwas_reference_frequency", 0),
    ("MAF", "__postgwas_reference_reason", 0),
    ("MAF", "__postgwas_reference_multiallelic", 0),
    ("MAF", "__postgwas_multiallelic_resolved_site", 0),
    ("MAF", "__postgwas_unresolved_site", 0),
    ("MAF", "__postgwas_external_eaf_assignment_reason", 0),
    ("MAF", "__postgwas_cached_order", 3),
    ("MAF", "__postgwas_reference_value__", 3),
    ("__postgwas_reference_value__", "__postgwas_reference_value___", 0),
])
def test_replacement_source_names_preserve_exact_fresh_schema_mapping_qc(tmp_path, monkeypatch, frequency_column, replacement_column, occupied):
    def modify(source, default, replacement):
        return source.with_columns(pl.lit(.42).alias(replacement_column + "_" * i) for i in range(occupied))
    parts = fixture(tmp_path, modify=modify, frequency_column=frequency_column)
    path = Path(parts[3]["1"]["user_eaf_file"])
    pl.read_csv(path, separator="\t").rename({"AF": replacement_column}).write_csv(path, separator="\t")
    parts[2]["eaf_replacement"]["column"] = replacement_column
    calls = {"read": 0, "join": 0}
    for name, key in (("read_reference_variant_table", "read"), ("allele_oriented_left_join", "join")):
        implementation = getattr(af, name)
        def observed(*args, _implementation=implementation, _key=key, **kwargs):
            calls[_key] += 1
            return _implementation(*args, **kwargs)
        monkeypatch.setattr(af, name, observed)
    decisions = validate(parts)
    result, qc, columns = adopt(parts, decisions)
    assert calls == {"read": 1, "join": 1}
    assert result[frequency_column].equals(parts[0][frequency_column])
    monkeypatch.setattr(af, "_consume_replacement_assignments", fresh_assignments)
    expected, expected_qc, expected_columns = adopt(parts, decisions)
    assert result.schema == expected.schema and result.equals(expected)
    assert columns == expected_columns
    ignored = {"external_merge_external_rows_loaded", "external_merge_reused_prepared_replacement_assignments"}
    assert {k: v for k, v in qc.items() if k not in ignored} == {
        k: v for k, v in expected_qc.items() if k not in ignored}


@pytest.mark.parametrize("orientations", [["direct"], ["direct", "swap"]])
@pytest.mark.parametrize("subset", [False, True])
@pytest.mark.parametrize("unresolved", [False, True])
def test_cached_frequency_metadata_matches_fresh_matching_mode(tmp_path, monkeypatch, orientations, subset, unresolved):
    def modify(source, default, replacement):
        with replacement.open("a") as handle:
            handle.write("1\t100\tA\tG\t%s\n" % ("NA" if unresolved else ".1"))
        return source
    parts = fixture(tmp_path, modify=modify, policy_overrides={"eaf.match_orientations": orientations})
    decisions = validate(parts)
    frame = pl.read_parquet(parts[5].paths["1"])
    if subset:
        frame = frame.slice(1, 7).reverse()
    result, qc, columns = adopt(parts, decisions, frame)
    monkeypatch.setattr(af, "_consume_replacement_assignments", fresh_assignments)
    expected, expected_qc, expected_columns = adopt(parts, decisions, frame)
    assert result.schema == expected.schema and result.equals(expected)
    assert columns == expected_columns
    ignored = {"external_merge_external_rows_loaded", "external_merge_reused_prepared_replacement_assignments"}
    assert {k: v for k, v in qc.items() if k not in ignored} == {
        k: v for k, v in expected_qc.items() if k not in ignored}
    for suffix in ("multiallelic_sites", "frequency_unresolved_sites", "frequency_reconstructed_rows", "frequency_unresolved_rows"):
        assert ("external_merge_reference_" + suffix in qc) is ("swap" in orientations)
