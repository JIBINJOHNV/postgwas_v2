"""Reference frequency attachment reuses evidence and preserves source columns."""

from unittest.mock import patch

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation import allele_frequency as af, strand
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.rejects import SOURCE_INPUT_ROW_COLUMN
from postgwas.modules.harmonisation.chromosome_partition import ChromosomePartitions


REFERENCE_MAPPING = {
    "chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF",
    "delimiter": "tab",
}
GUARD_KEY = "eaf_swap_guard_col"
ELIGIBILITY_KEY = "strand_reference_biallelic_col"


def _study(tmp_path, *, same_source):
    frame = pl.DataFrame({
        "CHR": ["1"] * 4, "POS": [10, 20, 30, 40],
        "EA": ["G" if same_source else "A", "T", "C", "C"], "OA": ["A" if same_source else "G", "C", "A", "A"],
        "MAF": [.3, .3, .2, .4], "strand_reference_biallelic": [.2, .3, .4, .5],
        "__postgwas_external_eaf_swap_guard": ["source one", "source two", "source three", "source four"],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN)
    columns = {
        "chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA",
        "eaf_col": "MAF", "beta_or_col": "strand_reference_biallelic",
        "output_folder": str(tmp_path / "output"), "gwas_outputname": "study",
    }
    return frame, columns


def _reference(path, *, extra=.1, reverse=False):
    first = "G\tA\t0.7" if reverse else "A\tG\t0.3"
    text = "CHROM\tPOS\tREF\tALT\tEUR\n1\t10\t" + first + "\n1\t20\tC\tT\t0.7\n1\t30\tA\tC\t0.8\n1\t40\tA\tC\t0.6\n"
    if extra is not None:
        text += "1\t10\tA\tC\t%s\n" % extra
    path.write_text(text)
    return path


def _decision(reference):
    return {
        "strand": "forward", "effect_type": "beta", "eaf_is_maf": True,
        "eaf_input_interpretation": "maf", "eaf_selected_source": "replacement_reference",
        "eaf_reference_validation": {"status": "replacement_selected"},
        "eaf_replacement": {
            "file": str(reference), "column": "EUR", "resolved_files": {"1": str(reference)},
        },
    }


def _approved(source, columns, comparison, replacement):
    policies = default_policies().with_overrides({"eaf.maf_reference_min_overlap": 2})
    maps = {"1": {"default_eaf_file": str(comparison), "default_eaf_reference_column": "EUR",
        "user_eaf_file": str(replacement)}}
    decisions, work = af.validate_study_allele_frequency(source, columns,
        dict(_decision(replacement), eaf_column="MAF"), maps, REFERENCE_MAPPING,
        policies=policies, eaf_replacement_file=str(replacement), eaf_replacement_column="EUR")
    decisions["eaf_replacement"]["resolved_files"] = {"1": str(replacement)}
    path = comparison.parent / "work.parquet"
    work.write_parquet(path)
    decisions = af.validate_replacement_allele_frequency(
        ChromosomePartitions({"1": str(path)}, {"1": work.height}), columns, decisions,
        maps, REFERENCE_MAPPING, policies=policies)
    return pl.read_parquet(path), decisions, policies


@pytest.mark.parametrize("prepared", [False, True])
def test_same_reference_replacement_reuses_aligned_frequencies_without_reread(
    tmp_path, prepared,
):
    reference = _reference(tmp_path / "reference.tsv")
    source, columns = _study(tmp_path, same_source=True)
    work, decisions, policies = _approved(source, columns, reference, reference)
    if not prepared:
        with pytest.raises(af.AlleleFrequencyError, match="replacement validation input columns"):
            af.harmonise_allele_frequency("1", source, columns,
                load_configuration().modules.harmonisation.output_layout.model_dump(), "\t",
                external_eaf_colmap=REFERENCE_MAPPING, default_eaf_file=str(reference),
                default_eaf_column="EUR", default_eaf_colmap=REFERENCE_MAPPING,
                study_decision=decisions, policies=policies)
        return
    layout = load_configuration().modules.harmonisation.output_layout.model_dump()
    with patch.object(af, "read_reference_variant_table", side_effect=AssertionError("replacement reread")):
        with patch.object(strand, "_load_reference", wraps=strand._load_reference) as reads:
            result, qc, resolved = af.harmonise_allele_frequency(
                "1", work, columns, layout, "\t", external_eaf_colmap=REFERENCE_MAPPING,
                default_eaf_file=str(reference), default_eaf_column="EUR",
                default_eaf_colmap=REFERENCE_MAPPING, study_decision=decisions,
                policies=policies,
            )
    assert reads.call_count == 0
    assert result.select(source.columns).equals(source)
    assert result[resolved["eaf_col"]].to_list() == pytest.approx([.3, .7, .8, .6])
    assert ELIGIBILITY_KEY not in resolved
    assert GUARD_KEY not in resolved
    assert result["__postgwas_external_eaf_swap_guard"].equals(source["__postgwas_external_eaf_swap_guard"])
    assert qc["external_merge_reused_aligned_strand_reference"] is True


def test_fresh_replacement_preserves_source_columns_without_export_guard_metadata(tmp_path):
    comparison = _reference(tmp_path / "comparison.tsv", extra=None, reverse=True)
    replacement = _reference(tmp_path / "replacement.tsv")
    source, columns = _study(tmp_path, same_source=False)
    # The resolved multiallelic row now contributes to diagnostics. Its original
    # GWAS EA is REF with B=1-(.3+.1)=.6, so align the unchanged input MAF and
    # default comparison to that same allele rather than relying on exclusion.
    source = source.with_columns(pl.when(pl.col("POS") == 10).then(.4)
        .otherwise(pl.col("MAF")).alias("MAF"))
    pl.read_csv(comparison, separator="\t").with_columns(
        pl.when(pl.col("POS") == 10).then(.6).otherwise(pl.col("EUR")).alias("EUR")
    ).write_csv(comparison, separator="\t")
    work, decisions, policies = _approved(source, columns, comparison, replacement)
    assert decisions["eaf_replacement_validation"]["comparable_variants"] == 4
    assert decisions["eaf_replacement_validation"]["correlation"] == pytest.approx(1.)
    layout = load_configuration().modules.harmonisation.output_layout.model_dump()
    with patch.object(af, "read_reference_variant_table", wraps=af.read_reference_variant_table) as reads:
        result, qc, resolved = af.harmonise_allele_frequency(
            "1", work, columns, layout, "\t", eaffile=str(replacement),
            external_eaf_colmap=REFERENCE_MAPPING, default_eaf_file=str(comparison),
            default_eaf_column="EUR", default_eaf_colmap=REFERENCE_MAPPING,
            study_decision=decisions, policies=policies,
        )
    assert reads.call_count == 0
    assert result.select(source.columns).equals(source)
    assert result[resolved["eaf_col"]].to_list() == pytest.approx([.6, .7, .8, .6])
    assert GUARD_KEY not in resolved
    assert ELIGIBILITY_KEY not in resolved
    assert "external_merge_frequency_swap_guard_column" not in qc


@pytest.mark.parametrize("same_source", [False, True])
def test_safe_biallelic_references_do_not_allocate_export_guard(tmp_path, same_source):
    replacement = _reference(tmp_path / "replacement.tsv", extra=None)
    comparison = replacement if same_source else _reference(
        tmp_path / "comparison.tsv", extra=None, reverse=True,
    )
    source, columns = _study(tmp_path, same_source=same_source)
    work, decisions, policies = _approved(source, columns, comparison, replacement)
    layout = load_configuration().modules.harmonisation.output_layout.model_dump()
    result, _, resolved = af.harmonise_allele_frequency(
        "1", work, columns, layout, "\t", eaffile=str(replacement),
        external_eaf_colmap=REFERENCE_MAPPING, default_eaf_file=str(comparison),
        default_eaf_column="EUR", default_eaf_colmap=REFERENCE_MAPPING,
        study_decision=decisions, policies=policies,
    )
    assert GUARD_KEY not in resolved
    assert result["__postgwas_external_eaf_swap_guard"].equals(source["__postgwas_external_eaf_swap_guard"])


def test_internal_study_frequency_does_not_add_external_export_guard_context(tmp_path):
    reference = _reference(tmp_path / "reference.tsv", extra=None)
    source, columns = _study(tmp_path, same_source=True)
    source = source.with_columns(pl.Series("MAF", [.3, .7, .8, .6]))
    layout = load_configuration().modules.harmonisation.output_layout.model_dump()
    result, _, resolved = af.harmonise_allele_frequency(
        "1", source, columns, layout, "\t", external_eaf_colmap=REFERENCE_MAPPING,
        default_eaf_file=str(reference), default_eaf_column="EUR",
        default_eaf_colmap=REFERENCE_MAPPING,
        study_decision={"strand": "forward", "effect_type": "beta", "eaf_is_maf": False,
                        "eaf_selected_source": "study_supplied",
                        "eaf_reference_validation": {"status": "not_required"}},
        policies=default_policies(),
    )
    assert result.select(source.columns).equals(source)
    assert resolved["eaf_col"] == "MAF"
    assert GUARD_KEY not in resolved
    assert ELIGIBILITY_KEY not in resolved


@pytest.mark.parametrize("reference_role", ["value", "chromosome"])
def test_explicit_reference_header_names_remain_usable_without_export_guard(tmp_path, reference_role):
    collision = "__postgwas_external_eaf_swap_guard"
    chromosome_column = collision if reference_role == "chromosome" else "CHROM"
    value_column = collision if reference_role == "value" else "EUR"
    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "%s\tPOS\tREF\tALT\t%s\n1\t10\tA\tG\t.3\n1\t10\tA\tC\t.1\n"
        % (chromosome_column, value_column)
    )
    reference_bytes = reference.read_bytes()
    source = pl.DataFrame({"CHR": ["1"], "POS": [10], "EA": ["G"], "OA": ["A"]})

    result, destination, stats = af.merge_external_allele_frequencies(
        source, str(reference), value_column,
        dict(REFERENCE_MAPPING, chr=chromosome_column),
        {"chr": "CHR", "pos": "POS", "ea": "EA", "oa": "OA"},
        policies=default_policies(), declared_alt_frequency=True,
    )

    assert result[destination].to_list() == [.3]
    assert "frequency_swap_guard_column" not in stats
    assert result.select(source.columns).equals(source)
    assert reference.read_bytes() == reference_bytes
