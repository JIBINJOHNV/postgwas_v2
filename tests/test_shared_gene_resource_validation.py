"""Static resource reuse must not reuse a scientific consumer's policies."""

from dataclasses import FrozenInstanceError
import gzip
from pathlib import Path
from unittest.mock import patch

import pytest

from postgwas.config import load_configuration
from postgwas.core import gene_coordinates, gene_sets, snp_sets
from postgwas.core.gene_sets import GeneSetFormat, validate_gene_set_source
from postgwas.core.gene_coordinates import read_gene_coordinates
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.plink import count_plink_samples
from postgwas.core.resource_preparation import ResourcePreparationError
from postgwas.modules.gcta_gene.errors import GctaGeneError
from postgwas.modules.gcta_gene.pathway_sets import scan_gmt_pathways, validate_gene_coordinate_reference
from postgwas.modules.gcta_gene.service import _prepare_analysis_set_list, _validate_gene_list
from postgwas.modules.magma.annotations import read_gene_locations
from postgwas.modules.magma.analysis import parse_gene_set_file
from postgwas.modules.magma.errors import MagmaError


def test_magma_gene_rows_reused_and_results_are_detached(tmp_path):
    path = tmp_path / "genes.loc"
    path.write_text("ENSG1 1 10 20 + GENE1\nENSG2 2 30 40 -\n")
    module = load_configuration().modules.magma
    with patch.object(gene_coordinates, "open_text", wraps=gene_coordinates.open_text) as reader:
        with InputValidationSession() as session:
            with session.scope("magma"):
                first = read_gene_locations(path, module)
            first.clear()
            with session.scope("n_magma"):
                second = read_gene_locations(path, module, "nMAGMA gene-location reference")
    assert reader.call_count == 1
    assert second == {"ENSG1": ("1", 10, 20, "+", "GENE1"), "ENSG2": ("2", 30, 40, "-", None)}
    record = next(record for record in session.records if record.role == "Gene coordinates")
    assert record.consumers == ("magma", "n_magma")


def test_gcta_gene_rows_reused_for_different_windows_and_chromosome_policies(tmp_path):
    path = tmp_path / "genes.txt"
    path.write_text("1 10 20 GENE1\n")
    module = load_configuration().modules.gcta_gene
    original = Path.open
    reads = []

    def counted(source, *args, **kwargs):
        if source == path:
            reads.append(source)
        return original(source, *args, **kwargs)

    with patch.object(Path, "open", counted), InputValidationSession():
        assert _validate_gene_list(path, module)["genes"] == 1
        zero, _ = validate_gene_coordinate_reference(path, {"1"}, "exact", 0, module.gene_annotation.columns)
        expanded, _ = validate_gene_coordinate_reference(path, {"1"}, "exact", 50, module.gene_annotation.columns)
        with pytest.raises(ResourcePreparationError, match="disallowed chromosome"):
            validate_gene_coordinate_reference(path, {"2"}, "exact", 0, module.gene_annotation.columns)
        with pytest.raises(GctaGeneError, match="do not share"):
            _validate_gene_list(path, module, {"2"})
    assert len(reads) == 1
    assert (zero["GENE1"].start, zero["GENE1"].end) == (10, 20)
    assert (expanded["GENE1"].start, expanded["GENE1"].end) == (1, 70)


@pytest.mark.parametrize("row", [
    "1 0 10 GENE1", "1 20 10 GENE1", "1 NA 20 GENE1", "1 1_0 20 GENE1",
    "1 ١٠ 20 GENE1", "1 10 20", "1 10 20 GENE1 extra",
    "1 10 20 GENE1\n1 30 40 GENE1",
])
def test_common_gene_reader_rejects_malformed_rows(tmp_path, row):
    path = tmp_path / "genes.txt"
    path.write_text(row + "\n")
    with pytest.raises(ValueError):
        read_gene_coordinates(path, column_roles=("chromosome", "start", "end", "gene"))


def test_gene_reader_preserves_whitespace_around_tab_delimited_coordinates(tmp_path):
    path = tmp_path / "genes.txt"
    path.write_text("1\t 10 \t +20 \tGENE1\n")
    rows = read_gene_coordinates(
        path, column_roles=("chromosome", "start", "end", "gene"), delimiter_pattern=r"\t",
    )
    assert (rows[0].start, rows[0].end) == (10, 20)


def test_gene_rows_immutable_and_changed_resource_rejected(tmp_path):
    path = tmp_path / "genes.txt"
    path.write_text("1 1 1 GENE1\n")
    with InputValidationSession():
        rows = read_gene_coordinates(path, column_roles=("chromosome", "start", "end", "gene"))
        with pytest.raises(FrozenInstanceError):
            rows[0].start = 2
        path.write_text("1 1 2 GENE1\n")
        with pytest.raises(ValueError, match="changed after pipeline preflight"):
            read_gene_coordinates(path, column_roles=("chromosome", "start", "end", "gene"))


def test_magma_consumer_chromosome_policy_not_hidden_by_cache(tmp_path):
    path = tmp_path / "genes.loc"
    path.write_text("ENSG1 1 10 20 +\n")
    module = load_configuration().modules.magma
    with InputValidationSession():
        read_gene_locations(path, module)
        module.input.invalid_chromosome_labels.append("1")
        with pytest.raises(MagmaError, match="invalid gene ID, chromosome"):
            read_gene_locations(path, module)


def test_fam_inventory_reused_across_blank_row_policies(tmp_path):
    path = tmp_path / "reference.fam"
    path.write_text("F I 0 0 1 -9\n\n")
    original = Path.open
    reads = []

    def counted(source, *args, **kwargs):
        if source == path:
            reads.append(source)
        return original(source, *args, **kwargs)

    with patch.object(Path, "open", counted), InputValidationSession():
        assert count_plink_samples(path, allow_blank_rows=True) == 1
        with pytest.raises(ValueError, match="expected 6"):
            count_plink_samples(path, allow_blank_rows=False)
    assert len(reads) == 1


def test_fastbat_static_counts_reused_for_actual_membership_filtering(tmp_path):
    path = tmp_path / "sets.txt"
    path.write_text("SET1\nrs1\nrs2\nEND\nSET2\nrs1\nrs3\nEND\n")
    destination = tmp_path / "selected.txt"
    with patch.object(snp_sets.sqlite3, "connect", wraps=snp_sets.sqlite3.connect) as database:
        with InputValidationSession():
            assert snp_sets.validate_fastbat_set_list(path) == {
                "input_sets": 2, "requested_set_variants": 4, "unique_requested_set_variants": 3,
            }
            _, metrics = _prepare_analysis_set_list(path, destination, {"rs1", "rs2"}, {"rs1"}, "error", "error", 10)
    assert database.call_count == 1
    assert database.call_args.args == ("",)
    assert destination.read_text().split() == ["SET1", "rs1", "END", "SET2", "rs1", "END"]
    assert metrics["requested_set_variants"] == 4
    assert metrics["unique_requested_set_variants"] == 3
    assert metrics["analysis_set_variants"] == 2


@pytest.mark.parametrize("text,message", [
    ("END\n", "without a set ID"),
    ("SET\nEND\n", "contains no variant"),
    ("SET\nrs1\n", "missing its terminating"),
    ("SET\nrs1\nrs1\nEND\n", "repeats variant ID"),
    ("SET\nrs1\nEND\nSET\nrs2\nEND\n", "repeats set ID"),
    ("bad name\nrs1\nEND\n", "must not contain whitespace"),
    ("SET\nrs1 rs2\nEND\n", "one variant ID"),
    ("\n", "contains no sets"),
])
def test_fastbat_invalid_static_structure_rejected(tmp_path, text, message):
    path = tmp_path / "sets.txt"
    path.write_text(text)
    with InputValidationSession():
        with pytest.raises(GctaGeneError, match=message):
            snp_sets.validate_fastbat_set_list(path, error_type=GctaGeneError)


def test_fastbat_mutation_during_consumption_rejected(tmp_path):
    path = tmp_path / "sets.txt"
    path.write_text("SET\nrs1\nEND\n")
    with pytest.raises(GctaGeneError, match="changed after pipeline preflight"):
        with snp_sets.open_fastbat_set_memberships(path, error_type=GctaGeneError) as (_, events):
            assert list(events) == [("SET", "rs1"), ("SET", None)]
            path.write_text("SET\nrs2\nEND\n")


def test_gmt_inventory_is_compact_and_static_checks_run_once(tmp_path):
    path = tmp_path / "sets.gmt"
    path.write_text("# resource comment\n\nSET1\tdesc\tG1\tG1\tG2\nSET2\t\tG2\n")
    original = gene_sets._iter_gene_set_rows
    with patch.object(gene_sets, "_iter_gene_set_rows", wraps=original) as reader:
        with InputValidationSession():
            source = scan_gmt_pathways(path, "deduplicate")
            reused = scan_gmt_pathways(path, "deduplicate")
            first = list(source.rows())
            second = list(reused.rows())
    assert sum(call.kwargs["validate"] for call in reader.call_args_list) == 1
    assert source.gene_sets == 2
    assert source.genes == frozenset({"G1", "G2"})
    assert first == second
    assert first[0].genes == ("G1", "G2")
    assert first[0].duplicate_gene_entries == 1
    assert not hasattr(source, "memberships")


def test_gmt_explicit_evidence_reused_in_direct_preparation_without_rescan(tmp_path):
    path = tmp_path / "sets.gmt"
    path.write_text("SET\tdesc\tG1\n")
    source = scan_gmt_pathways(path, "deduplicate")
    with patch.object(gene_sets, "_iter_gene_set_rows", side_effect=AssertionError("Repeated static scan")):
        assert scan_gmt_pathways(path, "deduplicate", evidence=source) is source
    with pytest.raises(ResourcePreparationError, match="does not match"):
        scan_gmt_pathways(path, "error", evidence=source)


@pytest.mark.parametrize("row,message", [
    ("SET\tdesc\tG1\tG1\n", "duplicate gene entries"),
    ("SET\tdesc\tG1\t\n", "empty gene identifier"),
    ("END\tdesc\tG1\n", "reserved"),
    ("SET NAME\tdesc\tG1\n", "whitespace-containing"),
    ("SET\tdesc\tG1\nSET\tdesc\tG2\n", "Duplicate gene-set"),
])
def test_gcta_gmt_policies_not_relaxed_by_magma_validation(tmp_path, row, message):
    path = tmp_path / "sets.gmt"
    path.write_text(row)
    permissive = GeneSetFormat("gmt", r"\s+", "#", True, "drop", "deduplicate", False)
    with InputValidationSession():
        if "Duplicate gene-set" not in message:
            validate_gene_set_source(path, permissive)
        with pytest.raises(ResourcePreparationError, match=message):
            scan_gmt_pathways(path, "error")


def test_gmt_changed_after_inventory_is_rejected_before_membership_consumption(tmp_path):
    path = tmp_path / "sets.gmt"
    path.write_text("SET\tdesc\tG1\n")
    source = scan_gmt_pathways(path, "error")
    path.write_text("SET\tdesc\tG2\n")
    with pytest.raises(ValueError, match="changed after pipeline preflight"):
        list(source.rows())


def test_magma_gmt_native_auto_and_membership_header_semantics_preserved(tmp_path):
    class Logger:
        def record(self, *_args, **_kwargs):
            pass

    module = load_configuration().modules.magma
    path = tmp_path / "sets.gmt"
    path.write_text("# comment\nSET NAME\t\tG1\t\tG1\tG2\n")
    with InputValidationSession():
        first, metadata = parse_gene_set_file(path, module, Logger(), return_metadata=True)
        second = parse_gene_set_file(path, module, Logger())
    assert first.equals(second)
    assert first["input_genes"].to_list() == ["G1,G2"]
    assert metadata["detected_format"] == "gmt"
    native = tmp_path / "sets.txt"
    native.write_text("SET\tG1 G2\n")
    _, metadata = parse_gene_set_file(native, module, Logger(), return_metadata=True)
    assert metadata["detected_format"] == "magma"
    mixed = tmp_path / "mixed.txt"
    mixed.write_text("SET1\tdesc\tG1\nSET2 G2\n")
    with pytest.raises(MagmaError, match="mixes GMT and native"):
        parse_gene_set_file(mixed, module, Logger())
    membership = tmp_path / "memberships.txt"
    membership.write_text("SET GENE\nSET1 G1\nSET1 G2\n")
    module.gene_sets.membership_has_header = True
    frame = parse_gene_set_file(membership, module, Logger(), input_format="membership")
    assert frame["input_genes"].to_list() == ["G1,G2"]


def test_membership_custom_columns_delimiter_and_deduplication_reused(tmp_path):
    path = tmp_path / "memberships.txt"
    path.write_text("GENE|SET\nG1|SET1\nG1|SET1\n# comment\n\nG2|SET2\nG3|SET1\n")
    options = GeneSetFormat(
        "membership", r"\|", "#", True, "reject", "deduplicate", False,
        membership_has_header=True, membership_set_column=1, membership_gene_column=0,
    )
    original = gene_sets._iter_gene_membership_rows
    with patch.object(gene_sets, "_iter_gene_membership_rows", wraps=original) as reader:
        with InputValidationSession():
            source = validate_gene_set_source(path, options)
            assert validate_gene_set_source(path, options) is source
            rows = list(source.rows())
    assert sum(call.kwargs["validate"] for call in reader.call_args_list) == 1
    assert [(row.name, row.genes) for row in rows] == [("SET1", ("G1", "G3")), ("SET2", ("G2",))]
    assert rows[0].duplicate_gene_entries == 1


@pytest.mark.parametrize("row", ["G1", "G1|", "|SET1"])
def test_membership_missing_columns_and_empty_identifiers_rejected(tmp_path, row):
    path = tmp_path / "memberships.txt"
    path.write_text(row + "\n")
    options = GeneSetFormat(
        "membership", r"\|", "#", True, "reject", "deduplicate", False,
        membership_set_column=1, membership_gene_column=0,
    )
    with pytest.raises(ValueError, match="Malformed|Empty"):
        validate_gene_set_source(path, options)


def test_compressed_gene_resources_preserve_magma_vs_gcta_openers(tmp_path):
    locations = tmp_path / "genes.loc.gz"
    with gzip.open(locations, "wt") as handle:
        handle.write("ENSG1 1 10 20 + GENE1\n")
    assert read_gene_locations(locations, load_configuration().modules.magma)["ENSG1"] == (
        "1", 10, 20, "+", "GENE1",
    )
    path = tmp_path / "sets.gmt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("SET\tdesc\tGENE1\n")
    options = GeneSetFormat("gmt", r"\s+", "#", True, "drop", "deduplicate", False)
    assert list(validate_gene_set_source(path, options).rows())[0].genes == ("GENE1",)
    with pytest.raises(ResourcePreparationError):
        scan_gmt_pathways(path, "deduplicate")


@pytest.mark.parametrize("input_format,contents", [
    ("gmt", "SET\tdesc\tG1\tG2\n"),
    ("magma", "SET G1 G2\n"),
    ("membership", "SET G1\nSET G2\n"),
])
def test_magma_first_use_and_cache_hit_each_read_membership_data_only_once(tmp_path, input_format, contents):
    class Logger:
        def record(self, *_args, **_kwargs):
            pass

    path = tmp_path / "gene_sets.txt"
    path.write_text(contents)
    module = load_configuration().modules.magma
    with patch.object(gene_sets, "open_text", wraps=gene_sets.open_text) as opener:
        with patch.object(gene_sets, "_iter_gene_set_rows", wraps=gene_sets._iter_gene_set_rows) as reader:
            with InputValidationSession():
                first = parse_gene_set_file(path, module, Logger(), input_format=input_format)
                assert opener.call_count == 1
                second = parse_gene_set_file(path, module, Logger(), input_format=input_format)
                assert opener.call_count == 2
    assert [call.kwargs["validate"] for call in reader.call_args_list] == [True, False]
    assert first.equals(second)


def test_gene_set_first_pass_consumer_must_finish_before_publishing(tmp_path):
    path = tmp_path / "sets.gmt"
    path.write_text("SET\tdesc\tG1\n")
    options = GeneSetFormat("gmt", r"\s+", "#", True, "drop", "deduplicate", False)

    def mutate_during_read(_row):
        path.write_text("SET\tdesc\tG2\n")

    with pytest.raises(MagmaError, match="changed after pipeline preflight"):
        validate_gene_set_source(path, options, error_type=MagmaError, on_row=mutate_during_read)


@pytest.mark.parametrize("resource", ["gene_coordinates", "snp_sets", "fam"])
def test_cached_resource_mutation_retains_callers_domain_error(tmp_path, resource):
    path = tmp_path / "resource.txt"
    if resource == "gene_coordinates":
        path.write_text("1 10 20 GENE1\n")

        def validate():
            return read_gene_coordinates(
                path, column_roles=("chromosome", "start", "end", "gene"), error_type=GctaGeneError,
            )
    elif resource == "snp_sets":
        path.write_text("SET\nrs1\nEND\n")

        def validate():
            return snp_sets.validate_fastbat_set_list(path, error_type=GctaGeneError)
    else:
        path.write_text("F I 0 0 1 -9\n")

        def validate():
            return count_plink_samples(path, error_type=GctaGeneError)
    with InputValidationSession():
        validate()
        path.write_text(path.read_text() + "\n")
        with pytest.raises(GctaGeneError, match="changed after pipeline preflight"):
            validate()
