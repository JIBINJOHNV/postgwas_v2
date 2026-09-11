import argparse
import errno
import subprocess
import stat
import sys
from pathlib import Path

import pytest
import polars as pl
import yaml

from postgwas.config import load_configuration
from postgwas.core.resource_preparation import ResourcePreparationError
from postgwas.modules.gcta_gene import pathway_sets
from postgwas.modules.gcta_gene.pathway_sets import (
    ChromosomeMappingTask,
    GeneInterval,
    _map_bim_variants,
    _mapping_worker_plan,
    prepare_resource,
)


SCRIPT = (
    Path(__file__).parents[1]
    / "tools"
    / "resource_preparation"
    / "gcta_gene_list.py"
)
GMT_SCRIPT = (
    Path(__file__).parents[1]
    / "tools"
    / "resource_preparation"
    / "gmt_to_gcta_fastbat_sets.py"
)


def _command(source: Path, output: Path) -> list[str]:
    return [
        sys.executable,
        str(SCRIPT),
        "--input", str(source),
        "--output-directory", str(output),
        "--output-name", "genes.txt",
        "--resource-name", "Test GCTA genes",
        "--source-label", "test source",
        "--genome-build", "GRCh37",
        "--expected-source-columns", "6",
        "--chromosome-column", "2",
        "--start-column", "3",
        "--end-column", "4",
        "--gene-column", "6",
        "--allowed-chromosomes", "1", "X", "Y",
    ]


def _gmt_command(gmt: Path, genes: Path, bim: Path, output: Path) -> list[str]:
    return [
        sys.executable,
        str(GMT_SCRIPT),
        "--gmt", str(gmt),
        "--gene-list", str(genes),
        "--bim", str(bim),
        "--output-directory", str(output),
        "--output-name", "pathways.set",
        "--resource-name", "Test fastBAT pathways",
        "--genome-build", "GRCh37",
        "--gene-window-kb", "0",
        "--allowed-chromosomes", "1", "X",
        "--chromosome-label-policy", "exact",
        "--duplicate-gene-policy", "deduplicate",
        "--unmapped-gene-policy", "report",
        "--empty-pathway-policy", "omit",
    ]


def _prepare_pathway_resource(
    gmt: Path,
    genes: Path,
    bim: Path,
    output: Path,
    *,
    audit_level: str = "normalized",
    analyzable_variant_ids: set[str] | None = None,
    analysis_variant_source: Path | None = None,
    maximum_set_variants: int | None = None,
    minimum_free_disk_gb: float = 0,
    minimum_gene_id_overlap_fraction: float | None = None,
    excluded_gene_ids: set[str] | None = None,
    analysis_scope: dict | None = None,
    allowed_chromosomes: list[str] | None = None,
):
    configuration = load_configuration()
    set_annotation = configuration.modules.gcta_gene.set_annotation
    conversion = set_annotation.conversion
    names = conversion.output_names
    return prepare_resource(
        argparse.Namespace(
            gmt=gmt,
            gene_list=genes,
            bim=bim,
            output_directory=output,
            output_name="pathways.set",
            resource_name="Test fastBAT pathways",
            genome_build="GRCh37",
            gene_window_kb=0,
            gene_columns=list(
                configuration.modules.gcta_gene.gene_annotation.columns
            ),
            allowed_chromosomes=allowed_chromosomes or ["1", "X"],
            chromosome_label_policy="exact",
            duplicate_gene_policy="deduplicate",
            unmapped_gene_policy="report",
            empty_pathway_policy="omit",
            pathway_mapping_name=names.pathway_mapping,
            pathway_gene_mapping_name=names.pathway_gene_mapping,
            gene_variant_mapping_name=names.gene_variant_mapping,
            expanded_mapping_name=names.expanded_mapping,
            unmapped_genes_name=names.unmapped_genes,
            manifest_name=names.manifest,
            readme_name=names.readme,
            checksums_name=names.checksums,
        ),
        mapping_workers=1,
        mapping_memory_gb=1,
        worker_memory_multiplier=conversion.parallelism.worker_memory_multiplier,
        minimum_gene_id_overlap_fraction=(
            minimum_gene_id_overlap_fraction
            if minimum_gene_id_overlap_fraction is not None
            else conversion.minimum_gene_id_overlap_fraction
        ),
        analyzable_variant_ids=analyzable_variant_ids,
        analysis_variant_source=analysis_variant_source,
        excluded_gene_ids=excluded_gene_ids,
        analysis_scope=analysis_scope,
        maximum_set_variants=(
            maximum_set_variants
            if maximum_set_variants is not None
            else set_annotation.maximum_set_variants
        ),
        oversized_set_policy=set_annotation.oversized_set_policy,
        audit_level=audit_level,
        audit_format=conversion.audit.format,
        audit_compression=conversion.audit.compression,
        audit_batch_rows=2,
        minimum_free_disk_gb=minimum_free_disk_gb,
        disk_estimation_safety_factor=1,
    )


def test_coordinate_reference_threshold_is_inclusive_and_fails_before_mapping(
    tmp_path,
):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text(
        "PATH_A\tdescription\tGENE1\tGENE2\n"
        "PATH_B\tdescription\tGENE2\n",
        encoding="utf-8",
    )
    genes = tmp_path / "genes.txt"
    genes.write_text(
        "1\t100\t200\tGENE1\n1\t300\t400\tGENE_COORDINATE_ONLY\n",
        encoding="utf-8",
    )
    bim = tmp_path / "reference.bim"
    bim.write_text("1\trs1\t0\t150\tA\tG\n", encoding="utf-8")

    manifest = _prepare_pathway_resource(
        gmt,
        genes,
        bim,
        tmp_path / "passing",
        minimum_gene_id_overlap_fraction=0.5,
    )

    assert manifest["validation"]["input_unique_genes"] == 2
    assert manifest["validation"]["genes_with_coordinates"] == 1
    assert manifest["validation"]["gmt_gene_mappability_fraction"] == 0.5
    assert manifest["validation"]["coordinate_reference_coverage_fraction"] == 0.5
    assert (
        manifest["validation"]["coordinate_reference_genes_absent_from_gmt"]
        == 1
    )
    assert manifest["validation"]["unmapped_unique_genes"] == 1
    assert manifest["validation"]["pathways_with_complete_gene_id_mapping"] == 0
    assert manifest["validation"]["pathways_with_partial_gene_id_mapping"] == 1
    assert manifest["validation"]["pathways_without_gene_id_mapping"] == 1
    mapping = pl.read_csv(
        tmp_path / "passing" / "pathway_mapping.tsv", separator="\t"
    )
    assert mapping["gene_id_mappability_fraction"].to_list() == [0.5, 0.0]
    assert mapping["analyzable_gene_fraction"].to_list() == [0.5, 0.0]

    failing_output = tmp_path / "failing"
    with pytest.raises(
        ResourcePreparationError,
        match=(
            "Coordinate-reference coverage is 1/2.*configured minimum of "
            "50.01%.*GMT gene mappability is 1/2.*reported only"
        ),
    ):
        _prepare_pathway_resource(
            gmt,
            genes,
            bim,
            failing_output,
            minimum_gene_id_overlap_fraction=0.5001,
        )
    assert not failing_output.exists()


def test_focused_gmt_fails_low_coordinate_reference_coverage(
    tmp_path,
):
    gmt = tmp_path / "focused.gmt"
    gmt.write_text("FOCUSED\tdescription\tGENE1\n", encoding="utf-8")
    genes = tmp_path / "genes.txt"
    genes.write_text(
        "1\t100\t200\tGENE1\n"
        "1\t300\t400\tGENE2\n"
        "1\t500\t600\tGENE3\n",
        encoding="utf-8",
    )
    bim = tmp_path / "reference.bim"
    bim.write_text("1\trs1\t0\t150\tA\tG\n", encoding="utf-8")

    output = tmp_path / "sets"
    with pytest.raises(
        ResourcePreparationError,
        match=(
            "Coordinate-reference coverage is 1/3.*below the configured "
            "minimum of 50.00%.*GMT gene mappability is 1/1"
        ),
    ):
        _prepare_pathway_resource(
            gmt,
            genes,
            bim,
            output,
            minimum_gene_id_overlap_fraction=0.5,
        )
    assert not output.exists()


def test_low_gmt_mappability_passes_when_coordinate_reference_coverage_passes(
    tmp_path,
):
    gmt = tmp_path / "broad.gmt"
    gmt.write_text(
        "PATH_A\tdescription\tGENE1\tMISSING1\tMISSING2\tMISSING3\n"
        "PATH_B\tdescription\tGENE2\n",
        encoding="utf-8",
    )
    genes = tmp_path / "genes.txt"
    genes.write_text(
        "1\t100\t200\tGENE1\n1\t300\t400\tGENE2\n",
        encoding="utf-8",
    )
    bim = tmp_path / "reference.bim"
    bim.write_text(
        "1\trs1\t0\t150\tA\tG\n1\trs2\t0\t350\tC\tT\n",
        encoding="utf-8",
    )

    manifest = _prepare_pathway_resource(
        gmt,
        genes,
        bim,
        tmp_path / "sets",
        minimum_gene_id_overlap_fraction=0.5,
    )

    validation = manifest["validation"]
    assert validation["gmt_gene_mappability_fraction"] == pytest.approx(0.4)
    assert validation["coordinate_reference_coverage_fraction"] == 1
    assert validation["coordinate_reference_genes_absent_from_gmt"] == 0


def test_pathway_scope_exclusions_are_distinct_from_unmapped_genes(tmp_path):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text(
        "PATH_A\tdescription\tGENE1\tGENE_MHC\tGENE_MISSING\n",
        encoding="utf-8",
    )
    genes = tmp_path / "genes.txt"
    genes.write_text(
        "1\t100\t200\tGENE1\n6\t29000000\t29100000\tGENE_MHC\n",
        encoding="utf-8",
    )
    bim = tmp_path / "reference.bim"
    bim.write_text(
        "1\trs1\t0\t150\tA\tG\n6\trs_mhc\t0\t29050000\tA\tG\n",
        encoding="utf-8",
    )
    analysis_scope = {
        "mhc_policy": "exclude_genes",
        "exclude_mhc_snps": False,
        "exclude_mhc_genes": True,
        "mhc_region": {"chromosome": "6", "start": 28477797, "end": 33448354},
        "mhc_source": "genome_resources",
        "exclude_chromosomes": [],
    }

    manifest = _prepare_pathway_resource(
        gmt,
        genes,
        bim,
        tmp_path / "sets",
        minimum_gene_id_overlap_fraction=0.5,
        excluded_gene_ids={"GENE_MHC"},
        analysis_scope=analysis_scope,
        allowed_chromosomes=["1", "6"],
    )

    assert manifest["schema_version"] == "gcta_fastbat_pathway_resource.v6"
    assert manifest["validation"]["genes_with_coordinates"] == 2
    assert manifest["validation"]["genes_eligible_for_analysis"] == 1
    assert manifest["validation"]["genes_excluded_by_analysis_scope"] == 1
    assert manifest["validation"]["unmapped_unique_genes"] == 1
    assert manifest["policies"]["analysis_scope"] == analysis_scope
    assert (tmp_path / "sets" / "pathways.set").read_text(
        encoding="utf-8"
    ) == "PATH_A\nrs1\nEND\n\n"
    assert (tmp_path / "sets" / "unmapped_genes.tsv").read_text(
        encoding="utf-8"
    ).splitlines() == ["pathway\tgene", "PATH_A\tGENE_MISSING"]


def test_resource_script_writes_exact_headerless_projection_and_provenance(tmp_path):
    source = tmp_path / "source.loc"
    source.write_text(
        "79501 1 69091 70008 + OR4F5\n9085 Y 25622115 25634745 + CDY1\n",
        encoding="utf-8",
    )
    output = tmp_path / "resource"

    result = subprocess.run(
        _command(source, output), capture_output=True, text=True, check=False,
    )

    assert result.returncode == 0, result.stderr
    assert (output / "genes.txt").read_text(encoding="utf-8") == (
        "1\t69091\t70008\tOR4F5\nY\t25622115\t25634745\tCDY1\n"
    )
    manifest = yaml.safe_load(
        (output / "resource_manifest.yaml").read_text(encoding="utf-8")
    )
    assert manifest["resource"]["methods"] == ["fastbat_gene", "mbat_combo"]
    assert manifest["resource"]["genome_build"] == "GRCh37"
    assert manifest["source"]["columns"] == {
        "chromosome": 2, "start": 3, "end": 4, "gene": 6,
    }
    assert manifest["validation"]["output_rows"] == 2
    assert manifest["validation"]["excluded_rows"] == 0
    assert (output / "README.md").is_file()
    assert (output / "SHA256SUMS").is_file()
    assert stat.S_IMODE(output.stat().st_mode) == stat.S_IMODE(tmp_path.stat().st_mode)


def test_resource_script_rejects_duplicate_gene_ids_without_output(tmp_path):
    source = tmp_path / "source.loc"
    source.write_text(
        "1 1 10 20 + GENE1\n2 X 30 40 - GENE1\n",
        encoding="utf-8",
    )
    output = tmp_path / "resource"

    result = subprocess.run(
        _command(source, output), capture_output=True, text=True, check=False,
    )

    assert result.returncode == 2
    assert "repeats gene identifier 'GENE1'" in result.stderr
    assert not output.exists()


def test_resource_script_rejects_undeclared_chromosome(tmp_path):
    source = tmp_path / "source.loc"
    source.write_text("1 2 10 20 + GENE1\n", encoding="utf-8")

    result = subprocess.run(
        _command(source, tmp_path / "resource"),
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 2
    assert "disallowed chromosome '2'" in result.stderr


def test_resource_script_refuses_existing_output_directory(tmp_path):
    source = tmp_path / "source.loc"
    source.write_text("1 1 10 20 + GENE1\n", encoding="utf-8")
    output = tmp_path / "resource"
    output.mkdir()
    marker = output / "curated.txt"
    marker.write_text("preserve me\n", encoding="utf-8")

    result = subprocess.run(
        _command(source, output), capture_output=True, text=True, check=False,
    )

    assert result.returncode == 2
    assert "output directory already exists" in result.stderr
    assert marker.read_text(encoding="utf-8") == "preserve me\n"


def test_gmt_converter_copies_exact_bim_ids_and_reports_mapping(tmp_path):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text(
        "PATH_A\tdescription A\tGENE1\tGENE2\tGENE2\n"
        "PATH_B\tdescription B\tGENE2\tMISSING\n"
        "PATH_EMPTY\tdescription C\tMISSING\n",
        encoding="utf-8",
    )
    genes = tmp_path / "genes.txt"
    genes.write_text(
        "1\t100\t200\tGENE1\n1\t180\t300\tGENE2\nX\t10\t20\tOTHER\n",
        encoding="utf-8",
    )
    bim = tmp_path / "reference.bim"
    bim.write_text(
        "1\trs1\t0\t150\tA\tG\n"
        "1\t1:190:A:G\t0\t190\tA\tG\n"
        "1\tcustom_id\t0\t250\tG\tT\n"
        "X\trsX\t0\t15\tA\tC\n",
        encoding="utf-8",
    )
    output = tmp_path / "sets"

    result = subprocess.run(
        _gmt_command(gmt, genes, bim, output),
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert (output / "pathways.set").read_text(encoding="utf-8") == (
        "PATH_A\nrs1\n1:190:A:G\ncustom_id\nEND\n\n"
        "PATH_B\n1:190:A:G\ncustom_id\nEND\n\n"
    )
    pathway_gene = pl.read_parquet(output / "pathway_gene_mapping.parquet")
    gene_variant = pl.read_parquet(output / "gene_variant_mapping.parquet")
    assert ("PATH_A", "GENE1") in set(pathway_gene.iter_rows())
    assert ("PATH_A", "GENE2") in set(pathway_gene.iter_rows())
    assert ("PATH_B", "GENE2") in set(pathway_gene.iter_rows())
    assert ("GENE1", "1:190:A:G") in set(gene_variant.iter_rows())
    assert ("GENE2", "1:190:A:G") in set(gene_variant.iter_rows())
    assert ("GENE2", "custom_id") in set(gene_variant.iter_rows())
    assert not (output / "pathway_gene_variant_mapping.parquet").exists()
    assert "PATH_EMPTY\tMISSING" in (
        output / "unmapped_genes.tsv"
    ).read_text(encoding="utf-8")
    manifest = yaml.safe_load(
        (output / "resource_manifest.yaml").read_text(encoding="utf-8")
    )
    assert manifest["policies"]["variant_identifier_policy"] == (
        "copy PLINK BIM column 2 verbatim; never infer rsIDs"
    )
    assert manifest["validation"]["pathways_written"] == 2
    assert manifest["validation"]["pathways_omitted"] == 1
    assert manifest["validation"]["unique_variants_written"] == 3
    assert manifest["validation"][
        "overlapping_gene_variant_memberships_collapsed"
    ] == 1
    assert manifest["generation"]["mapping_algorithm"] == (
        "chromosome_bisect_compact_variant_index_cache"
    )
    assert manifest["policies"]["audit"]["level"] == "normalized"
    assert manifest["policies"]["maximum_set_variants_applied"] is False
    assert manifest["validation"]["expanded_audit_rows"] == 0
    assert manifest["disk_preflight"]["passed"] is True


def test_analyzable_intersection_precedes_pathway_size_validation(tmp_path):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("PATH_A\tdescription\tGENE1\n", encoding="utf-8")
    genes = tmp_path / "genes.txt"
    genes.write_text("1\t100\t200\tGENE1\n", encoding="utf-8")
    bim = tmp_path / "reference.bim"
    bim.write_text(
        "1\tused\t0\t150\tA\tG\n"
        "1\treference_only\t0\t160\tC\tT\n",
        encoding="utf-8",
    )
    analysis_source = tmp_path / "input.ma"
    analysis_source.write_text("used\n", encoding="utf-8")
    output = tmp_path / "sets"

    manifest = _prepare_pathway_resource(
        gmt,
        genes,
        bim,
        output,
        analyzable_variant_ids={"used"},
        analysis_variant_source=analysis_source,
        maximum_set_variants=1,
    )

    assert (output / "pathways.set").read_text(encoding="utf-8").split() == [
        "PATH_A", "used", "END",
    ]
    assert manifest["validation"]["pathways_omitted_oversized"] == 0
    assert manifest["validation"]["analyzable_bim_variants_cached"] == 1
    assert manifest["policies"]["variant_universe"] == "gwas_bim_intersection"
    assert manifest["policies"]["maximum_set_variants_applied"] is True


def test_standalone_reference_resource_defers_gwas_specific_set_limit(tmp_path):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("PATH_A\tdescription\tGENE1\n", encoding="utf-8")
    genes = tmp_path / "genes.txt"
    genes.write_text("1\t100\t200\tGENE1\n", encoding="utf-8")
    bim = tmp_path / "reference.bim"
    bim.write_text(
        "1\trs1\t0\t150\tA\tG\n1\trs2\t0\t160\tC\tT\n",
        encoding="utf-8",
    )
    output = tmp_path / "sets"

    manifest = _prepare_pathway_resource(
        gmt, genes, bim, output, maximum_set_variants=1,
    )

    assert (output / "pathways.set").read_text(encoding="utf-8").split() == [
        "PATH_A", "rs1", "rs2", "END",
    ]
    assert manifest["validation"]["pathways_omitted_oversized"] == 0
    assert manifest["policies"]["variant_universe"] == "plink_bim_reference"
    assert manifest["policies"]["maximum_set_variants_applied"] is False


def test_summary_audit_omits_relationship_tables(tmp_path):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("PATH_A\tdescription\tGENE1\n", encoding="utf-8")
    genes = tmp_path / "genes.txt"
    genes.write_text("1\t100\t200\tGENE1\n", encoding="utf-8")
    bim = tmp_path / "reference.bim"
    bim.write_text("1\trs1\t0\t150\tA\tG\n", encoding="utf-8")
    output = tmp_path / "sets"

    manifest = _prepare_pathway_resource(
        gmt, genes, bim, output, audit_level="summary",
    )

    assert manifest["outputs"]["pathway_gene_mapping"] is None
    assert manifest["outputs"]["gene_variant_mapping"] is None
    assert manifest["outputs"]["expanded_mapping"] is None
    assert not list(output.glob("*.parquet"))


def test_expanded_audit_is_explicit_bounded_parquet_output(tmp_path):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text(
        "PATH_A\tdescription\tGENE1\tGENE2\n",
        encoding="utf-8",
    )
    genes = tmp_path / "genes.txt"
    genes.write_text(
        "1\t100\t200\tGENE1\n1\t150\t250\tGENE2\n",
        encoding="utf-8",
    )
    bim = tmp_path / "reference.bim"
    bim.write_text(
        "1\trs1\t0\t175\tA\tG\n1\trs2\t0\t225\tC\tT\n",
        encoding="utf-8",
    )
    output = tmp_path / "sets"

    manifest = _prepare_pathway_resource(
        gmt, genes, bim, output, audit_level="expanded",
    )

    expanded = pl.read_parquet(output / "pathway_gene_variant_mapping.parquet")
    assert set(expanded.iter_rows()) == {
        ("PATH_A", "GENE1", "rs1"),
        ("PATH_A", "GENE2", "rs1"),
        ("PATH_A", "GENE2", "rs2"),
    }
    assert manifest["validation"]["expanded_audit_rows"] == 3
    assert (output / "pathway_gene_mapping.parquet").is_file()
    assert (output / "gene_variant_mapping.parquet").is_file()


def test_disk_preflight_fails_before_output_writing(tmp_path, monkeypatch):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("PATH_A\tdescription\tGENE1\n", encoding="utf-8")
    genes = tmp_path / "genes.txt"
    genes.write_text("1\t100\t200\tGENE1\n", encoding="utf-8")
    bim = tmp_path / "reference.bim"
    bim.write_text("1\trs1\t0\t150\tA\tG\n", encoding="utf-8")
    output = tmp_path / "sets"
    monkeypatch.setattr(
        pathway_sets.shutil,
        "disk_usage",
        lambda _path: argparse.Namespace(total=1, used=0, free=1),
    )

    with pytest.raises(ResourcePreparationError, match="insufficient disk space"):
        _prepare_pathway_resource(gmt, genes, bim, output)

    assert not output.exists()
    assert not list(tmp_path.glob(".sets.*"))


def test_disk_full_error_reports_failed_file_and_safe_restart(tmp_path, monkeypatch):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("PATH_A\tdescription\tGENE1\n", encoding="utf-8")
    genes = tmp_path / "genes.txt"
    genes.write_text("1\t100\t200\tGENE1\n", encoding="utf-8")
    bim = tmp_path / "reference.bim"
    bim.write_text("1\trs1\t0\t150\tA\tG\n", encoding="utf-8")
    output = tmp_path / "sets"

    def disk_full(_self, *_values):
        raise OSError(errno.ENOSPC, "No space left on device")

    monkeypatch.setattr(pathway_sets._ParquetRowWriter, "append", disk_full)

    with pytest.raises(ResourcePreparationError) as error:
        _prepare_pathway_resource(gmt, genes, bim, output)

    message = str(error.value)
    assert "failed output=" in message
    assert "bytes written=" in message
    assert "audit level=normalized" in message
    assert "restart from the safe preparation boundary" in message
    assert not output.exists()
    assert not list(tmp_path.glob(".sets.*"))


def test_gmt_converter_uses_bounded_parallel_chromosome_caches(
    tmp_path, monkeypatch,
):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text(
        "PATH_A\tdescription\tGENE1\tGENEX\n", encoding="utf-8",
    )
    genes = tmp_path / "genes.txt"
    genes.write_text(
        "1\t100\t200\tGENE1\nX\t10\t20\tGENEX\n", encoding="utf-8",
    )
    bim = tmp_path / "reference.bim"
    bim.write_text(
        "1\trs1\t0\t150\tA\tG\nX\trsX\t0\t15\tA\tC\n",
        encoding="utf-8",
    )
    output = tmp_path / "sets"
    names = (
        load_configuration()
        .modules.gcta_gene.set_annotation.conversion.output_names
        .model_dump(mode="python")
    )
    defaults = load_configuration()
    set_annotation = defaults.modules.gcta_gene.set_annotation
    conversion = set_annotation.conversion

    class SynchronousExecutor:
        def __init__(self, max_workers):
            self.max_workers = max_workers

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        @staticmethod
        def map(function, tasks):
            return map(function, tasks)

    monkeypatch.setattr(pathway_sets, "ProcessPoolExecutor", SynchronousExecutor)
    manifest = prepare_resource(
        argparse.Namespace(
            gmt=gmt,
            gene_list=genes,
            bim=bim,
            output_directory=output,
            output_name="pathways.set",
            resource_name="Test fastBAT pathways",
            genome_build="GRCh37",
            gene_window_kb=0,
            gene_columns=list(
                defaults.modules.gcta_gene.gene_annotation.columns
            ),
            allowed_chromosomes=["1", "X"],
            chromosome_label_policy="exact",
            duplicate_gene_policy="deduplicate",
            unmapped_gene_policy="report",
            empty_pathway_policy="omit",
            pathway_mapping_name=names["pathway_mapping"],
            pathway_gene_mapping_name=names["pathway_gene_mapping"],
            gene_variant_mapping_name=names["gene_variant_mapping"],
            expanded_mapping_name=names["expanded_mapping"],
            unmapped_genes_name=names["unmapped_genes"],
            manifest_name=names["manifest"],
            readme_name=names["readme"],
            checksums_name=names["checksums"],
        ),
        mapping_workers=2,
        mapping_memory_gb=1,
        worker_memory_multiplier=8,
        minimum_gene_id_overlap_fraction=(
            conversion.minimum_gene_id_overlap_fraction
        ),
        maximum_set_variants=set_annotation.maximum_set_variants,
        oversized_set_policy=set_annotation.oversized_set_policy,
        audit_level=conversion.audit.level,
        audit_format=conversion.audit.format,
        audit_compression=conversion.audit.compression,
        audit_batch_rows=conversion.audit.batch_rows,
        minimum_free_disk_gb=0,
        disk_estimation_safety_factor=1,
    )

    assert (output / "pathways.set").read_text(encoding="utf-8") == (
        "PATH_A\nrs1\nrsX\nEND\n\n"
    )
    assert manifest["generation"]["requested_mapping_workers"] == 2
    assert manifest["generation"]["mapping_workers"] == 2


def test_mapping_worker_plan_respects_thread_and_memory_limits(tmp_path):
    tasks = []
    for chromosome in ("1", "2", "3"):
        cache = tmp_path / (chromosome + ".cache")
        cache.write_bytes(b"x" * 1024)
        tasks.append(ChromosomeMappingTask(
            chromosome=chromosome,
            bim_cache=cache,
            intervals=(),
            result_cache=tmp_path / (chromosome + ".result"),
        ))
    two_worker_memory_gb = (2 * 1024 * 8) / float(1024**3)

    workers, metrics = _mapping_worker_plan(
        tasks,
        requested_workers=3,
        memory_gb=two_worker_memory_gb,
        worker_memory_multiplier=8,
    )

    assert workers == 2
    assert metrics["requested_mapping_workers"] == 3
    assert metrics["mapping_workers"] == 2
    assert metrics["estimated_worker_memory_bytes"] == 8192


def test_mapping_worker_plan_rejects_insufficient_memory(tmp_path):
    cache = tmp_path / "1.cache"
    cache.write_bytes(b"x" * 1024)
    task = ChromosomeMappingTask(
        chromosome="1",
        bim_cache=cache,
        intervals=(),
        result_cache=tmp_path / "1.result",
    )

    with pytest.raises(
        ResourcePreparationError,
        match="configured mapping memory is insufficient",
    ):
        _mapping_worker_plan(
            [task],
            requested_workers=1,
            memory_gb=1024 / float(1024**3),
            worker_memory_multiplier=8,
        )


def test_parallel_mapping_reports_worker_failure(tmp_path, monkeypatch):
    bim = tmp_path / "reference.bim"
    bim.write_text(
        "1\trs1\t0\t100\tA\tG\n2\trs2\t0\t100\tC\tT\n",
        encoding="utf-8",
    )
    workspace = tmp_path / "mapping"
    workspace.mkdir()

    class FailingExecutor:
        def __init__(self, max_workers):
            self.max_workers = max_workers

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        @staticmethod
        def map(_function, _tasks):
            raise RuntimeError("worker crashed")

    monkeypatch.setattr(pathway_sets, "ProcessPoolExecutor", FailingExecutor)

    with pytest.raises(
        ResourcePreparationError,
        match=r"BIM chromosome mapping failed with 2 worker\(s\): worker crashed",
    ):
        _map_bim_variants(
            bim,
            {
                "1": [GeneInterval("1", 50, 150, "GENE1")],
                "2": [GeneInterval("2", 50, 150, "GENE2")],
            },
            {"1", "2"},
            "exact",
            workspace,
            expected_variant_total=2,
            requested_workers=2,
            memory_gb=1,
            worker_memory_multiplier=8,
            bim_ids_prevalidated=False,
        )


def test_gmt_converter_rejects_duplicate_bim_ids_without_published_output(tmp_path):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("PATH_A\tdescription\tGENE1\n", encoding="utf-8")
    genes = tmp_path / "genes.txt"
    genes.write_text("1\t100\t200\tGENE1\n", encoding="utf-8")
    bim = tmp_path / "reference.bim"
    bim.write_text(
        "1\trs1\t0\t150\tA\tG\n1\trs1\t0\t160\tA\tC\n",
        encoding="utf-8",
    )
    output = tmp_path / "sets"

    result = subprocess.run(
        _gmt_command(gmt, genes, bim, output),
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 2
    assert "repeats variant ID 'rs1'" in result.stderr
    assert not output.exists()
    assert not list(tmp_path.glob(".sets.*"))


def test_gmt_converter_rejects_reserved_end_bim_identifier(tmp_path):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("PATH_A\tdescription\tGENE1\n", encoding="utf-8")
    genes = tmp_path / "genes.txt"
    genes.write_text("1\t100\t200\tGENE1\n", encoding="utf-8")
    bim = tmp_path / "reference.bim"
    bim.write_text("1\tEND\t0\t150\tA\tG\n", encoding="utf-8")
    output = tmp_path / "sets"

    result = subprocess.run(
        _gmt_command(gmt, genes, bim, output),
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 2
    assert "reserved END variant ID" in result.stderr
    assert not output.exists()
