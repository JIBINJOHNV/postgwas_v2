import subprocess
import stat
import sys
from pathlib import Path

import yaml


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
    mapping = (output / "gene_variant_mapping.tsv").read_text(encoding="utf-8")
    assert "PATH_A\tGENE1\t1:190:A:G" in mapping
    assert "PATH_A\tGENE2\t1:190:A:G" in mapping
    assert "PATH_B\tGENE2\tcustom_id" in mapping
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
