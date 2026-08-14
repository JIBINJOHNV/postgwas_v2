import argparse
import gzip
from pathlib import Path
import subprocess

import polars as pl
import pytest
import yaml

from postgwas.config import load_configuration
from postgwas.core.errors import MissingRequiredArgumentsError
from postgwas.modules.ld_clumping.cli import build_parser
from postgwas.modules.ld_clumping import ld_prune_standard
from postgwas.modules.ld_clumping.service import (
    LDClumpingError,
    LDClumpingPreflight,
    _validate_reference_contract,
    resolve_ld_clumping_configuration,
    run_ld_clump_direct,
    validate_ld_clumping_configuration,
)


def _vcf(path: Path) -> Path:
    path.write_bytes(b"vcf")
    return path


def _manifest(module, **overrides):
    document = {
        "format_version": module.reference.format_version,
        "genome_build": module.genome_build.value,
        "populations": [module.population.value],
        "orientation": module.reference.orientation,
        "file_pattern": module.reference.file_pattern,
        "columns": module.reference.columns,
        "window_kb": module.window_kb,
        "minimum_r2": min(module.clump_r2, module.lead_r2),
        "allele_order_preserved": True,
        "plink_version": "PLINK v1.90",
    }
    document.update(overrides)
    return document


def test_module_yaml_controls_methods_and_thresholds_without_cli_shadowing(tmp_path):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    run_config = tmp_path / "ld_clumping.yaml"
    run_config.write_text(
        "methods: [region]\n"
        "lead_pvalue: 1.0e-6\n"
        "candidate_pvalue: 0.01\n"
        "clump_r2: 0.3\n"
        "lead_r2: 0.05\n",
        encoding="utf-8",
    )
    args = build_parser().parse_args([
        "--run-config", str(run_config),
        "--vcf", str(vcf),
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path / "output"),
    ])

    module = resolve_ld_clumping_configuration(args).modules.ld_clumping

    assert module.methods == ["region"]
    assert module.lead_pvalue == 1e-6
    assert module.candidate_pvalue == 0.01
    assert module.clump_r2 == 0.3
    assert module.lead_r2 == 0.05


def test_standard_method_reports_exact_missing_ld_folder_before_external_tools(tmp_path):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--clumping-methods", "standard",
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path / "output"),
    ])

    with pytest.raises(MissingRequiredArgumentsError) as captured:
        validate_ld_clumping_configuration(args)

    assert str(captured.value) == (
        "Required argument not provided: --ld-folder. Provide --ld-folder VALUE "
        "or set modules.ld_clumping.reference.directory in the run configuration."
    )


def test_reference_manifest_requires_matching_build_and_symmetric_indexing(tmp_path):
    module = load_configuration().modules.ld_clumping
    manifest_path = tmp_path / module.reference.manifest_filename
    manifest_path.write_text(
        yaml.safe_dump(_manifest(module), sort_keys=False), encoding="utf-8",
    )

    observed = _validate_reference_contract(module, tmp_path)

    assert observed.orientation == "symmetric_first_endpoint"
    manifest_path.write_text(
        yaml.safe_dump(_manifest(module, genome_build="GRCh38"), sort_keys=False),
        encoding="utf-8",
    )
    with pytest.raises(LDClumpingError, match="genome_build=GRCh38"):
        _validate_reference_contract(module, tmp_path)
    manifest_path.write_text(
        yaml.safe_dump(
            _manifest(module, plink_version="PLINK v2.00"), sort_keys=False,
        ),
        encoding="utf-8",
    )
    with pytest.raises(LDClumpingError, match="PLINK 1.9"):
        _validate_reference_contract(module, tmp_path)


def test_preflight_rejects_vcf_and_ld_clumping_build_mismatch(monkeypatch, tmp_path):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--clumping-methods", "region",
        "--genome-build", "GRCh37",
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path / "output"),
    ])
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.resolve_executable",
        lambda value, *args, **kwargs: str(value),
    )
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.read_vcf_header",
        lambda *args, **kwargs: (
            "##genome_build=GRCh38\n"
            "##FORMAT=<ID=ES,Number=1,Type=Float,Description=\"effect\">\n"
            "##FORMAT=<ID=SE,Number=1,Type=Float,Description=\"se\">\n"
            "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"af\">\n"
            "##FORMAT=<ID=LP,Number=1,Type=Float,Description=\"lp\">\n"
            "##INFO=<ID=EUR_LDblock,Number=1,Type=String,Description=\"ld\">\n"
            "#CHROM\n"
        ),
    )

    with pytest.raises(LDClumpingError, match="GRCh38.*does not match.*GRCh37"):
        validate_ld_clumping_configuration(args)


def test_selected_method_failure_is_fatal_and_changed_outputs_are_quarantined(
    monkeypatch, tmp_path,
):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    output = tmp_path / "output"
    reference = tmp_path / "reference"
    reference.mkdir()
    configuration = load_configuration(cli_overrides={
        "modules.ld_clumping.inputs.vcf": str(vcf),
        "modules.ld_clumping.inputs.dataset_id": "STUDY",
        "modules.ld_clumping.output_directory": str(output),
        "modules.ld_clumping.reference.directory": str(reference),
    })
    preflight = LDClumpingPreflight(
        configuration=configuration,
        vcf=vcf,
        output_directory=output,
        dataset_id="STUDY",
        bcftools="bcftools",
        tabix="tabix",
        reference_directory=reference,
        reference_manifest=None,
    )
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.validate_ld_clumping_configuration",
        lambda *args, **kwargs: preflight,
    )

    def region(**kwargs):
        path = output / "STUDY_LDpruned_EUR.tsv"
        path.write_text("region output\n", encoding="utf-8")
        return {"ldpruned_file": str(path)}

    def standard(**kwargs):
        path = output / "STUDY_formatted.tsv"
        path.write_text("incomplete standard output\n", encoding="utf-8")
        raise RuntimeError("standard failed")

    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.ld_clump_by_regions", region,
    )
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.ld_clump_standard", standard,
    )

    with pytest.raises(RuntimeError, match="standard failed"):
        run_ld_clump_direct(argparse.Namespace(), configuration=configuration)

    assert not (output / "STUDY_LDpruned_EUR.tsv").exists()
    assert not (output / "STUDY_formatted.tsv").exists()
    assert list(output.glob("STUDY_LDpruned_EUR.tsv.partial.*"))
    assert list(output.glob("STUDY_formatted.tsv.partial.*"))


def test_symmetric_ld_query_assigns_partner_when_index_is_first_endpoint(monkeypatch):
    index_id = "1_100_A_G"
    partner_id = "1_200_C_T"

    monkeypatch.setattr(
        ld_prune_standard.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0,
            stdout=f"1\t100\t{index_id}\t1\t200\t{partner_id}\t0.8\n",
            stderr="",
        ),
    )

    partners = ld_prune_standard.get_ld_partners(
        "1", 100, index_id, "EUR_chr1.ld.gz", 0.6,
        tabix_bin="tabix", window_kb=250, missing_index_action="error",
    )

    assert partners["uniq_id"].to_list() == [index_id, partner_id]


def test_fully_resolved_ld_query_does_not_reload_configuration(monkeypatch):
    monkeypatch.setattr(
        ld_prune_standard,
        "load_configuration",
        lambda: (_ for _ in ()).throw(AssertionError("configuration reloaded")),
    )
    monkeypatch.setattr(
        ld_prune_standard.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0,
            stdout="1\t100\t1_100_A_G\t1\t200\t1_200_C_T\t0.8\n",
            stderr="",
        ),
    )

    partners = ld_prune_standard.get_ld_partners(
        "1", 100, "1_100_A_G", "EUR_chr1.ld.gz", 0.6,
        tabix_bin="tabix", window_kb=250, missing_index_action="error",
    )

    assert partners["uniq_id"].to_list() == ["1_100_A_G", "1_200_C_T"]


def test_missing_index_defaults_to_error_instead_of_false_independence(monkeypatch):
    monkeypatch.setattr(
        ld_prune_standard.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0, stdout="", stderr="",
        ),
    )

    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match="No LD rows were returned.*action=error",
    ):
        ld_prune_standard.get_ld_partners(
            "1", 100, "1_100_A_G", "EUR_chr1.ld.gz", 0.6,
            tabix_bin="tabix", window_kb=250, missing_index_action="error",
        )


def test_configured_window_excludes_distant_ld_pairs(monkeypatch):
    index_id = "1_100_A_G"
    distant_id = "1_400000_C_T"
    monkeypatch.setattr(
        ld_prune_standard.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0,
            stdout=f"1\t100\t{index_id}\t1\t400000\t{distant_id}\t0.9\n",
            stderr="",
        ),
    )

    partners = ld_prune_standard.get_ld_partners(
        "1", 100, index_id, "EUR_chr1.ld.gz", 0.6,
        tabix_bin="tabix", window_kb=250, missing_index_action="error",
    )

    assert partners["uniq_id"].to_list() == [index_id]


def test_no_significant_chromosome_does_not_require_an_ld_file(tmp_path):
    gwas = ld_prune_standard.add_canonical_ids(pl.DataFrame({
        "chrcol": ["X"],
        "poscol": [100],
        "neacol": ["A"],
        "eacol": ["G"],
        "rsIDcol": ["rs1"],
        "pcol": [0.5],
        "becol": [-0.2],
        "secol": [0.1],
        "eafcol": [0.3],
    }))

    result = ld_prune_standard.process_chromosome(
        "X", gwas, tmp_path, "EUR", 5e-8, 0.6, 0.1, 250000,
    )

    assert result["progress"]["status"] == "no_significant_variants"


def test_canonical_matching_id_does_not_flip_alleles_or_effects():
    result = ld_prune_standard.add_canonical_ids(pl.DataFrame({
        "chrcol": ["chr01"],
        "poscol": [100],
        "neacol": ["G"],
        "eacol": ["A"],
        "rsIDcol": ["rs1"],
        "pcol": [1e-9],
        "becol": [-0.2],
        "secol": [0.1],
        "eafcol": [0.3],
    }))

    row = result.row(0, named=True)
    assert row["uniq_id"] == "1_100_A_G"
    assert row["chrcol"] == "1"
    assert (row["neacol"], row["eacol"], row["becol"], row["eafcol"]) == (
        "G", "A", -0.2, 0.3,
    )


def test_conflicting_effect_allele_orientation_is_not_silently_collapsed():
    gwas = pl.DataFrame({
        "chrcol": ["1", "chr1"],
        "poscol": [100, 100],
        "neacol": ["G", "A"],
        "eacol": ["A", "G"],
        "rsIDcol": ["rs1", "rs1"],
        "pcol": [1e-9, 1e-9],
        "becol": [-0.2, 0.2],
        "secol": [0.1, 0.1],
        "eafcol": [0.3, 0.7],
    })

    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match="conflicting effect-allele orientation.*flipping is prohibited",
    ):
        ld_prune_standard.add_canonical_ids(gwas)


def test_reference_preparation_always_preserves_allele_order_and_symmetrizes():
    script = Path("tools/resource_preparation/ld_file_preparation.sh")
    text = script.read_text(encoding="utf-8")

    assert "--keep-allele-order" in text
    assert "orientation: symmetric_first_endpoint" in text
    assert "print $4, $5, $6, $1, $2, $3, $7" in text
    subprocess.run(["bash", "-n", str(script)], check=True)


def test_reference_preparation_rejects_invalid_ranges_before_plink(tmp_path):
    script = Path("tools/resource_preparation/ld_file_preparation.sh")
    fileset = tmp_path / "reference"
    for suffix in ("bed", "bim", "fam"):
        (tmp_path / f"reference.{suffix}").write_text("input\n", encoding="utf-8")

    completed = subprocess.run(
        [
            "bash", str(script),
            "--bfile", str(fileset),
            "--output-dir", str(tmp_path / "output"),
            "--population", "EUR",
            "--genome-build", "GRCh37",
            "--chromosomes", "1",
            "--plink", "printf",
            "--bgzip", "printf",
            "--tabix", "printf",
            "--window-kb", "250",
            "--ld-window-variants", "99999",
            "--minimum-r2", "1.1",
            "--minimum-maf", "0.00001",
            "--threads", "8",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 2
    assert "Invalid --minimum-r2" in completed.stderr
    assert not (tmp_path / "output").exists()


def test_reference_preparation_executes_the_exact_allele_safe_plink_contract(
    tmp_path,
):
    script = Path("tools/resource_preparation/ld_file_preparation.sh")
    fileset = tmp_path / "reference"
    for suffix in ("bed", "bim", "fam"):
        (tmp_path / f"reference.{suffix}").write_text("input\n", encoding="utf-8")

    fake_plink = tmp_path / "plink"
    fake_plink.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "if [[ ${1:-} == --version ]]; then echo 'PLINK v1.90b-test'; exit 0; fi\n"
        "out=''\n"
        "args=(\"$@\")\n"
        "for ((i=0; i < ${#args[@]}; i++)); do\n"
        "  if [[ ${args[$i]} == --out ]]; then out=${args[$((i + 1))]}; fi\n"
        "done\n"
        "printf '%s\\n' \"$@\" > \"${out}.args\"\n"
        "printf 'CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\\n1 100 1_100_A_G 1 200 1_200_C_T 0.8\\n' | gzip -c > \"${out}.ld.gz\"\n",
        encoding="utf-8",
    )
    fake_bgzip = tmp_path / "bgzip"
    fake_bgzip.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\ngzip -c\n",
        encoding="utf-8",
    )
    fake_tabix = tmp_path / "tabix"
    fake_tabix.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\nlast=${!#}\nprintf 'index\\n' > \"${last}.tbi\"\n",
        encoding="utf-8",
    )
    for executable in (fake_plink, fake_bgzip, fake_tabix):
        executable.chmod(0o755)

    output = tmp_path / "output"
    subprocess.run(
        [
            "bash", str(script),
            "--bfile", str(fileset),
            "--output-dir", str(output),
            "--population", "EUR",
            "--genome-build", "GRCh37",
            "--chromosomes", "1",
            "--plink", str(fake_plink),
            "--bgzip", str(fake_bgzip),
            "--tabix", str(fake_tabix),
            "--window-kb", "250",
            "--ld-window-variants", "99999",
            "--minimum-r2", "0.05",
            "--minimum-maf", "0.00001",
            "--threads", "8",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    arguments = (output / "EUR_chr1.args").read_text(encoding="utf-8").splitlines()
    assert arguments == [
        "--bfile", str(fileset),
        "--chr", "1",
        "--keep-allele-order",
        "--r2", "gz",
        "--ld-window", "99999",
        "--ld-window-kb", "250",
        "--ld-window-r2", "0.05",
        "--maf", "0.00001",
        "--threads", "8",
        "--out", str(output / "EUR_chr1"),
    ]
    with gzip.open(output / "EUR_chr1.ld.gz", "rt", encoding="utf-8") as handle:
        pairs = [line.rstrip("\n").split("\t") for line in handle]
    assert pairs == [
        ["1", "100", "1_100_A_G", "1", "200", "1_200_C_T", "0.8"],
        ["1", "200", "1_200_C_T", "1", "100", "1_100_A_G", "0.8"],
    ]
    manifest = yaml.safe_load(
        (output / "ld_reference.yaml").read_text(encoding="utf-8")
    )
    assert manifest["orientation"] == "symmetric_first_endpoint"
    assert manifest["allele_order_preserved"] is True
    assert manifest["plink_version"] == "PLINK v1.90b-test"
    assert (output / "EUR_chr1.ld.gz.tbi").is_file()
