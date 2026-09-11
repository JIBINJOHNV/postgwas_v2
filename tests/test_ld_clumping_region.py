import math
from pathlib import Path

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.ld_clumping.ld_prune_region import ld_clump_by_regions


def test_region_clumping_queries_only_the_selected_population_block(
    monkeypatch, tmp_path,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf")
    observed = {}

    def fake_extract(vcf_path, table_path, dataset_id, columns, *args, **kwargs):
        observed["columns"] = columns
        pl.DataFrame({
            "CHR": ["1"],
            "BP": [15],
            "ID": ["rs15"],
            "REF": ["A"],
            "ALT": ["G"],
            "BETA": [0.2],
            "SE": [0.1],
            "AF": [0.2],
            "LP": [9.0],
            "EUR_LDblock": ["EUR-1_10_20"],
        }).write_csv(table_path, separator="\t")
        return str(table_path)

    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.ld_prune_region.extract_vcf_table",
        fake_extract,
    )

    result = ld_clump_by_regions(
        str(vcf), str(tmp_path), "STUDY", population="EUR", bcftools="bcftools",
    )

    assert observed["columns"]["EUR_LDblock"] == "%INFO/EUR_LDblock"
    assert observed["columns"]["ID"] == "%ID"
    assert "AFR_LDblock" not in observed["columns"]
    assert "EAS_LDblock" not in observed["columns"]
    significant = pl.read_csv(result["ldpruned_sig_file"], separator="\t")
    assert significant.select("SNP", "START", "END").row(0) == (
        "1_15_A_G", 10, 20,
    )
    assert result["significant_blocks"] == 1


def test_region_clumping_preserves_x_after_autosomal_schema_inference(
    monkeypatch, tmp_path,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf")
    configuration = load_configuration().modules.ld_clumping.model_copy(deep=True)
    configuration.table.infer_schema_length = 1

    def fake_extract(_vcf, table_path, *_args, **_kwargs):
        pl.DataFrame({
            "CHR": ["1", "X"],
            "BP": [15, 25],
            "ID": ["rs15", "rsX25"],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
            "BETA": ["0.2", "-0.3"],
            "SE": ["0.1", "0.2"],
            "AF": ["0.2", "0.4"],
            "LP": ["9.0", "8.0"],
            "EUR_LDblock": ["EUR-1_10_20", "EUR-X_20_30"],
        }).write_csv(table_path, separator="\t")
        return str(table_path)

    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.ld_prune_region.extract_vcf_table",
        fake_extract,
    )

    result = ld_clump_by_regions(
        str(vcf),
        str(tmp_path),
        "STUDY",
        population="EUR",
        bcftools="bcftools",
        configuration=configuration,
    )

    significant = pl.read_csv(
        result["ldpruned_sig_file"],
        separator="\t",
        schema_overrides={"CHR": pl.String, "BETA": pl.String},
    )
    assert significant["CHR"].to_list() == ["1", "X"]
    assert significant["SNP"].to_list() == ["1_15_A_G", "X_25_C_T"]
    assert significant["BETA"].to_list() == ["0.2", "-0.3"]


@pytest.mark.parametrize("compressor", ["pigz", None])
def test_region_raw_table_is_a_readable_gzip_file(compressor, monkeypatch, tmp_path):
    """The raw table must decompress.

    polars writes bytes through its own writer, so a text-mode gzip wrapper
    truncates the deflate stream and appends the CSV uncompressed. The file
    still starts with a valid gzip header, so the corruption only surfaces when
    something tries to read it back.
    """
    import gzip

    from postgwas.config import load_configuration
    from postgwas.modules.ld_clumping import ld_prune_region

    configuration = load_configuration(cli_overrides={
        "modules.ld_clumping.remove_mhc": False,
    }).modules.ld_clumping.model_copy(deep=True)
    # Cover the external compressor and the pure-Python fallback.
    configuration.table.compressor = compressor

    def fake_extract(vcf, table_path, dataset_id, columns, bcftools, **kwargs):
        pl.DataFrame({
            "CHR": ["1"] * 500,
            "BP": [str(100 + index) for index in range(500)],
            "ID": ["rs%s" % index for index in range(500)],
            "REF": ["A"] * 500,
            "ALT": ["G"] * 500,
            "BETA": ["0.2"] * 500,
            "SE": ["0.1"] * 500,
            "AF": ["0.3"] * 500,
            "LP": ["9.0"] * 500,
            "EUR_LDblock": ["EUR_1_1_100000"] * 500,
        }).write_csv(table_path, separator="\t")
        return str(table_path)

    monkeypatch.setattr(ld_prune_region, "extract_vcf_table", fake_extract)
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf")

    result = ld_prune_region.ld_clump_by_regions(
        str(vcf), str(tmp_path / "out"), "STUDY",
        population="EUR", bcftools="bcftools", threads=2,
        configuration=configuration,
    )

    with gzip.open(result["raw_file"], "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert len(text.splitlines()) == 501
    assert text.startswith("CHR\t")


def test_region_reports_every_significant_variant_without_ldblock_annotation(
    monkeypatch, tmp_path,
):
    configuration = load_configuration().modules.ld_clumping
    lead_lp_threshold = -math.log10(configuration.lead_pvalue)
    records = []

    class CaptureLogger:
        def record(self, marker, subject, **values):
            records.append((marker, subject, values))

    def fake_extract(_vcf, table_path, *_args, **_kwargs):
        pl.DataFrame({
            "CHR": ["1", "1", "1", "2", "6"],
            "BP": [100, 200, 300, 400, 30_000_000],
            "ID": [
                "rsBlock", "rsOutside", "rsNotSignificant", "rsBlock2",
                "rsMHC",
            ],
            "REF": ["A", "C", "G", "T", "A"],
            "ALT": ["G", "T", "A", "C", "C"],
            "BETA": ["0.2", "-0.3", "0.1", "0.4", "0.5"],
            "SE": ["0.1", "0.1", "0.2", "0.1", "0.1"],
            "AF": ["0.2", "0.3", "0.4", "0.1", "0.2"],
            "LP": [
                "9.0", str(lead_lp_threshold), "6.0", "7.5", "10.0",
            ],
            "EUR_LDblock": [
                "EUR_1_50_150", None, None, "EUR_2_350_450", None,
            ],
        }).write_csv(table_path, separator="\t")
        return str(table_path)

    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.ld_prune_region.extract_vcf_table",
        fake_extract,
    )
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf")

    result = ld_clump_by_regions(
        str(vcf),
        str(tmp_path / "out"),
        "STUDY",
        population="EUR",
        bcftools="bcftools",
        include_report_tables=True,
        configuration=configuration,
        logger=CaptureLogger(),
    )

    assert result["genome_wide_significant_outside_ld_regions"] == 1
    audit = pl.read_csv(
        result["significant_outside_ld_regions_file"],
        separator="\t",
        schema_overrides={"CHR": pl.String},
    )
    assert audit.height == 1
    outside = audit.select(
        "CHR", "BP", "ID", "uniq_id", "LP", "exclusion_reason",
    ).row(0)
    assert outside[:4] == (
        "1", 200, "rsOutside", "1_200_C_T",
    )
    assert outside[4] == pytest.approx(lead_lp_threshold)
    assert outside[5] == "missing_EUR_LDblock_annotation"
    assert result["_significant_outside_ld_region_details"][
        "total_rows"
    ] == 1
    log = Path(result["log_file"]).read_text(encoding="utf-8")
    assert "genome_wide_significant_outside_ld_regions=1" in log
    assert "input_variant_id=rsOutside" in log
    assert "rsNotSignificant" not in log
    assert "rsMHC" not in log
    assert any(
        marker == "OUTPUT"
        and subject == "region_clumping"
        and values["genome_wide_significant_outside_ld_regions"] == 1
        for marker, subject, values in records
    )
