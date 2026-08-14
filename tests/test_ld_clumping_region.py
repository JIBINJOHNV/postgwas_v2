import polars as pl

from postgwas.modules.ld_clumping.ld_prune_region import ld_clump_by_regions


def test_region_clumping_queries_only_the_selected_population_block(monkeypatch, tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf")
    observed = {}

    def fake_extract(vcf_path, table_path, dataset_id, columns, *args, **kwargs):
        observed["columns"] = columns
        pl.DataFrame({
            "CHR": ["1"],
            "BP": [15],
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
    assert "AFR_LDblock" not in observed["columns"]
    assert "EAS_LDblock" not in observed["columns"]
    significant = pl.read_csv(result["ldpruned_sig_file"], separator="\t")
    assert significant.select("SNP", "START", "END").row(0) == (
        "1_15_A_G", 10, 20,
    )
    assert result["significant_blocks"] == 1
