import gzip
from pathlib import Path

import polars as pl

from postgwas.modules.ld_clumping.ld_prune_region import ld_clump_by_regions


def test_region_clumping_queries_only_the_selected_population_block(monkeypatch, tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf")
    raw_output = tmp_path / "STUDY_vcf.tsv.gz"
    observed = {}

    def fake_run(command, **kwargs):
        observed["command"] = command
        with gzip.open(raw_output, "wt", encoding="utf-8") as handle:
            handle.write(
                "SNP\tCHR\tBP\tREF\tALT\tBETA\tSE\tNEF\tAF\tAFR\tEAS\t"
                "EUR\tSAS\tLP\tNC\tNCO\tEUR_LDblock\tBCSQ\n"
            )
            handle.write(
                "1:15:A:G\t1\t15\tA\tG\t0.2\t0.1\t1000\t0.2\t0.1\t0.1\t"
                "0.2\t0.1\t9\t500\t500\tEUR-1_10_20\t.\n"
            )

    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.ld_prune_region.subprocess.run",
        fake_run,
    )

    result = ld_clump_by_regions(
        str(vcf), str(tmp_path), "STUDY", population="EUR", bcftools="bcftools",
    )

    command = observed["command"]
    assert "set -o pipefail" in command
    assert "%INFO/EUR_LDblock" in command
    assert "%INFO/AFR_LDblock" not in command
    assert "%INFO/EAS_LDblock" not in command
    significant = pl.read_csv(result["ldpruned_sig_file"], separator="\t")
    assert significant.select("SNP", "START", "END").row(0) == (
        "1_15_A_G", 10, 20,
    )
