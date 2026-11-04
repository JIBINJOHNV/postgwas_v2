"""Complete-site AF attachment retains the existing adapter input contract."""

import json
from pathlib import Path
import subprocess
import sys

import polars as pl
import pysam
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation import gwas2vcf_export as export_module
from postgwas.modules.harmonisation.gwas2vcf_export import export_gwas2vcf_input
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.shared.allele_join import allele_oriented_left_join


def _complete_site_export_fixture(ea="A", oa="G", extra=.1):
    mapping = {
        "chr_col": "CHR", "pos_col": "POS", "snp_id_col": "SNP",
        "ea_col": "EA", "oa_col": "OA", "eaf_col": "EAF",
        "beta_col": "BETA", "se_col": "SE", "imp_z_col": "Z",
        "pval_col": "P", "ncontrol_col": "N", "neff_col": "NEFF",
    }
    source = pl.DataFrame({
        "CHR": ["1", "1"], "POS": [10, 20], "SNP": ["variant1", "variant2"],
        "EA": [ea, ea], "OA": [oa, oa],
        "BETA": [.2, .3], "SE": [.1, .1], "Z": [2., 3.],
        "P": [.0455, .0027], "N": [1000, 1000], "NEFF": [1000, 1000],
        "RAW_MAF": [.3, .3], "strand_action": ["forward", "forward"],
    })
    reference = pl.DataFrame({
        "CHROM": ["1"] * 4, "POS": [10, 10, 20, 20],
        "REF": ["A"] * 4, "ALT": ["G", "C", "G", "C"],
        "AF": [.3, extra, .3, extra],
    })
    joined, _, _ = allele_oriented_left_join(
        source, reference,
        study_columns={"chr": "CHR", "pos": "POS", "ea": "EA", "oa": "OA"},
        reference_columns={"chr": "CHROM", "pos": "POS", "ea": "ALT", "oa": "REF"},
        policies=default_policies(), value_column="AF", output_column="EAF",
        duplicate_exact_action="keep_one", duplicate_non_identical_action="discard_all",
        swapped_value="one_minus", prefer_non_null_value=True,
        screen_frequency=False, require_complete_frequency_site=True,
    )
    assert joined.select(source.columns).equals(source)
    return joined, mapping


def _export(frame, mapping, output):
    config = load_configuration().modules.harmonisation
    return export_gwas2vcf_input(
        frame, mapping, str(output), "study", "1", genome_build="GRCh37",
        layout=config.output_layout.model_dump(),
        input_config=config.gwas2vcf_input.model_dump(),
    )


def test_complete_reference_frequency_exports_without_extra_fasta_or_liftover_input(tmp_path):
    frame, mapping = _complete_site_export_fixture()
    original = frame.clone()
    output = tmp_path / "export"

    summary = _export(frame, mapping, output)

    exported = pl.read_csv(output / "study_chr1_vcf_input.tsv", separator="\t")
    positions = json.loads((output / "study_chr1.dict").read_text(encoding="utf-8"))
    assert frame.equals(original)
    assert exported["EAF"].to_list() == pytest.approx([.6, .6])
    assert exported["BETA"].to_list() == [.2, .3]
    assert exported["SE"].to_list() == [.1, .1]
    assert exported["Z"].to_list() == [2., 3.]
    assert exported["P"].to_list() == [.0455, .0027]
    assert "RAW_MAF" not in exported.columns
    assert "eaf_swap_guard_col" not in positions
    assert "eaf_swap_guard_col" not in summary["key"].to_list()


@pytest.mark.parametrize("ea, oa, fasta_base, extra, expected_ref, expected_alt, expected_af, sign", [
    ("G", "A", "A", .1, "A", "G", .3, 1),
    ("A", "G", "G", .1, "G", "A", .6, 1),
    ("A", "G", "A", 0., "A", "G", .3, -1),
])
def test_native_adapter_keeps_supported_reference_frequencies_and_allele_aware_statistics(
    tmp_path, ea, oa, fasta_base, extra, expected_ref, expected_alt, expected_af, sign,
):
    frame, mapping = _complete_site_export_fixture(ea, oa, extra)
    source = frame.clone()
    output = tmp_path / "export"
    _export(frame, mapping, output)
    fasta = tmp_path / "source.fa"
    fasta.write_text(">1\n" + fasta_base * 100 + "\n", encoding="utf-8")
    pysam.faidx(str(fasta))
    script = Path(export_module.__file__).parent / "adapters" / "gwas2vcf" / "main.py"
    destination = output / "native.vcf"
    completed = subprocess.run(
        [sys.executable, str(script), "--data", str(output / "study_chr1_vcf_input.tsv"),
         "--ref", str(fasta), "--out", str(destination), "--id", "STUDY",
         "--json", str(output / "study_chr1.dict")],
        capture_output=True, text=True, check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    destination = Path(str(destination) + ".gz")
    assert Path(str(destination) + ".tbi").is_file()
    assert frame.equals(source)
    with pysam.VariantFile(str(destination)) as vcf:
        records = list(vcf.fetch("1"))
    assert len(records) == 2
    for index, record in enumerate(records):
        assert record.pos == [10, 20][index]
        assert record.ref == expected_ref
        assert record.alts == (expected_alt,)
        sample = record.samples["STUDY"]
        assert "GT" not in sample
        assert sample["AF"] == pytest.approx((expected_af,))
        assert record.info["AF"] == pytest.approx((expected_af,))
        assert sample["ES"] == pytest.approx(([.2, .3][index] * sign,))
        assert sample["SE"] == pytest.approx((.1,))
        assert sample["EZ"] == pytest.approx(([2., 3.][index] * sign,))
        # Native text rounds LP to six significant digits.
        assert sample["LP"] == pytest.approx(([1.34199, 2.56864][index],))
        assert sample["SS"] == (1000,)
        assert sample["NEF"] == (1000,)
