"""Missing-token contracts for summary-statistic fields and identifiers."""

import polars as pl

from postgwas.core.values import missing_tokens, optional_text
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.summary_statistics_io import (
    load_summary_statistics_table,
)
from postgwas.modules.harmonisation.variant_identifiers import (
    harmonise_variant_identifiers,
)


MAPPING = {
    "snp_id_col": "SNP",
    "chr_col": "CHR",
    "pos_col": "POS",
    "ea_col": "EA",
    "oa_col": "OA",
}


def test_shared_missing_tokens_support_normalized_and_exact_contexts():
    assert "-" in missing_tokens()
    assert "NULL" in missing_tokens([" null "])
    assert missing_tokens(
        [" null "], include_standard=False, exact=True,
    ) == (" null ",)
    assert optional_text(" null ") is None


def test_identifier_placeholders_are_rebuilt_without_erasing_deletion_allele(
    tmp_path,
):
    source = tmp_path / "study.tsv"
    source.write_text(
        "SNP\tCHR\tPOS\tEA\tOA\n"
        "-\t1\t100\tA\tG\n"
        "null\t1\t101\tC\tT\n"
        "NoNe\t1\t102\tG\tA\n"
        "N/A\t1\t103\tT\tC\n"
        "not_applicable\t1\t104\tA\t-\n",
        encoding="utf-8",
    )
    policies = default_policies()

    frame, mapping = load_summary_statistics_table(
        str(source), dict(MAPPING), policies=policies,
    )

    assert frame["SNP"].to_list() == ["-", "null", "NoNe", "N/A", "not_applicable"]
    assert frame["OA"].to_list()[-1] == "-"

    harmonised, mapping = harmonise_variant_identifiers(
        "1", frame, mapping, policies=policies,
    )

    assert mapping["snp_id_col"] == "SNP"
    assert harmonised["SNP"].to_list() == [
        "1_100_A_G",
        "1_101_C_T",
        "1_102_G_A",
        "1_103_T_C",
        "1_104_A_-",
    ]
    assert harmonised["OA"].to_list()[-1] == "-"


def test_configured_reader_token_is_also_a_missing_identifier(tmp_path):
    source = tmp_path / "study.tsv"
    source.write_text(
        "SNP\tCHR\tPOS\tEA\tOA\nMISSING_ID\t2\t200\tA\tC\n",
        encoding="utf-8",
    )
    configured = list(default_policies().get("input.null_values")) + [
        "MISSING_ID",
    ]
    policies = default_policies().with_overrides({
        "input.null_values": configured,
    })

    frame, mapping = load_summary_statistics_table(
        str(source), dict(MAPPING), policies=policies,
    )
    assert frame["SNP"].to_list() == [None]

    harmonised, _mapping = harmonise_variant_identifiers(
        "2", frame, mapping, policies=policies,
    )
    assert harmonised["SNP"].to_list() == ["2_200_A_C"]


def test_identifier_matching_remains_vectorized():
    frame = pl.DataFrame({
        "SNP": ["rs1", "-"],
        "CHR": ["1", "1"],
        "POS": [100, 101],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
    })

    harmonised, _mapping = harmonise_variant_identifiers(
        "1", frame, dict(MAPPING), policies=default_policies(),
    )
    assert harmonised["SNP"].to_list() == ["rs1", "1_101_C_T"]
