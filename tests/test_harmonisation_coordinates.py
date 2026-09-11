"""Scientific regression tests for chromosome normalization and retention."""

from types import SimpleNamespace

import polars as pl

from postgwas.modules.harmonisation.concordance.analysis import (
    INPUT_SOURCE_ROW_COLUMN,
    PARTITION_COLUMN,
    input_chromosome_expression,
    reference_chromosome_expression,
    vcf_chromosome_expression,
)
from postgwas.modules.harmonisation.concordance.service import (
    _stage_indexed_csv_batches,
)
from postgwas.modules.harmonisation.coordinates import (
    harmonise_coordinates_and_alleles,
)
from postgwas.modules.harmonisation.policies import (
    default_policies,
)
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    SOURCE_INPUT_ROW_COLUMN,
)
from postgwas.modules.harmonisation.summary_statistics_io import (
    read_summary_statistics,
)


def test_packaged_chromosome_scope_is_autosomes_and_x_only():
    policies = default_policies()
    expected = [str(chromosome) for chromosome in range(1, 23)] + ["X"]

    assert policies.get("chromosome.allowed") == expected
    assert policies.get("chromosome.allowed_after_split") == expected
    assert policies.get("chromosome.strip_leading_zero") is True
    assert policies.get("chromosome.rename_map") == {
        "23": "X",
        "24": "Y",
        "25": "X",
        "26": "MT",
        "XY": "X",
        "PAR1": "X",
        "PAR2": "X",
        "M": "MT",
    }


def test_position_zero_is_rejected_as_invalid_with_original_row_provenance(tmp_path):
    source = pl.DataFrame({
        "variant": ["zero", "one"],
        "CHR": ["1", "1"],
        "POS": [0, 1],
        "EA": ["A", "A"],
        "OA": ["G", "G"],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    rejects = RejectCollector(
        source_snapshot=source,
        logger=None,
        out_path=str(tmp_path / "rejected.tsv"),
        delimiter="\t",
        compress=False,
    )

    result, _mapping, qc = harmonise_coordinates_and_alleles(
        chromosome="All_Chrs",
        df=source,
        sample_column_dict={
            "chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA",
        },
        policies=default_policies(),
        rejects=rejects,
        return_qc_info=True,
    )

    assert result["variant"].to_list() == ["one"]
    assert qc["position_min_value"] == 1
    assert qc["dropped_below_min_pos"] == 1
    assert qc["removed_by_reason"]["invalid_position"] == 1
    rejected = rejects.frame()
    assert rejected["variant"].to_list() == ["zero"]
    assert rejected["POS"].to_list() == ["0"]
    assert rejected["reject_reason"].to_list() == ["invalid_position"]
    assert rejected["reject_detail"].to_list() == [
        "position below position.min_value (1)"
    ]


def test_plink_par_is_retained_on_x_and_unsupported_chromosomes_are_counted():
    labels = [
        "1", "23", "25", "XY", "PAR1", "PAR2",
        "24", "26", "M", "MT", "Y", "GL000225.1", None,
    ]
    frame = pl.DataFrame({
        "CHR": labels,
        "POS": list(range(100, 100 + len(labels))),
        "EA": ["a"] * len(labels),
        "OA": ["g"] * len(labels),
    })

    result, _, qc = harmonise_coordinates_and_alleles(
        chromosome="All_Chrs",
        df=frame,
        sample_column_dict={
            "chr_col": "CHR",
            "pos_col": "POS",
            "ea_col": "EA",
            "oa_col": "OA",
        },
        policies=default_policies(),
        return_qc_info=True,
    )

    assert result["CHR"].to_list() == ["1", "X", "X", "X", "X", "X"]
    assert result["EA"].to_list() == ["A"] * 6
    assert result["OA"].to_list() == ["G"] * 6
    assert qc["par_variants_mapped_to_x"] == 4
    assert qc["dropped_null_chr"] == 1
    assert qc["dropped_null_pos"] == 0
    assert qc["dropped_y"] == 2
    assert qc["dropped_mt"] == 3
    assert qc["dropped_other_chromosomes"] == {"GL000225.1": 1}
    assert qc["dropped_unsupported_chromosomes"] == 6
    assert qc["removed_by_reason"]["invalid_chromosome"] == 1
    assert qc["removed_by_reason"]["invalid_position"] == 0
    assert qc["removed_by_reason"]["unsupported_chromosome"] == 6
    assert qc["initial_variants"] == 13
    assert qc["final_variants"] == 6
    assert qc["removed_total"] == 7


def test_input_accounting_separates_unsupported_chromosomes(tmp_path):
    source = tmp_path / "study.tsv"
    source.write_text(
        "CHR\tPOS\tEA\tOA\tSNP\tBETA\tP\n"
        "1\t100\tA\tG\tgood\t0.1\t0.1\n"
        "Y\t101\tA\tG\ty\t0.1\t0.1\n"
        "MT\t102\tA\tG\tmt\t0.1\t0.1\n"
        "2\tbad\tA\tG\tbad_position\t0.1\t0.1\n"
        "3\t103\tN\tG\tbad_allele\t0.1\t0.1\n",
        encoding="utf-8",
    )

    result = read_summary_statistics(
        sumstat_file=str(source),
        output_dir=str(tmp_path),
        sample_column_dict={
            "gwas_outputname": "study",
            "chr_col": "CHR",
            "pos_col": "POS",
            "snp_id_col": "SNP",
            "ea_col": "EA",
            "oa_col": "OA",
            "beta_or_col": "BETA",
            "pval_col": "P",
        },
        output_layout={
            "rejected_directory": "rejected",
            "input_reject": "rejected/{dataset_id}_input.tsv",
            "duplicates": "{dataset_id}_duplicates.tsv",
        },
        table_delimiter="\t",
        policies=default_policies(),
        input_line_count=5,
    )

    retained = result[0]
    invalid_coordinates = result[5]
    unsupported_chromosomes = result[6]
    invalid_alleles = result[7]
    duplicates = result[8]
    missing_mandatory = result[9]
    rejected = pl.read_csv(
        tmp_path / "rejected" / "study_input.tsv.gz",
        separator="\t",
        infer_schema_length=0,
    )

    assert retained["SNP"].to_list() == ["good"]
    assert (
        invalid_coordinates,
        unsupported_chromosomes,
        invalid_alleles,
        duplicates,
        missing_mandatory,
    ) == (1, 2, 1, 0, 0)
    assert 5 == retained.height + sum((1, 2, 1, 0, 0))
    assert rejected.group_by("reject_reason").len().sort("reject_reason").rows() == [
        ("invalid_position", 1),
        ("non_standard_allele", 1),
        ("unsupported_chromosome", 2),
    ]


def test_drop_mt_false_does_not_add_mt_to_an_autosomes_x_scope():
    policies = default_policies().with_overrides({
        "chromosome": {"drop_mt": False},
    })
    result, _, qc = harmonise_coordinates_and_alleles(
        chromosome="All_Chrs",
        df=pl.DataFrame({
            "CHR": ["X", "MT"],
            "POS": [100, 200],
            "EA": ["A", "C"],
            "OA": ["G", "T"],
        }),
        sample_column_dict={
            "chr_col": "CHR",
            "pos_col": "POS",
            "ea_col": "EA",
            "oa_col": "OA",
        },
        policies=policies,
        return_qc_info=True,
    )

    assert result["CHR"].to_list() == ["X"]
    assert qc["allowed_chromosomes"] == [
        str(chromosome) for chromosome in range(1, 23)
    ] + ["X"]
    assert qc["dropped_mt"] == 1


def test_concordance_uses_the_same_configured_chromosome_aliases():
    policies = default_policies()
    labels = [
        "23", "25", "XY", "PAR1", "PAR2", "24", "26", "M",
        "chr01", "001", "chr023", "023.0",
    ]
    expected = [
        "X", "X", "X", "X", "X", "Y", "MT", "MT",
        "1", "1", "X", "X",
    ]
    frame = pl.DataFrame({
        "CHR": labels,
        "POS": list(range(1, len(labels) + 1)),
        "REF_CHR": labels,
        "CHROM": labels,
    })
    row = SimpleNamespace(
        chromosome_column="CHR",
        position_column="POS",
        chromosome_position_column=None,
    )

    normalized = frame.select(
        input_chromosome_expression(row, policies).alias("input"),
        reference_chromosome_expression("REF_CHR", policies).alias("reference"),
        vcf_chromosome_expression(policies).alias("vcf"),
    )

    assert normalized["input"].to_list() == expected
    assert normalized["reference"].to_list() == expected
    assert normalized["vcf"].to_list() == expected


def test_leading_zero_policy_is_shared_by_input_reference_and_vcf():
    policies = default_policies().with_overrides({
        "chromosome": {"strip_leading_zero": False},
    })
    frame = pl.DataFrame({
        "CHR": ["01"], "POS": [1], "REF_CHR": ["01"], "CHROM": ["01"],
    })
    row = SimpleNamespace(
        chromosome_column="CHR",
        position_column="POS",
        chromosome_position_column=None,
    )

    normalized = frame.select(
        input_chromosome_expression(row, policies).alias("input"),
        reference_chromosome_expression("REF_CHR", policies).alias("reference"),
        vcf_chromosome_expression(policies).alias("vcf"),
    )

    assert normalized.row(0) == ("01", "01", "01")


def test_concordance_batch_staging_preserves_source_rows_and_columns(tmp_path):
    source = tmp_path / "input.tsv"
    destination = tmp_path / "input.parquet"
    frame = pl.DataFrame({
        "CHR": [1, 2, 23, 25] * 5,
        "POS": list(range(100, 120)),
    })
    frame.write_csv(source, separator="\t")
    policies = default_policies()
    row = SimpleNamespace(
        chromosome_column="CHR",
        position_column="POS",
        chromosome_position_column=None,
    )

    rows = _stage_indexed_csv_batches(
        source,
        destination,
        separator="\t",
        columns=["CHR", "POS"],
        partition_expression=input_chromosome_expression(row, policies),
        null_values=list(policies.get("input.null_values")),
        schema_inference_rows=int(policies.get("input.schema_inference_rows")),
        batch_rows=5,
        compression="zstd",
        row_index_name=INPUT_SOURCE_ROW_COLUMN,
        comment_prefix=None,
    )
    staged = pl.read_parquet(destination)

    assert rows == frame.height
    assert staged[INPUT_SOURCE_ROW_COLUMN].to_list() == list(range(1, 21))
    assert staged["CHR"].to_list() == [str(value) for value in frame["CHR"]]
    assert staged["POS"].to_list() == [str(value) for value in frame["POS"]]
    assert staged[PARTITION_COLUMN].to_list() == [
        "1", "2", "X", "X",
    ] * 5
