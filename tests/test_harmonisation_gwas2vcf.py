"""Tests for the bundled GWAS-to-VCF process boundary."""

import json
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pytest
import polars as pl
import pysam

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.cli import _engine_defaults
from postgwas.modules.harmonisation.effect_validation import (
    final_completeness_check,
)
from postgwas.modules.harmonisation.gwas2vcf_runner import run_gwas2vcf
from postgwas.modules.harmonisation.gwas2vcf_export import (
    export_gwas2vcf_input,
)
from postgwas.modules.harmonisation.imputation_quality import (
    harmonise_imputation_quality,
)
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.adapters.gwas2vcf.gwas import normalize_alleles
from postgwas.modules.harmonisation.adapters.gwas2vcf.main import main as adapter_main
from postgwas.modules.harmonisation.vcf_processing import (
    MissingExecutableError,
    annotate_and_liftover_vcf,
    require_binaries,
)


def _inputs(root: Path):
    for name in ("reference.fa", "dbsnp.vcf.gz", "adapter.py"):
        (root / name).write_text("fixture\n", encoding="utf-8")
    return root / "reference.fa", root / "dbsnp.vcf.gz", root / "adapter.py"


def _harmonisation_defaults():
    return _engine_defaults(load_configuration())


def test_adapter_uses_current_python_without_passing_dbsnp():
    with TemporaryDirectory() as directory:
        root = Path(directory)
        fasta, dbsnp, adapter = _inputs(root)
        defaults = _harmonisation_defaults()
        with patch("subprocess.run", return_value=subprocess.CompletedProcess([], 0)) as run:
            command, exit_code = run_gwas2vcf(
                "study", "1", "GRCh37", str(root), str(fasta),
                main_script_path=str(adapter),
                output_layout=defaults["output_layout"],
                python_bin=sys.executable,
            )

        argv = run.call_args.args[0]
        assert exit_code == 0
        assert argv[0] == str(Path(sys.executable).resolve())
        assert "--dbsnp" not in argv
        assert str(dbsnp) not in command
        assert argv[argv.index("--out") + 1].endswith("study_chr1_GRCh37.vcf")


def test_adapter_failure_raises_with_captured_diagnostic():
    with TemporaryDirectory() as directory:
        root = Path(directory)
        fasta, _dbsnp, adapter = _inputs(root)
        defaults = _harmonisation_defaults()

        def fail(_argv, **kwargs):
            kwargs["stdout"].write("ModuleNotFoundError: missing dependency\n")
            return subprocess.CompletedProcess([], 7)

        with patch("subprocess.run", side_effect=fail):
            with pytest.raises(RuntimeError, match="missing dependency"):
                run_gwas2vcf(
                    "study", "2", "GRCh37", str(root), str(fasta),
                    main_script_path=str(adapter),
                    output_layout=defaults["output_layout"],
                    python_bin=sys.executable,
                )

        transcript = (
            root / "logs" / "adapters" / "gwas2vcf" / "study_chr2.log"
        ).read_text(
            encoding="utf-8"
        )
        assert transcript.startswith("tool=gwas2vcf\ncommand=")
        assert "missing dependency" in transcript
        assert "command=" in transcript
        assert "exit_code=7" in transcript


def test_export_preserves_rows_and_summarises_only_adapter_mappings():
    columns = {
        "chr_col": "CHR", "pos_col": "BP", "snp_id_col": "SNP",
        "ea_col": "A1", "oa_col": "A2", "eaf_col": "EAF",
        "beta_col": "BETA", "se_col": "SE", "imp_z_col": "Z",
        "pval_col": "P", "ncontrol_col": "N",
        # Source/provenance mappings can point to final exported columns after
        # harmonisation, but they are not separate GWAS-to-VCF inputs.
        "eafcolumn": "EAF", "beta_or_col": "BETA",
    }
    frame = pl.DataFrame(
        {
            "CHR": ["1", "1"], "BP": [1, 2], "SNP": ["a", "b"],
            "A1": ["A", "C"], "A2": ["G", "T"], "EAF": [0.2, 0.3],
            "BETA": [0.1, 0.0], "SE": [0.2, 0.0], "Z": [0.5, 0.0],
            "P": [0.6, 1.0], "N": [100, 100],
            "strand_action": ["forward", "forward_swapped"],
        }
    )
    with TemporaryDirectory() as directory:
        defaults = _harmonisation_defaults()
        assert defaults["gwas2vcf_input"]["audit_columns"] == ["strand_action"]
        with pytest.raises(ValueError, match="configured audit columns are absent"):
            export_gwas2vcf_input(
                frame.drop("strand_action"), columns, directory, "study", "1",
                genome_build="GRCh37", layout=defaults["output_layout"],
                input_config=defaults["gwas2vcf_input"],
            )
        summary = export_gwas2vcf_input(
            frame, columns, directory, "study", "1", genome_build="GRCh37",
            layout=defaults["output_layout"],
            input_config=defaults["gwas2vcf_input"],
        )
        exported = pl.read_csv(Path(directory) / "study_chr1_vcf_input.tsv", separator="\t")
        mapping = json.loads(
            (Path(directory) / "study_chr1.dict").read_text(encoding="utf-8")
        )
    assert exported.height == frame.height
    assert exported.columns[-1] == "strand_action"
    assert exported["strand_action"].to_list() == ["forward", "forward_swapped"]
    assert "strand_action" not in mapping
    assert summary["num_rows"][0] == frame.height
    assert summary["key"].to_list() == list(
        defaults["gwas2vcf_input"]["required_column_keys"]
    )
    assert "eafcolumn" not in summary["key"].to_list()
    assert "beta_or_col" not in summary["key"].to_list()
    assert summary.filter(pl.col("key") == "eaf_col").height == 1
    assert "tsv_path" not in summary.columns
    assert "dict_path" not in summary.columns


def test_fixed_info_survives_completeness_check_and_reaches_adapter_export(tmp_path):
    columns = {
        "chr_col": "CHR", "pos_col": "BP", "snp_id_col": "SNP",
        "ea_col": "A1", "oa_col": "A2", "eaf_col": "EAF",
        "beta_col": "BETA", "se_col": "SE", "imp_z_col": "Z",
        "pval_col": "P", "ncontrol_col": "NEFF", "imp_info_col": None,
        "fixed_info": 0.9, "fixed_info_column": "__postgwas_fixed_info",
        "info_source": "fixed_cli",
    }
    frame = pl.DataFrame({
        "CHR": ["1"], "BP": [10], "SNP": ["rs1"],
        "A1": ["A"], "A2": ["G"], "EAF": [0.2],
        "BETA": [0.1], "SE": [0.05], "Z": [2.0],
        "P": [0.0455], "NEFF": [1000], "__temporary_qc_flag": [True],
        "strand_action": ["forward"],
    })
    policies = default_policies()

    frame, _info_qc, columns = harmonise_imputation_quality(
        "1", frame, columns, policies=policies,
    )
    frame, completeness_qc, columns = final_completeness_check(
        "1", frame, columns, policies=policies,
    )
    defaults = _harmonisation_defaults()
    export_gwas2vcf_input(
        frame, columns, str(tmp_path), "study", "1", genome_build="GRCh37",
        layout=defaults["output_layout"],
        input_config=defaults["gwas2vcf_input"],
    )

    exported = pl.read_csv(
        tmp_path / "study_chr1_vcf_input.tsv", separator="\t",
    )
    assert "__temporary_qc_flag" not in frame.columns
    assert frame["__postgwas_fixed_info"].to_list() == [0.9]
    assert exported["__postgwas_fixed_info"].to_list() == [0.9]
    assert completeness_qc["internal_columns_dropped"] == [
        "__temporary_qc_flag",
    ]
    assert completeness_qc["internal_columns_retained_for_export"] == [
        "__postgwas_fixed_info",
    ]


@pytest.mark.parametrize(
    ("reference", "start", "stop", "alleles", "expected"),
    [
        ("CCATGG", 2, 4, ("AT", "AG"), (3, 4, ("T", "G"))),
        ("CAAAG", 2, 3, ("A", "AA"), (0, 0, ("C", "CA"))),
        ("GACCT", 1, 3, ("AC", "TC"), (1, 2, ("A", "T"))),
    ],
)
def test_pure_python_allele_normalization(reference, start, stop, alleles, expected):
    assert normalize_alleles(reference, start, stop, alleles) == expected


def test_adapter_cli_needs_no_dbsnp_and_preserves_study_id_in_format(tmp_path):
    reference = tmp_path / "reference.fa"
    reference.write_text(">1\nAAAA\n", encoding="utf-8")
    pysam.faidx(str(reference))
    input_table = tmp_path / "study.tsv"
    input_table.write_text(
        "CHR\tPOS\tEA\tOA\tBETA\tSE\tP\tID\tstrand_action\n"
        "1\t2\tG\tA\t0.1\t0.05\t0.01\tstudy:variant=1\tforward\n",
        encoding="utf-8",
    )
    output = tmp_path / "study.vcf"
    mapping = tmp_path / "study.json"
    mapping.write_text(
        json.dumps({
            "chr_col": 0, "pos_col": 1, "ea_col": 2, "oa_col": 3,
            "beta_col": 4, "se_col": 5, "pval_col": 6, "snp_col": 7,
            "delimiter": "\t", "header": True, "build": "GRCh37",
        }),
        encoding="utf-8",
    )

    with patch.object(sys, "argv", [
        "gwas2vcf", "--data", str(input_table), "--ref", str(reference),
        "--out", str(output), "--id", "STUDY", "--json", str(mapping),
    ]):
        adapter_main()

    with pysam.VariantFile(str(output) + ".gz") as vcf:
        record = next(iter(vcf))
        assert tuple(vcf.header.samples) == ("STUDY",)
        assert record.id == "study_variant_1"
        assert record.samples["STUDY"]["ID"] == "study_variant_1"


def test_runtime_preflight_names_missing_bcftools_plugin():
    plugin_list = subprocess.CompletedProcess([], 0, stdout="fill-tags\nfixref\n")
    with patch("shutil.which", side_effect=lambda name: "/tools/" + name), patch(
        "subprocess.run", return_value=plugin_list
    ):
        with pytest.raises(MissingExecutableError, match="liftover"):
            require_binaries(
                {"bcftools": "bcftools"}, plugins=("liftover",),
            )


def test_missing_id_format_reaches_bcftools_with_escaped_separators(tmp_path):
    defaults = _harmonisation_defaults()
    layout = defaults["output_layout"]
    vcf_config = defaults["vcf_processing"]
    (tmp_path / "study_chr1_GRCh37.vcf.gz").write_bytes(b"fixture")

    with patch(
        "postgwas.modules.harmonisation.vcf_processing._run_bcftools_step"
    ) as run_step, patch(
        "postgwas.modules.harmonisation.vcf_processing.get_vcf_variant_count",
        return_value=1,
    ):
        annotate_and_liftover_vcf(
            str(tmp_path), "study", "1", "af.vcf.gz", "dbsnp.vcf.gz",
            "GRCh37.fa", "GRCh38.fa", "genes.gff3.gz", "chain.gz",
            "GRCh37", 1, layout,
            {"bash": "/bin/bash", "bcftools": "bcftools", "tabix": "tabix"},
            vcf_config,
        )

    id_call = next(
        call for call in run_step.call_args_list if call.args[2] == "STEP1_ID"
    )
    script = id_call.args[0][2]
    assert "--pair-logic exact -Ou" in script
    assert "--set-id '+%CHROM\\_%POS\\_%REF\\_%ALT'" in script
    assert script.count("bcftools annotate") == 2
    assert "bcftools view" not in script
