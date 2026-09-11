"""Tests for the bundled GWAS-to-VCF process boundary."""

import ast
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
from postgwas.modules.harmonisation.p_values import harmonise_p_values
from postgwas.modules.harmonisation.adapters.gwas2vcf import gwas as gwas_module
from postgwas.modules.harmonisation.adapters.gwas2vcf import main as adapter_main_module
from postgwas.modules.harmonisation.adapters.gwas2vcf import vcf as vcf_module
from postgwas.modules.harmonisation.adapters.gwas2vcf.gwas import (
    Gwas,
    ReferenceAlleleMismatchError,
    ReferenceSequenceFetchError,
    normalize_alleles,
)
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
        raw_vcf = root / "study_chr1_GRCh37.vcf.gz"

        def succeed(argv, **_kwargs):
            if argv[1:3] == ["index", "-n"]:
                return subprocess.CompletedProcess(argv, 0, stdout="1\n", stderr="")
            raw_vcf.write_text("raw-vcf-fixture\n", encoding="utf-8")
            Path(str(raw_vcf) + ".tbi").write_text("index\n", encoding="utf-8")
            return subprocess.CompletedProcess(argv, 0)

        with patch("subprocess.run", side_effect=succeed) as run:
            command, exit_code = run_gwas2vcf(
                "study", "1", "GRCh37", str(root), str(fasta),
                main_script_path=str(adapter),
                output_layout=defaults["output_layout"],
                python_bin=sys.executable,
                bcftools_bin=sys.executable,
                expected_variants=1,
            )

        argv = run.call_args_list[0].args[0]
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
                    bcftools_bin=sys.executable,
                    expected_variants=1,
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


def test_adapter_success_without_raw_vcf_is_rejected():
    with TemporaryDirectory() as directory:
        root = Path(directory)
        fasta, _dbsnp, adapter = _inputs(root)
        defaults = _harmonisation_defaults()
        with patch(
            "subprocess.run", return_value=subprocess.CompletedProcess([], 0)
        ):
            with pytest.raises(
                RuntimeError, match="expected output is missing or empty"
            ):
                run_gwas2vcf(
                    "study", "3", "GRCh37", str(root), str(fasta),
                    main_script_path=str(adapter),
                    output_layout=defaults["output_layout"],
                    python_bin=sys.executable,
                    bcftools_bin=sys.executable,
                    expected_variants=1,
                )


def test_adapter_record_count_must_equal_exported_rows():
    with TemporaryDirectory() as directory:
        root = Path(directory)
        fasta, _dbsnp, adapter = _inputs(root)
        defaults = _harmonisation_defaults()
        raw_vcf = root / "study_chr4_GRCh37.vcf.gz"

        def partial_output(argv, **_kwargs):
            if argv[1:3] == ["index", "-n"]:
                return subprocess.CompletedProcess(argv, 0, stdout="1\n", stderr="")
            raw_vcf.write_text("raw-vcf-fixture\n", encoding="utf-8")
            Path(str(raw_vcf) + ".tbi").write_text("index\n", encoding="utf-8")
            return subprocess.CompletedProcess(argv, 0)

        with patch("subprocess.run", side_effect=partial_output):
            with pytest.raises(
                RuntimeError,
                match=r"2 harmonised row\(s\).*contains 1 record\(s\)",
            ):
                run_gwas2vcf(
                    "study", "4", "GRCh37", str(root), str(fasta),
                    main_script_path=str(adapter),
                    output_layout=defaults["output_layout"],
                    python_bin=sys.executable,
                    bcftools_bin=sys.executable,
                    expected_variants=2,
                )


def test_adapter_json_validation_exits_nonzero(monkeypatch, tmp_path):
    malformed_json = tmp_path / "malformed.json"
    malformed_json.write_text("{not-json", encoding="utf-8")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gwas2vcf",
            "--ref", str(tmp_path / "reference.fa"),
            "--json", str(malformed_json),
        ],
    )
    with pytest.raises(SystemExit) as exc_info:
        adapter_main()
    assert exc_info.value.code == 1


def _gwas_record(chromosome: str, position: int, reference: str) -> Gwas:
    return Gwas(
        chromosome, position, reference, "G", 0.1, 0.2, 1.0, 1000,
        0.2, "variant", None, 1000, 1000, 0.9, 0.5,
    )


def test_fasta_access_error_is_not_treated_as_allele_mismatch(tmp_path):
    fasta_path = tmp_path / "reference.fa"
    fasta_path.write_text(">1\nACGTACGT\n", encoding="utf-8")
    pysam.faidx(str(fasta_path))

    with pysam.FastaFile(str(fasta_path)) as fasta:
        with pytest.raises(ReferenceAlleleMismatchError):
            _gwas_record("1", 1, "T").check_reference_allele(fasta)
        with pytest.raises(ReferenceSequenceFetchError, match="chromosome naming"):
            _gwas_record("chr1", 1, "A").check_reference_allele(fasta)
        with pytest.raises(ReferenceSequenceFetchError, match="outside"):
            _gwas_record("1", 20, "A").check_reference_allele(fasta)


def test_adapter_runtime_validation_does_not_use_assert_statements():
    for module in (gwas_module, vcf_module):
        source_path = Path(module.__file__)
        syntax_tree = ast.parse(source_path.read_text(encoding="utf-8"))
        assert not any(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))


def test_adapter_validation_never_uses_zero_status_exit():
    source_path = Path(adapter_main_module.__file__)
    syntax_tree = ast.parse(source_path.read_text(encoding="utf-8"))
    exit_calls = [
        node for node in ast.walk(syntax_tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "sys"
        and node.func.attr == "exit"
    ]
    assert exit_calls
    assert all(
        len(node.args) == 1
        and isinstance(node.args[0], ast.Constant)
        and node.args[0].value != 0
        for node in exit_calls
    )


def test_export_preserves_rows_and_summarises_only_adapter_mappings():
    columns = {
        "chr_col": "CHR", "pos_col": "BP", "snp_id_col": "SNP",
        "ea_col": "A1", "oa_col": "A2", "eaf_col": "EAF",
        "beta_col": "BETA", "se_col": "SE", "imp_z_col": "Z",
        "pval_col": "PVAL", "ncontrol_col": "N",
        # Source/provenance mappings can point to final exported columns after
        # harmonisation, but they are not separate GWAS-to-VCF inputs.
        "eafcolumn": "EAF", "beta_or_col": "BETA",
    }
    frame = pl.DataFrame(
        {
            "CHR": ["1", "1"], "BP": [1, 2], "SNP": ["a", "b"],
            "A1": ["A", "C"], "A2": ["G", "T"], "EAF": [0.2, 0.3],
            "BETA": [0.1, 0.0], "SE": [0.2, 0.0], "Z": [0.5, 0.0],
            "PVAL": [0.6, 1.0], "N": [100, 100],
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
        with pytest.raises(ValueError, match="neither reference-resolved"):
            export_gwas2vcf_input(
                frame.with_columns(pl.lit("disabled").alias("strand_action")),
                columns, directory, "study", "1", genome_build="GRCh37",
                layout=defaults["output_layout"],
                input_config=defaults["gwas2vcf_input"],
            )
        retained_frame = frame.with_columns(pl.Series(
            "strand_action",
            [
                "reference_unmatched_retained",
                "reference_unmatched_retained_reverse_complement",
            ],
        ))
        export_gwas2vcf_input(
            retained_frame, columns, directory, "study", "1",
            genome_build="GRCh37", layout=defaults["output_layout"],
            input_config=defaults["gwas2vcf_input"],
        )
        retained_export = pl.read_csv(
            Path(directory) / "study_chr1_vcf_input.tsv", separator="\t",
        )
        assert retained_export["strand_action"].to_list() == (
            retained_frame["strand_action"].to_list()
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


def test_extreme_raw_pvalue_survives_final_gate_and_adapter_lp_export(tmp_path):
    columns = {
        "chr_col": "CHR", "pos_col": "BP", "snp_id_col": "SNP",
        "ea_col": "A1", "oa_col": "A2", "eaf_col": "EAF",
        "beta_col": "BETA", "se_col": "SE", "imp_z_col": "Z",
        "pval_col": "P", "ncontrol_col": "N",
    }
    frame = pl.DataFrame({
        "CHR": ["1"], "BP": [2], "SNP": ["rs1"],
        "A1": ["G"], "A2": ["A"], "EAF": [0.2],
        "BETA": [0.1], "SE": [0.002334], "Z": [42.826406],
        "P": ["1e-400"], "N": [1000], "strand_action": ["forward"],
    })
    policies = default_policies()
    frame, _p_qc, columns = harmonise_p_values(
        "1", frame, columns, decision="raw", policies=policies,
    )
    frame, completeness_qc, columns = final_completeness_check(
        "1", frame, columns, policies=policies,
    )
    assert frame.height == 1
    assert completeness_qc["missing_by_field"]["pval"] == 0

    defaults = _harmonisation_defaults()
    export_gwas2vcf_input(
        frame, columns, str(tmp_path), "study", "1", genome_build="GRCh37",
        layout=defaults["output_layout"],
        input_config=defaults["gwas2vcf_input"],
    )
    input_table = tmp_path / "study_chr1_vcf_input.tsv"
    header, record = input_table.read_text(encoding="utf-8").splitlines()
    values = dict(zip(header.split("\t"), record.split("\t")))
    assert values[columns["pval_col"]].lower() == "1e-400"

    reference = tmp_path / "reference.fa"
    reference.write_text(">1\nAAAA\n", encoding="utf-8")
    pysam.faidx(str(reference))
    output = tmp_path / "study.vcf"
    mapping = tmp_path / "study_chr1.dict"
    with patch.object(sys, "argv", [
        "gwas2vcf", "--data", str(input_table), "--ref", str(reference),
        "--out", str(output), "--id", "STUDY", "--json", str(mapping),
    ]):
        adapter_main()

    with pysam.VariantFile(str(output) + ".gz") as vcf:
        vcf_record = next(iter(vcf))
        assert vcf_record.samples["STUDY"]["LP"] == pytest.approx((400.0,))


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
        "1\t2\tG\tA\t0.1\t0.05\t0.01\tstudy:variant=1\t"
        "reference_unmatched_retained\n",
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
        assert record.samples["STUDY"]["LP"] == pytest.approx((2.0,))


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

    header = (
        "##INFO=<ID=AF,Number=A,Type=Float>\n"
        "##INFO=<ID=AFR,Number=A,Type=Float>\n"
        "##INFO=<ID=EAS,Number=A,Type=Float>\n"
        "##INFO=<ID=EUR,Number=A,Type=Float>\n"
        "##INFO=<ID=SAS,Number=A,Type=Float>\n"
        "##FORMAT=<ID=AF,Number=A,Type=Float>\n"
        "##FORMAT=<ID=ES,Number=A,Type=Float>\n"
        "##FORMAT=<ID=EZ,Number=A,Type=Float>\n"
    )
    liftover_summary = (
        "Lines   total/swapped/reference added/rejected:\t1/0/0/0\n"
        "Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:"
        "\t1/0/0/0/0/0/0\n"
        "Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:"
        "\t1/0/0/0/0/0/0\n"
    )

    def completed_step(_command, _transcript, step, _chromosome, logger=None):
        return liftover_summary if step == "STEP4_LIFTOVER" else ""

    with patch(
        "postgwas.modules.harmonisation.vcf_processing._run_bcftools_step",
        side_effect=completed_step,
    ) as run_step, patch(
        "postgwas.modules.harmonisation.vcf_processing.get_vcf_variant_count",
        return_value=1,
    ) as count_variants, patch(
        "postgwas.modules.harmonisation.vcf_processing.read_vcf_header",
        return_value=header,
    ):
        annotate_and_liftover_vcf(
            str(tmp_path), "study", "1", "af.vcf.gz", "dbsnp.vcf.gz",
            "GRCh37.fa", "GRCh38.fa", "genes.gff3.gz", "chain.gz",
            "GRCh37", 1, layout,
            {"bash": "/bin/bash", "bcftools": "bcftools", "tabix": "tabix"},
            vcf_config,
        )
        annotate_and_liftover_vcf(
            str(tmp_path), "study", "1", "af.vcf.gz", "dbsnp.vcf.gz",
            "GRCh37.fa", "GRCh38.fa", "genes.gff3.gz", "chain.gz",
            "GRCh37", 1, layout,
            {"bash": "/bin/bash", "bcftools": "bcftools", "tabix": "tabix"},
            vcf_config,
            policies=default_policies().with_overrides({
                "vcf.liftover_swap": "keep",
            }),
        )

    raw_vcf = tmp_path / "study_chr1_GRCh37.vcf.gz"
    assert count_variants.call_args_list[0].args[0] == str(raw_vcf)
    assert not list(tmp_path.glob("*_ORIGINAL.vcf.gz*"))

    id_call = next(
        call for call in run_step.call_args_list if call.args[2] == "STEP1_ID"
    )
    script = id_call.args[0][2]
    assert "--pair-logic exact -Ou" in script
    assert "--set-id '+%CHROM\\_%POS\\_%REF\\_%ALT'" in script
    assert script.count("bcftools annotate") == 2
    assert "bcftools view" not in script

    liftover_calls = [
        call for call in run_step.call_args_list
        if call.args[2] == "STEP4_LIFTOVER"
    ]
    assert len(liftover_calls) == 2
    liftover_script = liftover_calls[0].args[0][2]
    assert (
        "--af-tags INFO/AF,FMT/AF,INFO/AFR,INFO/EAS,INFO/EUR,INFO/SAS"
        in liftover_script
    )
    assert "--es-tags FMT/ES,FMT/EZ" in liftover_script
    assert "INFO/SWAP==1 || INFO/SWAP==-1" in liftover_script
    first_norm = liftover_script.index("bcftools norm")
    sort = liftover_script.index("bcftools sort")
    deduplicate = liftover_script.index("bcftools norm", first_norm + 1)
    assert first_norm < sort < deduplicate
    assert "-d exact" not in liftover_script[first_norm:sort]
    assert "-d exact" in liftover_script[deduplicate:]

    keep_script = liftover_calls[1].args[0][2]
    assert "--es-tags FMT/ES,FMT/EZ" in keep_script
    assert "INFO/SWAP==1 || INFO/SWAP==-1" not in keep_script
