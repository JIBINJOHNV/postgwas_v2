import argparse
from argparse import Namespace
from io import StringIO
from pathlib import Path
import re

import pytest
import polars as pl
from rich.cells import cell_len
from rich.console import Console
import yaml

from postgwas.config import load_configuration, load_module_configuration
from postgwas.core.contracts import RunContext
from postgwas.core.errors import ConfigurationError
from postgwas.core.preflight import PipelinePreflightEvidence
from postgwas.core.resource_preparation import ResourcePreparationError
from postgwas.core.ui import format_cli_default
from postgwas.modules.gcta_gene import cli as gcta_gene_cli
from postgwas.modules.gcta_gene.adapters import build_gcta_command
from postgwas.modules.gcta_gene.cli import build_parser, get_gcta_gene_parser
from postgwas.modules.gcta_gene.errors import GctaGeneError
from postgwas.modules.gcta_gene.reporting import (
    build_gcta_scientific_summary,
    gcta_ma_input_outcome_fields,
)
from postgwas.modules.gcta_gene.results import normalize_gcta_results
from postgwas.modules.gcta_gene.service import (
    GctaGenePipelineResources,
    _gcta_configuration_digest,
    _prepare_analysis_set_list,
    _resolved_configuration,
    preflight_gcta_gene_pipeline,
    run_gcta_gene_direct,
)
from postgwas.modules.gcta_gene.stages import (
    build_gcta_gene_pipeline_progress_plan,
)
from postgwas.pipeline import runners
from postgwas.pipeline.executor import execute_pipeline
from postgwas.pipeline.planner import PipelinePlan
from postgwas.pipeline.runners import run_gcta_gene_runner
from postgwas.pipeline.registry import REGISTRY
from preflight_support import pipeline_input_vcf_evidence


def _write_reference(tmp_path, variant="rs1"):
    prefix = tmp_path / "reference"
    Path(str(prefix) + ".bed").write_bytes(bytes((0x6C, 0x1B, 0x01, 0x00)))
    Path(str(prefix) + ".fam").write_text("F1 I1 0 0 0 -9\n", encoding="utf-8")
    Path(str(prefix) + ".bim").write_text(
        "1\t%s\t0\t100\tG\tA\n" % variant, encoding="utf-8",
    )
    return prefix


def _write_gene_list(tmp_path):
    path = tmp_path / "genes.txt"
    path.write_text("1\t50\t150\tENSG000001\n", encoding="utf-8")
    return path


def _write_set_list(tmp_path, variant="rs1"):
    path = tmp_path / "sets.txt"
    path.write_text("SET_1\n%s\nEND\n" % variant, encoding="utf-8")
    return path


def _write_input(tmp_path, method="fastbat_gene", variant="rs1"):
    path = tmp_path / "input.ma"
    path.write_text(
        "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\n"
        "%s\tG\tA\t0.2\t0.1\t0.05\t0.01\t1000\n" % variant,
        encoding="utf-8",
    )
    return path


def _write_fake_gcta(tmp_path, version="1.94.1", version_exit_code=0):
    path = tmp_path / "fake_gcta"
    path.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        " print('GCTA v%s')\n"
        " raise SystemExit(%d)\n"
        "prefix = sys.argv[sys.argv.index('--out') + 1]\n"
        "if '--fastBAT' in sys.argv:\n"
        " if '--fastBAT-seg' in sys.argv:\n"
        "  Path(prefix + '.log').write_text("
        "'Running fastBAT analysis at genomic segments ...\\n'"
        "+ '1 of 1 sets.\\r')\n"
        "  Path(prefix + '.seg.fastbat').write_text("
        "'Chr Start End No.SNPs SNP_start SNP_end Chisq(Obs) Pvalue TopSNP.Pvalue TopSNP\\n'"
        "+ '1 1 100000 1 rs1 rs1 4.0 0.01 0.01 rs1\\n')\n"
        " elif '--fastBAT-set-list' in sys.argv:\n"
        "  Path(prefix + '.log').write_text("
        "'Running fastBAT analysis ...\\n1 of 1 sets.\\r')\n"
        "  Path(prefix + '.fastbat').write_text("
        "'Set No.SNPs SNP_start SNP_end Chisq(Obs) Pvalue TopSNP.Pvalue TopSNP\\n'"
        "+ 'SET_1 1 rs1 rs1 4.0 0.01 0.01 rs1\\n')\n"
        " else:\n"
        "  Path(prefix + '.log').write_text("
        "'1 genes have been mapped to SNP data.\\n'"
        "+ 'Running fastBAT analysis for genes ...\\n'"
        "+ '1 of 1 genes.\\r')\n"
        "  Path(prefix + '.gene.fastbat').write_text("
        "'Gene Chr Start End No.SNPs SNP_start SNP_end Chisq(Obs) Pvalue TopSNP.Pvalue TopSNP\\n'"
        "+ 'ENSG000001 1 50 150 1 rs1 rs1 4.0 0.01 0.01 rs1\\n')\n"
        "else:\n"
        " Path(prefix + '.log').write_text("
        "'Running mBAT-combo analysis for 1 gene(s) ...\\n')\n"
        " Path(prefix + '.gene.assoc.mbat').write_text("
        "'Gene Chr Start End No.SNPs SNP_start SNP_end TopSNP TopSNP_Pvalue No.Eigenvalues Chisq_mBAT P_mBATcombo P_mBAT Chisq_fastBAT P_fastBAT\\n'"
        "+ 'ENSG000001 1 50 150 1 rs1 rs1 rs1 0.01 1 4.0 0.02 0.03 4.0 0.01\\n')\n"
        % (version, version_exit_code),
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _write_incremental_mbat_gcta(tmp_path, *, fail=False):
    path = tmp_path / ("incremental_gcta_fail" if fail else "incremental_gcta")
    path.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "import sys\n"
        "import time\n"
        "if '--version' in sys.argv:\n"
        " print('GCTA v1.94.1')\n"
        " raise SystemExit(0)\n"
        "prefix = sys.argv[sys.argv.index('--out') + 1]\n"
        "Path(prefix + '.log').write_text("
        "'Running mBAT-combo analysis for 3 gene(s) ...\\n')\n"
        "header = ('Gene Chr Start End No.SNPs SNP_start SNP_end TopSNP '"
        "'TopSNP_Pvalue No.Eigenvalues Chisq_mBAT P_mBATcombo P_mBAT '"
        "'Chisq_fastBAT P_fastBAT\\n')\n"
        "rows = ["
        "f'GENE{i} 1 {i * 100} {i * 100 + 50} 1 rs1 rs1 rs1 0.01 1 4.0 '"
        "f'0.02 0.03 4.0 0.01\\n' for i in range(1, 4)]\n"
        "with Path(prefix + '.gene.assoc.mbat').open('w') as output:\n"
        " output.write(header)\n"
        " output.flush()\n"
        " for index, row in enumerate(rows):\n"
        "  output.write(row)\n"
        "  output.flush()\n"
        "  time.sleep(0.06)\n"
        "  if %r and index == 0:\n"
        "   raise SystemExit(7)\n"
        % fail,
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _write_incremental_fastbat_gcta(tmp_path):
    path = tmp_path / "incremental_fastbat_gcta"
    path.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "import sys\n"
        "import time\n"
        "if '--version' in sys.argv:\n"
        " print('GCTA v1.94.1')\n"
        " raise SystemExit(0)\n"
        "prefix = sys.argv[sys.argv.index('--out') + 1]\n"
        "with Path(prefix + '.log').open('w') as log:\n"
        " log.write('Reading PLINK BED file in SNP-major format ...\\n')\n"
        " log.flush()\n"
        " time.sleep(0.04)\n"
        " log.write('Matching the GWAS meta-analysis results to the genotype data ...\\n')\n"
        " log.flush()\n"
        " time.sleep(0.04)\n"
        " log.write('Mapping the physical positions of genes to SNP data ...\\n')\n"
        " log.flush()\n"
        " time.sleep(0.04)\n"
        " log.write('Running fastBAT analysis for genes ...\\n')\n"
        " log.flush()\n"
        " for completed in (100, 200, 300):\n"
        "  log.write(f'{completed} of 300 genes.\\r')\n"
        "  log.flush()\n"
        "  time.sleep(0.06)\n"
        "Path(prefix + '.gene.fastbat').write_text("
        "'Gene Chr Start End No.SNPs SNP_start SNP_end Chisq(Obs) Pvalue '"
        "'TopSNP.Pvalue TopSNP\\n'"
        "+ ''.join("
        "f'GENE{i} 1 {i * 100} {i * 100 + 50} 1 rs1 rs1 4.0 0.01 0.01 rs1\\n' "
        "for i in range(1, 4)))\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _set_progress_refresh(args, tmp_path, seconds):
    module = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    run_config = tmp_path / "gcta_progress_run.yaml"
    run_config.write_text(
        yaml.safe_dump({
            "logging": {"progress_refresh_seconds": seconds},
            "modules": {"gcta_gene": module},
        }),
        encoding="utf-8",
    )
    args.run_config = str(run_config)


def _configure_three_variant_gmt(args, tmp_path):
    Path(str(tmp_path / "reference") + ".bim").write_text(
        "1\trs1\t0\t100\tG\tA\n"
        "1\trs2\t0\t110\tG\tA\n"
        "1\trs3\t0\t120\tG\tA\n",
        encoding="utf-8",
    )
    (tmp_path / "genes.txt").write_text(
        "1\t50\t150\tGENE1\n", encoding="utf-8",
    )
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("SET_1\tdescription\tGENE1\n", encoding="utf-8")
    config = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    config["set_annotation"] = {"file": None, "gmt_file": str(gmt)}
    Path(args.run_config).write_text(yaml.safe_dump(config), encoding="utf-8")
    _set_progress_refresh(args, tmp_path, 0.01)


def _write_fake_bcftools(
    tmp_path, *, rows=(("1", 100, "rs1", "A", "G"),),
):
    path = tmp_path / "fake_bcftools"
    payload = "\n".join("\t".join(map(str, row)) for row in rows)
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if '-l' in sys.argv:\n"
        " print('STUDY')\n"
        "else:\n"
        " print(%r)\n" % payload,
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _write_config(tmp_path, method, reference, genes, sets):
    path = tmp_path / "gcta_gene.yaml"
    path.write_text(
        "enabled: true\n"
        "method: %s\n"
        "genome_build: GRCh37\n"
        "reference:\n"
        "  prefix: '%s'\n"
        "  population: EUR\n"
        "gene_annotation:\n"
        "  file: '%s'\n"
        "set_annotation:\n"
        "  file: '%s'\n"
        % (method, reference, genes, sets),
        encoding="utf-8",
    )
    return path


def _args(tmp_path, method="fastbat_gene", *, version="1.94.1", variant="rs1"):
    reference = _write_reference(tmp_path)
    genes = _write_gene_list(tmp_path)
    sets = _write_set_list(tmp_path, variant)
    return Namespace(
        gcta_input_file=str(_write_input(tmp_path, method, variant)),
        dataset_id="STUDY",
        output_directory=str(tmp_path / "out"),
        run_config=str(_write_config(tmp_path, method, reference, genes, sets)),
        gcta=str(_write_fake_gcta(tmp_path, version)),
        threads=2,
        dry_run=False,
        resume=False,
        overwrite=False,
    )


def _pipeline_context(args, variant_id_type="rsid"):
    formatted_input = args.gcta_input_file
    del args.gcta_input_file
    return RunContext({
        "formatter": {
            "gcta_gene": {
                "summary_statistics_input_file": formatted_input,
                "variant_id_type": variant_id_type,
            },
        },
    })


def test_help_exposes_scientific_compatibility_and_configuration():
    parser = build_parser()
    help_text = parser.format_help()
    normalized_help = " ".join(help_text.split())
    group_destinations = {
        group.title: [action.dest for action in group._group_actions]
        for group in parser._action_groups
    }
    visible_group_titles = [
        group.title
        for group in parser._action_groups
        if group.description or any(
            action.help is not argparse.SUPPRESS
            for action in group._group_actions
        )
    ]
    assert visible_group_titles[:4] == [
        "options",
        "Output",
        "Input files",
        "Reference files",
    ]
    assert group_destinations["Input files"] == [
        "gcta_input_file",
        "fastbat_set_list",
        "gmt",
    ]
    assert group_destinations["Reference files"] == [
        "gcta_reference_prefix",
        "gene_list",
    ]
    assert group_destinations["GCTA reference settings"] == [
        "genome_build",
        "gcta_reference_population",
    ]
    assert group_destinations["GCTA variant matching"] == [
        "gcta_chromosome_label_policy",
        "gcta_allow_strand_complement",
    ]
    assert group_destinations["GCTA analysis scope"] == [
        "mhc_policy",
        "mhc_chrom",
        "mhc_start",
        "mhc_end",
        "exclude_chromosomes",
    ]
    assert group_destinations["GCTA pathway conversion"] == [
        "gmt_chromosome_label_policy",
        "gmt_duplicate_gene_policy",
        "gmt_unmapped_gene_policy",
        "gcta_minimum_gene_id_overlap",
        "gmt_empty_pathway_policy",
        "fastbat_oversized_set_policy",
    ]
    assert group_destinations["GCTA result reporting"] == [
        "gcta_top_results",
        "gcta_nominal_alpha",
        "gcta_reporting_alpha",
        "gcta_fdr_alpha",
        "gcta_p_value_digits",
    ]
    assert "GCTA binary" not in group_destinations
    assert "GCTA inputs" not in group_destinations
    assert "GCTA association inputs:" not in help_text
    assert "GCTA association settings:" not in help_text
    assert (
        "--method METHOD"
        in help_text
    )
    assert (
        "Available options: fastbat_gene, fastbat_segment, fastbat_set, mbat_combo"
        in normalized_help
    )
    for description in (
        "fastbat_gene: Test SNPs assigned to each gene and its configured window",
        "fastbat_segment: Test consecutive fixed-size genomic segments",
        "fastbat_set: Test custom SNP sets supplied through --fastbat-set-list",
        "mbat_combo: Test each gene by combining signed mBAT and unsigned fastBAT",
    ):
        assert description in normalized_help
        assert any(
            line.lstrip().startswith(description.split(":", 1)[0] + ":")
            for line in help_text.splitlines()
        )
    assert "--fastbat-set-list" in help_text
    assert "--gmt" in help_text
    assert "mutually exclusive with --fastbat-set-list" in normalized_help
    assert "--fastbat-segment-size-kb" in help_text
    assert "--mhc-policy {include,exclude_snps,exclude_genes,exclude_both}" in help_text
    assert "--exclude-chromosomes CHROM [CHROM ...]" in help_text
    assert "(Default: exclude_both)" in normalized_help
    assert "(Default: Y MT)" in normalized_help
    build_actions = [
        action for action in parser._actions if "--genome-build" in action.option_strings
    ]
    assert len(build_actions) == 1
    help_lines = help_text.splitlines()
    build_line = next(
        line
        for line in help_lines
        if line.lstrip().startswith("--genome-build BUILD")
    )
    build_description_column = build_line.index("Required:")
    population_index = next(
        index
        for index, line in enumerate(help_lines)
        if line.lstrip().startswith("--gcta-reference-population CODE")
    )
    population_description = next(
        line
        for line in help_lines[population_index:population_index + 2]
        if "Required:" in line
    )
    assert build_line.lstrip().startswith("--genome-build BUILD")
    assert population_description.index("Required:") == build_description_column
    assert "all coordinate-bearing inputs must match" in normalized_help
    assert "--gcta-reference-genome-build" not in help_text
    assert "--gene-list-genome-build" not in help_text
    assert "--gcta-reference-population" in help_text
    assert "--vcf" not in help_text
    assert "--bcftools" not in help_text
    assert "--gcta PATH" not in help_text
    with pytest.raises(SystemExit):
        parser.parse_args(["--gcta", "gcta64"])
    assert "--gcta-coordinate-fallback" not in help_text
    assert "--gcta-minimum-reference-overlap" not in help_text
    pipeline_help = get_gcta_gene_parser(direct_controls=False).format_help()
    assert "--gcta-coordinate-fallback" not in pipeline_help
    assert "--gcta-minimum-reference-overlap" in pipeline_help
    assert "postgwas config export" in help_text
    for label in (
        "Run gene-based fastBAT:",
        "Run fixed-segment fastBAT:",
        "Run custom-set fastBAT from a prepared set list:",
        "Convert a GMT file and run custom-set fastBAT:",
        "Run gene-based mBAT-combo:",
        "Run entirely from exported settings:",
    ):
        assert label in help_text
    for argument in (
        "--method fastbat_gene",
        "--method fastbat_segment",
        "--fastbat-set-list pathways.fastbat.set",
        "--gmt pathways.gmt",
        "--method mbat_combo",
        "--run-config gcta_gene.yaml",
    ):
        assert argument in help_text


def test_fastbat_set_gmt_pipeline_plan_matches_the_thirteen_stage_contract():
    module = load_configuration().modules.gcta_gene.model_copy(deep=True)
    module.method = "fastbat_set"
    module.set_annotation.file = None
    module.set_annotation.gmt_file = "pathways.gmt"

    plan = build_gcta_gene_pipeline_progress_plan(module)

    assert plan["variant"] == "fastbat_set_gmt"
    assert plan["modules"] == {"formatter": (1, 6), "gcta_gene": (7, 13)}
    assert plan["stages"] == (
        "Validate the input harmonised GWAS-VCF",
        "Validate the PLINK LD-reference files and determine the BIM variant-ID format",
        "Validate the GCTA gene-coordinate file",
        "Validate the original pathway GMT file",
        "Compare pathway genes with the gene-coordinate reference and enforce the configured minimum coordinate-reference coverage",
        "Create and validate the BIM-compatible GCTA .ma input",
        "Map analyzable BIM variants to matched gene intervals",
        "Create candidate pathway-to-variant memberships",
        "Apply duplicate-membership, empty-set and set-size policies and write the final fastBAT sets",
        "Run GCTA fastBAT-set",
        "Validate the raw GCTA results",
        "Add and validate nominal, Bonferroni and Benjamini-Hochberg FDR results",
        "Validate and publish all outputs and display the scientific summary",
    )


def test_gcta_formatter_outcome_reports_every_ma_source_and_transformation(
    tmp_path,
):
    configuration = load_configuration()
    fields = gcta_ma_input_outcome_fields(
        configuration,
        tmp_path / "study_gcta.ma",
        {"variants": 100},
        formatter_result={
            "rows_in": 100,
            "rows_out": 99,
            "rows_excluded": 1,
            "p_values_bounded": 2,
            "variant_id_type": "unique",
        },
        variant_id_type="unique",
    )
    values = [(field[1], field[2]) for field in fields if len(field) == 3]

    assert ("SNP", "CHROM + POS + REF + ALT → CHROM_POS_REF_ALT") in values
    assert ("A1", "ALT → effect allele") in values
    assert ("A2", "REF → other allele") in values
    assert ("freq", "FORMAT/AF → A1 effect-allele frequency") in values
    assert (
        "BETA", "FORMAT/ES → copied without numerical transformation",
    ) in values
    assert (
        "SE", "FORMAT/SE → copied without numerical transformation",
    ) in values
    assert ("P", "FORMAT/LP → raw P = 10⁻ᴸᴾ") in values
    assert ("N", "FORMAT/SS → per-variant total N") in values
    assert ("Variants extracted", 100) in values
    assert ("Variants written", 99) in values
    assert ("Variants excluded during formatting", 1) in values
    assert ("Converted P values bounded at configured minimum", 2) in values

    direct_fields = gcta_ma_input_outcome_fields(
        configuration,
        tmp_path / "supplied.ma",
        {"variants": 100},
    )
    assert (
        "analysis", "Extract and transform GWAS-VCF fields",
    ) not in direct_fields


@pytest.mark.parametrize(
    ("method", "set_source", "expected_variant", "expected_total"),
    (
        ("fastbat_gene", None, "fastbat_gene", 8),
        ("fastbat_segment", None, "fastbat_segment", 7),
        ("fastbat_set", "prepared", "fastbat_set_prepared", 9),
        ("fastbat_set", "gmt", "fastbat_set_gmt", 13),
        ("mbat_combo", None, "mbat_combo", 8),
    ),
)
def test_gcta_pipeline_plans_reuse_common_validation_and_result_stages(
    method, set_source, expected_variant, expected_total,
):
    module = load_configuration().modules.gcta_gene.model_copy(deep=True)
    module.method = method
    module.set_annotation.file = (
        "pathways.fastbat.set" if set_source == "prepared" else None
    )
    module.set_annotation.gmt_file = (
        "pathways.gmt" if set_source == "gmt" else None
    )

    plan = build_gcta_gene_pipeline_progress_plan(module)

    assert plan["variant"] == expected_variant
    assert len(plan["stages"]) == expected_total
    assert plan["stages"][:2] == (
        "Validate the input harmonised GWAS-VCF",
        "Validate the PLINK LD-reference files and determine the BIM variant-ID format",
    )
    assert plan["stages"][-3:] == (
        "Validate the raw GCTA results",
        "Add and validate nominal, Bonferroni and Benjamini-Hochberg FDR results",
        "Validate and publish all outputs and display the scientific summary",
    )


def test_configurable_cli_values_do_not_own_defaults():
    parser = build_parser()
    for destination in (
        "gcta_input_file", "gcta_gene_method", "gcta_reference_prefix",
        "gene_list", "fastbat_set_list", "gmt", "genome_build",
        "gcta_reference_population", "gene_window_kb", "fastbat_segment_size_kb",
        "gmt_chromosome_label_policy", "gmt_duplicate_gene_policy",
        "gmt_unmapped_gene_policy", "gcta_minimum_gene_id_overlap",
        "gmt_empty_pathway_policy",
        "fastbat_oversized_set_policy",
        "gcta_reference_maf_min", "fastbat_ld_cutoff", "mbat_svd_gamma",
        "frequency_difference_max", "print_component_p_values", "write_snpset",
        "gcta_top_results", "gcta_nominal_alpha", "gcta_reporting_alpha",
        "gcta_fdr_alpha", "gcta_p_value_digits",
        "run_config", "resume", "overwrite",
        "dry_run", "dataset_id",
        "output_directory", "threads", "memory_gb", "seed",
        "gcta_chromosome_label_policy", "gcta_allow_strand_complement",
        "mhc_policy", "mhc_chrom", "mhc_start", "mhc_end",
        "exclude_chromosomes",
    ):
        action = next(item for item in parser._actions if item.dest == destination)
        assert action.default == argparse.SUPPRESS


def test_cli_help_defaults_are_read_from_canonical_configuration(monkeypatch):
    defaults = load_configuration()
    module = defaults.modules.gcta_gene
    module.method = "mbat_combo"
    module.gene_window_kb = 73
    module.fastbat_ld_cutoff = 0.83
    module.print_component_p_values = False
    module.reporting.top_result_count = 7
    module.reporting.nominal_alpha = 0.02
    module.reporting.familywise_alpha = 0.01
    module.reporting.fdr_alpha = 0.03
    module.reporting.p_value_significant_digits = 6
    module.set_annotation.conversion.chromosome_label_policy = "strip_chr_prefix"
    module.set_annotation.maximum_set_variants = 12_345
    module.mhc.policy = "exclude_snps"
    module.chromosomes.exclude = ["Y"]
    module.reference.required_extensions = [".bim", ".fam", ".bed"]
    monkeypatch.setattr(gcta_gene_cli, "load_configuration", lambda: defaults)

    parser = gcta_gene_cli.get_gcta_gene_parser(direct_controls=True)
    expected_by_destination = {
        "gcta_gene_method": "mbat_combo",
        "gene_window_kb": 73,
        "fastbat_ld_cutoff": 0.83,
        "print_component_p_values": False,
        "gcta_top_results": 7,
        "gcta_nominal_alpha": 0.02,
        "gcta_reporting_alpha": 0.01,
        "gcta_fdr_alpha": 0.03,
        "gcta_p_value_digits": 6,
        "gmt_chromosome_label_policy": "strip_chr_prefix",
        "gcta_minimum_gene_id_overlap": 0.5,
        "fastbat_oversized_set_policy": "omit",
        "mhc_policy": "exclude_snps",
        "exclude_chromosomes": "Y",
    }
    for destination, expected in expected_by_destination.items():
        action = next(item for item in parser._actions if item.dest == destination)
        assert format_cli_default(expected) in action.help
    rendered_help = " ".join(parser.format_help().split())
    assert "configured companion files: .bim, .fam, .bed" in rendered_help
    assert "configured GCTA limit of 12345 variants" in rendered_help


def test_output_help_uses_canonical_run_defaults_without_argparse_defaults():
    defaults = load_configuration()
    parser = build_parser()
    expected_by_destination = {
        "dataset_id": defaults.run.dataset_id,
        "output_directory": defaults.run.output_directory,
    }
    for destination, expected in expected_by_destination.items():
        action = next(item for item in parser._actions if item.dest == destination)
        assert format_cli_default(expected) in action.help
        assert action.default == argparse.SUPPRESS


def test_mapping_worker_memory_policy_does_not_invalidate_scientific_digest():
    baseline = load_configuration()
    changed = baseline.model_copy(deep=True)
    parallelism = (
        changed.modules.gcta_gene.set_annotation.conversion.parallelism
    )
    parallelism.worker_memory_multiplier = 9

    assert _gcta_configuration_digest(changed) == _gcta_configuration_digest(
        baseline
    )


def test_pathway_audit_and_disk_defaults_are_schema_validated_and_non_scientific():
    baseline = load_configuration()
    conversion = baseline.modules.gcta_gene.set_annotation.conversion

    assert conversion.audit.level == "normalized"
    assert conversion.audit.format == "parquet"
    assert conversion.audit.compression == "zstd"
    assert conversion.audit.batch_rows > 0
    assert conversion.disk.minimum_free_gb >= 0
    assert conversion.disk.estimation_safety_factor >= 1

    changed = baseline.model_copy(deep=True)
    changed_conversion = changed.modules.gcta_gene.set_annotation.conversion
    changed_conversion.audit.level = "summary"
    changed_conversion.disk.minimum_free_gb += 1
    assert _gcta_configuration_digest(changed) == _gcta_configuration_digest(
        baseline
    )


def test_maximum_set_variants_remains_in_scientific_digest():
    baseline = load_configuration()
    changed = baseline.model_copy(deep=True)
    changed.modules.gcta_gene.set_annotation.maximum_set_variants -= 1

    assert _gcta_configuration_digest(changed) != _gcta_configuration_digest(
        baseline
    )


def test_missing_direct_inputs_are_reported_together_logged_and_show_help(
    tmp_path, monkeypatch, capsys,
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.resolve_executable",
        lambda *args, **kwargs: pytest.fail(
            "runtime validation must follow required-argument validation"
        ),
    )

    assert gcta_gene_cli.main([]) == 1

    captured = capsys.readouterr()
    terminal = " ".join((captured.out + captured.err).split())
    required_options = (
        "--gcta-input-file",
        "--gcta-reference-prefix",
        "--genome-build",
        "--gcta-reference-population",
        "--gene-list",
    )
    for option in required_options:
        assert "Required argument not provided: %s." % option in terminal
    assert "usage: postgwas gcta_gene" in terminal.lower()
    assert "Run gene-based fastBAT:" in terminal

    service_log = (
        tmp_path / "results/logs/postgwas_fastbat_gene_gcta_gene.log"
    )
    log_text = " ".join(
        re.sub(
            r"(?:^|\n)\[[^]]+\]\s+RUN\s+(?:FAILED\s+)?",
            " ",
            service_log.read_text(encoding="utf-8"),
        ).split()
    )
    for option in required_options:
        assert "Required argument not provided: %s." % option in log_text
    assert not (tmp_path / "results/raw").exists()


def test_pipeline_runner_defers_complete_cli_validation_to_service(
    tmp_path, monkeypatch,
):
    resources = GctaGenePipelineResources(
        configuration=load_configuration(),
        reference_prefix=tmp_path / "reference",
        reference_paths=(),
        executable="gcta64",
        version="1.94.1",
        identifier_observation={},
        analysis_scope={},
        gene_list=None,
        gene_annotation_metrics=None,
        gmt=None,
        pathway_gene_preflight=None,
        set_list=None,
        set_source_metrics=None,
        file_identities=(),
    )
    context = RunContext(
        {"formatter": {"gcta_gene": {}}},
        validations={
            "gcta_gene": PipelinePreflightEvidence(
                module="gcta_gene",
                input_vcf={},
                resources=resources,
            ),
        },
    )
    args = Namespace(
        output_directory=str(tmp_path / "pipeline"),
        _step_num="02",
        gcta_gene_method="fastbat_set",
    )
    captured = {}
    expected = object()

    def reject_partial_configuration(*_args, **_kwargs):
        raise AssertionError("runner performed partial GCTA configuration validation")

    def fake_service(service_args, service_context, *, pipeline_resources):
        captured["has_input_override"] = hasattr(service_args, "gcta_input_file")
        captured["output_directory"] = service_args.output_directory
        captured["context"] = service_context
        captured["pipeline_resources"] = pipeline_resources
        return expected

    monkeypatch.setattr(
        "postgwas.pipeline.runners.load_module_configuration",
        reject_partial_configuration,
    )
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.run_gcta_gene_direct",
        fake_service,
    )

    result = run_gcta_gene_runner(args, context)

    assert result is expected
    assert captured == {
        "has_input_override": False,
        "output_directory": str(tmp_path / "pipeline" / "02_gcta_gene"),
        "context": context,
        "pipeline_resources": resources,
    }
    assert args.output_directory == str(tmp_path / "pipeline")


def test_pipeline_preflight_validates_external_resources_before_formatter(
    tmp_path,
):
    args = _args(tmp_path)
    del args.gcta_input_file

    evidence = preflight_gcta_gene_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )

    assert REGISTRY.get("gcta_gene").preflight.endswith(
        ":preflight_gcta_gene_pipeline"
    )
    assert isinstance(evidence, PipelinePreflightEvidence)
    assert isinstance(evidence.resources, GctaGenePipelineResources)
    assert evidence.resources.version == "1.94.1"
    assert evidence.resources.identifier_observation["variant_id_type"] == "rsid"
    assert evidence.resources.gene_annotation_metrics["genes"] == 1
    assert len(evidence.resources.file_identities) == 5
    assert "formatter-created GCTA .ma" in " ".join(evidence.deferred_checks)
    assert not (tmp_path / "out").exists()


def test_pipeline_preflight_rejects_malformed_gene_resource_before_output(
    tmp_path,
):
    args = _args(tmp_path)
    del args.gcta_input_file
    (tmp_path / "genes.txt").write_text(
        "1\t50\t150\tENSG000001\n1\t60\t160\tENSG000001\n",
        encoding="utf-8",
    )

    with pytest.raises(GctaGeneError, match="duplicate gene ID"):
        preflight_gcta_gene_pipeline(
            args,
            preflight_evidence=pipeline_input_vcf_evidence(),
        )

    assert not (tmp_path / "out").exists()


def test_pipeline_execution_reuses_static_gcta_gene_preflight(
    tmp_path, monkeypatch,
):
    args = _args(tmp_path)
    evidence = preflight_gcta_gene_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )
    context = _pipeline_context(args)
    context.publish_validation("gcta_gene", evidence)
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.validate_gcta_reference_files",
        lambda *_args, **_kwargs: pytest.fail(
            "PLINK resource validation was repeated"
        ),
    )
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.require_supported_gcta",
        lambda *_args, **_kwargs: pytest.fail("GCTA version probe was repeated"),
    )
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service._validate_gene_list",
        lambda *_args, **_kwargs: pytest.fail(
            "gene-coordinate validation was repeated"
        ),
    )

    result = run_gcta_gene_direct(
        args,
        context,
        pipeline_resources=evidence.resources,
    )

    assert result.metrics["tested_units"] == 1


def test_pipeline_rejects_resource_changed_after_gcta_gene_preflight(tmp_path):
    args = _args(tmp_path)
    evidence = preflight_gcta_gene_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )
    context = _pipeline_context(args)
    context.publish_validation("gcta_gene", evidence)
    with (tmp_path / "genes.txt").open("a", encoding="utf-8") as handle:
        handle.write("1\t200\t300\tENSG000002\n")

    with pytest.raises(GctaGeneError, match="changed after pipeline preflight"):
        run_gcta_gene_direct(
            args,
            context,
            pipeline_resources=evidence.resources,
        )

    assert not (tmp_path / "out").exists()


def test_fastbat_set_gmt_pipeline_reports_all_thirteen_stages_and_log_summaries(
    tmp_path, monkeypatch,
):
    args = _args(tmp_path, "fastbat_set")
    formatter_input = Path(args.gcta_input_file)
    del args.gcta_input_file
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("SET_1\tdescription\tENSG000001\n", encoding="utf-8")
    args.gmt = str(gmt)
    args.modules = ["gcta_gene"]
    args.format = None
    args.vcf = str(tmp_path / "study.vcf.gz")
    Path(args.vcf).write_text("vcf\n", encoding="utf-8")
    args.bcftools = "bcftools"
    args._step_num = "01"
    module_values = yaml.safe_load(
        Path(args.run_config).read_text(encoding="utf-8")
    )
    Path(args.run_config).write_text(
        yaml.safe_dump({"modules": {"gcta_gene": module_values}}),
        encoding="utf-8",
    )
    configuration = _resolved_configuration(args)
    plan = build_gcta_gene_pipeline_progress_plan(
        configuration.modules.gcta_gene,
    )
    stream = StringIO()
    progress_console = Console(
        file=stream,
        force_terminal=False,
        color_system=None,
        width=140,
    )
    monkeypatch.setattr(
        "postgwas.pipeline.executor.console", progress_console,
    )

    input_evidence = pipeline_input_vcf_evidence()
    input_evidence["input_vcf"]["harmonised"]["postgwas_dataset_id"] = "STUDY"
    gcta_evidence = preflight_gcta_gene_pipeline(
        args,
        preflight_evidence=input_evidence,
    )
    current_evidence = input_evidence["input_vcf"]
    monkeypatch.setattr(
        runners,
        "_validate_current_pipeline_vcf",
        lambda *_args, **_kwargs: current_evidence,
    )
    monkeypatch.setattr(runners, "check_and_resolve_binaries", lambda *_: None)

    def configure_identifiers(namespace, *_args):
        namespace.variant_id_types = {"gcta_gene": "rsid"}
        namespace.variant_id_observations = {
            "gcta_gene": {
                "variant_id_type": "rsid",
                "variants": 1,
                "consumers": ["GCTA gene"],
            },
        }

    monkeypatch.setattr(
        runners,
        "configure_reference_variant_identifiers",
        configure_identifiers,
    )

    def fake_formatter(formatter_args, context):
        output_file = Path(formatter_args.output_directory) / formatter_input.name
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_bytes(formatter_input.read_bytes())
        result = {
            "gcta_gene": {
                "summary_statistics_input_file": str(output_file),
                "variant_id_type": "rsid",
                "rows_in": 1,
                "rows_out": 1,
                "rows_excluded": 0,
                "p_values_bounded": 0,
            },
        }
        context["formatter"] = result
        return result

    monkeypatch.setattr(
        "postgwas.modules.formatting.service.run_formatter_direct",
        fake_formatter,
    )
    context = execute_pipeline(
        args,
        PipelinePlan(
            ("gcta_gene",),
            ("formatter", "gcta_gene"),
            ("formatter", "gcta_gene"),
        ),
        configuration,
        preflight_evidence={
            **input_evidence,
            "current_vcf": current_evidence,
            "gcta_gene": gcta_evidence,
        },
    )
    result = context["gcta_gene"]

    assert result.metrics["tested_units"] == 1
    output = stream.getvalue()
    normalized_output = " ".join(output.split())
    for number, title in enumerate(plan["stages"], 1):
        assert "Completed %d/13 · %s" % (number, title) in output
    for label in (
        "Input summary-statistics VCF",
        "PLINK LD reference",
        "Required companion files",
        "BIM structural validation",
        "GCTA gene-coordinate reference",
        "Coordinate columns",
        "Original pathway GMT",
        "Detected format",
        "Pathway/gene-coordinate compatibility",
        "Required minimum coordinate-reference coverage",
        "Prepare GCTA summary statistics from the harmonised GWAS-VCF",
        "Extract and transform GWAS-VCF fields",
        "Variants extracted",
        "Write the GCTA .ma input",
        "Required columns and non-missing values",
        "Scientific numeric ranges",
        "Compare the final .ma input with PLINK BIM",
        "Compatible shared allele pairs",
    ):
        assert label in normalized_output
    for provenance in (
        "ID → SNP",
        "ALT → effect allele",
        "REF → other allele",
        "FORMAT/AF → A1 effect-allele frequency",
        "FORMAT/ES → copied without numerical transformation",
        "FORMAT/SE → copied without numerical transformation",
        "FORMAT/LP → raw P = 10⁻ᴸᴾ",
        "FORMAT/SS → per-variant total N",
    ):
        assert provenance in normalized_output
    assert "All 13 stages completed" in output
    log = (
        tmp_path
        / "out"
        / "02_gcta_gene"
        / "logs"
        / "STUDY_fastbat_set_gcta_gene.log"
    ).read_text(encoding="utf-8")
    normalized_log = " ".join(
        re.sub(r"^\[[^]]+\] RUN\s*", "", line).strip()
        for line in log.splitlines()
    )
    for number, title in enumerate(plan["stages"], 1):
        assert "step=%d" % number in log
        assert title in normalized_log
    assert "gcta_input_validation_summary" in log
    assert "input_rewritten=false" in log
    for label in (
        "Required companion files",
        "BIM structural validation",
        "Coordinate columns",
        "Detected format",
        "Required minimum coordinate-reference coverage",
        "Variants extracted",
        "Required columns and non-missing values",
        "Compatible shared allele pairs",
    ):
        assert label in normalized_log
    for provenance in (
        "ID → SNP",
        "ALT → effect allele",
        "REF → other allele",
        "FORMAT/AF → A1 effect-allele frequency",
        "FORMAT/ES → copied without numerical transformation",
        "FORMAT/SE → copied without numerical transformation",
        "FORMAT/LP → raw P = 10⁻ᴸᴾ",
        "FORMAT/SS → per-variant total N",
    ):
        assert provenance in normalized_log
    assert "status=VALIDATED_PENDING_PIPELINE_CHECKPOINT" in log
    assert "gcta_pipeline_stage step=13 total=13" in log
    assert "status=COMPLETED" in log
    assert not hasattr(args, "_pipeline_progress_completion")
    assert not hasattr(args, "_pipeline_progress_failure")


def test_pipeline_cli_path_values_resolve_as_canonical_configuration(tmp_path):
    reference = _write_reference(tmp_path)
    gene_list = _write_gene_list(tmp_path)
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("SET_1\tdescription\tENSG000001\n", encoding="utf-8")
    output = tmp_path / "pipeline"
    args = Namespace(
        gcta_gene_method="fastbat_set",
        gcta_reference_prefix=str(reference),
        gene_list=gene_list,
        gmt=gmt,
        genome_build="GRCh37",
        gcta_reference_population="EUR",
        gcta_minimum_gene_id_overlap=0.75,
        dataset_id="STUDY",
        output_directory=output,
    )

    configuration = _resolved_configuration(args)
    module = configuration.modules.gcta_gene

    assert module.method == "fastbat_set"
    assert module.gene_annotation.file == str(gene_list)
    assert module.set_annotation.gmt_file == str(gmt)
    assert (
        module.set_annotation.conversion.minimum_gene_id_overlap_fraction
        == 0.75
    )
    assert configuration.run.dataset_id == "STUDY"
    assert configuration.run.output_directory == output


def test_gcta_scope_cli_overrides_canonical_yaml(tmp_path):
    reference = _write_reference(tmp_path)
    gene_list = _write_gene_list(tmp_path)
    configuration = _resolved_configuration(Namespace(
        gcta_reference_prefix=str(reference),
        gene_list=gene_list,
        genome_build="GRCh37",
        gcta_reference_population="EUR",
        mhc_policy="exclude_genes",
        mhc_chrom="chr6",
        mhc_start=100,
        mhc_end=200,
        exclude_chromosomes=["MT"],
    ))

    module = configuration.modules.gcta_gene
    assert module.mhc.policy == "exclude_genes"
    assert module.mhc.region_override.model_dump() == {
        "chromosome": "chr6", "start": 100, "end": 200,
    }
    assert module.chromosomes.exclude == ["MT"]


@pytest.mark.parametrize(
    "method",
    [
        "fastbat_set",
        "mbat_combo",
    ],
)
def test_pipeline_service_selects_shared_formatted_input(tmp_path, method):
    args = _args(tmp_path, method)
    context = _pipeline_context(args)

    result = run_gcta_gene_direct(args, context)

    assert result.metrics["method"] == method
    assert context["gcta_gene"] is result


def test_fastbat_command_uses_only_official_gene_test_options():
    module = load_module_configuration("gcta_gene")
    command = build_gcta_command(
        "gcta64", "study.tsv", "reference", "genes.txt", "out/STUDY", 4, module,
    )
    assert command == [
        "gcta64", "--bfile", "reference", "--maf", "0.01",
        "--diff-freq", "0.2", "--thread-num", "4",
        "--out", "out/STUDY", "--fastBAT", "study.tsv", "--fastBAT-ld-cutoff",
        "0.9", "--fastBAT-gene-list", "genes.txt", "--fastBAT-wind", "50",
    ]
    assert "--mBAT-combo" not in command


def test_gcta_command_uses_documented_exact_id_exclusion_file():
    module = load_module_configuration("gcta_gene")
    command = build_gcta_command(
        "gcta64",
        "study.tsv",
        "reference",
        "genes.txt",
        "out/STUDY",
        4,
        module,
        exclude_variants_file="excluded.snplist",
    )

    assert command[command.index("--exclude") + 1] == "excluded.snplist"
    assert command.index("--exclude") < command.index("--fastBAT")


def test_fastbat_segment_command_uses_official_segment_option(tmp_path):
    config = tmp_path / "segment.yaml"
    config.write_text("method: fastbat_segment\n", encoding="utf-8")
    module = load_module_configuration("gcta_gene", config)
    command = build_gcta_command(
        "gcta64", "study.tsv", "reference", None, "out/STUDY", 4, module,
    )
    assert "--fastBAT-seg" in command
    assert command[command.index("--fastBAT-seg") + 1] == "100"
    assert "--fastBAT-gene-list" not in command
    assert "--fastBAT-set-list" not in command


def test_fastbat_set_command_uses_official_set_and_snpset_options(tmp_path):
    config = tmp_path / "set.yaml"
    config.write_text(
        "method: fastbat_set\nwrite_snpset: true\n", encoding="utf-8",
    )
    module = load_module_configuration("gcta_gene", config)
    command = build_gcta_command(
        "gcta64", "study.tsv", "reference", "sets.txt", "out/STUDY", 4,
        module,
    )
    assert "--fastBAT-set-list" in command
    assert "--fastBAT-write-snpset" in command
    assert "--fastBAT-gene-list" not in command


def test_mbat_combo_command_includes_every_documented_method_option(tmp_path):
    config = tmp_path / "mbat.yaml"
    config.write_text(
        "method: mbat_combo\nwrite_snpset: true\n", encoding="utf-8",
    )
    module = load_module_configuration("gcta_gene", config)
    command = build_gcta_command(
        "gcta64", "study.ma", "reference", "genes.txt", "out/STUDY", 4,
        module,
    )
    assert command == [
        "gcta64", "--bfile", "reference", "--maf", "0.01",
        "--diff-freq", "0.2", "--thread-num", "4",
        "--out", "out/STUDY", "--mBAT-combo", "study.ma", "--mBAT-gene-list",
        "genes.txt", "--mBAT-wind", "50", "--mBAT-svd-gamma", "0.9",
        "--fastBAT-ld-cutoff", "0.9",
        "--mBAT-print-all-p", "--mBAT-write-snpset",
    ]


@pytest.mark.parametrize(
    "method", ["fastbat_gene", "fastbat_segment", "fastbat_set", "mbat_combo"],
)
@pytest.mark.parametrize("threshold", [0.0001, 0.1, 1.0])
def test_all_gcta_methods_pass_the_configured_frequency_threshold(method, threshold):
    module = load_module_configuration(
        "gcta_gene",
        cli_overrides={"method": method, "frequency_difference_max": threshold},
    )
    command = build_gcta_command(
        "gcta64", "study.ma", "reference", "annotation.txt", "out/STUDY", 2,
        module,
    )

    assert command.count("--diff-freq") == 1
    assert command[command.index("--diff-freq") + 1] == str(threshold)


@pytest.mark.parametrize(
    "method", ["fastbat_gene", "fastbat_segment", "fastbat_set", "mbat_combo"],
)
def test_all_gcta_methods_expose_and_record_frequency_threshold(tmp_path, method):
    parser = get_gcta_gene_parser(direct_controls=False)
    parser._postgwas_help_arguments = ("--method", method)
    gcta_gene_cli._customize_gcta_gene_pipeline_help(parser)
    assert "--frequency-difference-max FRACTION" in parser.format_help()

    args = _args(tmp_path, method)
    args.frequency_difference_max = 0.1
    configuration = _resolved_configuration(args)
    baseline = configuration.model_copy(deep=True)
    baseline.modules.gcta_gene.frequency_difference_max = 0.2
    assert _gcta_configuration_digest(configuration) != _gcta_configuration_digest(
        baseline
    )
    result = run_gcta_gene_direct(args)
    report = result.artifacts["html_report"].path.read_text(encoding="utf-8")
    assert "<td>Maximum frequency difference</td><td>0.1</td>" in report
    log = tmp_path / "out" / "logs" / ("STUDY_%s_gcta_gene.log" % method)
    assert "frequency_difference_max=0.1" in log.read_text(encoding="utf-8")


def test_gcta_input_accepts_distinct_nonempty_indel_alleles(tmp_path):
    args = _args(tmp_path, "fastbat_gene")
    Path(args.gcta_input_file).write_text(
        "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\n"
        "rs1\tCAAATAAAT\tC\t0.2\t0.1\t0.05\t0.01\t1000\n",
        encoding="utf-8",
    )
    Path(str(tmp_path / "reference") + ".bim").write_text(
        "1\trs1\t0\t100\tC\tCAAATAAAT\n", encoding="utf-8",
    )

    result = run_gcta_gene_direct(args)

    assert result.metrics["compatible_allele_pairs"] == 1


def test_mbat_combo_runs_once_and_reports_internal_fastbat(tmp_path, capsys):
    args = _args(tmp_path, "mbat_combo")
    result = run_gcta_gene_direct(args)

    assert result.metrics["method"] == "mbat_combo"
    assert result.metrics["tested_units"] == 1
    assert result.metrics["minimum_p_value"] == pytest.approx(0.02)
    top = result.metrics["scientific_summary"]["top_associations"][0]
    assert top["p_value"] == pytest.approx(0.02)
    assert top["component_p_values"] == {
        "P_mBAT": pytest.approx(0.03),
        "P_fastBAT": pytest.approx(0.01),
    }
    raw = result.artifacts["raw_results"].path
    normalized = result.artifacts["normalized_results"].path
    assert raw.name == "STUDY.gene.assoc.mbat" and raw.is_file()
    assert "P_fastBAT" in normalized.read_text(encoding="utf-8").splitlines()[0]
    html_report = result.artifacts["html_report"].path
    report_text = html_report.read_text(encoding="utf-8")
    assert html_report.name == "STUDY_mbat_combo_gcta_gene_report.html"
    assert "GCTA mBAT-combo gene-association report" in report_text
    assert "Complete association results" in report_text
    assert "P_mBATcombo" in report_text
    assert "P_mBAT" in report_text
    assert "P_fastBAT" in report_text
    log = tmp_path / "out" / "logs" / "STUDY_mbat_combo_gcta_gene.log"
    log_text = log.read_text(encoding="utf-8")
    assert log_text.count("--mBAT-combo") == 1
    assert "--mBAT-print-all-p" in log_text
    assert "--fastBAT," not in log_text
    assert "gcta_scientific_summary" in log_text
    assert "gcta_top_association" in log_text
    terminal = capsys.readouterr().out
    assert "GCTA mBAT-combo gene-association summary" in terminal
    assert "P_mBATcombo 0.02" in terminal
    assert "P_fdr_bh 0.02" in terminal
    assert "P_bonferroni 0.02" in terminal
    assert "P_mBAT 0.03" in terminal
    assert "P_fastBAT 0.01" in terminal
    assert "GCTA mbat combo measured progress" in terminal
    assert "Progress 0/1 · Test mapped genes · 0%" in terminal
    assert (
        "Completed 1/1 · Validate 1 result row for genes"
        in terminal
    )
    assert "gcta_progress_plan" in log_text
    normalized_log = " ".join(
        re.sub(
            r"(?:^|\n)\[[^]]+\]\s+RUN\s+(?:FAILED\s+)?",
            " ",
            log_text,
        ).split()
    )
    assert "percentage=100 status=VALIDATED" in normalized_log
    assert not list((tmp_path / "out" / ".partial").glob("*"))
    resolved = yaml.safe_load(
        (tmp_path / "out" / "run_metadata" / "resolved_config.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert list(resolved["modules"]) == ["formatting", "gcta_gene"]
    assert set(resolved["resources"]) == {
        "executables", "genomes", "populations",
    }
    assert list(resolved["resources"]["genomes"]) == ["GRCh37"]
    assert list(resolved["resources"]["populations"]) == ["EUR"]


def test_mbat_progress_tracks_incrementally_appended_gene_rows(tmp_path, capsys):
    args = _args(tmp_path, "mbat_combo")
    args.gcta = str(_write_incremental_mbat_gcta(tmp_path))
    _set_progress_refresh(args, tmp_path, 0.01)

    result = run_gcta_gene_direct(args)

    assert result.metrics["tested_units"] == 3
    terminal = capsys.readouterr().out
    assert "Progress 1/3 · Test mapped genes · 33%" in terminal
    assert "Progress 2/3 · Test mapped genes · 66%" in terminal
    assert (
        "Completed 3/3 · Validate 3 result rows for genes"
        in terminal
    )
    assert terminal.index("Progress 1/3") < terminal.index("Progress 2/3")
    assert terminal.index("Progress 2/3") < terminal.index("Completed 3/3")
    log = tmp_path / "out/logs/STUDY_mbat_combo_gcta_gene.log"
    log_text = log.read_text(encoding="utf-8")
    normalized_log = " ".join(
        re.sub(
            r"(?:^|\n)\[[^]]+\]\s+RUN\s+(?:FAILED\s+)?",
            " ",
            log_text,
        ).split()
    )
    assert "completed=1 total=3 percentage=33 status=RUNNING" in normalized_log
    assert "completed=2 total=3 percentage=66 status=RUNNING" in normalized_log
    assert (
        "completed=3 total=3 validated_result_units=3 percentage=100 "
        "status=VALIDATED"
    ) in normalized_log


def test_fastbat_progress_uses_native_compute_counter_before_results_exist(
    tmp_path, capsys,
):
    args = _args(tmp_path, "fastbat_gene")
    args.gcta = str(_write_incremental_fastbat_gcta(tmp_path))
    _set_progress_refresh(args, tmp_path, 0.01)

    result = run_gcta_gene_direct(args)

    assert result.metrics["tested_units"] == 3
    terminal = capsys.readouterr().out
    phases = (
        "Started 0/? · Load and filter LD-reference and GWAS variants",
        "Current 0/? · Match GWAS variants to the LD reference",
        "Current 0/? · Map gene positions to LD-reference variants",
        "Current 0/? · Evaluate annotated genes with fastBAT",
    )
    for phase in phases:
        assert phase in terminal
    assert [terminal.index(phase) for phase in phases] == sorted(
        terminal.index(phase) for phase in phases
    )
    assert "Progress 100/300 · Evaluate annotated genes with fastBAT · 33%" in terminal
    assert "Progress 200/300 · Evaluate annotated genes with fastBAT · 66%" in terminal
    assert (
        "Completed 300/300 · Validate 3 result rows for genes"
        in terminal
    )
    assert terminal.index("Progress 100/300") < terminal.index("Progress 200/300")
    assert terminal.index("Progress 200/300") < terminal.index("Completed 300/300")
    log = tmp_path / "out/logs/STUDY_fastbat_gene_gcta_gene.log"
    normalized_log = " ".join(
        re.sub(
            r"(?:^|\n)\[[^]]+\]\s+RUN\s+(?:FAILED\s+)?",
            " ",
            log.read_text(encoding="utf-8"),
        ).split()
    )
    assert "gcta_execution_phase" in normalized_log
    for phase in ("load_variants", "match_variants", "map_genes", "analyse_units"):
        assert "phase=%s" % phase in normalized_log
    assert "completed=100 total=300 percentage=33 status=RUNNING" in normalized_log
    assert "completed=200 total=300 percentage=66 status=RUNNING" in normalized_log
    assert "validated_result_units=3 percentage=100 status=VALIDATED" in normalized_log


def test_mbat_progress_failure_stays_below_completion(tmp_path, capsys):
    args = _args(tmp_path, "mbat_combo")
    args.gcta = str(_write_incremental_mbat_gcta(tmp_path, fail=True))
    _set_progress_refresh(args, tmp_path, 0.01)

    with pytest.raises(GctaGeneError, match="exit status 7"):
        run_gcta_gene_direct(args)

    terminal = capsys.readouterr().out
    assert "Failed 1/3 · Test mapped genes" in terminal
    assert "Completed 3/3" not in terminal
    log = tmp_path / "out/logs/STUDY_mbat_combo_gcta_gene.log"
    log_text = log.read_text(encoding="utf-8")
    assert "completed=1 total=3 status=FAILED" in log_text
    assert "percentage=100 status=VALIDATED" not in log_text


def test_fastbat_run_validates_and_normalizes_results(tmp_path):
    result = run_gcta_gene_direct(_args(tmp_path))
    assert result.metrics["gcta_version"] == "1.94.1"
    assert result.metrics["overlapping_variants"] == 1
    assert result.metrics["tested_units"] == 1
    assert result.artifacts["raw_results"].path.name == "STUDY.gene.fastbat"
    assert result.artifacts["normalized_results"].path.is_file()
    assert result.artifacts["html_report"].path.is_file()
    assert result.artifacts["completion_manifest"].path.is_file()
    summary = result.metrics["scientific_summary"]
    assert summary["tested_units"] == 1
    assert summary["bonferroni_threshold"] == pytest.approx(0.05)
    assert summary["nominal_significant"] == 1
    assert summary["fdr_bh_significant"] == 1
    assert summary["bonferroni_significant"] == 1
    corrected = pl.read_csv(
        result.artifacts["normalized_results"].path,
        separator="\t",
    )
    assert corrected["P_bonferroni"].to_list() == pytest.approx([0.01])
    assert corrected["P_fdr_bh"].to_list() == pytest.approx([0.01])
    assert corrected["significant_nominal"].to_list() == [True]
    assert corrected["significant_bonferroni"].to_list() == [True]
    assert corrected["significant_fdr_bh"].to_list() == [True]
    report_text = result.artifacts["html_report"].path.read_text(encoding="utf-8")
    for heading in (
        "Results at a glance",
        "Multiple-testing results",
        "Input compatibility and analysis coverage",
        "Complete association results",
        "How to interpret",
        "Analysis parameters and provenance",
    ):
        assert heading in report_text
    assert "All 1 validated genes are included" in report_text
    assert "ENSG000001" in report_text
    assert 'data-page-size="50"' in report_text


def test_gcta_default_scope_excludes_mhc_and_configured_chromosomes(tmp_path):
    args = _args(tmp_path)
    Path(args.gcta_input_file).write_text(
        "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\n"
        "rs1\tG\tA\t0.2\t0.1\t0.05\t0.01\t1000\n"
        "rs_mhc\tG\tA\t0.2\t0.1\t0.05\t0.02\t1000\n"
        "rs_y\tG\tA\t0.2\t0.1\t0.05\t0.03\t1000\n",
        encoding="utf-8",
    )
    Path(str(tmp_path / "reference") + ".bim").write_text(
        "1\trs1\t0\t100\tG\tA\n"
        "6\trs_mhc\t0\t29000000\tG\tA\n"
        "Y\trs_y\t0\t100\tG\tA\n",
        encoding="utf-8",
    )
    (tmp_path / "genes.txt").write_text(
        "1\t50\t150\tENSG000001\n"
        "6\t28450000\t28470000\tENSG_MHC\n"
        "Y\t50\t150\tENSG_Y\n",
        encoding="utf-8",
    )

    result = run_gcta_gene_direct(args)

    assert result.metrics["analyzable_variants"] == 1
    assert result.metrics["excluded_input_mhc_variants"] == 1
    assert result.metrics["excluded_input_chromosome_variants"] == 1
    assert result.metrics["excluded_mhc_genes"] == 1
    assert result.metrics["excluded_chromosome_genes"] == 1
    assert result.artifacts["variant_exclusion_list"].path.read_text(
        encoding="utf-8"
    ).splitlines() == ["rs_mhc", "rs_y"]
    assert result.artifacts["scoped_gene_list"].path.read_text(
        encoding="utf-8"
    ) == "1\t50\t150\tENSG000001\n"
    assert result.artifacts["excluded_gene_list"].path.read_text(
        encoding="utf-8"
    ).splitlines() == [
        "6\t28450000\t28470000\tENSG_MHC",
        "Y\t50\t150\tENSG_Y",
    ]
    completion = yaml.safe_load(
        result.artifacts["completion_manifest"].path.read_text(encoding="utf-8")
    )
    assert completion["schema_version"] == 3
    assert completion["outputs"]["html_report"]["path"] == str(
        result.artifacts["html_report"].path.resolve()
    )
    assert completion["inputs"]["source_annotation"]["path"].endswith(
        "/genes.txt"
    )
    assert completion["inputs"]["analysis_annotation"]["path"] == str(
        result.artifacts["scoped_gene_list"].path.resolve()
    )
    assert completion["inputs"]["variant_exclusion_list"]["path"] == str(
        result.artifacts["variant_exclusion_list"].path.resolve()
    )
    summary = result.metrics["scientific_summary"]
    assert summary["resolved_variants"] == 1
    assert summary["excluded_scope_variants"] == 2
    assert summary["excluded_scope_genes"] == 2
    assert summary["expected_chromosomes"] == ["1"]
    assert summary["missing_result_chromosomes"] == []


@pytest.mark.parametrize("method", ["fastbat_segment", "fastbat_set"])
def test_gene_only_mhc_policy_rejects_inputs_without_gene_identities(
    tmp_path,
    method,
):
    args = _args(tmp_path, method)
    config_path = Path(args.run_config)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["mhc"] = {"policy": "exclude_genes"}
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(
        GctaGeneError,
        match="requires a coordinate-defined gene input",
    ):
        run_gcta_gene_direct(args)


def test_gcta_summary_uses_one_value_column_at_every_indent(tmp_path, capsys):
    run_gcta_gene_direct(_args(tmp_path))

    terminal = capsys.readouterr().out
    report = terminal.split(
        "GCTA fastBAT gene-association summary", 1,
    )[1]
    field_lines = [line for line in report.splitlines() if " :" in line]
    assert field_lines
    assert len({
        cell_len(line.split(" :", 1)[0]) for line in field_lines
    }) == 1


def test_scientific_summary_ranks_results_and_counts_bonferroni_hits(tmp_path):
    raw_result = tmp_path / "raw.tsv"
    raw_result.write_text(
        "Gene\tChr\tStart\tEnd\tNo.SNPs\tPvalue\n"
        "GENE_C\t1\t300\t350\t3\t0.2\n"
        "GENE_A\t1\t100\t150\t2\t0.001\n"
        "GENE_B\t1\t200\t250\t4\t0.01\n",
        encoding="utf-8",
    )
    module = load_module_configuration("gcta_gene")
    module.reporting.top_result_count = 2
    result = tmp_path / "results.tsv"
    normalize_gcta_results(raw_result, result, module)

    summary = build_gcta_scientific_summary(
        result,
        module,
        {
            "tested_units": 3,
            "variants": 100,
            "overlapping_variants": 90,
            "unresolved_variants_removed": 10,
            "genes": 3,
            "shared_chromosomes": ["1", "2"],
        },
    )

    assert [item["identifier"] for item in summary["top_associations"]] == [
        "GENE_A", "GENE_B",
    ]
    assert summary["bonferroni_threshold"] == pytest.approx(0.05 / 3)
    assert summary["nominal_significant"] == 2
    assert summary["fdr_bh_significant"] == 2
    assert summary["bonferroni_significant"] == 2
    assert summary["overlap_fraction"] == pytest.approx(0.9)
    assert summary["source_units"] == 3
    assert summary["result_chromosomes"] == ["1"]
    assert summary["missing_result_chromosomes"] == ["2"]
    assert "10 GWAS variants" in summary["warnings"][0]
    assert "shared chromosome(s): 2" in summary["warnings"][1]
    corrected = pl.read_csv(result, separator="\t")
    assert corrected["P_bonferroni"].to_list() == pytest.approx([
        0.6, 0.003, 0.03,
    ])
    assert corrected["P_fdr_bh"].to_list() == pytest.approx([
        0.2, 0.003, 0.015,
    ])


def test_completed_gcta_result_resumes_only_after_manifest_validation(
    tmp_path, monkeypatch, capsys,
):
    args = _args(tmp_path)
    first = run_gcta_gene_direct(args)
    capsys.readouterr()
    args.resume = True
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.run_gcta_command",
        lambda *_args, **_kwargs: pytest.fail("validated result was rerun"),
    )

    resumed = run_gcta_gene_direct(args)

    assert resumed.metrics["tested_units"] == first.metrics["tested_units"]
    assert resumed.metrics["resumed"] is True
    assert (
        "GCTA outputs validated; continuing from completed step"
        in capsys.readouterr().out
    )
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_gene_gcta_gene.log"
    assert "reason=validated_resume" in log.read_text(encoding="utf-8")


def test_reporting_changes_do_not_invalidate_completed_analysis(
    tmp_path, monkeypatch, capsys,
):
    args = _args(tmp_path)
    run_gcta_gene_direct(args)
    capsys.readouterr()
    args.resume = True
    args.gcta_top_results = 1
    args.gcta_nominal_alpha = 0.005
    args.gcta_reporting_alpha = 0.005
    args.gcta_fdr_alpha = 0.005
    args.gcta_p_value_digits = 6
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.run_gcta_command",
        lambda *_args, **_kwargs: pytest.fail("reporting change reran GCTA"),
    )

    resumed = run_gcta_gene_direct(args)

    assert resumed.metrics["resumed"] is True
    summary = resumed.metrics["scientific_summary"]
    assert summary["nominal_alpha"] == 0.005
    assert summary["familywise_alpha"] == 0.005
    assert summary["fdr_alpha"] == 0.005
    assert summary["nominal_significant"] == 0
    assert summary["bonferroni_significant"] == 0
    assert summary["fdr_bh_significant"] == 0
    corrected = pl.read_csv(
        resumed.artifacts["normalized_results"].path,
        separator="\t",
    )
    assert corrected["significant_nominal"].to_list() == [False]
    assert corrected["significant_bonferroni"].to_list() == [False]
    assert corrected["significant_fdr_bh"].to_list() == [False]
    terminal = capsys.readouterr().out
    assert "GCTA outputs validated; continuing from completed step" in terminal
    assert "GCTA fastBAT gene-association summary" in terminal
    assert "Nominally significant associations" in terminal
    assert "FDR-significant associations" in terminal
    assert "Bonferroni threshold" in terminal
    assert "Bonferroni-significant associations" in terminal
    assert "Raw p ≤" in terminal
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_gene_gcta_gene.log"
    assert "gcta_multiple_testing_corrections" in log.read_text(encoding="utf-8")


def test_html_presentation_change_rebuilds_report_without_rerunning_gcta(
    tmp_path, monkeypatch,
):
    args = _args(tmp_path)
    first = run_gcta_gene_direct(args)
    report = first.artifacts["html_report"].path
    assert 'data-page-size="50"' in report.read_text(encoding="utf-8")
    config_path = Path(args.run_config)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["html_report"] = {"page_size": 1}
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    args.resume = True
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.run_gcta_command",
        lambda *_args, **_kwargs: pytest.fail("HTML change reran GCTA"),
    )

    resumed = run_gcta_gene_direct(args)

    assert resumed.metrics["resumed"] is True
    assert 'data-page-size="1"' in report.read_text(encoding="utf-8")
    manifest = yaml.safe_load(
        resumed.artifacts["completion_manifest"].path.read_text(encoding="utf-8")
    )
    assert manifest["schema_version"] == 3
    assert manifest["outputs"]["html_report"]["path"] == str(report.resolve())


def test_legacy_completion_manifest_adds_html_without_rerunning_gcta(
    tmp_path, monkeypatch,
):
    args = _args(tmp_path)
    first = run_gcta_gene_direct(args)
    manifest_path = first.artifacts["completion_manifest"].path
    report_path = first.artifacts["html_report"].path
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = 2
    manifest.pop("report_configuration_sha256")
    manifest["outputs"].pop("html_report")
    manifest_path.write_text(yaml.safe_dump(manifest), encoding="utf-8")
    report_path.unlink()
    args.resume = True
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.run_gcta_command",
        lambda *_args, **_kwargs: pytest.fail("legacy report upgrade reran GCTA"),
    )

    resumed = run_gcta_gene_direct(args)

    assert resumed.metrics["resumed"] is True
    assert report_path.is_file()
    upgraded = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    assert upgraded["schema_version"] == 3
    assert upgraded["outputs"]["html_report"]["path"] == str(
        report_path.resolve()
    )


def test_legacy_completion_manifest_refuses_untracked_html_output(tmp_path):
    args = _args(tmp_path)
    first = run_gcta_gene_direct(args)
    manifest_path = first.artifacts["completion_manifest"].path
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = 2
    manifest.pop("report_configuration_sha256")
    manifest["outputs"].pop("html_report")
    manifest_path.write_text(yaml.safe_dump(manifest), encoding="utf-8")
    args.resume = True

    with pytest.raises(
        GctaGeneError,
        match="legacy GCTA completion manifest does not track",
    ):
        run_gcta_gene_direct(args)


def test_changed_html_path_refuses_untracked_destination(tmp_path):
    args = _args(tmp_path)
    run_gcta_gene_direct(args)
    config_path = Path(args.run_config)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["output_layout"] = {
        "html_report": "reports/{dataset_id}_{method}_replacement.html",
    }
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    replacement = (
        tmp_path / "out/reports/STUDY_fastbat_gene_replacement.html"
    )
    replacement.parent.mkdir(parents=True, exist_ok=True)
    replacement.write_text("user-owned report\n", encoding="utf-8")
    args.resume = True

    with pytest.raises(
        GctaGeneError,
        match="newly configured GCTA HTML-report path already exists",
    ):
        run_gcta_gene_direct(args)
    assert replacement.read_text(encoding="utf-8") == "user-owned report\n"


def test_completed_gcta_result_rejects_changed_html_report_on_resume(tmp_path):
    args = _args(tmp_path)
    result = run_gcta_gene_direct(args)
    result.artifacts["html_report"].path.write_text(
        "changed\n", encoding="utf-8",
    )
    args.resume = True

    with pytest.raises(GctaGeneError, match="HTML report changed"):
        run_gcta_gene_direct(args)


def test_failed_report_refresh_restores_completed_resume_outputs(
    tmp_path, monkeypatch,
):
    args = _args(tmp_path)
    first = run_gcta_gene_direct(args)
    normalized = first.artifacts["normalized_results"].path
    report = first.artifacts["html_report"].path
    manifest = first.artifacts["completion_manifest"].path
    original = {
        "normalized": normalized.read_bytes(),
        "report": report.read_bytes(),
        "manifest": manifest.read_bytes(),
    }
    args.resume = True
    args.gcta_nominal_alpha = 0.005
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.write_gcta_html_report",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            GctaGeneError("synthetic HTML failure")
        ),
    )

    with pytest.raises(GctaGeneError, match="synthetic HTML failure"):
        run_gcta_gene_direct(args)

    assert normalized.read_bytes() == original["normalized"]
    assert report.read_bytes() == original["report"]
    assert manifest.read_bytes() == original["manifest"]
    assert not list(report.parent.glob(".*.pre_refresh.*"))
    assert not list(normalized.parent.glob(".*.pre_correction.*"))


def test_completed_gcta_result_rejects_changed_output_on_resume(tmp_path):
    args = _args(tmp_path)
    result = run_gcta_gene_direct(args)
    result.artifacts["raw_results"].path.write_text(
        "changed\n", encoding="utf-8",
    )
    args.resume = True

    with pytest.raises(GctaGeneError, match="completion outputs changed"):
        run_gcta_gene_direct(args)


def test_pipeline_rejects_formatter_input_that_is_not_bim_compatible(tmp_path):
    args = _args(tmp_path, "fastbat_set", variant="rs1")
    Path(str(tmp_path / "reference") + ".bim").write_text(
        "1\t1_100_A_G\t0\t100\tG\tA\n", encoding="utf-8",
    )
    (tmp_path / "sets.txt").write_text(
        "SET_1\n1_100_A_G\nEND\n", encoding="utf-8",
    )
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"VCF")
    args.vcf = str(vcf)
    context = _pipeline_context(args)

    with pytest.raises(
        GctaGeneError,
        match="No GCTA input variant identifiers occur in the PLINK BIM",
    ):
        run_gcta_gene_direct(args, context)


def test_pipeline_vcf_reproduces_formatter_unique_id_when_vcf_has_rsid(tmp_path):
    unique_id = "1_100_A_G"
    args = _args(tmp_path, "fastbat_set", variant=unique_id)
    Path(str(tmp_path / "reference") + ".bim").write_text(
        "1\t%s\t0\t100\tG\tA\n" % unique_id, encoding="utf-8",
    )
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"VCF")
    args.vcf = str(vcf)
    input_path = Path(args.gcta_input_file)
    input_bytes = input_path.read_bytes()
    context = _pipeline_context(args, variant_id_type="unique")

    result = run_gcta_gene_direct(args, context)

    assert result.metrics["formatter_variant_id_type"] == "unique"
    assert result.metrics["exact_id_matches"] == 1
    assert result.metrics["variant_ids_replaced"] == 0
    assert result.metrics["formatter_input_unmodified"] is True
    assert input_path.read_bytes() == input_bytes
    assert "harmonised_input" not in result.artifacts
    log = tmp_path / "out/logs/STUDY_fastbat_set_gcta_gene.log"
    assert "formatter_variant_id_type=unique" in log.read_text(encoding="utf-8")


def test_pipeline_vcf_requires_formatter_identifier_provenance(tmp_path):
    args = _args(tmp_path)
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"VCF")
    args.vcf = str(vcf)
    args.bcftools = str(_write_fake_bcftools(tmp_path))
    context = _pipeline_context(args)
    del context["formatter"]["gcta_gene"]["variant_id_type"]

    with pytest.raises(
        GctaGeneError,
        match="formatter output is missing variant_id_type provenance",
    ):
        run_gcta_gene_direct(args, context)


def test_pipeline_exact_bim_matching_enforces_configured_minimum_overlap(tmp_path):
    args = _args(tmp_path)
    input_path = Path(args.gcta_input_file)
    input_path.write_text(
        input_path.read_text(encoding="utf-8")
        + "rs2\tC\tT\t0.3\t0.2\t0.1\t0.02\t1000\n",
        encoding="utf-8",
    )
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"VCF")
    args.vcf = str(vcf)
    args.bcftools = str(_write_fake_bcftools(
        tmp_path,
        rows=(("1", 100, "rs1", "A", "G"), ("1", 200, "rs2", "C", "T")),
    ))
    config_path = Path(args.run_config)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["variant_harmonisation"] = {"minimum_overlap_fraction": 0.75}
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    context = _pipeline_context(args)

    with pytest.raises(GctaGeneError, match="Only 1/2 unique GCTA input variants"):
        run_gcta_gene_direct(args, context)


def test_pipeline_rejects_formatter_input_with_bim_allele_conflict(tmp_path):
    args = _args(tmp_path)
    Path(str(tmp_path / "reference") + ".bim").write_text(
        "1\trs1\t0\t100\tC\tA\n", encoding="utf-8",
    )
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"VCF")
    args.vcf = str(vcf)
    context = _pipeline_context(args)

    with pytest.raises(GctaGeneError, match="1 overlapping.*cannot be matched"):
        run_gcta_gene_direct(args, context)


def test_direct_input_reports_partial_exact_id_overlap_without_rewriting(
    tmp_path, capsys,
):
    args = _args(tmp_path)
    input_path = Path(args.gcta_input_file)
    input_path.write_text(
        input_path.read_text(encoding="utf-8")
        + "rs_missing\tC\tT\t0.3\t0.2\t0.1\t0.02\t1000\n",
        encoding="utf-8",
    )
    config_path = Path(args.run_config)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["variant_harmonisation"] = {"minimum_overlap_fraction": 0.99}
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    original_input = input_path.read_bytes()

    result = run_gcta_gene_direct(args)

    assert input_path.read_bytes() == original_input
    assert result.metrics["reference_variants"] == 1
    assert result.metrics["overlapping_variants"] == 1
    assert result.metrics["input_variants_absent_from_reference"] == 1
    assert result.metrics["reference_variants_absent_from_input"] == 0
    assert result.metrics["variant_ids_replaced"] == 0
    assert result.metrics["unresolved_variants_removed"] == 0
    assert result.metrics["direct_input_unmodified"] is True
    assert "harmonised_input" not in result.artifacts
    assert not (tmp_path / "out" / "prepared_inputs").exists()
    log = tmp_path / "out/logs/STUDY_fastbat_gene_gcta_gene.log"
    log_text = log.read_text(encoding="utf-8")
    normalized_log = " ".join(log_text.split())
    assert "input_ids_absent_from_reference=1" in log_text
    assert "input_rewritten=false" in log_text
    assert "direct-input GWAS variant is absent" in normalized_log
    assert "scientific_warnings=1" in normalized_log
    for evidence in (
        "input_format=gcta_ma",
        "required_columns_validated=true",
        "scientific_ranges_validated=true",
        "unique_variant_identifiers_validated=true",
        "allele_pairs_validated=true",
        "declared_genome_build=GRCh37",
        "declared_population=EUR",
        "compatible_allele_pairs=1",
        "bim_structure_validated=true",
        "coordinate_structure_validated=true",
    ):
        assert evidence in log_text
    terminal = " ".join(capsys.readouterr().out.split())
    for label in (
        "Supplied GCTA summary-statistics input",
        "GCTA .ma file",
        "Input format",
        "Required columns",
        "Scientific numeric ranges",
        "Unique non-empty variant identifiers",
        "Non-empty distinct allele pairs",
        "PLINK LD reference",
        "Required companion files",
        "Summary-statistic unique IDs",
        "PLINK BIM unique IDs",
        "Exact IDs shared",
        "Summary IDs absent from BIM",
        "BIM IDs absent from summary",
        "Compatible shared allele pairs",
        "Declared genome build",
        "Declared population",
        "Build and population provenance",
        "GCTA gene-coordinate reference",
        "Coordinate columns",
        "Unique gene identifiers",
        "Positive ordered coordinates",
    ):
        assert label in terminal
    assert "1; will not be used by GCTA" in terminal
    assert "COMPLETED WITH SCIENTIFIC WARNINGS" in terminal


def test_direct_service_rejects_vcf_input(tmp_path):
    args = _args(tmp_path)
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"VCF")
    args.vcf = str(vcf)

    with pytest.raises(
        GctaGeneError,
        match="Direct gcta_gene does not accept --vcf",
    ):
        run_gcta_gene_direct(args)


def test_old_gcta_version_fails_and_is_logged(tmp_path):
    args = _args(tmp_path, version="1.93.3")
    with pytest.raises(GctaGeneError, match="too old"):
        run_gcta_gene_direct(args)
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_gene_gcta_gene.log"
    assert "FAILED" in log.read_text(encoding="utf-8")


def test_valid_version_banner_is_accepted_when_probe_exits_nonzero(tmp_path):
    args = _args(tmp_path)
    args.gcta = str(_write_fake_gcta(tmp_path, version_exit_code=1))

    result = run_gcta_gene_direct(args)

    assert result.metrics["gcta_version"] == "1.94.1"
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_gene_gcta_gene.log"
    assert "probe_exit_code=1" in log.read_text(encoding="utf-8")


def test_zero_reference_overlap_fails_before_analysis(tmp_path):
    args = _args(tmp_path, variant="rs_missing")
    with pytest.raises(
        GctaGeneError,
        match="No direct GCTA summary-statistic SNP IDs occur exactly",
    ):
        run_gcta_gene_direct(args)
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_gene_gcta_gene.log"
    assert "--fastBAT" not in log.read_text(encoding="utf-8")


def test_mbat_combo_rejects_any_incompatible_reference_allele_pair(tmp_path):
    args = _args(tmp_path, "mbat_combo")
    input_path = Path(args.gcta_input_file)
    input_path.write_text(
        input_path.read_text(encoding="utf-8")
        + "rs2\tC\tA\t0.3\t0.2\t0.1\t0.02\t1000\n",
        encoding="utf-8",
    )
    reference = Path(str(tmp_path / "reference") + ".bim")
    reference.write_text(
        reference.read_text(encoding="utf-8")
        + "1\trs2\t0\t110\tG\tA\n",
        encoding="utf-8",
    )

    with pytest.raises(GctaGeneError, match="1 overlapping.*cannot be matched"):
        run_gcta_gene_direct(args)


def test_duplicate_nested_genome_build_declarations_are_rejected(tmp_path):
    config = tmp_path / "invalid.yaml"
    config.write_text(
        "genome_build: GRCh37\nreference:\n  genome_build: GRCh37\n",
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="Extra inputs are not permitted"):
        load_module_configuration("gcta_gene", config)


@pytest.mark.parametrize(
    "text",
    (
        "fastbat_ld_cutoff: 1.1\n",
        "mbat_svd_gamma: 0\n",
        "frequency_difference_max: .nan\n",
        "frequency_difference_max: 0\n",
        "frequency_difference_max: -0.1\n",
        "frequency_difference_max: 1.1\n",
        "variant_harmonisation:\n  minimum_overlap_fraction: 0\n",
        "set_annotation:\n  conversion:\n    minimum_gene_id_overlap_fraction: 0\n",
        "reporting:\n  top_result_count: 0\n",
        "reporting:\n  nominal_alpha: 0\n",
        "reporting:\n  familywise_alpha: 1.1\n",
        "reporting:\n  fdr_alpha: .nan\n",
        "reporting:\n  p_value_significant_digits: 11\n",
        "reporting:\n  correction_columns:\n    fdr_bh_adjusted_p: P_bonferroni\n",
        "reporting:\n  chromosome_columns:\n    fastbat_gene: Missing\n",
    ),
)
def test_scientific_ranges_are_schema_validated(tmp_path, text):
    config = tmp_path / "invalid.yaml"
    config.write_text(text, encoding="utf-8")
    with pytest.raises(ConfigurationError):
        load_module_configuration("gcta_gene", config)


def test_gcta_html_report_configuration_is_schema_validated():
    module = load_module_configuration("gcta_gene")

    assert module.html_report.page_size == 50
    assert module.results.schemas["mbat_combo"].p_value_column in (
        module.html_report.columns["mbat_combo"]
    )
    with pytest.raises(ConfigurationError, match="page_size"):
        load_module_configuration(
            "gcta_gene", cli_overrides={"html_report.page_size": 0},
        )
    with pytest.raises(ConfigurationError, match="html_report"):
        load_module_configuration(
            "gcta_gene",
            cli_overrides={
                "output_layout.html_report": (
                    "reports/{dataset_id}_{method}_gcta_gene_report.txt"
                ),
            },
        )


def test_reference_and_gene_column_orders_are_configuration_driven(tmp_path):
    args = _args(tmp_path, "mbat_combo")
    reference = Path(str(tmp_path / "reference") + ".bim")
    reference.write_text("rs1\t1\tA\tG\t100\t0\n", encoding="utf-8")
    genes = tmp_path / "genes.txt"
    genes.write_text("ENSG000001\t150\t1\t50\n", encoding="utf-8")
    config = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    config["reference"]["bim_columns"] = [
        "variant_id", "chromosome", "allele2", "allele1", "position",
        "genetic_distance",
    ]
    config["gene_annotation"]["columns"] = ["gene", "end", "chromosome", "start"]
    Path(args.run_config).write_text(yaml.safe_dump(config), encoding="utf-8")

    result = run_gcta_gene_direct(args)

    assert result.metrics["compatible_allele_pairs"] == 1
    assert result.metrics["genes"] == 1


@pytest.mark.parametrize(
    ("method", "result_name", "unit_label"),
    (
        ("fastbat_segment", "STUDY.seg.fastbat", "segments"),
        ("fastbat_set", "STUDY.fastbat", "sets"),
    ),
)
def test_all_fastbat_analysis_modes_run_and_normalize(
    tmp_path, method, result_name, unit_label, capsys,
):
    result = run_gcta_gene_direct(_args(tmp_path, method))

    assert result.artifacts["raw_results"].path.name == result_name
    assert result.metrics["tested_units"] == 1
    assert result.metrics["unit_label"] == unit_label
    assert result.artifacts["normalized_results"].path.is_file()
    assert result.artifacts["html_report"].path.is_file()
    report_text = result.artifacts["html_report"].path.read_text(encoding="utf-8")
    assert "All 1 validated %s are included" % unit_label in report_text
    assert "significant_nominal" in report_text
    assert "P_bonferroni" in report_text
    assert "P_fdr_bh" in report_text
    terminal = capsys.readouterr().out
    assert "Top associated %s" % unit_label in terminal
    if method == "fastbat_set":
        for label in (
            "Prepared fastBAT set-list input",
            "Input format",
            "Unique set identifiers",
            "Within-set variant uniqueness",
            "END terminators and non-empty sets",
            "GWAS/BIM set-membership compatibility",
            "Unavailable variant memberships",
            "Final fastBAT set list",
        ):
            assert label in terminal


def test_fastbat_set_rejects_an_unterminated_set_before_execution(
    tmp_path, capsys,
):
    args = _args(tmp_path, "fastbat_set")
    (tmp_path / "sets.txt").write_text("SET_1\nrs1\n", encoding="utf-8")

    with pytest.raises(GctaGeneError, match="missing its terminating END"):
        run_gcta_gene_direct(args)

    log = tmp_path / "out" / "logs" / "STUDY_fastbat_set_gcta_gene.log"
    assert "--fastBAT" not in log.read_text(encoding="utf-8")
    terminal = capsys.readouterr().out
    assert (
        "Failed 3/7 · Validate and prepare the fastBAT set-list input"
        in terminal
    )
    assert "All 7 stages completed" not in terminal


def test_fastbat_set_converts_gmt_with_exact_bim_identifiers(tmp_path, capsys):
    args = _args(tmp_path, "fastbat_set", variant="custom_1_100_A_G")
    Path(str(tmp_path / "reference") + ".bim").write_text(
        "1\tcustom_1_100_A_G\t0\t100\tG\tA\n"
        "1\treference_only\t0\t110\tG\tA\n",
        encoding="utf-8",
    )
    genes = tmp_path / "genes.txt"
    genes.write_text("1\t50\t150\tGENE1\n", encoding="utf-8")
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("SET_1\tdescription\tGENE1\n", encoding="utf-8")
    config = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    config["set_annotation"] = {"file": None, "gmt_file": str(gmt)}
    Path(args.run_config).write_text(yaml.safe_dump(config), encoding="utf-8")

    result = run_gcta_gene_direct(args)

    prepared = tmp_path / "out" / "prepared_sets" / "STUDY"
    set_text = (prepared / "pathways.fastbat.set").read_text(encoding="utf-8")
    assert set_text.split() == ["SET_1", "custom_1_100_A_G", "END"]
    assert "reference_only" not in set_text
    assert result.artifacts["prepared_set_resource"].path == prepared
    assert result.artifacts["analysis_set_list"].path == (
        prepared / "pathways.fastbat.set"
    )
    assert result.metrics["matched_set_variants"] == 1
    summary = result.metrics["scientific_summary"]
    assert summary["gmt_gene_mappability_fraction"] == 1
    assert summary["coordinate_reference_coverage_fraction"] == 1
    assert summary["pathways_with_complete_gene_id_mapping"] == 1
    manifest = yaml.safe_load(
        (prepared / "resource_manifest.yaml").read_text(encoding="utf-8")
    )
    assert manifest["resource"]["genome_build"] == "GRCh37"
    assert manifest["policies"]["variant_identifier_policy"].startswith("copy PLINK BIM")
    terminal = capsys.readouterr().out
    normalized_terminal = " ".join(terminal.split())
    for title in (
        "Validate the GCTA summary-statistics input",
        "Validate the PLINK LD reference and compare it with the GWAS input",
        "Validate GMT and gene resources and prepare fastBAT sets",
        "Run GCTA fastbat_set",
        "Validate the raw GCTA results",
        "Add and validate multiple-testing corrections",
        "Publish validated GCTA outputs and build the scientific summary",
    ):
        assert title in terminal
    for label in (
        "Original pathway GMT",
        "Detected format",
        "Pathway records",
        "GCTA gene-coordinate reference",
        "Unique gene identifiers",
        "Pathway/gene-coordinate compatibility",
        "GMT gene mappability",
        "GMT gene mappability criterion",
        "Coordinate-reference genes absent from GMT",
        "Coordinate-reference coverage by GMT",
        "Required minimum coordinate-reference coverage",
        "Pathways with complete gene-ID mapping",
        "Pathways with partial gene-ID mapping",
        "Pathways without gene-ID mapping",
        "Final fastBAT set resource",
        "Final sets written",
        "Published resource validation",
    ):
        assert label in normalized_terminal
    service_log = (
        tmp_path / "out/logs/STUDY_fastbat_set_gcta_gene.log"
    ).read_text(encoding="utf-8")
    for evidence in (
        "pathway_structure_validated=true",
        "gene_coordinate_structure_validated=true",
        "gene_id_compatibility_validated=true",
        "published_set_resource_validated=true",
    ):
        assert evidence in service_log
    html_report = result.artifacts["html_report"].path.read_text(
        encoding="utf-8"
    )
    for label in (
        "GMT gene mappability",
        "GMT gene mappability criterion",
        "Coordinate-reference coverage by GMT",
        "Required minimum coordinate-reference coverage",
        "Pathways with complete gene-ID mapping",
    ):
        assert label in html_report
    for title in (
        "Validate inputs and scan GMT pathways",
        "Resolve genes and apply configured boundaries",
        "Map analyzable BIM variants to genes",
        "Validate pathway-to-variant memberships",
        "Preflight pathway sizes and disk requirements",
        "Write fastBAT sets and normalized audit tables",
        "Write provenance and checksums",
        "Publish validated pathway resource",
    ):
        assert title in terminal
    assert "All 8 stages completed" in terminal
    assert "All 7 stages completed" in terminal


def test_gmt_bim_mapping_has_separate_measured_variant_progress(
    tmp_path, monkeypatch, capsys,
):
    args = _args(tmp_path, "fastbat_set")
    _configure_three_variant_gmt(args, tmp_path)
    clock = [0.0]

    def advancing_clock():
        clock[0] += 0.02
        return clock[0]

    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.pathway_sets.time.monotonic",
        advancing_clock,
    )

    result = run_gcta_gene_direct(args)

    assert result.metrics["tested_units"] == 1
    terminal = capsys.readouterr().out
    assert "BIM-to-pathway mapping progress" in terminal
    assert (
        "Progress 1/3 · Map analyzable BIM variants to genes · 33%"
        in terminal
    )
    assert (
        "Progress 2/3 · Map analyzable BIM variants to genes · 66%"
        in terminal
    )
    assert "Completed 3/3 · Validate BIM-to-pathway mapping" in terminal
    assert terminal.index("Progress 1/3") < terminal.index("Progress 2/3")
    assert terminal.index("Progress 2/3") < terminal.index("Completed 3/3")


def test_gmt_bim_mapping_failure_preserves_last_measured_variant(
    tmp_path, monkeypatch, capsys,
):
    args = _args(tmp_path, "fastbat_set")
    _configure_three_variant_gmt(args, tmp_path)
    clock = [0.0]

    def advancing_clock():
        clock[0] += 0.02
        return clock[0]

    def fail_mapping(*_args, progress_callback=None, **_kwargs):
        progress_callback(1)
        raise ResourcePreparationError("synthetic BIM mapping failure")

    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.pathway_sets.time.monotonic",
        advancing_clock,
    )
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.pathway_sets._map_bim_variants",
        fail_mapping,
    )

    with pytest.raises(GctaGeneError, match="synthetic BIM mapping failure"):
        run_gcta_gene_direct(args)

    terminal = capsys.readouterr().out
    assert "Failed 1/3 · Map analyzable BIM variants to genes" in terminal
    assert "Completed 3/3 · Validate BIM-to-pathway mapping" not in terminal
    assert "Failed 3/8 · Map analyzable BIM variants to genes" in terminal


def test_fastbat_set_omits_zero_overlap_sets_from_analysis_input(tmp_path, capsys):
    args = _args(tmp_path, "fastbat_set")
    (tmp_path / "sets.txt").write_text(
        "EMPTY_SET\nrs_missing\nEND\n\nSET_1\nrs1\nEND\n",
        encoding="utf-8",
    )

    result = run_gcta_gene_direct(args)

    analysis_set = result.artifacts["analysis_set_list"].path
    assert analysis_set.read_text(encoding="utf-8").split() == [
        "SET_1", "rs1", "END",
    ]
    assert result.metrics["input_sets"] == 2
    assert result.metrics["sets"] == 1
    assert result.metrics["omitted_empty_sets"] == 1
    assert result.metrics["omitted_empty_set_examples"] == ["EMPTY_SET"]
    assert len(result.metrics["scientific_summary"]["warnings"]) == 1
    terminal = capsys.readouterr().out
    assert "COMPLETED WITH SCIENTIFIC WARNINGS" in terminal
    assert "1 custom sets had no variants shared" in terminal
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_set_gcta_gene.log"
    log_text = log.read_text(encoding="utf-8")
    assert "fastbat_set_intersection" in log_text
    assert "1 custom sets had no variants shared" in log_text


def test_fastbat_set_empty_overlap_error_policy_stops_before_gcta(tmp_path):
    args = _args(tmp_path, "fastbat_set")
    (tmp_path / "sets.txt").write_text(
        "EMPTY_SET\nrs_missing\nEND\n\nSET_1\nrs1\nEND\n",
        encoding="utf-8",
    )
    config_path = Path(args.run_config)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["set_annotation"]["conversion"] = {"empty_pathway_policy": "error"}
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(GctaGeneError, match="1 fastBAT sets have no variants"):
        run_gcta_gene_direct(args)

    log = tmp_path / "out" / "logs" / "STUDY_fastbat_set_gcta_gene.log"
    assert "--fastBAT" not in log.read_text(encoding="utf-8")


@pytest.mark.parametrize("policy", ["omit", "error"])
def test_fastbat_set_enforces_gcta_variant_limit_before_execution(tmp_path, policy):
    maximum_set_variants = (
        load_configuration().modules.gcta_gene.set_annotation.maximum_set_variants
    )
    variants = ["rs%d" % index for index in range(maximum_set_variants + 1)]
    source = tmp_path / "sets.txt"
    source.write_text(
        "OVERSIZED\n%s\nEND\n\nVALID\nrs_valid\nEND\n"
        % "\n".join(variants),
        encoding="utf-8",
    )
    identifiers = set(variants) | {"rs_valid"}
    destination = tmp_path / "analysis.set"

    if policy == "error":
        with pytest.raises(GctaGeneError, match="20,000-variant hard limit"):
            _prepare_analysis_set_list(
                source,
                destination,
                identifiers,
                identifiers,
                "omit",
                policy,
                maximum_set_variants,
            )
        return

    _, metrics = _prepare_analysis_set_list(
        source,
        destination,
        identifiers,
        identifiers,
        "omit",
        policy,
        maximum_set_variants,
    )
    assert destination.read_text(encoding="utf-8").split() == [
        "VALID", "rs_valid", "END",
    ]
    assert metrics["sets"] == 1
    assert metrics["omitted_oversized_sets"] == 1
    assert metrics["omitted_oversized_set_examples"] == [
        "OVERSIZED (20001)"
    ]
    assert metrics["analysis_set_variants"] == 1


def test_fastbat_set_rejects_gmt_and_prepared_set_together(tmp_path):
    args = _args(tmp_path, "fastbat_set")
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("SET_1\tdescription\tENSG000001\n", encoding="utf-8")
    config = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    config["set_annotation"]["gmt_file"] = str(gmt)
    Path(args.run_config).write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(ConfigurationError, match="mutually exclusive"):
        run_gcta_gene_direct(args)


def test_prepared_set_input_bypasses_gmt_conversion(tmp_path, monkeypatch):
    args = _args(tmp_path, "fastbat_set")
    gmt = tmp_path / "configured_pathways.gmt"
    gmt.write_text("SET_1\tdescription\tENSG000001\n", encoding="utf-8")
    config = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    config["set_annotation"] = {"file": None, "gmt_file": str(gmt)}
    Path(args.run_config).write_text(yaml.safe_dump(config), encoding="utf-8")
    args.fastbat_set_list = str(tmp_path / "sets.txt")

    def unexpected_conversion(_):
        raise AssertionError("GMT converter must not run for --fastbat-set-list")

    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.prepare_resource",
        unexpected_conversion,
    )
    result = run_gcta_gene_direct(args)

    assert result.metrics["tested_units"] == 1
    assert "prepared_set_resource" not in result.artifacts


def test_gmt_is_rejected_for_non_set_method(tmp_path):
    args = _args(tmp_path, "fastbat_gene")
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("SET_1\tdescription\tENSG000001\n", encoding="utf-8")
    config = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    config["set_annotation"] = {"file": None, "gmt_file": str(gmt)}
    Path(args.run_config).write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(GctaGeneError, match="valid only with --method fastbat_set"):
        run_gcta_gene_direct(args)


def test_dry_run_removes_empty_staging_directory(tmp_path):
    args = _args(tmp_path)
    args.dry_run = True
    result = run_gcta_gene_direct(args)
    assert result.metrics["dry_run"] is True
    partial = tmp_path / "out" / ".partial" / "STUDY_fastbat_gene"
    assert not partial.exists()
