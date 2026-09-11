"""Native, small-data regressions for association plotting and configuration."""

import argparse
import gzip
import math
import os
from pathlib import Path
import shutil
import subprocess

import pandas as pd
import pytest
import yaml

from postgwas.config import load_configuration
from postgwas.modules.manhattan import service
from postgwas.modules.manhattan.cli import build_parser


@pytest.fixture
def native_vcf(tmp_path):
    pysam = pytest.importorskip("pysam")
    if not all(shutil.which(tool) for tool in ("bcftools", "Rscript")):
        pytest.skip("native Manhattan tests require bcftools and Rscript")
    raw = tmp_path / "input.vcf"
    raw.write_text(
        '##fileformat=VCFv4.2\n##genome_build=GRCh37\n'
        '##postgwas_version="test"\n'
        '##postgwas_dataset_id="STUDY"\n'
        '##postgwas_vcf_status="harmonised"\n'
        '##contig=<ID=1,length=249250621>\n##contig=<ID=5,length=180915260>\n'
        '##FORMAT=<ID=LP,Number=A,Type=Float,Description="LP">\n'
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">\n'
        '##INFO=<ID=AS,Number=2,Type=Integer,Description="Allelic counts">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY\n'
        '1\t100\t1_100_A_G\tA\tG\t.\tPASS\tAS=5,5\tLP:AF\t3.25:0.2\n'
        '5\t200\t5_200_C_T\tC\tT\t.\tPASS\tAS=0,10\tLP:AF\t20:0.3\n'
        '5\t300\t5_300_G_A\tG\tA\t.\tPASS\tAS=2,8\tLP:AF\t.:0.4\n'
    )
    vcf = Path(str(raw) + ".gz")
    pysam.tabix_compress(str(raw), str(vcf))
    pysam.tabix_index(str(vcf), preset="vcf")
    return vcf


def plot_args(tmp_path, vcf, **values):
    return argparse.Namespace(
        vcf=str(vcf), dataset_id="STUDY", output_directory=str(tmp_path / "out"),
        threads=1, memory_gb=1, **values,
    )


@pytest.mark.parametrize("role", ("version", "dataset_id", "status"))
@pytest.mark.parametrize("value", (None, '""', '"unterminated', '<ID=invalid>'))
def test_direct_manhattan_rejects_invalid_provenance_before_work(
    monkeypatch, tmp_path, role, value,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"header-query-fixture")
    headers = load_configuration().modules.formatting.input_contract.provenance_headers.model_dump()
    metadata = {name: '"recorded"' for name in headers.values()}
    name = headers[role]
    if value is None:
        del metadata[name]
    else:
        metadata[name] = value
    header = "##fileformat=VCFv4.2\n" + "".join(
        f"##{key}={payload}\n" for key, payload in metadata.items()
    ) + "##genome_build=GRCh37\n##contig=<ID=1,length=249250621>\n#CHROM\n"

    def forbidden(*_args, **_kwargs):
        pytest.fail("Invalid provenance reached index/sample/R/plotting work")

    monkeypatch.setattr(service, "resolve_executable", lambda value, _label: value)
    monkeypatch.setattr(service.subprocess, "run", forbidden)
    monkeypatch.setattr("postgwas.core.vcf.read_vcf_header", lambda *_args, **_kwargs: header)
    monkeypatch.setattr("postgwas.core.vcf.count_indexed_vcf_records", forbidden)
    monkeypatch.setattr("postgwas.core.vcf.select_vcf_sample", forbidden)

    with pytest.raises(RuntimeError, match=name):
        service.run_assoc_plot_direct(plot_args(tmp_path, vcf))

    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("custom_headers", (False, True))
def test_manhattan_accepts_recorded_origin_with_configured_header_names(
    monkeypatch, tmp_path, custom_headers,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"header-query-fixture")
    overrides = {}
    if custom_headers:
        overrides["modules.formatting.input_contract.provenance_headers"] = {
            "version": "origin_version", "dataset_id": "origin_study",
            "status": "origin_status",
        }
    configuration = load_configuration(cli_overrides=overrides)
    headers = configuration.modules.formatting.input_contract.provenance_headers.model_dump()
    values = {"version": "older-release", "dataset_id": "STUDY", "status": "harmonised"}
    header = "##fileformat=VCFv4.2\n" + "".join(
        f'##{headers[role]}="{value}"\n' for role, value in values.items()
    ) + (
        "##genome_build=GRCh37\n##contig=<ID=1,length=249250621>\n"
        '##FORMAT=<ID=LP,Number=A,Type=Float,Description="LP">\n'
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY\n"
    )
    monkeypatch.setattr(service, "load_run_configuration_for_module", lambda *_args, **_kwargs: configuration)
    monkeypatch.setattr(service, "resolve_executable", lambda value, _label: value)
    monkeypatch.setattr("postgwas.core.vcf.read_vcf_header", lambda *_args, **_kwargs: header)
    monkeypatch.setattr("postgwas.core.vcf.count_indexed_vcf_records", lambda *_args, **_kwargs: 1)
    monkeypatch.setattr("postgwas.core.vcf.select_vcf_sample", lambda *_args, **_kwargs: "STUDY")

    def runtime(command, **_kwargs):
        assert command[-1] == "--check-runtime"
        return subprocess.CompletedProcess(command, 0, stdout="test runtime\n")

    monkeypatch.setattr(service.subprocess, "run", runtime)
    _, indexed, _, runtime_text = service.resolve_plot_inputs(plot_args(tmp_path, vcf))

    assert indexed.header == header
    assert indexed.variant_count == 1 and indexed.sample == "STUDY"
    assert indexed.genome_build == "GRCh37"
    assert runtime_text == "test runtime"
    assert not (tmp_path / "out").exists()


def test_native_first_row_offsets_and_loglog(native_vcf, tmp_path):
    args = plot_args(tmp_path, native_vcf, min_lp=3.1, loglog_pval=5.5, max_height=25.5)
    result = service.run_assoc_plot_direct(args)
    points = pd.read_csv(result["plot_data"], sep="\t")
    assert list(points.pos) == [100, 200]  # headerless VCF query must retain row 1
    assert list(points.raw_lp) == [3.25, 20]
    expected_offset = sum([249250621, 243199373, 198022430, 191154276]) + 4 * 20_000_000
    assert list(points.chrompos) == [100, expected_offset + 200]
    assert points.lp.iloc[0] == 3.25
    assert points.lp.iloc[1] == pytest.approx(5.5 * math.log10(20) / math.log10(5.5))
    log = Path(result["log"]).read_text()
    assert "--max-height=25.5" in log and "--min-lp=3.1" in log
    assert "--significance=5e-08" in log and "--suggestive=1e-05" in log
    assert "VALIDATED: 2 plotted records; 1 excluded" in log
    assert "--axis-label-rows=2" in log
    assert "Plotted 2 of 3 records; minimum LP: 3.1; minimum AF filter: 0." in log
    assert "Above LP=5.5, Y-axis spacing is log-log" in log
    assert "Chromosomes with no plotted records: 2, 3, 4, 6" in log
    assert Path(result["plot"]).read_bytes().startswith(b"%PDF-")


@pytest.mark.parametrize("output_format", ["png", "pdf"])
def test_yaml_values_and_explicit_cli_precedence(native_vcf, tmp_path, output_format):
    config = tmp_path / "run.yaml"
    config.write_text(yaml.safe_dump({"modules": {"manhattan": {
        "file_format": "png", "minimum_neglog10_p": 10.5,
        "width": 8.25, "significance_threshold": 1e-9,
    }}}))
    values = {"run_config": str(config)}
    if output_format == "pdf":
        values.update(pdf=str(tmp_path / "explicit.pdf"), min_lp=3.1)
    result = service.run_assoc_plot_direct(plot_args(tmp_path, native_vcf, **values))
    assert Path(result["plot"]).suffix == "." + output_format
    log = Path(result["log"]).read_text()
    assert "--width=8.25" in log and "--significance=1e-09" in log
    assert '"file_format": "%s"' % output_format in log
    points = pd.read_csv(result["plot_data"], sep="\t")
    assert len(points) == (2 if output_format == "pdf" else 1)


def test_allelic_shift_exact_binomial_is_capped(native_vcf, tmp_path):
    result = service.run_assoc_plot_direct(plot_args(tmp_path, native_vcf, allelic_shift=True, min_lp=0))
    points = pd.read_csv(result["plot_data"], sep="\t")
    assert points.raw_lp.iloc[0] == 0
    assert points.raw_lp.iloc[1] == pytest.approx(-math.log10(2 / 1024))


@pytest.mark.parametrize("values,match", [
    ({"genome_build": "GRCh38"}, "does not match"),
    ({"pheno": "WRONG"}, "does not match"),
    ({"pdf": "a.pdf", "png": "a.png"}, "mutually exclusive"),
])
def test_invalid_inputs_fail_before_output(native_vcf, tmp_path, values, match):
    with pytest.raises(Exception, match=match):
        service.run_assoc_plot_direct(plot_args(tmp_path, native_vcf, **values))
    assert not (tmp_path / "out").exists()


def test_no_points_is_failure_not_completion(native_vcf, tmp_path):
    with pytest.raises(RuntimeError, match="failed"):
        service.run_assoc_plot_direct(plot_args(tmp_path, native_vcf, min_lp=100))
    assert not list((tmp_path / "out").glob("*.pdf"))
    assert "FAILED:" in (tmp_path / "out" / "STUDY_manhattan.log").read_text()


def test_successful_external_exit_without_outputs_is_failure(native_vcf, tmp_path, monkeypatch):
    original = service.subprocess.run

    def fake_render(command, **kwargs):
        if command[0].endswith("Rscript") and "--check-runtime" not in command:
            return subprocess.CompletedProcess(command, 0)
        return original(command, **kwargs)

    monkeypatch.setattr(service.subprocess, "run", fake_render)
    with pytest.raises(RuntimeError, match="does not exist or is empty"):
        service.run_assoc_plot_direct(plot_args(tmp_path, native_vcf))


@pytest.mark.parametrize("values", [
    {"loglog_pvalue": 1}, {"loglog_pvalue": 0.5},
    {"minimum_af": 0.6}, {"significance_threshold": 0},
    {"maximum_height": 1}, {"allelic_shift": True, "minimum_af": 0.1},
    {"pvalue_field": "LP;invalid"}, {"png_dpi": 0},
    {"axis_label_rows": 0}, {"caption_font_size": 0},
])
def test_invalid_configuration(values):
    with pytest.raises(Exception):
        load_configuration(cli_overrides={"modules.manhattan." + key: value for key, value in values.items()})


def test_parser_defaults_do_not_override_yaml_and_output_is_exclusive():
    parser = build_parser()
    args = parser.parse_args([])
    for key in (*service._PLOT_OPTIONS, "genome_build", "png", "pdf"):
        assert not hasattr(args, key)
    with pytest.raises(SystemExit):
        parser.parse_args(["--pdf", "a.pdf", "--png", "a.png"])


def test_bundled_assembly_lengths_match_ucsc_reference_invariants():
    import re

    adapter = (Path(service.__file__).parent / "resources" / "assoc_plot.R").read_text()
    # Full primary-chromosome vectors from UCSC hg19/hg38 chrom.sizes, cited in
    # the module documentation; guards the upstream +1000-bp hg19 errors.
    expected = [
        [249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
         159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
         115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
         59128983, 63025520, 48129895, 51304566, 155270560],
        [248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
         159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
         114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
         58617616, 64444167, 46709983, 50818468, 156040895],
    ]
    observed = [[int(value.strip()) for value in match.split(",")]
                for match in re.findall(r"chrlen <- setNames\(c\(([\d, ]+)\)", adapter)]
    assert observed == expected


def test_pipeline_help_does_not_claim_unimplemented_qq_plots():
    from postgwas.pipeline.cli import _add_pipeline_selection_arguments
    from postgwas.pipeline.registry import REGISTRY

    parser = argparse.ArgumentParser()
    _add_pipeline_selection_arguments(parser)
    assert "QQ" not in parser.format_help()
    assert REGISTRY.get("manhattan").description == "Generate a Manhattan plot."


def test_every_direct_example_has_required_inputs_and_parses(tmp_path):
    import re
    import shlex
    from rich.text import Text
    from postgwas.config.cli import build_parser as config_parser

    parser = build_parser()
    input_path = tmp_path / "example.vcf.gz"
    input_path.write_text("parser-only existence fixture\n")
    plain = Text.from_markup(parser.epilog).plain
    commands = [" ".join(block.replace("\\\n", " ").split())
                for block in re.findall(r"  postgwas [\s\S]+?(?=\n\n|$)", plain)]
    assert len(commands) == 4
    parsed_modes = []
    for command in commands:
        tokens = shlex.split(command)
        if tokens[1] == "config":
            config = config_parser().parse_args(tokens[2:])
            assert config.module == "manhattan"
            continue
        tokens[tokens.index("--vcf") + 1] = str(input_path)
        tokens[tokens.index("--output-directory") + 1] = str(tmp_path / "results")
        args = parser.parse_args(tokens[2:])
        assert args.vcf and args.dataset_id == "STUDY" and args.output_directory
        parsed_modes.append((getattr(args, "csq", False), getattr(args, "allelic_shift", False)))
    assert parsed_modes == [(False, False), (True, False), (False, True)]


def test_manhattan_wiki_contract_matches_validated_runtime():
    root = Path(__file__).parents[1]
    text = (root / "docs/wiki/modules/manhattan.md").read_text()
    template = (root / "docs/templates/module-page.md").read_text()
    headings = lambda value: [line for line in value.splitlines() if line.startswith("## ")]
    assert headings(text) == headings(template)
    for obsolete in ("unknown_sample", "first VCF", "PNG wins", "assocplot_<timestamp>"):
        assert obsolete not in text
    for required in ("single-sample", "mutually exclusive", "_manhattan_points.tsv", "log10(LP)"):
        assert required in text


@pytest.mark.parametrize("old,new,match", [
    ("Number=A", "Number=2", "must be numeric"),
    ("1\t100\t1_100_A_G\tA\tG\t", "1\t100\t1_100_A_G\tA\tG,T\t", "failed"),
    ("3.25:0.2", "-1:0.2", "failed"),
])
def test_invalid_record_contract(native_vcf, tmp_path, old, new, match):
    import pysam

    with gzip.open(native_vcf, "rt") as handle:
        content = handle.read().replace(old, new)
    raw = tmp_path / "invalid.vcf"
    raw.write_text(content)
    compressed = Path(str(raw) + ".gz")
    pysam.tabix_compress(str(raw), str(compressed))
    pysam.tabix_index(str(compressed), preset="vcf")
    with pytest.raises(Exception, match=match):
        service.run_assoc_plot_direct(plot_args(tmp_path, compressed, min_lp=0))


def test_pipeline_preflight_resolves_plot_runtime_and_header(native_vcf, tmp_path):
    from postgwas.pipeline.preflight import preflight_manhattan
    from preflight_support import pipeline_input_vcf_evidence

    evidence = preflight_manhattan(
        plot_args(tmp_path, native_vcf),
        preflight_evidence=pipeline_input_vcf_evidence(),
    )
    assert evidence.module == "manhattan"
    assert "data.table" in evidence.resources["runtime"]
    assert not (tmp_path / "out").exists()


def test_pipeline_common_vcf_evidence_precedes_plot_runtime(monkeypatch):
    from postgwas.pipeline.preflight import preflight_manhattan

    def forbidden_runtime(*args, **kwargs):
        raise AssertionError("R runtime must not run before common VCF validation")

    monkeypatch.setattr(service, "resolve_plot_inputs", forbidden_runtime)
    with pytest.raises(RuntimeError, match="shared input-VCF validation"):
        preflight_manhattan(argparse.Namespace(), preflight_evidence=None)


def test_native_coding_highlight_keeps_unannotated_sites(native_vcf, tmp_path):
    import pysam

    bcftools = os.environ.get("POSTGWAS_TEST_SPLIT_VEP_BCFTOOLS")
    if not bcftools:
        pytest.skip("Set POSTGWAS_TEST_SPLIT_VEP_BCFTOOLS to a tested split-vep-enabled executable")
    with gzip.open(native_vcf, "rt") as handle:
        content = handle.read().replace(
            '#CHROM\t', '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. Format: Allele|Consequence">\n#CHROM\t',
        ).replace("AS=5,5\t", "AS=5,5;CSQ=G|missense_variant\t")
    raw = tmp_path / "annotated.vcf"
    raw.write_text(content)
    compressed = Path(str(raw) + ".gz")
    pysam.tabix_compress(str(raw), str(compressed))
    pysam.tabix_index(str(compressed), preset="vcf")
    config = tmp_path / "annotated.yaml"
    config.write_text(yaml.safe_dump({"modules": {"manhattan": {}}, "resources": {"executables": {"bcftools": bcftools}}}))
    result = service.run_assoc_plot_direct(plot_args(tmp_path, compressed, run_config=str(config), csq=True))
    points = pd.read_csv(result["plot_data"], sep="\t")
    assert list(points.pos) == [100, 200]
    assert list(points.coding) == [True, False]


@pytest.mark.parametrize("chromosome", ["chr5", "chrX", "MT"])
def test_configured_chromosome_normalization(native_vcf, tmp_path, chromosome):
    import pysam

    with gzip.open(native_vcf, "rt") as handle:
        content = handle.read().replace("ID=1,", "ID=chr1,").replace("\n1\t", "\nchr1\t")
    length = {"chr5": 180915260, "chrX": 155270560, "MT": 16569}[chromosome]
    content = content.replace("ID=5,length=180915260", "ID=%s,length=%s" % (chromosome, length)).replace("\n5\t", "\n%s\t" % chromosome)
    raw = tmp_path / "labels.vcf"
    raw.write_text(content)
    compressed = Path(str(raw) + ".gz")
    pysam.tabix_compress(str(raw), str(compressed))
    pysam.tabix_index(str(compressed), preset="vcf")
    if chromosome == "MT":
        with pytest.raises(RuntimeError, match="failed"):
            service.run_assoc_plot_direct(plot_args(tmp_path, compressed))
    else:
        result = service.run_assoc_plot_direct(plot_args(tmp_path, compressed))
        points = pd.read_csv(result["plot_data"], sep="\t")
        assert list(points.source_chrom) == ["chr1", chromosome]
        assert list(points.chrom.astype(str)) == ["1", chromosome[3:]]
        assert "labels changed: 2" in Path(result["log"]).read_text()
