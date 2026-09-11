"""Scientific-contract tests for population LD-block annotation."""

import argparse
import csv
import gzip
from io import StringIO
from pathlib import Path
import shutil
import subprocess

import pytest
from rich.cells import cell_len
from rich.console import Console

from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError, MissingRequiredArgumentsError
from postgwas.core.interval_validation import validate_vcf_bed_contigs
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.core.ui import StageProgress
from postgwas.modules.ld_annotation import annot_ldblock
from postgwas.modules.ld_annotation import reporting as annotation_reporting
from postgwas.modules.ld_annotation import service as annotation_service
from postgwas.modules.ld_annotation.reporting import LDReferenceSummary
from postgwas.pipeline.registry import REGISTRY


BED_FILENAME_TEMPLATE = "{genome_build}_{population}_ldetect.bed.gz"
INFO_FIELD_TEMPLATE = "{population}_LDblock"
INFO_DESCRIPTION_TEMPLATE = (
    "LD Block ID from Berisa & Pickrell 2016 ({population}, {genome_build})"
)
OUTPUT_FILENAME_TEMPLATE = "{dataset_id}_ldblock.vcf.gz"
SUMMARY_FILENAME_TEMPLATE = "{dataset_id}_ldblock_summary.csv"
HTML_REPORT_FILENAME_TEMPLATE = "{dataset_id}_ldblock_report.html"
PROVENANCE_HEADERS = (
    load_configuration().modules.formatting.input_contract.provenance_headers.model_dump()
)
STUDY_PROVENANCE = (
    '##postgwas_version="test"\n'
    '##postgwas_dataset_id="STUDY"\n'
    '##postgwas_vcf_status="harmonised"\n'
)


def _write_module_configuration(
    path: Path,
    *,
    vcf: Path | None,
    ld_region_dir: Path | None,
    dataset_id: str = "YAML_STUDY",
    output_directory: Path | None = None,
) -> None:
    output_directory = output_directory or path.parent / "yaml_output"
    path.write_text(
        "inputs:\n"
        "  vcf: %s\n"
        "  ld_region_dir: %s\n"
        "  dataset_id: %s\n"
        "output_directory: %s\n"
        "populations: [AFR]\n"
        % (
            "null" if vcf is None else vcf,
            "null" if ld_region_dir is None else ld_region_dir,
            dataset_id,
            output_directory,
        ),
        encoding="utf-8",
    )


def _fake_annotation_outputs(vcf: Path, population: str) -> dict:
    summary_file = vcf.with_name("STUDY_ldblock_summary.csv")
    html_report = vcf.with_name("STUDY_ldblock_report.html")
    return {
        "annotated_vcf": vcf,
        "annotated_index": Path(str(vcf) + ".tbi"),
        "summary_file": summary_file,
        "html_report": html_report,
        "summary": {
            "dataset_id": "STUDY",
            "genome_build": "GRCh38",
            "total_variants": 1,
            "any_population_annotated": 1,
            "all_populations_annotated": 1,
            "fully_unassigned_variants": 0,
            "unassigned_on_bed_contigs": 0,
            "unassigned_outside_bed_contigs": 0,
            "unassigned_by_contig": {},
            "bed_contigs": ["1"],
            "populations": {
                population: {
                    "annotated_variants": 1,
                    "unassigned_variants": 0,
                    "annotation_percent": 100.0,
                    "bed_blocks": 1,
                    "blocks_used": 1,
                    "empty_blocks": 0,
                    "duplicate_bed_labels": 0,
                    "unexpected_annotation_labels": 0,
                    "bed_file": "reference.bed.gz",
                }
            },
        },
    }


def test_ld_annotation_yaml_inputs_and_cli_precedence(tmp_path):
    yaml_vcf = tmp_path / "yaml.vcf.gz"
    yaml_beds = tmp_path / "yaml_beds"
    yaml_output = tmp_path / "yaml_output"
    run_config = tmp_path / "ld_annotation.yaml"
    _write_module_configuration(
        run_config,
        vcf=yaml_vcf,
        ld_region_dir=yaml_beds,
        output_directory=yaml_output,
    )

    configured = annotation_service.resolve_ld_annotation_configuration(
        argparse.Namespace(run_config=str(run_config))
    )
    module = configured.modules.ld_annotation
    assert module.inputs.vcf == yaml_vcf
    assert module.inputs.ld_region_dir == yaml_beds
    assert module.inputs.dataset_id == "YAML_STUDY"
    assert module.output_directory == yaml_output
    assert [population.value for population in module.populations] == ["AFR"]

    cli_vcf = tmp_path / "cli.vcf.gz"
    cli_beds = tmp_path / "cli_beds"
    cli_output = tmp_path / "cli_output"
    overridden = annotation_service.resolve_ld_annotation_configuration(
        argparse.Namespace(
            run_config=str(run_config),
            vcf=str(cli_vcf),
            ld_region_dir=str(cli_beds),
            dataset_id="CLI_STUDY",
            output_directory=str(cli_output),
            ld_block_populations=["EUR", "EAS"],
            threads=7,
        )
    )
    module = overridden.modules.ld_annotation
    assert module.inputs.vcf == cli_vcf
    assert module.inputs.ld_region_dir == cli_beds
    assert module.inputs.dataset_id == "CLI_STUDY"
    assert module.output_directory == cli_output
    assert [population.value for population in module.populations] == ["EUR", "EAS"]
    assert overridden.execution.threads == 7


def test_ld_annotation_run_uses_yaml_directory_without_cli_argument(
    monkeypatch,
    tmp_path,
    capsys,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    beds = tmp_path / "beds"
    beds.mkdir()
    output = tmp_path / "output"
    run_config = tmp_path / "ld_annotation.yaml"
    _write_module_configuration(
        run_config,
        vcf=vcf,
        ld_region_dir=beds,
        output_directory=output,
    )
    captured = {}

    def annotate(**kwargs):
        captured.update(kwargs)
        return _fake_annotation_outputs(
            output / "YAML_STUDY_ldblock.vcf.gz",
            "AFR",
        )

    monkeypatch.setattr(annotation_service, "annotate_ldblocks", annotate)
    outputs = annotation_service.run_annot_ldblock(
        argparse.Namespace(run_config=str(run_config))
    )

    assert captured["vcf_path"] == str(vcf.resolve())
    assert captured["ld_dir"] == str(beds.resolve())
    assert captured["output_directory"] == str(output.resolve())
    assert captured["dataset_id"] == "YAML_STUDY"
    assert captured["provenance_headers"] == PROVENANCE_HEADERS
    assert [population.value for population in captured["populations"]] == ["AFR"]
    assert captured["html_report_filename_template"] == (
        HTML_REPORT_FILENAME_TEMPLATE
    )
    assert isinstance(captured["stage_progress"], StageProgress)
    assert captured["stage_progress"].enabled is True
    assert outputs["annotated_vcf"].name == "YAML_STUDY_ldblock.vcf.gz"
    assert outputs["log_file"].name == "YAML_STUDY_ldblock.log"
    assert outputs["log_file"].is_file()
    terminal = capsys.readouterr().out
    assert "LD-block annotation summary" in terminal
    assert "AFR LD blocks" in terminal
    assert "Summary CSV" in terminal
    assert "Detailed HTML report" in terminal
    assert "Detailed log" in terminal
    assert ".csv" in terminal
    assert terminal.endswith("\n\n")
    terminal_fields = [
        line for line in terminal.splitlines() if " : " in line
    ]
    assert {
        cell_len(line.split(" : ", 1)[0]) for line in terminal_fields
    } == {captured["stage_progress"].outcome_separator_column}
    assert "ld_annotation_run" in outputs["log_file"].read_text(
        encoding="utf-8"
    )


def test_direct_cli_accepts_yaml_only_inputs(monkeypatch, tmp_path):
    from postgwas.modules.ld_annotation.cli import main

    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    beds = tmp_path / "beds"
    beds.mkdir()
    run_config = tmp_path / "ld_annotation.yaml"
    _write_module_configuration(
        run_config,
        vcf=vcf,
        ld_region_dir=beds,
    )
    called = False

    def annotate(**_kwargs):
        nonlocal called
        called = True
        return _fake_annotation_outputs(tmp_path / "annotated.vcf.gz", "AFR")

    monkeypatch.setattr(annotation_service, "annotate_ldblocks", annotate)
    monkeypatch.setattr(
        "sys.argv",
        ["postgwas annot_ldblock", "--run-config", str(run_config)],
    )
    assert main() == 0
    assert called


def test_direct_cli_reports_canonical_missing_input_paths(monkeypatch, capsys):
    from postgwas.modules.ld_annotation.cli import main

    monkeypatch.setattr(
        "sys.argv",
        ["postgwas annot_ldblock", "--threads", "1"],
    )
    with pytest.raises(SystemExit) as captured:
        main()
    assert captured.value.code == 2
    terminal = capsys.readouterr().err
    assert (
        "Required argument not provided: --vcf. Provide --vcf VALUE or set "
        "modules.ld_annotation.inputs.vcf in the run configuration."
    ) in terminal
    assert (
        "Required argument not provided: --ld-region-dir. Provide "
        "--ld-region-dir VALUE or set modules.ld_annotation.inputs.ld_region_dir "
        "in the run configuration."
    ) in terminal


def test_ld_annotation_reports_all_missing_resolved_inputs_before_execution(
    monkeypatch,
):
    called = False

    def annotate(**_kwargs):
        nonlocal called
        called = True

    monkeypatch.setattr(annotation_service, "annotate_ldblocks", annotate)
    with pytest.raises(MissingRequiredArgumentsError) as captured:
        annotation_service.run_annot_ldblock(argparse.Namespace())

    message = str(captured.value)
    assert (
        "Required argument not provided: --vcf. Provide --vcf VALUE or set "
        "modules.ld_annotation.inputs.vcf in the run configuration."
    ) in message
    assert (
        "Required argument not provided: --ld-region-dir. Provide "
        "--ld-region-dir VALUE or set modules.ld_annotation.inputs.ld_region_dir "
        "in the run configuration."
    ) in message
    assert not called


def test_invalid_yaml_ld_directory_fails_before_annotation(monkeypatch, tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    run_config = tmp_path / "ld_annotation.yaml"
    _write_module_configuration(
        run_config,
        vcf=vcf,
        ld_region_dir=tmp_path / "missing_beds",
    )
    called = False

    def annotate(**_kwargs):
        nonlocal called
        called = True

    monkeypatch.setattr(annotation_service, "annotate_ldblocks", annotate)
    with pytest.raises(ConfigurationError, match="LD-block directory does not exist"):
        annotation_service.run_annot_ldblock(
            argparse.Namespace(run_config=str(run_config))
        )
    assert not called


def test_pipeline_registry_declares_configurable_ld_directory():
    specification = REGISTRY.get("annot_ldblock")
    required = {
        option.dest: (option.flag, option.config_path)
        for option in specification.required_options
    }
    assert required["ld_region_dir"] == (
        "--ld-region-dir",
        "modules.ld_annotation.inputs.ld_region_dir",
    )
    assert specification.preflight == (
        "postgwas.modules.ld_annotation.service:preflight_ld_annotation"
    )


def _write_bed(path: Path, rows: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(rows)


def _stub_external_execution(monkeypatch, header: str, commands: list[list[str]]):
    monkeypatch.setattr(
        annot_ldblock,
        "resolve_executable",
        lambda value, _label: value,
    )
    monkeypatch.setattr(
        annot_ldblock,
        "read_vcf_header",
        lambda *_args, **_kwargs: header,
    )
    monkeypatch.setattr(
        annot_ldblock,
        "count_indexed_vcf_records",
        lambda *_args, **_kwargs: 1,
    )
    monkeypatch.setattr(
        annot_ldblock,
        "select_vcf_sample",
        lambda *_args, **_kwargs: "STUDY",
    )

    def run(pipeline, _purpose, **_kwargs):
        commands.extend(pipeline)
        destination = Path(pipeline[-1][pipeline[-1].index("-o") + 1])
        destination.write_bytes(b"annotated-vcf")
        Path(str(destination) + ".tbi").write_bytes(b"annotated-index")
        return tuple(
            subprocess.CompletedProcess(command, 0, stdout="", stderr="")
            for command in pipeline
        )

    monkeypatch.setattr(annot_ldblock, "run_checked_pipeline", run)
    monkeypatch.setattr(
        annot_ldblock,
        "_validate_annotated_output",
        lambda *_args, **_kwargs: 1,
    )

    def summarize(_vcf, **kwargs):
        populations = {}
        for population in kwargs["populations"]:
            reference = kwargs["references"][population]
            populations[population] = {
                "annotated_variants": 1,
                "unassigned_variants": 0,
                "annotation_percent": 100.0,
                "bed_blocks": reference.block_count,
                "unique_bed_labels": len(reference.labels),
                "blocks_used": min(1, len(reference.labels)),
                "empty_blocks": max(0, len(reference.labels) - 1),
                "duplicate_bed_labels": (
                    reference.block_count - len(reference.labels)
                ),
                "unexpected_annotation_labels": 0,
                "bed_file": str(reference.path),
            }
        return {
            "dataset_id": kwargs["dataset_id"],
            "genome_build": kwargs["genome_build"],
            "total_variants": 1,
            "any_population_annotated": 1,
            "all_populations_annotated": 1,
            "fully_unassigned_variants": 0,
            "unassigned_on_bed_contigs": 0,
            "unassigned_outside_bed_contigs": 0,
            "unassigned_by_contig": {},
            "bed_contigs": list(kwargs["vcf_contigs"]),
            "populations": populations,
        }

    monkeypatch.setattr(
        annot_ldblock,
        "calculate_ld_annotation_summary",
        summarize,
    )


def _captured_progress() -> tuple[StageProgress, StringIO]:
    stream = StringIO()
    console = Console(
        file=stream,
        force_terminal=False,
        color_system=None,
        width=120,
    )
    return (
        StageProgress(
            "LD-block annotation progress",
            enabled=True,
            outcome_label_width=42,
            console=console,
        ),
        stream,
    )


@pytest.mark.parametrize(
    ("present", "requested", "expected"),
    (
        (0, 1, "Contig absent from the requested BED file"),
        (1, 1, "Contig present in the requested BED file"),
        (0, 3, "Contig absent from all 3 requested BED files"),
        (2, 3, "Contig present in 2 of 3 requested BED files"),
        (3, 3, "Contig present in all 3 requested BED files"),
    ),
)
def test_bed_contig_availability_is_explicit(present, requested, expected):
    assert (
        annotation_reporting._bed_contig_availability(present, requested)
        == expected
    )


def test_annotation_summary_counts_variants_blocks_and_uncovered_contigs(
    monkeypatch,
    tmp_path,
):
    references = {
        "EUR": LDReferenceSummary(
            path=tmp_path / "GRCh38_EUR_ldetect.bed.gz",
            contigs=frozenset({"1"}),
            labels=frozenset({"EUR_A", "EUR_B"}),
            block_count=2,
        ),
        "AFR": LDReferenceSummary(
            path=tmp_path / "GRCh38_AFR_ldetect.bed.gz",
            contigs=frozenset({"1"}),
            labels=frozenset({"AFR_A", "AFR_B"}),
            block_count=2,
        ),
    }
    query_lines = (
        "1\tEUR_A\tAFR_A\n",
        "1\tEUR_B\t.\n",
        "1\t.\t.\n",
        "X\t.\t.\n",
    )

    def stream(arguments, _purpose, *, line_consumer, **_kwargs):
        assert arguments[:3] == ["bcftools", "query", "-f"]
        assert "%INFO/EUR_LDblock" in arguments[3]
        assert "%INFO/AFR_LDblock" in arguments[3]
        for line in query_lines:
            line_consumer(line)
        return len(query_lines)

    monkeypatch.setattr(
        annotation_reporting,
        "run_checked_line_command",
        stream,
    )
    summary = annotation_reporting.calculate_ld_annotation_summary(
        tmp_path / "annotated.vcf.gz",
        bcftools="bcftools",
        dataset_id="STUDY",
        genome_build="GRCh38",
        populations=("EUR", "AFR"),
        info_ids={"EUR": "EUR_LDblock", "AFR": "AFR_LDblock"},
        references=references,
        expected_variant_count=4,
        vcf_contigs=("1", "X"),
    )

    assert summary["total_variants"] == 4
    assert summary["any_population_annotated"] == 2
    assert summary["all_populations_annotated"] == 1
    assert summary["fully_unassigned_variants"] == 2
    assert summary["unassigned_on_bed_contigs"] == 1
    assert summary["unassigned_outside_bed_contigs"] == 1
    assert summary["unassigned_by_contig"] == {"1": 1, "X": 1}
    assert summary["populations"]["EUR"]["annotated_variants"] == 2
    assert summary["populations"]["EUR"]["blocks_used"] == 2
    assert summary["populations"]["EUR"]["empty_blocks"] == 0
    assert summary["populations"]["AFR"]["annotated_variants"] == 1
    assert summary["populations"]["AFR"]["blocks_used"] == 1
    assert summary["populations"]["AFR"]["empty_blocks"] == 1

    summary_file = annotation_reporting.write_ld_annotation_summary(
        summary,
        tmp_path / "STUDY_ldblock_summary.csv",
    )
    with summary_file.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["record_type"] for row in rows] == [
        "overall",
        "population",
        "population",
        "unassigned_contig",
        "unassigned_contig",
    ]
    assert rows[0]["any_population_annotated"] == "2"
    assert rows[0]["all_populations_annotated"] == "1"
    assert rows[0]["unassigned_on_bed_contigs"] == "1"
    assert rows[0]["unassigned_outside_bed_contigs"] == "1"
    assert rows[1]["bed_blocks"] == "2"
    assert rows[1]["unique_bed_labels"] == "2"
    rendered = annotation_reporting.render_ld_annotation_summary(
        summary,
        summary_file,
        label_width=30,
        separator_column=60,
        annotated_vcf=tmp_path / "STUDY_ldblock.vcf.gz",
        annotated_index=tmp_path / "STUDY_ldblock.vcf.gz.tbi",
        html_report=tmp_path / "STUDY_ldblock_report.html",
        log_file=tmp_path / "STUDY_ldblock.log",
    )
    assert "LD-block annotation summary" in rendered
    assert "Assigned in any population" in rendered
    assert "Assigned in every population" in rendered
    assert "Unassigned in every population" in rendered
    assert "EUR annotated" in rendered
    assert "EUR unassigned" in rendered
    assert "2 used / 2 total; 0 empty" in rendered
    assert "Contig X without LD blocks" in rendered
    assert summary_file.name in rendered
    assert "STUDY_ldblock_report.html" in rendered
    assert "STUDY_ldblock.log" in rendered
    assert "STUDY_ldblock.vcf.gz" in rendered
    assert "STUDY_ldblock.vcf.gz.tbi" in rendered
    assert str(summary_file.parent) not in rendered
    summary_line = next(
        line for line in rendered.splitlines() if "Summary CSV" in line
    )
    assert summary_line.endswith(summary_file.name)
    rendered_lines = rendered.splitlines()
    population_heading = next(
        index
        for index, line in enumerate(rendered_lines)
        if "Population coverage" in line
    )
    summary_output = next(
        index
        for index, line in enumerate(rendered_lines)
        if "Summary CSV" in line
    )
    saved_outputs = next(
        index
        for index, line in enumerate(rendered_lines)
        if "Saved outputs" in line
    )
    afr_population = next(
        index
        for index, line in enumerate(rendered_lines)
        if "AFR annotated" in line
    )
    assert rendered_lines[population_heading - 1] == ""
    assert rendered_lines[afr_population - 1] == ""
    assert rendered_lines[saved_outputs - 1] == ""
    assert summary_output == saved_outputs + 3
    assert rendered.endswith("\n\n")
    field_lines = [line for line in rendered_lines if " : " in line]
    assert {
        cell_len(line.split(" : ", 1)[0]) for line in field_lines
    } == {60}
    assert len({len(line.split(" : ", 1)[0]) for line in field_lines}) == 1
    assert {len(line) - len(line.lstrip()) for line in field_lines} == {6}

    html_file = annotation_reporting.write_ld_annotation_html_report(
        summary,
        tmp_path / "STUDY_ldblock_report.html",
        input_vcf=tmp_path / "STUDY<input>.vcf.gz",
        vcf_contigs=("1", "X"),
        info_ids={"EUR": "EUR_LDblock", "AFR": "AFR_LDblock"},
        references=references,
        threads=4,
        bcftools="/tools/bcftools",
        outputs={
            "annotated_vcf": tmp_path / "STUDY_ldblock.vcf.gz",
            "annotated_index": tmp_path / "STUDY_ldblock.vcf.gz.tbi",
            "summary_file": summary_file,
            "html_report": tmp_path / "STUDY_ldblock_report.html",
            "log_file": tmp_path / "STUDY_ldblock.log",
        },
    )
    html = html_file.read_text(encoding="utf-8")
    assert "LD-block annotation report" in html
    assert "LD annotation does not remove variants" in html
    assert "Population coverage and block use" in html
    assert "EUR_LDblock" in html
    assert "AFR_LDblock" in html
    assert "Population LD-block reference inventory" in html
    assert "BED interval rows" in html
    assert "Unique block labels" in html
    assert "Variants with no LD-block assignment in any requested population" in html
    assert "Variants with no LD-block assignment" in html
    assert "Share of all unassigned variants" in html
    assert "LD-block BED availability" in html
    assert "Contig present in all 2 requested BED files" in html
    assert "Contig absent from all 2 requested BED files" in html
    assert "Reference status" not in html
    assert "BED contigs" in html
    assert "Validated outputs" in html
    assert "STUDY&lt;input&gt;.vcf.gz" in html
    assert "STUDY<input>.vcf.gz" not in html


def test_annotation_selects_bed_build_from_vcf_header(monkeypatch, tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    beds = tmp_path / "beds"
    beds.mkdir()
    expected_bed = beds / "blocks_GRCh38_EUR.bed.gz"
    _write_bed(expected_bed, "1\t100\t200\tEUR_1_100_200\n")
    commands = []
    _stub_external_execution(
        monkeypatch,
        "##fileformat=VCFv4.2\n" + STUDY_PROVENANCE +
        "##genome_build=GRCh38\n"
        "##contig=<ID=1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        commands,
    )

    outputs = annot_ldblock.annotate_ldblocks(
        vcf_path=str(vcf),
        output_directory=str(tmp_path / "output"),
        ld_dir=str(beds),
        genome_build_header="##genome_build={build}",
        supported_genome_builds=("GRCh37", "GRCh38"),
        bed_filename_template="blocks_{genome_build}_{population}.bed.gz",
        info_field_template=INFO_FIELD_TEMPLATE,
        info_description_template=INFO_DESCRIPTION_TEMPLATE,
        output_filename_template=OUTPUT_FILENAME_TEMPLATE,
        summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
        html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
        bcftools_bin="bcftools",
        provenance_headers=PROVENANCE_HEADERS,
        populations=("EUR",),
        threads=2,
        dataset_id="STUDY",
    )

    annotate_command = next(command for command in commands if "annotate" in command)
    assert annotate_command[annotate_command.index("-a") + 1] == str(expected_bed)
    assert "blocks_GRCh37_EUR.bed.gz" not in " ".join(annotate_command)
    assert outputs["annotated_vcf"].is_file()
    assert outputs["annotated_index"].is_file()
    assert outputs["summary_file"].is_file()
    assert outputs["html_report"].is_file()
    assert "Population coverage and block use" in outputs[
        "html_report"
    ].read_text(encoding="utf-8")
    assert outputs["variant_count"] == 1


@pytest.mark.parametrize("role", ("version", "dataset_id", "status"))
@pytest.mark.parametrize("value", (None, '""', '"unterminated', '<ID=invalid>'))
def test_direct_annotation_rejects_invalid_provenance_before_work(
    monkeypatch, tmp_path, role, value,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"header-query-fixture")
    beds = tmp_path / "beds"
    beds.mkdir()
    output = tmp_path / "output"
    configuration = load_configuration(cli_overrides={
        "modules.ld_annotation.inputs.vcf": str(vcf),
        "modules.ld_annotation.inputs.ld_region_dir": str(beds),
        "modules.ld_annotation.inputs.dataset_id": "STUDY",
        "modules.ld_annotation.output_directory": str(output),
    })
    metadata = {name: '"recorded"' for name in PROVENANCE_HEADERS.values()}
    name = PROVENANCE_HEADERS[role]
    if value is None:
        del metadata[name]
    else:
        metadata[name] = value
    header = "##fileformat=VCFv4.2\n" + "".join(
        f"##{key}={payload}\n" for key, payload in metadata.items()
    ) + "##genome_build=GRCh38\n##contig=<ID=1,length=248956422>\n#CHROM\n"
    header_queries = []

    def read_header(*_args, **_kwargs):
        header_queries.append(True)
        return header

    def forbidden(*_args, **_kwargs):
        pytest.fail("Invalid provenance reached index/sample/reference/analysis work")

    monkeypatch.setattr(annot_ldblock, "resolve_executable", lambda value, _label: value)
    monkeypatch.setattr(annot_ldblock, "read_vcf_header", read_header)
    for operation in (
        "count_indexed_vcf_records", "select_vcf_sample",
        "validate_ld_block_references", "run_checked_pipeline",
    ):
        monkeypatch.setattr(annot_ldblock, operation, forbidden)

    with pytest.raises(RuntimeError, match=name):
        annotation_service.run_annot_ldblock(
            argparse.Namespace(), configuration=configuration,
        )

    assert header_queries == [True]
    log_path = output / configuration.modules.ld_annotation.canonical_log_filename_template.format(
        dataset_id="STUDY",
    )
    log_text = log_path.read_text(encoding="utf-8")
    assert "postgwas_vcf_provenance" in log_text and "FAILED" in log_text
    assert name in log_text and str(vcf) in log_text
    assert set(output.iterdir()) == {log_path}


def test_missing_vcf_build_fails_before_output_creation(monkeypatch, tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    output = tmp_path / "output"
    commands = []
    _stub_external_execution(
        monkeypatch,
        "##fileformat=VCFv4.2\n" + STUDY_PROVENANCE +
        "##contig=<ID=1,length=249250621>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        commands,
    )

    with pytest.raises(ValueError, match="exactly one ##genome_build=<build>"):
        annot_ldblock.annotate_ldblocks(
            vcf_path=str(vcf),
            output_directory=str(output),
            ld_dir=str(tmp_path / "beds"),
            genome_build_header="##genome_build={build}",
            supported_genome_builds=("GRCh37", "GRCh38"),
            bed_filename_template=BED_FILENAME_TEMPLATE,
            info_field_template=INFO_FIELD_TEMPLATE,
            info_description_template=INFO_DESCRIPTION_TEMPLATE,
            output_filename_template=OUTPUT_FILENAME_TEMPLATE,
            summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
            html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
            bcftools_bin="bcftools",
            provenance_headers=PROVENANCE_HEADERS,
            populations=("EUR",),
            threads=2,
            dataset_id="STUDY",
        )

    assert not output.exists()
    assert commands == []


@pytest.mark.parametrize(
    ("vcf_contig", "bed_contig"),
    (("1", "chr1"), ("chr1", "1")),
)
def test_chromosome_naming_mismatch_fails_before_output_creation(
    monkeypatch,
    tmp_path,
    vcf_contig,
    bed_contig,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    beds = tmp_path / "beds"
    beds.mkdir()
    bed = beds / "GRCh38_EUR_ldetect.bed.gz"
    _write_bed(bed, f"{bed_contig}\t100\t200\tEUR_1_100_200\n")
    output = tmp_path / "output"
    commands = []
    _stub_external_execution(
        monkeypatch,
        "##fileformat=VCFv4.2\n" + STUDY_PROVENANCE +
        "##genome_build=GRCh38\n"
        f"##contig=<ID={vcf_contig},length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        commands,
    )

    with pytest.raises(ValueError, match="Chromosome names are incompatible"):
        annot_ldblock.annotate_ldblocks(
            vcf_path=str(vcf),
            output_directory=str(output),
            ld_dir=str(beds),
            genome_build_header="##genome_build={build}",
            supported_genome_builds=("GRCh37", "GRCh38"),
            bed_filename_template=BED_FILENAME_TEMPLATE,
            info_field_template=INFO_FIELD_TEMPLATE,
            info_description_template=INFO_DESCRIPTION_TEMPLATE,
            output_filename_template=OUTPUT_FILENAME_TEMPLATE,
            summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
            html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
            bcftools_bin="bcftools",
            provenance_headers=PROVENANCE_HEADERS,
            populations=("EUR",),
            threads=2,
            dataset_id="STUDY",
        )

    assert not output.exists()
    assert commands == []


def test_partial_mixed_chromosome_naming_is_rejected():
    with pytest.raises(ValueError, match="1 versus chr1"):
        validate_vcf_bed_contigs(
            ["1", "2"],
            {"chr1", "2"},
            Path("blocks.bed.gz"),
        )


def test_matching_chr_prefixed_names_are_accepted():
    validate_vcf_bed_contigs(
        ["chr1"],
        {"chr1"},
        Path("blocks.bed.gz"),
    )


def test_all_population_beds_are_preflighted_before_annotation(monkeypatch, tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    beds = tmp_path / "beds"
    beds.mkdir()
    _write_bed(
        beds / "GRCh38_EUR_ldetect.bed.gz",
        "1\t100\t200\tEUR_1_100_200\n",
    )
    _write_bed(
        beds / "GRCh38_AFR_ldetect.bed.gz",
        "chr1\t100\t200\tAFR_1_100_200\n",
    )
    output = tmp_path / "output"
    commands = []
    _stub_external_execution(
        monkeypatch,
        "##fileformat=VCFv4.2\n" + STUDY_PROVENANCE +
        "##genome_build=GRCh38\n"
        "##contig=<ID=1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        commands,
    )

    with pytest.raises(ValueError, match="GRCh38_AFR_ldetect"):
        annot_ldblock.annotate_ldblocks(
            vcf_path=str(vcf),
            output_directory=str(output),
            ld_dir=str(beds),
            genome_build_header="##genome_build={build}",
            supported_genome_builds=("GRCh37", "GRCh38"),
            bed_filename_template=BED_FILENAME_TEMPLATE,
            info_field_template=INFO_FIELD_TEMPLATE,
            info_description_template=INFO_DESCRIPTION_TEMPLATE,
            output_filename_template=OUTPUT_FILENAME_TEMPLATE,
            summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
            html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
            bcftools_bin="bcftools",
            provenance_headers=PROVENANCE_HEADERS,
            populations=("EUR", "AFR"),
            threads=2,
            dataset_id="STUDY",
        )

    assert not output.exists()
    assert commands == []


def test_missing_vcf_contigs_fails_before_bed_or_output_access(monkeypatch, tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    output = tmp_path / "output"
    commands = []
    _stub_external_execution(
        monkeypatch,
        "##fileformat=VCFv4.2\n" + STUDY_PROVENANCE +
        "##genome_build=GRCh38\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        commands,
    )

    with pytest.raises(ValueError, match="no ##contig declarations"):
        annot_ldblock.annotate_ldblocks(
            vcf_path=str(vcf),
            output_directory=str(output),
            ld_dir=str(tmp_path / "beds"),
            genome_build_header="##genome_build={build}",
            supported_genome_builds=("GRCh37", "GRCh38"),
            bed_filename_template=BED_FILENAME_TEMPLATE,
            info_field_template=INFO_FIELD_TEMPLATE,
            info_description_template=INFO_DESCRIPTION_TEMPLATE,
            output_filename_template=OUTPUT_FILENAME_TEMPLATE,
            summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
            html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
            bcftools_bin="bcftools",
            provenance_headers=PROVENANCE_HEADERS,
            populations=("EUR",),
            threads=2,
            dataset_id="STUDY",
        )

    assert not output.exists()
    assert commands == []


def test_real_bcftools_uses_vcf_declared_build_for_bed_selection(tmp_path):
    bcftools = shutil.which("bcftools")
    bgzip = shutil.which("bgzip")
    if not all((bcftools, bgzip)):
        pytest.skip("bcftools and bgzip are required for integration test")

    plain_vcf = tmp_path / "study.vcf"
    plain_vcf.write_text(
        "##fileformat=VCFv4.2\n" + STUDY_PROVENANCE +
        "##genome_build=GRCh38\n"
        "##contig=<ID=1,length=248956422>\n"
        "##FORMAT=<ID=ES,Number=1,Type=Float,Description=\"Effect size\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY\n"
        "1\t101\t.\tA\tG\t.\tPASS\t.\tES\t0.1\n",
        encoding="utf-8",
    )
    vcf = tmp_path / "study.vcf.gz"
    with vcf.open("wb") as output_handle:
        subprocess.run(
            [bgzip, "-c", str(plain_vcf)],
            check=True,
            stdout=output_handle,
        )
    subprocess.run([bcftools, "index", "--tbi", str(vcf)], check=True)

    beds = tmp_path / "beds"
    beds.mkdir()
    plain_bed = beds / "GRCh38_EUR_ldetect.bed"
    plain_bed.write_text("1\t100\t200\tEUR_1_100_200\n", encoding="utf-8")
    bed = beds / "GRCh38_EUR_ldetect.bed.gz"
    with bed.open("wb") as output_handle:
        subprocess.run(
            [bgzip, "-c", str(plain_bed)],
            check=True,
            stdout=output_handle,
        )
    plain_afr_bed = beds / "GRCh38_AFR_ldetect.bed"
    plain_afr_bed.write_text("1\t100\t200\tAFR_1_100_200\n", encoding="utf-8")
    afr_bed = beds / "GRCh38_AFR_ldetect.bed.gz"
    with afr_bed.open("wb") as output_handle:
        subprocess.run(
            [bgzip, "-c", str(plain_afr_bed)],
            check=True,
            stdout=output_handle,
        )

    outputs = annot_ldblock.annotate_ldblocks(
        vcf_path=str(vcf),
        output_directory=str(tmp_path / "output"),
        ld_dir=str(beds),
        genome_build_header="##genome_build={build}",
        supported_genome_builds=("GRCh37", "GRCh38"),
        bed_filename_template=BED_FILENAME_TEMPLATE,
        info_field_template=INFO_FIELD_TEMPLATE,
        info_description_template=INFO_DESCRIPTION_TEMPLATE,
        output_filename_template=OUTPUT_FILENAME_TEMPLATE,
        summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
        html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
        bcftools_bin=bcftools,
        provenance_headers=PROVENANCE_HEADERS,
        populations=("EUR", "AFR"),
        threads=2,
        dataset_id="STUDY",
    )

    queried = subprocess.run(
        [
            bcftools,
            "query",
            "-f",
            "%INFO/EUR_LDblock\t%INFO/AFR_LDblock\\n",
            str(outputs["annotated_vcf"]),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert queried.stdout == "EUR_1_100_200\tAFR_1_100_200\n"
    assert outputs["annotated_index"].is_file()
    assert outputs["summary_file"].is_file()
    assert outputs["html_report"].is_file()
    assert outputs["variant_count"] == 1
    with outputs["summary_file"].open(encoding="utf-8", newline="") as handle:
        summary_rows = list(csv.DictReader(handle))
    population_rows = {
        row["population"]: row
        for row in summary_rows
        if row["record_type"] == "population"
    }
    assert set(population_rows) == {"EUR", "AFR"}
    assert all(row["annotated_variants"] == "1" for row in population_rows.values())
    assert all(row["bed_blocks"] == "1" for row in population_rows.values())
    assert all(row["blocks_used"] == "1" for row in population_rows.values())
    assert all(row["empty_blocks"] == "0" for row in population_rows.values())

    rerun = annot_ldblock.annotate_ldblocks(
        vcf_path=str(outputs["annotated_vcf"]),
        output_directory=str(tmp_path / "output"),
        ld_dir=str(beds),
        genome_build_header="##genome_build={build}",
        supported_genome_builds=("GRCh37", "GRCh38"),
        bed_filename_template=BED_FILENAME_TEMPLATE,
        info_field_template=INFO_FIELD_TEMPLATE,
        info_description_template=INFO_DESCRIPTION_TEMPLATE,
        output_filename_template=OUTPUT_FILENAME_TEMPLATE,
        summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
        html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
        bcftools_bin=bcftools,
        provenance_headers=PROVENANCE_HEADERS,
        populations=("EUR", "AFR"),
        threads=2,
        dataset_id="STUDY",
    )
    rerun_header = subprocess.run(
        [bcftools, "view", "--header-only", str(rerun["annotated_vcf"])],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    header_lines = rerun_header.splitlines()
    assert sum(
        line.startswith("##INFO=<ID=EUR_LDblock,") for line in header_lines
    ) == 1
    assert sum(
        line.startswith("##INFO=<ID=AFR_LDblock,") for line in header_lines
    ) == 1


def test_multiple_populations_use_one_uncompressed_pipeline(monkeypatch, tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    beds = tmp_path / "beds"
    beds.mkdir()
    for population in ("EUR", "AFR", "EAS"):
        _write_bed(
            beds / ("GRCh38_%s_ldetect.bed.gz" % population),
            "1\t100\t200\t%s_1_100_200\n" % population,
        )
    commands = []
    progress, stream = _captured_progress()
    _stub_external_execution(
        monkeypatch,
        "##fileformat=VCFv4.2\n" + STUDY_PROVENANCE +
        "##genome_build=GRCh38\n"
        "##contig=<ID=1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        commands,
    )

    annot_ldblock.annotate_ldblocks(
        vcf_path=str(vcf),
        output_directory=str(tmp_path / "output"),
        ld_dir=str(beds),
        genome_build_header="##genome_build={build}",
        supported_genome_builds=("GRCh37", "GRCh38"),
        bed_filename_template=BED_FILENAME_TEMPLATE,
        info_field_template=INFO_FIELD_TEMPLATE,
        info_description_template=INFO_DESCRIPTION_TEMPLATE,
        output_filename_template=OUTPUT_FILENAME_TEMPLATE,
        summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
        html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
        bcftools_bin="bcftools",
        provenance_headers=PROVENANCE_HEADERS,
        populations=("EUR", "AFR", "EAS"),
        threads=3,
        dataset_id="STUDY",
        stage_progress=progress,
    )

    assert len(commands) == 3
    assert commands[0][-2:] == ["-Ou", str(vcf.resolve())]
    assert commands[1][-2:] == ["-Ou", "-"]
    assert "-Oz" in commands[2]
    assert "--write-index=tbi" in commands[2]
    assert commands[2][-1] == "-"
    assert "--threads" not in commands[0]
    assert "--threads" not in commands[1]
    assert commands[2][commands[2].index("--threads") + 1] == "3"
    assert sum(command.count("-o") for command in commands) == 1
    rendered = stream.getvalue()
    assert "Required references" in rendered
    assert "3/3 available" in rendered
    assert "EUR reference" in rendered
    assert "AFR reference" in rendered
    assert "EAS reference" in rendered
    assert "Populations annotated" in rendered
    assert "EUR, AFR, EAS" in rendered
    assert "Coverage summary" in rendered
    assert "calculated for 3 populations" in rendered
    assert "Reports prepared" in rendered
    assert "Analysis artifacts published" in rendered
    assert "Assigned in any population" not in rendered
    assert "EUR coverage" not in rendered
    progress_fields = [
        line for line in rendered.splitlines() if " : " in line
    ]
    assert {
        cell_len(line.split(" : ", 1)[0]) for line in progress_fields
    } == {progress.outcome_separator_column}


def test_annotation_reports_six_validated_progress_stages(monkeypatch, tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    beds = tmp_path / "beds"
    beds.mkdir()
    _write_bed(
        beds / "GRCh38_EUR_ldetect.bed.gz",
        "1\t100\t200\tEUR_1_100_200\n",
    )
    commands = []
    _stub_external_execution(
        monkeypatch,
        "##fileformat=VCFv4.2\n" + STUDY_PROVENANCE +
        "##genome_build=GRCh38\n"
        "##contig=<ID=1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        commands,
    )
    progress, stream = _captured_progress()
    log_file = tmp_path / "annotation.log"
    logger = PipelineLogger(
        "STUDY",
        "run",
        str(tmp_path),
        log_path=str(log_file),
        stage_progress=progress,
    )

    try:
        annot_ldblock.annotate_ldblocks(
            vcf_path=str(vcf),
            output_directory=str(tmp_path / "output"),
            ld_dir=str(beds),
            genome_build_header="##genome_build={build}",
            supported_genome_builds=("GRCh37", "GRCh38"),
            bed_filename_template=BED_FILENAME_TEMPLATE,
            info_field_template=INFO_FIELD_TEMPLATE,
            info_description_template=INFO_DESCRIPTION_TEMPLATE,
            output_filename_template=OUTPUT_FILENAME_TEMPLATE,
            summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
            html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
            bcftools_bin="bcftools",
            provenance_headers=PROVENANCE_HEADERS,
            populations=("EUR",),
            threads=2,
            dataset_id="STUDY",
            logger=logger,
            stage_progress=progress,
        )
    finally:
        logger.close()

    rendered = stream.getvalue()
    expected = (
        "Validate input VCF metadata",
        "Validate population LD-block references",
        "Annotate population LD blocks",
        "Validate annotated VCF and index",
        "Calculate and save annotation reports",
        "Publish validated VCF, index, and reports",
    )
    assert "LD-block annotation progress" in rendered
    for number, title in enumerate(expected, start=1):
        assert "Completed %d/6 · %s" % (number, title) in rendered
    assert "All 6 stages completed" in rendered
    assert "Genome build" in rendered
    stage_one = rendered.split(
        "Completed 1/6 · Validate input VCF metadata", 1
    )[1].split("Completed 2/6", 1)[0]
    stage_two = rendered.split(
        "Completed 2/6 · Validate population LD-block references", 1
    )[1].split("Completed 3/6", 1)[0]
    assert "Total variants" in stage_one
    assert "VCF sample columns" in stage_one
    assert "Populations requested" not in stage_one
    assert "Populations requested" in stage_two
    assert "Required references" in rendered
    assert "INFO fields created" in rendered
    assert "VCF variants validated" in rendered
    assert "Coverage summary" in rendered
    assert "calculated for 1 population" in rendered
    assert "Reports prepared" in rendered
    assert "Analysis artifacts published" in rendered
    assert "4/4" in rendered
    assert "EUR coverage" not in rendered
    assert "Detailed HTML report" not in rendered
    log_text = log_file.read_text(encoding="utf-8")
    assert log_text.count("STEP") == 6
    assert log_text.count("stage_outcome") == 6
    assert "population_coverage" in log_text
    assert "EUR_LDblock" in log_text


def test_annotation_progress_stops_at_failed_reference_stage(monkeypatch, tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    beds = tmp_path / "beds"
    beds.mkdir()
    commands = []
    _stub_external_execution(
        monkeypatch,
        "##fileformat=VCFv4.2\n" + STUDY_PROVENANCE +
        "##genome_build=GRCh38\n"
        "##contig=<ID=1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        commands,
    )
    progress, stream = _captured_progress()

    with pytest.raises(FileNotFoundError, match="LD block BED file missing"):
        annot_ldblock.annotate_ldblocks(
            vcf_path=str(vcf),
            output_directory=str(tmp_path / "output"),
            ld_dir=str(beds),
            genome_build_header="##genome_build={build}",
            supported_genome_builds=("GRCh37", "GRCh38"),
            bed_filename_template=BED_FILENAME_TEMPLATE,
            info_field_template=INFO_FIELD_TEMPLATE,
            info_description_template=INFO_DESCRIPTION_TEMPLATE,
            output_filename_template=OUTPUT_FILENAME_TEMPLATE,
            summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
            html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
            bcftools_bin="bcftools",
            provenance_headers=PROVENANCE_HEADERS,
            populations=("EUR",),
            threads=2,
            dataset_id="STUDY",
            stage_progress=progress,
        )

    rendered = stream.getvalue()
    assert "Completed 1/6 · Validate input VCF metadata" in rendered
    assert "Failed 2/6 · Validate population LD-block references" in rendered
    assert "Completed 2/6" not in rendered
    assert "All 6 stages completed" not in rendered
    assert commands == []
    assert not (tmp_path / "output").exists()


def test_annotation_rejects_changed_variant_count_before_publication(
    monkeypatch,
    tmp_path,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    beds = tmp_path / "beds"
    beds.mkdir()
    _write_bed(
        beds / "GRCh38_EUR_ldetect.bed.gz",
        "1\t100\t200\tEUR_1_100_200\n",
    )
    output = tmp_path / "output"
    commands = []
    _stub_external_execution(
        monkeypatch,
        "##fileformat=VCFv4.2\n" + STUDY_PROVENANCE +
        "##genome_build=GRCh38\n"
        "##contig=<ID=1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        commands,
    )
    monkeypatch.setattr(
        annot_ldblock,
        "_validate_annotated_output",
        lambda *_args, **_kwargs: 2,
    )

    with pytest.raises(
        ValueError,
        match=(
            "input VCF has 1 variants but the annotated VCF has 2.*"
            "must retain every input variant"
        ),
    ):
        annot_ldblock.annotate_ldblocks(
            vcf_path=str(vcf),
            output_directory=str(output),
            ld_dir=str(beds),
            genome_build_header="##genome_build={build}",
            supported_genome_builds=("GRCh37", "GRCh38"),
            bed_filename_template=BED_FILENAME_TEMPLATE,
            info_field_template=INFO_FIELD_TEMPLATE,
            info_description_template=INFO_DESCRIPTION_TEMPLATE,
            output_filename_template=OUTPUT_FILENAME_TEMPLATE,
            summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
            html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
            bcftools_bin="bcftools",
            provenance_headers=PROVENANCE_HEADERS,
            populations=("EUR",),
            threads=2,
            dataset_id="STUDY",
        )

    assert not (output / "STUDY_ldblock.vcf.gz").exists()
    assert not (output / "STUDY_ldblock.vcf.gz.tbi").exists()
    assert not list(output.glob(".*.tmp.vcf.gz*"))


@pytest.mark.parametrize(
    "failure_point",
    ("pipeline", "validation", "summary", "html_report"),
)
def test_failure_before_publication_preserves_existing_output_set(
    monkeypatch,
    tmp_path,
    failure_point,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"input-vcf")
    beds = tmp_path / "beds"
    beds.mkdir()
    _write_bed(
        beds / "GRCh38_EUR_ldetect.bed.gz",
        "1\t100\t200\tEUR_1_100_200\n",
    )
    output = tmp_path / "output"
    output.mkdir()
    final_vcf = output / "STUDY_ldblock.vcf.gz"
    final_index = output / "STUDY_ldblock.vcf.gz.tbi"
    final_summary = output / "STUDY_ldblock_summary.csv"
    final_html = output / "STUDY_ldblock_report.html"
    final_vcf.write_bytes(b"previous-vcf")
    final_index.write_bytes(b"previous-index")
    final_summary.write_bytes(b"previous-summary")
    final_html.write_bytes(b"previous-html")
    commands = []
    _stub_external_execution(
        monkeypatch,
        "##fileformat=VCFv4.2\n" + STUDY_PROVENANCE +
        "##genome_build=GRCh38\n"
        "##contig=<ID=1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        commands,
    )
    if failure_point == "pipeline":
        def fail_pipeline(*_args, **_kwargs):
            raise RuntimeError("pipeline failed")

        monkeypatch.setattr(
            annot_ldblock,
            "run_checked_pipeline",
            fail_pipeline,
        )
    elif failure_point == "validation":
        def fail_validation(*_args, **_kwargs):
            raise ValueError("validation failed")

        monkeypatch.setattr(
            annot_ldblock,
            "_validate_annotated_output",
            fail_validation,
        )
    elif failure_point == "summary":
        def fail_summary(*_args, **_kwargs):
            raise ValueError("summary failed")

        monkeypatch.setattr(
            annot_ldblock,
            "calculate_ld_annotation_summary",
            fail_summary,
        )
    else:
        def fail_html_report(*_args, **_kwargs):
            raise ValueError("html report failed")

        monkeypatch.setattr(
            annot_ldblock,
            "write_ld_annotation_html_report",
            fail_html_report,
        )

    with pytest.raises((RuntimeError, ValueError), match="failed"):
        annot_ldblock.annotate_ldblocks(
            vcf_path=str(vcf),
            output_directory=str(output),
            ld_dir=str(beds),
            genome_build_header="##genome_build={build}",
            supported_genome_builds=("GRCh37", "GRCh38"),
            bed_filename_template=BED_FILENAME_TEMPLATE,
            info_field_template=INFO_FIELD_TEMPLATE,
            info_description_template=INFO_DESCRIPTION_TEMPLATE,
            output_filename_template=OUTPUT_FILENAME_TEMPLATE,
            summary_filename_template=SUMMARY_FILENAME_TEMPLATE,
            html_report_filename_template=HTML_REPORT_FILENAME_TEMPLATE,
            bcftools_bin="bcftools",
            provenance_headers=PROVENANCE_HEADERS,
            populations=("EUR",),
            threads=2,
            dataset_id="STUDY",
        )

    assert final_vcf.read_bytes() == b"previous-vcf"
    assert final_index.read_bytes() == b"previous-index"
    assert final_summary.read_bytes() == b"previous-summary"
    assert final_html.read_bytes() == b"previous-html"
    assert not list(output.glob(".*.tmp.vcf.gz*"))
    assert not list(output.glob(".*.tmp.csv"))
    assert not list(output.glob(".*.tmp.html"))
