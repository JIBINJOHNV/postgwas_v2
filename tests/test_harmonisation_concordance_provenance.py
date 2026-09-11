"""Standalone concordance proves study provenance before scientific I/O."""

import argparse
import json
from pathlib import Path
from unittest.mock import Mock

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.core import vcf as core_vcf
from postgwas.core.paths import configured_output_path
from postgwas.modules.harmonisation.concordance import cli, service
from postgwas.modules.harmonisation.sample_sheet import (
    HarmonisationSampleSheetRow,
    SAMPLE_SHEET_VERSION,
)


@pytest.fixture
def standalone(tmp_path, monkeypatch):
    config = load_configuration(cli_overrides={"run.output_directory": str(tmp_path / "results")})
    row = HarmonisationSampleSheetRow(
        config_version=SAMPLE_SHEET_VERSION, dataset_id="study",
        input_file=tmp_path / "input.tsv", chromosome_column="CHR", position_column="BP",
        effect_allele_column="EA", other_allele_column="OA",
        effect_allele_frequency_column="EAF", effect_column="BETA", effect_type="beta",
        standard_error_column="SE", z_score_column="Z", p_value_column="P", p_value_type="raw",
        delimiter="tab", control_count=1000,
    )
    args = argparse.Namespace(sample_sheet=str(tmp_path / "studies.csv"), dataset_id=row.dataset_id,
                              vcf=str(tmp_path / "study.vcf.gz"))
    # Only the mocked header query sees this synthetic primary-VCF fixture.
    Path(args.vcf).write_bytes(b"synthetic VCF fixture")
    monkeypatch.setattr(cli, "load_configuration", lambda *args, **kwargs: config)
    monkeypatch.setattr(cli, "load_harmonisation_sample_sheet", lambda path: [row])
    return config, row, args


def _header(config, *, changed_field=None, changed_value=None):
    values = {"version": "0.9.0", "dataset_id": "study", "status": "harmonised"}
    if changed_field is not None:
        values[changed_field] = changed_value
    names = config.modules.formatting.input_contract.provenance_headers.model_dump()
    return "##contig=<ID=1,assembly=GRCh37>\n" + "".join(
        f"##{names[field]}={json.dumps(value)}\n"
        for field, value in values.items() if value is not None
    )


def _output_paths(config, dataset_id):
    layout = config.modules.harmonisation.output_layout.root
    base = configured_output_path(config.run.output_directory, layout["dataset_directory"],
                                  dataset_id=dataset_id)
    return {
        key: configured_output_path(base, layout[key], dataset_id=dataset_id)
        for key in ("run_manifest", "concordance_log", "concordance_summary")
    }


@pytest.mark.parametrize("renamed_headers", [False, True])
@pytest.mark.parametrize("external_eaf", [False, True])
def test_standalone_reuses_valid_header_and_exempts_external_eaf(
    standalone, monkeypatch, renamed_headers, external_eaf,
):
    config, row, args = standalone
    if renamed_headers:
        contract = config.modules.formatting.input_contract
        contract.provenance_headers = type(contract.provenance_headers)(
            version="configured_origin_version", dataset_id="configured_origin_dataset",
            status="configured_origin_status",
        )
    study = pl.DataFrame({
        "CHR": [1, 1], "BP": [101, 102], "EA": ["A", "C"], "OA": ["G", "T"],
        "EAF": [0.2, 0.3], "BETA": [0.2, 0.3], "SE": [0.1, 0.1],
        "Z": [2.0, 3.0], "P": [0.05, 0.01],
    })
    if external_eaf:
        mapping = config.modules.harmonisation.external_eaf_mapping
        reference = row.input_file.parent / "reference.tsv"
        pl.DataFrame({
            mapping.chromosome: [1, 1], mapping.position: [101, 102],
            mapping.effect_allele: ["A", "C"], mapping.other_allele: ["G", "T"],
            "reference_frequency": [0.2, 0.3],
        }).write_csv(reference, separator="\t")
        row.effect_allele_frequency_column = None
        row.external_eaf_file = reference
        row.external_eaf_column = "reference_frequency"
        study = study.drop("EAF")
    study.write_csv(row.input_file, separator="\t")
    paths = _output_paths(config, row.dataset_id)
    paths["run_manifest"].parent.mkdir(parents=True, exist_ok=True)
    paths["run_manifest"].write_text(json.dumps({
        "dataset": {"genome_build": {"inferred_build": "GRCh37"}, "study_decisions": {
            "effect_type": "beta", "pvalue_type": "raw", "eaf_is_maf": False,
            "eaf_is_maf_source": "study_level_statistic", "strand": "forward",
        }},
    }))
    header = _header(config)
    header_query = Mock(return_value=header)
    monkeypatch.setattr(core_vcf, "run_checked_command", header_query)
    redundant_query = Mock(side_effect=AssertionError("Header must be reused after preflight"))
    monkeypatch.setattr(service, "run_checked_command", redundant_query)
    extracted = pl.DataFrame({
        "CHROM": [1, 1], "POS": [101, 102], "ID": ["v1", "v2"],
        "REF": ["G", "T"], "ALT": ["A", "C"], "ES": [0.2, 0.3],
        "SE": [0.1, 0.1], "EZ": [2.0, 3.0], "AF": [0.2, 0.3],
        "LP": [1.3010299956639813, 2.0],
    })

    def extract(_vcf, destination, _dataset, _columns, _bcftools, **kwargs):
        assert header_query.call_count == 1
        extracted.write_csv(destination, separator=kwargs["delimiter"])
        return str(destination)

    monkeypatch.setattr(service, "extract_vcf_table", extract)
    result = cli.run_standalone_validation(args)

    assert result["status"] == "PASS"
    assert result["summary"]["matched_variants"] == 2
    assert result["metrics"]["allele_frequency"]["concordant"] == 2
    header_query.assert_called_once()
    assert header_query.call_args.args[0] == [
        config.resources.executables.bcftools, "view", "--header-only", str(Path(args.vcf).resolve()),
    ]
    redundant_query.assert_not_called()
    log = Path(result["reports"]["log"]).read_text()
    assert "postgwas_vcf_provenance" in log and "PASSED" in log
    assert "0.9.0" in log


@pytest.mark.parametrize("field", ["version", "dataset_id", "status"])
@pytest.mark.parametrize("value", [None, "", "   "])
def test_standalone_rejects_missing_or_blank_provenance_before_work_and_finalizes_failure(
    standalone, monkeypatch, capsys, field, value,
):
    config, row, args = standalone
    header_query = Mock(return_value=_header(config, changed_field=field, changed_value=value))
    monkeypatch.setattr(core_vcf, "run_checked_command", header_query)
    work = {}
    for name in ("_read_run_manifest", "_stage_input", "_stage_duplicate_report", "_stage_strand_actions",
                 "_extract_vcf", "_stage_external_eaf", "_compare_staged_partitions"):
        work[name] = Mock(side_effect=AssertionError(f"{name} ran before provenance passed"))
        monkeypatch.setattr(service, name, work[name])
    monkeypatch.setattr("sys.argv", ["postgwas --validate", "--sample-sheet", args.sample_sheet,
                                    "--dataset-id", args.dataset_id, "--vcf", args.vcf])

    assert cli.main() == 1

    name = config.modules.formatting.input_contract.provenance_headers.model_dump()[field]
    stderr = capsys.readouterr().err
    assert name in stderr
    assert "PostGWAS-harmonised study VCFs only" in stderr
    assert "do not add provenance headers manually" in stderr
    header_query.assert_called_once()
    for operation in work.values():
        operation.assert_not_called()
    paths = _output_paths(config, row.dataset_id)
    log = paths["concordance_log"].read_text()
    summary = paths["concordance_summary"].read_text()
    assert name in log and "Validation failed" in log
    assert "Concordance validation finished" in log
    assert "\tstatus\tFAIL" in summary and name in summary
    assert not row.input_file.exists()
    assert list(paths["concordance_summary"].parent.iterdir()) == [paths["concordance_summary"]]
