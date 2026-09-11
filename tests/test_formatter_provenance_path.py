"""Formatter origin failures retain the known input path and audit logger."""

from pathlib import Path
from unittest.mock import Mock

import pytest

from postgwas.config import load_configuration
from postgwas.core.errors import FormattingError
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.preflight import PreflightLogBuffer
from postgwas.modules.formatting import table


@pytest.mark.parametrize("value", [None, '""'])
def test_formatter_provenance_failure_names_input_path_and_records_logger(tmp_path, monkeypatch, value):
    config = load_configuration().modules.formatting
    name = config.input_contract.provenance_headers.version
    fixture = (Path(__file__).parent / "fixtures/formatting/harmonised.vcf").read_text()
    header = "\n".join(line for line in fixture.splitlines() if not line.startswith("##" + name + "="))
    if value is not None:
        header += f"\n##{name}={value}\n"
    header_query = Mock(return_value=header)
    extraction = Mock(side_effect=AssertionError("Origin failure must precede record extraction"))
    monkeypatch.setattr(table, "read_vcf_header", header_query)
    monkeypatch.setattr(table, "extract_vcf_table", extraction)
    source = tmp_path / "foreign study.vcf.gz"
    projection = tmp_path / "projection.tsv"
    logger = PreflightLogBuffer()

    with InputValidationSession() as session:
        with pytest.raises(FormattingError) as caught:
            table.load_harmonised_vcf(source, projection, "bcftools", config, logger=logger)

    assert str(source) in str(caught.value) and name in str(caught.value)
    assert session.records[0].path == str(source.resolve())
    assert session.records[0].status == "failed"
    assert logger.events[0].arguments == ("VALIDATE", "postgwas_vcf_provenance")
    assert dict(logger.events[0].fields)["status"] == "FAILED"
    assert str(source) in dict(logger.events[0].fields)["message"]
    header_query.assert_called_once()
    extraction.assert_not_called()
    assert not projection.exists()
