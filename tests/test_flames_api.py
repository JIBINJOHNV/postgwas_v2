"""Local fixtures only: no annotation requests leave the test process."""

import importlib
import json
import sys

import pandas as pd
import pytest
import requests
from pydantic import ValidationError

from postgwas.config import load_module_configuration
from postgwas.config.models.common import GenomeBuild
from postgwas.config.models.modules.flames import FlamesAnnotationAPI
from postgwas.modules.flames import Query_api as api


def _settings():
    return load_module_configuration("flames").annotation_api.model_dump(mode="json")


def _response(data, status=200, headers=None, malformed=False):
    response = requests.Response()
    response.status_code = status
    response._content = b"not JSON" if malformed else json.dumps(data).encode()
    response.headers.update(headers or {})
    return response


def _vep(transcripts=None):
    value = {"assembly_name": "GRCh37", "seq_region_name": "1", "start": 10, "end": 10,
             "allele_string": "A/G", "most_severe_consequence": "intergenic_variant"}
    if transcripts is not None:
        value["transcript_consequences"] = transcripts
    return [value]


@pytest.fixture
def client(monkeypatch, tmp_path):
    monkeypatch.setattr(api.time, "sleep", lambda delay: None)
    settings = _settings()
    settings["retry_count"] = 2
    return api.AnnotationAPIClient(settings, tmp_path / "requests.jsonl")


def _mock(monkeypatch, responses):
    calls = []
    iterator = iter(responses)
    def get(url, **kwargs):
        calls.append((url, kwargs))
        item = next(iterator)
        if isinstance(item, Exception):
            raise item
        return item
    monkeypatch.setattr(api.requests, "get", get)
    return calls


def test_vep_valid_intergenic_is_not_missing_annotation(client, monkeypatch):
    value = _vep()
    calls = _mock(monkeypatch, [_response(value)])
    assert client.vep(1, 10, "A", "G", GenomeBuild.GRCH37) == value
    assert calls[0][0].startswith("https://grch37.rest.ensembl.org/")
    assert calls[0][1]["timeout"] == _settings()["timeout_seconds"]
    assert calls[0][1]["allow_redirects"] is False
    record = json.loads(client.audit_path.read_text())
    assert record["status"] == 200
    assert json.loads(record["response_body"]) == value
    assert len(record["response_sha256"]) == 64


def test_vep_throttle_retry_returns_result_and_preserves_attempts(client, monkeypatch):
    calls = _mock(monkeypatch, [_response({}, 429, {"Retry-After": "0"}), _response(_vep())])
    assert client.vep(1, 10, "A", "G", GenomeBuild.GRCH37) == _vep()
    assert len(calls) == 2
    assert [json.loads(line)["status"] for line in client.audit_path.read_text().splitlines()] == [429, 200]


def test_vep_400_retries_other_allele_but_never_fabricates(client, monkeypatch):
    calls = _mock(monkeypatch, [_response({}, 400), _response(_vep())])
    assert client.vep(1, 10, "A", "G", GenomeBuild.GRCH37) == _vep()
    assert [url.rsplit("/", 1)[1] for url, _ in calls] == ["A", "G"]


@pytest.mark.parametrize("status", [301, 403, 404, 429, 500])
def test_http_errors_fail_closed_with_bounded_retries(client, monkeypatch, status):
    calls = _mock(monkeypatch, [_response({}, status)] * 3)
    with pytest.raises(requests.HTTPError, match=str(status)):
        client.vep(1, 10, "A", "G", GenomeBuild.GRCH37)
    assert len(calls) == (3 if status == 429 else 1)
    assert client.audit_path.is_file()


@pytest.mark.parametrize("value", [None, [], {}, [{}], _vep([{}]), _vep([{"gene_id": "ENSG1", "impact": "UNKNOWN"}])])
def test_missing_or_malformed_vep_variant_is_not_intergenic(client, monkeypatch, value):
    _mock(monkeypatch, [_response(value)])
    with pytest.raises(RuntimeError):
        client.vep(1, 10, "A", "G", GenomeBuild.GRCH37)


@pytest.mark.parametrize("field,value", [("assembly_name", "GRCh38"), ("seq_region_name", "2"), ("start", 11), ("end", None), ("allele_string", "A/T"), ("transcript_consequences", None)])
def test_incompatible_vep_identity_or_schema(client, monkeypatch, field, value):
    decoded = _vep()
    decoded[0][field] = value
    _mock(monkeypatch, [_response(decoded)])
    with pytest.raises(RuntimeError):
        client.vep(1, 10, "A", "G", GenomeBuild.GRCH37)


@pytest.mark.parametrize("response", [requests.Timeout("fixture timeout"), _response({}, malformed=True)])
def test_transport_and_json_failures_leave_audit(client, monkeypatch, response):
    _mock(monkeypatch, [response])
    with pytest.raises(RuntimeError):
        client.vep(1, 10, "A", "G", GenomeBuild.GRCH37)
    assert client.audit_path.is_file()


def _cadd(score="12.3"):
    return [{"Chrom": "1", "Pos": "10", "Ref": "A", "Alt": "G", "PHRED": score}]


def test_cadd_valid_score_and_allele_order_are_preserved(client, monkeypatch):
    _mock(monkeypatch, [_response(_cadd())])
    assert client.cadd(1, 10, "G", "A", GenomeBuild.GRCH37) == 12.3


@pytest.mark.parametrize("value", [[], None, [{}], _cadd("nan"), _cadd(-1), _cadd() * 2])
def test_cadd_missing_invalid_or_ambiguous_score_never_becomes_zero(client, monkeypatch, value):
    _mock(monkeypatch, [_response(value)])
    with pytest.raises(RuntimeError):
        client.cadd(1, 10, "A", "G", GenomeBuild.GRCH37)


@pytest.mark.parametrize("field,value", [("Alt", "T"), ("Pos", "11"), ("Chrom", "2")])
def test_cadd_requires_matching_coordinates_and_alleles(client, monkeypatch, field, value):
    decoded = _cadd()
    decoded[0][field] = value
    _mock(monkeypatch, [_response(decoded)])
    with pytest.raises(RuntimeError):
        client.cadd(1, 10, "A", "G", GenomeBuild.GRCH37)


@pytest.mark.parametrize("endpoint", ["http://example.org", "https://user:secret@example.org", "https://example.org?token=secret"])
def test_config_rejects_insecure_or_authenticated_endpoints(endpoint):
    settings = _settings()
    settings["vep_servers"]["GRCh37"] = endpoint
    with pytest.raises(ValidationError):
        FlamesAnnotationAPI.model_validate(settings)


@pytest.mark.parametrize("field,value", [("timeout_seconds", 0), ("timeout_seconds", float("inf")), ("retry_count", -1), ("requests_per_second", 0)])
def test_config_rejects_invalid_request_limits(field, value):
    settings = _settings()
    settings[field] = value
    with pytest.raises(ValidationError):
        FlamesAnnotationAPI.model_validate(settings)


@pytest.fixture
def annotate(monkeypatch):
    monkeypatch.setitem(sys.modules, "Query_api", api)
    return importlib.import_module("postgwas.modules.flames.annotate")


def test_partial_vep_failure_aborts_locus_instead_of_zero_imputation(annotate, monkeypatch):
    responses = iter([_vep(), RuntimeError("fixture failed variant")])
    def query(*args):
        value = next(responses)
        if isinstance(value, Exception):
            raise value
        return value
    monkeypatch.setattr(api, "query_VEP", query)
    creds = pd.DataFrame({"chr": [1, 1], "pos": [10, 20], "a1": ["A", "A"], "a2": ["G", "G"], "prob1": [0.5, 0.5]})
    genes = pd.DataFrame({"ensg": ["ENSG1"]})
    with pytest.raises(RuntimeError, match="1:20:A:G"):
        annotate.get_VEP(creds, genes, "prob1", "GRCH37")
    assert "VEP_sum" not in genes


def test_valid_transcript_weights_and_intergenic_zero_are_unchanged(annotate, monkeypatch):
    responses = iter([_vep([{"gene_id": "ENSG1", "impact": "MODERATE"}]), _vep()])
    monkeypatch.setattr(api, "query_VEP", lambda *args: next(responses))
    creds = pd.DataFrame({"chr": [1, 1], "pos": [10, 20], "a1": ["A", "A"], "a2": ["G", "G"], "prob1": [0.5, 0.5]})
    genes = pd.DataFrame({"ensg": ["ENSG1", "ENSG2"]})
    result = annotate.get_VEP(creds, genes, "prob1", "GRCH37")
    assert result["VEP_sum"].tolist() == [0.3, 0.0]
    assert result["VEP_max"].tolist() == [0.3, 0.0]


def test_multiple_transcripts_and_variants_aggregate_once_per_gene(annotate, monkeypatch):
    responses = iter([
        _vep([{"gene_id": "ENSG1", "impact": "MODERATE"},
              {"gene_id": "ENSG1", "impact": "HIGH"},
              {"gene_id": "ENSG2", "impact": "LOW"}]),
        _vep([{"gene_id": "ENSG1", "impact": "MODIFIER"}]),
        _vep(),
    ])
    monkeypatch.setattr(api, "query_VEP", lambda *args: next(responses))
    creds = pd.DataFrame({"chr": [1, 1, 1], "pos": [10, 20, 30], "a1": ["A"] * 3,
                          "a2": ["G"] * 3, "prob1": [0.5, 0.3, 0.2]})
    genes = pd.DataFrame({"ensg": ["ENSG1", "ENSG2", "ENSG3"]})
    result = annotate.get_VEP(creds, genes, "prob1", "GRCH37")
    assert result.ensg.tolist() == genes.ensg.tolist()
    assert result.VEP_sum.tolist() == pytest.approx([0.53, 0.2, 0])
    assert result.VEP_max.tolist() == pytest.approx([0.5, 0.2, 0])
