"""Fail-closed, provenance-recorded clients for FLAMES annotation protocols."""

import hashlib
import json
import math
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path

import requests

from postgwas.config.models.common import GenomeBuild
from postgwas.config.models.modules.flames import FlamesAnnotationAPI


class AnnotationAPIClient:
    def __init__(self, settings, audit_path):
        self.settings = FlamesAnnotationAPI.model_validate(settings)
        self.audit_path = Path(audit_path)
        self.audit_path.parent.mkdir(parents=True, exist_ok=True)
        self.last_request = None

    def _record(self, **fields):
        with self.audit_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps({
                "timestamp": datetime.now(timezone.utc).isoformat(), **fields,
            }, allow_nan=False) + "\n")

    def request(self, url):
        """Retry only explicit throttling; preserve every response before parsing."""
        for attempt in range(self.settings.retry_count + 1):
            if self.last_request is not None:
                delay = 1 / self.settings.requests_per_second - (time.monotonic() - self.last_request)
                if delay > 0:
                    time.sleep(delay)
            self.last_request = time.monotonic()
            try:
                response = requests.get(
                    url, headers={"Accept": "application/json"},
                    timeout=self.settings.timeout_seconds, allow_redirects=False,
                )
            except requests.RequestException as exc:
                self._record(url=url, attempt=attempt + 1, error=str(exc))
                raise RuntimeError("Annotation request failed; see %s" % self.audit_path) from exc
            self._record(
                url=url, attempt=attempt + 1, status=response.status_code,
                response_sha256=hashlib.sha256(response.content).hexdigest(),
                response_body=response.text,
            )
            if response.status_code != 429 or attempt == self.settings.retry_count:
                if not 200 <= response.status_code < 300:
                    raise requests.HTTPError(
                        "Annotation HTTP %s; see %s" % (response.status_code, self.audit_path),
                        response=response,
                    )
                try:
                    return response.json()
                except ValueError as exc:
                    raise RuntimeError("Malformed annotation JSON; see %s" % self.audit_path) from exc
            retry_after = response.headers.get("Retry-After")
            delay = self.settings.retry_delay_seconds
            if retry_after is not None:
                try:
                    delay = float(retry_after)
                except ValueError:
                    try:
                        delay = max(0, parsedate_to_datetime(retry_after).timestamp() - time.time())
                    except (ValueError, TypeError, OverflowError) as exc:
                        raise RuntimeError("Invalid annotation Retry-After header") from exc
                if not math.isfinite(delay) or delay < 0:
                    raise RuntimeError("Invalid annotation Retry-After header")
            time.sleep(delay)

    def vep(self, chromosome, position, a1, a2, build):
        server = self.settings.vep_servers[build].rstrip("/")
        # VEP rejects a reference-identical alternate with HTTP 400. Preserve
        # FLAMES' explicit alternate-allele retry, never its fabricated result.
        try:
            decoded = self.request(f"{server}/vep/human/region/{chromosome}:{position}:{position}/{a1}")
        except requests.HTTPError as exc:
            if exc.response.status_code != 400:
                raise
            decoded = self.request(f"{server}/vep/human/region/{chromosome}:{position}:{position}/{a2}")
        if not isinstance(decoded, list) or len(decoded) != 1 or not isinstance(decoded[0], dict):
            raise RuntimeError("VEP must return exactly one variant response")
        variant = decoded[0]
        if (
            variant.get("assembly_name") != build.value
            or str(variant.get("seq_region_name")) != str(chromosome)
            or variant.get("start") != int(position)
            or variant.get("end") != int(position)
            or not isinstance(variant.get("allele_string"), str)
            or set(variant["allele_string"].split("/")) != {a1, a2}
            or not isinstance(variant.get("most_severe_consequence"), str)
            or not variant["most_severe_consequence"]
        ):
            raise RuntimeError("VEP response has missing or incompatible build, coordinates, alleles or consequence")
        transcripts = variant.get("transcript_consequences", [])
        if not isinstance(transcripts, list):
            raise RuntimeError("VEP transcript consequences must be a list when present")
        for transcript in transcripts:
            if (
                not isinstance(transcript, dict)
                or not isinstance(transcript.get("gene_id"), str)
                or not transcript["gene_id"]
                or transcript.get("impact") not in {"HIGH", "MODERATE", "LOW", "MODIFIER"}
            ):
                raise RuntimeError("VEP transcript consequence is missing a gene or valid impact")
        # Intergenic/regulatory variants may legitimately have no transcripts.
        return decoded

    def cadd(self, chromosome, position, a1, a2, build):
        server = self.settings.cadd_servers[build].rstrip("/")
        decoded = self.request(f"{server}/{chromosome}:{position}")
        if not isinstance(decoded, list):
            raise RuntimeError("CADD response must be a list of variant scores")
        matches = []
        for variant in decoded:
            if not isinstance(variant, dict) or not {"Ref", "Alt", "Chrom", "Pos", "PHRED"}.issubset(variant):
                raise RuntimeError("CADD response is missing variant identity or PHRED")
            if str(variant["Chrom"]) != str(chromosome) or str(variant["Pos"]) != str(position):
                raise RuntimeError("CADD response coordinates do not match the requested variant")
            if {variant["Ref"], variant["Alt"]} == {a1, a2}:
                try:
                    score = float(variant["PHRED"])
                except (ValueError, TypeError) as exc:
                    raise RuntimeError("CADD PHRED must be finite and nonnegative") from exc
                if not math.isfinite(score) or score < 0:
                    raise RuntimeError("CADD PHRED must be finite and nonnegative")
                matches.append(score)
        if len(matches) != 1:
            raise RuntimeError("CADD must return exactly one matching allele score; missing scores are not zero")
        return matches[0]


_client = None


def configure_api(settings, audit_path):
    global _client
    _client = AnnotationAPIClient(settings, audit_path)


def _configured_client(build):
    if _client is None:
        raise RuntimeError("FLAMES annotation API requires resolved configuration and an audit log")
    normalized = str(build).upper()
    matches = [item for item in GenomeBuild if item.value.upper() == normalized]
    if len(matches) != 1:
        raise ValueError("Unsupported annotation genome build: %s" % build)
    return _client, matches[0]


def query_VEP(chr, pos, a1, a2, build):
    client, genome_build = _configured_client(build)
    return client.vep(chr, pos, a1, a2, genome_build)


def query_CADD(chr, pos, a1, a2, build):
    client, genome_build = _configured_client(build)
    return client.cadd(chr, pos, a1, a2, genome_build)
