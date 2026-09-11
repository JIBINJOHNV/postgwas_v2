"""User-facing and scientific-evidence contracts for harmonisation HTML."""

from copy import deepcopy
from html.parser import HTMLParser
from pathlib import Path
import re
from urllib.parse import unquote

import pytest

from postgwas.modules.harmonisation.html_report import (
    render_dataset_report,
    render_run_report,
    write_dataset_report,
)
from postgwas.modules.harmonisation.report_steps import (
    chromosome_steps,
    dataset_steps,
    post_merge_steps,
)
from test_harmonisation_html_report import _manifest, _record


class _ReportHTML(HTMLParser):
    """Inspect semantic HTML without a browser or third-party dependency."""

    def __init__(self, document):
        super().__init__()
        self.elements = []
        self.headings = []
        self.words = []
        self._heading = None
        self._heading_words = []
        self._ignore = False
        self.feed(document)

    def handle_starttag(self, tag, attrs):
        self.elements.append((tag, dict(attrs)))
        if tag in {"script", "style"}:
            self._ignore = True
        if re.fullmatch(r"h[1-6]", tag):
            self._heading = tag
            self._heading_words = []

    def handle_data(self, data):
        if not self._ignore:
            self.words.append(data)
            if self._heading:
                self._heading_words.append(data)

    def handle_endtag(self, tag):
        if tag in {"script", "style"}:
            self._ignore = False
        if tag == self._heading:
            self.headings.append((tag, " ".join("".join(self._heading_words).split())))
            self._heading = None

    @property
    def text(self):
        return " ".join(" ".join(self.words).split())


def test_all_eight_dataset_sections_are_visible_in_analysis_order():
    report = _ReportHTML(render_dataset_report(_record(), _manifest()))
    expected = (
        "Results at a glance",
        "Inputs and initial checks",
        "Whole-study preparation",
        "Chromosome harmonisation",
        "Chromosome completion and variant accounting",
        "Merging and final dataset QC",
        "Input–VCF concordance",
        "Final interpretation and downloads",
    )
    headings = [text for tag, text in report.headings if tag == "h3"]
    main = [text for text in headings if re.match(r"^[1-8][. ]", text)]
    assert len(main) == len(expected)
    for number, (actual, title) in enumerate(zip(main, expected), 1):
        assert actual.startswith(str(number))
        assert title in actual


def test_report_has_accessible_headings_scoped_tables_and_resolving_navigation():
    report = _ReportHTML(render_dataset_report(_record(), _manifest()))
    assert sum(tag == "h1" for tag, _ in report.headings) == 1
    assert ("html", {"lang": "en"}) in report.elements
    headers = [attrs for tag, attrs in report.elements if tag == "th"]
    assert headers
    assert all(attrs.get("scope") in {"row", "col"} for attrs in headers)
    ids = [attrs["id"] for _, attrs in report.elements if "id" in attrs]
    assert len(ids) == len(set(ids))
    fragments = [
        unquote(attrs["href"][1:])
        for tag, attrs in report.elements
        if tag == "a" and attrs.get("href", "").startswith("#")
    ]
    assert fragments
    assert all(fragment in ids for fragment in fragments)


@pytest.mark.parametrize("with_chromosome", [False, True])
def test_all_workflow_steps_are_rendered_even_before_chromosome_processing(with_chromosome):
    manifest = _manifest() if with_chromosome else {}
    summary = (_manifest()["dataset"]["chromosome_summaries"]["1"]
               if with_chromosome else {})
    report = _ReportHTML(render_dataset_report(_record(), manifest))
    stage_headings = [text for tag, text in report.headings if tag in {"h4", "h5"}]
    sequences = (
        dataset_steps(manifest),
        chromosome_steps(summary),
        post_merge_steps(manifest),
    )
    assert [len(sequence) for sequence in sequences] == [8, 16, 5]
    for steps in sequences:
        positions = []
        for step in steps:
            matches = [
                index for index, title in enumerate(stage_headings)
                if title.startswith(f"{step['number']}.") and step["title"] in title
            ]
            assert len(matches) == 1, step["title"]
            positions.append(matches[0])
        assert positions == sorted(positions)


def test_missing_evidence_is_not_rendered_as_zero_or_scientific_success():
    report = _ReportHTML(render_dataset_report(
        {"dataset_id": "not_started", "status": "NOT_RUN"}, {},
    ))
    assert "Not recorded" in report.text
    assert "0 / 0" not in report.text
    assert "0.0000%" not in report.text
    assert "PASSED" not in report.text


@pytest.mark.parametrize("requested", [False, True])
def test_concordance_not_requested_differs_from_requested_without_results(requested):
    manifest = _manifest()
    manifest["report_context"] = {"concordance_requested": requested}
    report = _ReportHTML(render_dataset_report(_record(), manifest))
    match = report.text.split("7. Input–VCF concordance")[-1]
    match = match.split("8. Final interpretation and downloads")[0]
    if requested:
        assert "Not requested" not in match
        assert "not recorded" in match.lower() or "not completed" in match.lower()
    else:
        assert "Not requested" in match


def test_run_search_has_an_accessible_label_and_preserves_all_dataset_content():
    names = ("study_one", "study two", "study-two", "study&lt;three&gt;")
    records = [_record(name) for name in names]
    records[-1]["status"] = "FAILED"
    records[-1]["failure_reason"] = "Missing required reference"
    manifests = {name: _manifest(name) for name in names}
    report = _ReportHTML(render_run_report(records, manifests))
    search = [attrs for tag, attrs in report.elements if tag == "input" and attrs.get("type") == "search"]
    assert len(search) == 1
    labels = {attrs.get("for") for tag, attrs in report.elements if tag == "label"}
    assert search[0].get("aria-label") or search[0].get("id") in labels
    ids = [attrs["id"] for _, attrs in report.elements if "id" in attrs]
    assert len(ids) == len(set(ids)), "Distinct dataset names must not collide after anchor encoding"
    assert len([
        title for _, title in report.headings
        if title.startswith("8.") and "Final interpretation and downloads" in title
    ]) == len(names)
    assert "Missing required reference" in report.text
    assert all(name in report.text for name in names)


def test_external_reference_identity_and_user_controlled_html_are_escaped():
    record = _record('study"><img src=x onerror="alert(1)">')
    record["input_file"] = "javascript:alert(1)"
    record["failure_reason"] = "invalid <input> & missing reference"
    manifest = _manifest()
    manifest["report_context"] = {
        "input_mapping": {
            "effect_column": "OR<script>alert(2)</script>",
            "external_info_file": '/reference/panel" onclick="alert(3).tsv.gz',
        }
    }
    document = render_dataset_report(record, manifest)
    report = _ReportHTML(document)
    assert not any(tag in {"img", "input"} for tag, _ in report.elements)
    assert all(
        not str(value).lower().startswith(("javascript:", "data:"))
        for _, attrs in report.elements for key, value in attrs.items()
        if key in {"href", "src"}
    )
    assert not any(key.startswith("on") for _, attrs in report.elements for key in attrs)
    assert "invalid &lt;input&gt; &amp; missing reference" in document
    assert "&lt;script&gt;alert(2)&lt;/script&gt;" in document


def test_existing_output_links_are_relative_and_missing_intermediates_are_not_linked(tmp_path):
    outputs = tmp_path / "study" / "harmonisation"
    outputs.mkdir(parents=True)
    vcf = outputs / "study GRCh38.vcf.gz"
    vcf.write_bytes(b"fixture")
    missing = outputs / "study_gwas2vcf_GRCh38.vcf.gz"
    record, manifest = _record(), _manifest()
    manifest["merged_vcfs"] = {"grch38": str(vcf), "gwas2vcf": str(missing)}
    manifest["gwas2vcf_intermediate"] = {
        "path": str(missing),
        "retained": False,
        "removed_artifacts": [str(missing)],
        "retention_reason": "successful_run_default_cleanup",
    }
    destination = outputs / "study_report.html"
    report = _ReportHTML(write_dataset_report(destination, record, manifest).read_text())
    links = [attrs["href"] for tag, attrs in report.elements if tag == "a" and attrs.get("href")]
    assert any(unquote(href) == vcf.name for href in links)
    assert not any(missing.name in unquote(href) for href in links)
    assert "removed" in report.text.lower() and "policy" in report.text.lower()


def test_custom_qc_rule_criteria_and_virtual_subset_meaning_are_displayed():
    manifest = _manifest()
    manifest["report_context"] = {
        "qc_assessment": {
            "raw": {"num_records": 94},
            "qc_passed": {"num_records": 88},
            "rules": [{
                "key": "imputation_quality",
                "label": "Imputation quality",
                "criterion": "0.83 <= FORMAT/SI <= 0.97; missing values retained",
                "purpose": "Evaluate the recorded score policy",
                "decision": "exclude_from_virtual_subset",
                "failed_raw": 6,
                "failed_fraction_raw": 6 / 94,
                "unique_only_raw": 4,
                "overlap_raw": 2,
            }],
        },
    }
    report = _ReportHTML(render_dataset_report(_record(), manifest))
    assert "0.83 <= FORMAT/SI <= 0.97; missing values retained" in report.text
    assert "6 / 94 (6.3830% of raw VCF records)" in report.text
    assert "virtual" in report.text.lower()
    assert "not physically filtered" in report.text.lower() or "remains unfiltered" in report.text.lower()
    assert "overlap" in report.text.lower()


def test_concordance_uses_stratum_match_counts_and_includes_p_values():
    manifest = _manifest()
    manifest["concordance_validation"] = {
        "status": "WARNING",
        "summary": {
            "variant_types": {
                "snp": {
                    "input_unique_variants": 90,
                    "vcf_unique_variants": 89,
                    "exact_matched_variants": 87,
                    "variant_union": 92,
                    "input_only_variants": 3,
                    "vcf_only_variants": 2,
                },
                "indel": {
                    "input_unique_variants": 5,
                    "vcf_unique_variants": 5,
                    "exact_matched_variants": 4,
                    "variant_union": 6,
                    "input_only_variants": 1,
                    "vcf_only_variants": 1,
                },
            },
        },
        "metrics": {"pval": {"checked": 91, "concordant": 90, "mismatches": 1}},
        "metrics_by_variant_type": {
            "snp": {"pval": {"checked": 87, "concordant": 86, "mismatches": 1}},
            "indel": {"pval": {"checked": 4, "concordant": 4, "mismatches": 0}},
        },
    }
    report = _ReportHTML(render_dataset_report(_record(), manifest))
    assert "87 / 92 (94.5652% of input/VCF variant union)" in report.text
    assert "4 / 6 (66.6667% of input/VCF variant union)" in report.text
    assert "90 / 91 (98.9011% of checked values)" in report.text
    assert "86 / 87 (98.8506% of checked values)" in report.text
    assert "pval" in report.text.lower() or "p-value" in report.text.lower()


def test_supplied_z_and_reconstructed_statistics_are_distinguished():
    record, manifest = _record(), _manifest()
    manifest["report_context"] = {
        "input_mapping": {
            "effect_column": None,
            "standard_error_column": None,
            "z_score_column": "Z_STUDY",
            "effect_allele_frequency_column": "AF_STUDY",
        },
    }
    manifest["policies"] = {"validation.se_division_floor": 0.000123}
    manifest["dataset"]["chromosome_summaries"]["1"]["stage_qc"]["effect_from_z_qc"] = {
        "status": "complete",
        "beta_computed": True,
        "se_computed": True,
        "effect_estimate_scale": "standardized",
        "effect_estimate_method_label": "Z/EAF/Neff approximation",
    }
    report = _ReportHTML(render_dataset_report(record, manifest))
    assert "Z_STUDY" in report.text and "AF_STUDY" in report.text
    assert "BETA derived: yes" in report.text
    assert "SE derived: yes" in report.text
    assert "standardized" in report.text
    assert "Z/EAF/Neff approximation" in report.text
    assert "0.000123" in report.text


def test_renderer_reuses_evidence_without_mutation_or_scientific_file_reads(monkeypatch):
    record, manifest = _record(), _manifest()
    before = deepcopy((record, manifest))

    def forbid_read(*args, **kwargs):
        pytest.fail("HTML rendering must not reopen a scientific input or report")

    monkeypatch.setattr(Path, "read_text", forbid_read)
    monkeypatch.setattr(Path, "read_bytes", forbid_read)
    first = render_dataset_report(record, manifest)
    assert first == render_dataset_report(record, manifest)
    assert (record, manifest) == before
