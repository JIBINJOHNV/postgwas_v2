"""Progress, CSV, HTML, and terminal contracts for LD clumping."""

import argparse
import csv
from pathlib import Path

import pytest
from rich.cells import cell_len

from postgwas.config import load_configuration
from postgwas.modules.ld_clumping.reporting import (
    build_ld_clumping_summary,
    render_ld_clumping_html,
    render_ld_clumping_summary,
    write_ld_clumping_html_report,
    write_ld_clumping_summary_csv,
)
from postgwas.modules.ld_clumping.service import (
    LDClumpingInputValidation,
    LDClumpingPreflight,
    LDReferenceManifest,
    run_ld_clump_direct,
)


def _summary(tmp_path: Path, *, status: str = "partial_reference") -> dict:
    module = load_configuration().modules.ld_clumping
    return build_ld_clumping_summary(
        dataset_id="STUDY",
        input_vcf=tmp_path / "study.vcf.gz",
        reference_directory=tmp_path / "reference",
        module=module,
        methods=("region", "standard"),
        status=status,
        region_result={
            "status": "completed",
            "annotated_variants": 1200,
            "ld_blocks": 300,
            "significant_blocks": 8,
            "genome_wide_significant_outside_ld_regions": 1,
            "significant_outside_ld_regions_file": (
                tmp_path / "outside_ld_regions.tsv"
            ),
            "_significant_outside_ld_region_details": {
                "columns": [
                    "CHR", "BP", "SNP", "uniq_id", "ID", "REF", "ALT",
                    "BETA", "SE", "AF", "LP", "P_value", "LDblock",
                    "exclusion_reason",
                ],
                "rows": [[
                    "3", 303, "3_303_A_G", "3_303_A_G", "rsRegionOutside",
                    "A", "G", "0.2", "0.1", "0.3", 9.0, 1e-9, None,
                    "missing_EUR_LDblock_annotation",
                ]],
                "total_rows": 1,
                "output_file": tmp_path / "outside_ld_regions.tsv",
            },
        },
        standard_result={
            "status": status,
            "reference_coverage_status": (
                "partial" if status == "partial_reference" else "complete"
            ),
            "input_significant_variants": 12,
            "significant_variants": 10,
            "skipped_significant_variants": 2,
            "independent_significant_snps": 5,
            "lead_snps": 3,
            "genomic_risk_loci": 2,
            "reference_exclusions": 2,
            "reference_exclusions_within_reported_locus_boundaries": 1,
            "reference_exclusions_outside_reported_locus_boundaries": 1,
            "reference_exclusion_reasons": {
                "missing_chromosome_reference": 1,
                "allele_mismatch": 1,
            },
            "warnings": 2,
            "other_warnings": 0,
            "skipped_chromosomes": ["X"],
            "reference_ld_rows_before_maf": 100,
            "reference_ld_rows_retained": 90,
            "low_maf_partner_rows_excluded": 10,
            "low_maf_partner_variants_excluded": 4,
            "chromosome_results": [{
                "chrom": "1",
                "status": "completed",
                "input_significant": 11,
                "significant": 10,
                "independent": 5,
                "lead": 3,
                "loci": 2,
                "warnings": 1,
                "other_warnings": 0,
                "reference_exclusions": 1,
                "reference_exclusion_reasons": {"allele_mismatch": 1},
                "reference_ld_rows_retained": 90,
            }],
            "missing_reference_details": [{
                "chromosome": "X",
                "significant_variants": 1,
                "missing_resources": ["EUR_chrX.ld.gz"],
            }],
            "_reference_exclusion_details": [
                {
                    "chromosome": "1",
                    "position": 101,
                    "canonical_id": "1_101_A_G",
                    "input_variant_id": "rsInsideExcluded",
                    "reason": "allele_mismatch",
                    "reported_locus_boundary_status": (
                        "inside_reported_locus_boundary"
                    ),
                    "overlapping_genomic_loci": "1",
                    "overlapping_locus_chromosome": "1",
                    "overlapping_locus_start": 90,
                    "overlapping_locus_end": 120,
                    "locus_membership_interpretation": (
                        "coordinate overlap only; excluded index was not used "
                        "as an LD member or locus seed"
                    ),
                },
                {
                    "chromosome": "X",
                    "position": 201,
                    "canonical_id": "X_201_C_T",
                    "input_variant_id": "rsOutsideExcluded",
                    "reason": "missing_chromosome_reference",
                    "reported_locus_boundary_status": (
                        "outside_all_reported_locus_boundaries"
                    ),
                    "overlapping_genomic_loci": "",
                    "overlapping_locus_chromosome": "",
                    "overlapping_locus_start": None,
                    "overlapping_locus_end": None,
                    "locus_membership_interpretation": (
                        "outside all reported locus boundaries"
                    ),
                },
            ],
            "_report_tables": {
                "genomic_risk_loci": {
                    "columns": [
                        "Genomic_locus", "Lead_uniqID", "CHR", "START", "END",
                    ],
                    "rows": [
                        [1, "1_101_A_G", "1", 90, 120],
                        [2, "2_202_C_T", "2", 180, 240],
                    ],
                    "total_rows": 2,
                    "output_file": tmp_path / "loci.tsv",
                },
                "lead_snps": {
                    "columns": ["Genomic_locus", "lead_SNP_id", "l_rsid"],
                    "rows": [[1, "1_101_A_G", "rsLead"]],
                    "total_rows": 1,
                    "output_file": tmp_path / "leads.tsv",
                },
                "independent_significant_snps": {
                    "columns": [
                        "Genomic_locus", "ind_sig_SNP_id", "rsID",
                        "lead_SNP_id", "r2_with_Lead",
                    ],
                    "rows": [[1, "1_101_A_G", "rsIndependent", "1_101_A_G", 1.0]],
                    "total_rows": 1,
                    "output_file": tmp_path / "independent.tsv",
                },
            },
        },
        output_paths={
            "summary_csv": tmp_path / "summary.csv",
            "html_report": tmp_path / "report.html",
            "canonical_log": tmp_path / "run.log",
        },
        threads=11,
        memory_gb=57,
    )


def test_csv_and_html_reuse_complete_partial_reference_evidence(tmp_path):
    summary = _summary(tmp_path)
    csv_path = write_ld_clumping_summary_csv(
        summary, tmp_path / "summary.csv",
    )
    html_path = write_ld_clumping_html_report(
        summary, tmp_path / "report.html",
    )

    with csv_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["record_type"] for row in rows} == {
        "overall",
        "method",
        "chromosome",
        "missing_reference",
        "reference_exclusion",
        "output",
    }
    overall = rows[0]
    assert overall["status"] == "partial_reference"
    assert overall["input_significant_variants"] == "12"
    assert overall["analysed_significant_variants"] == "10"
    assert overall["independent_significant_snps"] == "5"
    assert overall["genomic_risk_loci"] == "2"
    assert overall[
        "genome_wide_significant_outside_ld_regions"
    ] == "1"
    assert overall[
        "reference_exclusions_within_reported_locus_boundaries"
    ] == "1"
    assert overall[
        "reference_exclusions_outside_reported_locus_boundaries"
    ] == "1"
    assert overall["lead_pvalue"] == "5e-08"
    missing = next(row for row in rows if row["record_type"] == "missing_reference")
    assert missing["chromosome"] == "X"
    assert missing["skipped_significant_variants"] == "1"
    assert "EUR_chrX.ld.gz" in missing["detail_value"]

    html = html_path.read_text(encoding="utf-8")
    assert "LD clumping report" in html
    assert "PARTIAL REFERENCE" in html
    assert "not a genome-wide result" in html
    assert "Chromosome-level standard clumping" in html
    assert "Alleles do not match LD reference" in html
    assert "Chromosome LD reference unavailable" in html
    assert "GWS excluded" in html
    assert "additional independent signals or loci may be missing" in html
    assert "EUR_chrX.ld.gz" in html
    assert "5e-08" in html
    assert "Genomic Risk Loci" in html
    assert "Lead SNPs" in html
    assert "Independent Significant SNPs" in html
    assert "1_101_A_G" in html
    assert "rsLead" in html
    assert "rsIndependent" in html
    assert "Search this table" in html
    assert "data-table-filter" in html
    assert "Complete TSV" in html
    assert "Excluded indexes and reported locus boundaries" in html
    assert "rsInsideExcluded" in html
    assert "rsOutsideExcluded" in html
    assert "Inside reported locus boundary" in html
    assert "Outside all reported locus boundaries" in html
    assert "Matching locus ID" in html
    assert "Matching locus boundary" in html
    assert "chr1:90–120" in html
    assert "Not applicable" in html
    assert "Coordinate overlap does not restore reference validation" in html
    assert "Genome-wide-significant variants outside annotated LD regions" in html
    assert "rsRegionOutside" in html
    assert "does not prove biological LD independence" in html


def test_terminal_summary_uses_one_configured_value_column(tmp_path):
    summary = _summary(tmp_path)
    summary["standard"]["input_significant_variants"] = 13
    summary["standard"]["skipped_significant_variants"] = 3
    summary["standard"]["reference_exclusions"] = 3
    summary["standard"]["reference_exclusion_reasons"][
        "missing_from_reference"
    ] = 1
    summary["standard"][
        "reference_exclusions_outside_reported_locus_boundaries"
    ] = 2
    rendered = render_ld_clumping_summary(summary, label_width=42)
    rendered_text = " ".join(rendered.split())
    fields = [line for line in rendered.splitlines() if " : " in line]

    assert len(fields) >= 12
    assert len({cell_len(line.split(" : ", 1)[0]) for line in fields}) == 1
    assert "partial reference · results are not genome-wide" in rendered
    assert "partial · chrX LD reference unavailable" in rendered
    assert "GWS variants supplied" in rendered
    assert "GWS variants passing LD-reference checks" in rendered
    assert "GWS variants excluded" in rendered
    assert "Chromosome LD reference unavailable" in rendered
    assert "Position not found in LD reference" in rendered
    assert "Alleles do not match LD reference" in rendered
    assert "3 GWS variants were excluded from LD clumping" in rendered_text
    assert (
        "additional independent signals or loci may be missing"
        in rendered_text
    )
    assert "Outside all locus boundaries" in rendered
    assert "Within reported locus boundaries" in rendered
    assert "coordinate overlap only" in rendered
    assert "Excluded variants relative to reported loci" in rendered
    assert "coordinate overlap only\n\n" in rendered
    assert "Chromosomes skipped" not in rendered
    lines = rendered.splitlines()
    excluded = next(line for line in lines if "GWS variants excluded" in line)
    reason_lines = [
        line for line in lines
        if any(label in line for label in (
            "Chromosome LD reference unavailable",
            "Position not found in LD reference",
            "Alleles do not match LD reference",
        ))
    ]
    assert all(
        len(line) - len(line.lstrip())
        == len(excluded) - len(excluded.lstrip()) + 2
        for line in reason_lines
    )
    assert "Significant outside LD regions" in rendered
    assert "1 · missing EUR_LDblock annotation" in rendered
    assert "Summary CSV" in rendered
    assert str(tmp_path) not in next(
        line for line in fields if "Summary CSV" in line
    )


def test_html_escapes_untrusted_result_metadata(tmp_path):
    summary = _summary(tmp_path, status="completed")
    summary["standard"]["reference_exclusion_reasons"] = {
        "<script>alert(1)</script>": 1,
    }
    summary["standard"]["_report_tables"]["lead_snps"]["rows"] = [[
        1, "<img src=x onerror=alert(2)>", "rsLead",
    ]]

    html = render_ld_clumping_html(summary)

    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html
    assert "<img src=x onerror=alert(2)>" not in html
    assert "&lt;img src=x onerror=alert(2)&gt;" in html


def test_region_only_html_omits_standard_reference_sections(tmp_path):
    summary = _summary(tmp_path, status="completed")
    summary["methods"] = ["region"]
    summary["standard"] = {}

    html = render_ld_clumping_html(summary)

    assert "selected methods: region" in html
    assert "standard LD-reference clumping was not selected" in html
    assert "Chromosome-level standard clumping" not in html
    assert "Reference LD filtering" not in html


def _region_configuration(tmp_path: Path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf")
    configuration = load_configuration(cli_overrides={
        "modules.ld_clumping.methods": ["region"],
        "modules.ld_clumping.inputs.vcf": str(vcf),
        "modules.ld_clumping.inputs.dataset_id": "STUDY",
        "modules.ld_clumping.output_directory": str(tmp_path / "output"),
    })
    preflight = LDClumpingPreflight(
        configuration=configuration,
        vcf=vcf,
        output_directory=tmp_path / "output",
        dataset_id="STUDY",
        bcftools="bcftools",
        tabix=None,
        reference_directory=None,
        reference_manifest=None,
    )
    return configuration, preflight


def _standard_configuration(tmp_path: Path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf")
    reference = tmp_path / "reference"
    reference.mkdir()
    configuration = load_configuration(cli_overrides={
        "modules.ld_clumping.methods": ["standard"],
        "modules.ld_clumping.inputs.vcf": str(vcf),
        "modules.ld_clumping.inputs.dataset_id": "STUDY",
        "modules.ld_clumping.output_directory": str(tmp_path / "output"),
        "modules.ld_clumping.reference.directory": str(reference),
    })
    preflight = LDClumpingPreflight(
        configuration=configuration,
        vcf=vcf,
        output_directory=tmp_path / "output",
        dataset_id="STUDY",
        bcftools="bcftools",
        tabix="tabix",
        reference_directory=reference,
        reference_manifest=LDReferenceManifest(
            format_version=configuration.modules.ld_clumping.reference.format_version,
            genome_build=configuration.modules.ld_clumping.genome_build,
            populations=[configuration.modules.ld_clumping.population],
            orientation=configuration.modules.ld_clumping.reference.orientation,
            file_pattern=configuration.modules.ld_clumping.reference.file_pattern,
            reverse_file_pattern=(
                configuration.modules.ld_clumping.reference.reverse_file_pattern
            ),
            variant_inventory_pattern=(
                configuration.modules.ld_clumping.reference.variant_inventory_pattern
            ),
            columns=configuration.modules.ld_clumping.reference.columns,
            variant_inventory_columns=(
                configuration.modules.ld_clumping.reference.variant_inventory_columns
            ),
            window_kb=configuration.modules.ld_clumping.window_kb,
            minimum_r2=min(
                configuration.modules.ld_clumping.clump_r2,
                configuration.modules.ld_clumping.lead_r2,
            ),
            minimum_maf=(
                configuration.modules.ld_clumping.minimum_reference_maf
            ),
            allele_order_preserved=True,
            plink_version="PLINK v1.90-test",
        ),
    )
    return configuration, preflight


def _mock_direct_validation(monkeypatch, preflight, *, variants=3):
    input_validation = LDClumpingInputValidation(
        configuration=preflight.configuration,
        vcf=preflight.vcf,
        output_directory=preflight.output_directory,
        dataset_id=preflight.dataset_id,
        bcftools=preflight.bcftools,
        genome_build=preflight.configuration.modules.ld_clumping.genome_build.value,
        sample=preflight.dataset_id,
        variant_count=variants,
        contigs=("1",),
        required_fields=("FORMAT/LP",),
    )
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service._validate_ld_clumping_input",
        lambda *args, **kwargs: input_validation,
    )
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service._validate_ld_clumping_references",
        lambda *args, **kwargs: preflight,
    )


def test_direct_service_shows_progress_and_publishes_reports(
    monkeypatch, tmp_path, capsys,
):
    configuration, preflight = _region_configuration(tmp_path)
    _mock_direct_validation(monkeypatch, preflight, variants=100)
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.ld_clump_by_regions",
        lambda **kwargs: {
            "status": "completed",
            "annotated_variants": 100,
            "ld_blocks": 20,
            "significant_blocks": 2,
            "genome_wide_significant_outside_ld_regions": 1,
            "_significant_outside_ld_region_details": {
                "columns": ["CHR", "BP", "ID"],
                "rows": [["1", 101, "rsOutsideRegion"]],
                "total_rows": 1,
                "output_file": tmp_path / "outside.tsv",
            },
            "output_files": {},
        },
    )

    result = run_ld_clump_direct(
        argparse.Namespace(), configuration=configuration,
    )

    assert Path(result["summary_csv"]).is_file()
    assert Path(result["html_report"]).is_file()
    terminal = capsys.readouterr().out
    assert "LD clumping analysis progress" in terminal
    assert "Completed 1/4 · Validate the input GWAS-VCF" in terminal
    assert "Input summary-statistics VCF" in terminal
    assert "Total variants" in terminal
    assert "100" in terminal
    assert "Genome build inferred from VCF" in terminal
    assert "VCF embedded dataset/sample" in terminal
    assert "VCF structural validation" in terminal
    assert "Completed 2/4 · Validate selected LD-reference inputs" in terminal
    assert "Annotated-region reference" in terminal
    assert "INFO/EUR_LDblock" in terminal
    validation_output = terminal.split(
        "Completed 1/4 · Validate the input GWAS-VCF", 1
    )[1].split("Completed 3/4", 1)[0]
    validation_fields = [
        line for line in validation_output.splitlines() if " : " in line
    ]
    assert len(validation_fields) >= 7
    assert len({
        cell_len(line.split(" : ", 1)[0]) for line in validation_fields
    }) == 1
    assert "Completed 4/4 · Save CSV summary and detailed HTML report" in terminal
    assert "All 4 stages completed" in terminal
    assert "LD clumping summary" in terminal
    assert "Summary CSV" in terminal
    assert "Detailed HTML report" in terminal
    assert "Significant outside LD regions" in terminal
    assert "rsOutsideRegion" in Path(result["html_report"]).read_text(
        encoding="utf-8"
    )
    assert "_significant_outside_ld_region_details" not in result[
        "ld_clump_region"
    ]
    log_text = Path(result["canonical_log"]).read_text(encoding="utf-8")
    assert "stage_outcome" in log_text
    assert "report_source=validated_in_memory_method_results" in log_text


def test_direct_service_embeds_standard_results_without_returning_transient_rows(
    monkeypatch, tmp_path, capsys,
):
    configuration, preflight = _standard_configuration(tmp_path)
    observed = {}
    _mock_direct_validation(monkeypatch, preflight, variants=1)

    def standard(**kwargs):
        observed.update(kwargs)
        return {
            "status": "completed",
            "reference_coverage_status": "complete",
            "input_significant_variants": 1,
            "significant_variants": 1,
            "skipped_significant_variants": 0,
            "independent_significant_snps": 1,
            "lead_snps": 1,
            "genomic_risk_loci": 1,
            "reference_exclusions": 1,
            "reference_exclusions_within_reported_locus_boundaries": 0,
            "reference_exclusions_outside_reported_locus_boundaries": 1,
            "chromosome_results": [],
            "output_files": {},
            "_reference_exclusion_details": [{
                "chromosome": "2",
                "position": 202,
                "canonical_id": "2_202_C_T",
                "input_variant_id": "rsExcluded",
                "reason": "missing_from_reference",
                "reported_locus_boundary_status": (
                    "outside_all_reported_locus_boundaries"
                ),
                "overlapping_genomic_loci": "",
                "locus_membership_interpretation": (
                    "outside all reported locus boundaries"
                ),
            }],
            "_report_tables": {
                "genomic_risk_loci": {
                    "columns": ["Genomic_locus", "Lead_uniqID"],
                    "rows": [[1, "1_101_A_G"]],
                    "total_rows": 1,
                    "output_file": tmp_path / "loci.tsv",
                },
                "lead_snps": {
                    "columns": ["lead_SNP_id", "l_rsid"],
                    "rows": [["1_101_A_G", "rsLead"]],
                    "total_rows": 1,
                    "output_file": tmp_path / "leads.tsv",
                },
                "independent_significant_snps": {
                    "columns": ["ind_sig_SNP_id", "rsID"],
                    "rows": [["1_101_A_G", "rsIndependent"]],
                    "total_rows": 1,
                    "output_file": tmp_path / "independent.tsv",
                },
            },
        }

    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.ld_clump_standard", standard,
    )

    result = run_ld_clump_direct(
        argparse.Namespace(), configuration=configuration,
    )

    assert observed["include_report_tables"] is True
    terminal = capsys.readouterr().out
    compact_terminal = " ".join(terminal.split())
    assert "Completed 1/4 · Validate the input GWAS-VCF" in terminal
    assert "Completed 2/4 · Validate selected LD-reference inputs" in terminal
    assert "FUMA-style indexed LD reference" in terminal
    assert "Reference directory" in terminal
    assert "Reference format" in terminal
    assert "version 2 · upper-triangle dual index" in compact_terminal
    assert "Manifest compatibility" in terminal
    assert "Chromosome resource checks" in terminal
    assert "checked after GWS chromosomes are identified" in compact_terminal
    assert "_report_tables" not in result["ld_clump_standard"]
    assert "_reference_exclusion_details" not in result["ld_clump_standard"]
    html = Path(result["html_report"]).read_text(encoding="utf-8")
    assert "Genomic Risk Loci" in html
    assert "1_101_A_G" in html
    assert "rsLead" in html
    assert "rsIndependent" in html
    assert "rsExcluded" in html
    assert "Outside all reported locus boundaries" in html
    log_text = Path(result["canonical_log"]).read_text(encoding="utf-8")
    assert "_report_tables" not in log_text
    assert "_reference_exclusion_details" not in log_text


def test_direct_service_progress_stops_below_completion_on_method_failure(
    monkeypatch, tmp_path, capsys,
):
    configuration, preflight = _region_configuration(tmp_path)
    _mock_direct_validation(monkeypatch, preflight)

    def fail_region(**_kwargs):
        raise RuntimeError("region failed")

    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.ld_clump_by_regions",
        fail_region,
    )

    with pytest.raises(RuntimeError, match="region failed"):
        run_ld_clump_direct(
            argparse.Namespace(), configuration=configuration,
        )

    terminal = capsys.readouterr().out
    assert "Failed 3/4 · Run annotated-region LD pruning" in terminal
    assert "All 4 stages completed" not in terminal
    output = preflight.output_directory
    assert not list(output.glob("reports/*.csv"))
    assert not list(output.glob("reports/*.html"))
