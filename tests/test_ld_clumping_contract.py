import argparse
import gzip
from pathlib import Path
import subprocess

import polars as pl
import pytest
import yaml

from postgwas.config import load_configuration
from postgwas.config.models.modules.ld_clumping import (
    LDClumpingReportingConfig,
)
from postgwas.core.errors import MissingRequiredArgumentsError
from postgwas.modules.ld_clumping.cli import build_parser
from postgwas.modules.ld_clumping import ld_prune_standard
from postgwas.modules.ld_clumping.service import (
    LDClumpingError,
    LDClumpingInputValidation,
    LDClumpingPreflight,
    LDReferenceManifest,
    _validate_reference_contract,
    resolve_ld_clumping_configuration,
    run_ld_clump_direct,
    validate_ld_clumping_configuration,
)


def _vcf(path: Path) -> Path:
    path.write_bytes(b"vcf")
    return path


def _fake_tabix(monkeypatch, rows, contigs=("1",)):
    """Stub tabix for both the contig listing and the region query."""

    def run(command, *args, **kwargs):
        stdout = "\n".join(contigs) + "\n" if "--list-chroms" in command else rows
        return subprocess.CompletedProcess(command, 0, stdout=stdout, stderr="")

    monkeypatch.setattr(ld_prune_standard.subprocess, "run", run)


def _manifest(module, **overrides):
    document = {
        "format_version": module.reference.format_version,
        "genome_build": module.genome_build.value,
        "populations": [module.population.value],
        "orientation": module.reference.orientation,
        "file_pattern": module.reference.file_pattern,
        "reverse_file_pattern": module.reference.reverse_file_pattern,
        "variant_inventory_pattern": module.reference.variant_inventory_pattern,
        "columns": module.reference.columns,
        "variant_inventory_columns": module.reference.variant_inventory_columns,
        "window_kb": module.window_kb,
        "minimum_r2": min(module.clump_r2, module.lead_r2),
        "minimum_maf": module.minimum_reference_maf,
        "allele_order_preserved": True,
        "plink_version": "PLINK v1.90",
    }
    document.update(overrides)
    return document


def _validated_manifest(configuration):
    module = configuration.modules.ld_clumping
    return LDReferenceManifest.model_validate(_manifest(module))


def _mock_direct_validation(monkeypatch, preflight, *, variants=1):
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


def test_module_yaml_controls_methods_and_thresholds_without_cli_shadowing(tmp_path):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    run_config = tmp_path / "ld_clumping.yaml"
    run_config.write_text(
        "methods: [region]\n"
        "lead_pvalue: 1.0e-6\n"
        "candidate_pvalue: 0.01\n"
        "clump_r2: 0.3\n"
        "lead_r2: 0.05\n",
        encoding="utf-8",
    )
    args = build_parser().parse_args([
        "--run-config", str(run_config),
        "--vcf", str(vcf),
        "--missing-chromosome-action", "error",
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path / "output"),
    ])

    module = resolve_ld_clumping_configuration(args).modules.ld_clumping

    assert module.methods == ["region"]
    assert module.lead_pvalue == 1e-6
    assert module.candidate_pvalue == 0.01
    assert module.clump_r2 == 0.3
    assert module.lead_r2 == 0.05
    assert module.missing_chromosome_action == "error"


def test_standard_method_reports_exact_missing_ld_folder_before_external_tools(tmp_path):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--clumping-methods", "standard",
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path / "output"),
    ])

    with pytest.raises(MissingRequiredArgumentsError) as captured:
        validate_ld_clumping_configuration(args)

    assert str(captured.value) == (
        "Required argument not provided: --ld-folder. Provide --ld-folder VALUE "
        "or set modules.ld_clumping.reference.directory in the run configuration."
    )


def test_cojo_method_reports_exact_missing_reference_prefix_before_tools(tmp_path):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--clumping-methods", "cojo-slct",
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path / "output"),
    ])

    with pytest.raises(MissingRequiredArgumentsError) as captured:
        validate_ld_clumping_configuration(args)

    assert str(captured.value) == (
        "Required argument not provided: --cojo-reference-prefix. Provide "
        "--cojo-reference-prefix VALUE or set "
        "modules.gcta_cojo.reference.prefix in the run configuration."
    )


def test_cojo_cli_keeps_gcta_defaults_and_separates_locus_distance(tmp_path):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--clumping-methods", "cojo-slct",
        "--cojo-reference-prefix", str(tmp_path / "reference"),
        "--cojo-merge-dist", "300000",
        "--genome-build", "GRCh38",
        "--population", "AFR",
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path / "output"),
    ])

    configuration = resolve_ld_clumping_configuration(args)
    ld_module = configuration.modules.ld_clumping
    cojo_module = configuration.modules.gcta_cojo

    assert ld_module.cojo.merge_distance_bp == 300000
    assert ld_module.merge_distance_bp == 250000
    assert cojo_module.mode == "slct"
    assert cojo_module.analysis.window_kb == 10000
    assert cojo_module.analysis.significance_threshold == 5e-8
    assert cojo_module.analysis.collinearity_cutoff == 0.9
    assert cojo_module.analysis.frequency_difference_max == 0.2
    assert cojo_module.genome_build == "GRCh38"
    assert cojo_module.reference.population == "AFR"


def test_reference_manifest_requires_matching_build_and_dual_indexing(tmp_path):
    module = load_configuration().modules.ld_clumping
    manifest_path = tmp_path / module.reference.manifest_filename
    manifest_path.write_text(
        yaml.safe_dump(_manifest(module), sort_keys=False), encoding="utf-8",
    )

    observed = _validate_reference_contract(module, tmp_path)

    assert observed.orientation == "upper_triangle_dual_index"
    manifest_path.write_text(
        yaml.safe_dump(_manifest(module, genome_build="GRCh38"), sort_keys=False),
        encoding="utf-8",
    )
    with pytest.raises(LDClumpingError, match="genome_build=GRCh38"):
        _validate_reference_contract(module, tmp_path)
    manifest_path.write_text(
        yaml.safe_dump(
            _manifest(module, plink_version="PLINK v2.00"), sort_keys=False,
        ),
        encoding="utf-8",
    )
    with pytest.raises(LDClumpingError, match="PLINK 1.9"):
        _validate_reference_contract(module, tmp_path)

    symmetric = module.model_copy(deep=True)
    symmetric.reference.orientation = "symmetric_first_endpoint"
    manifest_path.write_text(
        yaml.safe_dump(_manifest(symmetric), sort_keys=False), encoding="utf-8",
    )
    assert (
        _validate_reference_contract(symmetric, tmp_path).orientation
        == "symmetric_first_endpoint"
    )


def test_preflight_rejects_vcf_and_ld_clumping_build_mismatch(monkeypatch, tmp_path):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--clumping-methods", "region",
        "--genome-build", "GRCh37",
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path / "output"),
    ])
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.resolve_executable",
        lambda value, *args, **kwargs: str(value),
    )
    monkeypatch.setattr(
        "postgwas.core.vcf.read_vcf_header",
        lambda *args, **kwargs: (
            "##genome_build=GRCh38\n"
            '##postgwas_version="test"\n'
            '##postgwas_dataset_id="STUDY"\n'
            '##postgwas_vcf_status="harmonised"\n'
            "##FORMAT=<ID=ES,Number=1,Type=Float,Description=\"effect\">\n"
            "##FORMAT=<ID=SE,Number=1,Type=Float,Description=\"se\">\n"
            "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"af\">\n"
            "##FORMAT=<ID=LP,Number=1,Type=Float,Description=\"lp\">\n"
            "##INFO=<ID=EUR_LDblock,Number=1,Type=String,Description=\"ld\">\n"
            "#CHROM\n"
        ),
    )

    with pytest.raises(LDClumpingError, match="GRCh38.*does not match.*GRCh37"):
        validate_ld_clumping_configuration(args)


def test_preflight_records_separate_vcf_and_standard_reference_evidence(
    monkeypatch, tmp_path,
):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    reference = tmp_path / "reference"
    reference.mkdir()
    args = build_parser().parse_args([
        "--vcf", str(vcf),
        "--clumping-methods", "standard",
        "--ld-folder", str(reference),
        "--genome-build", "GRCh37",
        "--population", "EUR",
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path / "output"),
    ])
    configuration = resolve_ld_clumping_configuration(args)
    module = configuration.modules.ld_clumping
    (reference / module.reference.manifest_filename).write_text(
        yaml.safe_dump(_manifest(module), sort_keys=False),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.resolve_executable",
        lambda value, *args, **kwargs: str(value),
    )
    monkeypatch.setattr(
        "postgwas.core.vcf.read_vcf_header",
        lambda *args, **kwargs: (
            "##fileformat=VCFv4.2\n"
            "##genome_build=GRCh37\n"
            '##postgwas_version="test"\n'
            '##postgwas_dataset_id="STUDY"\n'
            '##postgwas_vcf_status="harmonised"\n'
            "##contig=<ID=1,length=249250621>\n"
            "##FORMAT=<ID=ES,Number=1,Type=Float,Description=\"effect\">\n"
            "##FORMAT=<ID=SE,Number=1,Type=Float,Description=\"se\">\n"
            "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"af\">\n"
            "##FORMAT=<ID=LP,Number=1,Type=Float,Description=\"lp\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsummary_stats\n"
        ),
    )
    monkeypatch.setattr(
        "postgwas.core.vcf.count_indexed_vcf_records",
        lambda *args, **kwargs: 9_441_540,
    )
    monkeypatch.setattr(
        "postgwas.core.vcf.select_vcf_sample",
        lambda *args, **kwargs: "summary_stats",
    )

    preflight = validate_ld_clumping_configuration(
        args, configuration=configuration,
    )

    assert preflight.input_variant_count == 9_441_540
    assert preflight.input_genome_build == "GRCh37"
    assert preflight.input_sample == "summary_stats"
    assert preflight.input_contigs == ("1",)
    assert set(preflight.input_required_fields) == {
        "FORMAT/ES", "FORMAT/SE", "FORMAT/AF", "FORMAT/LP",
    }
    assert preflight.reference_directory == reference.resolve()
    assert preflight.reference_manifest is not None
    assert preflight.reference_manifest.orientation == "upper_triangle_dual_index"


@pytest.mark.parametrize("role", ("version", "dataset_id", "status"))
@pytest.mark.parametrize("value", (None, '""', '"unterminated', '<ID=invalid>'))
def test_direct_clumping_rejects_invalid_provenance_before_work(
    monkeypatch, tmp_path, role, value,
):
    from postgwas.modules.ld_clumping import service

    vcf = _vcf(tmp_path / "study.vcf.gz")
    output = tmp_path / "output"
    configuration = load_configuration(cli_overrides={
        "modules.ld_clumping.inputs.vcf": str(vcf),
        "modules.ld_clumping.inputs.dataset_id": "STUDY",
        "modules.ld_clumping.output_directory": str(output),
        "modules.ld_clumping.methods": ["region"],
    })
    headers = configuration.modules.formatting.input_contract.provenance_headers.model_dump()
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
        pytest.fail("Invalid provenance reached index/sample/reference/analysis work")

    monkeypatch.setattr(service, "resolve_executable", lambda value, *_args, **_kwargs: value)
    monkeypatch.setattr("postgwas.core.vcf.read_vcf_header", lambda *_args, **_kwargs: header)
    monkeypatch.setattr("postgwas.core.vcf.count_indexed_vcf_records", forbidden)
    monkeypatch.setattr("postgwas.core.vcf.select_vcf_sample", forbidden)
    monkeypatch.setattr(service, "_validate_ld_clumping_references", forbidden)
    monkeypatch.setattr(service, "ld_clump_by_regions", forbidden)

    with pytest.raises(LDClumpingError, match=name):
        run_ld_clump_direct(argparse.Namespace(), configuration=configuration)

    log_path = output / configuration.modules.ld_clumping.output_layout.canonical_log.format(
        dataset_id="STUDY", population=configuration.modules.ld_clumping.population.value,
    )
    log_text = log_path.read_text(encoding="utf-8")
    assert "postgwas_vcf_provenance" in log_text and "FAILED" in log_text
    assert name in log_text and str(vcf) in log_text
    assert {path for path in output.rglob("*") if path.is_file()} == {log_path}


def test_standard_projection_preserves_x_and_textual_variant_ids(
    monkeypatch, tmp_path,
):
    configuration = load_configuration().modules.ld_clumping.model_copy(deep=True)
    configuration.table.infer_schema_length = 1
    vcf = _vcf(tmp_path / "study.vcf.gz")

    def fake_extract(_vcf, table_path, *_args, **_kwargs):
        pl.DataFrame({
            "chrcol": ["1", "X"],
            "poscol": [100, 200],
            "neacol": ["A", "C"],
            "eacol": ["G", "T"],
            "rsIDcol": ["000123", "rsX"],
            "lpcol": ["1.0", "2.0"],
            "becol": ["0.2", "-0.3"],
            "secol": ["0.1", "0.2"],
            "eafcol": ["0.2", "0.4"],
        }).write_csv(table_path, separator="\t")
        return str(table_path)

    monkeypatch.setattr(ld_prune_standard, "extract_vcf_table", fake_extract)

    output, returned = ld_prune_standard.vcf_to_standard_ldclump(
        str(vcf),
        str(tmp_path),
        "STUDY",
        bcftools_path="bcftools",
        configuration=configuration,
    )

    observed = pl.read_csv(
        output,
        separator="\t",
        schema_overrides={"chrcol": pl.String, "rsIDcol": pl.String},
    )
    assert observed["chrcol"].to_list() == ["1", "X"]
    assert observed["rsIDcol"].to_list() == ["000123", "rsX"]
    # The frame is returned as well as written so the caller never re-parses it,
    # and the lossless LP column survives the projection.
    assert returned.to_dicts() == observed.to_dicts()
    assert "lpcol" in observed.columns


def test_selected_method_failure_is_fatal_and_changed_outputs_are_quarantined(
    monkeypatch, tmp_path,
):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    output = tmp_path / "output"
    reference = tmp_path / "reference"
    reference.mkdir()
    configuration = load_configuration(cli_overrides={
        "modules.ld_clumping.inputs.vcf": str(vcf),
        "modules.ld_clumping.inputs.dataset_id": "STUDY",
        "modules.ld_clumping.output_directory": str(output),
        "modules.ld_clumping.reference.directory": str(reference),
    })
    preflight = LDClumpingPreflight(
        configuration=configuration,
        vcf=vcf,
        output_directory=output,
        dataset_id="STUDY",
        bcftools="bcftools",
        tabix="tabix",
        reference_directory=reference,
        reference_manifest=_validated_manifest(configuration),
    )
    _mock_direct_validation(monkeypatch, preflight)

    def region(**kwargs):
        path = output / "STUDY_LDpruned_EUR.tsv"
        path.write_text("region output\n", encoding="utf-8")
        return {"ldpruned_file": str(path)}

    def standard(**kwargs):
        path = output / "STUDY_EUR_formatted.tsv"
        path.write_text("incomplete standard output\n", encoding="utf-8")
        raise RuntimeError("standard failed")

    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.ld_clump_by_regions", region,
    )
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.ld_clump_standard", standard,
    )

    with pytest.raises(RuntimeError, match="standard failed"):
        run_ld_clump_direct(argparse.Namespace(), configuration=configuration)

    assert not (output / "STUDY_LDpruned_EUR.tsv").exists()
    assert not (output / "STUDY_EUR_formatted.tsv").exists()
    assert list(output.glob("STUDY_LDpruned_EUR.tsv.partial.*"))
    assert list(output.glob("STUDY_EUR_formatted.tsv.partial.*"))


def test_validated_failure_exclusion_audit_is_not_quarantined(
    monkeypatch, tmp_path,
):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    output = tmp_path / "output"
    reference = tmp_path / "reference"
    reference.mkdir()
    configuration = load_configuration(cli_overrides={
        "modules.ld_clumping.methods": ["standard"],
        "modules.ld_clumping.inputs.vcf": str(vcf),
        "modules.ld_clumping.inputs.dataset_id": "STUDY",
        "modules.ld_clumping.output_directory": str(output),
        "modules.ld_clumping.reference.directory": str(reference),
    })
    preflight = LDClumpingPreflight(
        configuration=configuration,
        vcf=vcf,
        output_directory=output,
        dataset_id="STUDY",
        bcftools="bcftools",
        tabix="tabix",
        reference_directory=reference,
        reference_manifest=_validated_manifest(configuration),
    )
    _mock_direct_validation(monkeypatch, preflight)

    def standard(**kwargs):
        output.mkdir(parents=True, exist_ok=True)
        formatted = output / "STUDY_EUR_formatted.tsv"
        audit = output / "STUDY_EUR_LD_Reference_Exclusions.tsv"
        formatted.write_text("incomplete\n", encoding="utf-8")
        audit.write_text("canonical_id\treason\n1_100_A_G\tmissing\n", encoding="utf-8")
        error = RuntimeError("all indexes excluded")
        error.validated_failure_outputs = (str(audit),)
        raise error

    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.ld_clump_standard", standard,
    )

    with pytest.raises(RuntimeError, match="all indexes excluded"):
        run_ld_clump_direct(argparse.Namespace(), configuration=configuration)

    assert (output / "STUDY_EUR_LD_Reference_Exclusions.tsv").is_file()
    assert not (output / "STUDY_EUR_formatted.tsv").exists()
    assert list(output.glob("STUDY_EUR_formatted.tsv.partial.*"))


def test_partial_reference_status_is_preserved_in_the_canonical_log(
    monkeypatch, tmp_path,
):
    vcf = _vcf(tmp_path / "study.vcf.gz")
    output = tmp_path / "output"
    reference = tmp_path / "reference"
    reference.mkdir()
    configuration = load_configuration(cli_overrides={
        "modules.ld_clumping.inputs.vcf": str(vcf),
        "modules.ld_clumping.inputs.dataset_id": "STUDY",
        "modules.ld_clumping.output_directory": str(output),
        "modules.ld_clumping.reference.directory": str(reference),
        "modules.ld_clumping.methods": ["standard"],
    })
    preflight = LDClumpingPreflight(
        configuration=configuration,
        vcf=vcf,
        output_directory=output,
        dataset_id="STUDY",
        bcftools="bcftools",
        tabix="tabix",
        reference_directory=reference,
        reference_manifest=_validated_manifest(configuration),
    )
    _mock_direct_validation(monkeypatch, preflight)
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.ld_clump_standard",
        lambda **kwargs: {
            "status": "partial_reference",
            "ldpruned_sig_file": None,
            "skipped_chromosomes": ["X"],
            "skipped_significant_variants": 112,
        },
    )

    result = run_ld_clump_direct(
        argparse.Namespace(), configuration=configuration,
    )

    assert result["status"] == "partial_reference"
    log_path = Path(result["canonical_log"])
    log_text = log_path.read_text(encoding="utf-8")
    assert "WARNING  ld_clumping_partial_reference" in log_text
    assert "skipped_significant_variants=112" in log_text
    assert "DONE     ld_clumping_completed" not in log_text


def _reference_data(index_id, rows=None, exclusions=()):
    rows = rows or []
    partners = (
        pl.DataFrame(
            rows,
            schema=["uniq_id", "ld_id", "ref_chr", "ref_pos", "r2"],
            orient="row",
        )
        if rows
        else ld_prune_standard._empty_ld_partners()
    )
    return ld_prune_standard.ChromosomeReferenceData(
        partners={index_id: partners} if rows else {},
        verified_indexes=frozenset([index_id]),
        exclusions=tuple(exclusions),
        ld_rows=len(rows),
    )


def test_verified_index_assigns_batched_partner():
    index_id = "1_100_A_G"
    partner_id = "1_200_C_T"
    reference_data = _reference_data(
        index_id,
        [(partner_id, partner_id, "1", 200, 0.8)],
    )

    partners = ld_prune_standard.get_ld_partners(
        "1", 100, index_id, "EUR_chr1.ld.gz", 0.6,
        window_kb=250, reference_data=reference_data,
    )

    assert partners["uniq_id"].to_list() == [index_id, partner_id]


def test_fully_resolved_ld_query_does_not_reload_configuration(monkeypatch):
    monkeypatch.setattr(
        ld_prune_standard,
        "load_configuration",
        lambda: (_ for _ in ()).throw(AssertionError("configuration reloaded")),
    )
    reference_data = _reference_data(
        "1_100_A_G",
        [("1_200_C_T", "1_200_C_T", "1", 200, 0.8)],
    )

    partners = ld_prune_standard.get_ld_partners(
        "1", 100, "1_100_A_G", "EUR_chr1.ld.gz", 0.6,
        window_kb=250, reference_data=reference_data,
    )

    assert partners["uniq_id"].to_list() == ["1_100_A_G", "1_200_C_T"]


def test_unverified_index_cannot_be_treated_as_independent():
    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match="excluded index variant reached LD clumping",
    ):
        ld_prune_standard.get_ld_partners(
            "1", 100, "1_100_A_G", "EUR_chr1.ld.gz", 0.6,
            window_kb=250,
            reference_data=ld_prune_standard.ChromosomeReferenceData(
                partners={},
                verified_indexes=frozenset(),
                exclusions=(),
                ld_rows=0,
            ),
        )


def test_configured_window_excludes_distant_ld_pairs():
    index_id = "1_100_A_G"
    distant_id = "1_400000_C_T"
    reference_data = _reference_data(
        index_id,
        [(distant_id, distant_id, "1", 400000, 0.9)],
    )

    partners = ld_prune_standard.get_ld_partners(
        "1", 100, index_id, "EUR_chr1.ld.gz", 0.6,
        window_kb=250, reference_data=reference_data,
    )

    assert partners["uniq_id"].to_list() == [index_id]


def test_window_reachable_annotations_preserve_fuma_membership_semantics():
    """The speed subset must not change candidate/reference classification."""
    index_id = "1_100_A_G"
    candidate_id = "1_200_C_T"
    high_p_id = "1_300_G_T"
    reference_only_id = "1_400_A_T"
    boundary_id = "1_1100_A_C"
    outside_id = "1_1101_C_G"
    gwas = ld_prune_standard.add_canonical_ids(pl.DataFrame({
        "chrcol": ["1"] * 5,
        "poscol": [100, 200, 300, 1100, 1101],
        "neacol": ["A", "C", "G", "A", "C"],
        "eacol": ["G", "T", "T", "C", "G"],
        "rsIDcol": ["rsIndex", "rsCandidate", "rsHighP", "rsBoundary", "rsOutside"],
        "pcol": [1e-9, 0.01, 0.5, 0.02, 0.01],
        "becol": [0.2, 0.1, -0.1, 0.3, -0.2],
        "secol": [0.1] * 5,
        "eafcol": [0.2, 0.3, 0.4, 0.25, 0.35],
    }))
    reference_data = _reference_data(
        index_id,
        [
            (candidate_id, candidate_id, "1", 200, 0.9),
            (high_p_id, high_p_id, "1", 300, 0.8),
            (reference_only_id, reference_only_id, "1", 400, 0.7),
            (boundary_id, boundary_id, "1", 1100, 0.65),
            (outside_id, outside_id, "1", 1101, 0.95),
        ],
    )
    verified = gwas.filter(pl.col("uniq_id") == index_id)
    reachable = ld_prune_standard._window_reachable_gwas_annotations(
        gwas, verified, 1,
    )

    assert boundary_id in reachable["uniq_id"].to_list()
    assert outside_id not in reachable["uniq_id"].to_list()

    arguments = {
        "ld_path": "EUR_chr1.ld.gz",
        "lead_p_threshold": 5e-8,
        "r2_clump_threshold": 0.6,
        "candidate_p_threshold": 0.05,
        "window_kb": 1,
        "reference_data": reference_data,
        "log": lambda _message: None,
    }
    full_members = ld_prune_standard.find_ind_sig_snps(gwas, **arguments)
    reachable_members = ld_prune_standard.find_ind_sig_snps(
        reachable, **arguments,
    )

    assert reachable_members.to_dicts() == full_members.to_dicts()
    assert reachable_members["uniq_id"].to_list() == [
        index_id,
        candidate_id,
        reference_only_id,
        boundary_id,
    ]
    reference_only = reachable_members.filter(
        pl.col("uniq_id") == reference_only_id
    ).row(0, named=True)
    assert reference_only["is_gwas_tagged"] is False
    assert reference_only["pcol"] is None
    assert high_p_id not in reachable_members["uniq_id"].to_list()


def test_no_significant_chromosome_does_not_require_an_ld_file(tmp_path):
    gwas = ld_prune_standard.add_canonical_ids(pl.DataFrame({
        "chrcol": ["X"],
        "poscol": [100],
        "neacol": ["A"],
        "eacol": ["G"],
        "rsIDcol": ["rs1"],
        "pcol": [0.5],
        "becol": [-0.2],
        "secol": [0.1],
        "eafcol": [0.3],
    }))

    result = ld_prune_standard.process_chromosome(
        "X", gwas, tmp_path, "EUR", 5e-8, 0.6, 0.1, 250000,
    )

    assert result["progress"]["status"] == "no_significant_variants"


def test_canonical_matching_id_does_not_flip_alleles_or_effects():
    result = ld_prune_standard.add_canonical_ids(pl.DataFrame({
        "chrcol": ["chr01"],
        "poscol": [100],
        "neacol": ["G"],
        "eacol": ["A"],
        "rsIDcol": ["rs1"],
        "pcol": [1e-9],
        "becol": [-0.2],
        "secol": [0.1],
        "eafcol": [0.3],
    }))

    row = result.row(0, named=True)
    assert row["uniq_id"] == "1_100_A_G"
    assert row["chrcol"] == "1"
    assert (row["neacol"], row["eacol"], row["becol"], row["eafcol"]) == (
        "G", "A", -0.2, 0.3,
    )


def test_conflicting_effect_allele_orientation_is_not_silently_collapsed():
    gwas = pl.DataFrame({
        "chrcol": ["1", "chr1"],
        "poscol": [100, 100],
        "neacol": ["G", "A"],
        "eacol": ["A", "G"],
        "rsIDcol": ["rs1", "rs1"],
        "pcol": [1e-9, 1e-9],
        "becol": [-0.2, 0.2],
        "secol": [0.1, 0.1],
        "eafcol": [0.3, 0.7],
    })

    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match="conflicting effect-allele orientation.*flipping is prohibited",
    ):
        ld_prune_standard.add_canonical_ids(gwas)


def test_reference_preparation_preserves_allele_order_and_supports_orientations():
    script = Path("tools/resource_preparation/ld_file_preparation.sh")
    text = script.read_text(encoding="utf-8")

    assert "--keep-allele-order" in text
    assert "upper_triangle_dual_index" in text
    assert "symmetric_first_endpoint" in text
    assert "print $4, $5, $6, $1, $2, $3, $7" in text
    assert "variant_inventory_columns" in text
    subprocess.run(["bash", "-n", str(script)], check=True)


def test_reference_preparation_rejects_invalid_ranges_before_plink(tmp_path):
    script = Path("tools/resource_preparation/ld_file_preparation.sh")
    fileset = tmp_path / "reference"
    for suffix in ("bed", "bim", "fam"):
        (tmp_path / f"reference.{suffix}").write_text("input\n", encoding="utf-8")

    completed = subprocess.run(
        [
            "bash", str(script),
            "--bfile", str(fileset),
            "--output-dir", str(tmp_path / "output"),
            "--population", "EUR",
            "--genome-build", "GRCh37",
            "--chromosomes", "1",
            "--plink", "printf",
            "--gzip", "gzip",
            "--bgzip", "printf",
            "--tabix", "printf",
            "--window-kb", "250",
            "--ld-window-variants", "99999",
            "--minimum-r2", "1.1",
            "--minimum-maf", "0.00001",
            "--threads", "8",
            "--chromosome-workers", "2",
            "--orientation", "upper_triangle_dual_index",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 2
    assert "Invalid --minimum-r2" in completed.stderr
    assert not (tmp_path / "output").exists()


def test_reference_preparation_executes_the_exact_allele_safe_plink_contract(
    tmp_path,
):
    script = Path("tools/resource_preparation/ld_file_preparation.sh")
    fileset = tmp_path / "reference"
    (tmp_path / "reference.bed").write_text("input\n", encoding="utf-8")
    (tmp_path / "reference.fam").write_text("input\n", encoding="utf-8")
    (tmp_path / "reference.bim").write_text(
        "1\t1_100_A_G\t0\t100\tA\tG\n"
        "1\t1_200_C_T\t0\t200\tC\tT\n"
        "2\t2_300_A_C\t0\t300\tA\tC\n"
        "2\t2_400_G_T\t0\t400\tG\tT\n",
        encoding="utf-8",
    )

    fake_plink = tmp_path / "plink"
    fake_plink.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "if [[ ${1:-} == --version ]]; then echo 'PLINK v1.90b-test'; exit 0; fi\n"
        "out=''\n"
        "chrom=''\n"
        "args=(\"$@\")\n"
        "for ((i=0; i < ${#args[@]}; i++)); do\n"
        "  if [[ ${args[$i]} == --out ]]; then out=${args[$((i + 1))]}; fi\n"
        "  if [[ ${args[$i]} == --chr ]]; then chrom=${args[$((i + 1))]}; fi\n"
        "done\n"
        "if printf '%s\\n' \"$@\" | grep -q -- --freq; then\n"
        "  printf '%s\\n' \"$@\" > \"${out}.args\"\n"
        "  printf ' CHR SNP A1 A2 MAF NCHROBS\\n 1 1_100_A_G A G 0.2 1006\\n 1 1_200_C_T C T 0.8 1006\\n 2 2_300_A_C A C 0.2 1006\\n 2 2_400_G_T G T 0.3 1006\\n' > \"${out}.frq\"\n"
        "  exit 0\n"
        "fi\n"
        "printf '%s\\n' \"$@\" > \"${out}.args\"\n"
        "if [[ $chrom == 1 ]]; then\n"
        "  printf 'CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\\n1 100 1_100_A_G 1 200 1_200_C_T 0.8\\n' | gzip -c > \"${out}.ld.gz\"\n"
        "else\n"
        "  printf 'CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\\n2 300 2_300_A_C 2 400 2_400_G_T 0.7\\n' | gzip -c > \"${out}.ld.gz\"\n"
        "fi\n",
        encoding="utf-8",
    )
    fake_bgzip = tmp_path / "bgzip"
    fake_bgzip.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\ngzip -c\n",
        encoding="utf-8",
    )
    fake_tabix = tmp_path / "tabix"
    fake_tabix.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\nlast=${!#}\nprintf 'index\\n' > \"${last}.tbi\"\n",
        encoding="utf-8",
    )
    for executable in (fake_plink, fake_bgzip, fake_tabix):
        executable.chmod(0o755)

    output = tmp_path / "output"
    command = [
        "bash", str(script),
        "--bfile", str(fileset),
        "--output-dir", str(output),
        "--population", "EUR",
        "--genome-build", "GRCh37",
        "--chromosomes", "1,2",
        "--plink", str(fake_plink),
        "--gzip", "gzip",
        "--bgzip", str(fake_bgzip),
        "--tabix", str(fake_tabix),
        "--window-kb", "250",
        "--ld-window-variants", "99999",
        "--minimum-r2", "0.05",
        "--minimum-maf", "0.00001",
        "--threads", "8",
        "--chromosome-workers", "2",
        "--orientation", "upper_triangle_dual_index",
    ]
    subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )

    arguments = (output / "EUR_chr1.args").read_text(encoding="utf-8").splitlines()
    assert arguments == [
        "--bfile", str(fileset),
        "--chr", "1",
        "--keep-allele-order",
        "--r2", "gz",
        "--ld-window", "99999",
        "--ld-window-kb", "250",
        "--ld-window-r2", "0.05",
        "--maf", "0.00001",
        "--threads", "4",
        "--out", str(output / "EUR_chr1"),
    ]
    inventory_arguments = (
        output / ".EUR.inventory_working.args"
    ).read_text(encoding="utf-8").splitlines()
    assert inventory_arguments == [
        "--bfile", str(fileset),
        "--keep-allele-order",
        "--maf", "0.00001",
        "--freq",
        "--threads", "8",
        "--out", str(output / ".EUR.inventory_working"),
    ]
    with gzip.open(output / "EUR_chr1.ld.gz", "rt", encoding="utf-8") as handle:
        forward_pairs = [line.rstrip("\n").split("\t") for line in handle]
    with gzip.open(
        output / "EUR_chr1.reverse.ld.gz", "rt", encoding="utf-8"
    ) as handle:
        reverse_pairs = [line.rstrip("\n").split("\t") for line in handle]
    assert forward_pairs == [
        ["1", "100", "1_100_A_G", "1", "200", "1_200_C_T", "0.8"],
    ]
    assert reverse_pairs == [
        ["1", "200", "1_200_C_T", "1", "100", "1_100_A_G", "0.8"],
    ]
    with gzip.open(
        output / "EUR_chr1.variants.tsv.gz", "rt", encoding="utf-8"
    ) as handle:
        inventory = [line.rstrip("\n").split("\t") for line in handle]
    assert inventory == [
        ["1", "100", "1_100_A_G", "1_100_A_G", "A", "G", "0.2"],
        ["1", "200", "1_200_C_T", "1_200_C_T", "C", "T", "0.2"],
    ]
    manifest = yaml.safe_load(
        (output / "ld_reference.yaml").read_text(encoding="utf-8")
    )
    assert manifest["format_version"] == 2
    assert manifest["orientation"] == "upper_triangle_dual_index"
    assert manifest["minimum_maf"] == 0.00001
    assert manifest["allele_order_preserved"] is True
    assert manifest["plink_version"] == "PLINK v1.90b-test"
    assert (output / "EUR_chr1.ld.gz.tbi").is_file()
    assert (output / "EUR_chr1.reverse.ld.gz.tbi").is_file()
    assert (output / "EUR_chr1.variants.tsv.gz.tbi").is_file()
    assert (output / "EUR_chr2.reverse.ld.gz.tbi").is_file()
    assert (output / "EUR_chr2.variants.tsv.gz.tbi").is_file()

    forward_before_migration = (output / "EUR_chr1.ld.gz").read_bytes()
    subprocess.run(
        [*command, "--reuse-existing-ld"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert (output / "EUR_chr1.ld.gz").read_bytes() == forward_before_migration


def test_block_lead_is_chosen_by_lp_when_p_underflows(monkeypatch, tmp_path):
    """Two variants in one block, both stronger than 1e-308.

    Recovering P as 10**-LP gives both exactly 0.0. Ordering on P therefore
    ties and falls back to the SNP string, which would pick the variant whose
    identifier sorts first rather than the stronger association.
    """
    from postgwas.modules.ld_clumping import ld_prune_region

    configuration = load_configuration(cli_overrides={
        "modules.ld_clumping.remove_mhc": False,
    }).modules.ld_clumping
    weaker, stronger = "1_100_A_G", "1_200_C_T"

    def fake_extract(vcf, table_path, dataset_id, columns, bcftools, **kwargs):
        pl.DataFrame({
            "CHR": ["1", "1"],
            "BP": ["100", "200"],
            "ID": ["rsWeaker", "rsStronger"],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
            "BETA": ["0.2", "0.3"],
            "SE": ["0.1", "0.1"],
            "AF": ["0.3", "0.4"],
            # 10 ** -400 and 10 ** -500 are both exactly 0.0 as floats.
            "LP": ["400.0", "500.0"],
            "EUR_LDblock": ["EUR_1_1_1000", "EUR_1_1_1000"],
        }).write_csv(table_path, separator="\t")
        return str(table_path)

    monkeypatch.setattr(ld_prune_region, "extract_vcf_table", fake_extract)

    result = ld_prune_region.ld_clump_by_regions(
        str(_vcf(tmp_path / "study.vcf.gz")),
        str(tmp_path / "out"),
        "STUDY",
        population="EUR",
        bcftools="bcftools",
        configuration=configuration,
    )

    pruned = pl.read_csv(result["ldpruned_file"], separator="\t")
    assert pruned.height == 1
    # Both recovered P values underflow, so P alone cannot discriminate.
    assert pruned["P_value"].to_list() == [0.0]
    assert pruned["LP"].to_list() == [500.0]
    assert pruned["SNP"].to_list() == [stronger]
    assert stronger > weaker, "the stronger variant must not also sort first"
    # The allele-order-independent key lets region output join standard output.
    assert pruned["uniq_id"].to_list() == ["1_200_C_T"]


def test_locus_summary_reports_lp_without_the_underflow_fallback():
    """A lead SNP stronger than 1e-308 must keep its real LP in the summary."""
    ind_sig = pl.DataFrame({
        "uniq_id": ["1_100_A_G"],
        "ind_sig_SNP_id": ["1_100_A_G"],
        "chrcol": ["1"],
        "poscol": [100],
        "pcol": [0.0],
        "lpcol": [412.5],
        "becol": [0.2],
        "secol": [0.1],
        "eafcol": [0.3],
        "rsIDcol": ["rs1"],
        "eacol": ["G"],
        "neacol": ["A"],
        "r2_with_IndSig": [1.0],
        "is_gwas_tagged": [True],
    })
    leads = pl.DataFrame({
        "uniq_id": ["1_100_A_G"],
        "lead_SNP_id": ["1_100_A_G"],
        "r2_with_Lead": [1.0],
    })

    summary, _, _, _, _ = ld_prune_standard.define_genomic_risk_loci(
        ind_sig, leads, 250000, 1, summary_pvalue_thresholds={"n_5e_8": 5e-8},
    )

    assert summary["LP"].to_list() == [412.5]
    assert summary["P_value"].to_list() == [0.0]


def test_lead_clumping_preserves_fuma_greedy_order_with_set_membership():
    variants = ["1_100_A_G", "1_200_C_T", "1_300_A_C"]
    heads = pl.DataFrame({
        "uniq_id": variants,
        "ind_sig_SNP_id": variants,
        "chrcol": ["1", "1", "1"],
        "poscol": [100, 200, 300],
        "pcol": [1e-12, 1e-10, 1e-9],
        "lpcol": [12.0, 10.0, 9.0],
    })
    reference_data = ld_prune_standard.ChromosomeReferenceData(
        partners={
            variants[0]: pl.DataFrame({
                "uniq_id": [variants[1], variants[2]],
                "ld_id": ["raw_b", "raw_c"],
                "ref_chr": ["1", "1"],
                "ref_pos": [200, 300],
                "r2": [0.2, 0.05],
            }),
        },
        verified_indexes=frozenset(variants),
        exclusions=(),
        ld_rows=2,
    )

    leads = ld_prune_standard.find_lead_snps(
        heads,
        "EUR_chr1.ld.gz",
        0.1,
        window_kb=250,
        reference_data=reference_data,
        log=lambda _message: None,
    )

    assert leads.select("uniq_id", "lead_SNP_id", "r2_with_Lead").rows() == [
        (variants[0], variants[0], 1.0),
        (variants[1], variants[0], 0.2),
        (variants[2], variants[2], 1.0),
    ]


def test_chr_prefixed_reference_is_queried_with_its_own_contig_label(monkeypatch):
    """A chr-prefixed reference must work for batched position retrieval."""
    observed_regions = []

    def run(command, *args, **kwargs):
        if "--list-chroms" in command:
            return subprocess.CompletedProcess(command, 0, stdout="chr1\n", stderr="")
        regions = Path(command[2]).read_text(encoding="utf-8")
        observed_regions.append(regions)
        return subprocess.CompletedProcess(command, 0, stdout="payload\n", stderr="")

    monkeypatch.setattr(ld_prune_standard.subprocess, "run", run)

    payload = ld_prune_standard._batch_tabix_query(
        "tabix", "EUR_chr1.ld.gz", "1", [200, 100, 100], contig_cache={},
    )

    assert payload == "payload\n"
    assert observed_regions == ["chr1\t100\t100\nchr1\t200\t200\n"]


def test_batched_tabix_table_streams_stdout_to_a_temporary_file(monkeypatch):
    query_used_file = []

    def run(command, *args, **kwargs):
        if "--list-chroms" in command:
            return subprocess.CompletedProcess(
                command, 0, stdout="1\n", stderr=""
            )
        query_used_file.append("stdout" in kwargs and "capture_output" not in kwargs)
        kwargs["stdout"].write(
            "1\t100\traw_index\t1\t200\traw_partner\t0.8\n"
        )
        return subprocess.CompletedProcess(command, 0, stdout=None, stderr="")

    monkeypatch.setattr(ld_prune_standard.subprocess, "run", run)

    table = ld_prune_standard._batch_tabix_table(
        "tabix",
        "EUR_chr1.ld.gz",
        "1",
        [100],
        contig_cache={},
        parser=ld_prune_standard._parse_ld_payload,
    )

    assert query_used_file == [True]
    assert table.select("snp_a", "snp_b", "r2").rows() == [
        ("raw_index", "raw_partner", 0.8)
    ]


def test_absent_reference_contig_names_the_mismatch(monkeypatch):
    """A genuine contig mismatch must not be reported as a missing index SNP."""
    _fake_tabix(monkeypatch, "", contigs=("1", "2"))

    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match="no contig matching chromosome 7.*different contig naming",
    ):
        ld_prune_standard._batch_tabix_query(
            "tabix", "EUR_chr7.ld.gz", "7", [100], contig_cache={},
        )


def test_upper_triangle_reference_batches_both_endpoints_and_skips_mismatch(
    monkeypatch,
):
    index_id = "1_100_A_G"
    reverse_partner = "1_50_C_T"
    forward_partner = "1_200_A_C"
    calls = []
    payloads = iter([
        # First inventory lookup validates index identity at positions 100/300.
        "1\t100\traw_index\t1_100_A_G\tA\tG\t0.2\n"
        "1\t300\traw_other\t1_300_A_T\tA\tT\t0.3\n",
        # Forward pairs: index is PLINK endpoint A.
        "1\t100\traw_index\t1\t200\traw_forward\t0.8\n",
        # Reverse sidecar: original endpoint B is now endpoint A.
        "1\t100\traw_index\t1\t50\traw_reverse\t0.9\n",
        # Inventory rows for every pair endpoint.
        "1\t50\traw_reverse\t1_50_C_T\tC\tT\t0.2\n"
        "1\t100\traw_index\t1_100_A_G\tA\tG\t0.2\n"
        "1\t200\traw_forward\t1_200_A_C\tA\tC\t0.3\n",
    ])

    def query(_tabix, path, _chrom, positions, **_kwargs):
        calls.append((str(path), tuple(sorted(set(positions)))))
        return next(payloads)

    monkeypatch.setattr(ld_prune_standard, "_batch_tabix_query", query)
    significant = pl.DataFrame({
        "chrcol": ["1", "1"],
        "poscol": [100, 300],
        "uniq_id": [index_id, "1_300_C_G"],
        "rsIDcol": ["rsIndex", "rsMismatch"],
    })

    reference_data = ld_prune_standard.prepare_chromosome_reference(
        significant,
        "1",
        "EUR_chr1.ld.gz",
        "EUR_chr1.variants.tsv.gz",
        reverse_ld_path="EUR_chr1.reverse.ld.gz",
        orientation="upper_triangle_dual_index",
        tabix_bin="tabix",
        minimum_reference_maf=0.00001,
        missing_index_action="warning_skip",
        log=lambda _message: None,
    )

    assert reference_data.verified_indexes == frozenset([index_id])
    assert [record["reason"] for record in reference_data.exclusions] == [
        "allele_mismatch"
    ]
    partners = ld_prune_standard.get_ld_partners(
        "1", 100, index_id, "EUR_chr1.ld.gz", 0.6,
        window_kb=250, reference_data=reference_data,
    )
    assert partners["uniq_id"].to_list() == [
        index_id,
        reverse_partner,
        forward_partner,
    ]
    assert len(calls) == 4


def test_reference_maf_filters_every_ld_partner_and_reports_counts(monkeypatch):
    payloads = iter([
        "1\t100\traw_index\t1_100_A_G\tA\tG\t0.2\n",
        "1\t100\traw_index\t1\t150\traw_low\t0.9\n"
        "1\t100\traw_index\t1\t200\traw_boundary\t0.8\n"
        "1\t100\traw_index\t1\t250\traw_high\t0.7\n",
        "",
        "1\t100\traw_index\t1_100_A_G\tA\tG\t0.2\n"
        "1\t150\traw_low\t1_150_C_T\tC\tT\t0.005\n"
        "1\t200\traw_boundary\t1_200_A_C\tA\tC\t0.01\n"
        "1\t250\traw_high\t1_250_G_T\tG\tT\t0.2\n",
    ])
    monkeypatch.setattr(
        ld_prune_standard,
        "_batch_tabix_query",
        lambda *_args, **_kwargs: next(payloads),
    )
    messages = []
    significant = pl.DataFrame({
        "chrcol": ["1"],
        "poscol": [100],
        "uniq_id": ["1_100_A_G"],
        "rsIDcol": ["rsIndex"],
    })

    reference_data = ld_prune_standard.prepare_chromosome_reference(
        significant,
        "1",
        "EUR_chr1.ld.gz",
        "EUR_chr1.variants.tsv.gz",
        reverse_ld_path="EUR_chr1.reverse.ld.gz",
        orientation="upper_triangle_dual_index",
        tabix_bin="tabix",
        minimum_reference_maf=0.01,
        missing_index_action="warning_skip",
        log=messages.append,
    )
    partners = ld_prune_standard.get_ld_partners(
        "1",
        100,
        "1_100_A_G",
        "EUR_chr1.ld.gz",
        0.6,
        window_kb=250,
        reference_data=reference_data,
    )

    assert partners["uniq_id"].to_list() == [
        "1_100_A_G",
        "1_200_A_C",
        "1_250_G_T",
    ]
    assert reference_data.ld_rows_before_maf == 3
    assert reference_data.ld_rows == 2
    assert reference_data.low_maf_partner_rows_excluded == 1
    assert reference_data.low_maf_partner_variants_excluded == 1
    assert any(
        "minimum_reference_maf=0.01" in message
        and "low_maf_partner_rows_excluded=1" in message
        for message in messages
    )


def test_verified_index_without_pair_rows_is_retained_as_self(monkeypatch):
    payloads = iter([
        "1\t100\traw_index\t1_100_A_G\tA\tG\t0.2\n",
        "",
        "",
    ])
    monkeypatch.setattr(
        ld_prune_standard,
        "_batch_tabix_query",
        lambda *_args, **_kwargs: next(payloads),
    )
    significant = pl.DataFrame({
        "chrcol": ["1"],
        "poscol": [100],
        "uniq_id": ["1_100_A_G"],
        "rsIDcol": ["rsIndex"],
    })
    reference_data = ld_prune_standard.prepare_chromosome_reference(
        significant,
        "1",
        "EUR_chr1.ld.gz",
        "EUR_chr1.variants.tsv.gz",
        reverse_ld_path="EUR_chr1.reverse.ld.gz",
        orientation="upper_triangle_dual_index",
        tabix_bin="tabix",
        minimum_reference_maf=0.00001,
        missing_index_action="warning_skip",
        log=lambda _message: None,
    )

    partners = ld_prune_standard.get_ld_partners(
        "1", 100, "1_100_A_G", "EUR_chr1.ld.gz", 0.6,
        window_kb=250, reference_data=reference_data,
    )
    assert partners.select("uniq_id", "r2").rows() == [("1_100_A_G", 1.0)]


def test_missing_and_low_maf_indexes_are_audited_and_strict_mode_fails(
    monkeypatch,
):
    inventory = "1\t200\traw_low\t1_200_C_T\tC\tT\t0.000001\n"
    monkeypatch.setattr(
        ld_prune_standard,
        "_batch_tabix_query",
        lambda *_args, **_kwargs: inventory,
    )
    significant = pl.DataFrame({
        "chrcol": ["1", "1"],
        "poscol": [100, 200],
        "uniq_id": ["1_100_A_G", "1_200_C_T"],
        "rsIDcol": ["rsMissing", "rsLow"],
    })

    reference_data = ld_prune_standard.prepare_chromosome_reference(
        significant,
        "1",
        "EUR_chr1.ld.gz",
        "EUR_chr1.variants.tsv.gz",
        reverse_ld_path="EUR_chr1.reverse.ld.gz",
        orientation="upper_triangle_dual_index",
        tabix_bin="tabix",
        minimum_reference_maf=0.00001,
        missing_index_action="warning_skip",
        log=lambda _message: None,
    )

    assert [record["reason"] for record in reference_data.exclusions] == [
        "missing_from_reference",
        "below_reference_maf",
    ]
    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match="Significant index variants are absent.*action=error",
    ):
        ld_prune_standard.prepare_chromosome_reference(
            significant,
            "1",
            "EUR_chr1.ld.gz",
            "EUR_chr1.variants.tsv.gz",
            reverse_ld_path="EUR_chr1.reverse.ld.gz",
            orientation="upper_triangle_dual_index",
            tabix_bin="tabix",
            minimum_reference_maf=0.00001,
            missing_index_action="error",
            log=lambda _message: None,
        )


def test_malformed_inventory_and_ld_rows_fail_the_reference_contract():
    malformed_inventory = pl.DataFrame({
        "chromosome": ["1"],
        "position": [100],
        "reference_id": ["raw"],
        "canonical_id": ["1_100_A_T"],
        "allele_1": ["A"],
        "allele_2": ["G"],
        "minor_allele_frequency": [0.2],
    })
    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match="canonical allele-ID contract",
    ):
        ld_prune_standard._validate_inventory(
            malformed_inventory,
            chromosome="1",
            reference_file="EUR_chr1.variants.tsv.gz",
        )

    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match="coordinate, or r2 contract",
    ):
        ld_prune_standard._parse_ld_payload(
            "1\t100\tindex\t1\t200\t.\tnan\n",
            reference_file="EUR_chr1.ld.gz",
            chromosome="1",
        )


def test_vectorised_inventory_validation_normalises_chromosome_and_allele_order():
    inventory = pl.DataFrame({
        "chromosome": ["chr01", "1"],
        "position": [100, 200],
        "reference_id": ["raw_reversed", "raw_boundary"],
        "canonical_id": ["1_100_A_G", "1_200_C_T"],
        "allele_1": ["G", "C"],
        "allele_2": ["A", "T"],
        "minor_allele_frequency": [0.2, 0.5],
    })

    observed = ld_prune_standard._validate_inventory(
        inventory,
        chromosome="1",
        reference_file="EUR_chr1.variants.tsv.gz",
    )

    assert observed is None


def test_excluded_indexes_are_classified_against_final_locus_boundaries():
    exclusions = [
        {
            "chromosome": chromosome,
            "position": position,
            "canonical_id": canonical_id,
            "input_variant_id": input_id,
            "reason": reason,
            "observed_reference_ids": "",
            "action": "warning_skip",
        }
        for chromosome, position, canonical_id, input_id, reason in (
            ("1", 100, "1_100_A_G", "rsInside", "missing_from_reference"),
            ("1", 250, "1_250_C_T", "rsOutside", "allele_mismatch"),
            ("X", 150, "X_150_A_C", "rsMissingChr", "missing_chromosome_reference"),
            ("2", 500, "2_500_G_T", "rsEnd", "below_reference_maf"),
        )
    ]
    loci = pl.DataFrame({
        "Genomic_locus": [7, 8],
        "CHR": ["1", "2"],
        "START": [100, 400],
        "END": [200, 500],
    })

    classified = ld_prune_standard._classify_reference_exclusions_by_locus(
        exclusions,
        loci,
    )

    assert [row["canonical_id"] for row in classified] == [
        "1_100_A_G", "1_250_C_T", "X_150_A_C", "2_500_G_T",
    ]
    assert [row["reported_locus_boundary_status"] for row in classified] == [
        "inside_reported_locus_boundary",
        "outside_all_reported_locus_boundaries",
        "outside_all_reported_locus_boundaries",
        "inside_reported_locus_boundary",
    ]
    assert [row["overlapping_genomic_loci"] for row in classified] == [
        "7", "", "", "8",
    ]
    assert [row["overlapping_locus_chromosome"] for row in classified] == [
        "1", "", "", "2",
    ]
    assert [row["overlapping_locus_start"] for row in classified] == [
        100, None, None, 400,
    ]
    assert [row["overlapping_locus_end"] for row in classified] == [
        200, None, None, 500,
    ]
    assert "not used as an LD member or locus seed" in classified[0][
        "locus_membership_interpretation"
    ]
    assert ld_prune_standard._reference_exclusion_locus_counts(classified) == {
        "reference_exclusions_within_reported_locus_boundaries": 2,
        "reference_exclusions_outside_reported_locus_boundaries": 2,
    }


def test_reference_exclusion_writer_preserves_late_matching_boundaries(
    tmp_path,
):
    """Null-only leading rows must not make boundary columns null-typed."""
    outside = [
        {
            "chromosome": "X",
            "position": index + 1,
            "canonical_id": f"X_{index + 1}_A_G",
            "input_variant_id": f"rsOutside{index + 1}",
            "reason": "missing_chromosome_reference",
            "observed_reference_ids": "",
            "action": "warning_skip",
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
        }
        for index in range(112)
    ]
    matching = {
        "chromosome": "1",
        "position": 2_373_168,
        "canonical_id": "1_2373168_A_G",
        "input_variant_id": "rsInside",
        "reason": "missing_from_reference",
        "observed_reference_ids": "",
        "action": "warning_skip",
        "reported_locus_boundary_status": "inside_reported_locus_boundary",
        "overlapping_genomic_loci": "1",
        "overlapping_locus_chromosome": "1",
        "overlapping_locus_start": 2_369_498,
        "overlapping_locus_end": 2_402_499,
        "locus_membership_interpretation": (
            "coordinate overlap only; excluded index was not used as an LD "
            "member or locus seed"
        ),
    }
    configuration = load_configuration().modules.ld_clumping

    output = ld_prune_standard._write_reference_exclusions(
        [*outside, matching],
        tmp_path,
        "STUDY",
        "EUR",
        configuration,
    )
    audit = pl.read_csv(
        output,
        separator="\t",
        infer_schema_length=None,
    )

    assert audit.columns == ld_prune_standard.LD_REFERENCE_EXCLUSION_COLUMNS
    assert audit.schema["overlapping_locus_start"] == pl.Int64
    assert audit.schema["overlapping_locus_end"] == pl.Int64
    assert audit.select(
        "overlapping_locus_chromosome",
        "overlapping_locus_start",
        "overlapping_locus_end",
    ).tail(1).row(0) == ("1", 2_369_498, 2_402_499)


def test_overlapping_reported_loci_cannot_produce_ambiguous_exclusion_audit():
    exclusion = {
        "chromosome": "1",
        "position": 175,
        "canonical_id": "1_175_A_G",
        "input_variant_id": "rsAmbiguous",
        "reason": "missing_from_reference",
        "observed_reference_ids": "",
        "action": "warning_skip",
    }
    overlapping_loci = pl.DataFrame({
        "Genomic_locus": [1, 2],
        "CHR": ["1", "1"],
        "START": [100, 150],
        "END": [200, 250],
    })

    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match="boundaries are invalid or overlap",
    ):
        ld_prune_standard._classify_reference_exclusions_by_locus(
            [exclusion],
            overlapping_loci,
        )


def test_missing_ld_endpoint_still_fails_after_vectorised_coverage_check(
    monkeypatch,
):
    payloads = iter([
        "1\t100\traw_index\t1_100_A_G\tA\tG\t0.2\n",
        "1\t100\traw_index\t1\t200\traw_missing\t0.8\n",
        "",
        "1\t100\traw_index\t1_100_A_G\tA\tG\t0.2\n",
    ])
    monkeypatch.setattr(
        ld_prune_standard,
        "_batch_tabix_query",
        lambda *_args, **_kwargs: next(payloads),
    )
    significant = pl.DataFrame({
        "chrcol": ["1"],
        "poscol": [100],
        "uniq_id": ["1_100_A_G"],
        "rsIDcol": ["rsIndex"],
    })

    with pytest.raises(
        ld_prune_standard.PipelineStageError,
        match=(
            "LD pair endpoints are absent.*missing_endpoint_ids=1.*"
            "examples=\\['raw_missing'\\]"
        ),
    ):
        ld_prune_standard.prepare_chromosome_reference(
            significant,
            "1",
            "EUR_chr1.ld.gz",
            "EUR_chr1.variants.tsv.gz",
            reverse_ld_path="EUR_chr1.reverse.ld.gz",
            orientation="upper_triangle_dual_index",
            tabix_bin="tabix",
            minimum_reference_maf=0.00001,
            missing_index_action="warning_skip",
            log=lambda _message: None,
        )


def test_scientific_settings_never_fall_back_to_packaged_defaults():
    """A resolved module configuration must not be mixed with global defaults."""
    with pytest.raises(
        Exception, match="scientific setting fall back to the packaged defaults",
    ):
        ld_prune_standard._resolve_standard_options(
            configuration=load_configuration().modules.ld_clumping,
            window_kb=None,
        )


def test_every_scientific_output_pattern_is_population_scoped():
    """Two populations in one directory must not overwrite each other."""
    layout = load_configuration().modules.ld_clumping.output_layout
    scientific = {
        name: pattern
        for name, pattern in layout.model_dump().items()
        if name not in layout.shared_cojo_directory_fields
    }
    assert scientific
    assert all("{population}" in pattern for pattern in scientific.values())
    assert layout.cojo_formatter_directory == "cojo/formatter"
    assert layout.cojo_root_directory == "cojo"
    assert layout.summary_csv.endswith(".csv")
    assert layout.html_report.endswith(".html")


def test_html_result_table_columns_are_complete_and_unique():
    columns = {
        "genomic_risk_loci": ["Genomic_locus"],
        "lead_snps": ["lead_SNP_id"],
        "independent_significant_snps": ["ind_sig_SNP_id"],
        "cojo_loci": ["COJO_locus"],
        "cojo_selected_signals": ["SNP"],
    }
    validated = LDClumpingReportingConfig(
        result_table_columns=columns,
    )
    assert validated.result_table_columns == columns

    with pytest.raises(ValueError, match="one or more unique column names"):
        LDClumpingReportingConfig(result_table_columns={
            **columns,
            "lead_snps": ["lead_SNP_id", "lead_SNP_id"],
        })

    with pytest.raises(ValueError, match="must define exactly these result tables"):
        LDClumpingReportingConfig(result_table_columns={
            "genomic_risk_loci": ["Genomic_locus"],
            "lead_snps": ["lead_SNP_id"],
        })

    with pytest.raises(ValueError, match="columns not produced"):
        LDClumpingReportingConfig(result_table_columns={
            **columns,
            "genomic_risk_loci": ["not_a_locus_column"],
        })
