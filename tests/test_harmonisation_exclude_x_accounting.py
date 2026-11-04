"""Explicit chromosome exclusions remain separate in every accounting view."""

from copy import deepcopy
from types import SimpleNamespace

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.modules.harmonisation import service
from postgwas.modules.harmonisation.cli import _engine_defaults
from postgwas.modules.harmonisation.html_report import _input_section, _overview_inputs
from postgwas.modules.harmonisation.qc_reporting import (
    harmonisation_qc_summary_lines,
    harmonisation_qc_takeaway_lines,
)
from postgwas.modules.harmonisation.report_steps import dataset_steps
from postgwas.modules.harmonisation.policies import load_policies
from postgwas.modules.harmonisation.rejects import SOURCE_INPUT_ROW_COLUMN
from postgwas.modules.qc_summary.assessment import assess_variant_table


EXCLUSION_METRIC = "total_variant_removed_excluded_chromosomes"


def _counts():
    return {
        "total_variant_infile": 12,
        "total_variant_read": 12,
        EXCLUSION_METRIC: 2,
        "total_variant_removed_unsupported_chromosomes": 1,
        "total_variant_removed_missing_values": 1,
        "total_variant_removed_null_coords": 1,
        "total_variant_removed_non_standard_alleles": 1,
        "total_variant_removed_duplicates": 1,
        "total_variant_remaining_for_harmonisation": 5,
        "total_variant_ready_snps": 5,
        "total_variant_ready_indels_or_other": 0,
    }


def test_input_terminal_separates_user_exclusions_from_invalid_chromosomes():
    text = service._input_validation_summary_block(
        input_variants=12, variants_read=12, excluded_chromosomes=2,
        missing_required=1, invalid_coordinates=1, unsupported_chromosomes=1,
        non_standard_alleles=1, duplicate_variants=1, ready_variants=5,
        ready_snps=5, ready_indels_or_other=0,
    )
    assert "User-excluded chromosomes" in text
    assert "2" in next(line for line in text.splitlines() if "User-excluded" in line)
    assert "1" in next(line for line in text.splitlines() if "Unsupported" in line)


def test_html_cards_and_fallback_qc_preserve_recorded_exclusion_count():
    counts = _counts()
    original = deepcopy(counts)
    manifest = {"pre_vcf_summary": counts}
    for text in (
        _input_section({}, manifest),
        " ".join(body for _title, body, _span in _overview_inputs({}, manifest)),
    ):
        assert "User-excluded chromosomes" in text
        assert "Unsupported chromosomes" in text
    card = dataset_steps(manifest)[2]
    assert card["details"]["Observed evidence"]["Input accounting"][EXCLUSION_METRIC] == 2
    assert card["metrics"]["User-excluded chromosomes"] == 2
    table = service._fallback_qc_frame(counts)
    assert table.loc[table.metric == EXCLUSION_METRIC, "GLOBAL_VALUE"].tolist() == [2]
    assert counts == original


def test_final_terminal_views_show_explicit_exclusion_once_per_view(tmp_path):
    table = tmp_path / "vcf_records.tsv"
    table.write_text(
        "CHROM\tPOS\tREF\tALT\tINFO_AF\tINFO_EXTERNAL_AF\t"
        "FORMAT_AF\tFORMAT_SI\tFORMAT_LP\tFORMAT_NEF\n"
        "1\t100\tA\tG\t0.2\t0.21\t0.2\t0.9\t8\t100\n"
        "1\t101\tC\tT\t0.3\t0.31\t0.3\t0.9\t8\t100\n"
    )
    configuration = load_configuration().modules.qc_summary
    assessment = assess_variant_table(
        table, configuration, "GRCh37", configuration.reference_af_column,
        configuration.vcf_fields.model_dump(),
        delimiter=configuration.table.delimiter,
        null_values=configuration.table.null_values,
    )
    assessment["raw_vcf"] = str(tmp_path / "study.vcf.gz")
    counts = _counts()
    detailed = "\n".join(harmonisation_qc_summary_lines(counts, assessment))
    takeaways = "\n".join(harmonisation_qc_takeaway_lines(
        counts, assessment, genome_build_info={}, study_decisions={},
        population_frequency={}, dataset_summary={}, merge_summary={},
        primary_outputs={}, final_status="OK",
    ))
    assert detailed.count("User-excluded chromosomes") == 1
    assert takeaways.count("User-excluded chromosomes") == 1
    assert "2" in next(line for line in takeaways.splitlines() if "User-excluded" in line)


@pytest.mark.parametrize("enabled,expected", [
    (False, ["2", "X", "1", "chrX", "23", "PAR1"]),
    (True, ["2", "1"]),
])
def test_startup_inventory_applies_canonical_exclusion_without_reordering(
    monkeypatch, enabled, expected,
):
    defaults = _engine_defaults(load_configuration())
    defaults["policies"]["chromosome"]["exclude_x"] = enabled
    seen = []

    def preflight(chromosomes, **kwargs):
        seen.append((list(chromosomes), kwargs["grch_version"]))
        return {}, {}

    monkeypatch.setattr(service, "preflight_harmonisation_resources", preflight)
    scope = SimpleNamespace(builds=["GRCh37", "GRCh38"],
                            chromosomes=["2", "X", "1", "chrX", "23", "PAR1"])
    service.preflight_configured_reference_scope(defaults, "/unused", scope)
    assert seen == [(expected, "GRCh37"), (expected, "GRCh38")]


def test_exact_resource_planning_never_resolves_excluded_chromosome_paths(monkeypatch):
    defaults = _engine_defaults(load_configuration())
    policies = load_policies({"chromosome": {"exclude_x": True}})
    resolved = []

    def resource_map(*, chromosome, **_kwargs):
        resolved.append(chromosome)
        return {"fixture": chromosome}

    monkeypatch.setattr(service, "build_resource_map", resource_map)
    monkeypatch.setattr(service, "validate_harmonisation_resource_maps",
                        lambda resources, **_kwargs: {"chromosomes": list(resources)})
    maps, summary = service.preflight_harmonisation_resources(
        ["2", "X", "1", "chrX", "PAR2"], grch_version="GRCh37", target_build="GRCh38",
        resource_folder="/unused", user_eaf_file=None, user_info_file=None,
        user_eaf_column=None, user_info_column=None,
        default_eaf_reference_source=defaults["default_eaf_reference_source"],
        default_comparison_af_file=defaults["default_comparison_af_file"],
        resource_layout=defaults["resource_layout"], dbsnp=defaults["default_dbsnp"],
        default_eaf_reference_column=defaults["default_eaf_reference_column"],
        default_comparison_af_column=defaults["default_comparison_af_column"],
        default_eaf_colmap=defaults["default_eaf_colmap"],
        external_eaf_colmap=defaults["external_eaf_colmap"],
        external_info_colmap=defaults["external_info_colmap"], executables={"bcftools": "unused"},
        vcf_config=defaults["vcf_processing"], policies=policies, require_default_eaf=True,
    )
    assert resolved == ["2", "1"]
    assert list(maps) == ["2", "1"]
    assert summary == {"chromosomes": ["2", "1"]}


def _preparation_arguments(tmp_path, *, only_x=False):
    module = load_configuration().modules.harmonisation
    source = tmp_path / "study.tsv"
    retained = [
        ["one", "1", "10", "G", "A", "-.1", ".1", ".2", ".1", ".9"],
        ["two", "1", "20", "G", "A", ".1", ".1", ".3", ".2", ".9"],
        ["three", "1", "30", "G", "A", ".2", ".1", ".4", ".7", ".9"],
        ["four", "1", "40", "G", "A", "-.2", ".1", ".5", ".8", ".9"],
    ]
    excluded = ["excluded", "chrX", ".", "bad", ".", "bad", "bad", "bad", "bad", "bad"]
    rows = [excluded] if only_x else [retained[0], excluded, *retained[1:]]
    source.write_text("\t".join([
        "SNP", "CHR", "POS", "EA", "OA", "BETA", "SE", "P", "AF", "INFO",
    ]) + "\n" + "\n".join("\t".join(row) for row in rows) + "\n")
    mapping = {
        "gwas_outputname": "study", "sumstat_file": str(source),
        "output_folder": str(tmp_path / "out"), "chr_col": "CHR", "pos_col": "POS",
        "snp_id_col": "SNP", "ea_col": "EA", "oa_col": "OA",
        "beta_or_col": "BETA", "se_col": "SE", "pval_col": "P", "eaf_col": "AF",
        "declared_effect_type": "beta", "declared_pvalue_type": "raw", "ncontrol": "1000",
        "imp_info_col": "INFO", "info_source": "internal",
    }
    return dict(
        sumstat_file=str(source), sample_column_dict=mapping,
        output_dir_path=tmp_path / "out", resource_folder=str(tmp_path),
        default_eaf_reference_source="fixture", default_comparison_af_file="fixture",
        dbsnp="fixture", resource_layout={}, output_layout=dict(module.output_layout.root),
        gwas2vcf_input=module.gwas2vcf_input.model_dump(), executables={},
        vcf_config=module.vcf_processing.model_dump(), user_eaf_file=None,
        user_eaf_column=None, user_info_file=None, user_info_column=None,
        default_eaf_reference_column="EUR", default_comparison_af_column="EUR",
        build_reference_files={}, build_check_colmap={}, default_eaf_colmap={},
        external_eaf_colmap={}, external_info_colmap={},
        external_reference_staging=module.external_reference_staging.model_dump(),
        policies=load_policies({"chromosome": {"exclude_x": True}}),
        input_line_count=len(rows),
    )


def test_real_preparation_excludes_x_before_every_scientific_decision(
    tmp_path, monkeypatch, capsys,
):
    arguments = _preparation_arguments(tmp_path)
    source = arguments["sumstat_file"]
    original_source = open(source).read()
    observations = {}

    def check_frame(name, frame):
        observations[name] = frame["SNP"].to_list()
        assert frame["CHR"].unique().to_list() == ["1"]
        assert frame[SOURCE_INPUT_ROW_COLUMN].to_list() == [1, 3, 4, 5]

    def build(frame, *_args, **_kwargs):
        check_frame("build", frame)
        return {"inferred_build": "GRCh37", "input_variants": 4,
                "testable_variants": 4, "untestable_variants": 0}

    def strand(*_args, **_kwargs):
        observations["strand"] = True
        return {"strand": "forward", "strand_mode": "auto", "forward": 4,
                "reverse": 0, "ambiguous": 0, "informative": 4,
                "dominant_fraction": 1., "minimum_informative_variants": 1,
                "consensus_threshold": .99, "reason": "known-forward fixture"}

    def resources(chromosomes, **_kwargs):
        observations["resources"] = list(chromosomes)
        assert list(chromosomes) == ["1"]
        maps = {"1": {"default_eaf_file": str(tmp_path / "unused.tsv"),
                "default_eaf_reference_column": "EUR", "user_eaf_file": None,
                "user_info_file": None}}
        summary = dict.fromkeys(("unique_resource_files_checked", "indexed_vcfs_checked",
            "frequency_vcf_headers_checked", "fasta_indexes_checked", "table_headers_checked",
            "structured_text_files_checked"), 0)
        summary.update(chromosomes=["1"], file_identities=[], resources_by_chromosome=maps)
        return maps, summary

    monkeypatch.setattr(service, "infer_genome_build", build)
    monkeypatch.setattr(service, "resolve_strand_consensus", strand)
    monkeypatch.setattr(service, "preflight_harmonisation_resources", resources)
    for name in ("resolve_study_properties", "validate_study_allele_frequency", "write_chromosome_partitions"):
        actual = getattr(service, name)

        def observe(frame, *args, _actual=actual, _name=name, **kwargs):
            check_frame(_name, frame)
            return _actual(frame, *args, **kwargs)

        monkeypatch.setattr(service, name, observe)

    with PipelineLogger("study", "dataset", str(tmp_path / "logs")) as logger:
        prepared = service._prepare_dataset_for_chromosome_processing(**arguments, logger=logger)

    assert set(observations) == {"build", "strand", "resources", "resolve_study_properties",
                                 "validate_study_allele_frequency", "write_chromosome_partitions"}
    assert prepared.rows_read == prepared.input_variants == 5
    assert prepared.removed_excluded_chromosomes == 1
    assert prepared.removed_unsupported_chromosomes == 0
    assert prepared.ready_variants == 4
    assert prepared.partition_rows == {"1": 4}
    assert prepared.field_completeness["input_rows"] == 4
    assert prepared.field_completeness["chromosome_exclusion"] == {
        "enabled": True, "excluded_chromosomes": ["X"], "excluded_rows": 1,
        "input_rows": 5, "retained_rows": 4,
    }
    partition = pl.read_parquet(prepared.chromosome_files["1"])
    assert partition["SNP"].to_list() == ["one", "two", "three", "four"]
    assert partition["BETA"].to_list() == [-.1, .1, .2, -.2]
    assert open(source).read() == original_source
    assert "User-excluded chromosomes" in capsys.readouterr().out


def test_x_only_preparation_fails_before_scientific_decisions_and_flushes_audit(
    tmp_path, monkeypatch,
):
    arguments = _preparation_arguments(tmp_path, only_x=True)
    for name in ("infer_genome_build", "resolve_strand_consensus", "resolve_study_properties",
                 "preflight_harmonisation_resources", "validate_study_allele_frequency",
                 "write_chromosome_partitions"):
        monkeypatch.setattr(service, name, lambda *args, _name=name, **kwargs:
                            pytest.fail(_name + " must not run for excluded-only data"))
    with PipelineLogger("study", "dataset", str(tmp_path / "logs")) as logger:
        with pytest.raises(ValueError, match="No variants remain after --exclude-chr-x"):
            service._prepare_dataset_for_chromosome_processing(**arguments, logger=logger)
    log = next((tmp_path / "logs").glob("*dataset.log")).read_text()
    assert "FAILED" in log
    assert "No variants remain after --exclude-chr-x" in log
    assert not list((tmp_path / "out").rglob("*.parquet"))
    reject_files = list((tmp_path / "out").rglob("*input*.tsv*"))
    assert len(reject_files) == 1
    rejects = pl.read_csv(reject_files[0], separator="\t")
    assert rejects.height == 1
    assert rejects["reject_reason"].to_list() == ["chromosome_excluded_by_user"]
