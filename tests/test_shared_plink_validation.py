"""Shared BIM/PLINK validation preserves method policies and pipeline reuse."""

from argparse import Namespace
from unittest.mock import patch

import pytest

from postgwas.config import load_configuration
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.plink import (
    count_plink_samples,
    validate_plink_bed_dimensions,
    validate_plink_bundle_dimensions,
    validate_plink_files,
)
from postgwas.core import variant_identifiers
from postgwas.modules.formatting.reference_identifiers import (
    BimIdentifierRequirement,
    configure_reference_variant_identifiers,
    configure_required_variant_identifier_type,
)


@pytest.fixture
def configuration():
    return load_configuration()


def _bim(tmp_path, rows, name="reference.bim"):
    path = tmp_path / name
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _inspect(path, configuration):
    formatting = configuration.modules.formatting
    return variant_identifiers.inspect_bim_identifier_type(
        path,
        column_roles=configuration.modules.magma.input.bim_columns,
        delimiter_pattern=configuration.modules.magma.input.table_delimiter_pattern,
        rsid_pattern=formatting.variant_identifiers.rsid_pattern,
        unique_id_template=formatting.variant_identifiers.unique_id_template,
        chromosome_prefix_pattern=formatting.chromosome_labels.prefix_pattern,
        chromosome_aliases=formatting.chromosome_labels.aliases,
    )


def _requirement(path, configuration, *, target="magma", consumer="MAGMA"):
    return BimIdentifierRequirement(
        consumer=consumer,
        formatter_target=target,
        bim_file=path,
        column_roles=configuration.modules.magma.input.bim_columns,
        delimiter_pattern=configuration.modules.magma.input.table_delimiter_pattern,
    )


@pytest.mark.parametrize("rows,expected", [
    (["1 rs1 0 100 G A", "1 rs2 0 200 T C"], "rsid"),
    (["1 1_100_A_G 0 100 G A", "1 1_200_C_T 0 200 C T"], "unique"),
])
def test_bim_convention_includes_both_allele_orders(tmp_path, configuration, rows, expected):
    with InputValidationSession() as session:
        result = _inspect(_bim(tmp_path, rows), configuration)
    assert result.identifier_type == expected
    assert result.variants == 2
    assert result.missing_ids == result.duplicate_ids == result.duplicate_rows == 0
    record = next(record for record in session.records if record.role == "PLINK BIM")
    assert record.metrics["identifier_type"] == expected
    assert "does not establish" in record.message


@pytest.mark.parametrize("identifier", ["1_101_A_G", "2_100_A_G", "1_100_A_C"])
def test_unique_id_must_agree_with_bim_coordinates_and_allele_pair(
    tmp_path, configuration, identifier,
):
    path = _bim(tmp_path, ["1 %s 0 100 G A" % identifier])
    with InputValidationSession() as session:
        with pytest.raises(ValueError, match="mixed or unsupported"):
            _inspect(path, configuration)
    assert any(record.metrics.get("other_ids") == 1 for record in session.records)


@pytest.mark.parametrize("position", ["1_00", "100.0", "1e2", "+100", "١٠٠", "-1", "NA"])
def test_bim_position_rejects_python_literals_and_nondecimal_text(tmp_path, configuration, position):
    path = _bim(tmp_path, ["1 rs1 0 %s G A" % position])
    with pytest.raises(ValueError, match="expected ASCII decimal digits"):
        _inspect(path, configuration)


@pytest.mark.parametrize("position", ["1", "000100", "2147483646"])
def test_bim_position_accepts_positive_ascii_decimal_text(tmp_path, configuration, position):
    path = _bim(tmp_path, ["1 rs1 0 %s G A" % position])
    assert _inspect(path, configuration).variants == 1


def test_bim_position_zero_remains_invalid(tmp_path, configuration):
    path = _bim(tmp_path, ["1 rs1 0 0 G A"])
    with pytest.raises(ValueError, match="non-positive position"):
        _inspect(path, configuration)


@pytest.mark.parametrize("token", [".", "-", "NA"])
def test_missing_identifiers_are_reported_not_accepted(tmp_path, configuration, token):
    path = _bim(tmp_path, ["1 rs1 0 100 G A", "1 %s 0 200 T C" % token])
    with InputValidationSession() as session:
        with pytest.raises(ValueError, match="missing=1"):
            _inspect(path, configuration)
    assert any(record.metrics.get("missing_ids") == 1 for record in session.records)


def test_duplicate_counts_do_not_replace_consumer_specific_rejection(tmp_path, configuration):
    path = _bim(tmp_path, ["1 rs1 0 100 G A", "1 rs1 0 100 T A", "1 rs1 0 100 C A"])
    with InputValidationSession() as session:
        result = _inspect(path, configuration)
    assert result.rsids == 3
    assert result.duplicate_ids == 1
    assert result.duplicate_rows == 2
    assert any(record.role == "PLINK BIM" and record.status == "warning" for record in session.records)


def test_same_bim_is_scanned_once_for_different_formatter_targets(tmp_path, configuration):
    path = _bim(tmp_path, ["1 rs1 0 100 G A"])
    args = Namespace()
    original = variant_identifiers._scan_bim_identifier_type
    with patch.object(variant_identifiers, "_scan_bim_identifier_type", wraps=original) as scanner:
        with InputValidationSession() as session:
            for target in ("magma", "gcta_gene"):
                with session.scope(target):
                    configure_reference_variant_identifiers(
                        args, configuration.modules.formatting,
                        [_requirement(path, configuration, target=target, consumer=target)],
                    )
        assert scanner.call_count == 1
    assert args.variant_id_types == {"magma": "rsid", "gcta_gene": "rsid"}
    assert session.records[0].consumers == ("magma", "gcta_gene")


def test_direct_calls_do_not_reuse_a_path_only_observation(tmp_path, configuration):
    path = _bim(tmp_path, ["1 rs1 0 100 G A"])
    args = Namespace()
    requirement = _requirement(path, configuration)
    configure_reference_variant_identifiers(args, configuration.modules.formatting, [requirement])
    path.write_text("1 1_100_A_G 0 100 G A\n", encoding="utf-8")
    with pytest.raises(ValueError, match="different identifier conventions"):
        configure_reference_variant_identifiers(args, configuration.modules.formatting, [requirement])


def test_changed_bim_cannot_reuse_session_evidence(tmp_path, configuration):
    path = _bim(tmp_path, ["1 rs1 0 100 G A"])
    with InputValidationSession():
        _inspect(path, configuration)
        path.write_text("1 rs1 0 100 G A\n1 rs2 0 200 G A\n", encoding="utf-8")
        with pytest.raises(RuntimeError, match="changed after pipeline preflight"):
            _inspect(path, configuration)


@pytest.mark.parametrize("same_call", [True, False])
def test_consumers_sharing_a_target_cannot_request_different_conventions(
    tmp_path, configuration, same_call,
):
    rsids = _bim(tmp_path, ["1 rs1 0 100 G A"], "rsids.bim")
    unique = _bim(tmp_path, ["1 1_100_A_G 0 100 G A"], "unique.bim")
    requirements = [_requirement(rsids, configuration), _requirement(unique, configuration, consumer="Other")]
    args = Namespace()
    with pytest.raises(ValueError, match="different identifier conventions"):
        if same_call:
            configure_reference_variant_identifiers(args, configuration.modules.formatting, requirements)
        else:
            for requirement in requirements:
                configure_reference_variant_identifiers(args, configuration.modules.formatting, [requirement])


def test_protocol_requirement_cannot_be_overwritten_by_bim_inference(tmp_path, configuration):
    args = Namespace()
    configure_required_variant_identifier_type(
        args, configuration.modules.formatting,
        consumer="Protocol", formatter_target="gcta_gene", required_type="unique",
    )
    path = _bim(tmp_path, ["1 rs1 0 100 G A"])
    with pytest.raises(ValueError, match="different identifier conventions"):
        configure_reference_variant_identifiers(
            args, configuration.modules.formatting,
            [_requirement(path, configuration, target="gcta_gene")],
        )


def test_companion_checks_do_not_require_genotypes_for_bim_only_consumer(tmp_path):
    path = _bim(tmp_path, ["1 rs1 0 100 G A"])
    assert validate_plink_files(path.with_suffix(""), [".bim"]) == {"bim": path}
    with InputValidationSession() as session:
        with pytest.raises(ValueError, match=r"reference\.bed.*reference\.fam"):
            validate_plink_files(path.with_suffix(""), [".bed", ".bim", ".fam"])
    assert sum(record.status == "failed" for record in session.records) == 2


@pytest.mark.parametrize("samples", [1, 4, 5, 8])
def test_bed_dimensions_use_ceil_samples_divided_by_four(tmp_path, samples):
    path = tmp_path / "reference.bed"
    expected = 3 + 2 * ((samples + 3) // 4)
    path.write_bytes(bytes((0x6C, 0x1B, 0x01)) + bytes(expected - 3))
    assert validate_plink_bed_dimensions(path, variants=2, samples=samples) == expected


@pytest.mark.parametrize("contents,message", [
    (b"bad", "SNP-major binary header"),
    (bytes((0x6C, 0x1B, 0x01)), "inconsistent with BIM/FAM"),
    (bytes((0x6C, 0x1B, 0x01, 0, 0)), "inconsistent with BIM/FAM"),
])
def test_bad_bed_bundle_is_rejected(tmp_path, contents, message):
    path = tmp_path / "reference.bed"
    path.write_bytes(contents)
    with InputValidationSession() as session:
        with pytest.raises(ValueError, match=message):
            validate_plink_bed_dimensions(path, variants=1, samples=1)
    assert any(record.status == "failed" for record in session.records)


def test_fam_sample_count_and_existing_blank_row_policies(tmp_path):
    path = tmp_path / "reference.fam"
    path.write_text("family sample 0 0 1 -9\n\n", encoding="utf-8")
    assert count_plink_samples(path, allow_blank_rows=True) == 1
    with pytest.raises(ValueError, match="expected 6"):
        count_plink_samples(path)


def test_complete_bundle_checks_and_cache(tmp_path):
    path = _bim(tmp_path, ["1 rs1 0 100 G A"])
    path.with_suffix(".bed").write_bytes(bytes((0x6C, 0x1B, 0x01, 0)))
    path.with_suffix(".fam").write_text("family sample 0 0 1 -9\n", encoding="utf-8")
    with InputValidationSession() as session:
        files = validate_plink_files(path.with_suffix(""), [".bed", ".bim", ".fam"])
        samples = count_plink_samples(files["fam"])
        validate_plink_bed_dimensions(files["bed"], variants=1, samples=samples)
    assert sum(record.role == "PLINK companion" for record in session.records) == 3
    assert any(record.role == "PLINK FAM" and record.metrics["samples"] == 1 for record in session.records)
    assert any(record.role == "PLINK BED" and record.metrics["bed_bytes"] == 4 for record in session.records)
    assert all(record.status == "passed" for record in session.records)


def test_bundle_reuses_fam_count_when_consumer_already_validated_it(tmp_path):
    path = _bim(tmp_path, ["1 rs1 0 100 G A"])
    path.with_suffix(".bed").write_bytes(bytes((0x6C, 0x1B, 0x01, 0)))
    path.with_suffix(".fam").write_text("family sample 0 0 1 -9\n", encoding="utf-8")
    files = validate_plink_files(path.with_suffix(""), [".bed", ".bim", ".fam"])
    samples = count_plink_samples(files["fam"])
    with patch("postgwas.core.plink.count_plink_samples", side_effect=AssertionError("Repeated FAM scan")):
        result = validate_plink_bundle_dimensions(files, variants=1, samples=samples)
    assert result == {"variants": 1, "samples": 1, "bed_bytes": 4}


@pytest.mark.parametrize("engine", ["susie", "finemap"])
def test_fine_mapping_pipeline_selects_identifier_type_before_formatter(
    tmp_path, monkeypatch, engine,
):
    from postgwas.modules.fine_mapping.service import preflight_fine_mapping_pipeline
    from preflight_support import pipeline_input_vcf_evidence

    path = _bim(tmp_path, ["1 1_100_A_G 0 100 G A"])
    path.with_suffix(".bed").write_bytes(bytes((0x6C, 0x1B, 0x01, 0)))
    path.with_suffix(".fam").write_text("family sample 0 0 1 -9\n", encoding="utf-8")
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._validate_runtime_tools", lambda _args: ((), None),
    )
    args = Namespace(finemap_method=engine, finemap_ld_reference=path.with_suffix(""))
    with InputValidationSession() as session:
        preflight_fine_mapping_pipeline(args, preflight_evidence=pipeline_input_vcf_evidence())
    assert args.variant_id_types[engine] == "unique"
    assert any(record.role == "PLINK BED" and record.metrics["variants"] == 1 for record in session.records)


def test_fine_mapping_pipeline_rejects_bad_bed_before_generated_inputs(tmp_path, monkeypatch):
    from postgwas.modules.fine_mapping.service import preflight_fine_mapping_pipeline
    from preflight_support import pipeline_input_vcf_evidence

    path = _bim(tmp_path, ["1 rs1 0 100 G A"])
    path.with_suffix(".bed").write_bytes(b"BED")
    path.with_suffix(".fam").write_text("family sample 0 0 1 -9\n", encoding="utf-8")
    monkeypatch.setattr(
        "postgwas.modules.fine_mapping.preflight._validate_runtime_tools", lambda _args: ((), None),
    )
    args = Namespace(finemap_method="susie", finemap_ld_reference=path.with_suffix(""))
    with pytest.raises(ValueError, match="SNP-major binary header"):
        preflight_fine_mapping_pipeline(args, preflight_evidence=pipeline_input_vcf_evidence())
