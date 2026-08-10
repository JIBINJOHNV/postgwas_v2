"""Regression contracts for complete and validated chromosome VCF merges."""

from pathlib import Path
from unittest.mock import patch

import pytest

from postgwas.config import load_configuration
from postgwas.core.paths import configured_output_path
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.vcf_processing import (
    LiftoverFailureRateError,
    VcfMergeError,
    annotate_and_liftover_vcf,
    concat_vcfs_by_build,
)


def _configuration():
    config = load_configuration().modules.harmonisation
    return dict(config.output_layout.root), config.vcf_processing.model_dump()


def _arguments(tmp_path, policies=None):
    output_layout, vcf_config = _configuration()
    return {
        "output_dir": str(tmp_path),
        "gwas_outputname": "study",
        "grch_version": "GRCh37",
        "expected_chromosomes": ["1"],
        "threads": 2,
        "output_layout": output_layout,
        "executables": {"bcftools": "bcftools", "tabix": "tabix"},
        "vcf_config": vcf_config,
        "policies": policies or default_policies(),
    }


def _required_inputs(tmp_path):
    output_layout, _vcf_config = _configuration()
    paths = [
        configured_output_path(
            tmp_path, output_layout[pattern], dataset_id="study", chromosome="1",
            build="GRCh37", target_build="GRCh38",
        )
        for pattern in (
            "chromosome_annotated_vcf",
            "chromosome_lifted_vcf",
            "chromosome_raw_vcf",
        )
    ]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"v" * 200)
        Path(str(path) + ".tbi").write_bytes(b"index")
    return paths


def _annotation_arguments(tmp_path):
    output_layout, vcf_config = _configuration()
    input_vcf = configured_output_path(
        tmp_path,
        output_layout["chromosome_raw_vcf"],
        dataset_id="study",
        chromosome="1",
        build="GRCh37",
        target_build="GRCh38",
    )
    input_vcf.parent.mkdir(parents=True, exist_ok=True)
    input_vcf.write_bytes(b"input-vcf")
    return {
        "output_dir": str(tmp_path),
        "gwas_outputname": "study",
        "chromosome": "1",
        "external_eaf_file": "external-af.vcf.gz",
        "default_dbsnp_file": "dbsnp.vcf.gz",
        "genome_fasta_file": "GRCh37.fa",
        "target_genome_fasta_file": "GRCh38.fa",
        "gff_file": "genes.gff3.gz",
        "chain_file": "GRCh37_to_GRCh38.chain",
        "grch_version": "GRCh37",
        "threads": 1,
        "output_layout": output_layout,
        "executables": {
            "bash": "bash", "bcftools": "bcftools", "tabix": "tabix",
        },
        "vcf_config": vcf_config,
        "policies": default_policies(),
    }


def test_merge_failure_default_is_fail():
    assert default_policies().get("vcf.on_merge_failure") == "fail"
    assert default_policies().get("vcf.concat_max_attempts") == 2
    assert default_policies().get("vcf.liftover_fail_fraction") == 0.5
    _output_layout, vcf_config = _configuration()
    assert vcf_config["required_merge_groups"] == [
        "input_build", "target_build", "raw_gwas2vcf",
    ]


def test_no_chromosome_inputs_fail_all_required_output_groups(tmp_path):
    with patch(
        "postgwas.modules.harmonisation.vcf_processing.run_checked_command"
    ) as command:
        with pytest.raises(VcfMergeError) as error:
            concat_vcfs_by_build(**_arguments(tmp_path))

    command.assert_not_called()
    message = str(error.value)
    assert "GRCh37" in message
    assert "GRCh38" in message
    assert "gwas2vcf_GRCh37" in message
    assert "chromosome 1 input is missing" in message


def test_zero_survivor_liftover_fails_independently_of_fraction(tmp_path):
    arguments = _annotation_arguments(tmp_path)
    reject_vcf = configured_output_path(
        tmp_path,
        arguments["output_layout"]["chromosome_not_lifted_vcf"],
        dataset_id="study",
        chromosome="1",
        build="GRCh37",
        target_build="GRCh38",
    )
    reject_vcf.write_bytes(b"rejected-vcf")
    counts = [100, 100, 100, 100, 100, 0, 100]

    with patch(
        "postgwas.modules.harmonisation.vcf_processing._run_bcftools_step"
    ), patch(
        "postgwas.modules.harmonisation.vcf_processing.get_vcf_variant_count",
        side_effect=counts,
    ), pytest.raises(LiftoverFailureRateError, match="zero variants") as error:
        annotate_and_liftover_vcf(**arguments)

    assert "100 annotated variants entered liftover" in str(error.value)
    assert "reported rejected count: 100" in str(error.value)
    assert "regardless of vcf.liftover_fail_fraction" in str(error.value)


def test_zero_record_required_chromosome_vcf_fails_validation(tmp_path):
    inputs = _required_inputs(tmp_path)

    def command(arguments, _purpose, **_kwargs):
        if len(arguments) > 2 and arguments[1:3] == ["index", "-n"]:
            return "0\n"
        raise AssertionError("unexpected command: %r" % arguments)

    with patch(
        "postgwas.modules.harmonisation.vcf_processing.run_checked_command",
        side_effect=command,
    ), pytest.raises(VcfMergeError, match="zero records in a required VCF"):
        concat_vcfs_by_build(**_arguments(tmp_path))

    assert all(path.is_file() for path in inputs)


def test_zero_record_required_merged_vcf_fails_validation(tmp_path):
    inputs = _required_inputs(tmp_path)

    def command(arguments, _purpose, **_kwargs):
        if len(arguments) > 1 and arguments[1] == "concat":
            output = Path(arguments[arguments.index("--output") + 1])
            output.write_bytes(b"v" * 200)
            Path(str(output) + ".tbi").write_bytes(b"index")
            return ""
        if len(arguments) > 2 and arguments[1:3] == ["index", "-n"]:
            return "0\n" if "merged" in str(arguments[-1]) else "1\n"
        raise AssertionError("unexpected command: %r" % arguments)

    with patch(
        "postgwas.modules.harmonisation.vcf_processing.run_checked_command",
        side_effect=command,
    ), pytest.raises(VcfMergeError, match="required merged VCF contains zero records"):
        concat_vcfs_by_build(**_arguments(tmp_path))

    assert all(path.is_file() for path in inputs)
    assert not list(tmp_path.glob("study_*_merged.vcf.gz"))


def test_continue_records_missing_required_groups_without_hiding_them(tmp_path):
    policies = default_policies().with_overrides({"vcf.on_merge_failure": "continue"})

    result = concat_vcfs_by_build(**_arguments(tmp_path, policies=policies))

    assert set(result["merge_failures"]) == {
        "GRCh37", "GRCh38", "gwas2vcf_GRCh37",
    }
    assert set(result["merge_failure_details"]) == set(result["merge_failures"])
    assert set(result["required_merge_failures"]) == set(result["merge_failures"])
    assert result["optional_merge_failures"] == []
    assert result["merge_status"] == "PARTIAL"
    assert result["grch37"] is None
    assert result["grch38"] is None
    assert result["gwas2vcf"] is None
    assert "notlifted_GRCh38" not in result["merge_failures"]


def test_continue_marks_successfully_merged_subset_as_partial(tmp_path):
    _required_inputs(tmp_path)
    policies = default_policies().with_overrides({"vcf.on_merge_failure": "continue"})
    arguments = _arguments(tmp_path, policies=policies)
    arguments["expected_chromosomes"] = ["1", "2"]

    def command(command_arguments, _purpose, **_kwargs):
        if len(command_arguments) > 1 and command_arguments[1] == "concat":
            output = Path(command_arguments[command_arguments.index("--output") + 1])
            output.write_bytes(b"v" * 200)
            Path(str(output) + ".tbi").write_bytes(b"index")
            return ""
        if len(command_arguments) > 2 and command_arguments[1:3] == ["index", "-n"]:
            return "1\n"
        raise AssertionError("unexpected command: %r" % command_arguments)

    with patch(
        "postgwas.modules.harmonisation.vcf_processing.run_checked_command",
        side_effect=command,
    ):
        result = concat_vcfs_by_build(**arguments)

    assert result["merge_status"] == "PARTIAL"
    assert set(result["required_merge_failures"]) == {
        "GRCh37", "GRCh38", "gwas2vcf_GRCh37",
    }
    assert all(Path(result[key]).is_file() for key in ("grch37", "grch38", "gwas2vcf"))
    assert all(
        "chromosome 2 input is missing" in result["merge_failure_details"][label][0]
        for label in result["required_merge_failures"]
    )


def test_transient_concat_failure_is_retried_and_can_recover(tmp_path):
    _required_inputs(tmp_path)
    attempts = {}

    def command(arguments, _purpose, **_kwargs):
        if len(arguments) > 1 and arguments[1] == "concat":
            output = Path(arguments[arguments.index("--output") + 1])
            attempts[output.name] = attempts.get(output.name, 0) + 1
            if output.name == "study_GRCh37_merged.vcf.gz" and attempts[output.name] == 1:
                raise RuntimeError("temporary concat failure")
            output.write_bytes(b"v" * 200)
            Path(str(output) + ".tbi").write_bytes(b"index")
            return ""
        if len(arguments) > 2 and arguments[1:3] == ["index", "-n"]:
            return "1\n"
        raise AssertionError("unexpected command: %r" % arguments)

    with patch(
        "postgwas.modules.harmonisation.vcf_processing.run_checked_command",
        side_effect=command,
    ):
        result = concat_vcfs_by_build(**_arguments(tmp_path))

    assert attempts["study_GRCh37_merged.vcf.gz"] == 2
    assert result["merge_attempts"]["GRCh37"] == 2
    assert result["merge_failures"] == []
    assert result["merge_status"] == "OK"


def test_optional_merge_failure_does_not_downgrade_required_outputs(tmp_path):
    _required_inputs(tmp_path)
    output_layout, _vcf_config = _configuration()
    optional_input = configured_output_path(
        tmp_path,
        output_layout["chromosome_not_lifted_vcf"],
        dataset_id="study",
        chromosome="1",
        build="GRCh37",
        target_build="GRCh38",
    )
    optional_input.write_bytes(b"v" * 200)
    Path(str(optional_input) + ".tbi").write_bytes(b"index")

    def command(arguments, _purpose, **_kwargs):
        if len(arguments) > 1 and arguments[1] == "concat":
            output = Path(arguments[arguments.index("--output") + 1])
            if output.name == "study_notlifted_GRCh38_merged.vcf.gz":
                raise RuntimeError("optional concat failure")
            output.write_bytes(b"v" * 200)
            Path(str(output) + ".tbi").write_bytes(b"index")
            return ""
        if len(arguments) > 2 and arguments[1:3] == ["index", "-n"]:
            return "1\n"
        raise AssertionError("unexpected command: %r" % arguments)

    with patch(
        "postgwas.modules.harmonisation.vcf_processing.run_checked_command",
        side_effect=command,
    ):
        result = concat_vcfs_by_build(**_arguments(tmp_path))

    assert result["required_merge_failures"] == []
    assert result["optional_merge_failures"] == ["notlifted_GRCh38"]
    assert result["merge_status"] == "OK"
    assert result["merge_attempts"]["notlifted_GRCh38"] == 2


def test_missing_merged_index_is_a_failure_and_inputs_are_retained(tmp_path):
    inputs = _required_inputs(tmp_path)

    def command(arguments, _purpose, **kwargs):
        if len(arguments) > 2 and arguments[1:3] == ["index", "-n"]:
            return "1\n"
        if len(arguments) > 1 and arguments[1] == "concat":
            output = Path(arguments[arguments.index("--output") + 1])
            output.write_bytes(b"v" * 200)
            missing = [
                str(path) for path in kwargs.get("expected_outputs", [])
                if not Path(path).is_file() or Path(path).stat().st_size == 0
            ]
            if missing:
                raise RuntimeError("expected output is missing or empty: %s" % missing)
        raise AssertionError("unexpected command: %r" % arguments)

    with patch(
        "postgwas.modules.harmonisation.vcf_processing.run_checked_command",
        side_effect=command,
    ):
        with pytest.raises(VcfMergeError, match="expected output is missing or empty"):
            concat_vcfs_by_build(**_arguments(tmp_path))

    assert all(path.is_file() for path in inputs)
    assert not list(tmp_path.glob("study_*_merged.vcf.gz"))


def test_success_requires_readable_merged_vcf_and_index(tmp_path):
    _required_inputs(tmp_path)

    def command(arguments, _purpose, **_kwargs):
        if len(arguments) > 1 and arguments[1] == "concat":
            output = Path(arguments[arguments.index("--output") + 1])
            output.write_bytes(b"v" * 200)
            Path(str(output) + ".tbi").write_bytes(b"index")
            return ""
        if len(arguments) > 2 and arguments[1:3] == ["index", "-n"]:
            return "1\n"
        raise AssertionError("unexpected command: %r" % arguments)

    with patch(
        "postgwas.modules.harmonisation.vcf_processing.run_checked_command",
        side_effect=command,
    ):
        result = concat_vcfs_by_build(**_arguments(tmp_path))

    assert result["merge_failures"] == []
    assert result["merge_failure_details"] == {}
    for key in ("grch37", "grch38", "gwas2vcf"):
        assert Path(result[key]).is_file()
        assert Path(result[key] + ".tbi").is_file()
