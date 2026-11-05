"""Detailed study evidence survives the real visible/hidden transcript path."""

import json
import os
from pathlib import Path
import subprocess
import sys
import textwrap

import pytest

from postgwas.core.pipeline_logging import _PREFIX_WIDTH


@pytest.mark.parametrize("show_screen", [False, True])
def test_detailed_report_reaches_canonical_and_dataset_logs_in_both_display_modes(
    tmp_path, show_screen,
):
    script = textwrap.dedent("""
        import json
        from pathlib import Path
        import sys
        import polars as pl
        from postgwas.core.pipeline_logging import PipelineLogger
        from postgwas.core.screen_logging import ScreenSettings, record_screen
        from postgwas.modules.harmonisation.allele_frequency import validate_study_allele_frequency
        from postgwas.modules.harmonisation.cli import _DatasetScreenRouter
        from postgwas.modules.harmonisation.policies import default_policies
        from postgwas.modules.harmonisation.service import _announce, _study_decisions_block
        from postgwas.modules.harmonisation.study_properties import resolve_study_properties

        output = Path(sys.argv[1])
        show_screen = sys.argv[2] == 'True'
        policies = default_policies().with_overrides({
            'effect.type': 'beta', 'pvalue.type': 'raw',
        })
        frame = pl.DataFrame({
            'CHR': ['1', '1', '1'], 'POS': [100, 200, 300],
            'EA': ['A', 'A', 'T'], 'OA': ['G', 'C', 'C'],
            'BETA': [-.2, .1, .3], 'P': [.1, .5, 1.], 'AF': [.2, .6, .8],
        })
        original = frame.clone()
        columns = {
            'chr_col': 'CHR', 'pos_col': 'POS', 'ea_col': 'EA', 'oa_col': 'OA',
            'beta_or_col': 'BETA', 'pval_col': 'P', 'eaf_col': 'AF',
        }
        reference_mapping = {
            'chr': 'CHROM', 'pos': 'POS', 'a1': 'ALT', 'a2': 'REF', 'delimiter': 'tab',
        }
        # If reporting introduced a comparison/read for this accepted branch,
        # the real preflight would fail on this nonexistent source.
        resources = {'1': {
            'default_eaf_file': str(output / 'not_needed_reference.tsv'),
            'default_eaf_reference_column': 'EUR',
        }}
        logger = PipelineLogger(
            'study', 'dataset', str(output), policies=policies, screen_level='ERROR',
        )
        try:
            with record_screen(ScreenSettings(show_screen, output / 'screen.log')):
                with _DatasetScreenRouter(
                    {'study': output / 'study_screen_report.txt'}, display=show_screen,
                ) as router:
                    router.select('study')
                    decisions = resolve_study_properties(
                        frame, columns, policies=policies, logger=logger,
                    )
                    decisions, prepared = validate_study_allele_frequency(
                        frame, columns, decisions, resources, reference_mapping,
                        policies=policies, logger=logger,
                    )
                    _announce(logger, _study_decisions_block(
                        decisions, policies=policies, resource_maps=resources,
                    ))
        finally:
            logger.close()
        (output / 'observed.json').write_text(json.dumps({
            'original_unchanged': frame.equals(original),
            'prepared_unchanged': prepared.equals(original),
            'effect_type': decisions['effect_type'],
            'effect_source': decisions['effect_type_source'],
            'pvalue_type': decisions['pvalue_type'],
            'frequency_source': decisions['eaf_selected_source'],
            'frequency_validation': decisions['eaf_reference_validation'],
            'frequency_values': prepared['AF'].to_list(),
        }))
    """)
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", POLARS_MAX_THREADS="4")
    completed = subprocess.run(
        [sys.executable, "-B", "-c", script, str(tmp_path), str(show_screen)],
        capture_output=True, text=True, timeout=60, check=False, env=environment,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    title = "Study-wide decisions and evidence"
    if show_screen:
        assert completed.stdout.count(title) == 1
    else:
        assert completed.stdout == ""

    observed = json.loads((tmp_path / "observed.json").read_text())
    assert observed == {
        "original_unchanged": True,
        "prepared_unchanged": True,
        "effect_type": "beta",
        "effect_source": "policy",
        "pvalue_type": "raw",
        "frequency_source": "study_supplied",
        "frequency_validation": {
            "status": "not_required",
            "decision_reason": "study_frequency_not_predominantly_minor",
        },
        "frequency_values": [.2, .6, .8],
    }
    paths = (
        tmp_path / "screen.log", tmp_path / "study_screen_report.txt",
        tmp_path / "study_dataset.log",
    )
    for path in paths:
        raw = path.read_text(encoding="utf-8")
        assert "\x1b" not in raw
        assert raw.count(title) == 1
        if path.name == "study_dataset.log":
            raw = "\n".join(line[_PREFIX_WIDTH:] for line in raw.splitlines())
        report = " ".join(raw.split()).split(title, 1)[1]
        for section in (
            "Effect type", "Standard-error scale", "P-value type",
            "BETA/Z direction for missing-SE recovery", "Allele frequency",
            "Strand consensus",
        ):
            assert section in report
        assert "YAML policy declaration" in report
        assert "Not required; input was not MAF-like" in report
        assert "Declared study EAF retained" in report
        assert "Detailed evidence" in report
        assert "Direct EAF correlation" not in report
        assert "Comparable variants" not in report
        assert "Folded-MAF correlation" not in report
        assert "Comparison reference" not in report


@pytest.mark.parametrize("show_screen,dataset_count,late_failure,width", [
    (False, 1, False, 76), (True, 1, False, 140),
    (False, 1, True, 140), (True, 1, True, 76),
    (False, 2, False, 76), (True, 2, False, 140),
    (False, 2, True, 140), (True, 2, True, 76),
])
def test_real_preparation_reports_once_per_dataset_across_transcript_modes(
        tmp_path, show_screen, dataset_count, late_failure, width):
    script = textwrap.dedent("""
        import json
        from pathlib import Path
        import sys
        import polars as pl
        from _pytest.monkeypatch import MonkeyPatch
        from postgwas.core.screen_logging import ScreenSettings, record_screen
        from postgwas.modules.harmonisation.cli import _DatasetScreenRouter
        from test_harmonisation_study_report_integration import _prepare

        output = Path(sys.argv[1])
        show_screen = sys.argv[2] == 'True'
        dataset_count = int(sys.argv[3])
        late_failure = sys.argv[4] == 'True'
        names = ['study_one', 'study_two'][:dataset_count]
        observed = {}
        with record_screen(ScreenSettings(show_screen, output / 'screen.log')):
            with _DatasetScreenRouter({name: output / (name + '_screen_report.txt')
                                      for name in names}, display=show_screen) as router:
                for index, name in enumerate(names):
                    router.select(name)
                    fail = late_failure and index == dataset_count - 1
                    with MonkeyPatch.context() as patch:
                        try:
                            prepared = _prepare(output / name, patch, dataset_id=name,
                                automatic_inference=True,
                                failure_stage='partition' if fail else None)
                        except OSError as error:
                            decisions = getattr(error, 'study_decisions', {})
                            build = getattr(error, 'genome_build_info', {}) or {}
                            observed[name] = {'status': 'failed', 'message': str(error),
                                'effect_type': decisions.get('effect_type'),
                                'pvalue_type': decisions.get('pvalue_type'),
                                'info_type': decisions.get('info_score_type'),
                                'build': build.get('inferred_build'),
                                'build_chromosomes': sorted(build.get('per_chromosome', {}))}
                        else:
                            frame = pl.read_parquet(prepared.chromosome_files['1'])
                            observed[name] = {'status': 'prepared',
                                'effect_type': prepared.study_decisions['effect_type'],
                                'pvalue_type': prepared.study_decisions['pvalue_type'],
                                'info_type': prepared.study_decisions['info_score_type'],
                                'build': prepared.genome_build_info['inferred_build'],
                                'build_chromosomes': sorted(prepared.genome_build_info['per_chromosome']),
                                'frequency_values': frame['AF'].to_list(),
                                'beta_values': frame['BETA'].to_list()}
        (output / 'observed.json').write_text(json.dumps(observed))
    """)
    tests_root = str(Path(__file__).resolve().parent)
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", POLARS_MAX_THREADS="4",
        COLUMNS=str(width), PYTHONPATH=tests_root + os.pathsep + os.environ.get("PYTHONPATH", ""))
    completed = subprocess.run([sys.executable, "-B", "-c", script, str(tmp_path),
        str(show_screen), str(dataset_count), str(late_failure)], capture_output=True,
        text=True, timeout=60, check=False, env=environment)
    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    title = "Study-wide decisions and evidence"
    if show_screen:
        assert completed.stdout.count(title) == dataset_count
    else:
        assert completed.stdout == ""
    canonical = (tmp_path / "screen.log").read_text()
    assert canonical.count(title) == dataset_count
    observed = json.loads((tmp_path / "observed.json").read_text())
    for index, name in enumerate(["study_one", "study_two"][:dataset_count]):
        fail = late_failure and index == dataset_count - 1
        outcome = observed[name]
        assert outcome["status"] == ("failed" if fail else "prepared")
        assert outcome["effect_type"] == "beta" and outcome["pvalue_type"] == "raw"
        assert outcome["info_type"] == "standard_info" and outcome["build"] == "GRCh37"
        assert outcome["build_chromosomes"] == ["1"]
        if fail:
            assert outcome["message"] == "partition fixture failure"
        else:
            assert outcome["frequency_values"] == [.1, .2, .7, .8]
            assert outcome["beta_values"] == [-.1, .1, .2, -.2]
        for path in (tmp_path / (name + "_screen_report.txt"),
                     tmp_path / name / "logs" / (name + "_dataset.log")):
            raw = path.read_text()
            assert "\x1b" not in raw
            assert raw.count(title) == 1
            assert raw.count("Selected build") == 1
            assert sum(line.rstrip().endswith("Input validation") for line in raw.splitlines()) == 1
            assert "Chromosome 1 build evidence" not in raw
            assert "Effect type automatically inferred" not in raw
            assert "P-value type automatically inferred" not in raw
            if path.name.endswith("_dataset.log"):
                raw = "\n".join(line[_PREFIX_WIDTH:] for line in raw.splitlines())
            report = " ".join(raw.split()).split(title, 1)[1]
            assert "effect_type" in report and "pvalue.type" in report
            assert "warning" in report.lower() and "INFO" in report
            assert ("Dataset stopped" in report) is fail
            if fail:
                assert "partition fixture failure" in report
                assert "accepted for assignment" not in report
