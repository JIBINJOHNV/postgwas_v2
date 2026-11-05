"""Structural startup and immutable per-study preparation contracts."""
import csv
import json
from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock

import pytest

from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.paths import configured_output_path
from postgwas.core.input_validation import InputValidationSession, capture_preflight_file_identities
from postgwas.modules.harmonisation import cli, service, resource_preflight as resources
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.sample_sheet import load_harmonisation_sample_sheet
from postgwas.modules.harmonisation.shared.runtime import prepared_work_fingerprint
from test_harmonisation_resource_preflight import _resource_maps, _mock_bcftools, MAPPING, VCF_CONFIG

# Existing local compressed-input tests import this path; the portable batch
# fixtures below generate their own inputs and never read this local manifest.
FIXTURE = Path(__file__).parent / 'data/harmonisation/manifest_v2.csv'


def successful_engine_result(**kwargs):
    """Persist the required report contract for a mocked successful engine."""
    columns=kwargs['sample_column_dict']
    manifest=configured_output_path(columns['output_folder'],
        kwargs['default_cfg']['output_layout']['run_manifest'],dataset_id=columns['gwas_outputname'])
    manifest.parent.mkdir(parents=True,exist_ok=True)
    assessment=manifest.parent/'fixture_assessment.json'
    assessment.write_text(json.dumps({'raw':{},'qc_passed':{}}))
    manifest.write_text(json.dumps({'status':'OK','qc_assessment':{'reports':{'json':str(assessment)}}}))
    return {'status':'OK','manifest':str(manifest)}


def sheet(tmp_path, count=3):
    """Create independent valid studies without repository or private data."""
    baseline = {
        'config_version': 2, 'trait_type': 'case_control',
        'chromosome_column': 'CHR', 'position_column': 'POS',
        'effect_allele_column': 'EA', 'other_allele_column': 'OA',
        'effect_allele_frequency_column': 'EAF', 'effect_column': 'BETA',
        'effect_type': 'beta', 'standard_error_column': 'SE',
        'p_value_column': 'P', 'p_value_type': 'raw',
        'control_count': 2000, 'control_count_column': None,
        'case_count': 1000, 'imputation_info_column': 'INFO',
        'z_score_column': None, 'external_eaf_file': None,
        'external_eaf_column': None, 'delimiter': 'tab',
    }
    records = []
    for index in range(count):
        source = tmp_path / ('study%d.tsv' % index)
        source.write_text('CHR\tPOS\tEA\tOA\tEAF\tBETA\tSE\tP\tINFO\n'
            '1\t101\tC\tA\t0.2\t0.1\t0.1\t0.3173105079\t0.95\n'
            '1\t202\tG\tA\t0.7\t-0.2\t0.1\t0.0455002639\t0.9\n', encoding='utf-8')
        records.append({**baseline, 'dataset_id': 'study%d' % index, 'input_file': str(source)})
    return records


def write_sheet(path, records):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=records[0])
        writer.writeheader()
        writer.writerows(records)
    return path


def config(tmp_path):
    root = tmp_path / 'references'
    root.mkdir()
    return load_configuration(cli_overrides={
        'resources.root': str(root), 'run.output_directory': str(tmp_path / 'results'),
        'execution.retries': 2,
    })


def report(cfg):
    path = cfg.run.output_directory / cfg.modules.harmonisation.runtime.top_metadata_directory / cfg.modules.harmonisation.runtime.preflight_report_file
    with path.open() as handle:
        return list(csv.DictReader(handle))


def test_all_schema_rows_and_independent_fields_are_reported(tmp_path):
    records = sheet(tmp_path)
    records[0].update(control_count='bad', control_count_column=None, effect_column=None, z_score_column=None)
    records[2].update(effect_allele_frequency_column=None, external_eaf_file=None, external_eaf_column=None)
    records[2]['input_file'] = str(tmp_path / 'missing.tsv')
    path = write_sheet(tmp_path / 'studies.csv', records)
    issues = []
    rows = load_harmonisation_sample_sheet(path, issues=issues)
    assert [row.dataset_id for row in rows] == ['study1']
    text = '\n'.join(item['problem'] for item in issues)
    for expected in ('row 2', 'row 4', 'positive whole number', 'effect_column', 'no allele-frequency source', 'missing.tsv'):
        assert expected in text
    with pytest.raises(ConfigurationError, match='row 4'):
        load_harmonisation_sample_sheet(path)


def test_malformed_row_and_blank_lines_keep_actual_line_numbers(tmp_path):
    records = sheet(tmp_path)
    records[2]['config_version'] = 99
    path = write_sheet(tmp_path / 'studies.csv', records)
    lines = path.read_text().splitlines()
    path.write_text('\n'.join([lines[0], '', lines[1], 'broken,row', lines[2], lines[3]])+'\n')
    issues = []
    load_harmonisation_sample_sheet(path, issues=issues)
    assert {item['row'] for item in issues} == {4, 6}


def test_bad_csv_dispatches_only_valid_rows_with_failed_row_evidence(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, 'preflight_configured_reference_scope', lambda *a, **k: None)
    cfg = config(tmp_path)
    records = sheet(tmp_path, 4)
    records[0]['config_version'] = 99
    path = write_sheet(tmp_path / 'studies.csv', records)
    events = []
    def engine(**kwargs):
        name = kwargs['sample_column_dict']['gwas_outputname']
        events.append(name)
        raise AssertionError('Malformed CSV must stop before scientific preparation')
    monkeypatch.setattr(cli, 'run_harmonisation_pipeline', engine)
    launch = Mock(return_value={'dispatched': True})
    monkeypatch.setattr(cli, '_run_validated_rows', launch)
    assert cli.run_harmonisation(Namespace(sample_sheet=path), configuration=cfg, resolved_executables={}) == {'dispatched': True}
    assert events == []
    launch.assert_called_once()
    assert [r.dataset_id for r in launch.call_args.args[2]] == ['study1', 'study2', 'study3']
    assert set(launch.call_args.kwargs['preflight_failures']) == {'study0'}
    assert (cfg.run.output_directory / 'run_metadata' / cfg.modules.harmonisation.runtime.resolved_config_file).is_file()
    rows = report(cfg)
    assert {row['dataset_id'] for row in rows if row['status']=='failed'} == {'study0'}
    assert next(row for row in rows if row['dataset_id']=='study2')['status']=='passed'
    log = (cfg.run.output_directory / 'run_metadata' / cfg.logging.filename).read_text()
    assert 'study0' in log


def test_startup_dispatches_rows_without_preparing_all_studies(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, 'preflight_configured_reference_scope', lambda *a, **k: None)
    cfg = config(tmp_path)
    path = write_sheet(tmp_path / 'studies.csv', sheet(tmp_path))
    events = []
    def engine(**kwargs):
        name = kwargs['sample_column_dict']['gwas_outputname']
        events.append(name)
        raise AssertionError('Startup must not read GWAS content')
    def execute(*args, **kwargs):
        assert events == []
        assert 'prepared_runs' not in kwargs
        assert [row.dataset_id for row in args[2]] == ['study0', 'study1', 'study2']
        return {'ok': True}
    monkeypatch.setattr(cli, 'run_harmonisation_pipeline', engine)
    monkeypatch.setattr(cli, '_run_validated_rows', execute)
    assert cli.run_harmonisation(Namespace(sample_sheet=path), configuration=cfg, resolved_executables={}) == {'ok': True}
    assert all(row['status']=='passed' for row in report(cfg))


def test_preparation_rejects_changed_inputs_and_config(tmp_path):
    path = tmp_path / 'study.tsv'
    path.write_text('initial\n')
    identities = capture_preflight_file_identities([path])
    prepared = service.PreparedHarmonisation(None, service.configuration_digest([{'study': 'one'}, {}]), identities)
    prepared.validate({'study': 'one'}, {})
    with pytest.raises(service.PipelineError, match='configuration changed'):
        prepared.validate({'study': 'two'}, {})
    path.write_text('changed\n')
    with pytest.raises(service.PipelineError, match='changed after'):
        prepared.validate({'study': 'one'}, {})


def validate_maps(maps):
    return resources.validate_harmonisation_resource_maps(
        maps, bcftools='bcftools', vcf_config=VCF_CONFIG,
        grch_version="GRCh37",
        default_eaf_colmap=MAPPING, external_eaf_colmap=MAPPING,
        external_info_colmap=MAPPING, policies=default_policies(), require_default_eaf=True,
    )


@pytest.mark.parametrize('alias', ['same', 'symlink', 'hardlink'])
def test_same_file_maf_conflict_reports_both_chromosomes(tmp_path, monkeypatch, alias):
    maps = _resource_maps(tmp_path)
    for chrom, values in maps.items():
        original = Path(values['default_eaf_file'])
        target = tmp_path / ('alias'+chrom+'.tsv')
        if alias == 'symlink':
            target.symlink_to(original)
        elif alias == 'hardlink':
            target.hardlink_to(original)
        else:
            target = original
        values['user_eaf_file'] = str(target)
        values['user_eaf_column'] = 'EUR'
    monkeypatch.setattr(resources, 'run_checked_command', _mock_bcftools)
    with pytest.raises(resources.ResourcePreflightError) as error:
        validate_maps(maps)
    issues = [item for item in error.value.issues if 'same file' in item['problem']]
    assert {chrom for item in issues for chrom in item['chromosomes']} == {'1', '2'}
    assert all('default_eaf.source' in item['problem'] for item in issues)


@pytest.mark.parametrize("alias", [False, True])
def test_same_file_non_maf_and_exact_cutoff_remain_valid(tmp_path, monkeypatch, alias):
    maps = _resource_maps(tmp_path)
    for chrom, values in maps.items():
        path = Path(values['default_eaf_file'])
        # Exactly 95%, not greater than 95%: keep declared EAF semantics.
        path.write_text('CHROM\tPOS\tREF\tALT\tEUR\n' + ''.join(
            '%s\t%d\tA\tG\t%s\n' % (chrom, index+1, .2 if index<19 else .8)
            for index in range(20)
        ))
        target = path
        if alias:
            target = tmp_path / ('same_reference_'+chrom+'.tsv')
            target.symlink_to(path)
        values.update(user_eaf_file=str(target), user_eaf_column='EUR')
    monkeypatch.setattr(resources, 'run_checked_command', _mock_bcftools)
    assert validate_maps(maps)['status'] == 'passed'


def test_shared_reference_checks_reused_but_changed_index_rechecked(tmp_path, monkeypatch):
    maps = _resource_maps(tmp_path)
    calls = Mock(side_effect=_mock_bcftools)
    monkeypatch.setattr(resources, 'run_checked_command', calls)
    with InputValidationSession(pipeline=False):
        validate_maps(maps)
        initial = calls.call_count
        validate_maps(maps)
        assert calls.call_count == initial
        index = Path(maps['1']['default_comparison_af_file']+'.tbi')
        index.write_text('changed index\n')
        validate_maps(maps)
        assert calls.call_count == initial+1


def test_config_header_and_resource_problems_are_combined(tmp_path):
    rows = load_harmonisation_sample_sheet(write_sheet(tmp_path / 'studies.csv', sheet(tmp_path, 1)))
    cfg = config(tmp_path)
    defaults = cli._engine_defaults(cfg, resolved_executables={})
    engine = cli.to_harmonisation_input(rows[0], resource_directory=cfg.resources.root,
        output_directory=cfg.run.output_directory, output_layout=defaults['output_layout'])
    engine['ea_col'] = 'MISSING_ALLELE'
    engine['eaffile'] = str(tmp_path / 'missing_eaf.tsv')
    with pytest.raises(service.ConfigError) as error:
        service.validate_harmonisation_inputs(engine, defaults)
    for value in ('MISSING_ALLELE', 'missing_eaf.tsv', 'GRCh37_check_file', 'GRCh38_check_file'):
        assert value in str(error.value)


def test_generated_work_uses_content_fingerprint_not_metadata(tmp_path):
    import os
    path = tmp_path / 'work.parquet'
    path.write_bytes(b'original prepared bytes')
    contract = service.configuration_digest([{}, {}])
    with InputValidationSession(pipeline=False):
        prepared = service.PreparedHarmonisation(None, contract, (), (prepared_work_fingerprint(path, error_type=service.PipelineError),))
        original = path.stat()
        os.utime(path, ns=(original.st_atime_ns, original.st_mtime_ns + 1000000000))
        prepared.validate({}, {})
        path.write_bytes(b'changed prepared bytes!')
        os.utime(path, ns=(original.st_atime_ns, original.st_mtime_ns))
        with pytest.raises(service.PipelineError, match='work file changed'):
            prepared.validate({}, {})


def test_invalid_path_does_not_hide_later_row_errors(tmp_path):
    records = sheet(tmp_path)
    invalid = tmp_path / 'symlink_loop'
    invalid.symlink_to(invalid)
    records[0]['input_file'] = str(invalid)
    records[2]['config_version'] = 99
    path = write_sheet(tmp_path / 'paths.csv', records)
    issues = []
    rows = load_harmonisation_sample_sheet(path, issues=issues)
    assert [row.dataset_id for row in rows] == ['study1']
    assert {item['row'] for item in issues} == {2, 4}
