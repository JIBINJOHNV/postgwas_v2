"""Scientific/preflight failure stops one sample, never subsequent valid samples."""
from argparse import Namespace
import csv
from pathlib import Path

import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation import cli
from postgwas.modules.harmonisation.service import PipelineError
from postgwas.modules.harmonisation.study_properties import StudyPropertyError
from postgwas.modules.harmonisation.sample_sheet import load_harmonisation_sample_sheet
from test_harmonisation_batch_preflight import successful_engine_result, sheet as sample_records, write_sheet


def setup(tmp_path, monkeypatch):
    resources = tmp_path / 'resources'; resources.mkdir()
    config = load_configuration(cli_overrides={'resources.root': str(resources),
        'run.output_directory': str(tmp_path / 'out'), 'execution.retries': 2,
        'modules.harmonisation.concordance_validation.enabled': False})
    records=sample_records(tmp_path,2)
    for record,name in zip(records,('bad','good')):record['dataset_id']=name
    fixture=write_sheet(tmp_path/'setup.csv',records)
    rows=load_harmonisation_sample_sheet(fixture)
    monkeypatch.setattr(cli, 'preflight_configured_reference_scope', lambda *a, **k: None)
    return config, rows


@pytest.mark.parametrize('error_type', [StudyPropertyError, PipelineError])
def test_preparation_failure_is_not_retried_and_next_sample_runs(tmp_path, monkeypatch, error_type):
    config, rows = setup(tmp_path, monkeypatch)
    calls = []
    def engine(**kwargs):
        name = kwargs['sample_column_dict']['gwas_outputname']
        calls.append((name, bool(kwargs.get('preparation_only'))))
        if name == 'bad':
            raise error_type('whole-study decision rejected')
        return successful_engine_result(**kwargs)
    monkeypatch.setattr(cli, 'run_harmonisation_pipeline', engine)
    with pytest.raises(PipelineError, match='batch is incomplete'):
        cli._run_validated_rows(Namespace(), config, rows)
    assert calls == [('bad', True), ('good', True), ('good', False)]
    with (tmp_path/'out/run_metadata/harmonisation_run_summary.csv').open() as handle:
        results = {r['dataset_id']: r for r in csv.DictReader(handle)}
    assert results['bad']['status'] == 'FAILED'
    assert results['bad']['failure_stage'] == 'dataset preparation'
    assert results['good']['status'] == 'OK'


def test_invalid_sample_sheet_row_does_not_block_valid_sample(tmp_path, monkeypatch):
    config, rows = setup(tmp_path, monkeypatch)
    fixture = tmp_path/'setup.csv'
    with fixture.open() as handle:
        reader=csv.DictReader(handle); template=next(reader); header=reader.fieldnames
    records=[]
    for name in ('bad','good'):
        record={**template,'dataset_id':name,'input_file':str(tmp_path/(name+'.tsv'))}
        if name == 'bad':record['effect_allele_frequency_column']=''
        records.append(record)
    sheet=tmp_path/'samples.csv'
    with sheet.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=header);writer.writeheader();writer.writerows(records)
    calls=[]
    def engine(**kwargs):
        calls.append(kwargs['sample_column_dict']['gwas_outputname'])
        return successful_engine_result(**kwargs)
    monkeypatch.setattr(cli,'run_harmonisation_pipeline',engine)
    with pytest.raises(PipelineError):
        cli.run_harmonisation(Namespace(sample_sheet=sheet),configuration=config,resolved_executables={})
    assert calls == ['good','good']
    with (tmp_path/'out/run_metadata/harmonisation_run_summary.csv').open() as handle:
        results={r['dataset_id']:r for r in csv.DictReader(handle)}
    assert results['bad']['status']=='PREFLIGHT_FAILED'
    assert results['good']['status']=='OK'


def test_report_collision_cannot_truncate_declared_input(tmp_path, monkeypatch):
    config, rows = setup(tmp_path, monkeypatch)
    target = cli._screen_report_paths(config, rows)['bad']
    target.parent.mkdir(parents=True)
    target.write_text('scientific input must survive')
    rows[0] = rows[0].model_copy(update={'input_file': str(target)})
    with pytest.raises(cli.ConfigurationError, match='overwrite a declared scientific input'):
        cli._run_validated_rows(Namespace(), config, rows)
    assert target.read_text() == 'scientific input must survive'


def report_target(config, kind):
    root=Path(config.run.output_directory)
    metadata=cli._runtime_path(config,root,'top_metadata_directory')
    if kind in ('run_summary_file','preflight_report_file','resolved_config_file','run_html_report'):
        return cli._runtime_path(config,metadata,kind)
    if kind=='top_log':return cli._top_run_log(config,root)
    if kind=='screen':return cli._screen_report_paths(config,[],extra_dataset_ids=['bad'])['bad']
    if kind=='html':return cli._html_report_paths(config,[],extra_dataset_ids=['bad'])[0]['bad']
    metadata=cli._runtime_path(config,root,'metadata_directory',dataset_id='bad')
    if kind=='dataset_log':return metadata/config.logging.filename
    return cli._runtime_path(config,metadata,kind.removeprefix('dataset_'))


@pytest.mark.parametrize('kind',['run_summary_file','preflight_report_file','resolved_config_file',
    'run_html_report','top_log','screen','html','dataset_log','dataset_command_file',
    'dataset_sample_sheet_row_file','dataset_resolved_config_file','dataset_supplied_config_file'])
@pytest.mark.parametrize('invalid',[True,False])
def test_public_reporting_never_overwrites_declared_input(tmp_path,monkeypatch,kind,invalid):
    from test_harmonisation_batch_preflight import sheet,write_sheet
    config,_=setup(tmp_path,monkeypatch)
    records=sheet(tmp_path,1);records[0]['dataset_id']='bad'
    if invalid:records[0]['effect_type']='typo'
    target=report_target(config,kind);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text('original scientific source\n')
    records[0]['input_file']=str(target)
    manifest=write_sheet(tmp_path/'samples.csv',records)
    monkeypatch.setattr(cli,'run_harmonisation_pipeline',lambda **kw:pytest.fail('unsafe output must stop before engine'))
    with pytest.raises(cli.ConfigurationError,match='overwrite a declared scientific input'):
        cli.run_harmonisation(Namespace(sample_sheet=manifest),configuration=config,resolved_executables={})
    assert target.read_text()=='original scientific source\n'


@pytest.mark.parametrize('alias',['symlink','hardlink'])
@pytest.mark.parametrize('contents',['external source sentinel\n',''])
def test_reporting_preserves_external_template_and_filesystem_aliases(tmp_path,monkeypatch,alias,contents):
    from test_harmonisation_batch_preflight import sheet,write_sheet
    config,_=setup(tmp_path,monkeypatch)
    target=report_target(config,'run_summary_file');target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(contents)
    external=tmp_path/'GRCh37_external_chr1.tsv'
    if alias=='symlink':external.symlink_to(target)
    else:external.hardlink_to(target)
    records=sheet(tmp_path,1);records[0].update(dataset_id='bad',effect_type='typo',
        external_eaf_file=str(tmp_path/'{build}_external_chr{chromosome}.tsv'),external_eaf_column='EUR')
    manifest=write_sheet(tmp_path/'samples.csv',records)
    with pytest.raises(cli.ConfigurationError,match='overwrite a declared scientific input'):
        cli.run_harmonisation(Namespace(sample_sheet=manifest),configuration=config,resolved_executables={})
    assert target.read_text()==contents


def test_failure_finalization_cannot_overwrite_an_input(tmp_path,monkeypatch):
    config,rows=setup(tmp_path,monkeypatch)
    target=report_target(config,'run_summary_file');target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text('failure finalization sentinel\n')
    rows=[rows[0].model_copy(update={'input_file':str(target)})]
    with pytest.raises(cli.ConfigurationError,match='overwrite a declared scientific input'):
        cli._write_preflight_failure_run_summary(config,rows,ValueError('shared setup failed'))
    assert target.read_text()=='failure finalization sentinel\n'


@pytest.mark.parametrize('source_kind',['replacement','replacement_tilde','reference','reference_index'])
def test_reporting_cannot_overwrite_configured_resource_aliases(tmp_path,monkeypatch,source_kind):
    from test_harmonisation_batch_preflight import sheet,write_sheet
    config,_=setup(tmp_path,monkeypatch)
    target=report_target(config,'run_summary_file');target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text('resource sentinel\n')
    module=config.modules.harmonisation
    if source_kind.startswith('replacement'):
        resource=tmp_path/'GRCh37_replacement_chr1.tsv'
        if source_kind=='replacement_tilde':
            original=Path.expanduser
            monkeypatch.setattr(Path,'expanduser',lambda path: Path(str(tmp_path)+'/'+str(path)[2:])
                if str(path).startswith('~/') else original(path))
        config.modules.harmonisation=module.model_copy(update={
            'eaf_replacement_file':('~/{build}_replacement_chr{chromosome}.tsv' if source_kind=='replacement_tilde'
                else str(tmp_path/'{build}_replacement_chr{chromosome}.tsv')),
            'eaf_replacement_column':'EUR'})
    else:
        pattern=module.resource_layout.default_eaf if source_kind=='reference' else module.resource_layout.dbsnp
        resource=Path(config.resources.root)/pattern.format(build='GRCh37',chromosome='1',
            source=module.default_eaf.source if source_kind=='reference' else module.reference.dbsnp_source)
        if source_kind=='reference_index':resource=Path(str(resource)+'.csi')
    resource.parent.mkdir(parents=True,exist_ok=True);resource.hardlink_to(target)
    records=sheet(tmp_path,1);records[0].update(dataset_id='bad',effect_type='typo')
    manifest=write_sheet(tmp_path/'samples.csv',records)
    with pytest.raises(cli.ConfigurationError,match='overwrite a declared scientific input'):
        cli.run_harmonisation(Namespace(sample_sheet=manifest),configuration=config,resolved_executables={})
    assert target.read_text()=='resource sentinel\n'
