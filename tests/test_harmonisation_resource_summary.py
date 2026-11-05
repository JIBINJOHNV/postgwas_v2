"""Typed resource summaries retain exact file audit and expand problems."""
from copy import deepcopy
from pathlib import Path
import pytest
from postgwas.config import load_configuration
from postgwas.core.input_validation import InputValidationSession, record_file_validation
from postgwas.core.validation_reporting import FileValidationDisplay
from postgwas.modules.harmonisation import resource_preflight as rp
from postgwas.modules.harmonisation.policies import default_policies
from test_harmonisation_resource_preflight import _resource_maps, _mock_bcftools, MAPPING, VCF_CONFIG


def validate(maps):
    return rp.validate_harmonisation_resource_maps(maps,grch_version='GRCh37',bcftools='bcftools',
        vcf_config=VCF_CONFIG,default_eaf_colmap=MAPPING,external_eaf_colmap=MAPPING,
        external_info_colmap=MAPPING,policies=default_policies(),require_default_eaf=True,
        resource_patterns={'dbsnp_file':'dbSNP_chr{chromosome}.vcf.gz',
            'default_comparison_af_file':'population_chr{chromosome}.vcf.gz'})


def test_resources_show_one_typed_table_and_unchanged_audit(tmp_path,monkeypatch,capsys):
    maps=_resource_maps(tmp_path)
    monkeypatch.setattr(rp,'run_checked_command',_mock_bcftools)
    with InputValidationSession() as session:
        display=FileValidationDisplay(session,load_configuration())
        result=validate(maps)
        before=deepcopy([r.to_dict() for r in session.records])
        display.flush()
        text=capsys.readouterr().out
        assert text.count('Resource summary')==1
        assert 'dbSNP VCF' in text and 'Population AF VCF' in text
        assert 'dbSNP_chr{chromosome}.vcf.gz' in text
        assert 'population_chr{chromosome}.vcf.gz' in text
        assert 'CSI' in text and 'TBI' in text and 'FAI' in text
        assert '2/2' in text
        assert 'column EUR' in ' '.join(text.replace('|',' ').split())
        assert 'CHECKS PASSED' not in text
        assert [r.to_dict() for r in session.records]==before
        assert result['indexed_vcfs_checked']==4
        display.flush()
        assert capsys.readouterr().out==''
        # A later recorded fact updates the group, not a per-chromosome card.
        record_file_validation(maps['1']['dbsnp_file'],'dbSNP',checks=('extra index check',),metrics={'variants':10})
        display.flush()
        update=capsys.readouterr().out
        assert 'Resource summary' in update and 'extra index check' in ' '.join(update.replace('|', ' ').split())
        assert 'File 1 ·' not in update
        assert 'Variants' in update and '10' in update


def test_problem_is_expanded_with_exact_path_and_partial_count(tmp_path,monkeypatch,capsys):
    maps=_resource_maps(tmp_path)
    monkeypatch.setattr(rp,'run_checked_command',_mock_bcftools)
    missing=Path(maps['2']['dbsnp_file']);missing.unlink()
    with InputValidationSession() as session:
        display=FileValidationDisplay(session,load_configuration())
        with pytest.raises(rp.ResourcePreflightError):validate(maps)
        display.flush()
        text=capsys.readouterr().out
        assert 'Resource summary' in text and '1/2' in text
        assert str(missing) in ''.join(line.strip() for line in text.splitlines()) and 'does not exist' in text
        assert any(r.status=='failed' and r.path==str(missing) for r in session.records)


def test_shared_roles_merge_build_inventories_and_expand_later_warning(tmp_path,capsys):
    with InputValidationSession() as session:
        display=FileValidationDisplay(session,load_configuration())
        paths=[tmp_path/'build37.fa',tmp_path/'build38.fa']
        for path in paths:
            record_file_validation(path,'FASTA',checks=('fasta index',))
            rp.register_resource_summaries({'1':{'genome_fasta_file':str(path)}},('genome_fasta_file',),
                patterns={'genome_fasta_file':'{build}_chr{chromosome}.fa'})
        rp.register_resource_summaries({'1':{'target_fasta':str(paths[0])}},('target_fasta',),
            patterns={'target_fasta':'{build}_chr{chromosome}.fa'})
        display.flush();text=capsys.readouterr().out
        assert text.count('Source-build FASTA')==1 and '2/2' in text
        assert 'Target-build FASTA' in text and '1/1' in text
        record_file_validation(paths[0],'FASTA',status='warning',message='index may be stale',checks=('index timestamp',))
        display.flush();text=capsys.readouterr().out
        assert 'index may be stale' in text and '0/1' in text and '1/2' in text
        assert str(paths[0]) in ''.join(line.strip() for line in text.splitlines())
