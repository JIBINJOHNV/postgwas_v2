"""Parent reports saved decisions, including evidence before a later failure."""
from copy import deepcopy
from io import StringIO
import pytest
from postgwas.core.pipeline_logging import PipelineLogger, print_chromosome_summary


def render(summary, status='OK'):
    return print_chromosome_summary(1,'study',1,status,1,'',summary,stream=StringIO(),width=120)


def test_completed_stage_evidence_survives_later_failure(tmp_path):
    logger = PipelineLogger('study',1,str(tmp_path))
    with logger.step(4,16,'Allele frequency','eaf',rows_in=100) as ctx:
        ctx.extra.update(final_variants=90, missing_eaf_count=3, eaf_degenerate_count=2,
            eaf_degenerate_action='reject', eaf_degenerate_removed=2,
            external_merge_usable_match_fraction=.9,external_merge_min_match_fraction=.8,
            post_orientation_duplicates={'duplicate_rows':5,'rows_removed':4},
            post_orientation_duplicate_policy={'conflicting_action':'fail_dataset','selection_order':['completeness','input_order']})
        ctx.rows_out=90
        logger.decide('Frequency direction',{'source':'external'},'Keep effect-allele frequency')
        logger.qc('SPECIAL_CHECK','Recorded policy kept rows',100,90,outcome='kept',matched=7)
    with pytest.raises(ValueError):
        with logger.step(5,16,'Sample size','size',rows_in=90):
            raise ValueError('sample counts rejected')
    summary=logger.summary();logger.close();before=deepcopy(summary)
    text=render(summary,'FAILED')
    assert 'Allele frequency' in text
    assert 'Keep effect-allele frequency' in text
    assert '90.00%' in text and '80.00%' in text
    assert 'Degenerate frequency' in text and 'reject' in text
    assert 'fail_dataset' in text and 'completeness' in text
    assert 'SPECIAL_CHECK' in text and 'matched 7' in text
    assert 'failed' in text and 'sample counts rejected' in text
    assert 'Imputation quality' not in text
    assert summary==before


@pytest.mark.parametrize('status',['ran','off','skipped','not recorded'])
def test_both_concordance_checks_have_explicit_status(status):
    summary={'stage_qc':{'effect_statistics_validation_qc':{
        'beta_se_z_concordance':{'status':status,'action':'reject','checked':1000,'discordant':7,'removed':7},
        'z_pval_concordance':{'status':'off','action':'warn'},'removed_total':7}}}
    text=render(summary)
    assert 'BETA/SE/Z status' in text and status in text
    assert 'Z/P status' in text and 'off' in text
    assert 'BETA, SE and Z' not in text
    assert 'BETA/SE/Z disagreements' in text and '7' in text


def test_duplicate_policy_is_saved_without_changing_scientific_rows(tmp_path):
    import polars as pl
    from postgwas.config import load_configuration
    from postgwas.modules.harmonisation.cli import _engine_defaults
    from postgwas.modules.harmonisation.allele_frequency import _resolve_reference_aligned_duplicates
    from postgwas.modules.harmonisation.policies import default_policies
    mapping={'chr_col':'CHR','pos_col':'POS','ea_col':'EA','oa_col':'OA','snp_id_col':'SNP','gwas_outputname':'study'}
    frame=pl.DataFrame({'CHR':['1'],'POS':[101],'EA':['G'],'OA':['A'],'SNP':['rs1']})
    logger=PipelineLogger('study',1,str(tmp_path/'logs'))
    policies=default_policies()
    with logger.step(4,16,'Allele frequency','eaf',rows_in=1) as ctx:
        retained, stats=_resolve_reference_aligned_duplicates(chromosome='1',df=frame,
            sample_column_dict={**mapping,'output_folder':str(tmp_path)},
            output_layout=_engine_defaults(load_configuration())['output_layout'],
            output_delimiter='\t',policies=policies,logger=logger,ctx=ctx,rejects=None)
    summary=logger.summary();logger.close()
    saved=summary['steps'][0]['extra']['post_orientation_duplicate_policy']
    assert saved['conflicting_action']==policies.get('duplicates.conflicting_action')
    assert saved['selection_order']==policies.get('duplicates.selection_order')
    assert retained.equals(frame) and stats['rows_removed']==0
