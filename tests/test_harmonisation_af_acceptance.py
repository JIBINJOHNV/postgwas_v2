"""Whole-study AF acceptance: configured gates, boundaries and preserved evidence."""
import polars as pl
import pytest
import yaml

from postgwas.config import load_run_configuration_for_module
from postgwas.config.exporter import render_module_configuration
from postgwas.modules.harmonisation import allele_frequency as af
from postgwas.modules.harmonisation.allele_frequency import StudyAlleleFrequencyError
from postgwas.modules.harmonisation.policies import default_policies, load_policies, PolicyError
from postgwas.modules.harmonisation.rejects import SOURCE_INPUT_ROW_COLUMN
from postgwas.modules.harmonisation.shared.runtime import NullStepContext
from postgwas.modules.harmonisation.study_properties import _decide_eaf_is_maf
from postgwas.modules.harmonisation.study_reporting import render_study_decisions

COLUMNS = {'chr_col': 'CHR', 'pos_col': 'POS', 'ea_col': 'EA', 'oa_col': 'OA', 'eaf_col': 'AF'}
MAPPING = {'chr': 'CHROM', 'pos': 'POS', 'a1': 'ALT', 'a2': 'REF', 'delimiter': 'tab'}


def prepare(frame, policies=None, reference=None, **options):
    policies = policies or default_policies()
    decisions = {}
    _decide_eaf_is_maf(frame, COLUMNS, policies, NullStepContext(), decisions, {})
    return af.validate_study_allele_frequency(frame, COLUMNS, decisions, {
        chrom: {'default_eaf_file': str(reference), 'default_eaf_reference_column': 'EUR'}
        for chrom in frame['CHR'].unique()}, MAPPING, policies=policies, **options)


def fixture(tmp_path):
    reference = tmp_path / 'reference.tsv'
    reference.write_text('CHROM\tPOS\tREF\tALT\tEUR\n' + ''.join(
        f'1\t{i}\tG\tA\t{value}\n' for i, value in enumerate([.8, .7, .6, .65, .75], 1)))
    frame = pl.DataFrame({'CHR': ['1'] * 5, 'POS': list(range(1, 6)), 'EA': ['A'] * 5,
        'OA': ['G'] * 5, 'AF': [.2, .3, .4, .35, .25]}).with_row_index(SOURCE_INPUT_ROW_COLUMN)
    return reference, frame, default_policies().with_overrides({'eaf.maf_reference_min_overlap': 2})


@pytest.fixture
def decisions():
    return {'eaf_column': 'AF', 'eaf_is_maf_initial': True, 'eaf_selected_source': 'study_supplied',
        'eaf_reference_validation': {'comparable_variants': 1100, 'minimum_overlap': 1000,
            'correlation_method': 'pearson', 'direct_eaf_correlation': .98, 'minimum_correlation': .90,
            'mean_eaf_error': .02, 'maximum_mean_absolute_difference': .10,
            'reference_effect_allele_minor_fraction': .99}}


def render(decisions):
    return render_study_decisions(decisions)


def normal(text):
    return ' '.join(text.split())


def _paired_comparison(study, reference=None, *, total=None, overrides=None, whole_study=True):
    from postgwas.modules.harmonisation.allele_frequency import _compare_frequency_pairs
    frame = pl.DataFrame({'study': study, 'reference': reference if reference is not None else study})
    policies = default_policies().with_overrides(overrides or {})
    return _compare_frequency_pairs(frame, 'study', 'reference', policies, NullStepContext(),
        accept_direct_eaf=whole_study, study_variant_count=total if total is not None else frame.height)


def _variable_minor_values(count):
    return [.1 + .3 * i / max(count - 1, 1) for i in range(count)]


def test_direct_correlation_uses_shared_threshold():
    import math
    reference = _variable_minor_values(1000)
    supplied = [value + .08 * math.sin(i) for i, value in enumerate(reference)]
    decision, evidence = _paired_comparison(supplied, reference)
    assert .80 < evidence['direct_eaf_correlation'] < .90
    assert evidence['minimum_maf_correlation'] == .90
    assert evidence['direct_eaf_acceptance_reason'] == 'correlation_below_minimum'
    assert decision == 'inconclusive'


def test_direct_correlation_threshold_is_inclusive():
    import math
    reference = _variable_minor_values(1000)
    supplied = [value + .08 * math.sin(i) for i, value in enumerate(reference)]
    measured = pl.DataFrame({'a': supplied, 'b': reference}).select(pl.corr('a', 'b')).item()
    decision, evidence = _paired_comparison(supplied, reference,
        overrides={'eaf.maf_reference_min_correlation': measured})
    assert evidence['direct_eaf_correlation'] == measured
    assert decision == 'eaf'
    assert evidence['direct_eaf_accepted'] is True


@pytest.mark.parametrize('minor_count,expected', [(900, 'inconclusive'), (901, 'eaf')])
def test_reference_minor_fraction_requires_strictly_more_than_threshold(minor_count, expected):
    values = _variable_minor_values(minor_count) + [.6] * (1000 - minor_count)
    decision, evidence = _paired_comparison(values)
    assert evidence['reference_effect_allele_minor_fraction'] == minor_count / 1000
    assert evidence['minimum_reference_effect_allele_minor_fraction'] == .90
    assert decision == expected
    if expected == 'inconclusive':
        assert evidence['direct_eaf_acceptance_reason'] == 'reference_minor_fraction_not_above_minimum'
        assert evidence['direct_eaf_accepted'] is False


def test_reference_half_frequency_counts_in_minor_fraction():
    values = _variable_minor_values(900) + [.5] + [.6] * 99
    decision, evidence = _paired_comparison(values)
    assert evidence['reference_effect_allele_minor_fraction'] == .901
    assert decision == 'eaf'


@pytest.mark.parametrize('pairs,total,expected', [(500, 1000, 'eaf'), (499, 1000, 'inconclusive'),
    (500, 1001, 'inconclusive'), (501, 1001, 'eaf')])
def test_overlap_fraction_equality_and_odd_totals(pairs, total, expected):
    decision, evidence = _paired_comparison(_variable_minor_values(pairs), total=total,
        overrides={'eaf.maf_reference_min_overlap': 2})
    assert evidence['study_variant_count'] == total
    assert evidence['observed_overlap_fraction'] == pairs / total
    assert evidence['minimum_overlap_fraction'] == .50
    assert decision == expected
    if expected == 'inconclusive':
        assert evidence['direct_eaf_acceptance_reason'] == 'insufficient_overlap_fraction'


def test_overlap_fraction_override_avoids_ceiling_roundoff():
    decision, evidence = _paired_comparison(_variable_minor_values(14), total=100,
        overrides={'eaf.maf_reference_min_overlap': 2, 'eaf.maf_reference_min_overlap_fraction': .14})
    assert evidence['observed_overlap_fraction'] == .14
    assert decision == 'eaf'


@pytest.mark.parametrize('pairs,total,expected', [(999, 999, 'inconclusive'), (1000, 1999, 'eaf')])
def test_absolute_overlap_floor_is_retained(pairs, total, expected):
    decision, evidence = _paired_comparison(_variable_minor_values(pairs), total=total)
    assert evidence['minimum_overlap'] == 1000
    assert decision == expected
    if expected == 'inconclusive':
        assert evidence['direct_eaf_acceptance_reason'] == 'insufficient_overlap'


@pytest.mark.parametrize('total', [None, 0, -1, True, 2.5, float('nan')])
def test_whole_study_comparison_requires_positive_total(total):
    from postgwas.modules.harmonisation.allele_frequency import _compare_frequency_pairs
    frame = pl.DataFrame({'s': [.1, .2], 'r': [.1, .2]})
    with pytest.raises(StudyAlleleFrequencyError, match='positive.*study.*count|study.*count.*positive'):
        _compare_frequency_pairs(frame, 's', 'r', default_policies(), NullStepContext(),
            accept_direct_eaf=True, study_variant_count=total)


def test_whole_study_overlap_uses_all_input_rows(tmp_path):
    # Six usable pairs over two chromosomes, plus missing/invalid, palindrome,
    # unmatched, split-site and ambiguously oriented diagnostic exclusions.
    reference = tmp_path / 'coverage.tsv'
    reference.write_text('CHROM\tPOS\tREF\tALT\tEUR\n' + ''.join(
        f'{chrom}\t{pos}\tG\tA\t{af}\n' for chrom, pos, af in
        [('1', 1, .1), ('1', 2, .2), ('1', 3, .3), ('2', 1, .15), ('2', 2, .25), ('2', 3, .35),
         ('1', 4, .1), ('1', 5, .2), ('1', 7, .2), ('1', 7, .1)])
        + '1\t6\tT\tA\t0.2\n1\t8\tC\tG\t0.2\n')
    frame = pl.DataFrame({'CHR': ['1'] * 3 + ['2'] * 3 + ['1'] * 6,
        'POS': [1, 2, 3, 1, 2, 3, 4, 5, 6, 99, 7, 8],
        'EA': ['A'] * 11 + ['G'], 'OA': ['G'] * 8 + ['T', 'G', 'G', 'C'],
        'AF': [.1, .2, .3, .15, .25, .35, None, -1., .2, .2, .2, .2]})
    # Different ALT at the same position makes site 7 ineligible.
    text = reference.read_text().replace('1\t7\tG\tA\t0.1', '1\t7\tG\tC\t0.1')
    reference.write_text(text)
    frame = frame.with_row_index(SOURCE_INPUT_ROW_COLUMN)
    pol = default_policies().with_overrides({'eaf.maf_reference_min_overlap': 2})
    result, prepared = prepare(frame, pol, reference)
    evidence = result['eaf_reference_validation']
    assert evidence['study_variant_count'] == 12
    assert evidence['comparable_variants'] == 6
    assert evidence['observed_overlap_fraction'] == .5
    assert prepared.select(frame.columns).equals(frame)
    bad = frame.with_columns(pl.when(pl.col(SOURCE_INPUT_ROW_COLUMN) == 0).then(None).otherwise(pl.col('AF')).alias('AF'))
    with pytest.raises(StudyAlleleFrequencyError, match='insufficient_overlap_fraction') as caught:
        prepare(bad, pol, reference)
    assert caught.value.study_decisions['eaf_reference_validation']['study_variant_count'] == 12
    assert caught.value.study_decisions['eaf_reference_validation']['comparable_variants'] == 5


def test_whole_study_pair_count_cannot_exceed_denominator():
    with pytest.raises(StudyAlleleFrequencyError, match='exceed'):
        _paired_comparison([.1, .2, .3], total=2)


@pytest.mark.parametrize('scope', ['all', 'common'])
@pytest.mark.parametrize('style', ['full', 'minimal', 'values'])
def test_whole_study_af_policies_export_and_reload(tmp_path, scope, style):
    expected = {'maf_reference_min_overlap': 1000, 'maf_reference_min_overlap_fraction': .50,
        'maf_reference_min_correlation': .90, 'maf_reference_max_mean_absolute_difference': .10,
        'maf_reference_error_margin': .02, 'study_reference_min_minor_fraction': .90}
    text = render_module_configuration('harmonisation', scope=scope, style=style)
    target = tmp_path / 'defaults.yaml'
    target.write_text(text)
    exported = yaml.safe_load(text)
    eaf = (exported['modules']['harmonisation'] if scope == 'common' else exported)['policies']['eaf']
    for key, value in expected.items():
        assert eaf[key] == value
    reloaded = load_run_configuration_for_module('harmonisation', target)
    policies = load_policies(reloaded.modules.harmonisation.policies)
    assert all(policies.get('eaf.' + key) == value for key, value in expected.items())
    source = tmp_path / 'overrides.yaml'
    source.write_text(yaml.safe_dump({'policies': {'eaf': {
        'maf_reference_min_overlap_fraction': .14, 'study_reference_min_minor_fraction': .85,
        'maf_reference_min_correlation': .93, 'maf_reference_error_margin': .02}}}))
    target.write_text(render_module_configuration('harmonisation', config_file=source, scope=scope, style=style))
    policies = load_policies(load_run_configuration_for_module('harmonisation', target).modules.harmonisation.policies)
    assert policies.get('eaf.maf_reference_min_overlap_fraction') == .14
    assert policies.get('eaf.study_reference_min_minor_fraction') == .85
    assert policies.get('eaf.maf_reference_min_correlation') == .93
    assert policies.get('eaf.maf_reference_error_margin') == .02


@pytest.mark.parametrize('key', ['eaf.maf_reference_min_overlap_fraction', 'eaf.study_reference_min_minor_fraction',
    'eaf.maf_reference_min_correlation'])
@pytest.mark.parametrize('value', [-.01, 1.01, float('nan'), float('inf'), -float('inf'), True, False, 90])
def test_af_proportion_policies_reject_invalid_values(key, value):
    with pytest.raises(PolicyError, match=key):
        default_policies().with_overrides({key: value})


@pytest.mark.parametrize('value', [.5, 1000.5, True, 1])
def test_af_absolute_overlap_keeps_integer_semantics(value):
    with pytest.raises(PolicyError, match='eaf.maf_reference_min_overlap'):
        default_policies().with_overrides({'eaf.maf_reference_min_overlap': value})


@pytest.mark.parametrize('method', ['pearson', 'spearman'])
@pytest.mark.parametrize('constant_column', ['study', 'reference'])
@pytest.mark.parametrize('folded_only', [False, True])
def test_constant_pair_vectors_have_unavailable_correlations(method, constant_column, folded_only):
    constant = [.25, .75] * 500 if folded_only else [.2] * 1000
    variable = [.1 + .3 * i / 999 for i in range(1000)]
    frame = pl.DataFrame({constant_column: constant,
        ('reference' if constant_column == 'study' else 'study'): variable})
    policies = af.resolve_policies(None).with_overrides({'eaf.maf_reference_correlation_method': method})
    _, evidence = af._compare_frequency_pairs(frame, 'study', 'reference', policies, NullStepContext())
    assert evidence['maf_correlation'] is None
    if not folded_only:
        assert evidence['direct_eaf_correlation'] is None


def test_fraction_gate_alone_cannot_select_replacement(tmp_path):
    reference = tmp_path / 'minor_fraction.tsv'
    reference.write_text('CHROM\tPOS\tREF\tALT\tEUR\n' + ''.join(
        f'1\t{i}\tG\tA\t{value}\n' for i, value in enumerate([.1, .2, .3, .4, .51], 1)))
    frame = pl.DataFrame({'CHR': ['1'] * 5, 'POS': list(range(1, 6)), 'EA': ['A'] * 5,
        'OA': ['G'] * 5, 'AF': [.1, .2, .3, .4, .49]}).with_row_index(SOURCE_INPUT_ROW_COLUMN)
    original = frame.clone()
    policies = default_policies().with_overrides({'eaf.maf_reference_min_overlap': 2})
    with pytest.raises(af.StudyAlleleFrequencyError, match='reference_minor_fraction_not_above_minimum') as caught:
        prepare(frame, policies, reference, eaf_replacement_file='/unused.tsv', eaf_replacement_column='EUR')
    evidence = caught.value.study_decisions['eaf_reference_validation']
    assert evidence['maf_reference_decision'] == 'inconclusive'
    assert evidence['reference_effect_allele_minor_fraction'] == .8
    assert evidence['direct_eaf_correlation'] > .90
    assert evidence['mean_eaf_error'] < .10
    assert 'eaf_replacement' not in caught.value.study_decisions
    assert frame.equals(original)


def test_insufficient_coverage_cannot_select_supported_maf_replacement(tmp_path):
    reference, frame, policies = fixture(tmp_path)
    frame = pl.concat([frame, frame.with_columns((pl.col('POS') + 10000).alias('POS'),
        (pl.col(SOURCE_INPUT_ROW_COLUMN) + 5).alias(SOURCE_INPUT_ROW_COLUMN))])
    policies = policies.with_overrides({'eaf.maf_reference_min_overlap_fraction': .51})
    with pytest.raises(af.StudyAlleleFrequencyError, match='insufficient_overlap_fraction') as caught:
        prepare(frame, policies, reference, eaf_replacement_file=str(reference), eaf_replacement_column='EUR')
    evidence = caught.value.study_decisions['eaf_reference_validation']
    assert evidence['observed_overlap_fraction'] == .5
    assert evidence['maf_reference_decision'] == 'inconclusive'
    assert 'eaf_replacement' not in caught.value.study_decisions


@pytest.mark.parametrize('status,reason', [('passed', 'allele_aligned_frequency_agreement'),
    ('failed', 'insufficient_overlap_fraction'), ('failed', 'correlation_below_minimum'),
    ('failed', 'reference_minor_fraction_not_above_minimum'), ('replacement_selected', 'maf_interpretation_supported')])
def test_frequency_report_shows_whole_study_gates_and_exact_operators(decisions, status, reason):
    decisions['eaf_reference_validation'].update(status=status, decision_reason=reason,
        direct_eaf_accepted=status == 'passed', study_variant_count=2000,
        observed_overlap_fraction=.55, minimum_overlap_fraction=.50,
        minimum_reference_effect_allele_minor_fraction=.90)
    text = normal(render(decisions))
    assert '1,100 / 2,000 (55.00%); required >= 50.00%' in text
    assert '0.98; required >= 0.9' in text
    assert '99.00%; required > 90.00%' in text
    assert 'required >= 0.8' not in text
    if status == 'failed':
        assert reason.replace('_', ' ') in text
        assert 'Study EAF retained' not in text


@pytest.mark.parametrize('key,value', [
    ('maf_reference_min_overlap_fraction', .60),
    ('study_reference_min_minor_fraction', .92),
    ('maf_reference_min_correlation', .93),
])
def test_af_policy_changes_scientific_checkpoint_identity_and_restarts(tmp_path, key, value):
    from postgwas.core.completion import configuration_digest
    from postgwas.core.direct_execution import _checkpoint_configuration, run_direct_with_checkpoint

    output = tmp_path / 'results'
    run = tmp_path / 'run.yaml'
    original = load_run_configuration_for_module('harmonisation')
    run.write_text(yaml.safe_dump({'policies': {'eaf': {key: value}}}))
    changed = load_run_configuration_for_module('harmonisation', run)
    # Empty arguments isolate resolved policy identity from the YAML file or CLI text.
    first = _checkpoint_configuration(original, 'harmonisation', 'harmonisation', ())
    second = _checkpoint_configuration(changed, 'harmonisation', 'harmonisation', ())
    assert configuration_digest(first) != configuration_digest(second)
    assert load_policies(changed.modules.harmonisation.policies).get('eaf.' + key) == value
    calls = []

    def invoke(config):
        def operation():
            observed = load_policies(config.modules.harmonisation.policies).get('eaf.' + key)
            calls.append(observed)
            output.mkdir(exist_ok=True)
            (output / 'scientific_output.txt').write_text(str(observed))
            return 0
        return run_direct_with_checkpoint(operation, module_name='harmonisation',
            public_command='harmonisation', arguments=(), configuration=config,
            output_directory=output)

    default = load_policies(original.modules.harmonisation.policies).get('eaf.' + key)
    assert original.run.resume is True and original.run.overwrite is False
    assert invoke(original) == invoke(original) == 0
    assert calls == [default]
    assert invoke(changed) == invoke(changed) == 0
    assert calls == [default, value]
    assert (output / 'scientific_output.txt').read_text() == str(value)
    assert changed.run.overwrite is False


def test_ordinary_folded_classifier_shares_correlation_but_not_whole_study_fraction_gates():
    import math
    reference = [.6 + .3 * i / 999 for i in range(1000)]
    supplied = [1 - value + .08 * math.sin(i) for i, value in enumerate(reference)]
    overrides = {'eaf.maf_reference_min_overlap_fraction': 1.,
        'eaf.study_reference_min_minor_fraction': 1.}
    decision, evidence = _paired_comparison(supplied, reference, total=10000,
        whole_study=False, overrides=overrides)
    assert .80 < evidence['maf_correlation'] < .90
    assert evidence['minimum_maf_correlation'] == .90
    assert evidence['minimum_error_margin'] == .02
    assert evidence['decision_reason'] == 'maf_correlation_below_minimum'
    assert decision == 'inconclusive'
    permissive, permissive_evidence = _paired_comparison(supplied, reference, total=10000,
        whole_study=False, overrides={**overrides, 'eaf.maf_reference_min_correlation': .80})
    assert permissive == 'maf'
    assert permissive_evidence['reference_suitability_status'] == 'suitable'
    assert permissive_evidence['reference_effect_allele_minor_fraction'] == 0
    assert permissive_evidence['absolute_error_difference'] > .02
    boundary, boundary_evidence = _paired_comparison(supplied, reference,
        whole_study=False, overrides={**overrides,
            'eaf.maf_reference_min_correlation': evidence['maf_correlation']})
    assert boundary_evidence['maf_correlation'] == boundary_evidence['minimum_maf_correlation']
    assert boundary == 'maf'
    for key in ('study_variant_count', 'minimum_overlap_fraction',
                'minimum_reference_effect_allele_minor_fraction'):
        assert key not in evidence and key not in permissive_evidence
    assert default_policies().get('eaf.external_min_match_fraction') == .80


@pytest.mark.parametrize('major_reference,expected,error', [(.625, 'inconclusive', .02), (.75, 'maf', .03)])
def test_ordinary_directional_margin_is_strict_and_remains_point_zero_two(major_reference, expected, error):
    # Two major-reference rows give total directional error .5/25=.02 or .75/25=.03.
    reference = _variable_minor_values(23) + [.625, major_reference]
    supplied = [min(value, 1 - value) for value in reference]
    decision, evidence = _paired_comparison(supplied, reference, whole_study=False,
        overrides={'eaf.maf_reference_min_overlap': 2})
    assert evidence['mean_maf_error'] == 0
    assert evidence['mean_eaf_error'] == error
    assert evidence['minimum_error_margin'] == .02
    assert evidence['reference_suitability_status'] == 'suitable'
    assert decision == expected


def test_true_maf_control_can_pass_whole_study_acceptance_heuristic():
    # This counterexample documents a heuristic limitation, not semantic identification:
    # every supplied value is true folded MAF, despite passing all configured gates.
    reference = [.01 + .48 * i / 9399 for i in range(9400)]
    reference += [.55 + .10 * i / 599 for i in range(600)]
    supplied = [min(value, 1 - value) for value in reference]
    decision, evidence = _paired_comparison(supplied, reference)
    assert len(supplied) == 10000
    assert evidence['direct_eaf_correlation'] == pytest.approx(.9523683946, abs=1e-10)
    assert evidence['reference_effect_allele_minor_fraction'] == .94
    assert evidence['mean_eaf_error'] == pytest.approx(.012, abs=1e-14)
    assert evidence['observed_overlap_fraction'] == 1.
    assert evidence['comparable_variants'] >= evidence['minimum_overlap']
    assert decision == 'eaf'
    assert evidence['direct_eaf_accepted'] is True
    assert evidence['reference_suitability_status'] == 'not_required'


def test_accepted_whole_study_matches_reused_and_worker_swaps_eaf_once(tmp_path, monkeypatch):
    from postgwas.config import load_configuration
    from postgwas.modules.harmonisation import strand

    reference = tmp_path / 'reference.tsv'
    reference.write_text('CHROM\tPOS\tREF\tALT\tEUR\n'
        '1\t1\tG\tA\t0.1\n1\t2\tA\tG\t0.8\n'
        '2\t1\tG\tA\t0.3\n2\t2\tA\tG\t0.6\n')
    original = pl.DataFrame({'CHR': ['1', '1', '2', '2'], 'POS': [1, 2, 1, 2],
        'EA': ['A'] * 4, 'OA': ['G'] * 4, 'AF': [.1, .2, .3, .4],
        'BETA': [.2] * 4, 'Z': [2.] * 4}).with_row_index(SOURCE_INPUT_ROW_COLUMN)
    policies = default_policies().with_overrides({'eaf.maf_reference_min_overlap': 2})
    decision, prepared = prepare(original, policies, reference)
    assert decision['eaf_reference_validation']['status'] == 'passed'
    assert decision['eaf_is_maf_initial'] is True and decision['eaf_is_maf'] is False
    assert prepared.select(original.columns).equals(original)
    decision.update(strand='forward', effect_type='beta')
    path = tmp_path / 'prepared.parquet'
    prepared.write_parquet(path)
    work = pl.read_parquet(path)

    def forbidden(*args, **kwargs):
        raise AssertionError('Accepted study must reuse prepared matches and AF decision')

    monkeypatch.setattr(strand, '_load_reference', forbidden)
    monkeypatch.setattr(strand, '_reference_orientation_candidates', forbidden)
    monkeypatch.setattr(af, '_compare_frequency_pairs', forbidden)
    layout = load_configuration().modules.harmonisation.output_layout.model_dump()
    for chrom in ('1', '2'):
        result, qc, columns = af.harmonise_allele_frequency(chrom,
            work.filter(pl.col('CHR') == chrom), dict(COLUMNS, beta_or_col='BETA', imp_z_col='Z',
                output_folder=str(tmp_path), gwas_outputname='study'), layout, '\t',
            external_eaf_colmap=MAPPING, default_eaf_file=str(reference),
            default_eaf_column='EUR', default_eaf_colmap=MAPPING,
            policies=policies, study_decision=decision)
        result = result.sort('POS')
        assert qc['strand_orientation']['reference_matches_reused'] is True
        assert qc['eaf_provenance'] == 'study_supplied'
        assert qc['final_status'] == 'success' and result.height == 2
        assert result[columns['eaf_col']].to_list() == pytest.approx(
            [.1, .8] if chrom == '1' else [.3, .6])
        assert result['BETA'].to_list() == [.2, -.2]
        assert result['Z'].to_list() == [2., -2.]
    assert original['AF'].to_list() == [.1, .2, .3, .4]


def _public_config_cli(*arguments):
    import subprocess
    import sys
    return subprocess.run([sys.executable, '-B', '-m', 'postgwas', 'config', *arguments],
        capture_output=True, text=True, check=False)


@pytest.mark.parametrize('scope,module_arguments', [('common', []), ('all', ['--module', 'harmonisation'])])
def test_public_config_validate_rejects_both_invalid_af_proportions(tmp_path, scope, module_arguments):
    target = tmp_path / 'invalid.yaml'
    exported = _public_config_cli('export', '--module', 'harmonisation', '--scope', scope,
        '--style', 'values', '--output', str(target))
    assert exported.returncode == 0, exported.stderr
    document = yaml.safe_load(target.read_text())
    module = document['modules']['harmonisation'] if scope == 'common' else document
    module['policies']['eaf'].update(maf_reference_min_overlap_fraction=90,
        study_reference_min_minor_fraction=True)
    target.write_text(yaml.safe_dump(document))
    original = target.read_bytes()
    result = _public_config_cli('validate', '--config', str(target), *module_arguments)
    assert result.returncode == 2
    assert 'Traceback' not in result.stderr
    assert 'Configuration is valid.' not in result.stdout
    assert 'eaf.maf_reference_min_overlap_fraction' in result.stderr
    assert 'eaf.study_reference_min_minor_fraction' in result.stderr
    assert '2 problems with the policy settings' in result.stderr
    assert target.read_bytes() == original


@pytest.mark.parametrize('scope,module_arguments', [('common', []), ('all', ['--module', 'harmonisation'])])
@pytest.mark.parametrize('dotted', [False, True])
def test_public_config_validate_preserves_valid_af_override_shape_and_exports(tmp_path, scope, module_arguments, dotted):
    from copy import deepcopy
    from postgwas.config import load_configuration
    overrides = {'maf_reference_min_overlap_fraction': '0.14', 'study_reference_min_minor_fraction': '0.85',
        'maf_reference_min_correlation': '0.93'}
    policies = {'eaf.' + key: value for key, value in overrides.items()} if dotted else {'eaf': overrides}
    source = tmp_path / 'valid.yaml'
    module = {'policies': deepcopy(policies)}
    source.write_text(yaml.safe_dump({'modules': {'harmonisation': module}} if scope == 'common' else module))
    original = source.read_bytes()
    result = _public_config_cli('validate', '--config', str(source), *module_arguments)
    assert result.returncode == 0, result.stderr
    assert 'Configuration is valid.' in result.stdout
    config = load_run_configuration_for_module('harmonisation', source)
    assert config.modules.harmonisation.policies == policies
    assert config.modules.harmonisation.model_dump()['policies'] == policies
    assert load_configuration().modules.harmonisation.policies == {}
    destination = tmp_path / 'export.yaml'
    exported = _public_config_cli('export', '--module', 'harmonisation', '--scope', scope,
        '--style', 'values', '--run-config', str(source), '--output', str(destination))
    assert exported.returncode == 0, exported.stderr
    reloaded = load_run_configuration_for_module('harmonisation', destination)
    assert load_policies(reloaded.modules.harmonisation.policies).resolved_dict() == load_policies(policies).resolved_dict()
    assert source.read_bytes() == original


def test_public_config_validation_does_not_require_polars():
    import subprocess
    import sys
    script = '''
import importlib.abc
import runpy
import sys
class NoPolars(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] == 'polars':
            raise AssertionError('Configuration validation must not import Polars')
sys.meta_path.insert(0, NoPolars())
sys.argv = ['postgwas', 'config', 'validate']
runpy.run_module('postgwas', run_name='__main__')
'''
    result = subprocess.run([sys.executable, '-B', '-c', script],
        capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert 'Configuration is valid.' in result.stdout
